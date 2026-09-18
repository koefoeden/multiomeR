// Kenward-Roger correction specialized to one donor random slope and one
// fixed-slope contrast. Equations follow pbkrtest 0.5.5 vcovAdj/.KR_adjust.
// Accepts independent fitted variance components; all batch inputs are immutable.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
Rcpp::DataFrame peak_gene_KR_batch_cpp(
    const Eigen::Map<MatrixXd> X, const Eigen::Map<MatrixXd> Z,
    Rcpp::List covariance, Rcpp::NumericVector coefficient,
    Rcpp::NumericVector slope_variance, Rcpp::NumericVector residual_variance) {
  const int n = X.rows(), p = X.cols(), m = covariance.size();
  if (n != Z.rows() || n <= p || p < 1 || !X.allFinite() || !Z.allFinite() || coefficient.size() != m ||
      slope_variance.size() != m || residual_variance.size() != m)
    Rcpp::stop("Incompatible batch dimensions");
  const MatrixXd G = Z * Z.transpose();
  const MatrixXd identity = MatrixXd::Identity(n, n);
  Rcpp::NumericVector pv(m, NA_REAL), df(m, NA_REAL), variance(m, NA_REAL);
  Rcpp::CharacterVector diagnostic(m, "");
  for (int k = 0; k < m; ++k) {
    Rcpp::checkUserInterrupt();
    try {
      MatrixXd Phi = Rcpp::as<MatrixXd>(covariance[k]);
      if (Phi.rows() != p || Phi.cols() != p || !Phi.allFinite() ||
          !std::isfinite(coefficient[k]) || !std::isfinite(slope_variance[k]) ||
          slope_variance[k] < 0 || !std::isfinite(residual_variance[k]) || residual_variance[k] <= 0)
        throw std::runtime_error("Invalid fitted covariance or variance");
      MatrixXd Sigma = slope_variance[k] * G + residual_variance[k] * identity;
      Eigen::LLT<MatrixXd> chol(Sigma);
      if (chol.info() != Eigen::Success) throw std::runtime_error("Singular observation covariance");
      MatrixXd inv = chol.solve(identity), T = inv * X;
      MatrixXd H[2] = {G * inv, inv};
      MatrixXd O[2] = {H[0] * X, T};
      MatrixXd P[2], Q[2][2];
      for (int i = 0; i < 2; ++i) {
        P[i] = -O[i].transpose() * T;
        P[i] = P[i].selfadjointView<Eigen::Upper>();
        for (int j = i; j < 2; ++j) Q[i][j] = O[i].transpose() * inv * O[j];
      }
      MatrixXd information(2, 2);
      for (int i = 0; i < 2; ++i) for (int j = i; j < 2; ++j) {
        information(i,j) = information(j,i) = (H[i].transpose().cwiseProduct(H[j])).sum()
          - 2 * Phi.cwiseProduct(Q[i][j]).sum()
          + (Phi * P[i]).cwiseProduct(P[j] * Phi).sum();
      }
      Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(information);
      if (eigen.info() != Eigen::Success) throw std::runtime_error("Variance information eigendecomposition failed");
      MatrixXd W;
      if (eigen.eigenvalues().cwiseAbs().minCoeff() > 1e-10) {
        Eigen::PartialPivLU<MatrixXd> lu(information);
        if (lu.rcond() < std::numeric_limits<double>::epsilon())
          throw std::runtime_error("Numerically singular variance information");
        W = 2 * lu.inverse();
      }
      else {
        // Match MASS::ginv's default relative singular-value tolerance.
        Eigen::JacobiSVD<MatrixXd> svd(information, Eigen::ComputeFullU | Eigen::ComputeFullV);
        VectorXd values = svd.singularValues();
        const double tolerance = std::sqrt(std::numeric_limits<double>::epsilon()) * values.maxCoeff();
        for (int i = 0; i < values.size(); ++i) values[i] = values[i] > tolerance ? 1 / values[i] : 0;
        W = 2 * svd.matrixV() * values.asDiagonal() * svd.matrixU().transpose();
      }
      MatrixXd U = W(0,1) * (Q[0][1] - P[0] * Phi * P[1]);
      U = (U + U.transpose()).eval();
      for (int i = 0; i < 2; ++i) U += W(i,i) * (Q[i][i] - P[i] * Phi * P[i]);
      MatrixXd adjusted = Phi + 2 * Phi * U * Phi;
      const int s = p - 1;
      MatrixXd theta = MatrixXd::Zero(p,p); theta(s,s) = 1 / Phi(s,s);
      MatrixXd ui[2] = {theta * Phi * P[0] * Phi, theta * Phi * P[1] * Phi};
      double A1 = 0, A2 = 0;
      for (int i = 0; i < 2; ++i) for (int j = i; j < 2; ++j) {
        double weight = (i == j ? 1 : 2) * W(i,j);
        A1 += weight * ui[i].trace() * ui[j].trace();
        A2 += weight * ui[i].cwiseProduct(ui[j].transpose()).sum();
      }
      double B = (A1 + 6 * A2) / 2, g = (2 * A1 - 5 * A2) / (3 * A2);
      double c1 = g / (3 + 2 * (1-g)), c2 = (1-g) / (3 + 2 * (1-g)), c3 = (3-g) / (3 + 2 * (1-g));
      double V0 = 1+c1*B, V1 = 1-c2*B, V2 = 1-c3*B;
      if (std::abs(V0) < 1e-10) V0 = 0;
      double ratio = (1-A2) / V1;
      double rho = ratio * ratio * V0 / V2;
      double ddf = 4 + 3 / (rho-1);
      double scaling = std::abs(ddf-2) < 0.01 ? 1 : ddf * (1-A2) / (ddf-2);
      double F = scaling * coefficient[k] * coefficient[k] / adjusted(s,s);
      if (!std::isfinite(ddf) || ddf <= 0 || !std::isfinite(F) || F < 0 || adjusted(s,s) <= 0)
        throw std::runtime_error("Invalid Kenward-Roger approximation");
      df[k] = ddf; variance[k] = adjusted(s,s); pv[k] = R::pf(F, 1, ddf, false, false);
    } catch (const std::exception& e) { diagnostic[k] = e.what(); }
  }
  return Rcpp::DataFrame::create(Rcpp::Named("p")=pv, Rcpp::Named("df")=df,
    Rcpp::Named("adjusted_variance")=variance, Rcpp::Named("diagnostic")=diagnostic);
}
