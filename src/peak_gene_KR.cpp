// Kenward-Roger correction specialized to one donor random slope and one
// fixed-slope contrast, in donor space. Equations follow pbkrtest 0.5.5
// vcovAdj/.KR_adjust. Z has one nonzero per row (the aggregate's donor column), so
// Sigma = tau2 ZZ' + sigma2 I is block diagonal with identity-plus-rank-one blocks.
// Using Phi^-1 = X' Sigma^-1 X = (X'X - B' diag(a) B) / sigma2 with B = Z'X, every
// trace reduces to M = B Phi B' (q by q), and the one-dimensional contrast needs only
// quadratic forms in phi = Phi e_last.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;

// [[Rcpp::export]]
Rcpp::DataFrame peak_gene_KR_batch_cpp(
    const Eigen::Map<MatrixXd> X, const Eigen::Map<MatrixXd> Z,
    Rcpp::List covariance, Rcpp::NumericVector coefficient,
    Rcpp::NumericVector slope_variance, Rcpp::NumericVector residual_variance) {
  const int n = X.rows(), p = X.cols(), q = Z.cols(), m = covariance.size();
  if (n != Z.rows() || n <= p || p < 1 || !X.allFinite() || !Z.allFinite() || coefficient.size() != m ||
      slope_variance.size() != m || residual_variance.size() != m)
    Rcpp::stop("Incompatible batch dimensions");
  MatrixXd B = MatrixXd::Zero(q, p);
  ArrayXd s = ArrayXd::Zero(q);
  {
    Eigen::VectorXi column = Eigen::VectorXi::Constant(n, -1);
    VectorXd z = VectorXd::Zero(n);
    for (int k = 0; k < q; ++k) for (int i = 0; i < n; ++i) if (Z(i,k) != 0) {
      if (column[i] >= 0) Rcpp::stop("Each row of Z must have at most one nonzero entry");
      column[i] = k; z[i] = Z(i,k);
    }
    for (int i = 0; i < n; ++i) if (column[i] >= 0) s[column[i]] += z[i] * z[i];
    for (int j = 0; j < p; ++j) for (int i = 0; i < n; ++i) if (column[i] >= 0) B(column[i], j) += z[i] * X(i,j);
  }
  const int last = p - 1;
  Rcpp::NumericVector pv(m, NA_REAL), df(m, NA_REAL), variance(m, NA_REAL);
  Rcpp::CharacterVector diagnostic(m, "");
  for (int k = 0; k < m; ++k) {
    Rcpp::checkUserInterrupt();
    try {
      const MatrixXd Phi = Rcpp::as<MatrixXd>(covariance[k]);
      const double tau2 = slope_variance[k], sigma2 = residual_variance[k];
      if (Phi.rows() != p || Phi.cols() != p || !Phi.allFinite() || !std::isfinite(coefficient[k]) ||
          !std::isfinite(tau2) || tau2 < 0 || !std::isfinite(sigma2) || sigma2 <= 0)
        throw std::runtime_error("Invalid fitted covariance or variance");
      const double s2 = sigma2, s4 = s2 * s2;
      const ArrayXd a = tau2 / (sigma2 + tau2 * s), as = a * s, d1 = (1 - as) / sigma2;
      // T'T = X' Sigma^-2 X = (X'X + B' diag(c) B) / sigma^4, and X'X = sigma2 Phi^-1 + B' diag(a) B.
      const ArrayXd e = a + (-2 * a + a * as);
      const ArrayXd g = s * (1 - as) / sigma2;
      const MatrixXd M = B * Phi * B.transpose();
      const ArrayXd Mdiag = M.diagonal().array();
      const MatrixXd M2 = M.array().square().matrix();
      // tr(Phi Q_ij) and tr(Phi P_i Phi P_j), with C = diag(d1) B.
      const double trQ00 = (g * d1.square() * Mdiag).sum();
      const double trQ01 = (d1.cube() * Mdiag).sum();
      const double trPhiTtT = (s2 * p + (e * Mdiag).sum()) / s4;
      const double trQ11 = (trPhiTtT - (a * d1.square() * Mdiag).sum()) / s2;
      const MatrixXd CPhiC = d1.matrix().asDiagonal() * M * d1.matrix().asDiagonal();
      const double trPP00 = CPhiC.squaredNorm();
      const double trPP01 = (s2 * (d1.square() * Mdiag).sum() +
        d1.square().matrix().dot(M2 * e.matrix())) / s4;
      const double trPP11 = (s4 * p + 2 * s2 * (e * Mdiag).sum() + e.matrix().dot(M2 * e.matrix())) / (s4 * s4);
      double trace[2][2];
      trace[1][1] = (n - (2 * as - as.square()).sum()) / s4;
      trace[0][1] = ((1 - as).square() * s).sum() / s4;
      trace[0][0] = g.matrix().squaredNorm();
      const double trQ[2][2] = {{trQ00, trQ01}, {trQ01, trQ11}};
      const double trPP[2][2] = {{trPP00, trPP01}, {trPP01, trPP11}};
      MatrixXd information(2, 2);
      for (int i = 0; i < 2; ++i) for (int j = i; j < 2; ++j)
        information(i,j) = information(j,i) = trace[i][j] - 2 * trQ[i][j] + trPP[i][j];
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
        Eigen::JacobiSVD<MatrixXd> svd(information, Eigen::ComputeFullU | Eigen::ComputeFullV);
        VectorXd values = svd.singularValues();
        const double tolerance = std::sqrt(std::numeric_limits<double>::epsilon()) * values.maxCoeff();
        for (int i = 0; i < values.size(); ++i) values[i] = values[i] > tolerance ? 1 / values[i] : 0;
        W = 2 * svd.matrixV() * values.asDiagonal() * svd.matrixU().transpose();
      }
      // Quadratic forms in phi: h = B phi, P_0 phi = -B' u0, P_1 phi = -(sigma2 e_last + B' (e h)) / sigma^4.
      const VectorXd phi = Phi.col(last);
      const double Phi_ll = phi[last];
      const ArrayXd h = (B * phi).array(), u0 = d1.square() * h, eh = e * h;
      const double phiTtTphi = (s2 * Phi_ll + (e * h.square()).sum()) / s4;
      const double phiQ[2][2] = {
        {(g * d1.square() * h.square()).sum(), (d1.cube() * h.square()).sum()},
        {(d1.cube() * h.square()).sum(), (phiTtTphi - (a * d1.square() * h.square()).sum()) / s2}};
      const VectorXd Mu0 = M * u0.matrix(), Meh = M * eh.matrix();
      const double phiPPP[2][2] = {
        {u0.matrix().dot(Mu0), (s2 * (h * u0).sum() + u0.matrix().dot(Meh)) / s4},
        {(s2 * (h * u0).sum() + u0.matrix().dot(Meh)) / s4,
         (s4 * Phi_ll + 2 * s2 * (h * eh).sum() + eh.matrix().dot(Meh)) / (s4 * s4)}};
      double phiUphi = 2 * W(0,1) * (phiQ[0][1] - phiPPP[0][1]);
      for (int i = 0; i < 2; ++i) phiUphi += W(i,i) * (phiQ[i][i] - phiPPP[i][i]);
      const double adjusted = Phi_ll + 2 * phiUphi;
      // For a one-dimensional contrast, A1 = A2 = sum_ij W_ij t_i t_j with t_i = phi' P_i phi / Phi_ll.
      const double t[2] = {-(d1.square() * h.square()).sum() / Phi_ll, -phiTtTphi / Phi_ll};
      double A1 = 0;
      for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) A1 += W(i,j) * t[i] * t[j];
      const double A2 = A1;
      double Bk = (A1 + 6 * A2) / 2, gk = (2 * A1 - 5 * A2) / (3 * A2);
      double c1 = gk / (3 + 2 * (1-gk)), c2 = (1-gk) / (3 + 2 * (1-gk)), c3 = (3-gk) / (3 + 2 * (1-gk));
      double V0 = 1+c1*Bk, V1 = 1-c2*Bk, V2 = 1-c3*Bk;
      if (std::abs(V0) < 1e-10) V0 = 0;
      double ratio = (1-A2) / V1;
      double rho = ratio * ratio * V0 / V2;
      double ddf = 4 + 3 / (rho-1);
      double scaling = std::abs(ddf-2) < 0.01 ? 1 : ddf * (1-A2) / (ddf-2);
      double F = scaling * coefficient[k] * coefficient[k] / adjusted;
      if (!std::isfinite(ddf) || ddf <= 0 || !std::isfinite(F) || F < 0 || adjusted <= 0)
        throw std::runtime_error("Invalid Kenward-Roger approximation");
      df[k] = ddf; variance[k] = adjusted; pv[k] = R::pf(F, 1, ddf, false, false);
    } catch (const std::exception& e) { diagnostic[k] = e.what(); }
  }
  return Rcpp::DataFrame::create(Rcpp::Named("p")=pv, Rcpp::Named("df")=df,
    Rcpp::Named("adjusted_variance")=variance, Rcpp::Named("diagnostic")=diagnostic);
}
