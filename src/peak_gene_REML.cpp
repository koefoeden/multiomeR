// Profiled REML fit of one random donor slope and its Kenward-Roger test, in the donor-rank subspace.
// Z has at most one nonzero per row (the aggregate's donor column), so Z'Z is diagonal
// and Z'X, Z'Y are donor sums. The residualized random-effect design W = (I - H)Z
// enters only through its q-by-q Gram matrix W'W = Z'Z - Z'X (X'X)^-1 X'Z, whose
// eigenvalues are the nonzero eigenvalues of the residual covariance K'ZZ'K.
// The fitted fixed-effect covariance Phi enters the Kenward-Roger correction only through
// M = B Phi B', h = B Phi e_last and Phi_ll, with B = Z'X; by the Woodbury identity all three
// follow from the eigendecomposition of W'W, so no p-by-p matrix is factorized per gene.
#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Kenward-Roger correction specialized to one donor random slope and one fixed-slope
// contrast, in donor space. Equations follow pbkrtest 0.5.5 vcovAdj/.KR_adjust. Sigma =
// tau2 ZZ' + sigma2 I is block diagonal with identity-plus-rank-one blocks. Using
// Phi^-1 = X' Sigma^-1 X = (X'X - B' diag(a) B) / sigma2, every trace reduces to M, and the
// one-dimensional contrast needs only quadratic forms in phi = Phi e_last. Returns {p, df}.
static std::pair<double, double> kenward_roger(const ArrayXd& s, const MatrixXd& M, const ArrayXd& h,
    double Phi_ll, int n, int p, double coefficient, double tau2, double sigma2) {
  if (!M.allFinite() || !h.allFinite() || !std::isfinite(Phi_ll) || !std::isfinite(coefficient) ||
      !std::isfinite(tau2) || tau2 < 0 || !std::isfinite(sigma2) || sigma2 <= 0)
    throw std::runtime_error("Invalid fitted covariance or variance");
  const double s2 = sigma2, s4 = s2 * s2;
  const ArrayXd a = tau2 / (sigma2 + tau2 * s), as = a * s, d1 = (1 - as) / sigma2;
  // T'T = X' Sigma^-2 X = (X'X + B' diag(c) B) / sigma^4, and X'X = sigma2 Phi^-1 + B' diag(a) B.
  const ArrayXd e = a + (-2 * a + a * as);
  const ArrayXd g = s * (1 - as) / sigma2;
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
  // Quadratic forms in phi: P_0 phi = -B' u0, P_1 phi = -(sigma2 e_last + B' (e h)) / sigma^4.
  const ArrayXd u0 = d1.square() * h, eh = e * h;
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
  double F = scaling * coefficient * coefficient / adjusted;
  if (!std::isfinite(ddf) || ddf <= 0 || !std::isfinite(F) || F < 0 || adjusted <= 0)
    throw std::runtime_error("Invalid Kenward-Roger approximation");
  return {R::pf(F, 1, ddf, false, false), ddf};
}

static Rcpp::List fit_REML_KR_batch(const Eigen::Map<MatrixXd> X,
    const Eigen::Map<MatrixXd> Z, const Eigen::Map<MatrixXd> Y) {
  const int n=X.rows(), p=X.cols(), q=Z.cols(), d=n-p, m=Y.cols();
  if (d<1 || p<1 || Z.rows()!=n || Y.rows()!=n || !X.allFinite() || !Z.allFinite() || !Y.allFinite())
    Rcpp::stop("Invalid model dimensions or nonfinite input");
  Eigen::VectorXi column=Eigen::VectorXi::Constant(n,-1);
  VectorXd z=VectorXd::Zero(n);
  for(int k=0;k<q;++k) for(int i=0;i<n;++i) if(Z(i,k)!=0) {
    if(column[i]>=0) Rcpp::stop("Each row of Z must have at most one nonzero entry");
    column[i]=k; z[i]=Z(i,k);
  }
  auto donor_sums=[&](const MatrixXd& M) {
    MatrixXd S=MatrixXd::Zero(q,M.cols());
    for(int j=0;j<M.cols();++j) for(int i=0;i<n;++i) if(column[i]>=0) S(column[i],j)+=z[i]*M(i,j);
    return S;
  };
  VectorXd s=VectorXd::Zero(q);
  for(int i=0;i<n;++i) if(column[i]>=0) s[column[i]]+=z[i]*z[i];
  const MatrixXd B=donor_sums(X), ZtY=donor_sums(Y);
  MatrixXd XtX=MatrixXd::Zero(p,p);
  XtX.selfadjointView<Eigen::Lower>().rankUpdate(X.transpose());
  XtX=XtX.selfadjointView<Eigen::Lower>();
  Rcpp::NumericVector coefficient(m,NA_REAL),sigma2(m,NA_REAL),tau2(m,NA_REAL),ratio(m,NA_REAL),pvalue(m,NA_REAL),df(m,NA_REAL);
  Rcpp::CharacterVector diagnostic(m,"");
  auto result=[&]() {
    return Rcpp::List::create(Rcpp::Named("coefficient")=coefficient,
      Rcpp::Named("residual_variance")=sigma2,Rcpp::Named("slope_variance")=tau2,
      Rcpp::Named("variance_ratio")=ratio,Rcpp::Named("p")=pvalue,Rcpp::Named("df")=df,
      Rcpp::Named("diagnostic")=diagnostic);
  };
  Eigen::LLT<MatrixXd> xchol(XtX);
  const VectorXd Ldiag=xchol.matrixLLT().diagonal();
  // L_ii^2 / (X'X)_ii is one minus the R^2 of column i on the preceding columns; the
  // caller already screens the peak column at 1e-12. A deficient design marks this
  // peak's genes instead of stopping the branch.
  if(xchol.info()!=Eigen::Success ||
      !(Ldiag.array().square()>1e-13*XtX.diagonal().array()).all()) {
    for(int j=0;j<m;++j) diagnostic[j]="Fixed-effect design is rank deficient";
    return result();
  }
  const MatrixXd XtY=X.transpose()*Y, beta_ols=xchol.solve(XtY);
  const MatrixXd LB=xchol.matrixL().solve(B.transpose());
  const MatrixXd G=LB.transpose()*LB;
  MatrixXd gram=-G;
  gram.diagonal()+=s;
  Eigen::SelfAdjointEigenSolver<MatrixXd> eigen(gram);
  const double tol=std::numeric_limits<double>::epsilon()*n*std::max(1.0,eigen.eigenvalues().maxCoeff());
  std::vector<int> kept;
  for(int i=0;i<q;++i) if(eigen.eigenvalues()[i]>=tol) kept.push_back(i);
  const int r=kept.size();
  VectorXd eigenvalues(r);
  MatrixXd V(q,r);
  for(int i=0;i<r;++i) {eigenvalues[i]=eigen.eigenvalues()[kept[i]]; V.col(i)=eigen.eigenvectors().col(kept[i]);}
  // Coordinates of the OLS residual R in the orthonormal basis W V Lambda^-1/2; W'R = Z'R.
  const MatrixXd ZtR=ZtY-B*beta_ols;
  const MatrixXd projections=(1/eigenvalues.array().sqrt()).matrix().asDiagonal()*(V.transpose()*ZtR);
  // The same basis vectors, applied as W c = Z c - X (X'X)^-1 Z'X' c.
  const MatrixXd C=V*(1/eigenvalues.array().sqrt()).matrix().asDiagonal()*projections;
  const MatrixXd fitted_beta=beta_ols-xchol.solve(B.transpose()*C);
  // Phi / sigma2 = (X'X)^-1 + lambda (X'X)^-1 B' (I + lambda W'W)^-1 B (X'X)^-1, with
  // (I + lambda W'W)^-1 = U diag(1 / (1 + lambda Lambda)) U' over all eigenvectors U of W'W.
  const MatrixXd U=eigen.eigenvectors(), GU=G*U, XtX_inverse_Bt=xchol.solve(B.transpose());
  const ArrayXd all_eigenvalues=eigen.eigenvalues().array().max(0.0);
  const VectorXd phi0=xchol.solve(VectorXd::Unit(p,p-1)), h0=B*phi0, Uh0=U.transpose()*h0;
  for(int j=0;j<m;++j) {
    Rcpp::checkUserInterrupt();
    // The residual outside the donor subspace, computed directly to avoid cancellation.
    VectorXd outside=Y.col(j)-X*fitted_beta.col(j);
    for(int i=0;i<n;++i) if(column[i]>=0) outside[i]-=z[i]*C(column[i],j);
    const double remainder=outside.squaredNorm();
    VectorXd squares=projections.col(j).array().square();
    const double total=remainder+squares.sum();
    if(total<=1e-24*std::max(1.0,Y.col(j).squaredNorm())) {
      diagnostic[j]="No residual response variation";continue;
    }
    auto objective=[&](double t) {
      double lambda=std::expm1(t), rss=remainder, ld=0;
      for(int i=0;i<r;++i) {double v=1+lambda*eigenvalues[i];rss+=squares[i]/v;ld+=std::log(v);}
      return d*std::log(rss/d)+ld;
    };
    const int maximum=32;
    double grid[maximum+1];
    for(int i=0;i<=maximum;++i)grid[i]=objective(i);
    double best_t=0, best=grid[0];
    const double golden=0.6180339887498948482;
    auto refine=[&](double left,double right) {
      double a=right-golden*(right-left),b=left+golden*(right-left),fa=objective(a),fb=objective(b);
      for(int iter=0;iter<100 && right-left>1e-9;++iter) {
        if(fa<fb) {right=b;b=a;fb=fa;a=right-golden*(right-left);fa=objective(a);}
        else {left=a;a=b;fa=fb;b=left+golden*(right-left);fb=objective(b);}
      }
      return (left+right)/2;
    };
    for(int i=1;i<maximum;++i) if(grid[i]<=grid[i-1] && grid[i]<=grid[i+1]) {
      double t=refine(i-1,i+1),f=objective(t);
      if(f<best) {best=f;best_t=t;}
    }
    if(grid[0]<=grid[1]) {
      double t=refine(0,1),f=objective(t);
      if(f<best-1e-10) {best=f;best_t=t;}
    }
    if(grid[maximum]<best) {diagnostic[j]="Variance ratio exceeds search range";continue;}
    double lambda=std::expm1(best_t);
    double rss=remainder;
    for(int i=0;i<r;++i)rss+=squares[i]/(1+lambda*eigenvalues[i]);
    sigma2[j]=rss/d;tau2[j]=lambda*sigma2[j];ratio[j]=lambda;
    // Sigma^-1 sigma^2 = I - Z diag(w) Z' with w = lambda / (1 + lambda s).
    const ArrayXd shrinkage=1/(1+lambda*all_eigenvalues);
    const VectorXd u=lambda*(U*(shrinkage*Uh0.array()).matrix());
    const VectorXd weighted_ZtY=(lambda/(1+lambda*s.array())*ZtY.col(j).array()).matrix();
    coefficient[j]=(phi0+XtX_inverse_Bt*u).dot(XtY.col(j)-B.transpose()*weighted_ZtY);
    const MatrixXd scaled_GU=GU*shrinkage.sqrt().matrix().asDiagonal();
    const MatrixXd M=sigma2[j]*(G+lambda*scaled_GU*scaled_GU.transpose());
    try {
      const std::pair<double,double> test=kenward_roger(s.array(), M, sigma2[j]*(h0+G*u).array(),
        sigma2[j]*(phi0[p-1]+h0.dot(u)), n, p, coefficient[j], tau2[j], sigma2[j]);
      pvalue[j]=test.first; df[j]=test.second;
    } catch (const std::exception& e) { diagnostic[j]=e.what(); }
  }
  return result();
}

extern "C" SEXP peak_gene_REML_KR_batch_cpp(SEXP X, SEXP Z, SEXP Y) {
  BEGIN_RCPP
  return fit_REML_KR_batch(Rcpp::as<Eigen::Map<MatrixXd>>(X), Rcpp::as<Eigen::Map<MatrixXd>>(Z),
    Rcpp::as<Eigen::Map<MatrixXd>>(Y));
  END_RCPP
}
