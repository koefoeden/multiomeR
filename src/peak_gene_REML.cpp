// Profiled REML fitter for one random donor slope, in the donor-rank subspace.
// Z has at most one nonzero per row (the aggregate's donor column), so Z'Z is diagonal
// and Z'X, Z'Y are donor sums. The residualized random-effect design W = (I - H)Z
// enters only through its q-by-q Gram matrix W'W = Z'Z - Z'X (X'X)^-1 X'Z, whose
// eigenvalues are the nonzero eigenvalues of the residual covariance K'ZZ'K.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <limits>
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
Rcpp::List peak_gene_REML_batch_cpp(const Eigen::Map<MatrixXd> X,
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
  Eigen::LLT<MatrixXd> xchol(XtX);
  const VectorXd Ldiag=xchol.matrixLLT().diagonal();
  if(xchol.info()!=Eigen::Success || Ldiag.minCoeff()<=std::sqrt(std::numeric_limits<double>::epsilon())*Ldiag.maxCoeff())
    Rcpp::stop("Fixed-effect design is rank deficient");
  const double logdetX=2*Ldiag.array().log().sum();
  const MatrixXd XtY=X.transpose()*Y, beta_ols=xchol.solve(XtY);
  const MatrixXd LB=xchol.matrixL().solve(B.transpose());
  MatrixXd gram=-LB.transpose()*LB;
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
  Rcpp::List covariance(m);
  Rcpp::NumericMatrix beta(p,m);
  std::fill(beta.begin(),beta.end(),NA_REAL);
  Rcpp::NumericVector sigma2(m,NA_REAL),tau2(m,NA_REAL),ratio(m,NA_REAL),deviance(m,NA_REAL);
  Rcpp::CharacterVector diagnostic(m,"");
  Rcpp::IntegerVector evaluations(m);
  for(int j=0;j<m;++j) {
    Rcpp::checkUserInterrupt();
    covariance[j]=Rcpp::NumericMatrix(p,p);
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
      ++evaluations[j];
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
    // Sigma^-1 sigma^2 = I - Z diag(w) Z' with w = lambda / (1 + lambda s).
    VectorXd weights=s.unaryExpr([&](double e){return lambda/(1+lambda*e);});
    MatrixXd information=XtX-B.transpose()*weights.asDiagonal()*B;
    Eigen::LLT<MatrixXd> chol(information);
    if(chol.info()!=Eigen::Success || information.partialPivLu().rcond()<1e-12) {
      diagnostic[j]="Ill-conditioned fitted fixed-effect information";continue;
    }
    VectorXd estimate=chol.solve(XtY.col(j)-B.transpose()*weights.asDiagonal()*ZtY.col(j));
    double rss=remainder;
    for(int i=0;i<r;++i)rss+=squares[i]/(1+lambda*eigenvalues[i]);
    sigma2[j]=rss/d;tau2[j]=lambda*sigma2[j];ratio[j]=lambda;
    deviance[j]=best+logdetX+d*(1+std::log(2*3.14159265358979323846));
    covariance[j]=Rcpp::wrap((sigma2[j]*chol.solve(MatrixXd::Identity(p,p))).eval());
    for(int i=0;i<p;++i)beta(i,j)=estimate[i];
  }
  return Rcpp::List::create(Rcpp::Named("beta")=beta,Rcpp::Named("covariance")=covariance,
    Rcpp::Named("residual_variance")=sigma2,Rcpp::Named("slope_variance")=tau2,
    Rcpp::Named("variance_ratio")=ratio,Rcpp::Named("REML_deviance")=deviance,
    Rcpp::Named("diagnostic")=diagnostic,Rcpp::Named("evaluations")=evaluations);
}
