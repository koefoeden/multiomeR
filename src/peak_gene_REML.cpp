// Profiled REML fitter for one random donor slope.
// Shared immutable spectral decomposition; independent scalar optimization per gene.
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <limits>
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
Rcpp::List peak_gene_REML_batch_cpp(const Eigen::Map<MatrixXd> X,
    const Eigen::Map<MatrixXd> Z, const Eigen::Map<MatrixXd> Y) {
  const int n=X.rows(), p=X.cols(), d=n-p, m=Y.cols();
  if (d<1 || p<1 || Z.rows()!=n || Y.rows()!=n || !X.allFinite() || !Z.allFinite() || !Y.allFinite())
    Rcpp::stop("Invalid model dimensions or nonfinite input");
  Eigen::ColPivHouseholderQR<MatrixXd> qr(X);
  if(qr.rank()!=p) Rcpp::stop("Fixed-effect design is rank deficient");
  MatrixXd Q=qr.householderQ()*MatrixXd::Identity(n,n);
  MatrixXd K=Q.rightCols(d), projected=K.transpose()*Z;
  MatrixXd S=projected*projected.transpose();
  Eigen::SelfAdjointEigenSolver<MatrixXd> spectrum(S);
  if(spectrum.info()!=Eigen::Success) Rcpp::stop("Residual covariance decomposition failed");
  VectorXd eigenvalues=spectrum.eigenvalues();
  const double tol=std::numeric_limits<double>::epsilon()*n*std::max(1.0,eigenvalues.maxCoeff());
  for(int i=0;i<d;++i) if(eigenvalues[i]<tol) eigenvalues[i]=0;
  MatrixXd residuals=spectrum.eigenvectors().transpose()*K.transpose()*Y;
  MatrixXd XtX=X.transpose()*X;
  Eigen::LLT<MatrixXd> xchol(XtX);
  MatrixXd L=xchol.matrixL();
  double logdetX=2*L.diagonal().array().log().sum();
  MatrixXd ZtZ=Z.transpose()*Z;
  Eigen::SelfAdjointEigenSolver<MatrixXd> zs(ZtZ);
  MatrixXd U=Z*zs.eigenvectors();
  VectorXd ze=zs.eigenvalues().cwiseMax(0);
  Rcpp::List covariance(m);
  Rcpp::NumericMatrix beta(p,m);
  std::fill(beta.begin(),beta.end(),NA_REAL);
  Rcpp::NumericVector sigma2(m,NA_REAL),tau2(m,NA_REAL),ratio(m,NA_REAL),deviance(m,NA_REAL);
  Rcpp::CharacterVector diagnostic(m,"");
  Rcpp::IntegerVector evaluations(m);
  for(int j=0;j<m;++j) {
    Rcpp::checkUserInterrupt();
    covariance[j]=Rcpp::NumericMatrix(p,p);
    VectorXd squares=residuals.col(j).array().square();
    if(squares.sum()<=1e-24*std::max(1.0,Y.col(j).squaredNorm())) {
      diagnostic[j]="No residual response variation";continue;
    }
    auto objective=[&](double t) {
      ++evaluations[j];
      double lambda=std::expm1(t), rss=0, ld=0;
      for(int i=0;i<d;++i) {double v=1+lambda*eigenvalues[i];rss+=squares[i]/v;ld+=std::log(v);}
      return d*std::log(rss/d)+ld;
    };
    // Scan the profile for interior minima, including the zero-variance edge.
    // The finite upper limit is explicit and never accepted as a valid optimum.
    const int maximum=32;
    double grid[maximum+1];
    for(int i=0;i<=maximum;++i)grid[i]=objective(i);
    double best_t=0, best=grid[0];
    for(int i=1;i<maximum;++i) if(grid[i]<=grid[i-1] && grid[i]<=grid[i+1]) {
      double left=i-1,right=i+1;
      const double golden=0.6180339887498948482;
      double a=right-golden*(right-left),b=left+golden*(right-left);
      double fa=objective(a),fb=objective(b);
      for(int iter=0;iter<100 && right-left>1e-9;++iter) {
        if(fa<fb) {right=b;b=a;fb=fa;a=right-golden*(right-left);fa=objective(a);}
        else {left=a;a=b;fa=fb;b=left+golden*(right-left);fb=objective(b);}
      }
      double t=(left+right)/2,f=objective(t);
      if(f<best) {best=f;best_t=t;}
    }
    // Also search a minimum between zero and the first grid point.
    if(grid[0]<=grid[1]) {
      double left=0,right=1;
      const double golden=0.6180339887498948482;
      double a=right-golden,b=golden,fa=objective(a),fb=objective(b);
      for(int iter=0;iter<100 && right-left>1e-9;++iter) {
        if(fa<fb) {right=b;b=a;fb=fa;a=right-golden*(right-left);fa=objective(a);}
        else {left=a;a=b;fa=fb;b=left+golden*(right-left);fb=objective(b);}
      }
      double t=(left+right)/2,f=objective(t);
      if(f<best-1e-10) {best=f;best_t=t;}
    }
    if(grid[maximum]<best) {diagnostic[j]="Variance ratio exceeds search range";continue;}
    double lambda=std::expm1(best_t);
    VectorXd weights=ze.unaryExpr([&](double e){return lambda/(1+lambda*e);});
    MatrixXd WX=X-U*weights.asDiagonal()*(U.transpose()*X);
    VectorXd Wy=Y.col(j)-U*weights.asDiagonal()*(U.transpose()*Y.col(j));
    MatrixXd information=X.transpose()*WX;
    Eigen::LLT<MatrixXd> chol(information);
    if(chol.info()!=Eigen::Success || information.partialPivLu().rcond()<1e-12) {
      diagnostic[j]="Ill-conditioned fitted fixed-effect information";continue;
    }
    VectorXd estimate=chol.solve(X.transpose()*Wy);
    double rss=0;
    for(int i=0;i<d;++i)rss+=squares[i]/(1+lambda*eigenvalues[i]);
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
