#include <Rcpp.h>
#if defined(_OPENMP)
  #include <omp.h>
#endif

#include "graycode.h"
#include "correlation.h"
#include <string>

using namespace Eigen;
using namespace Rcpp;
using namespace std;

//' Perform Bayesian Linear Model Averaging over all of the possible linear models where vy is the response
//' and the covariates are in mX.
//'
//' @importFrom Rcpp evalCpp
//' @useDynLib blma
//'
//' @param vy Vector of responses
//' @param mX Covariate matrix
//' @param prior -- the choice of mixture $g$-prior used to perform Bayesian model averaging. The choices
//' available include:
//'   \itemize{
//'     \item{"BIC"}{-- the Bayesian information criterion obtained by using the cake prior
//'     of Ormerod et al. (2017).}
//'
//'     \item{"ZE"}{-- special case of the prior structure described by Maruyama and George (2011).}
//'
//'     \item{"liang_g1"}{-- the mixture \eqn{g}-prior of Liang et al. (2008) with prior hyperparameter
//'     \eqn{a=3} evaluated directly using Equation (10) of Greenaway and Ormerod (2018) where the Gaussian
//'     hypergeometric function is evaluated using the {gsl} library. Note: this option can lead to numerical problems and is only
//''    meant to be used for comparative purposes.}
//'
//'     \item{"liang_g2"}{-- the mixture \eqn{g}-prior of Liang et al. (2008) with prior hyperparameter
//'      \eqn{a=3} evaluated directly using Equation (11)  of Greenaway and Ormerod (2018).}
//'
//'     \item{"liang_g_n_appell"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008) with prior
//'      hyperparameter \eqn{a=3} evaluated using the {appell R} package.}
//'
//'     \item{"liang_g_approx"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008) with prior hyperparameter
//'      \eqn{a=3} using the approximation Equation (15)  of Greenaway and Ormerod (2018) for model with more
//'       than two covariates and numerical quadrature (see below) for models with one or two covariates.}
//'
//'     \item{"liang_g_n_quad"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008) with prior hyperparameter
//'      \eqn{a=3} evaluated using a composite trapezoid rule.}
//'
//'     \item{"robust_bayarri1"}{-- the robust prior of Bayarri et al. (2012) using default prior hyper
//'     parameter choices evaluated directly using Equation (18)  of Greenaway and Ormerod (2018) with the
//'     {gsl} library.}
//'
//'     \item{"robust_bayarri2"}{-- the robust prior of Bayarri et al. (2012) using default prior hyper
//'     parameter choices evaluated directly using Equation (19) of Greenaway and Ormerod (2018).}
//' }
//' @param modelprior The model prior to use. The choices of model prior are "uniform", "beta-binomial" or
//' "bernoulli". The choice of model prior dictates the meaning of the parameter modelpriorvec.
//' @param modelpriorvec If modelprior is "uniform", then the modelpriorvec is ignored and can be null.
//'
//' If
//' the modelprior is "beta-binomial" then modelpriorvec should be length 2 with the first element containing
//' the alpha hyperparameter for the beta prior and the second element containing the beta hyperparameter for
//' the beta prior.
//'
//' If modelprior is "bernoulli", then modelpriorvec must be of the same length as the number
//' of columns in mX. Each element i of modelpriorvec contains the prior probability of the the ith covariate
//' being included in the model.
//' @param cores The number of cores to use. Defaults to 1
//' @return A list containing
//' \describe{
//' \item{vR2}{the vector of correlations for each model}
//' \item{vp_gamma}{the vector of number of covariates for each model}
//' \item{vlogp}{the vector of logs of the likelihoods of each model}
//' \item{vinclusion_prob}{the vector of inclusion probabilities for each of the covariates}
//' }
//' @examples
//' mD <- MASS::UScrime
//' notlog <- c(2,ncol(MASS::UScrime))
//' mD[,-notlog] <- log(mD[,-notlog])
//'
//' for (j in 1:ncol(mD)) {
//'   mD[,j] <- (mD[,j] - mean(mD[,j]))/sd(mD[,j])
//' }
//'
//' varnames <- c(
//'   "log(AGE)",
//'   "S",
//'   "log(ED)",
//'   "log(Ex0)",
//'   "log(Ex1)",
//'   "log(LF)",
//'   "log(M)",
//'   "log(N)",
//'   "log(NW)",
//'   "log(U1)",
//'   "log(U2)",
//'   "log(W)",
//'   "log(X)",
//'   "log(prison)",
//'   "log(time)")
//'
//' y.t <- mD$y
//' X.f <- data.matrix(cbind(mD[1:15]))
//' colnames(X.f) <- varnames
//' blma_result <- blma(y.t, X.f, "maruyama")
//' @references
//' Bayarri, M. J., Berger, J. O., Forte, A., Garcia-Donato, G., 2012. Criteria for Bayesian
//' model choice with application to variable selection. Annals of Statistics 40 (3), 1550-
//' 1577.
//'
//' Greenaway, M. J., J. T. Ormerod (2018) Numerical aspects of Bayesian linear models averaging using mixture
//' g-priors.
//'
//' Liang, F., Paulo, R., Molina, G., Clyde, M. a., Berger, J. O., 2008. Mixtures of g priors for
//' Bayesian variable selection. Journal of the American Statistical Association 103 (481),
//' 410-423.
//'
//' Ormerod, J. T., Stewart, M., Yu, W., Romanes, S. E., 2017. Bayesian hypothesis tests
//' with diffuse priors: Can we have our cake and eat it too?
//' @export
// [[Rcpp::export]]
List blma(NumericVector vy, NumericMatrix mX, std::string prior,
          std::string modelprior = "uniform", Nullable<NumericVector> modelpriorvec = R_NilValue,
          const unsigned int cores = 1) {
  Map<VectorXd> vy_m = as< Map<VectorXd> >(vy);
  Map<MatrixXd> mX_m = as< Map<MatrixXd> >(mX);
  NumericVector modelpriorvec_r(0);
  const int intercept_col = 1;
  const bool bNatural_Order = false;
  const bool bIntercept = false;
  const bool bCentre = false;
  if (modelpriorvec.isNotNull()) {
    modelpriorvec_r = modelpriorvec.get();
  }
  Map<VectorXd> modelpriorvec_m = as< Map<VectorXd> >(modelpriorvec_r);
  List result = blma_cpp(vy_m, mX_m, prior, modelprior, modelpriorvec_m, intercept_col - 1, bNatural_Order,
                          bIntercept, bCentre, cores);
  return result;
}

//' Perform Bayesian Linear Model Averaging over all of the possible linear models where vy is the response,
//' covariates that may be included are in mZ and covariates which are always included are in mX.
//'
//' @param vy The vector of responses
//' @param mX The matrix of fixed covariates which will be included in every model
//' @param mZ The matrix of varying covariates, which may or may not be included in each model
//' @param prior -- the choice of mixture $g$-prior used to perform Bayesian model averaging. The choices
//' available include:
//'   \itemize{
//'     \item{"BIC"}{-- the Bayesian information criterion obtained by using the cake prior
//'     of Ormerod et al. (2017).}
//'
//'     \item{"ZE"}{-- special case of the prior structure described by Maruyama and George (2011).}
//'
//'     \item{"liang_g1"}{-- the mixture \eqn{g}-prior of Liang et al. (2008) with prior hyperparameter
//'     \eqn{a=3} evaluated directly using Equation (10) of Greenaway and Ormerod (2018) where the Gaussian
//'     hypergeometric function is evaluated using the {gsl} library. Note: this option can lead to numerical problems and is only
//''    meant to be used for comparative purposes.}
//'
//'     \item{"liang_g2"}{-- the mixture \eqn{g}-prior of Liang et al. (2008) with prior hyperparameter
//'      \eqn{a=3} evaluated directly using Equation (11)  of Greenaway and Ormerod (2018).}
//'
//'     \item{"liang_g_n_appell"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008) with prior
//'      hyperparameter \eqn{a=3} evaluated using the {appell R} package.}
//'
//'     \item{"liang_g_approx"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008) with prior hyperparameter
//'      \eqn{a=3} using the approximation Equation (15)  of Greenaway and Ormerod (2018) for model with more
//'       than two covariates and numerical quadrature (see below) for models with one or two covariates.}
//'
//'     \item{"liang_g_n_quad"}{-- the mixture \eqn{g/n}-prior of Liang et al. (2008) with prior hyperparameter
//'      \eqn{a=3} evaluated using a composite trapezoid rule.}
//'
//'     \item{"robust_bayarri1"}{-- the robust prior of Bayarri et al. (2012) using default prior hyper
//'     parameter choices evaluated directly using Equation (18)  of Greenaway and Ormerod (2018) with the
//'     {gsl} library.}
//'
//'     \item{"robust_bayarri2"}{-- the robust prior of Bayarri et al. (2012) using default prior hyper
//'     parameter choices evaluated directly using Equation (19) of Greenaway and Ormerod (2018).}
//' }
//' @param modelprior The model prior to use. The choices of model prior are "uniform", "beta-binomial" or
//' "bernoulli". The choice of model prior dictates the meaning of the parameter modelpriorvec.
//' @param modelpriorvec If modelprior is "uniform", then the modelpriorvec is ignored and can be null.
//'
//' If the modelprior is "beta-binomial" then modelpriorvec should be length 2 with the first element containing
//' the alpha hyperparameter for the beta prior and the second element containing the beta hyperparameter for
//' the beta prior.
//'
//' If modelprior is "bernoulli", then modelpriorvec must be of the same length as the number
//' of columns in mX. Each element i of modelpriorvec contains the prior probability of the the ith covariate
//' being included in the model.
//'
//' @param cores The number of cores to use. Defaults to 1
//'
//' @return A list containing
//' \describe{
//' \item{vR2}{the vector of correlations for each model}
//' \item{vp_gamma}{the vector of number of covariates for each model}
//' \item{vlogp}{the vector of logs of the likelihoods of each model}
//' \item{vinclusion_prob}{the vector of inclusion probabilities for each of the covariates}
//' }
//' @examples
//' mD <- MASS::UScrime
//' notlog <- c(2,ncol(MASS::UScrime))
//' mD[,-notlog] <- log(mD[,-notlog])
//'
//' for (j in 1:ncol(mD)) {
//'   mD[,j] <- (mD[,j] - mean(mD[,j]))/sd(mD[,j])
//' }
//'
//' varnames <- c(
//'   "log(AGE)",
//'   "S",
//'   "log(ED)",
//'   "log(Ex0)",
//'   "log(Ex1)",
//'   "log(LF)",
//'   "log(M)",
//'   "log(N)",
//'   "log(NW)",
//'   "log(U1)",
//'   "log(U2)",
//'   "log(W)",
//'   "log(X)",
//'   "log(prison)",
//'   "log(time)")
//'
//' y.t <- mD$y
//' X.f <- data.matrix(cbind(mD[, 1:10]))
//' colnames(X.f) <- varnames[1:10]
//' Z.f <- data.matrix(cbind(mD[, 11:15]))
//' blma_result <- blma_fixed(y.t, X.f, Z.f, "maruyama")
//' @references
//' Bayarri, M. J., Berger, J. O., Forte, A., Garcia-Donato, G., 2012. Criteria for Bayesian
//' model choice with application to variable selection. Annals of Statistics 40 (3), 1550-
//' 1577.
//'
//' Greenaway, M. J., J. T. Ormerod (2018) Numerical aspects of Bayesian linear models averaging using mixture
//' g-priors.
//'
//' Liang, F., Paulo, R., Molina, G., Clyde, M. a., Berger, J. O., 2008. Mixtures of g priors for
//' Bayesian variable selection. Journal of the American Statistical Association 103 (481),
//' 410-423.
//'
//' Ormerod, J. T., Stewart, M., Yu, W., Romanes, S. E., 2017. Bayesian hypothesis tests
//' with diffuse priors: Can we have our cake and eat it too?
//' @export
// [[Rcpp::export]]
List blma_fixed(NumericVector vy, NumericMatrix mX, NumericMatrix mZ, std::string prior,
                std::string modelprior = "uniform", Nullable<NumericVector> modelpriorvec = R_NilValue,
                const unsigned int cores = 1) {
  Map<VectorXd> vy_m = as< Map<VectorXd> >(vy);
  Map<MatrixXd> mX_m = as< Map<MatrixXd> >(mX);
  Map<MatrixXd> mZ_m = as< Map<MatrixXd> >(mZ);
  const int intercept_col = 1;
  const bool bNatural_Order = false;
  const bool bIntercept = false;
  const bool bCentre = false;
  NumericVector modelpriorvec_r(0);
  if (modelpriorvec.isNotNull()) {
    modelpriorvec_r = modelpriorvec.get();
  }
  Map<VectorXd> modelpriorvec_m = as< Map<VectorXd> >(modelpriorvec_r);
  List result = blma_fixed_cpp(vy_m, mX_m, mZ_m, prior, modelprior, modelpriorvec_m, intercept_col - 1,
                                bNatural_Order, bIntercept, bCentre, cores);
  return result;
}

//' Return the graycode matrix
//'
//' @param varying The number of covariates varying in the graycode matrix
//' @param fixed The number of fixed covariates in the graycode matrix. These covariates will always be included
//' @return The graycode matrix. The number of fixed columns will be included in the lower indexed columns
//' as 1s, while the higher indexed columns will varying depending on whether each covariate in the varying
//' set of covariates is included or not.
//' @export
// [[Rcpp::export]]
IntegerMatrix graycode(unsigned int varying, unsigned int fixed = 0) {
  Graycode gray(fixed, varying);
  MatrixXi result = gray.to_MatrixXi();
  return wrap(result);
}
