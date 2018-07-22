#include <Rcpp.h>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <sstream>
#include <unordered_map>
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#include <boost/functional/hash.hpp>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "graycode.h"
#include "correlation.h"
#include "priors.h"
#include "pva.h"
#include "sampler.h"

// #define DEBUG

using namespace std;
using namespace Rcpp;

//' sampler
//' @usage sampler_new(iterations, vy, mX, prior = "BIC", modelprior = "uniform",
//' modelpriorvec_in=NULL)
//' @param iterations The number of iterations to run the MCMC sampler for
//' @param vy_in Vector of responses
//' @param mX_in The matrix of covariates which may or may not be included in each model
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
//' @param modelpriorvec_in If modelprior is "uniform", then the modelpriorvec is ignored and can be null.
//'
//' If
//' the modelprior is "beta-binomial" then modelpriorvec should be length 2 with the first element containing
//' the alpha hyperparameter for the beta prior and the second element containing the beta hyperparameter for
//' the beta prior.
//'
//' If modelprior is "bernoulli", then modelpriorvec must be of the same length as the number
//' of columns in mX. Each element i of modelpriorvec contains the prior probability of the the ith covariate
//' being included in the model.
//' @param cores The number of cores to use. Defaults to 1.
//' @return The object returned is a list containing:
//' \itemize{
//'   \item{"mGamma"}{-- An iterations by p binary matrix containing the sampled models.}
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
//' K <- 100
//' p <- ncol(X.f)
//' sampler_result <- sampler_new(10000, y.t, X.f, prior = "BIC", modelprior = "uniform",
//'                               modelpriorvec_in=NULL)
//'
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
List sampler_new(const int iterations,
            const NumericVector vy_in,
            const NumericMatrix mX_in,
            const std::string prior,
            const std::string modelprior,
            const Nullable<NumericVector> modelpriorvec_in = R_NilValue,
            const int cores = 1L)
{
  // Try using the parallelisation in Eigen. This is an inherently serial algorithm,
  // and I don't think OpenMP is going to help us here.
  #ifdef _OPENMP
    Eigen::initParallel();
    Eigen::setNbThreads(cores);
  #endif

  VectorXd vy(vy_in.length());   // = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(vy_in);
  for (auto i = 0; i < vy_in.length(); i++) vy[i] = vy_in[i];
                                 // = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mX_in);
  MatrixXd mX(mX_in.nrow(), mX_in.ncol());
  for (auto i = 0; i < mX_in.nrow(); i++)
    for (auto j = 0; j < mX_in.ncol(); j++)
      mX(i, j) = mX_in(i, j);
  const auto n = mX.rows();
  const auto p = mX.cols();
  // Normalise vy and mX
  normalise(vy, mX);

  NumericVector modelpriorvec_r(0);
  if (!modelpriorvec_in.isNull()) {
    modelpriorvec_r = modelpriorvec_in.get();
  }
  VectorXd modelpriorvec(modelpriorvec_r.length());
  for (auto i = 0; i < modelpriorvec_r.length(); i++)
    modelpriorvec(i) = modelpriorvec_r(i);
  MatrixXd mGamma(iterations, mX_in.ncol());
  const MatrixXd mXTX = mX.transpose() * mX;
  const MatrixXd mXTy = mX.transpose() * vy;
  dbitset gamma_curr(p); // The current model
  VectorXd log_probs(iterations);
  MatrixXd mXTX_inv_curr;
  double sigma2_curr;

  log_prob_fn log_prob;
  set_log_prob(prior, log_prob);

  // Initialise gamma
  for (auto i = 0; i < p; i++) {
    if (R::runif(0., 1.) < 0.5) {
      gamma_curr[i] = false;
    } else {
      gamma_curr[i] = true;
    }
  }
  mGamma.row(0) = gamma_to_row(gamma_curr);

  // Initialise mXTX_inv and sigma2
  auto p_gamma_curr = gamma_curr.count();
  if (p_gamma_curr == 0) {
    stringstream ss;
    ss << "gamma has no bits set" << endl;
    throw domain_error(ss.str());
  }
  MatrixXd mX_gamma(n, p_gamma_curr);
  get_cols(mX, gamma_curr, mX_gamma);
  MatrixXd mX_gamma_Ty(p_gamma_curr, 1);
  get_rows(mXTy, gamma_curr, mX_gamma_Ty);
  mXTX_inv_curr = (mX_gamma.transpose() * mX_gamma).inverse();
  sigma2_curr = 1. - (mX_gamma_Ty.transpose() * mXTX_inv_curr * mX_gamma_Ty).value() / n;
  #ifdef DEBUG
  Rcpp::Rcout << "sigma2" << sigma2 << std::endl;
  #endif

  // Generate sample gammas

  double log_p_curr = calculate_log_prob(n, p, 1. - sigma2_curr, p_gamma_curr, gamma_curr, log_prob, modelprior, modelpriorvec);
  for (auto i = 0; i < iterations - 1; i++) {
    Rcpp::checkUserInterrupt();
    #ifdef DEBUG
        Rcpp::Rcout << "Iteration " << i << std::endl;
    #endif
    p_gamma_curr = gamma_curr.count();
    MatrixXd mX_gamma(n, p_gamma_curr);
    get_cols(mX, gamma_curr, mX_gamma);
    MatrixXd mX_gamma_Ty(p_gamma_curr, 1);
    get_rows(mXTy, gamma_curr, mX_gamma_Ty);
    mXTX_inv_curr = (mX_gamma.transpose() * mX_gamma).inverse();
    sigma2_curr = 1. - (mX_gamma_Ty.transpose() * mXTX_inv_curr * mX_gamma_Ty).value() / n;
    log_p_curr = calculate_log_prob(n, p, 1. - sigma2_curr, p_gamma_curr, gamma_curr, log_prob, modelprior, modelpriorvec);
    // Try to alter model covariates
    for (auto j = 0; j < p; j++) {
      dbitset gamma_prop = gamma_curr; // The next model we will be considering
      gamma_prop[j] = !gamma_curr[j];

      p_gamma_curr = gamma_curr.count();
      auto p_gamma_prop = gamma_prop.count();

      // Update or downdate mXTX_inv
      MatrixXd mXTX_inv_prop(p_gamma_prop, p_gamma_prop);
      // calculate_mXTX_inv_prop(gamma_curr, gamma_prop, j, mXTX, mXTX_inv_curr, mXTX_inv_prop, bUpdate);
      MatrixXd mX_gamma_prop(n, p_gamma_prop);
      get_cols(mX, gamma_prop, mX_gamma_prop);
      mXTX_inv_prop = (mX_gamma_prop.transpose() * mX_gamma_prop).inverse();

      // Calculate sigma2_prop
      // double sigma2_prop = calculate_sigma2_prop(n, p_gamma_prop, mX, gamma_prop, vy, mXTX_inv_prop);
      MatrixXd mX_gamma_prop_Ty = mX_gamma_prop.transpose() * vy;
      double sigma2_prop = 1. - (mX_gamma_prop_Ty.transpose() * mXTX_inv_prop * mX_gamma_prop_Ty).value() / n;

      double log_p_prop;
      log_p_prop = calculate_log_prob(n, p, 1. - sigma2_prop, p_gamma_prop, gamma_prop, log_prob, modelprior, modelpriorvec);
      double r = 1. / (1. + exp(log_p_prop - log_p_curr));
      // double r = exp(log_p_1) / (exp(log_p_0) + exp(log_p_1));
      auto prob = R::runif(0., 1.);
      if (r < prob) {
        gamma_curr = gamma_prop;
	log_p_curr = log_p_prop;
        sigma2_curr = sigma2_prop;
        mXTX_inv_curr = mXTX_inv_prop;
      }
    }
    mGamma.row(i + 1) = gamma_to_row(gamma_curr);
  }

  return List::create(Named("mGamma") = mGamma);
}
