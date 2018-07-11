#include <Rcpp.h>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
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
#include "cva.h"

// #define DEBUG

using namespace std;
using namespace Rcpp;


template <typename Derived1>
void gamma_to_row(const dbitset& gamma, Eigen::MatrixBase<Derived1>& row)
{
  auto p = gamma.size();
  for (auto i = 0; i < p; i++) {
    row(0, i) = gamma[i] ? 1. : 0.;
  }
}

//' mcmcmc
//'
//' @param n     Number of MCMCMC samples
//' @param vy_in Vector of responses
//' @param mX_in The matrix of covariates which may or may not be included in each model
//' @param mGamma_in Matrix of initial models, a K by p logical matrix
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
//'   \item{"mGamma"}{-- An n by p binary matrix containing the final population of models}
//'
//'   \item{"vgamma.hat"}{-- The most probable model found by cva}
//'
//'   \item{"vlogp"}{-- The null-based Bayes factor for each model in the population}
//'
//'   \item{"posterior_model_probabilities"}{-- The estimated posterior model parameters for each model in
//'   the population.}
//'
//'   \item{"posterior_inclusion_probabilities"}{-- The estimated variable inclusion probabilities for each
//'   model in the population.}
//'
//'   \item{"vR2"}{-- The fitted R-squared values for each model in the population.}
//'
//'   \item{"vp"}{-- The model size for each model in the population.}
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
//' n_iterations <- 100
//' p <- ncol(X.f)
//' initial_gamma <- matrix(rbinom(K * p, 1, .5), K, p)
//' mcmcmc_result <- mcmcmc(n_iterations, y.t, X.f, prior = "BIC", modelprior = "uniform",
//'                   modelpriorvec_in=NULL)
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
List mcmcmc(const int iterations,
            const NumericVector vy_in,
            const NumericMatrix mX_in,
            const std::string prior,
            const std::string modelprior,
            const Nullable<NumericVector> modelpriorvec_in = R_NilValue,
            const int cores = 1L)
{
  VectorXd vy(vy_in.length());   // = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(vy_in);
  for (auto i = 0; i < vy_in.length(); i++) vy[i] = vy_in[i];
                                 // = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mX_in);
  MatrixXd mX(mX_in.nrow(), mX_in.ncol());
  for (auto i = 0; i < mX_in.nrow(); i++)
    for (auto j = 0; j < mX_in.ncol(); j++)
      mX(i, j) = mX_in(i, j);
  const auto n = mX.rows();
  const auto p = mX.cols();
  vy = sqrt(n) * (vy.array() - vy.mean()) / vy.norm();
  for (auto j = 0; j < p; j++) {
    mX.col(j) = sqrt(n) * (mX.col(j).array() - mX.col(j).mean()) / mX.col(j).norm();
  }
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
  dbitset gamma(p); // The current model
  VectorXd log_probs(iterations);
  MatrixXd mXTX_inv;
  double sigma2;

  log_prob_fn log_prob;
  set_log_prob(prior, log_prob);

  // Initialise gamma
  for (auto i = 0; i < p; i++) {
    if (R::runif(0., 1.) < 0.5) {
      gamma[i] = false;
    } else {
      gamma[i] = true;
    }
  }

  // Initialise mXTX_inv and sigma2
  auto p_gamma = gamma.count();
  if (p_gamma == 0) {
    stringstream ss;
    ss << "gamma[" << k + 1 << "] has no bits set" << endl;
    throw domain_error(ss.str());
  }
  MatrixXd mX_gamma(n, p_gamma);
  get_cols(mX, gamma, mX_gamma);
  MatrixXd mX_gamma_Ty(p_gamma, 1);
  get_rows(mXTy, gamma, mX_gamma_Ty);
  mXTX_inv = (mX_gamma.transpose() * mX_gamma).inverse();
  sigma2 = 1. - (mX_gamma_Ty.transpose() * mXTX_inv * mX_gamma_Ty).value() / n;
  #ifdef DEBUG
  Rcpp::Rcout << "sigma2" << sigma2 << std::endl;
  #endif
  calculate_log_probabilities(gamma, sigma2, n, log_probs, log_prob, modelprior, modelpriorvec);

  // Generate sample gammas
  auto iteration = 1;
  #ifdef DEBUG
  Rcpp::Rcout << "Iteration " << iteration << std::endl;
  #endif

  for (auto i = 0; i < iterations - 1; i++) {
    // Try to alter model covariates
    for (auto j = 0; j < p; j++) {
      dbitset gamma_prime = gamma; // The next model we will be considering
      gamma_prime[j] = !gamma_prime[j];

      auto p_gamma = gamma.count();
      auto p_gamma_prime = gamma_prime.count();
      if ((p_gamma_prime == 0) || (p_gamma_prime >= n - 1))
        continue;
      bool bUpdate = !gamma_current[j];

      #ifdef DEBUG
      Rcpp::Rcout << bUpdate ? "Updating" : "Downdating" << j << std::endl;
      #endif

      // Update or downdate mXTX_inv
      MatrixXd mXTX_inv_prime(p_gamma_prime, p_gamma_prime);
      calculate_mXTX_inv_prime(gamma, gamma_prime, j, mXTX, mXTX_inv, mXTX_inv_prime, bUpdate);

      // Calculate sigma2_prime
      double sigma2_prime = calculate_sigma2_prime(n, p_gamma_prime, mX, gamma_prime, vy, mXTX_inv_prime);

      #ifdef DEBUG
      Rcpp::Rcout << "sigma2 " << sigma2 << std::endl;
      Rcpp::Rcout << "sigma2_prime " << sigma2_prime << std::endl;
      #endif
      double log_p_gamma;
      double log_p_gamma_prime;
      log_p_gamma = calculate_log_prob(n, p, 1. - sigma2, p_gamma, gamma, log_prob, modelprior, modelpriorvec);
      log_p_gamma_prime = calculate_log_prob(n, p, 1. - sigma2_prime, p_gamma_prime, gamma_prime, log_prob, modelprior, modelpriorvec);
      #ifdef DEBUG
      Rcpp::Rcout << "log_p_gamma " << log_p_gamma;
      Rcpp::Rcout << " log_p_gamma_prime " << log_p_gamma_prime;
      Rcpp::Rcout << " difference " << log_p_gamma_prime - log_p_gamma << std::endl;
      #endif
      double log_p_0, log_p_1;
      if (bUpdate) {
        log_p_0 = log_p_gamma;
        log_p_1 = log_p_gamma_prime;
      } else {
        log_p_0 = log_p_gamma;
        log_p_1 = log_p_gamma_prime;
      }
      double r = exp(log_p_1);
      #ifdef DEBUG
        // Do the probabilities sum to 1?
      #endif
      if (R::runif(0., 1.) < r) {
        gamma[j] = true;
        #ifdef DEBUG
        Rcpp::Rcout << "Keep update" << std::endl;
        #endif
        if (bUpdate) {
          sigma2 = sigma2_prime;
          mXTX_inv = mXTX_inv_prime;
        }
      } else {
        gamma[j] = false;
        #ifdef DEBUG
        Rcpp::Rcout << "Don't keep update" << std::endl;
        #endif
        if (!bUpdate) {
          sigma2 = sigma2_prime;
          mXTX_inv = mXTX_inv_prime;
        }
      }
    }
    gamma_to_row(gamma, mGamma.row(i + 1));
  }
  calculate_log_probabilities(gamma, sigma2, n, log_probs, log_prob, modelprior, modelpriorvec);
  
  auto M = log_probs.array().maxCoeff(); // Maximum log-likelihood
  // Rcpp::Rcout << "M " << M << std::endl;
  VectorXd vmodel_prob = (log_probs.array() - M).array().exp() / (log_probs.array() - M).array().exp().sum();
  VectorXd vinclusion_prob = mGamma_prime.transpose() * vmodel_prob;

  uint max_idx;
  auto max_prob = vmodel_prob.maxCoeff(&max_idx);
  VectorXd vgamma_hat = mGamma_prime.row(max_idx);

  List result = List::create(Named("mGamma") = mGamma_prime,
                             Named("vgamma.hat") = vgamma_hat,
                             Named("vlogp") = log_probs,
                             Named("posterior_model_probabilities") = vmodel_prob,
                             Named("posterior_inclusion_probabilities") = vinclusion_prob,
                             Named("vR2") = 1. - sigma2.array(),
                             Named("vp") = vp_gamma);
  return result;
}
