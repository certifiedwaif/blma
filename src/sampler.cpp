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
#include "cva.h"

// #define DEBUG

using namespace std;
using namespace Rcpp;


VectorXd gamma_to_row(const dbitset& gamma)
{
  auto p = gamma.size();
  VectorXd v(p);
  for (auto j = 0; j < p; j++) {
    v(j) = gamma[j] ? 1. : 0.;
  }
  return v;
}


//' sampler
//'
//' @param iterations
//' @param vy_in
//' @param mX_in
//' @param prior
//' @param modelprior
//' @param modelpriorvec_in
//' @param cores
//' @return The object returned is a list containing:
//' @export
// [[Rcpp::export]]
List sampler(const int iterations,
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
  mGamma.row(0) = gamma_to_row(gamma);

  // Initialise mXTX_inv and sigma2
  auto p_gamma = gamma.count();
  if (p_gamma == 0) {
    stringstream ss;
    ss << "gamma has no bits set" << endl;
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

  // Generate sample gammas
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
      bool bUpdate = !gamma[j];

      #ifdef DEBUG
      if (bUpdate) {
        Rcpp::Rcout << "Updating " << j << std::endl;
      } else {
        Rcpp::Rcout << "Downdating " << j << std::endl;
      }
      #endif

      // Update or downdate mXTX_inv
      MatrixXd mXTX_inv_prime(p_gamma_prime, p_gamma_prime);
      calculate_mXTX_inv_prime(gamma, gamma_prime, j, mXTX, mXTX_inv, mXTX_inv_prime, bUpdate);
      // MatrixXd mX_gamma_prime(n, p_gamma_prime);
      // get_cols(mX, gamma_prime, mX_gamma_prime);
      // MatrixXd mX_gamma_prime_Ty(p_gamma_prime, 1);
      // get_rows(mXTy, gamma_prime, mX_gamma_prime_Ty);
      // mXTX_inv_prime = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();

      // Calculate sigma2_prime
      double sigma2_prime = calculate_sigma2_prime(n, p_gamma_prime, mX, gamma_prime, vy, mXTX_inv_prime);
      // double sigma2_prime = 1. - (mX_gamma_prime_Ty.transpose() * mXTX_inv_prime * mX_gamma_prime_Ty).value() / n;

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
        log_p_0 = log_p_gamma_prime;
        log_p_1 = log_p_gamma;
      }
      double r = 1. / (1. + exp(log_p_0 - log_p_1));
      #ifdef DEBUG
        // Do the probabilities sum to 1?
        Rcpp::Rcout << "r " << r << std::endl;
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
    mGamma.row(i + 1) = gamma_to_row(gamma);
  }

  return List::create(Named("mGamma") = mGamma);
}
