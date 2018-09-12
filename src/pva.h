// pva.h

#pragma once

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

#include "graycode.h"
#include "correlation.h"
#include "priors.h"

// #define DEBUG

using namespace std;
using namespace Rcpp;

double calculate_log_prob(const uint n, const uint p, const double R2, const uint p_gamma,
                          const dbitset& gamma,
                          const std::function<double (const int n, const int p, double vR2, int vp_gamma)> log_prob,
                          const std::string modelprior, const VectorXd& modelpriorvec);
void calculate_log_probabilities(const vector< dbitset >& gamma, const VectorXd& sigma2,
                                  const int n,
                                  VectorXd& log_probs,
                                  const std::function<double (const int n, const int p, double vR2, int vp_gamma)> log_prob,
                                  const std::string& modelprior, const VectorXd& modelpriorvec);
void calculate_weights(const VectorXd& sigma2,
                        const VectorXd& log_probs,
                        VectorXd& w);
double calculate_entropy(const VectorXd& w);
double calculate_w_dot_prob(const VectorXd& w,
                            const VectorXd& log_probs);
void gamma_to_MatrixXd(const vector< dbitset >& gamma, VectorXd& m);
void calculate_mXTX_inv_prime(const dbitset& gamma, const dbitset& gamma_prime, int j,
                              const MatrixXd& mXTX, const MatrixXd& mXTX_inv,
                              MatrixXd& mXTX_inv_prime,  bool bUpdate);
double calculate_sigma2_prime(const uint n, const uint p_gamma_prime,
                              const MatrixXd& mX, const dbitset& gamma_prime,
                              const VectorXd& vy,
                              const MatrixXd& mXTX_inv_prime);
List pva(const NumericVector vy_in, const NumericMatrix mX_in,
         const NumericMatrix mGamma_in,
         const std::string prior,
         const std::string modelprior, const Nullable<NumericVector> modelpriorvec_in,
         const bool bUnique,
         const double lambda,
         const int cores);
