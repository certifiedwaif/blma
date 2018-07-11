// cva.h
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

typedef log_prob_fn std::function<double (const int N, const int p, double vR2, int vp_gamma)>;

void set_log_prob(const string prior, std::function<double (const int n, const int p, double vR2, int vp_gamma)>& log_prob);
double calculate_log_prob(const uint n, const uint p, const double R2, const uint p_gamma,
                          const dbitset& gamma,
                          const std::function<double (const int n, const int p, double vR2, int vp_gamma)> log_prob,
                          const std::string modelprior, const VectorXd& modelpriorvec);
template <typename Derived1, typename Derived2>
void calculate_log_probabilities(const vector< dbitset >& gamma, const Eigen::MatrixBase<Derived1>& sigma2,
                                  const int n,
                                  Eigen::MatrixBase<Derived2>& log_probs,
                                  const std::function<double (const int n, const int p, double vR2, int vp_gamma)> log_prob,
                                  const std::string& modelprior, const VectorXd& modelpriorvec);
template <typename Derived1, typename Derived2, typename Derived3>
void calculate_weights(const Eigen::MatrixBase<Derived1>& sigma2,
                        const Eigen::MatrixBase<Derived2>& log_probs,
                        Eigen::MatrixBase<Derived3>& w);
template <typename Derived1>
double calculate_entropy(const Eigen::MatrixBase<Derived1>& w);
template <typename Derived1, typename Derived2>
double calculate_w_dot_prob(const Eigen::MatrixBase<Derived1>& w,
                            const Eigen::MatrixBase<Derived2>& log_probs);
template <typename Derived1>
void gamma_to_MatrixXd(const vector< dbitset >& gamma, Eigen::MatrixBase<Derived1>& m);
template <typename Derived1, typename Derived2, typename Derived3>
void calculate_mXTX_inv_prime(const dbitset& gamma, const dbitset& gamma_prime, int j,
                              const Eigen::MatrixBase<Derived1>& mXTX, const Eigen::MatrixBase<Derived2>& mXTX_inv,
                              Eigen::MatrixBase<Derived3>& mXTX_inv_prime,  bool bUpdate);
template <typename Derived1, typename Derived2, typename Derived3>
double calculate_sigma2_prime(const uint n, const uint p_gamma_prime,
                              const Eigen::MatrixBase<Derived1>& mX, const dbitset& gamma_prime,
                              const Eigen::MatrixBase<Derived2>& vy,
                              const Eigen::MatrixBase<Derived3>& mXTX_inv_prime);
List cva(const NumericVector vy_in, const NumericMatrix mX_in,
         const NumericMatrix mGamma_in,
         const std::string prior,
         const std::string modelprior, const Nullable<NumericVector> modelpriorvec_in = R_NilValue,
         const bool bUnique = true,
         const double lambda = 1.,
         const int cores = 1L);
