// correlation.cpp
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]

#include "priors.h"
#include "correlation.h"
#include "graycode.h"

#include <sys/time.h>
#include <unistd.h>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <boost/tokenizer.hpp>

using namespace boost;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::RowVectorXi;
using Eigen::MatrixXi;
using namespace std;
using namespace Rcpp;

const bool NUMERIC_FIX = false;
// #define DEBUG


vector<uint>& get_indices_from_dbitset(const dbitset& gamma, vector<uint>& v)
{
  for (size_t i = 0; i < gamma.size(); i++) {
    if (gamma[i]) {
      v.push_back(i);
    }
  }

  return v;
}


// Get columns
// Need a corresponding get rows
// And get rows and columns
MatrixXd& get_cols(const MatrixXd& m1, const dbitset& gamma, MatrixXd& m2)
{
  // Special case of get_rows_and_cols
  vector<uint> columns;
  columns = get_indices_from_dbitset(gamma, columns);

  for (size_t i = 0; i < columns.size(); i++) {
    m2.col(i) = m1.col(columns[i]);
  }

  return m2;
}


MatrixXd& get_rows(const MatrixXd& m1, const dbitset& gamma, MatrixXd& m2)
{
  // Special case of get_rows_and_cols
  vector<uint> rows;
  rows = get_indices_from_dbitset(gamma, rows);
  // MatrixXd m2(rows.size(), m1.cols());

  for (size_t i = 0; i < rows.size(); i++) {
    m2.row(i) = m1.row(rows[i]);
  }

  return m2;
}


MatrixXd& get_rows_and_cols(const MatrixXd& m1, const dbitset& rows_bs,
const dbitset& cols_bs, MatrixXd& m2)
{
  vector<uint> row_indices;
  row_indices = get_indices_from_dbitset(rows_bs, row_indices);
  vector<uint> col_indices;
  col_indices = get_indices_from_dbitset(cols_bs, col_indices);

  // Matrices are stored in column major order, so the intent is to access contiguous memory in
  // sequence by looping down each row in the inner loop.
  for (size_t j = 0; j < col_indices.size(); j++) {
    for (size_t i = 0; i < row_indices.size(); i++) {
      m2(i, j) = m1(row_indices[i], col_indices[j]);
    }
  }

  return m2;
}


double var(const VectorXd& v)
{
  const auto n = v.size();
  return (v.array() - v.mean()).array().square().sum() / (n - 1);
}


double sd(const VectorXd& v)
{
	return (sqrt(var(v)));
}


Normed normalise(VectorXd& vy, MatrixXd& mX)
{
  const auto n = vy.size();
  const auto p = mX.cols();
  VectorXd mu_mX(p);
  VectorXd sigma2_mX(p);
  #ifdef DEBUG
  Rcpp::Rcout << "vy " << vy.head(10) << std::endl;
  Rcpp::Rcout << "mX " << mX.topRows(10) << std::endl;
  #endif
  Normed normed;
  normed.vy = vy;
  normed.mX = mX;

  // Normalise vy and mX
  auto mu_vy = vy.mean();
  auto sigma2_vy = (n - 1) * var(vy) / n;
  normed.vy = (vy.array() - mu_vy) / sqrt(sigma2_vy);
  for (auto i = 0; i < p; i++) {
    mu_mX(i) = mX.col(i).mean();
    sigma2_mX(i) = (n - 1) * var(mX.col(i)) / n;
    normed.mX.col(i) = (mX.col(i).array() - mu_mX(i)) / sqrt(sigma2_mX(i));
  }
  #ifdef DEBUG
  Rcpp::Rcout << "vy " << vy.head(10) << std::endl;
  Rcpp::Rcout << "mX " << mX.topRows(10) << std::endl;
  Rcpp::Rcout << "mu_vy " << mu_vy << " sigma2_vy " << sigma2_vy << std::endl;
  Rcpp::Rcout << "mu_mX " << mu_mX << " sigma2_mX " << sigma2_mX << std::endl;
  #endif
  return normed;
}


void show_matrix_difference(ostream& os, const MatrixXd& m1, const MatrixXd& m2, const double epsilon = 1e-8)
{
  // Check that m1 and m2 have the same dimensions.
  if (!(m1.rows() == m2.rows() && m1.cols() == m2.cols())) {
    os << "Dimensions of m1 and m2 do not match, cannot display difference." << endl;
    return;
  }

  // Iterate through the elements of m1 and m2, looking for and reporting differences
  for (auto i = 0; i < m1.rows(); i++) {
    for (auto j = 0; j < m1.cols(); j++) {
      if (abs(m1(i, j) - m2(i, j)) > epsilon) {
        os << "Row " << i << ", column " << j << " m1 " << m1(i, j) << " m2 " << m2(i, j);
        os <<  " difference " << m1(i, j) - m2(i, j);
        os << " relative difference " << (m1(i, j) - m2(i, j)) / m1(i, j) << endl;
      }
    }
  }
}


// Perform the rank one update on (X_gamma^T X_gamma)^{-1}
MatrixXd& rank_one_update(const dbitset& gamma, const uint col_abs, const uint min_idx,
  const uint fixed,
const MatrixXd& mXTX, const MatrixXd& mA, MatrixXd& mA_prime, bool& bLow)
{
  auto p = mXTX.cols();
  auto p_gamma_prime = mA_prime.cols();
  auto p_gamma = mA.cols();

  if (gamma.count() == 0) {
    // No need to do a rank-one update. Just construct a 1 by 1 matrix.
    // TODO: This situation should never arise. But somehow it does. Find out
    // how.
    mA_prime.resize(1, 1);
    mA_prime << 1. / mXTX(col_abs, col_abs);
    bLow = false;
    return mA_prime;
  }
  // Construct mA_prime
  // b = 1 / (x^T x - x^T X_gamma A X_gamma^T x)
  auto xTx = mXTX(col_abs, col_abs);
  dbitset col_bs(p);
  col_bs[col_abs] = true;
  MatrixXd X_gamma_T_x(p_gamma, 1);
  #ifdef DEBUG
  Rcpp::Rcout << "gamma " << gamma << endl;
  Rcpp::Rcout << "col_bs " << col_bs << endl;
  #endif
  X_gamma_T_x = get_rows_and_cols(mXTX, gamma, col_bs, X_gamma_T_x);

  // const double b = 1 / (x^T x - x^T X_gamma A X_gamma^T x).value();
  const auto b = 1 / (xTx - (X_gamma_T_x.transpose() * mA * X_gamma_T_x).value());
  // b is supposed to be positive definite.
  #ifdef DEBUG
  Rcpp::Rcout << "b " << b << endl;
  #endif
  const auto epsilon = 1e-12;
  if (b > epsilon) {
    // Do rank one update
    // Matrix m1 = A + b A X_gamma^T x x^T X_gamma A
    // The relative column index
    auto col = col_abs - min_idx + fixed;
    MatrixXd X_gamma_x(p_gamma, p_gamma);
    MatrixXd A_X_gamma_T_x = mA * X_gamma_T_x;
    // Re-arrange.
    MatrixXd b_X_gamma_T_x_A = b * A_X_gamma_T_x;
    if (col == 0) {
      mA_prime << b, -b * A_X_gamma_T_x.transpose(),
        -b * A_X_gamma_T_x, mA + b * A_X_gamma_T_x * A_X_gamma_T_x.transpose();
    }
    else if (0 < col && col < p_gamma_prime - 1) {
      MatrixXd m1(p_gamma, p_gamma);
      m1 = mA + b * A_X_gamma_T_x * A_X_gamma_T_x.transpose();
      // mA_prime << 1, 2, 3,
      //            4, 5, 6,
      //            7, 8, 9;
      mA_prime.topLeftCorner(col, col) = m1.topLeftCorner(col, col);
      mA_prime.row(col).leftCols(col) = -b_X_gamma_T_x_A.topRows(col).transpose();
      mA_prime.col(col).bottomRows(p_gamma_prime - (col + 1)) = -b_X_gamma_T_x_A.bottomRows(p_gamma_prime - (col + 1));
      mA_prime(col, col) = b;
      mA_prime.bottomLeftCorner(p_gamma_prime - (col + 1), col) = m1.bottomLeftCorner(p_gamma_prime - (col + 1), col);
      mA_prime.bottomRightCorner(p_gamma_prime - (col + 1), p_gamma_prime - (col + 1)) = m1.bottomRightCorner(p_gamma_prime - (col + 1), p_gamma_prime - (col + 1));

      // Should take advantage of the symmetry of mA_prime. For now, just fill in upper triangular entries.
      #ifdef DEBUG
      Rcpp::Rcout << "mA_prime " << endl << mA_prime << endl;
      #endif
      for (auto j = 0; j < p_gamma_prime; j++) {
        for (auto i = 0; i < j; i++) {
          mA_prime(i, j) = mA_prime(j, i);
        }
      }
      #ifdef DEBUG
      Rcpp::Rcout << "mA_prime " << mA_prime << endl;
      #endif
    } else                                   // col == p_gamma_prime
    {
      mA_prime << mA + b * A_X_gamma_T_x * A_X_gamma_T_x.transpose(), -b_X_gamma_T_x_A,
        -b_X_gamma_T_x_A.transpose(), b;
    }
    bLow = false;
  }
  else {
    // Signal that a rank one update was impossible so that the calling code can perform a full inversion.
    bLow = true;
    // bLow = false;
  }

  return mA_prime;
}


// Perform the rank one downdate on (X_gamma^T X_gamma)^{-1}
MatrixXd& rank_one_downdate(const uint col_abs, const uint min_idx, const uint fixed,
const MatrixXd& mA, MatrixXd& mA_prime)
{
  auto p_gamma_prime = mA_prime.cols();
  auto p_gamma = mA.cols();
  // The relative column index
  auto col = col_abs - min_idx + fixed;

  // Need to deal with three cases
  if (col == 0) {
    const MatrixXd mA_11 = mA.bottomRightCorner(p_gamma_prime, p_gamma_prime);
    const VectorXd va_12 = mA.col(0).tail(p_gamma_prime);
    // const MatrixXd va_12 = mA.block(0, 0, p_gamma_prime, 1);
    auto a_22 = mA(0, 0);
    mA_prime = mA_11 - (va_12 * va_12.transpose()) / a_22;
  }
  else if (1 <= col && col <= p_gamma - 1) {
    // 1 2 3
    // 4 5 6
    // 7 8 9
    MatrixXd mA_11(p_gamma_prime, p_gamma_prime);
    VectorXd va_12(p_gamma_prime);
    mA_11.topLeftCorner(col, col) = mA.topLeftCorner(col, col);
    mA_11.bottomRows(p_gamma_prime - col).leftCols(col) = mA.bottomRows(p_gamma_prime - col).leftCols(col);
    mA_11.bottomRightCorner(p_gamma_prime - col, p_gamma_prime - col) = mA.bottomRightCorner(p_gamma_prime - col, p_gamma_prime - col);
    va_12.head(col) = mA.col(col).head(col);
    va_12.tail(p_gamma_prime - col) = mA.col(col).tail(p_gamma_prime - col);
    auto a_22 = mA(col, col);
    mA_prime = mA_11 - (va_12 * va_12.transpose()) / a_22;
  } else                                     // col == p_gamma
  {
    const MatrixXd mA_11 = mA.topLeftCorner(p_gamma_prime, p_gamma_prime);
    const VectorXd va_12 = mA.col(p_gamma - 1).head(p_gamma - 1);
    auto a_22 = mA(p_gamma - 1, p_gamma - 1);
    mA_prime = mA_11 - (va_12 * va_12.transpose()) / a_22;
  }

  // Should take advantage of the symmetry of mA_prime. For now, just fill in upper triangular entries.
  for (auto j = 0; j < p_gamma_prime; j++) {
    for (auto i = 0; i < j; i++) {
      mA_prime(i, j) = mA_prime(j, i);
    }
  }

  return mA_prime;
}


void update_mA_prime(bool bUpdate, const dbitset& gamma,
const uint col, const uint min_idx, const uint fixed,
const MatrixXd& mXTX, const MatrixXd& mA, MatrixXd& mA_prime,
bool& bLow)
{
  if (bUpdate) {
    // Rank one update of mA_prime from mA
    #ifdef DEBUG
    Rcpp::Rcout << "Updating " << col << endl;
    #endif
    mA_prime = rank_one_update(gamma, col, min_idx, fixed, mXTX, mA, mA_prime, bLow);
  }
  else {
    // Rank one downdate
    #ifdef DEBUG
    Rcpp::Rcout << "Downdating " << col << endl;
    #endif
    mA_prime = rank_one_downdate(col, min_idx, fixed, mA, mA_prime);
  }
}


double logp2(int n, double R2, int p)
{
  auto a = -3./4.;
  auto b = (n - 5.) / 2. - p / 2. - a;
  auto result = lgamma(p / 2. + a + 1.) + lgamma(a + b + 2.) - lgamma((n - 1.) / 2) - lgamma(a + 1.);
  result += -(b + 1.) * log(1. - R2);

  return result;
}


bool check_model_prior_parameters(const std::string modelprior, const VectorXd& modelpriorvec,
                                  const MatrixXd& mX, std::string& invalid_reason)
{
  if (modelprior == "uniform") {
    return true;
  } else if (modelprior == "beta-binomial") {
    if (modelpriorvec.size() != 2) {
      invalid_reason = "modelpriorvec was not of size 2";
      return false;
    }

    if (modelpriorvec(0) <= 0.|| modelpriorvec(1) <= 0.) {
      invalid_reason = "elements of modelpriorvec should be positive";
      return false;
    }

    return true;
  } else if (modelprior == "bernoulli") {
    if (modelpriorvec.size() != mX.cols()) {
      invalid_reason = "modelpriorvec was not of the same size as the number of columns in mX";
      return false;
    }

    if (abs(modelpriorvec.sum() - 1.) > 1e-5) {
      invalid_reason = "modelpriorvec does not sum to 1";
      return false;
    }

    return true;
  } else {
    invalid_reason = "modelprior unknown";
    return false;
  }
}


void calculate_probabilities(const std::string prior, const std::string modelprior, const VectorXd& modelpriorvec,
                             const int n, const int p, const VectorXd& vR2_all,
                             const VectorXi& vpgamma_all, const Graycode& graycode,
                             VectorXd& vlogp_all,
                             VectorXd& vinclusion_prob)
{
  log_prob_fn log_prob;
  set_log_prob(prior, log_prob);

  auto nmodels = vR2_all.size();
  #pragma omp parallel for
  for (auto i = 0; i < nmodels; i++) {
    vlogp_all(i) = log_prob(n, vpgamma_all(i), vR2_all(i));
    if (modelprior == "beta-binomial") {
      double alpha = modelpriorvec(0);
      double beta = modelpriorvec(1);
      vlogp_all(i) += ::Rf_lbeta(alpha + vpgamma_all(i), beta + p - vpgamma_all(i)) - ::Rf_lbeta(alpha, beta);
    }
    if (modelprior == "bernoulli") {
      for (auto j = 0; j < p; j++) {
        auto gamma = graycode[i][j] ? 1. : 0.;
        if (modelpriorvec(j) == 0. || modelpriorvec(j) == 1.)
          continue;
        vlogp_all(i) += gamma * log(modelpriorvec(j)) + (1 - gamma) * log(1. - modelpriorvec(j));
      }
    }
  }
  auto M = vlogp_all.array().maxCoeff(); // Maximum log-likelihood
  // Rcpp::Rcout << "M " << M << std::endl;
  VectorXd vmodel_prob = (vlogp_all.array() - M).array().exp() / (vlogp_all.array() - M).array().exp().sum();
  // If we have memory problems, this can be recoded so that we don't have to construct the entire mGamma
  // matrix
  // MatrixXd mGamma = graycode.to_MatrixXi().cast<double>();
  // vinclusion_prob = mGamma.transpose() * vmodel_prob;
  vinclusion_prob = VectorXd::Zero(p);
  for (int i = 0; i < nmodels; i++) {
    for (int j = 0; j < p; j++) {
      auto gamma = graycode[i][j] ? 1. : 0.;
      vinclusion_prob(j) += gamma * vmodel_prob(i);
    }
  }
}


// Calculate the correlations for every subset of the covariates in mX
List all_correlations_main(const Graycode& graycode, VectorXd vy, MatrixXd mX, std::string prior,
  std::string modelprior, VectorXd modelpriorvec,
  const uint fixed, const uint intercept_col, const uint max_iterations, const bool bNatural_Order = false,
  const bool bIntercept = false,
  const bool bCentre = true, uint cores = 1L)
{
  #ifdef _OPENMP
    // Eigen::initParallel();
    // omp_set_num_threads(cores);
    // Eigen::setNbThreads(cores);
  #endif

  #ifdef DEBUG
    // None of R's functions, such as Rcpp::checkUserInterrupt() or
    // Rcpp::Rcout, are threadsafe. If you're debugging and calling these
    // functions in multiple threads, R will crash.
    omp_set_num_threads(1);
  #endif

  const uint n = mX.rows();                  // The number of observations
  const uint p = mX.cols();                  // The number of covariates
  VectorXd vR2_all(max_iterations);          // Vector of correlations for all models
  VectorXi vpgamma_all(max_iterations);      // Vector of number of covariates included in each model
  bool bmA_set = false;                      // Whether mA has been set yet
  bool bUpdate;                              // True for an update, false for a downdate
  uint diff_idx;                             // The covariate which is changing
  uint min_idx;                              // The minimum bit which is set in gamma_prime
  dbitset gamma(p);                          // The model gamma
  dbitset gamma_prime(p);                    // The model gamma_prime
  uint p_gamma_prime;                        // The number of columns in the matrix mX_gamma_prime
  uint p_gamma;                              // The number of columns in the matrix mX_gamma
  vector<MatrixXd> vec_mA(p);
  vector<MatrixXd> vec_mX_gamma(p);
  vector<MatrixXd> vec_m1(p);
  const MatrixXd mXTX = mX.transpose() * mX;
  const MatrixXd mXTy = mX.transpose() * vy;
  const double yTy = vy.squaredNorm();

  // Check that modelprior parameters are correct
  std::string invalid_reason;
  if (!check_model_prior_parameters(modelprior, modelpriorvec, mX, invalid_reason)) {
    stringstream ss;
    ss << "modelprior parameters are invalid - " << invalid_reason;
    Rcpp::stop(ss.str());
  }

  // Pre-allocate memory
  for (uint i = 0; i < p; i++) {
    vec_mA[i].resize(i + 1, i + 1);
    vec_mX_gamma[i].resize(n, i + 1);
    vec_m1[i].resize(i + 1, 1);
  }

  // if (bCentre) {
  	Normed normed = normalise(vy, mX);
    vy = normed.vy;
    mX = normed.mX;
  // }

  vpgamma_all(0) = 0;
  vR2_all(0) = 0.;

  // Loop through models, updating and downdating mA as necessary
  Rcpp::checkUserInterrupt();
  #pragma omp parallel for\
    firstprivate(gamma, gamma_prime, bmA_set, vec_mX_gamma, vec_mA, vec_m1)\
    private(diff_idx, min_idx, p_gamma_prime, p_gamma, bUpdate)\
      shared(mX, vR2_all, vpgamma_all, graycode)\
      default(none)
  for (uint idx = 1; idx < max_iterations; idx++) {
    #ifdef DEBUG
    Rcpp::Rcout << endl << "Iteration " << idx << endl;
    #endif
    // By properties of Greycode, only one element can be different. And it's either one higher or
    // one lower.
    // Check if update or downdate, and for which variable
    gamma = gamma_prime;
    gamma_prime = graycode[idx];

    #ifdef DEBUG
    Rcpp::Rcout << "Previous gamma: " << gamma << endl;
    Rcpp::Rcout << "Current gamma:  " << gamma_prime << endl;
    #endif

    graycode.change(gamma_prime, gamma, bUpdate, diff_idx, min_idx, p_gamma_prime);

    #ifdef DEBUG
    Rcpp::Rcout << "min_idx " << min_idx << " diff_idx " << diff_idx;
    Rcpp::Rcout << " bUpdate " << (bUpdate ? "true" : "false") << endl;
    #endif

    // Get mX matrix for gamma
    MatrixXd& mA = vec_mA[p_gamma - 1];
    MatrixXd& mA_prime = vec_mA[p_gamma_prime - 1];
    MatrixXd& mX_gamma_prime = vec_mX_gamma[p_gamma_prime - 1];

    // Only needed when bmA_set is false.
    // If we haven't previously calculated this inverse, calculate it the first time.
    if (!bmA_set) {
      // Calculate full inverse mA, O(p^3)
      mX_gamma_prime = get_cols(mX, gamma_prime, mX_gamma_prime);
      mA_prime = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
      bmA_set = true;
    }
    else {
      bool bLow;
      #ifdef DEBUG
      Rcpp::Rcout << "mA_prime before update " << mA_prime << endl;
      #endif
      update_mA_prime(bUpdate, gamma,
        diff_idx, min_idx, fixed,
        mXTX,   mA, mA_prime,
        bLow);
      if (bLow) {
        mX_gamma_prime = get_cols(mX, gamma_prime, mX_gamma_prime);
        mA_prime = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
      }
      #ifdef DEBUG
      Rcpp::Rcout << "mA_prime after update " << mA_prime << endl;
      #endif

      #ifdef DEBUG
      // Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
      mX_gamma_prime = get_cols(mX, gamma_prime, mX_gamma_prime);
      MatrixXd identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime;
      if (!identity_prime.isApprox(MatrixXd::Identity(p_gamma_prime, p_gamma_prime)) && NUMERIC_FIX) {
        Rcpp::Rcout << "(mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime" << endl;
        Rcpp::Rcout << identity_prime << endl;
        Rcpp::Rcout << "Iterative calculation of inverse is wrong, recalculating ..." << endl;
        // This inverse is nonsense. Do a full inversion.
        MatrixXd mA_prime_full = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
        show_matrix_difference(Rcpp::Rcout, mA_prime, mA_prime_full);
        // TODO: Find the differences between mA_prime_full and mA_prime
        identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime_full;
        mA_prime = mA_prime_full;
      }

      // Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
      Rcpp::Rcout << "(mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime" << endl;
      Rcpp::Rcout << identity_prime << endl;
      #endif
    }

    double R2;
    double numerator;
    // Can pre-compute, using vXTy
    // VectorXd v1(p_gamma_prime);
    // v1 = vy.transpose() * mX_gamma_prime;
    MatrixXd& m1 = vec_m1[p_gamma_prime - 1];
    m1 = get_rows(mXTy, gamma_prime, m1);
    numerator = (m1.transpose() * mA_prime * m1).value();
    R2 = numerator / yTy;
    #ifdef DEBUG
    Rcpp::Rcout << "m1 " << m1 << endl;
    Rcpp::Rcout << "mA_prime " << mA_prime << endl;
    Rcpp::Rcout << "Numerator " << numerator << " denominator " << yTy;
    Rcpp::Rcout << " R2 " << R2 << endl;
    #endif
    vR2_all(idx) = R2;                       // Calculate correlation
    vpgamma_all(idx) = p_gamma_prime;

    p_gamma = p_gamma_prime;
    // FIXME: How do you do this in a thread-safe way?
    // Rcpp::checkUserInterrupt();
  }

  double R2_full = (mXTy.transpose() * mXTX.inverse() * mXTy).value() / yTy; // R2 of full model
  auto p_star = 0; // Count of significant coefficients
  VectorXd vbeta_hat = (mXTX.inverse()) * mXTy;
  VectorXd se_beta_hat = (1 - R2_full) * mXTX.inverse().diagonal();
  VectorXd t_beta = vbeta_hat.array() / se_beta_hat.array();
  double threshold = ::Rf_qt(0.975, n - p, 1, 0);
  for (uint i = 1; i < p; i++) {
    if (abs(t_beta(i)) > threshold)
      p_star++;
  }

  // auto M = 0.;
  VectorXd vlogp_all(max_iterations);        // Vector of model likelihoods
  VectorXd vinclusion_prob(p);               // Vector of variable inclusion likelihoods
  calculate_probabilities(prior, modelprior, modelpriorvec,
                          n, p, vR2_all, vpgamma_all, graycode,
                          vlogp_all, vinclusion_prob);

  // Rcpp::Rcout << "vmodel_prob " << vmodel_prob << std::endl;
  // Rcpp::Rcout << "vinclusion_prob " << vinclusion_prob << std::endl;

  if (!bNatural_Order) {
    return List::create(Named("vR2") = vR2_all,
                        Named("vp_gamma") = vpgamma_all,
                        Named("vlogp") = vlogp_all,
                        Named("vinclusion_prob") = vinclusion_prob);
  } else {
    VectorXd vR2(max_iterations);
    VectorXi vp_gamma(max_iterations);
    VectorXd vlogp(max_iterations);
    for (uint i = 1; i < max_iterations; i++) {
      vR2(i) = vR2_all(graycode.gray_to_binary(i));
      vp_gamma(i) = vpgamma_all(graycode.gray_to_binary(i));
      vlogp(i) = vlogp_all(graycode.gray_to_binary(i));
    }

    return List::create(Named("vR2") = vR2,
                        Named("vp_gamma") = vp_gamma,
                        Named("vlogp") = vlogp,
                        Named("vinclusion_prob") = vinclusion_prob);
  }
}

// [[Rcpp:export]]
List blma_cpp(VectorXd vy, MatrixXd mX, std::string prior, std::string modelprior,
              VectorXd modelpriorvec,
              const uint intercept_col,
              const bool bNatural_Order, const bool bIntercept, const bool bCentre,
              const uint cores)
{
  const uint p = mX.cols();
  const uint fixed = 0;
  const uint max_iterations = 1 << p;

  Graycode graycode(p);
  return all_correlations_main(graycode, vy, mX, prior, modelprior, modelpriorvec,
                                fixed, intercept_col, max_iterations, bNatural_Order,
                                bIntercept, bCentre, cores);
}


// Calculate the correlations for every subset of the covariates in mX
// [[Rcpp:export]]
List blma_fixed_cpp(VectorXd vy, MatrixXd mX, MatrixXd mZ, std::string prior,
                    std::string modelprior, VectorXd modelpriorvec,
                    const uint intercept_col, const bool bNatural_Order, const bool bIntercept,
                    const bool bCentre, const uint cores)
{
  const uint n = mX.rows();
  const uint p1 = mX.cols();
  const uint p2 = mZ.cols();
  MatrixXd mC(n, p1 + p2);
  const uint max_iterations = 1 << p2;

  // mC << mX, mZ;
  mC.leftCols(p1) = mX;
  mC.rightCols(p2) = mZ;
  Graycode graycode(p1, p2);
  return all_correlations_main(graycode, vy, mC, prior, modelprior, modelpriorvec,
                                p1, intercept_col, max_iterations, bNatural_Order,
                                bIntercept, bCentre), cores;
}


VectorXd one_correlation(VectorXd vy, MatrixXd mX, MatrixXd mZ)
{
  const uint n = mX.rows();
  const uint p = mX.cols();
  const uint m = mZ.cols();

  MatrixXd mC(n, p + m);
  mC << mX, mZ;

  MatrixXd m1(p + m, p + m);                 // Symmetric
  VectorXd m2(p + m);
  m1 << mX.transpose() * mX, mX.transpose() * mZ,
    mZ.transpose() * mX, mZ.transpose() * mZ;
  m2 << mX.transpose() * vy,
    mZ.transpose() * vy;
  VectorXd R2 = (vy.transpose() * mC * m1.inverse() * m2) / vy.squaredNorm();

  return R2;
}
