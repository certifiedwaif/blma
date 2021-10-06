// correlation.hpp

#pragma once

#include <string>
#include <Rcpp.h>
#include <RcppEigen.h>
#define EIGEN_USE_BLAS
#include "normalised.h"
#include "graycode.h"
#include "priors.h"
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

// "If you want work with matrices, you should use C++ with Eigen or Armadillo. It's pretty fast." - Hadley Wickham,
// completely unprompted.

using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::string;

MatrixXd parseCSVfile_double(string infilename);
List blma_fixed_cpp(VectorXd vy, MatrixXd mX, MatrixXd mZ, std::string prior, std::string modelprior, VectorXd modelpriorvec,
                    const int intercept_col, const bool bNatural_Order = false, const bool bIntercept = false,
                    const bool bCentre = true, const int cores = 1L);
List blma_cpp(VectorXd vy, MatrixXd mX, std::string prior, std::string modelprior, VectorXd modelpriorvec,
              const int intercept_col, const bool bNatural_Order = false, const bool bIntercept = false,
              const bool bCentre = true, const int cores = 1L);
MatrixXd& get_cols(const MatrixXd& m1, const dbitset& gamma, MatrixXd& m2);
MatrixXd& get_rows(const MatrixXd& m1, const dbitset& gamma, MatrixXd& m2);
MatrixXd& get_rows_and_cols(const MatrixXd& m1, const dbitset& rows_bs,
const dbitset& cols_bs, MatrixXd& m2);
MatrixXd& rank_one_update(const dbitset& gamma, const int col_abs, const int min_idx,
  const int fixed,
const MatrixXd& mXTX, const MatrixXd& mA, MatrixXd& mA_prime, bool& bLow);
MatrixXd& rank_one_downdate(const int col_abs, const int min_idx, const int fixed,
const MatrixXd& mA, MatrixXd& mA_prime);
double var(const VectorXd& v);
double sd(const VectorXd& v);
Normed normalise(VectorXd& vy, MatrixXd& mX);
