#pragma once

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

// #define DEBUG

using namespace std;
using namespace Rcpp;


VectorXd gamma_to_row(const dbitset& gamma);
List sampler(const int iterations,
            const NumericVector vy_in,
            const NumericMatrix mX_in,
            const std::string prior,
            const std::string modelprior,
            const Nullable<NumericVector> modelpriorvec_in,
            const int cores);
