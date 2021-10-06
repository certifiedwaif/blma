// normalised.h
#pragma once

#include <Rcpp.h>
#include <RcppEigen.h>
#define EIGEN_USE_BLAS

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]


using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::MatrixXd;

struct Normed {
    VectorXd vy;
    MatrixXd mX;

    Normed(VectorXd vy, MatrixXd mX): vy(vy), mX(mX) {}

    // Copy constructor
    Normed(const Normed& other) : vy(other.vy), mX(other.mX) {}

    // Move constructor
    Normed(const Normed&& other) : vy(other.vy), mX(other.mX) {}
};

Normed normalise(VectorXd& vy, MatrixXd& mX);
