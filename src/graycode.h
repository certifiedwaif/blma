// graycode.h

// [[Rcpp::depends(BH)]]

#pragma once

#include <Rcpp.h>
#include <boost/dynamic_bitset.hpp>
#include <Eigen/Dense>

using namespace boost;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::RowVectorXi;
using Eigen::MatrixXi;
using namespace std;

typedef dynamic_bitset<> dbitset;

struct Graycode {
  Graycode(int _p);
  Graycode(int _fixed, int _varying);
  const int fixed;
  const int varying;
  const int size;

  int binary_to_gray(const int num) const;
  int gray_to_binary(const int num) const;
  VectorXd binary_to_vec(const int num);
  VectorXd gray_vec(const int i);
  MatrixXi to_MatrixXi() const;
  dbitset operator[](const int idx) const;
  void change(const dbitset& gamma_prime, const dbitset& gamma,
              bool& update, int& diff_idx, int& min_idx,
              int& bits_set) const;
};
