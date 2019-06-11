// graycode.cpp

#include "graycode.h"

Graycode::Graycode(int _p) : fixed(0), varying(_p), size(fixed + varying)
{ 
}

Graycode::Graycode(int _fixed, int _varying) : fixed(_fixed), varying(_varying), size(fixed + varying)
{
}


int Graycode::binary_to_gray(int num) const
{
  	return (num >> 1) ^ num;
}


/*
   The purpose of this function is to convert a reflected binary
   Gray code number to a binary number.
   */
int Graycode::gray_to_binary(int num) const
{
  	int mask;
  	for (mask = num >> 1; mask != 0; mask = mask >> 1) {
    	num = num ^ mask;
  	}
  	return num;
}


VectorXd Graycode::binary_to_vec(int num)
{
  	VectorXd result(size);
  	for (int i = 0; i < size; i++) {
    	result[(size - 1) - i] = num & 1;
    	num >>= 1;
  	}
  	return(result);
}


VectorXd Graycode::gray_vec(int i)
{
  	return binary_to_vec(binary_to_gray(i)).transpose();
}


MatrixXi Graycode::to_MatrixXi() const
{
  	int rows = 1 << varying;
  	MatrixXi result(rows, size);
#pragma omp parallel for
  	for (int i = 0; i < rows; i++) {
    	dbitset bs = (*this)[i];
    	for (int j = 0; j < size; j++) {
      		result(i, j) = bs[j] ? 1 : 0;
    	}
  	}
  	return(result);
}


dbitset Graycode::operator[](const int idx) const
{
  	dbitset bs_varying(varying, idx);
  	bs_varying =  bs_varying ^ (bs_varying >> 1);
  	if (fixed != 0) {
    	dbitset bs(size);

    	for (int i = 0; i < fixed; i++) {
      		bs[i] = true;
    	}

    	for (int i = 0; i < varying; i++) {
      		bs[i + fixed] = bs_varying[i];
    	}

    	return bs;
  	} else {
    	return bs_varying;
  	}
}


//' Given a new dbitset, find out whether we are updating or downdating, which bit has changed in the
//' new dbitset, what the minimum bit is which is set and which bits are set.
//'
//' @param[in] gamma_prime dbitset The new bitset we are changing to
//' @param[in] gamma dbiset The old bitset that we are changing from
//' @param[out] update bool A flag which is true if we are updating, and false if we are downdating
//' @param[out] diff_idx int Which bit has changed from gamma to gamma_prime
//' @param[out] min_idx The minimum index of bit which is set
//' @param[out] bits_set A count of how many bits are set in gamma_prime
void Graycode::change(const dbitset& gamma_prime, const dbitset& gamma,
        bool& update, int& diff_idx, int& min_idx, int& bits_set) const
{

  	// Find the LSB of the varying bitset.
  	// min_idx = min(gamma.find_first(), gamma_prime.find_first());
  	// min_idx = min(gamma.find_next(fixed), gamma_prime.find_next(fixed));
  	for (auto idx = fixed; idx < size; idx++) {
    	if (gamma[idx] || gamma_prime[idx]) {
      		min_idx = idx;
      		break;
    	}
  	}

  	// Find bit that has changed.
  	// #ifdef DEBUG
  	// Rcpp::Rcout << "gamma_prime ^ prime " << (gamma_prime ^ gamma) << endl;
  	// #endif
  	// diff_idx = (gamma_prime ^ gamma).find_next(fixed - 1);
  	for (auto idx = fixed; idx < size; idx++) {
    	if (gamma[idx] != gamma_prime[idx]) {
      		diff_idx = idx;
      		break;
    	}
  	}

  	// Has it been set, or unset?
  	update = gamma_prime[diff_idx];

  	bits_set = gamma_prime.count();

}
