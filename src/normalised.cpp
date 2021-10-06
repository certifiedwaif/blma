// normalised.cpp
#include "normalised.h"

double var(const VectorXd& v)
{
  	const auto n = v.size();
  	VectorXd normalised = v.array() - v.mean();
  	return normalised.squaredNorm() / (n - 1);
}


double sd(const VectorXd& v)
{
	return (sqrt(var(v)));
}


Normed normalise(VectorXd& vy, MatrixXd& mX)
{
  	const auto n = vy.size();
  	const auto p = mX.cols();
#ifdef DEBUG
  	Rcpp::Rcout << "n " << n << std::endl;
  	Rcpp::Rcout << "p " << p << std::epdl;
#endif
  	VectorXd mu_mX(p);
  	VectorXd sigma2_mX(p);
#ifdef DEBUG
  	Rcpp::Rcout << "vy " << vy.head(10) << std::endl;
  	Rcpp::Rcout << "mX " << mX.topRows(10) << std::endl;
#endif
  	Normed normed(vy, mX);

  	// Normalise vy and mX
  	auto mu_vy = normed.vy.mean();
  	auto sigma2_vy = (n - 1) * var(normed.vy) / n;
  	normed.vy = (normed.vy.array() - mu_vy) / sqrt(sigma2_vy);
#ifdef DEBUG
  	Rcpp::Rcout << normed.vy.mean() << std::endl;
  	Rcpp::Rcout << normed.vy.squaredNorm() << std::endl;
#endif
  	for (auto i = 0; i < p; i++) {
    	mu_mX(i) = normed.mX.col(i).mean();
    	sigma2_mX(i) = (n - 1) * var(normed.mX.col(i)) / n;
    	normed.mX.col(i) = (normed.mX.col(i).array() - mu_mX(i)) / sqrt(sigma2_mX(i));
#ifdef DEBUG
		Rcpp::Rcout << normed.mX.col(i).mean() << std::endl;
		Rcpp::Rcout << normed.mX.col(i).squaredNorm() << std::endl;
#endif
  	}
#ifdef DEBUG
  	Rcpp::Rcout << "vy " << vy.head(10) << std::endl;
  	Rcpp::Rcout << "mX " << mX.topRows(10) << std::endl;
  	Rcpp::Rcout << "mu_vy " << mu_vy << " sigma2_vy " << sigma2_vy << std::endl;
  	Rcpp::Rcout << "mu_mX " << mu_mX << " sigma2_mX " << sigma2_mX << std::endl;
#endif
  	return normed;
}

