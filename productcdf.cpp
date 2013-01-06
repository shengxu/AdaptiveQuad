#include <cmath>

#include "parameters.h"
#include "util.h"
#include "productcdf.h"

using namespace std;

static inline double getcdf(double x, double sqalpha) {
	double cdf;
	if (x <= 0) {
		cdf = 0.5 - 0.5*std::erf(-sqalpha*x);
	} else {
		cdf = 0.5 + 0.5*std::erf(sqalpha*x);
	}
	return cdf;
}

static inline double getcdf_MB(double x, double alpha, double sqalpha) {
	return std::erf(sqalpha*x) - 2/CONST::SQRT_PI*sqalpha*x*std::exp(-alpha*std::pow(x, 2.0));
}

void CDFmuvt::setcdf() {
//	double pdf_p, pdf;
	double sqalpha = std::sqrt(alpha);
	
//	cdf.push_back(0);
//	pdf_p = getpdf(grid[0], sqalpha);
//	for (unsigned int i = 0; i < grid.size(); i++) {
//		pdf = getpdf(grid[i], sqalpha);
//		cdf.push_back(0.5*(pdf_p + pdf)*(grid[i] - grid[i-1]));
//		pdf_p = pdf;
	for (auto it = grid.begin(); it != grid.end(); it++) {
		cdf.push_back(getcdf(*it, sqalpha));
	}
}

double CDFmuvt::getx(double cdfin) const {
//	const std::vector<double>::iterator it = std::lower_bound(cdf.begin(), cdf.end(), cdfin);
//	unsigned int ind = it - cdf.begin();
//	return grid[ind] + (cdfin - cdf[ind])*(cdf[ind+1] - cdf[ind])/(grid[ind+1] - grid[ind]);
	return interp(cdf.begin(), grid.begin(), cdfin, cdf.size());
}


void CDFMB::setcdf() {
	double sqalpha = std::sqrt(alpha);
	for (auto it = grid.begin(); it != grid.end(); it++) {
		cdf.push_back(getcdf_MB(*it, alpha, sqalpha));
	}
}
