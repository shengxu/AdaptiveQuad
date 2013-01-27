#ifndef UTIL_H
#define UTIL_H

#include <cmath>

#include "parameters.h"
#include "interp.h"
#include "productcdf.h"

//double interp(double *x, double *y, double x_int, int N);
//template <typename T1, typename T2>
//T2 interp(const T1 *x, const T2 *y, T1 x_int, int N);


double adaprandomint(double E, const CDFmuvt & cdf, const std::vector<double> &xs_E, const std::vector<double> &xs_sig, double eps);


// E in eV, v in m/s
inline double EtoV(const double E) {
	return std::sqrt(2*E*1.e-6/CONST::M_NEUT);
}

inline double VtoE(const double v) {
	return 0.5*CONST::M_NEUT*1.e6*v*v;
}

#endif  // UTIL_H
