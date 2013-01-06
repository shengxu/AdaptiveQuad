//double interp(double *x, double *y, double x_int, int N);
//template <typename T1, typename T2>
//T2 interp(const T1 *x, const T2 *y, T1 x_int, int N);
#include "interp.h"
#include "productcdf.h"

double adaprandomint(double E, const CDFmuvt & cdf, const std::vector<double> &xs_E, const std::vector<double> &xs_sig, double eps);
