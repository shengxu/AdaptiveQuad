#ifndef INTERP_H
#define INTERP_H

#include <algorithm>

template <typename T1, typename T2>
T2 interp(const T1 *x, const T2 *y, T1 x_int, int N) {
	T2 y_int;
	if (x_int < x[0]) {
		return y[0]+(y[1]-y[0])/(x[1]-x[0])*(x_int-x[0]);
	} else if (x_int == x[0]) {
		return y[0];
	} else if (x_int > x[N-1]) {
		return y[N-1]+(y[N-1]-y[N-2])/(x[N-1]-x[N-2])*(x_int-x[N-1]);
	} else if (x_int == x[N-1]) {
		return y[N-1];  
	} else {
		const T2 *p_int = std::lower_bound(x, x + N, x_int);
		unsigned int ind = p_int - x;
		if (x_int != x[ind]) {  // x_int falls in the interval
			y_int = y[ind]+(y[ind+1]-y[ind])/(x[ind+1]-x[ind])*(x_int-x[ind]);
		} else {          // x_int falls on a grid point
			y_int = y[ind];
		}
		return y_int;
	} 
} 

#endif // INTERP_H
