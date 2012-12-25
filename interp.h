#ifndef INTERP_H
#define INTERP_H

#include <algorithm>

#include "parameters.h"

//input:
// x: velocity array
// y: xs array
// x_int: input velocity
//output:
// y_int: xs corresponding to x_int 
// only supports random access iterator
template <class T1_Iterator, class T2_Iterator>
typename std::iterator_traits<T2_Iterator>::value_type interp(T1_Iterator x, T2_Iterator y, typename std::iterator_traits<T1_Iterator>::value_type x_int, int N) {
	typename std::iterator_traits<T2_Iterator>::value_type y_int;
	if (x_int < x[0]) {
//		y_int = y[0]+(y[1]-y[0])/(x[1]-x[0])*(x_int-x[0]);
		// specific for 1/v xs
		if (x_int < 0.1*PARAM::vmin) {
			y_int = x[0]*y[0]/(0.1*PARAM::vmin);
		} else {
			y_int = x[0]*y[0]/x_int;
		}
	} else if (x_int == x[0]) {
		y_int = y[0];
	} else if (x_int > x[N-1]) {
		y_int = y[N-1]+(y[N-1]-y[N-2])/(x[N-1]-x[N-2])*(x_int-x[N-1]);
	} else if (x_int == x[N-1]) {
		y_int = y[N-1];  
	} else {
		const T2_Iterator p_int = std::lower_bound(x, x + N, x_int);
		unsigned int ind = p_int - x;
		if (x_int != x[ind]) {  // x_int falls in the interval
			y_int = y[ind]+(y[ind+1]-y[ind])/(x[ind+1]-x[ind])*(x_int-x[ind]);
		} else {          // x_int falls on a grid point
			y_int = y[ind];
		}
	} 
	return y_int;
} 

#endif // INTERP_H
