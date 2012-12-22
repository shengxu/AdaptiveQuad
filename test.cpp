#include <iostream>
#include <cmath>

#include "util.h"
#include "xsdata.h"

using namespace std;

#include <algorithm>

//template <typename T1, typename T2>
//T2 interp(const T1 *x, const T2 *y, T1 x_int, int N) {
//	T2 y_int;
//	if (x_int < x[0]) {
//		return y[0]+(y[1]-y[0])/(x[1]-x[0])*(x_int-x[0]);
//	} else if (x_int == x[0]) {
//		return y[0];
//	} else if (x_int > x[N-1]) {
//		return y[N-1]+(y[N-1]-y[N-2])/(x[N-1]-x[N-2])*(x_int-x[N-1]);
//	} else if (x_int == x[N-1]) {
//		return y[N-1];  
//	} else {
//		const T2 *p_int = std::lower_bound(x, x + N, x_int);
//		unsigned int ind = p_int - x;
//		if (x_int != x[ind]) {  // x_int falls in the interval
//			y_int = y[ind]+(y[ind+1]-y[ind])/(x[ind+1]-x[ind])*(x_int-x[ind]);
//		} else {          // x_int falls on a grid point
//			y_int = y[ind];
//		}
//		return y_int;
//	} 
//} 

//double interp(double *x, double *y, double x_int, int N);
//template <typename T1, typename T2>
//T2 interp(const T1 *x, const T2 *y, T1 x_int, int N);

int main(void) {
	double x[6] = {0, 1, 2, 3, 4, 5};
	double y[6] = {0, 2, 4, 6, 8, 10};
	
	double x0[6] = {3.5, 4, -0.1, 5.1, 0, 5};
	
	for (int i=0; i < 6; i++) {
		double x_int = interp(x, y, x0[i], 6);
		cout<<"x0["<<i<<"], x_int: "<<x0[i]<<" "<<x_int<<endl;
	}
	
	cout<<"PI is "<<M_PI<<endl;
	cout<<"sqrt(2) is "<<pow(2, 0.5)<<endl;
	
	int r;
	string xsfile("sample_xs.txt");
	if ( (r = readxs(xsfile, xs_E, xs_sig)) < 0) {
		cout<<"Error in reading xs file!"<<endl;
		return r;
	}
	gridEtoV(xs_E, xs_v);

	cout<<"xs vector size: "<<xs_E.size()<<endl;
	cout<<"xs_E       "<<"xs_v      "<<"xs_sig       "<<endl;
	for (int i=0; i<xs_sig.size(); i++) {
		cout<<xs_E[i]<<"  "<<xs_v[i]<<"  "<<xs_sig[i]<<endl;
	}
	
	return 0;
}
