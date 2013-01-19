#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <random>

#include "parameters.h"
#include "util.h"
#include "xsdata.h"
#include "productcdf.h"

using namespace std;

// initialize PARAM::T here
namespace PARAM {
	double T = 300;
}
	
#include <algorithm>
#include <map>

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


namespace {
	// build erf table
	const int NERF = 12001;
	const double ULIMIT = 6.;
	const double INTERF = 2*ULIMIT/(NERF-1);
	double *erftb = new double[NERF];
}

void seterftb () {
	for (int i=0; i<NERF; i++) {
		erftb[i] =  erf(INTERF * i - ULIMIT);
	}
}

inline double myerf(double x) {
	int ind;
	double result;
	
	if (x >= -ULIMIT && x <= ULIMIT) {
		ind = (int) ((x + ULIMIT) / INTERF);	
		result = erftb[ind] + (x + ULIMIT - ind * INTERF) * (erftb[ind+1] - erftb[ind]);
	} else if (x < -ULIMIT) {
		result = -1;
	} else {
		result = 1;
	}
	
#ifdef DEBUG
	cout<<"x = "<<x<<", ind = "<<ind<<", result = "<<result<<", ULIMIT = "<<ULIMIT<<endl;
#endif	
	
	return result;
}

void delerftv() {
	delete [] erftb;
}

int main(void) {

//	double x[6] = {0, 1, 2, 3, 4, 5};
//	double y[6] = {0, 2, 4, 6, 8, 10};
//	
//	double x0[6] = {3.5, 4, -0.1, 5.1, 0, 5};
//	
//	for (int i=0; i < 6; i++) {
//		double x_int = interp(x, y, x0[i], 6);
//		cout<<"x0["<<i<<"], x_int: "<<x0[i]<<" "<<x_int<<endl;
//	}
//	
//	cout<<"PI is "<<M_PI<<endl;
//	cout<<"sqrt(2) is "<<pow(2, 0.5)<<endl;
//	
////	int r;
////	isotope U238;
////	string xsfile("sample_xs.txt");
////	if ( (r = U238.readxs(xsfile, U238.xs_E, U238.xs_sig)) < 0) {
////		cout<<"Error in reading xs file!"<<endl;
////		return r;
////	}
////	U238.gridEtoV(U238.xs_E, U238.xs_v);

////	cout<<"xs vector size: "<<U238.xs_E.size()<<endl;
////	cout<<"xs_E       "<<"xs_v      "<<"xs_sig       "<<endl;
////	for (int i=0; i<U238.xs_sig.size(); i++) {
////		cout<<U238.xs_E[i]<<"  "<<U238.xs_v[i]<<"  "<<U238.xs_sig[i]<<endl;
////	}
////	
////	
////	cout<<"vmin: "<<PARAM::vmin<<", vmax: "<<PARAM::vmax<<endl;

//	double alpha = 238*CONST::M_NUCLEON/2/CONST::K_BOLTZMANN/PARAM::T;
//	CDFmuvt cdf238(alpha);
//	CDFMB cdfMB(alpha);
//	const double LIMIT = 500;
//	const double INTVL = 1;
//	for (int i=0; i < 2*LIMIT/INTVL; i++) {
//		cdf238.grid.push_back(-LIMIT + i*INTVL);
//		cdfMB.grid.push_back(-LIMIT + i*INTVL);
//	}
//	cdf238.setcdf();
//	cdfMB.setcdf();
//	ofstream outfile("sample_cdf.out");
////	ofstream outfile4("sample_cdfMB.out");
//	for (auto it = cdf238.cdf.begin(); it != cdf238.cdf.end(); it++) {
//		outfile<<*it<<endl;
//	}
//	
//	default_random_engine e;
//	uniform_real_distribution<double> u(0,1);
//	ofstream outfile2("muvt_dist.out");	
//	ofstream outfile3("mu_vt_dist.out");
//	
//	for (auto i = 0; i < 100000; i++) {
////		cout<<u(e)<<" ";
//		outfile2<<cdf238.getx(u(e))<<endl;
//	}
//	
//	for (auto i = 0; i < 100000; i++) {
//		double vt = cdfMB.getx(u(e));
//		double mu = 2*u(e) - 1;
//		outfile3<<vt*mu<<endl;
//	}
	
	seterftb();
	
	ofstream outfile("erftable_test.out");	
	outfile<<setprecision(15);
		
	default_random_engine e;
	normal_distribution<double> u(0, 3);
	multimap<double, double> erfout;
	for (auto i = 0; i < 100000; i++) {
		double rnd = u(e);
		erfout.insert(pair<double, double>(rnd, myerf(rnd)));
	}
	
	for (multimap<double, double>::iterator it=erfout.begin(); it!=erfout.end(); ++it)
    outfile << (*it).first << "   " << (*it).second << '\n';
    
	delerftv();
	return 0;
}
