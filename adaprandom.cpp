#include <iostream>
#include <iomanip>
#include <vector>
#include <random>

#include "parameters.h"
#include "xsdata.h"
#include "productcdf.h"
#include "util.h"

using namespace std;

namespace {
	int NMAX = 25;
}

static inline double E2v(double E) {
	return sqrt(2*E*1.e-6/CONST::M_NEUT);
}
// E: eV
double adaprandomint(double E, const CDFmuvt & cdf, const std::vector<double> &xs_E, const std::vector<double> &xs_sig, double eps) {

	size_t cnt = 0;
	bool stat = false;
//	double Erel = 0;
	double xs_brdn = 0;
	double xs_brdn_p;
	
	double v = E2v(E);
	
	#ifdef DEBUG
		cout<<"input: v = "<<v<<endl;
	#endif
	
	// initial evaluation
	int L0 = 2;
	int N0 = std::pow(2, L0);
	double intv = 1.0/N0;

	for (auto i = 1; i < N0; i++) {
		double vtmu = cdf.getx(i*intv);
		double Erel = E*(1 + 2*vtmu/v);
		double vrel = E2v(Erel);
		double xs_rel = interp(xs_E.begin(), xs_sig.begin(), Erel, xs_E.size());
		xs_brdn += vrel*xs_rel;
		#ifdef DEBUG
			cout<<"Erel = "<<Erel<<", xs_rel = "<<xs_rel<<endl;
		#endif
		cnt++;
	}
	
	xs_brdn /= (N0 - 1);
	xs_brdn_p = xs_brdn;
	intv /= 2;
	L0++;
	while (1) {
		xs_brdn *= (N0 - 1);
		for (auto i = 0; i < N0; i++) {
			double vtmu = cdf.getx((2*i + 1)*intv);
			double Erel = E*(1 + 2*vtmu/v);
			double vrel = E2v(Erel);
			double xs_rel = interp(xs_E.begin(), xs_sig.begin(), Erel, xs_E.size());				
			xs_brdn += vrel*xs_rel;
			#ifdef DEBUG
				cout<<"Erel = "<<Erel<<", xs_rel = "<<xs_rel<<endl;
			#endif
			cnt++;
		}
		if (L0 >= NMAX) {
			cout<<cnt<<"       ";
			return xs_brdn/v;
		} else {
			xs_brdn /= (2*N0 - 1);
			if (abs((xs_brdn - xs_brdn_p)/(xs_brdn + PARAM::eps)) < eps) {
				if (stat) {
					cout<<cnt<<"       ";
					return xs_brdn/v;
				}  else {
					stat = true;
				}
			}
		}
		xs_brdn_p = xs_brdn;
		N0 *= 2;
		intv /= 2;
		L0++;
	}
}
