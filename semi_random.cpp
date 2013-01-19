#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>

#include "util.h"
#include "parameters.h"
#include "xsdata.h"

using namespace std;

// initialize PARAM::T here
namespace PARAM {
	double T = 300;
}


namespace {
	// build erf table
	const int NERF = 12001;
	const double ULIMIT = 6.;
	const double INTERF = 2*ULIMIT/(NERF-1);
	double *erftb = new double[NERF];
}

inline void seterftb () {
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

int main(int argc, char **argv) {
	
	seterftb();
//	cout<<"INTERF = "<<INTERF<<endl;
//	cout<<"INTERF * 5 - ULIMIT = "<<INTERF * 5 - ULIMIT<<endl;
//	for (int i=580; i < 620; i++) {
//		cout<<"i, erftb[i]: "<<i<<" "<<erftb[i]<<endl;	
//	}
	
	int r;
	isotope U238;
	string xsfile(argv[1]);
	// string xsfile("pendf_0K_102");
	
	
	// set the precision of output
	cout<<setprecision(15);

	if ( (r = U238.readxs(xsfile, U238.xs_E, U238.xs_sig)) < 0) {
		cout<<"Error in reading xs file!"<<endl;
		return r;
	}
	U238.gridEtoV(U238.xs_E, U238.xs_v);
	
	
	int A = 238;
	double alpha = CONST::M_NUCLEON*A/(2.*CONST::K_BOLTZMANN*PARAM::T);  // alpha, as used in cullen's method
	double sqalpha = sqrt(alpha);   // square root of alpha
	double delv = 4./sqrt(alpha);
	double vpsq = 1/alpha;  // square of most probable velocity
//	cout<<"vpsq = "<<vpsq<<endl;
	double xs_brdn;
	double cdf_p, cdf;
	double muvt, xs_ave, vave;
	
//	cout<<"xs_E       "<<"xs_v      "<<"xs_sig       "<<endl;

//	for (unsigned int i = 800000; i < U238.xs_v.size(); i++) {
//		U238.xs_sig[i] *= 10;
//	}
	
	for (unsigned int i = 0; i < U238.xs_v.size(); i += 10) {
//	for (unsigned int i = 140; i < 2000; i += 10) {
//	for (unsigned int i = 640; i == 640; i += 10) {
//		cout<<"i = "<<i<<endl;
//		double El = 0.5*CONST::M_NEUT*pow(U238.xs_v[i] - delv,2), Eu = 0.5*CONST::M_NEUT*pow(U238.xs_v[i] + delv,2);
		int indl = lower_bound(U238.xs_v.begin(), U238.xs_v.end(), U238.xs_v[i] - delv) - U238.xs_v.begin();
		int indu = upper_bound(U238.xs_v.begin(), U238.xs_v.end(), U238.xs_v[i] + delv) - U238.xs_v.begin();
//		cout<<"indl = "<<indl<<", indu = "<<indu<<endl;
		xs_brdn= 0;
		cdf_p = 0;
		double vT_over_v = vpsq/pow(U238.xs_v[i], 2);
		for (int j=indl; j <= indu; j++) {
//			muvt = U238.xs_v[i]*(sqrt(U238.xs_E[j]/U238.xs_E[i]) - 1);
			muvt = U238.xs_v[i]*(2./3. - 8./(9.*U238.xs_E[j]/U238.xs_E[i] + 3.));
//			muvt = 0.5*U238.xs_v[i]*(U238.xs_E[j]/U238.xs_E[i] - vT_over_v - 1);
//			muvt = 0.5*U238.xs_v[i]*(sqrt(2*U238.xs_E[j]/U238.xs_E[i] - 1) - 1);
			cdf = 0.5*(1 + myerf(sqalpha*muvt));
//			xs_brdn += U238.xs_sig[j] * (cdf - cdf_p) * U238.xs_v[j]/U238.xs_v[i];
//			xs_brdn += 0.5*(U238.xs_sig[j] + U238.xs_sig[j-1])* (cdf - cdf_p)* 0.5*(U238.xs_v[j] + U238.xs_v[j-1])/U238.xs_v[i];
			xs_ave = 0.5*(U238.xs_sig[j] + U238.xs_sig[j-1]);
//			vave = sqrt((U238.xs_E[j] + U238.xs_E[j-1])*1.e-6/CONST::M_NEUT);
//			cout<<"vave = "<<vave<<", U238.xs_v[i] = "<<U238.xs_v[i]<<endl;
//			xs_brdn += xs_ave * (cdf - cdf_p) * vave/U238.xs_v[i];
			xs_brdn += xs_ave * (cdf - cdf_p);
			cdf_p = cdf;
		}
//		cout<<"v = "<<v<<", delv = "<<delv<<", sig(v) = "<<U238.xs_sig[i]<<endl;

		cout<<U238.xs_E[i]<<"  "<<U238.xs_v[i]<<"  "<<xs_brdn<<"  "<<indu - indl + 1<<endl;
//		cout<<"Broadened xs at 0.025 eV at 300K: "<<xs_brdn<<endl;	
	}
	
	delete [] erftb;

}
