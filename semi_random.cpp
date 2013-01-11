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

int main(int argc, char **argv) {
	
	int r;
	isotope U238;
	string xsfile(argv[1]);
	// string xsfile("pendf_0K_102");
	if ( (r = U238.readxs(xsfile, U238.xs_E, U238.xs_sig)) < 0) {
		cout<<"Error in reading xs file!"<<endl;
		return r;
	}
	U238.gridEtoV(U238.xs_E, U238.xs_v);
	
	int A = 238;
	double alpha = CONST::M_NUCLEON*A/(2.*CONST::K_BOLTZMANN*PARAM::T);
	double sqalpha = sqrt(alpha);
	double delv = 4./sqrt(alpha);
	double xs_brdn;
	double cdf_p, cdf;
	double muvt;
	
//	cout<<"xs_E       "<<"xs_v      "<<"xs_sig       "<<endl;
//	for (unsigned int i=0; i < U238.xs_v.size(); i += 10) {
	for (unsigned int i = 140; i < U238.xs_v.size(); i += 10) {
//	for (unsigned int i = 640; i == 640; i += 10) {
//		cout<<"i = "<<i<<endl;
//		double El = 0.5*CONST::M_NEUT*pow(U238.xs_v[i] - delv,2), Eu = 0.5*CONST::M_NEUT*pow(U238.xs_v[i] + delv,2);
		int indl = lower_bound(U238.xs_v.begin(), U238.xs_v.end(), U238.xs_v[i] - delv) - U238.xs_v.begin();
		int indu = upper_bound(U238.xs_v.begin(), U238.xs_v.end(), U238.xs_v[i] + delv) - U238.xs_v.begin();
//		cout<<"indl = "<<indl<<", indu = "<<indu<<endl;
		xs_brdn= 0;
		cdf_p = 0;
		for (int j=indl; j <= indu; j++) {
			muvt = U238.xs_v[i]*(sqrt(U238.xs_E[j]/U238.xs_E[i]) - 1);
			cdf = 0.5*(1 + erf(sqalpha*muvt));
			xs_brdn += U238.xs_sig[j] * (cdf - cdf_p);
			cdf_p = cdf;
		}
//		cout<<"v = "<<v<<", delv = "<<delv<<", sig(v) = "<<U238.xs_sig[i]<<endl;

		cout<<U238.xs_E[i]<<"  "<<U238.xs_v[i]<<"  "<<xs_brdn<<endl;
//		cout<<"Broadened xs at 0.025 eV at 300K: "<<xs_brdn<<endl;	
	}

}
