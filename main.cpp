#include <iostream>
#include <iomanip>
#include <vector>

#include "adapsimpsint.h"

using namespace std;

// initialize PARAM::T here
namespace PARAM {
	double T = 300;
}

int main() {
	Fun f;
	AdapSimps simpson1(f);
	cout<<"Adaptive Simpson Int:"<<setprecision(7)<<simpson1(0, 2, 1e-7, 20)<<endl;	
	
	int r;
	isotope U238;
	string xsfile("pendf_0K_102");
	if ( (r = U238.readxs(xsfile, U238.xs_E, U238.xs_sig)) < 0) {
		cout<<"Error in reading xs file!"<<endl;
		return r;
	}
	U238.gridEtoV(U238.xs_E, U238.xs_v);

	cout<<"xs vector size: "<<U238.xs_E.size()<<endl;
	cout<<"xs_E       "<<"xs_v      "<<"xs_sig       "<<endl;
//	for (int i=0; i<U238.xs_sig.size(); i++) {
	for (int i=0; i<10; i++) {
		cout<<U238.xs_E[i]<<"  "<<U238.xs_v[i]<<"  "<<U238.xs_sig[i]<<endl;
	}
	
	int A = 238;
	double v = 2.2e3;
	sigv f2(A, v, U238.xs_E, U238.xs_v, U238.xs_sig);
	AdapSimps simpson2(f2);
	double alpha = CONST::M_NUCLEON*A/(2.*CONST::K_BOLTZMANN*PARAM::T);
	double delv = 4./sqrt(alpha);
	double xs_brdn = 0;
	cout<<"xs_E       "<<"xs_v      "<<"xs_sig       "<<endl;
	for (int i=0; i < U238.xs_v.size(); i += 100) {
		v = U238.xs_v[i];
		f2.v = v;
		xs_brdn = simpson2(v - delv, v + delv, 1e-2, 80);
		f2.v = -v;
		xs_brdn -= simpson2(0, delv, 1e-2, 80);
		cout<<U238.xs_E[i]<<"  "<<U238.xs_v[i]<<"  "<<xs_brdn<<endl;
//		cout<<"Broadened xs at 0.025 eV at 300K: "<<xs_brdn<<endl;	
	}
}
