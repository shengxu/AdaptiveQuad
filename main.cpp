#include <iostream>
#include <iomanip>
#include <vector>

#include "adapsimpsint.h"

using namespace std;

// initialize PARAM::T here
namespace PARAM {
	double T = 300;
}

int main(int argc, char **argv) {
	Fun f;
	AdapSimps simpson1(f);
	cout<<"Adaptive Simpson Int:"<<setprecision(15)<<simpson1(0, 2, 1e-7, 20)<<endl;	
	
	int r;
	isotope U238;
	string xsfile(argv[1]);
	// string xsfile("pendf_0K_102");
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
//	f2.setalpha(PARAM::T);
	AdapSimps simpson2(f2);
	double alpha = CONST::M_NUCLEON*A/(2.*CONST::K_BOLTZMANN*PARAM::T);
	f2.alpha = alpha;
	double delv = 4./sqrt(alpha);
	double xs_brdn = 0;
	cout<<"xs_E       "<<"xs_v      "<<"xs_sig       "<<endl;
	for (unsigned int i=0; i < U238.xs_v.size(); i += 10) {
//	for (unsigned int i=380180; i == 380180; i += 10) {
//		cout<<"i = "<<i<<endl;
		v = U238.xs_v[i];
		f2.v = v;
//		cout<<"v = "<<v<<", delv = "<<delv<<", sig(v) = "<<U238.xs_sig[i]<<endl;
		xs_brdn = simpson2(max(v - delv, 0.0), v + delv, 1e-4, 20);
		if (v < delv) {
			f2.v = -v;
			xs_brdn -= simpson2(0, delv - v, 1e-4, 20);
		} else {
			cout<<"0      ";		
		}
		cout<<U238.xs_E[i]<<"  "<<U238.xs_v[i]<<"  "<<xs_brdn<<endl;
//		cout<<"Broadened xs at 0.025 eV at 300K: "<<xs_brdn<<endl;	
	}

//		int i = 530930;
//		v = U238.xs_v[i];
//		f2.v = v;
//		xs_brdn = simpson2(max(v - delv, 0.0), v + delv, 1e-3, 10);
//		cout<<"xs_brdn_1 = "<<xs_brdn<<endl;
//		if (v < delv) {
//			f2.v = -v;
//			xs_brdn -= simpson2(0, delv - v, 1e-3, 10);
//			cout<<"xs_brdn_2 = "<<xs_brdn<<endl;
//		} else {
//			cout<<"0      ";		
//		}
//		cout<<U238.xs_E[i]<<"  "<<U238.xs_v[i]<<"  "<<xs_brdn<<endl;
////		cout<<"Broadened xs at 0.025 eV at 300K: "<<xs_brdn<<endl;	
}
