#include <iostream>
#include <iomanip>
#include <vector>

#include "adapsimpsint.h"
#include "util.h"

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
	double alpha = A*CONST::M_NUCLEON/2/CONST::K_BOLTZMANN/PARAM::T;
	CDFmuvt cdf238(alpha);
//	CDFMB cdfMB(alpha);
	const double LIMIT = 500;
	const double INTVL = 1;
	for (int i=0; i < 2*LIMIT/INTVL; i++) {
		cdf238.grid.push_back(-LIMIT + i*INTVL);
//		cdfMB.grid.push_back(-LIMIT + i*INTVL);
	}
	cdf238.setcdf();
	
	cout<<setprecision(15);
	
	double xs_brdn = 0;
//	cout<<"xs_E       "<<"xs_v      "<<"xs_sig       "<<endl;
	for (unsigned int i=100; i < U238.xs_E.size(); i += 10) {
//	for (unsigned int i=100; i <= 20000; i += 10) {
//	for (unsigned int i=19960; i <= 19960; i += 10) {
//	for (unsigned int i=0; i <= 100; i += 10) {
		xs_brdn = adaprandomint(U238.xs_E[i], cdf238, U238.xs_E, U238.xs_sig, 1e-3);
		cout<<U238.xs_E[i]<<"  "<<U238.xs_v[i]<<"  "<<xs_brdn<<endl;
//		cout<<"Broadened xs at 0.025 eV at 300K: "<<xs_brdn<<endl;	
	}	
}
