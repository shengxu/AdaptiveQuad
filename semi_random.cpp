#include <iostream>
#include <fstream>
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
	const int NERF = 120001;
	const double ULIMIT = 6.;
	const double INTERF = 2*ULIMIT/(NERF-1);
	double *erftb = new double[NERF];
}

void seterftb() {
	for (int i=0; i<NERF; i++) {
		erftb[i] =  erf(INTERF * i - ULIMIT);
	}
}

void delerftb() {
	delete [] erftb;
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

#ifdef MYERF	
	seterftb();
#endif
//	cout<<"INTERF = "<<INTERF<<endl;
//	cout<<"INTERF * 5 - ULIMIT = "<<INTERF * 5 - ULIMIT<<endl;
//	for (int i=580; i < 620; i++) {
//		cout<<"i, erftb[i]: "<<i<<" "<<erftb[i]<<endl;	
//	}

	PARAM::T = atof(argv[4]);
		
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
	vector<double> xs_E_ref, xs_sig_ref, xs_sig_ave;
	// input: energy width, xs differencd
	U238.refinemesh(atof(argv[2]), atof(argv[3]), xs_E_ref, xs_sig_ref, xs_sig_ave);	
//	cout<<"input parameters: "<<argv[2]<<" "<<argv[3]<<endl;
	U238.xs_E = xs_E_ref;
	U238.xs_sig = xs_sig_ref;
	U238.gridEtoV(U238.xs_E, U238.xs_v);
	
	// sequence of energy points to broaden
	// uniform in lethargy simulating neutron slowing down with H2O as moderator
	double Ebegin = 1.95e4;  // upper bound to evaluate
	double Eend = 1.;    // lower bound
	int Npoints = 100000;  // number of equal-lethargy points
	double ksi = log(Ebegin/Eend)/Npoints;
	vector<double> Eseq;
	while (Ebegin >= Eend) {
		Eseq.push_back(Ebegin);
		Ebegin *= exp(-ksi);
	}

	ofstream outfile;
	outfile.open("fordebug.out", ios::out | ios::app);
	outfile<<setprecision(15);
//	outfile<<"size of xs_E: "<<U238.xs_E.size()<<", size of xs_sig_ave: "<<xs_sig_ave.size()<<endl;
//	for (auto i=0; i < xs_sig_ave.size(); i++) {
//		outfile<<U238.xs_E[i]<<"  "<<U238.xs_sig[i]<<"  "<<xs_sig_ave[i]<<endl;
//	}
	outfile<<atof(argv[2])<<"  "<<atof(argv[3])<<"  "<<U238.xs_E.size()<<endl;
	outfile.close();

	int A = 238;
	double alpha = CONST::M_NUCLEON*A/(2.*CONST::K_BOLTZMANN*PARAM::T);  // alpha, as used in cullen's method
//	int A = 236;
//	double alpha = CONST::M_NEUT*A/(2.*CONST::K_BOLTZMANN*PARAM::T);
	
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

#ifdef DEBUG
		ofstream outfile2("fordebug.out");
#endif	

	for (auto i=0; i < Eseq.size(); i++) {
		double Etmp = Eseq[i];
		double vtmp = EtoV(Etmp);

	for (auto i=0; i<Eseq.size(); i++) {
		double Etmp = Eseq[i];
		double vtmp = EtoV(Eseq[i]);
//	for (unsigned int i = 0; i < U238.xs_v.size(); i++) {	
//		double vtmp = U238.xs_v[i];
//		double Etmp = U238.xs_E[i];
//		double El = 0.5*CONST::M_NEUT*pow(U238.xs_v[i] - delv,2), Eu = 0.5*CONST::M_NEUT*pow(U238.xs_v[i] + delv,2);
		int indl, indu;
		vector<double>::iterator tmp;
		// lower_bound returns the first element that is not less than input value
		if ( (tmp = lower_bound(U238.xs_v.begin(), U238.xs_v.end(), vtmp - delv)) != U238.xs_v.begin() ) {
			indl = tmp - U238.xs_v.begin() -1;
		} else {
			indl = tmp - U238.xs_v.begin();
		}
		indu = upper_bound(U238.xs_v.begin(), U238.xs_v.end(), vtmp + delv) - U238.xs_v.begin();
//		cout<<"indl = "<<indl<<", indu = "<<indu<<endl;
		xs_brdn= 0;

#ifdef DEBUG
			outfile2<<"vtmp = "<<vtmp<<", delv = "<<delv<<", lower bound:"<<U238.xs_v[indl]<<", upper bound:"<<U238.xs_v[indu]<<endl;
			outfile2<<"lower bound: "<<VtoE(vtmp-delv)<<", upper bound: "<<VtoE(vtmp+delv)<<endl;
#endif
		
//		muvt = vtmp*(2./3. - 8./(9.*U238.xs_E[indl]/Etmp + 3.));
		muvt = vtmp*(U238.xs_E[indl]/Etmp - 1)/(U238.xs_E[indl]/Etmp + 1);
//		muvt = vtmp*(sqrt(U238.xs_E[indl]/Etmp) - 1);


#ifdef MYERF
		cdf_p = 0.5*(1 + myerf(sqalpha*muvt));
#else
		cdf_p = 0.5*(1 + erf(sqalpha*muvt));
#endif
		for (int j=indl+1; j <= indu; j++) {
//			muvt = U238.xs_v[i]*(sqrt(U238.xs_E[j]/U238.xs_E[i]) - 1);

//			muvt = vtmp*(2./3. - 8./(9.*U238.xs_E[j]/Etmp + 3.));
			muvt = vtmp*(U238.xs_E[j]/Etmp - 1)/(U238.xs_E[j]/Etmp + 1);
//			muvt = vtmp*(sqrt(U238.xs_E[j]/Etmp) - 1);

//			muvt = 0.5*U238.xs_v[i]*(U238.xs_E[j]/U238.xs_E[i] - vT_over_v - 1);
//			muvt = 0.5*U238.xs_v[i]*(sqrt(2*U238.xs_E[j]/U238.xs_E[i] - 1) - 1);
#ifdef MYERF
			cdf = 0.5*(1 + myerf(sqalpha*muvt));
#else
			cdf = 0.5*(1 + erf(sqalpha*muvt));
#endif			

//#ifdef DEBUG
//			outfile2<<U238.xs_E[j-1]<<" -- "<<U238.xs_E[j]<<"  "<<U238.xs_sig[j-1]<<"  "<<xs_sig_ave[j-1]<<endl;
//			outfile2<<cdf_p<<"  "<<cdf<<endl;
//#endif


//			xs_brdn += U238.xs_sig[j] * (cdf - cdf_p) * U238.xs_v[j]/U238.xs_v[i];
//			xs_brdn += 0.5*(U238.xs_sig[j] + U238.xs_sig[j-1])* (cdf - cdf_p)* 0.5*(U238.xs_v[j] + U238.xs_v[j-1])/U238.xs_v[i];
//			xs_ave = 0.5*(U238.xs_sig[j] + U238.xs_sig[j-1]);
			xs_ave = xs_sig_ave[j-1];
//			vave = sqrt((U238.xs_E[j] + U238.xs_E[j-1])*1.e-6/CONST::M_NEUT);
//			cout<<"vave = "<<vave<<", U238.xs_v[i] = "<<U238.xs_v[i]<<endl;
//			xs_brdn += xs_ave * (cdf - cdf_p) * vave/U238.xs_v[i];
			xs_brdn += xs_ave * (cdf - cdf_p);
			cdf_p = cdf;
		}
//		cout<<"v = "<<v<<", delv = "<<delv<<", sig(v) = "<<U238.xs_sig[i]<<endl;


		cout<<Etmp<<"  "<<vtmp<<"  "<<xs_brdn<<"  "<<indu - indl + 1<<endl;
//		cout<<"Broadened xs at 0.025 eV at 300K: "<<xs_brdn<<endl;	
	}

#ifdef DEBUG	
	outfile2.close();
#endif
	
#ifdef MYERF	
	delerftb();
#endif	

}
