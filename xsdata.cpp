#include <iostream>
#include <fstream>
#include <cmath>

#include "xsdata.h"
#include "parameters.h"

using namespace std;

//vector<double> xs_E;
//vector<double> xs_v;
//vector<double> xs_sig;

typedef struct {
	double E;
	double sig;
} XSpoint;

inline void findinterp(const double sig_b, const XSpoint& left, const XSpoint& right, const double percent, XSpoint& result) {
	if(right.sig >= sig_b) {
		result.sig = sig_b*(1+percent);
	} else {
		result.sig = sig_b*(1-percent);
	}
	
	result.E = left.E + (right.E - left.E)/(right.sig - left.sig)*(result.sig - left.sig);
}

int isotope::readxs(const string &filename, vector<double> &xs_E, vector<double> &xs_sig) {
	double E, sig;
	
	ifstream infile(filename.c_str());
	
	//Always test the file open.
	if(!infile) {
		cerr<<"Error opening output file"<<endl;
//		system("pause");
		return -1;
	}
	
	// Start reading data
//	while (!infile.eof()){
//		infile >> E;
//		infile >> sig;
	while (infile>>E>>sig){
		xs_E.push_back(E);
		xs_sig.push_back(sig);
	}
	
	infile.close();
	return 0;
}

void isotope::gridEtoV(vector<double> &xs_E, vector<double> &xs_v) {
	for (auto it=xs_E.begin(); it<xs_E.end(); it++) {
		xs_v.push_back(sqrt(2*(*it)*1.e-6/CONST::M_NEUT));
	}
}

void isotope::refinemesh(const double delE, const double relxs, vector<double> &xs_E_ref, vector<double> &xs_sig_ref) {
//	xs_E_ref.push_back(xs_E[0]);
//	xs_sig_ref.push_back(xs_sig[0]);
//	int j = 0;   // to record last index pushed
//	double E_ave = xs_E[0];        // to record average energy
//	double sig_ave = xs_sig[0];    // to record average xs
//	
//	for (int i = 1; i < xs_sig.size(); i++) {
//		if (abs(xs_sig[i] - xs_sig[j])/xs_sig[j] <= relxs 
//				&& abs(xs_E[i] - xs_E[j]) <= delE) {
//			E_ave += xs_E[i];
//			sig_ave += xs_sig[i];
//		} else {
//			int k = i - j;
//			sig_ave /= k;
//			xs_sig_ref.push_back(sig_ave);
//			E_ave /= k;
//			xs_E_ref.push_back(E_ave);
//			j = i;
//			E_ave = xs_E[i];
//			sig_ave = xs_sig[i];
//		}
//	}

	xs_E_ref.push_back(xs_E[0]);
	xs_sig_ref.push_back(xs_sig[0]);
	int j = 0;
	
	for (unsigned int i = 1; i < xs_sig.size(); i++) {
		if (abs(xs_sig[i] - xs_sig[j])/xs_sig[j] <= relxs
				&& abs(xs_E[i] - xs_E[j]) <= delE*sqrt(xs_E[j]/6.67) ) {
			i++;
		} else {
			xs_E_ref.push_back(xs_E[i]);
			xs_sig_ref.push_back(xs_sig[i]);
			j = i;
		}
	}
}
	
void isotope::refinemesh(const double delE, const double relxs, vector<double> &xs_E_ref, vector<double> &xs_sig_ref, vector<double> &xs_sig_ave) {
//	xs_E_ref.push_back(xs_E[0]);
//	xs_sig_ref.push_back(xs_sig[0]);
//	int j = 0;   // to record last index pushed
//	double E_ave = xs_E[0];        // to record average energy
//	double sig_ave = xs_sig[0];    // to record average xs
//	
//	for (int i = 1; i < xs_sig.size(); i++) {
//		if (abs(xs_sig[i] - xs_sig[j])/xs_sig[j] <= relxs 
//				&& abs(xs_E[i] - xs_E[j]) <= delE) {
//			E_ave += xs_E[i];
//			sig_ave += xs_sig[i];
//		} else {
//			int k = i - j;
//			sig_ave /= k;
//			xs_sig_ref.push_back(sig_ave);
//			E_ave /= k;
//			xs_E_ref.push_back(E_ave);
//			j = i;
//			E_ave = xs_E[i];
//			sig_ave = xs_sig[i];
//		}
//	}

	XSpoint tmp;
	// first element
	double El = xs_E[0], sigl = xs_sig[0];
	double Er, sigr;
	
	xs_E_ref.push_back(xs_E[0]);
	xs_sig_ref.push_back(xs_sig[0]);
	unsigned int i;
	double xs_ave = 0;
	bool stat = false;  // indicate left boundary
	
	for (i = 1; i < xs_sig.size(); i++) {
		bool stat1 = abs(xs_sig[i] - sigl)/sigl <= relxs;
		double delEnew = delE*sqrt(El/6.67);
		bool stat2 = abs(xs_E[i] - El) <= delEnew;
		if (stat1 && stat2) {
			if (stat) {
				xs_ave += 0.5*(xs_sig[i-1] + xs_sig[i])*(xs_E[i] - xs_E[i-1]);
			} else {
				xs_ave += 0.5*(sigl + xs_sig[i])*(xs_E[i] - El);
				stat = true;
			}
		} else {
			if (!stat2) { // energy criterion violated
//				findinterp(sigl, {xs_E[i-1], xs_sig[i-1]}, {xs_E[i], xs_sig[i]}, relxs, tmp);
//				Er = tmp.E; sigr = tmp.sig;
				Er = El + delEnew;
				sigr = xs_sig[i-1] + (xs_sig[i] -xs_sig[i-1])/(xs_E[i] - xs_E[i-1])*(Er - xs_E[i-1]);
				xs_E_ref.push_back(Er);
				xs_sig_ref.push_back(sigr);
				if (stat) {
					xs_ave += 0.5*(xs_sig[i-1] + sigr)*(Er - xs_E[i-1]);
				} else {
					xs_ave += 0.5*(sigl + sigr)*(Er - El);
				}
				xs_sig_ave.push_back(xs_ave/(Er - El));
				El = Er; sigl = sigr;
				stat = false;
				xs_ave = 0;
				i--;
			} else {   // ignore if xs criterion violated
				Er = xs_E[i]; sigr = xs_sig[i];
				xs_E_ref.push_back(Er);
				xs_sig_ref.push_back(sigr);
				if (stat) {
					xs_ave += 0.5*(xs_sig[i-1] + sigr)*(Er - xs_E[i-1]);
				} else {
					xs_ave += 0.5*(sigl + sigr)*(Er - El);
				}
				xs_sig_ave.push_back(xs_ave/(Er - El));
				El = Er; sigl = sigr;
				stat = false;
				xs_ave = 0;
			}
			
			
		}		
	}	
	
	// last element			
	xs_E_ref.push_back(xs_E[i]);
	xs_sig_ref.push_back(xs_sig[i]);
	xs_sig_ave.push_back(xs_ave/(xs_E[i] - El));
	
}
