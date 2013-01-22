#include <iostream>
#include <fstream>
#include <cmath>

#include "xsdata.h"
#include "parameters.h"

using namespace std;

//vector<double> xs_E;
//vector<double> xs_v;
//vector<double> xs_sig;

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
				&& abs(xs_E[i] - xs_E[j]) <= delE) {
			i++;
		} else {
			xs_E_ref.push_back(xs_E[i]);
			xs_sig_ref.push_back(xs_sig[i]);
			j = i;
		}
	}
	
}
