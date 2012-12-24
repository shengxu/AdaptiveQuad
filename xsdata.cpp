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
