#ifndef XSDATA_H
#define XSDATA_H

#include <string>
#include <vector>

#include "parameters.h"

//extern vector<double> xs_E;
//extern vector<double> xs_v;
//extern vector<double> xs_sig;

class isotope {
	public:
		int readxs(const std::string &filename, std::vector<double> &xs_E, std::vector<double> &xs_sig);
		void gridEtoV(std::vector<double> &xs_E, std::vector<double> &xs_v);
		void refinemesh(const double delE, const double relxs, std::vector<double> &xs_E_ref, std::vector<double> &xs_sig_ref);
		void refinemesh(const double delE, const double relxs, std::vector<double> &xs_E_ref, std::vector<double> &xs_sig_ref, std::vector<double> &xs_sig_ave);
		
		std::vector<double> xs_E;
		std::vector<double> xs_v;
		std::vector<double> xs_sig;
};

#endif
