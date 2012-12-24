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
		
		std::vector<double> xs_E;
		std::vector<double> xs_v;
		std::vector<double> xs_sig;
};

#endif
