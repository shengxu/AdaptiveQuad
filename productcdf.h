#ifndef PRODUCTPDF_H
#define PRODUCTPDF_H

#include <vector>
#include "parameters.h"

class CDFmuvt {
	public:
		double alpha;
		std::vector<double> grid;
		std::vector<double> cdf;
		
		CDFmuvt(const double _alpha) : alpha(_alpha) {}
		void setcdf();
		double getx(double cdfin) const;
};


class CDFMB : public CDFmuvt {
	public:
		CDFMB(const double _alpha) : CDFmuvt(_alpha) {}
		void setcdf();
};

#endif  // PRODUCTPDF_H
