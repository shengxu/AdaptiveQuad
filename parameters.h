#ifndef PARAMETERS_H
#define PARAMETERS_H

//// constants
//namespace CONST {
//	//real(8), parameter:: ZERO = 0.0_8
//	//real(8), parameter:: ONE = 1.0_8
//	extern double const PI;  //= 3.1415926535898;
//	extern double const K_BOLTZMANN;  //= 8.617342e-11;     // unit: MeV/K
//	extern double const M_NEUT;  //= 939.565378/std::pow(3.0e+8, 2);  //unit: MeV/((m/s)^2)
//	extern double const M_NUCLEON;  //= 931.494061/std::pow(3.0e+8, 2);  //unit: MeV/((m/s)^2)
//	extern double const C;  //= 3.e8;  //unit: m/s
//}

// constants
namespace CONST {
	//real(8), parameter:: ZERO = 0.0_8
	//real(8), parameter:: ONE = 1.0_8
	const double PI = 3.1415926535898;
	const double SQRT_PI = std::sqrt(PI);
	const double K_BOLTZMANN = 8.617342e-11;     // unit: MeV/K
	const double M_NEUT = 939.565378/std::pow(3.0e+8, 2);  //unit: MeV/((m/s)^2)
	const double M_NUCLEON = 931.494061/std::pow(3.0e+8, 2);  //unit: MeV/((m/s)^2)
	const double C = 3.e8;  //unit: m/s
}

//namespace PARAM {
//	extern double const eps;  //= 1.0e-16;     // for devision of small numbers
//	// set max and min energy
//	extern double const emin;  //= 1e-11;
//	extern double const emax;  //= 20.0;
//	extern double const vmin;  //= std::sqrt(2.*1e-11/CONST::M_NEUT);
//	extern double const vmax;  //= std::sqrt(2.*20./CONST::M_NEUT);

//	//upper bound of velocity in M-B distribution, relative to Vmp
//	extern double const upper_bound_v;  //= 20.0;

//	//TODO: needs to be initialized somewhere else
//	extern double T;  //= 300;
//}

namespace PARAM {
	const double eps = 1.0e-16;     // for devision of small numbers
	// set max and min energy
	const double emin = 1e-11;
	const double emax = 20.0;
	const double vmin = std::sqrt(2.*1e-11/CONST::M_NEUT);
	const double vmax = std::sqrt(2.*20./CONST::M_NEUT);

	//upper bound of velocity in M-B distribution, relative to Vmp
	const double upper_bound_v = 20.0;

	//TODO: needs to be initialized somewhere else
//	double T = 300;
	extern double T;
}

#endif //PARAMETERS_H
