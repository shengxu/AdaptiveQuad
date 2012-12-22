#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

using namespace std;
namespace CONST {
	//real(8), parameter:: ZERO = 0.0_8
	//real(8), parameter:: ONE = 1.0_8
	double PI = 3.1415926535898;
	double K_BOLTZMANN = 8.617342e-11;     // unit: MeV/K
	double M_NEUT = 939.565378/pow(3.0e+8, 2);  //unit: MeV/((m/s)^2)
	double M_NUCLEON = 931.494061/pow(3.0e+8, 2);  //unit: MeV/((m/s)^2)
	double C = 3.e8;  //unit: m/s
	double EPS = 1.0e-16;     // for devision of small numbers
}

#endif  // CONSTANTS_H
