#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include "constants.h"
using namespace std;

// set max and min energy
double emin = 1e-11;
double emax = 20.0;
double vmin = sqrt(2.*1e-11/CONST::M_NEUT);
double vmax = sqrt(2.*20./CONST::M_NEUT);

//upper bound of velocity in M-B distribution, relative to Vmp
double upper_bound_v = 20.0;

//TODO: needs to be initialized somewhere else
double T = 300;

#endif //PARAMETERS_H
