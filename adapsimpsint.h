#ifndef ADAPSIMPSONINT_H
#define ADAPSIMPSONINT_H

#include <iostream>
#include <cmath>
#include <vector>
#include "util.h"
#include "parameters.h"
#include "xsdata.h"

class F	
{
public:
	virtual double operator ()(double x) const=0;
};

class Fun:public F {
public:
	double operator()(double x) const {
	      return log(1.0+x)/(1.0+x*x);
	}
};

class sigv:public F {
public:
	int A;
	double v;
	double alpha;
	std::vector<double> xs_E;
	std::vector<double> xs_v;
	std::vector<double> xs_sig;
	
	// constructor
	sigv(int _A, double _v, std::vector<double> &_xs_E, std::vector<double> &_xs_v, std::vector<double> &_xs_sig) : A(_A), v(_v), xs_E(_xs_E), xs_v(_xs_v), xs_sig(_xs_sig) {}	
	
	void setalpha(double T) {
		alpha = CONST::M_NUCLEON*A/(2.*CONST::K_BOLTZMANN*T);
	}
	
	double operator()(double V) const {
		return (double) std::pow(alpha/CONST::PI, 0.5)*std::pow(V/v, 2)*std::exp(-alpha*std::pow(V-v,2))*interp(xs_v.begin(), xs_sig.begin(), V, xs_v.size());
	}
};


class Integ {
public:
	virtual double operator ()(double a,double b,double eps, int N) const=0;
};

class AdapSimps:public Integ {
public:
	AdapSimps(const F &pf) : f(pf) {}
	double operator ()(double a, double b,double eps, int N) const;
private:
	const F &f;
};	

struct info {
	double a;
	double F[3];  // [F(a), F(c), F(b)]
	double h;
	double TOL;
	double S;
	int L;
};

#endif // ADAPSIMPSONINT_H
