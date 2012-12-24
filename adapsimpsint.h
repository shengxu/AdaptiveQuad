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
	
	// constructor
	sigv(int _A, double _v) : A(_A), v(_v) {}
	
	double operator()(double V) const {
		double alpha = CONST::M_NUCLEON*A/(2.*CONST::K_BOLTZMANN*T);
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
