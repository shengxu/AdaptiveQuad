//Trapzint.cpp
#include <vector>
#include <cmath>
#include <assert.h>
#include "adapsimpsint.h"

using namespace std;

double AdapSimps::operator ()(double a,double b,double eps, int N) const {
	vector<struct info> infostack;
	struct info tmp;
	struct info tmp2;
	double result = 0;
	int cnt = 0;

	tmp.a = a;
	tmp.h = (b-a)/2;
	tmp.F[0] = f(a);
	tmp.F[1] = f(a + tmp.h);
	tmp.F[2] = f(b);
	tmp.TOL = 10*eps;
	tmp.S = tmp.h*(tmp.F[0] + 4*tmp.F[1] + tmp.F[2])/3;
	tmp.L = 1;

	infostack.push_back(tmp);

	while (!infostack.empty()) {
		++cnt;
		tmp = infostack.back();
		#ifdef DEUBG
			cout<<"tmp.L = "<<tmp.L<<", tmp.TOL = "<<tmp.TOL<<endl;
		#endif
		infostack.pop_back();
		if (cnt == 1) {
			assert(infostack.empty());
		}

		double FD = f(tmp.a + tmp.h/2);
		double FE = f(tmp.a + 3*tmp.h/2);
		double S1 = tmp.h*(tmp.F[0] + 4*FD + tmp.F[1])/6;
		double S2 = tmp.h*(tmp.F[1] + 4*FE + tmp.F[2])/6;

		#ifdef DEUBG
			cout<<"level = "<<tmp.L<<", S1, S2, and tmp.S: "<<S1<<" "<<S2<<" "<<tmp.S<<endl;		
		#endif
		if (abs(S1 + S2 - tmp.S) < tmp.TOL) {
			result += (S1 + S2);
			#ifdef DEUBG
				cout<<"level = "<<tmp.L<<", result = "<< result<<endl;
			#endif
		} else {
			if (tmp.L < N) {
				tmp2.a = tmp.a + tmp.h;
				tmp2.F[0] = tmp.F[1];
				tmp2.F[1] = FE;
				tmp2.F[2] = tmp.F[2];
				tmp2.h = tmp.h/2;
				tmp2.TOL = tmp.TOL/2;
				tmp2.S = S2;
				tmp2.L = tmp.L + 1;
				infostack.push_back(tmp2);
				tmp2.a = tmp.a;
				tmp2.F[0] = tmp.F[0];
				tmp2.F[1] = FD;
				tmp2.F[2] = tmp.F[1];
	//			tmp2.h = tmp.h/2;
	//			tmp2.TOL = tmp.TOL/2;
				tmp2.S = S1;
	//			tmp2.L = tmp.L + 1;
				infostack.push_back(tmp2);
			} else {
				cout<<"Level exceeded!"<<endl;
				result += (S1 + S2);
				#ifdef DEUBG
					cout<<"S1: "<<S1<<", S2: "<<S2<<", tmp.S: "<<tmp.S<<", TOL: "<<tmp.TOL<<endl;
				#endif
//				return 0;
			}
		}
	}

	return result;
}
//End of file Trapzint.cpp
