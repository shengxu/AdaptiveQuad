//Trapzint.cpp
#include <vector>
#include <cmath>
#include <assert.h>
#include "adapsimpsint.h"

using namespace std;

double AdapSimps::operator ()(double a,double b,double eps, int N2) const {

	if (a == b) 
		return 0;
		
	vector<struct info> infostack;
	struct info tmp;
	struct info tmp2;
	double result = 0;
	int cnt = 0;

	#ifdef DEBUG
		cout<<"a = "<<a<<", b = "<<b<<", intv = "<<b-a<<endl;
	#endif

	tmp.a = a;
	tmp.h = (b-a)/2;
	tmp.F[0] = f(a);
	tmp.F[1] = f(a + tmp.h);
	tmp.F[2] = f(b);
	tmp.S = tmp.h*(tmp.F[0] + 4*tmp.F[1] + tmp.F[2])/3;
	//	tmp.TOL = 10*eps;	
	double TOL = tmp.S*eps;  // here use eps as upperbound for relative error
	#ifdef DEBUG
		cout<<"tmp.S = "<<tmp.S<<"F[0:2]:"<<tmp.F[0]<<", "<<tmp.F[1]<<", "<<tmp.F[2]<<endl;
	#endif
	
//	infostack.push_back(tmp);
	
	int N1 = 3;
	int Nint = pow(2, N1-1);
	double intv = (b-a)/Nint;
	for (int i=0; i < Nint; i++) {	
		++cnt;
		tmp.a = a + i*intv;
		tmp.h = intv/2;
		tmp.F[0] = f(tmp.a);
		tmp.F[1] = f(tmp.a + tmp.h);
		tmp.F[2] = f(tmp.a + intv);
		tmp.S = tmp.h*(tmp.F[0] + 4*tmp.F[1] + tmp.F[2])/3;
		//	tmp.TOL = 10*eps;	
		tmp.TOL = TOL/Nint;  // here use eps as upperbound for relative error
		tmp.L = 1;
		#ifdef DEBUG
			cout<<"a + tmp.h = "<<a + (b-a)/2<<"tmp.a + intv = "<<tmp.a + intv<<endl;
			cout<<"i = "<<i<<", tmp.S = "<<tmp.S<<"F[0:2]:"<<tmp.F[0]<<", "<<tmp.F[1]<<", "<<tmp.F[2]<<endl;
		#endif

		infostack.push_back(tmp);		
	}


	while (!infostack.empty()) {
		++cnt;
		tmp = infostack.back();
		#ifdef DEBUG
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

		if (abs(S1 + S2 - tmp.S) < tmp.TOL) {
			result += (S1 + S2);
		#ifdef DEBUG
			cout<<"level = "<<tmp.L<<", S1, S2, and tmp.S: "<<S1<<" "<<S2<<" "<<tmp.S<<" Passed the test!"<<endl;		
		#endif			
		} else {
		#ifdef DEBUG
			cout<<"level = "<<tmp.L<<", S1, S2, and tmp.S: "<<S1<<" "<<S2<<" "<<tmp.S<<" Didn't pass the test...."<<endl;		
		#endif		
			if (tmp.L < N2) {
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
				result += (S1 + S2);
				#ifdef DEBUG
					cout<<"Level exceeded!"<<endl;
					cout<<"S1: "<<S1<<", S2: "<<S2<<", tmp.S: "<<tmp.S<<", result: "<<result<<", TOL: "<<tmp.TOL<<endl;
					cout<<"a = "<<tmp.a<<", b = "<<tmp.a + 2*tmp.h<<", intv = "<<2*tmp.h<<endl;
				#endif
//				return 0;
			}
		}
	}

	cout<<cnt<<"    ";
	return result;
}
//End of file Trapzint.cpp
