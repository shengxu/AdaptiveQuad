#include <iostream>
#include <iomanip>
#include <vector>

#include "adapsimpsint.h"

using namespace std;

int main() {
	Fun f;
	AdapSimps simpson1(f);
	cout<<"Adaptive Simpson Int:"<<setprecision(7)<<simpson1(0, 2, 1e-7, 20)<<endl;	
}
