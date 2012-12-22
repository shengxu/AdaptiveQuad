#include <iostream.h>
#include <iomanip>
#include <vector.h>

#include "adapsimpsint.h"
#include "xsdata.h"

using namespace std;

int main() {
	Fun f;
	Trapz trapz1(f);
	cout<<"TRAPZ Int:"<<setprecision(7)<<trapz1(0,2,1e-7,20)<<endl;	
//ŒÆËã²¢Êä³ö»ý·Öœá¹û
}
