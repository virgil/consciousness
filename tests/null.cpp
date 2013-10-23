#include <iostream>
#include <assert.h>
#include <float.h>
#include <limits.h>

using namespace std;

int main()
{

//	double x = NULL;
	double x = numeric_limits<double>::quiet_NaN();
	double y = numeric_limits<double>::quiet_NaN();;
	
	cout << "x=" << x << "\t y=" << y << endl;
	cout << "------------------------------------" << endl;
	cout << "x+y=" << x+y << endl;
	cout << "x-y=" << x-y << endl;
	cout << "x*y=" << x*y << endl;
	cout << "x/y=" << x/y << endl;

	cout << "y-x=" << y-x << endl;	
	cout << "y/x=" << y/x << endl;	

	
	cout << "x<y  \t=" << (x<y) << endl;
	cout << "x==y \t=" << (x==y) << endl;	
	cout << "x>y  \t=" << (x>y) << endl;
	cout << "x>=0  \t=" << (x>=0) << endl;
	
	if( x == 0.0 )
		cout << "x == 0.0" << endl;
	else
		cout << "x != 0.0" << endl;


	return 0;
}
