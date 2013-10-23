// this program checks that newly allocated memory is filled with zeros.

#include <iostream>
#include <assert.h>

const unsigned int MAX_SIZE = 100;


using namespace std;

int main()
{

	unsigned int* blah = new unsigned int[MAX_SIZE];
	double* blah2 = new double[MAX_SIZE];
	
	
	for( int i=0; i<MAX_SIZE; i++ ) {
		assert( blah[i] == 0 );
		cout << i << " was 0!" << endl;
		
		assert( blah2[i] == 0.0 );
		cout << i << " was 0.0!" << endl;
        
        assert( blah[i] == 0 );
        assert( blah2[i] == 0.0 );
	}
	
    cout << "All allocated entries were zero." << endl;
    
	return 0;
}