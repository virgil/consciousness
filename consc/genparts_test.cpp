#include <assert.h>
#include <string>
#include <iostream>
#include <sstream>
#include <string>

#include "t_partition.h"
#include "PartitionEnumerator.h"

using namespace std;

long double BELL_NUMBERS[32] = { 
	0.0, 1.0, 2.0, 5.0, 15.0, 52.0, 203.0, 877.0, 4140.0, 21147.0, 115975.0, 678570.0, 4213597.0, 27644437.0, 190899322.0, 
	1382958545.0, 10480142147.0, 82864869804.0, 682076806159.0, 5832742205057.0, 51724158235372.0, 474869816156751.0, 4506715738447323.0, 
	44152005855084346.0, 445958869294805289.0, 4638590332229999353.0, 49631246523618756274.0, 545717047936059989389.0, 6160539404599934652455.0, 
	71339801938860275191172.0, 846749014511809332450147.0, 10293358946226376485095653.0 };


int main(void)
{
	unsigned int network_size = 5;
	PartitionEnumerator pgen( network_size );
	const unsigned int* result;

	// current partition
	unsigned int incP[network_size];

	t_partition P( network_size );
	
	//	vector<unsigned int> previous;
	
	int count = 0;
	
	while( 1 )
	{

		result = pgen.nextPartition();
		unsigned int numparts = result[0];
		memcpy(incP, result + 1, sizeof(unsigned int) * numparts);
		
//		vector<unsigned int> temp;
//		for( int i=0; i<numparts; i++ ) {
//			temp.push_back( (unsigned int) incP[i] );
//		}
		
//		if( previous == temp ) {
//			break;
//		}
		

		P.update( numparts, incP );
		count += 1;			
		
		cout << count << ". \t" << P.BAREstr() << endl;
		
//		cout << count << endl;

		
		if( P.size() == network_size )
			break;



	}
	
	if( count == BELL_NUMBERS[network_size] )
	{	
		cout << "Count is correct at " << count << endl;
	}
	cout << "count=" << count << endl;
	cout << "Bell[" << network_size << "] = " << BELL_NUMBERS[network_size] << endl;
	
	
	return 0;

}

