//////////////////////////////////////////////////////////////////////////////////////////
// gen_allpartitions.cpp -- outputs partitions in various formats
//////////////////////////////////////////////////////////////////////////////////////////
// Original recursive Algorithm by Arend Hintze, modified and optimized by Virgil Griffith
//////////////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <math.h>			// for ceil()
#include <assert.h>
#include "consc/helpers.h"
#include "consc/t_partition.h"
#include "consc/PartitionEnumerator.h"
#include "consc/t_consc.h"
//#include <string>

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////
// Only one of these can be on.
// FLAG__PRINT_IN_BINARY_FORMAT is fastest
//////////////////////////////////////////////////////////////////////////////////////////

// can also be "MASK", "BARE", or "BOTH"
const string FLAG__PRINT_IN = "BARE";
bool PRINT_PARTITIONS = FALSE;
const bool SKIP_TOTAL_PARTITION=TRUE;

//////////////////////////////////////////////////////////////////////////////////////////
const bool FLAG__CHECK_COUNTS = TRUE;
const bool FLAG__PRINT_NUMBERS = TRUE;
//////////////////////////////////////////////////////////////////////////////////////////
// Global variables that the user sets
//////////////////////////////////////////////////////////////////////////////////////////
unsigned int N=0;					// Required argument
unsigned int MAX_NUM_PARTS=0;		// optional. default: MAX_NUM_PARTS = N
unsigned int SUBSET=0;				// optional. default: no subset == 2**(N-1) -- meaning all N nodes are 1.
//////////////////////////////////////////////////////////////////////////////////////////

// we use long double because long long wasn't big enough
long double BELL_NUMBERS[32] = { 
0.0, 1.0, 2.0, 5.0, 15.0, 52.0, 203.0, 877.0, 4140.0, 21147.0, 115975.0, 678570.0, 4213597.0, 27644437.0, 190899322.0, 
1382958545.0, 10480142147.0, 82864869804.0, 682076806159.0, 5832742205057.0, 51724158235372.0, 474869816156751.0, 4506715738447323.0, 
44152005855084346.0, 445958869294805289.0, 4638590332229999353.0, 49631246523618756274.0, 545717047936059989389.0, 6160539404599934652455.0, 
71339801938860275191172.0, 846749014511809332450147.0, 10293358946226376485095653.0 };

// contains the number of partitions printed of size >= index
//	vector<long double> NUMBER_PRINTED;
vector<unsigned long> NUMBER_PRINTED;


//////////////////////////////////////////////////////////////////////////////////////////
// Function prototypes
//////////////////////////////////////////////////////////////////////////////////////////
void print_partition( t_partition* P );
void print_usage();
//string binary_repr( unsigned int number, int padding=1 );
//////////////////////////////////////////////////////////////////////////////////////////

void print_usage()
// print how to use
{
	cout << endl << "USAGE: \n ./gen_partitions N [K] [S]" << endl;
	cout << "\tN - the number of nodes. N >= 1" << endl;
	cout << "\tK - the maximum number of parts to divide the network into. 1<=K<=N.  Default=N" << endl;
	cout << "\tFor just bipartions, K=2.  For bi and tri partitions, K=3.  For all partitions (default), K=N." << endl;	
	cout << "\tS - N-digit binary string representing the nodes to compute the partitions within. If specified, N is ignored." << endl;
//	cout << "\tP - Print the partitions?  {0,1}" << endl;
	
	exit(-1);
}

string new_binary_repr( const unsigned int digit, unsigned int numplaces )
{
	// if numplaces isn't defined, assume whole network
	if( numplaces == 0 )
        numplaces = N;
	
	string z="";
	for( int i=0; i<numplaces; i++ )
	{
		int bit = ((digit >> i)&1);
		
		if( bit )
			z.insert(0, "1" );
		else
			z.insert(0, "0" );
	}
	
	return z;
}


bool get_next_partition( PartitionEnumerator* pgen, t_partition* P, unsigned int* nodemap, unsigned int numnodes )
{
	
	const unsigned int* result = pgen->nextPartition();
	unsigned int numparts = result[0];
	
	// if numparts==0, then we've exhausted all of the partitions
	if( ! numparts ) {
		return false;
	}

//	if( numparts == numnodes ) {
//		P->update( numparts, result+1 );
//		return true;
//	}

	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// remap result -> incoming_P
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int* incoming_P = (unsigned int*) calloc(numparts, sizeof(unsigned int));
	remapNodes( incoming_P, result, nodemap, numnodes );
	
	// copy incoming_P into P; free incoming_P
	P->update( numparts, incoming_P );
	free( incoming_P );
	
	return true;	
}

unsigned int numbits_on( unsigned int subset )
{
	// clear the old part
	unsigned int z = 0;
	
	do	{
		if( (subset&1) )
			z += 1;
		subset >>= 1;
		
	}while( subset );
	
	return z;
}


int main(int argc, char* argv[])
{	
	// only one of these can be on.
	assert( FLAG__PRINT_IN == "BARE" || FLAG__PRINT_IN == "MASK" || FLAG__PRINT_IN == "BOTH" );

    string subset_string = "";
    
	switch( argc )
	{
		//specified N, MAX_NUM_PARTS, SUBSET, and PRINT_PARTITIONS
		case 5:
			PRINT_PARTITIONS = (bool) atoi(argv[4]);

		//specified N, MAX_NUM_PARTS, and SUBSET	
		case 4:
            subset_string = string(argv[3]);

		//specified N and MAX_NUM_PARTS			
		case 3:
			MAX_NUM_PARTS = atoi(argv[2]);

		//specified only N		
		case 2:
			N = atoi(argv[1]);
			break;

		default:
			print_usage();
	}


	///////////////////////////////////////////////////////////////////////////////////////////
	// set defaults for SUBSET and MAX_NUM_PARTS
	///////////////////////////////////////////////////////////////////////////////////////////
	if( SUBSET == 0 ) 
		SUBSET = (1 << N)-1;

    if( subset_string != "" && subset_string.length() != N ) {
        cerr << "Length of binary string didn't match N.  binary_string='" << subset_string << "'" << endl;
        print_usage();
    }
    
    // if the subset string was defined...
    if( subset_string.length() == N ) {
        SUBSET = 0;
        for( int i=0; i<subset_string.length(); i++ ) {
            
            // if this bit isn't zero, make it one.
            if( subset_string[i] != '0' )
                SUBSET |= 1;
            
            // shift subset over by one.
            SUBSET <<= 1;
        }

        // now shift it back so the last digit is in the 1's place.
        SUBSET >>= 1;
        
//        cout << "subset_bitmask=" << subset_string << endl;
    }

	///////////////////////////////////////////////////////////////////////////////////////////
	// set subset_string to binary_repr of SUBSET
	///////////////////////////////////////////////////////////////////////////////////////////
    subset_string = new_binary_repr( SUBSET, N );

    
	if( numbits_on(SUBSET) < N )
		N = numbits_on(SUBSET);
	
	if( MAX_NUM_PARTS == 0  || MAX_NUM_PARTS > N )
		MAX_NUM_PARTS = N;

    
	
	///////////////////////////////////////////////////////////////////////////////////////////
	// sanity checks
	///////////////////////////////////////////////////////////////////////////////////////////
	assert( 1 <= N );
	assert( N <= 32 );
		
	assert( 1 <= SUBSET );
	
//	unsigned int highest = (1<<N)-1;
//	cout << "N=" << N << endl;
//	cout << "highest=" << highest << endl;
//	cout << "SUBSET=" << SUBSET << endl;
	
	assert( numbits_on(SUBSET) <= N );
	
	assert( 1 <= MAX_NUM_PARTS );
	assert( MAX_NUM_PARTS <= N );
	

	// this is MAX_NUM_PARTS+1 because the vector is zero-indexed
	if( FLAG__CHECK_COUNTS )
		NUMBER_PRINTED.resize( MAX_NUM_PARTS+1 );

	///////////////////////////////////////////////////////////////////////////////////////////
	

	///////////////////////////////////////////////////////////////////////////////////////////
	// Initialize the partition enumerator, the partition, and the nodemap
	///////////////////////////////////////////////////////////////////////////////////////////
    cerr << "N=" << N << endl;
	PartitionEnumerator pgen( N, SKIP_TOTAL_PARTITION );
//	PartitionEnumerator pgen( 12, SKIP_TOTAL_PARTITION );    

	
//	// define a partition with the maximum number of parts
//	t_partition P( MAX_NUM_PARTS, SUBSET );
	t_partition P( N, SUBSET );
	
	// initialize nodemap to be as large as it could ever be, numnodes
	// Now to create the array that maps from result -> incoming partition
	unsigned int nodemap[N];
	memset( nodemap, 0, sizeof(unsigned int) * N );
	subset2array( SUBSET, nodemap );
	
//	for( int i=0; i<N; i++ ) {
//		cout << "nodemap[" << i << "] = " << nodemap[i] << endl;
//	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// END partitions generator
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// now read the parts from GENERATOR...
	while( get_next_partition(&pgen, &P, nodemap, N) )
	{
		if( P.size() > MAX_NUM_PARTS )
			continue;
		
		if( PRINT_PARTITIONS )
			print_partition( &P );

		// increment all NUMBER_PRINTED >= P.size();
		if( FLAG__CHECK_COUNTS )
		{
			for( int i=P.size(); i<=MAX_NUM_PARTS; i++ ) {
				NUMBER_PRINTED[i] += 1;
			}
		}
		
	}
		
	
	if( FLAG__CHECK_COUNTS )
	{
		
		// if we skipped the total partition, add one to each
		cout << "- nodes=" << N << " \t max_parts=" << MAX_NUM_PARTS << " \t #partitions=" << NUMBER_PRINTED.back() << endl;		
		cout << "- subset=" << subset_string << "\t mask=" << SUBSET << endl;		

		// sanity check for a few max_num_parts sizes
		assert( NUMBER_PRINTED[0] == 0 );
		
		unsigned long number_expected = 0;
		
		if( MAX_NUM_PARTS >= 1 ) {								// K = 1
			number_expected = 1;
			if( SKIP_TOTAL_PARTITION ) number_expected--;

			assert( NUMBER_PRINTED[1] == number_expected );
			cout << "+ Passed K =1 check.";
			cout << "\t\t #partitions=" << NUMBER_PRINTED[1] << endl;
		}	
		if( MAX_NUM_PARTS >= 2 ) {								// K = 2
			number_expected = (1 << (N-1) );
			if( SKIP_TOTAL_PARTITION ) number_expected--;
			assert( NUMBER_PRINTED[2] == number_expected );
			cout << "+ Passed K<=2 check.";
			cout << "\t\t #partitions=" << NUMBER_PRINTED[2] << endl;		
		}
		if( MAX_NUM_PARTS >= 3) {								// K = 3
			number_expected = ceil(pow(3,N-1) / 2.0);
			if( SKIP_TOTAL_PARTITION ) number_expected--;

			assert( NUMBER_PRINTED[3] == number_expected );
			cout << "+ Passed K<=3 check.";
			cout << "\t\t #partitions=" << NUMBER_PRINTED[3] << endl;		
		}
		if( MAX_NUM_PARTS >= 6 && MAX_NUM_PARTS >= (N-2) ) {	// K=N-2	
			// K = N-2
			number_expected = BELL_NUMBERS[N] - 1 -(N*(N-1)/2.0);
			if( SKIP_TOTAL_PARTITION ) number_expected--;

			assert( NUMBER_PRINTED[N-2] == number_expected );
			cout << "+ Passed K<=" << N-2 << " check.";
			cout << "\t\t #partitions=" << NUMBER_PRINTED[N-2] << endl;
		}
		if( MAX_NUM_PARTS >= 5 && MAX_NUM_PARTS >= (N-1) ) {		// K = N-1
			number_expected = BELL_NUMBERS[N]-1;
			if( SKIP_TOTAL_PARTITION ) number_expected--;

			assert( NUMBER_PRINTED[N-1] == number_expected );

			cout << "+ Passed K<=" << N-1 << " check.";
			cout << "\t\t #partitions=" << NUMBER_PRINTED[N-1] << endl;		
		}
		if( MAX_NUM_PARTS >= 4 && MAX_NUM_PARTS == N )				// K = N
		{
			number_expected = BELL_NUMBERS[N];
			if( SKIP_TOTAL_PARTITION ) number_expected--;	

			assert( NUMBER_PRINTED[N] == number_expected );
			cout << "+ Passed K<=" << N << " check.";
			cout << "\t\t #partitions=" << NUMBER_PRINTED[N] << endl;		
		}

		assert( NUMBER_PRINTED.back() <= BELL_NUMBERS[N] );
	
	}

	return 0;
}

void print_partition( t_partition* P )
{
	static unsigned long count = 1;
	
	if( FLAG__PRINT_NUMBERS )
		cout << count << ".\t";
	
	if( FLAG__PRINT_IN == "MASK" || FLAG__PRINT_IN == "BOTH" ) {
		cout << P->MASKstr();
	}
	
	if( FLAG__PRINT_IN == "BOTH" )
		cout << "\t\t";

	if( FLAG__PRINT_IN == "BARE" || FLAG__PRINT_IN == "BOTH" ) {
		cout << P->BAREstr();
	}

	cout << " \t K=" << P->size() << endl;
	
//	P->sort();
//	cout << "\t\t P: " << P->BAREstr() << endl;
	
	count += 1;
}

