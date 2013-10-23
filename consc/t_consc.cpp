#ifndef T_CONSC_CPP
#define T_CONSC_CPP

///////////////////////////////////////////////////////////////////////////////////////////
// Main file for consciousness
// by Virgil Griffith, originally derived from algorithms by Arend Hintze
//////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <algorithm>
#include <boost/foreach.hpp>
//#include <boost/array.hpp>    // for a better fixed-length array that supports STL functions like .size(), etc.
#include <limits.h>		// for min, max values of integers, doubles, floats, etc.

// for openMP
//#include <omp.h>
//#pragma omp

//#include <boost/unordered_map.hpp>
//#include <bitset>

#include "t_state.h"

//#define NDEBUG
#include <assert.h>
#include "helpers.h"

using namespace std;



////////////////////////////////////////////////////////////////
// Configuration Options
////////////////////////////////////////////////////////////////
const unsigned short OUTPUT_PRECISION = 4;

////////////////////////////////////////////////////////////////
//const string FLAG__MIP_METHOD = "ATOMIC";
const string FLAG__MIP_METHOD = "TONONI";
//const string FLAG__MIP_METHOD = "TONONI_PROD_M0";
//const string FLAG__MIP_METHOD = "TONONI_PROD_M1";
//const string FLAG__MIP_METHOD = "NONE";
////////////////////////////////////////////////////////////////

const string FLAG__STATE_EI = "UNITS";
//const string FLAG__STATE_EI = "WIRES";   // <ei(x1/P)> using the ``wire perturbations'' from Balduzzi.

#define DO_CAPITALIZED_ASSERTS 1

const bool FLAG__CACHE_M1=true;
const bool FLAG__CACHE_M0_GIVEN_M1=true;


// output a line for the atomic partition no matter what?
const bool SAVE_ATOMIC_PARTITION=true;

// Don't get any more MIPs than this.
const unsigned int MAXIMUM_NUMBER_OF_MIPS_TO_STORE = 100;
const unsigned int MAXIMUM_NUMBER_OF_MAIN_COMPLEXES_TO_STORE = 20;

// What format to print the partitions in?  BARE or JSON
const string FLAG__PRINT_PARTITIONS_FORMAT = "BARE";

// no trailing slash here
string DEFAULT_TEMP_DIRECTORY="/tmp";

const unsigned int DECIMAL_PRECISION=5;
const long double PRECISION = 1.0/pow(10.0, (double) DECIMAL_PRECISION+1);




////////////////////////////////////////////////////////////////
#if DO_CAPITALIZED_ASSERTS
	#define ASSERT(x) assert(x)		//defined as plain-ol assert()
#else
	#define ASSERT(x)				//defined as nothing
#endif
////////////////////////////////////////////////////////////////

#ifndef foreach
	#define foreach BOOST_FOREACH
#endif

#ifndef reverse_foreach
	#define reverse_foreach BOOST_REVERSE_FOREACH
#endif

#define restrict __restrict__

#include "t_consc.h"
#include "t_partition.h"
#include "PartitionEnumerator.h"
#include "helpers.h"
#include "t_subset.h"
#include "t_consc_part2.cpp"


bool check_use_total_partition();

PartitionEnumerator PARTITIONS_GENERATOR( 1, (check_use_total_partition()==false) );


/////////////////////////////////////////////////////////////////
// Do we allow the total partition? 
/////////////////////////////////////////////////////////////////
bool check_use_total_partition()
// returns TRUE if FLAG__MIP_METHOD disallows the total partition
{
	const char *MIP_methods_that_allow_total_partition[] = { "NONE" };
	
	foreach( const char* str, MIP_methods_that_allow_total_partition )
		if( strcmp(FLAG__MIP_METHOD.c_str(),str) == 0 )
		   return true;
		   
   return false;
	
}
/////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
//sanity check constructor
////////////////////////////////////////////////////////////////
t_consciousness::t_consciousness()
{
	

	// ASSERT()s are ON or OFF ?
#if DO_CAPITALIZED_ASSERTS
	cerr << "- CAPITALIZED ASSERTS ARE ON" << endl;
#endif
	
	if( FLAG__CACHE_M1==false && FLAG__CACHE_M0_GIVEN_M1==false )
		cerr << "- CACHEs: ** NONE **" << endl;
	else
	{
		cerr << "- CACHEs: ";
		
		if( FLAG__CACHE_M1 )
			cerr << "H[M1]; ";
		if( FLAG__CACHE_M0_GIVEN_M1 )
			cerr << "H[M0|M1]; ";
		cerr << "\b\b " << endl;
//		cerr << endl;
	}
	
	assert( (FLAG__PRINT_PARTITIONS_FORMAT == "BARE") ^ (FLAG__PRINT_PARTITIONS_FORMAT == "JSON") );
	
	assert( 1 <= DECIMAL_PRECISION );
	
	// if DEFAULT_TEMP_DIRECTORY ends with a '/', remove it.
	if( DEFAULT_TEMP_DIRECTORY[DEFAULT_TEMP_DIRECTORY.length()-1] == '/' ) {
		DEFAULT_TEMP_DIRECTORY.erase( DEFAULT_TEMP_DIRECTORY.length()-1 );
	}


	////////////////////////////////////////////////////////////////////////
	// set our fine default values
	////////////////////////////////////////////////////////////////////////
	this->neuron_fires_only_at_threshold = false;
	this->states = NULL;
	this->max_num_parts = 0;
	this->FULL_MASK = 0;
	this->numunits = 0;
	this->numstates = 0;	
	this->H_X1 = DBL_NONE;	
	this->H_X0_GIVEN_X1 = DBL_NONE;
	this->NUM_X1_STATES = DBL_NONE;
	this->H_M1_cache = NULL;
	this->H_M0_GIVEN_M1_cache = NULL;
	this->prob_s1_given_mu0__Vnodes = NULL;
	
	this->use_total_partition = check_use_total_partition();

	cerr << "- MIP_METHOD=" << FLAG__MIP_METHOD;
	cerr << "\t Use Total Partition: " << use_total_partition << endl;
	
		
//	For testing whether CAPITALIZED ASSERTs are ON
//	ASSERT( 0 == 1 );
}

// destructor
t_consciousness::~t_consciousness()
{
	if( x1_states != NULL ) {
		free( x1_states );
		x1_states = NULL;
	}
	
	if( num_x0s_to_this_x1 != NULL ) {
		free( num_x0s_to_this_x1 );
		num_x0s_to_this_x1 = NULL;
	}
	
	if( H_M1_cache != NULL ) {
		delete [] H_M1_cache;
		H_M1_cache = NULL;
	}

	if( H_M0_GIVEN_M1_cache != NULL ) {
		delete [] H_M0_GIVEN_M1_cache;
		H_M0_GIVEN_M1_cache = NULL;
	}
	
	
	// this is the cache for prob_s1_given_mu0
	if( prob_s1_given_mu0__Vnodes != NULL ) {
		delete [] prob_s1_given_mu0__Vnodes;
		prob_s1_given_mu0__Vnodes = NULL;
	}
}


double t_consciousness::H_M0( const t_subset& restrict S )
// returns H[M0] for a given subset M
{
	double z = S.numnodes();
	return z;
}

double t_consciousness::H_M0( const unsigned int S )
// returns H[M0] for a given subset M
{
	return (double) numunits_in_mask(S);
}

// entropy of a partition at X0
// TODO: Probably want to remove this H_M0_normalized function
double t_consciousness::H_M0_normalized( const t_partition& restrict P )
{
	const double summ = numunits;
	double z = 0.0;

	assert( 2 <= P.size() );
	
	for( int i=0; i<P.size(); i++ )
	{
		const double this_prob = P.size(i) / summ;
		ASSERT( this_prob > 0.0 );
		z -= this_prob * log2( this_prob );
	}
	
	// With K symbols the maximum value of the entropy is log2( K )	
	// In this case K = #symbols = #parts = P.size()
//	z /= log2( P.size() );
	
	assert( 0.0 < z );
	return z;
}

// entropy of a partition at X1
// TODO: Probably want to remove this H_M1_normalized function
double t_consciousness::H_M1_normalized( const t_partition& restrict P, bool environment )
// NOTE: THIS CAN RETURN ZERO (and be correct)!!!!!!
{
	double summ=0.0, z=0.0;

	// make an array to store the entropy of each part
	// define a block of memory to hold the number of times a mu1 occurs
	double* restrict M1_entropies = (double*) calloc(P.size(), sizeof(double));
	
	// calculate the summ of the entropies
	for( int i=0; i<P.size(); i++ ) {
			M1_entropies[i] = H_M1( P[i] );
		summ += M1_entropies[i];
	}
	
//	cerr << "entropies: ";
	// calculate the entropy of the summ of the M1_entropies
	for( int i=0; i<P.size(); i++ ) {

//		cerr << M1_entropies[i] << " ";
		if( fequals(M1_entropies[i],0.0) )
			continue;

		const double this_prob = M1_entropies[i] / summ;
		
		z -= this_prob * log2( this_prob );
	}
	
//	cerr << endl;
	
	// With K symbols the maximum value of the entropy is log2( K )
	// In this case K = #symbols = #parts = P.size()
	assert( 2 <= P.size() );
//	z /= log2( P.size() );

	// free the array
	free( M1_entropies );

	assert( 0.0 < z );
	return z;
}

double t_consciousness::H_S0_GIVEN_M1( const unsigned int S0mask, const unsigned int M1mask )
// returns H[S0|M1] via some algebra
// H[S0|M1] = H(S0) + H(M1|S0) - H(M1)
{
	return H_A0_B1( S0mask, M1mask ) - H_M1( M1mask );
	
//	return H_M0(S0mask) + H_M1_GIVEN_S0(M1mask, S0mask) - H_M1(M1mask);
}

double t_consciousness::H_M1_GIVEN_M0( const unsigned int partmask )
// returns H[M1|M0] via some algebra
// H[M1|M0] = H(M1) + H(M0|M1) - H(M0)
{
	return H_M1( partmask ) + H_M0_GIVEN_M1__ELEMENTS( partmask ) - H_M0( partmask );
}



t_psi_result t_consciousness::psi( const t_state x1 )
// computes the psi measure for a given state x1
{
    assert( is_valid_mu1( x1 ) );
    assert( x1.size() == this->numunits );
    
    t_psi_result z;
    z.ei = this->ei( x1 );
    z.x1 = x1;
//    z.lowerbound = psi_lowerbound( x1 );
//    z.upperbound = psi_upperbound( x1 );
    
    if( fequals(z.lowerbound,z.upperbound) )
        z.lowerbound = min(z.lowerbound,z.upperbound);
    
    assert( z.lowerbound <= z.upperbound );
    
    return z;
}





bool t_consciousness::break_MIPscore_ties( const t_partition& restrict P0, const t_partition& restrict P1 )
// returns TRUE if P1 > P0
{
	
	// If P1 has MORE PARTS than P0, use P1. 
	return ( P1.size() > P0.size() );
	
}


double t_consciousness::prob_mu0_given_mu1__anticond( const unsigned int mu0, const unsigned int mu1, const unsigned int partmask, const unsigned int partsize )
// returns p(mu0|mu1 \anticond \widetilde{M} )
// this is defined in terms of Bayes rule.
{
	assert( (mu0|partmask) == partmask );	
	assert( is_valid_mu1(mu1, partmask) );
	
	const double top = prob_mu1_given_mu0__anticond( mu1, mu0, partmask, partsize );
	
	if( top == 0.0 )
		return 0.0;
	
	double bottom=0.0;
	for( uint x0=0; x0<=partmask; x0++ )
	{
		// if this x0 has bits on outside of the part, skip it.
		if( (x0 | partmask) != partmask )
			continue;

		// remove all bits outside of the partmask.
		const unsigned int m0 = x0 & partmask;
		
		bottom += prob_mu1_given_mu0__anticond( mu1, m0, partmask, partsize );
	}

//	cout << "top=" << top << "\tbottom=" << bottom << endl;
	
	double z = top / bottom;
	
	assert( 0.0 <= z );
	assert( z <= 1.0 );
	return z;
}


/////////////////////////////////////////////////////////////////////////////
// Functions for computing INFORMATION FLOW via Ay, Polani (2006)
/////////////////////////////////////////////////////////////////////////////
double t_consciousness::I_A0_ARROW_B1( const unsigned int Amask, const unsigned int Bmask )
// computes I( A0 -> B1 ) without imposing (preventing flow through) any subset S0
// I(A0 -> B1) = \sum_{a0} p(a0) \sum_{b1} p(b1 | \hat{a0}) * log2( \frac{p(b1|\hat{a0})} {\sum_{a`0} p(a`0) * p(b1|\hat{a`0} }  )
// = p(a0) \sum_{a0} \sum_{b1} p(b1 | \hat{a0}) * log2( \frac{p(b1|\hat{a0})} {p(a`0) * \sum_{a`0} p(b1|\hat{a`0} } )
// = 2^{-|A|} * \sum_{a0} \sum_{b1} p(b1 | \hat{a0}) * log2( \frac{p(b1|\hat{a0})} {2^{-|A|} * \sum_{a`0} p(b1|\hat{a`0} } )
// ---------------------------------------------------------------------------
// this can be misleading because it IGNORES info-flow due solely to INTERACTIONS with other subset(s).
// For example: If B = A xor C
// then I(A0 -> B1)=0 even though I(A0 -> B1|\hat{C0})=1
// ---------------------------------------------------------------------------
// For multiple timesteps it gets I(A0 -> B1) gets worse because it can be positive when A has no *direct* influence on B.
// As such, it's unclear how useful this function is.  It may need to be changed.
// Typically, this function will only be called by I_A0_ARROW_B1_IMPOSING_S0()
{
	// assert Amask is valid
	assert( Amask && Amask <= FULL_MASK );

	// assert Bmask is valid
	assert( Bmask && Bmask <= FULL_MASK );
	
	// setting Smask=s0=0 computes I(A0 -> B1)
	double z = I_A0_ARROW_B1_IMPOSING_s0( Amask, Bmask, 0, 0 );

	// remove any -0.0s
	if( z == -0.0 )
		z = 0.0;

	// check bounds
	assert( 0.0 <= z );
	
	return z;
}

double t_consciousness::I_A0_ARROW_B1_IMPOSING_S0( const unsigned int Amask, const unsigned int Bmask, const unsigned int Smask )
// computes I( A0 -> B1 | \hat{S0} )
// = \sum_{s0} p(s0) I( A0 -> B1 | \hat{s0} )
// = p(s0) \sum_{s0} I( A0 -> B1 | \hat{s0} )
// = 1.0/2^{|S|} \sum_{s0} I( A0 -> B1 | \hat{s0} )
// ----------------------------------------------------
// imposing a variable S0 can *INCREASE* or DECREASE I(A0 -> B1)
// Note that Smask is ALLOWED TO BE BLANK
{	
	// if Smask is BLANK, use I(A0 -> B1)
	if( Smask == 0 )
		return I_A0_ARROW_B1( Amask, Bmask );
	
	
	// assert Amask and Smask are valid
	assert( Amask && Amask <= FULL_MASK );
	assert( Smask && Smask <= FULL_MASK );
	
	
//	cout << "A & S = " << binary_repr( (Amask & Smask) ) << endl;
	// Amask and Smask should not overlap
	assert( (Amask & Smask) == 0 );

	// assert Bmask is valid
	assert( Bmask && Bmask <= FULL_MASK );
	
	// p(s0) = 1.0 / 2^{|S|}
	const double prob_s0 = 1.0 / (1 << numunits_in_mask(Smask));
								  
	double summ=0.0;
	
	//foreach s0 state...
	for( int x0=0; x0<=Smask; x0++ )
	{
		// if this x0 has anything ON outside of Smask, skip it
		if( (x0 | Smask) > Smask )
			continue;
		
		const unsigned int s0 = x0;
		
		summ += I_A0_ARROW_B1_IMPOSING_s0( Amask, Bmask, Smask, s0 );		
	}
	
	double z = prob_s0 * summ;
	
	// remove any -0.0's
	if( z == -0.0 )
		z = 0.0;
	
	// assert lower bound
	assert( 0.0 <= z );
	
	return z;
}

double t_consciousness::I_A0_ARROW_B1_IMPOSING_s0( const unsigned int Amask, const unsigned int Bmask, const unsigned int Smask, const unsigned int s0 )
// computes I( A0 -> B1 | \hat s0 )
// should only be called by I_A0_ARROW_B1_IMPOSING_S0() or I_A0_ARROW_B1()
// Note that Smask and s0 ARE ALLOWED TO BE ZERO (be blank)
{
	// assert Amask is valid
	assert( Amask && Amask <= FULL_MASK );
	
	// assert Smask is valid
	assert( Smask <= FULL_MASK );
	
	// assert s0 is subset-or-equal to Smask
	assert( (s0 | Smask) == Smask );
	
	// TODO: This should be made more efficient by using the old Amask_size
	// p(a0) = 1.0 / 2^{|A|}
	const double prob_a0 = 1.0 / (1 << numunits_in_mask(Amask));
	const double prob_a0prime = prob_a0;
	
	const unsigned int ASmask = Amask | Smask;
    // assert ASmask is valid
	assert( ASmask && ASmask <= FULL_MASK );

	// this is the big loop for \sum_{a0} \sum_{b1}
	double summ_over_a0b1=0.0;
	
	//foreach a0 state...
	for( int possible_a0=0; possible_a0<=Amask; possible_a0++ )
	{
		// if there is anything ON outside of Amask, skip it
		if( (possible_a0 | Amask) > Amask )
			continue;
		
		const unsigned int a0 = possible_a0;
		const unsigned int a0s0 = a0 | s0;
		
		//foreach b1 state...
		for( int possible_b1=0; possible_b1<=Bmask; possible_b1++ )
		{
			// if counter has anything ON outside of Bmask, skip it.
			if( (possible_b1|Bmask) > Bmask )
				continue;

            if( ! is_valid_mu1(possible_b1,Bmask) )
                continue;			
            
			const unsigned int b1 = possible_b1;
            
			const double prob_b1_imposing_a0s0 = prob_b1_IMPOSING_a0s0( Bmask, b1, ASmask, a0s0 );
			
			// avoid computing the logterm is p(b1|\hat{a0},\hat{s0})=0.0
			if( prob_b1_imposing_a0s0 == 0.0 )
				continue;

			
			// Now compute the summ over \sum_{a0`} p(b1| \hat{a`0} \hat{s0} )			
			double summ_over_a0prime=0.0;
			for( int possible_a0prime=0; possible_a0prime<=Amask; possible_a0prime++ ) {

				// if there is anything ON outside of Amask, skip it
				if( (possible_a0prime | Amask) > Amask )
					continue;

				const unsigned int s0a0prime = s0 | possible_a0prime;
				
				summ_over_a0prime += prob_b1_IMPOSING_a0s0( Bmask, b1, ASmask, s0a0prime );
			}
			
			const double logterm = prob_b1_imposing_a0s0 / (prob_a0prime * summ_over_a0prime);
			
			// add this to the summ over b1
			summ_over_a0b1 += prob_b1_imposing_a0s0 * log2( logterm );
		}		
	}

	const double z = prob_a0 * summ_over_a0b1;
	
	assert( 0.0 <= z );
		
	return z;
}

double t_consciousness::prob_b1_IMPOSING_a0s0( const unsigned int Bmask, const unsigned int b1, const unsigned int ASmask, const unsigned int a0s0 )
// should only be called by I_A0_ARROW_B1_IMPOSING_s0
// computes p( b1 | \hat{a0}, \hat{s0} )
{
	
	// assert b1 is valid
	assert( Bmask && Bmask <= FULL_MASK );
	assert( is_valid_mu1(b1,Bmask) );

	
	// assert a0s0 is valid
	assert( (a0s0 | ASmask) == ASmask );
	assert( ASmask && ASmask <= FULL_MASK );
	
	// todo: save this result
	const unsigned int Bsize = numunits_in_mask( Bmask );
	const unsigned int ASsize = numunits_in_mask( ASmask );
	
	// if ASmask is the FULL_MASK, then we can skip this multiplication step.
	if( ASmask == FULL_MASK )
	{
		const unsigned int x1 = states[a0s0];

		// if x1 matches b1: return 1.0. else return 0.0;
		return ((x1 & Bmask) == b1) ? 1.0 : 0.0;
	}

	
	double z=1.0;
	unsigned int* restrict Bnodes = new unsigned int[Bsize];
	subset2array(Bmask, Bnodes);
	
	for( int i=0; i<Bsize; i++ ) {
		const uint n_index = Bnodes[i];
		// only keep the bit at location n_index
		const uint n1 = (b1 & (1<<n_index));
		
		z *= prob_n1_GIVEN_a0s0( n1, n_index, a0s0, ASmask, ASsize );
		
		// if the probability becomes 0.0, break out of the loop.  We must BREAK
		// instead of return 0.0 because we still have to delete [] the nodes.
		if( z == 0.0 )
			break;
	}
	
	delete [] Bnodes;
	
	// check bounds
	assert( 0.0 <= z );
	assert( z <= 1.0 );
	
	return z;
}

double t_consciousness::prob_n1_GIVEN_a0s0( const unsigned int n1, const unsigned int n_index, const unsigned int a0s0, 
										    const unsigned int ASmask, const unsigned int ASsize )
// should only be called by prob_b1_IMPOSING_a0s0
// computes p( n1 | a0, s0 ) via summing over everything outside of a0,s0
{
	//check n1 is valid (is zero or has n_index set to 1
	assert( n1==0 || n1 == (1<<n_index) );
	
	//check n_index is valid
	assert( n_index < numunits );

	//assert ASmask is valid (exists and is <= FULL_MASK)
	assert( ASmask && ASmask <= FULL_MASK );
	
	//assert a0s0 is a subset of ASmask
	assertmsg( (a0s0 | ASmask) == ASmask, "Error: s0 has a entry ON outside of S0mask." );
	
	//check ASsize is valid
	assert( ASsize );
//	cout << "ASsize=" << ASsize << "\t\t numunits=" << numunits << endl;
//	cout.flush();
	assert( ASsize <= numunits );
	
	
	// Z is defined as EVERYTHING OUTSIDE OF ASmask
	const unsigned int Zmask = (FULL_MASK & (~ASmask) );
	const unsigned int Zsize = numunits_in_mask(Zmask);
	
	assert( Zmask );
	
	// anticond_mask is everything outside of S0mask
	//const unsigned int anticond_mask = (FULL_MASK & (~S0mask));
	unsigned int num_matches=0;
	
	//foreach x0 state...
	for( int x0=0; x0<=FULL_MASK; x0++ )
	{
		// skip all x0 states that dont match a0s0.
		if( (x0 & ASmask) != a0s0 )
			continue;
		
		const unsigned int x1 = states[x0];
		
		// for this x1 only keep the state of the node at location n_index
		//const unsigned int temp_n1 = (x1 >> n_index)&1;		// OLD BAD WAY
		const unsigned int temp_n1 = (x1 & (1<<n_index));
		
		num_matches += int(temp_n1 == n1);
	}
	
	
	// divide by the number of states we tried -- 2**(SIZE_OF_ANTICOND_MASK)
	const double states_tried = (double) (1 << Zsize);
	
	// p(n1|a0s0) = num_matches / states_tried
	const double z = ((double) num_matches) / states_tried;
	
	// check the bounds
	assert( 0.0 <= z );
	assert( z <= 1.0 );
		
	return z;
}


double t_consciousness::prob_s1_given_mu0( const register unsigned int s1, const register bitmask Smask, const register unsigned int mu0, const bitmask MUmask, const unsigned int MUsize )
// calculates p(s1|mu0) for a specified s1, mu0
// the user must pass MUsize, the number of bits ON in MUmask, for efficiency reasons
{
	static unsigned int old_Vmask=UINT_MAX;
	// for doing the fast log2
	static const int MultiplyDeBruijnBitPosition[32] = 
	{
		0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
		8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
	};
	

	// verify the masks are valid
	// note that these will SOMETIMES BE THE FULL MASK
	assert( MUmask <= FULL_MASK );
	assert( Smask <= FULL_MASK );
	
	// check that mu0 is a subset of MUmask
	// check that s1 is a subset of Smask
	assert( (mu0 | MUmask) == MUmask );
	assert( is_valid_mu1(s1,Smask) );
	assert( MUsize > 0 );
	
	
	// If MUmask is FULL_MASK, go ahead and return that one.
	if( MUmask == FULL_MASK )
	{
		if( (states[mu0] & Smask) == s1 )
			return 1.0;
		else
			return 0.0;
	}
	unsigned int mu0s1_matches=0;
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Nov 9, 2010	
	// COMPLEX OPTIMIZATION: Looping ONLY over the x0s that match mu0
	// We create the x0 by weaving state mu0 and v0.  v0 is in Vmask, and mu0 is in MUmask
	// This is faster because:	(1) Only loop through 2^|V| states instead of 2^|X| states
	// Major parts of this optimization:
	// -- caching the Vnodes (prob_s1_given_mu0__Vnodes)
	// -- knowing the bit within v0 that's changed using fast_log2( (v0-1) ^ v0 )
	// -- Setting the v0 bits in ONLY TWO STEPS
	//		x0 = ( (x0 >> x0_pivot) | 1 ) << x0_pivot;
	// -- Setting the mu0 bits in ONE STEP
	//		x0 |= mu0;
	////////////////////////////////////////////////////////////////////////////////////////////////
	
	// get all of the bits OUTSIDE of MUmask
	const unsigned int Vmask = FULL_MASK & (~MUmask);
	const unsigned int Vsize = numunits - MUsize;
	const unsigned int Vstates = (1 << Vsize);
	
	// just a shorter alias to the class variable
	unsigned int* Vnodes = this->prob_s1_given_mu0__Vnodes;

	if( Vmask != old_Vmask )
	{
		mask2array( Vmask, Vnodes );
		
		// set the current Vmask to old_Vmask
		old_Vmask = Vmask;
	}


//	 cout << "-----------------mask=" << binary_repr( MUmask, numunits ) << "-----------------------------" << endl;	
//	 cout << "------------------mu0=" << binary_repr( mu0, numunits ) << "-----------------------------" << endl;
	 
//	 cout << "Vnodes: \t";
//	 for( int i=0; i<Vsize; i++ ) {
//	 cout << Vnodes[i] << " ";
//	 } cout << "\t\t Vsize=" << Vsize << endl;
	
	
	// We'll be Bit-shifting x0 a bunch, so lets make it a register.
	register unsigned int x0=mu0;

	// We have to do this one OUTSIDE OF THE LOOP because fast_log2() breaks on it.
	if( ((states[mu0])&Smask) == s1 )
		mu0s1_matches += 1;
	
	// foreach v0 state...
	for( unsigned int v0=1; v0<Vstates; v0++ )
	{
		////////////////////////////////////////////////////////////////////////////////////////////////
		// Nov 9, 2010	
		// OPTIMIZATION: We compute log2 for integers really fast.
		// We do this to get the bit that changed between (v0-1) and v0.
		// This method uses a table defined at the top of the function, and it ONLY WORKS because (v0-1)^(v0) is
		// ONE LESS than a POWER OF TWO.
		// from: http://www-graphics.stanford.edu/~seander/bithacks.html
		// aka: const unsigned int leftmost_bit_that_changed = fast_log2( (v0-1) ^ v0 );		
		////////////////////////////////////////////////////////////////////////////////////////////////
		const unsigned int leftmost_bit_that_changed = MultiplyDeBruijnBitPosition[(uint32_t)( ((v0-1) ^ v0) * 0x07C4ACDDU) >> 27];

//		cout << binary_repr( v0-1, Vsize ) << " (" << v0-1 << ") ^ " << binary_repr( v0, Vsize ) << " (" << v0 << ") = " << binary_repr( (v0-1) ^ v0, Vsize );

		assert( leftmost_bit_that_changed <= (Vsize-1) );
		
		const unsigned int x0_pivot = Vnodes[ leftmost_bit_that_changed ];

//		cout << "\t guess_bit=" << guess_bit << "\t x0_pivot=" << x0_pivot << endl;
		
		// set everything to the RIGHT of x0_pivot to ZERO; put a 1 at x0_pivot; then shift everything back in place.
		// then re-apply the mu0 mask
		x0 = (( (x0 >> x0_pivot) | 1 ) << x0_pivot) | mu0;
		
		// this should now be true
		assert( (x0&MUmask) == mu0 );
		
		if( ((states[x0])&Smask) == s1 )
			mu0s1_matches += 1;
	}
	//	cout << "=====================done.========================" << endl << endl;
	
		
	// if there are no intersections, then the probability is zero.
	if( mu0s1_matches == 0 )
		return 0.0;
	
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Nov 5, 2010
	// OPTIMIZATION: efficiently calculating the total number of times mu0 will match.
	// This is faster because:	(1) can calculate #mu0_matches in one step;
	//							(2) can collapse the if statements together
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Full details are on page 124 of Virgil's Laboratory Notebook.
	// number_of_mu0_matches = 2^|X| / 2^|MUmask| = 2^numunits / 2^MUsize = 2^(numunits - MUsize)
	////////////////////////////////////////////////////////////////////////////////////////////////
	const double mu0_matches = (double) (1 << Vsize);	
	
	// this should be true
	assert( mu0s1_matches <= mu0_matches );
	
	const double z = (double) mu0s1_matches / mu0_matches;
	
	// check bounds for z
	assert( 0.0 <= z );
	assert( z <= 1.0 );
	
	return z;
}


double t_consciousness::prob_s1_given_mu0__slower( const unsigned int s1, const bitmask Smask, const unsigned int mu0, const bitmask MUmask, const unsigned int MUsize )
// calculates p(s1|mu0) for a specified s1, mu0
// the user must pass MUsize, the number of bits ON in MUmask, for efficiency reasons
// this is the less optimized version of prob_s1_given_mu0.  It is ~25% slower than the standard prob_s1_given_mu0()
{
	// verify the masks are valid
	// note that these will SOMETIMES BE THE FULL MASK
	assert( MUmask <= FULL_MASK );
	assert( Smask <= FULL_MASK );
	
	// check that mu0 is a subset of MUmask
	// check that s1 is a subset of Smask
	assert( (mu0 | MUmask) == MUmask );
	assert( is_valid_mu1(s1,Smask) );
	assert( MUsize > 0 );

	// If MUmask is FULL_MASK, go ahead and return that one.
	if( MUmask == FULL_MASK )
	{
		if( (states[mu0] & Smask) == s1 )
			return 1.0;
		else
			return 0.0;
	}	
	
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Nov 6, 2010	
	// OPTIMIZATION: Incrementing by the LOWEST bit OUTSIDE of MUmask
	// This is faster because:	(1) If the right-most bit is INSIDE MUmask, then we can increment further than 1 per x0.
	////////////////////////////////////////////////////////////////////////////////////////////////
	
	// get all of the bits OUTSIDE of MUmask
	const unsigned int not_MUmask = FULL_MASK & (~MUmask);
	unsigned int i_bit = 1;

	// if the i_bit is ON outside of the MUmask, STOP
	// and use the i_bit for x0_inc
	while( (not_MUmask & i_bit) == 0 )
	{	
		// Move the '1' bit to the left.
		i_bit <<= 1;			
	}

	unsigned int mu0s1_matches=0;
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Nov 5, 2010	
	// OPTIMIZATION: Skipping some initial and final x0 states that will never match mu0
	// This is faster because:	(1) can skip some initial iterations of the loop for free
	//							(2) can skip some final iterations of the loop for free
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Full details are on page 127 of Virgil's Laboratory Notebook.
	// The lowest value of x0 that bitmask matches mu0 is mu0 itself. Thus can start iterating x0 at mu0 instead of 0.
	// The highest value of x0 that bitmask matches mu0 is all 1s outside MUmask, with mu0 inside MUmask
	////////////////////////////////////////////////////////////////////////////////////////////////
	
	// these expressions are the simplest expressions for the bounds I could come up with
	const unsigned int x0_start = mu0;
	const unsigned int x0_stop = (FULL_MASK & (~MUmask)) | mu0;
	const unsigned int x0_inc = i_bit;
//	const unsigned int x0_start=0, x0_stop=FULL_MASK, x0_inc=1;
	
	for( unsigned int x0=x0_start; x0<=x0_stop; x0 += x0_inc )	
	{
		// if this x0 state matches mu0 AND matches s1...
		// After trying both, checking mu0 first is faster.
		if( ((x0&MUmask) == mu0) && (((states[x0])&Smask) == s1) )
				mu0s1_matches += 1;
	}

	// if there are no intersections, then the probability is zero.
	if( mu0s1_matches == 0 )
		return 0.0;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// Nov 5, 2010
	// OPTIMIZATION: efficiently calculating the total number of times mu0 will match.
	// This is faster because:	(1) can calculate #mu0_matches in one step;
	//							(2) can collapse the if statements together
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Full details are on page 124 of Virgil's Laboratory Notebook.
	// number_of_mu0_matches = 2^|X| / 2^|MUmask| = 2^numunits / 2^MUsize = 2^(numunits - MUsize)
	////////////////////////////////////////////////////////////////////////////////////////////////
	const double mu0_matches = (double) (1 << (numunits - MUsize));	
	
	// this should be true
	assert( mu0s1_matches <= mu0_matches );
	
	const double z = (double) mu0s1_matches / mu0_matches;

	// check bounds for z
	assert( 0.0 <= z );
	assert( z <= 1.0 );
	
	return z;
}

double t_consciousness::prob_mu0_given_s1( const t_state& mu0, const t_state& s1 )
// just an alias to the bitmask version
{
    return this->prob_mu0_given_s1( mu0.value, mu0.mask, s1.value, s1.mask, mu0.size() );
}

double t_consciousness::prob_mu0_given_s1( const unsigned int mu0, const bitmask MUmask, const unsigned int s1, const bitmask Smask, const bitmask MUsize )
// calculates p(mu0|s1) from p(s1|mu0)
// November 7, 2010 -- OPTIMIZATION: Takes MUsize as input so we don't have to calculate it everytime we call this function.  This saved ~10% of compute time.
{
	// assert the masks are valid
	assert( MUmask <= FULL_MASK );
	assert( Smask <= FULL_MASK );

	// assert the states are valid
	assert( is_valid_mu1(s1,Smask) );
	assert( (mu0|MUmask) == MUmask );
	
	const double top = prob_s1_given_mu0( s1, Smask, mu0, MUmask, MUsize );	
	
	if( top == 0.0 )
		return 0.0;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// Nov 5, 2010	
	// OPTIMIZATION: efficiently calculating: \sum_{\mu_0 \in M_0} p( s1|\mu0 )
	// \sum_{\mu_0 \in M_0} p( s1|\mu0 ) = 2^|MUmask| * p(S1 = s1)
	// This is faster because:	(1) no call to the expensive function p(s1|mu0)
	//							(2) no for loop over the expensive function
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Full details are on page 126 of Virgil's Laboratory Notebook.
	// \sum_{\mu_0 \in M_0} p( s1|\mu0 ) = bottom = 2^|MUmask| * p(S1 = s1)
	////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	const double this_prob_s1 = prob_s1( s1, Smask );
	const double bottom = this_prob_s1 * (1 << MUsize);
	
	const double z = top / bottom;
	
	assert( 0.0 <= z );
	assert( z <= 1.0 );
	
	return z;
}

double t_consciousness::prob_mu1_given_s0__anticond( const unsigned int mu1, const unsigned int s0, const unsigned int M1mask, const unsigned int S0mask, const unsigned int anticond_mask, const unsigned int M1size, const unsigned int S0size, const unsigned int anticond_size )
// returns p(mu1|s0 \anticond \anticond_mask ) = p(mu1|s0 \anticond \widetilde{S})
// this is defined in terms of Bayes rule.
{
	// assert S0 and anticond are disjoint
	assert( (S0mask & anticond_mask) == 0 );
	// assert S0 and anticond collectively pave the system
	assert( (S0mask | anticond_mask) == FULL_MASK );
	
	double z=1.0;
	unsigned int* restrict nodes = new unsigned int[M1size];
	subset2array(M1mask, nodes);
	
	for( int i=0; i<M1size; i++ ) {
		const uint n_index = nodes[i];
		// for this mu1 only keep the state of the node at location n_index
		const unsigned int n1 = (mu1 & (1<<n_index));   // Keep only the [0,1] state at value n_index.  Everything rest goes to zero.
		
		
		//z *= prob_n1_given_mu0__anticond( n1, n_index, mu0, partmask, partsize );
		z *= prob_n1_given_s0__anticond( n1, n_index, s0, S0mask, anticond_size );		
		
		// if the probability becomes 0.0, break out of the loop.  We must BREAK
		// instead of return 0.0 because we still have to delete [] the nodes.
		if( z == 0.0 )
			break;
	}
	
	delete [] nodes;
	
	return z;
	
	
}


double t_consciousness::prob_mu1_given_mu0__anticond( const unsigned int mu1, const unsigned int mu0, const unsigned int partmask, const unsigned int partsize )
// returns p(mu1|mu0 \anticond \widetilde{M} )
{

	assert( (mu0|partmask) == partmask );	
//	assert( is_valid_mu1(mu1,partmask) );
	
	assert( numunits_in_mask(partmask) == partsize );
	
	double z=1.0;
	const uint anticond_size = numunits - partsize;
	unsigned int* restrict nodes = new unsigned int[partsize];
	subset2array(partmask, nodes);
	
	for( int i=0; i<partsize; i++ ) {
		const uint n_index = nodes[i];
//		const uint n1 = (mu1 >> n_index)&1;
		const uint n1 = (1 << n_index) & mu1;
		
		assert( (n1 == 0) || n1 == (1 << n_index) );
		
		z *= prob_n1_given_s0__anticond( n1, n_index, mu0, partmask, anticond_size );

		// if the probability becomes 0.0, break out of the loop.  We must BREAK
		// instead of return 0.0 because we still have to delete [] the nodes.
		if( z == 0.0 )
			break;
	}

	delete [] nodes;
	
	return z;
}



double t_consciousness::prob_s1( const unsigned int s1, const bitmask Smask )
// returns p( S1 = s1 )
{
	static unsigned int old_s1=0, old_Smask=0;
	static double old_z=-1.0;
	
	// if we're requesting the same value, return from the cache.
	if( s1==old_s1 && Smask==old_Smask )
		return old_z;
	
	// verify partmask
	assert( Smask <= FULL_MASK );
	// verify s1 is a subset of partmask
	assert( (s1 | Smask) == Smask );
	
	unsigned int num_x0s_going_to_s1=0;

	// Foreach output state x1...	
	for( int x1_index=0; x1_index<NUM_X1_STATES; x1_index++ )
	{
		const unsigned int x1 = x1_states[x1_index];

		// if this x1 matches s1, then add it's #x0s to s1
		if( (x1&Smask) == s1 )
			num_x0s_going_to_s1 += this->num_x0s_to_this_x1[x1_index];
	}
	
	const double z = (double) num_x0s_going_to_s1 / (double) this->numstates;
	
	old_s1 = s1;
	old_Smask = Smask;
	old_z = z;
	
    assert( 0.0 <= z );
    assert( z <= 1.0 );
    
	return z;
}





double t_consciousness::prob_n1_given_s0__anticond( const unsigned int n1, const unsigned int n_index, const unsigned int s0, const unsigned int S0mask, const unsigned int anticond_size )
// using a generalized version of perturbing the wires from Balduzzi returns p( n1| s0 \anticond ~S0mask )
// n_index is the index WITHIN the whole network.    0 <= n_index < numunits
// s0 is NOT shifted to the right.  It occupies the same bits as S0mask
// anticond_size is there solely for efficiency purposes
{

//	cerr << "i'm here!!!" << endl;
//	cerr.flush();

//	assertmsg( (s0 & (~S0mask)) == 0, "Error: s0 has a entry ON outside of S0mask." );
	assertmsg( (s0 | S0mask) == S0mask, "Error: s0 has a entry ON outside of S0mask." );
	
	
	assert( 0 < anticond_size );
	assert( n_index < numunits );
	
	// assert n1 is blank or has a 1 at n_index
	assert( n1==0 || n1 == (1<<n_index) );

	// anticond_mask is everything outside of S0mask
	const unsigned int anticond_mask = (FULL_MASK & (~S0mask));
	unsigned int num_matches=0;
	
	// assert that anticond_mask makes sense
	assertmsg( numunits_in_mask(anticond_mask) == anticond_size, "derived anticond_mask and anticond_size didn't agree." );
	
	//foreach x0 state...
	for( uint x0=0; x0<numstates; x0++ )
	{
		// skip all x0 states that dont match s0.
		if( (x0 & S0mask) != s0 )
			continue;
		
		const uint x1 = states[x0];

		// get the n_index'th node 
		const unsigned int temp_n1 = x1 & (1<<n_index);
		
		num_matches += int(temp_n1 == n1);
	}
	
	
	// divide by the number of states we tried -- 2**(SIZE_OF_ANTICOND_MASK)
	const double states_tried = (double) (1 << anticond_size);
	
	// z = p(n1|mu0 \anticond A )
	const double z = ((double) num_matches) / states_tried;
	
	return z;	
}



double t_consciousness::H_M1( const unsigned int S )
// returns H[M1] for a given subset M
{
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Create the composite mask for caching purposes
	////////////////////////////////////////////////////////////////////////////////////////////////
	// if the cache_part_entropies is not NONE, then it's a valid cache
	if( FLAG__CACHE_M1 && H_M1_cache[S] >= 0.0 )
		return H_M1_cache[S];
	////////////////////////////////////////////////////////////////////////////////////////////////
	
	unsigned int numnodes = numunits_in_mask(S);

	ASSERT( 1 <= numnodes );
	ASSERT( numnodes <= numunits );
	
	//unsigned int numnodes = S.numnodes();		//numunits_in_subset( part_mask );
	unsigned int d = 1 << numnodes;		//d=2**part.size()

	unsigned int* restrict nodes = new unsigned int[numnodes];
	mask2array( S, nodes );

	
	// define a block of memory to hold the number of times a mu1 occurs
	unsigned int* restrict mu1_instances = (unsigned int*) calloc(d, sizeof(unsigned int));	
	assert( mu1_instances != NULL );

/*	
	// the bitmask for the MUstate
	unsigned int MUmask = (1 << numnodes) -1;

	// the bits that will be UNCHANGED.  The MU bits already INSIDE MU.
	unsigned int MUINmask = S & MUmask;
		
	// the bitmask for the bits that will transfer from outside MU to inside MU (SOURCE bits)
	unsigned int XOUTmask = S & ~MUmask;

	// the bitmask for the destination bits that will be received from outside MU (DESTINATION bits)	
	unsigned int MUDESTmask = MUmask & ~MUINmask;	

	cout << "XOUTmask=" << XOUTmask << endl;
	cout << "MUDESTmask=" << MUDESTmask << endl;
	cout.flush();
	unsigned int XOUTsize = numunits_in_mask(XOUTmask);
	cout << "XOUTsize=" << XOUTsize << endl;
	cout.flush();
	
	unsigned int* restrict XOUTnodes = new unsigned int[XOUTsize];
	unsigned int* restrict MUDESTnodes = new unsigned int[XOUTsize];
	
	mask2array( XOUTmask, XOUTnodes );
	mask2array( MUDESTmask, MUDESTnodes );
*/

	//foreach x0...
	unsigned int temp_mu1=0;
	unsigned int old_M1_mask=UINT_MAX, cur_M1_mask=0;

	for(uint temp_x0=0; temp_x0<numstates; temp_x0++ )
	{
		// remove the bits in x1 that aren't part of the Xstate
		cur_M1_mask = states[temp_x0] & S;

		if( cur_M1_mask != old_M1_mask ) {
			Xstate2MUstate( nodes, numnodes, cur_M1_mask, temp_mu1 );
//			temp_mu1 = Xstate2MUstate_exchg( XOUTnodes, MUDESTnodes, XOUTsize, MUINmask, cur_M1_mask );
			old_M1_mask = cur_M1_mask;
		}
		//printf("%i	%i\n",temp_mu0,temp_mu1);
		mu1_instances[temp_mu1] += 1;
	}
/*	
	// don't need the nodes anymore.
	delete [] XOUTnodes;
	delete [] MUDESTnodes;
*/
	
	const double Dnumstates = (double) this->numstates;
	double z=0.0, sum_prob=0.0;
	
	for(int mu1=0; mu1<d; mu1++)
	{
		if( mu1_instances[mu1] == 0 )
			continue;
		
		const double prob_mu1 = (double) mu1_instances[mu1] / Dnumstates;
	
		//Make each row a conditional probability by dividing by #mu1's/
		//			//Now run over both to return the information using the probabilities of p(mu0,mu1) of p(mu0|mu1)
		//			z -= ( mu0mu1_instances[mu1][mu0] / Dnumstates ) * log2( (mu0mu1_instances[mu1][mu0] / Dmu1_instances) );
		z -= ( prob_mu1 ) * log2( prob_mu1 );
		
		sum_prob += prob_mu1;
		
	}
	
	assert( fequals(sum_prob, 1.0) );
	delete [] nodes;
	
	// free the memory in both.
	free( mu1_instances );
	
	if( z < 0.0 ) {
		if( myceil(z) < 0.0 ) {
			cerr << "function: H_M1()" << endl;
			cerr << "Warning: Had negative entropy in entropy_of_part() that didn't round to 0.0 (" << myceil(z) << ").  Maybe decimal precision is too high?" << endl;
			assert( 0.0 <= myceil(z) );
		}
		z = 0.0;
	}

	
	// Store the part entropy in the cache, if we're doing that.
	if( FLAG__CACHE_M1 ) {
		
		// if this is a brand new subset, flush the old cache.
		ASSERT( this->H_M1_cache[S] != this->H_M1_cache[S] );
		
		// set the part entropy
		this->H_M1_cache[ S ] = z;
	}
	
	// DO NOT DO ANY ROUNDING HERE!	
	ASSERT( 0.0 <= z );
	
	return z;
}


double t_consciousness::H_A0_B1( const unsigned int A0mask, const unsigned int B1mask )
// calculate the JOINT ENTROPY of a variables A0 and B1.
// aka: H( A0, B1 )
{
	// take care of the null cases.
	if( A0mask == 0 )	return H_M1( B1mask );
	if( B1mask == 0 )	return H_M0( A0mask );
	
	assert( A0mask && A0mask <= FULL_MASK );
	assert( B1mask && B1mask <= FULL_MASK );

	const unsigned int Asize = numunits_in_mask( A0mask );
	const unsigned int Bsize = numunits_in_mask( B1mask );
		
	const unsigned int A0states = 1 << Asize;		//d=2**part.size()
	const unsigned int B1states = 1 << Bsize;		//d=2**part.size()	

	const double Dnumstates = this->numstates;
	
	unsigned int* restrict a0b1_instances = (unsigned int*) calloc(A0states * B1states, sizeof(unsigned int));
//	double* restrict a0b1_instances = (double*) calloc(A0states * B1states, sizeof(double));
	assert( a0b1_instances != NULL );
	
	
	unsigned int total=A0states*B1states;
//	cout << "A0states=" << A0states << "  \t\t B1states=" << B1states << endl;
//	for( int i=0; i<total;i++ ) {
//		cout << "i=" << i << " \t = " << a0b1_instances[i] << endl;
//	}
	
	unsigned int* restrict Anodes = new unsigned int[Asize];
	unsigned int* restrict Bnodes = new unsigned int[Bsize];
	
	mask2array( A0mask, Anodes );
	mask2array( B1mask, Bnodes );
	
	// create the big matrix
	unsigned int count=0;
	for( int x0=0; x0<=FULL_MASK; x0++ )
	{
		const unsigned int x1 = states[x0];
		
		unsigned int a0=0, b1=0;
		
		// make the a0 state
		Xstate2MUstate( Anodes, Asize, x0, a0 );
		// make the b1 state
		Xstate2MUstate( Bnodes, Bsize, x1, b1 );

		assert( a0 <= A0mask );
		assert( b1 <= B1mask );

		// increment this entry in the matrix
//		a0b1_instances[ (a0*A0states)+ b1 ] += (1.0/Dnumstates);
		a0b1_instances[ (a0*B1states)+ b1 ] += 1;
//		cout << "\t [" << (a0*B1states)+ b1 << "] += 1 = " << a0b1_instances[ (a0*B1states)+ b1 ] << endl;
//		cout << "adding " << 1.0/Dnumstates << endl;
		count += 1;
	}
	
//	unsigned int t=A0states*B1states;
//	for( int i=0; i<=t;i++ )
//		cout << "i=" << i << " \t = " << a0b1_instances[i] << endl;

	
//	cout << "count=" << count << endl;
	assert( count == this->numstates );

	// don't need the nodes anymore.
	delete [] Anodes;
	delete [] Bnodes;
	
//	cout << "=====================" << endl;
	double z=0.0, summ_prob=0.0;
//	double summ_prob=0.0;
	for( int a0=0; a0<A0states; a0++ ) {
		for( int b1=0; b1<B1states; b1++ ) {
			
			if( a0b1_instances[ (a0*B1states) + b1 ] == 0 )
				continue;
			
			const double this_prob = ((double) a0b1_instances[ (a0*B1states) + b1 ]) / Dnumstates;
			
//			cout << "sum_prob=" << summ_prob << " \t +" << this_prob << endl;
			summ_prob += this_prob;
			
			z -= this_prob * log2(this_prob);

		}
	}
	
	// don't need the matrix anymore
	free( a0b1_instances );

	// remove any -0.0's, etc.
	if( fequals(z,0.0) )
		z = 0.0;

	// sanity checks
//	cout << "z=" << z << "\t summ_prob=" << summ_prob << endl;
	assert( fequals(summ_prob,1.0) );
	assert( 0.0 <= z );
	assert( z <= numunits );
	
	
	return z;
}
double t_consciousness::I_A0_B1( const unsigned int A0, const unsigned int B1 )
// the mutual information of
// I(A_0:B_1) = H(B_1) - H(B_1|A_0)
{

//	cerr << "doing term1..." << endl;
//	double term1 = H_M1( B1 );
//	assert( 0.0 <= term1 );
//	double term2 = H_M1_GIVEN_S0( B1, A0 );
//	assert( 0.0 <= term2 );
	
//	double z = term1 - term2;

	if( A0 == FULL_MASK && B1 == FULL_MASK )
		return H_X1;

	// assert A0 and B1 are valid
	assert( A0 && A0 <= FULL_MASK );
	assert( B1 && B1 <= FULL_MASK );	
	
    double z = H_M0(A0) - H_M0_GIVEN_S1(A0,B1);
    
	if( fequals(z,0.0) )
		z = 0.0;
	
	assert( 0.0 <= z );
	
	return z;
}



double t_consciousness::I_A0_B1_GIVEN_C1( const unsigned int A0, const unsigned int B1, const unsigned int C1 )
// computes conditional mutual information I( A0 : B1 | C1 )
{
	// if C1 is emptry, goto the simpler version of this function
	if( C1 == 0 )	return I_A0_B1( A0, B1 );
	
	// assert A0, B1, C1 are valid
	assert( A0 && A0 <= FULL_MASK );
	assert( B1 && B1 <= FULL_MASK );
	assert( C1 && C1 <= FULL_MASK );	
		
	// assert no overlap
	assert( (B1 & C1) == 0 );

	double z = H_A0_B1(A0,C1) + H_M1(B1|C1) - H_A0_B1(A0, B1|C1)  - H_M1(C1);
	
	if( fequals(z,0.0) )
		z = 0.0;
	
	assert( 0.0 <= z );
    
    return z;
}

double t_consciousness::I_A0_B1_GIVEN_C0_D1( unsigned int A0, unsigned int B1, unsigned int C0, unsigned int D1 )
// compute the conditional mutual information condition on C_0 and D_1
// I(X;Y|Z) = H(X|Z) − H(X|Y,Z)
// I(A;B1|C0,D1) = H(A0|C0,D1) − H(A0|B1,C0,D1)
//               = H(A0|D1) - H(A0|B1,D1)                 

{
	// if C0 or D1 isn't present then we use the simpler forms of this function.
	if( C0 == 0 && D1 == 0 ) return I_A0_B1( A0, B1 );
	else if( C0 == 0 )		 return I_A0_B1_GIVEN_C1( A0, B1, D1 );
	else if( D1 == 0 )		 return I_A0_B1_GIVEN_C0( A0, B1, C0 );

	assert( A0 && A0 <= FULL_MASK );
	assert( B1 && B1 <= FULL_MASK );
	assert( C0 && C0 <= FULL_MASK );
	assert( D1 && D1 <= FULL_MASK );
	
	// assert no overlap between A0,C0 and B1, D1
	assert( (A0 & C0) == 0 );
	assert( (B1 & D1) == 0 );
	
	
//	double z = H_M0_GIVEN_S1( A0, D1 ) - H_M0_GIVEN_S1( A0, B1 | D1 );
	// H(X,Z) + H(Y,Z) − H(X,Y,Z) − H(Z)
	const unsigned int X = A0, Y = B1;
	double z = H_A0_B1(X|C0,D1) + H_A0_B1(C0,Y|D1) - H_A0_B1(X|C0,Y|D1 ) - H_A0_B1(C0,D1);
	
	if( fequals(z,0.0) )
		z = 0.0;
	
	assert( 0.0 <= z );

	return z;
}

double t_consciousness::I_A0_B1_GIVEN_C0( unsigned int A0, unsigned int B1, unsigned int C0 )
// the mutual information of
// I(A_0:B_1|C_0) = H(B_1|C_0) - H(B_1|A_0,C_0)
// I(A0:B1|C0) = H(A0|C0) - H(A0|B1,C0)
//             = H(A0) - H(A0|B1)
{
	assert( A0 && A0 <= FULL_MASK );
	assert( B1 && B1 <= FULL_MASK );
	assert( C0 && C0 <= FULL_MASK );
	
	// assert no overlap between A0 and C0
	assert( (A0 & C0) == 0 );
		
	double z1 = H_M1_GIVEN_S0( B1, C0 ) - H_M1_GIVEN_S0( B1, A0 | C0 );
//	double z2 = H_M0( A0 ) - H_M0_GIVEN_S1( A0, B1 );
//	double z3 = numunits_in_mask(A0) - H_M0_GIVEN_S1( A0, B1 );	
	
//	cout << "z1=" << z1 << "\tz2=" << z2 << "\tz3=" << z3 << endl;
	
	double z = z1;
	
	if( fequals(z,0.0) )
		z = 0.0;
	
	if( z < 0.0 ) {
//		cerr << "both_mask=" << binary_repr(both) << endl;
//		cerr << "H[B1|C0]=" << term1 << "\tH[B1|A0,C0]=" << term2 << "\tdiff=" << term1-term2 << endl;
		cerr << "A=" << binary_repr(A0) << "\tB=" << binary_repr(B1) << "\tC=" << binary_repr(C0) << endl;
		cerr << "z=" << z << endl;
		exit(-1);
	}
	
	assert( 0.0 <= z );
	
	return z;
	
}



inline unsigned int t_consciousness::cachekey( unsigned int S0, unsigned int M1 )
// returns the cache key for H[M1|S0]
{
	// Cache key has format: S0,M1
	assert( numunits <= 16 );
	
	assert( M1 >= 1 );
	assert( M1 <= FULL_MASK );
	
	// S0 can be zero when calculating H(M1) (with no conditioning)
	assert( S0 >= 0 );
	assert( S0 < FULL_MASK );
	
//	unsigned long z = (S0 << numunits) | M1;
	unsigned int z = (S0 << numunits) | M1;	
	
	if( S0 )
		assert( z > FULL_MASK );

	return z;
}

double t_consciousness::H_M0_GIVEN_S1( const unsigned int M0, const unsigned int S1 )
// !!! WARNING: This function may not be correct. !!!
{

	// if M and S are the same, return H[M0|M1]
	if( M0 == S1 )
		return H_M0_GIVEN_M1( M0 );

	const unsigned int M0size=numunits_in_mask(M0), S1size=numunits_in_mask(S1);
	
	
	const unsigned int M0d = 1 << M0size;		//d=2**part.size()
	const unsigned int S1d = 1 << S1size;		//d=2**part.size()	
	
	// now define a 2-dimensional array as a single contiguous block of memory
	// this requires us to create a big-ass 1-d array but we PRETEND it is a 2-d array
	// However, we must now perform subscript calculations manually, accessing the i,jth element via:
	// myarray[i * ncolumns + j]
	
	unsigned int* restrict mu0mu1_instances = (unsigned int*) calloc(M0d * S1d, sizeof(unsigned int));
	assert( mu0mu1_instances != NULL );	
	
	// the bitmask for the MUstate
	unsigned int MU0mask = M0d -1;
	unsigned int MU1mask = S1d -1;
	
	unsigned int* restrict Mnodes = new unsigned int[M0size];
	unsigned int* restrict Snodes = new unsigned int[S1size];
	
	mask2array( M0, Mnodes );
	mask2array( S1, Snodes );

	
	unsigned int mu0=0, mu1=0;
	//	unsigned int old_X0OUT_mask=UINT_MAX, old_X1OUT_mask=UINT_MAX;
	unsigned int old_M0_mask = UINT_MAX, old_S1_mask=UINT_MAX;
	
	for(unsigned int temp_x0=0; temp_x0<numstates; temp_x0++ )
	{
		// clear the bits not in the part
		//		const unsigned int x1 = states[temp_x0] & M1;
		//		temp_x0 &= S0;
		//		mu0 = temp_x0 & S0;
		//		mu1 = states[temp_x0] & M1;
		
		
//		const unsigned int cur_M0_mask = temp_x0 & M0;
//		const unsigned int cur_S1_mask = states[temp_x0] & S1;

		const unsigned int cur_M0_mask = temp_x0;
		const unsigned int cur_S1_mask = states[temp_x0];

		if( cur_M0_mask != old_M0_mask ) {		
			
			//			Xstate2MUstate( Snodes, S0size, temp_x0, mu0 );
			Xstate2MUstate( Mnodes, M0size, cur_M0_mask, mu0 );			
			//			old_X0OUT_mask = cur_X0OUT_mask;
			old_M0_mask = cur_M0_mask;
		}
		
		if( cur_S1_mask != old_S1_mask ) {			
			//			Xstate2MUstate( Snodes, S0size, temp_x0, mu0 );
			Xstate2MUstate( Snodes, S1size, cur_S1_mask, mu1 );			
			//			old_X0OUT_mask = cur_X0OUT_mask;
			old_S1_mask = cur_S1_mask;
		}
		
		//		assert( mu0 <= cur_S0_mask );
		//		assert( mu1 <= cur_M1_mask );
		assert( mu0 <= MU0mask );
		assert( mu1 <= MU1mask );
		
		
		mu0mu1_instances[ (mu0*S1d)+mu1 ] += 1;
	}
	
	
	// don't need the nodes anymore.
	delete [] Mnodes;
	delete [] Snodes;
	
	
	const double Dnumstates = (double) this->numstates;
	double H_M1_S0=0.0;
	double sum_prob=0.0;
	
	// create H( M0 S1 )
	for( int mu0=0; mu0<M0d; mu0++ )
	{
		for(int mu1=0; mu1<S1d; mu1++)
		{
			if( mu0mu1_instances[(mu0*S1d)+mu1] == 0 )
				continue;
			
			const double prob_mu0mu1 = mu0mu1_instances[(mu0*S1d)+mu1] / Dnumstates;
			
			
			// H( XY ) = prob(mu0,mu1) * log2( p(mu1|mu0) )
			H_M1_S0 -= ( prob_mu0mu1 * log2( prob_mu0mu1 ) );
			sum_prob += prob_mu0mu1;
		}
		
	}
	
	
	//	cerr << "made it here [1]" << endl;
	//	cerr.flush();	
	
	// free the memory in both.
	free( mu0mu1_instances );
	
	assert( fequals(sum_prob, 1.0) );
	
	
	//	cerr << "made it here [2]" << endl;
	//	cerr.flush();	
	
	// H( X|Y ) = H(XY) - H(Y)
	//	cerr << "- H[M1,S0]=" << H_M1_S0 << endl;
	//	cerr.flush();
	double z = H_M1_S0 - H_M1( S1 );
	
	if( z < 0.0 ) {
		if( myceil(z) < 0.0 ) {
			cerr << "function: entropy_of_part__ELEMENTS" << endl;
			cerr << "Warning: Had negative entropy in entropy_of_part() that didn't round to 0.0 (" << myceil(z) << ").  Maybe decimal precision is too high?" << endl;
			assert( 0.0 <= myceil(z) );
		}
		z = 0.0;
	}
	
	
	// DO NOT DO ANY ROUNDING HERE!
	ASSERT( 0.0 <= z );
	ASSERT( z <= numunits );
	
	return z;
}



double t_consciousness::H_M1_GIVEN_S0( unsigned int M1, unsigned int S0 )
// Calculates H[M_1|S_0] for different values of M and S.
{
	////////////////////////////////////////////////////////////////////////////////////////////////
	// if the S0-state is the FULL MASK, then the entropy of M1 is zero.  Return zero.
	////////////////////////////////////////////////////////////////////////////////////////////////
	if( S0 == FULL_MASK )
		return 0.0;
	
	if( S0 == 0 )
		return H_M1( M1 );
	
	assert( S0 != FULL_MASK );
	assert( S0 > 0 );
	assert( M1 > 0 );
	
	
	return H_A0_B1( S0, M1 ) - H_M0(S0);
}

			
vector<t_phi_result> t_consciousness::bracket_maincomplexes()
// discovers the main complexes via enumerating all supersets.
{
	vector<t_phi_result> MCs;
	MCs.clear();
	
	// if there's only one x1 state, ei=0.0 and phi=0.0.  Period.
	if( this->NUM_X1_STATES <= 1 )
		return MCs;
	
	//foreach subset...calculate the MIP and the PHI
	//we start from the FULL_MASK because we can often skip steps that way
	for( uint subset=FULL_MASK; subset>=1; subset-- )
	{
		//OPTIMIZATION: only attempt all partitions to find the MIPs if it's possible
		//for this subset to have higher or equal PHI to the current highestphi
		// UPDATE: The optimization using bracket_ei() wasn't worth the time it took to mantain it.
		
		t_phi_result candidate = bracketphi( subset );

		// MC must have a MIP
		if( candidate.MIPs.empty() )
			continue;

		// MC can't have zero phi.
		if( fequals(candidate.min_phi, 0.0 ) )
			continue;

		
		bool removed_all=true;
		
		// going BACKWARDS so the indices don't mess up when removing entries...
		for( int i=MCs.size()-1; i>=0; i-- )
		{
			// if the candidate is a STRICT SUBSET of a MC....
			if( candidate.subset.is_subset( MCs[i].subset ) )
			{
				// if the candidate has GREATER OR EQUAL phi, remove the subset.
				if( candidate.min_phi >= MCs[i].min_phi || fequals(candidate.min_phi, MCs[i].min_phi) )
					MCs.erase( MCs.begin()+i );
				else
					removed_all=false;
			}

			// if found subsets.... and the candidate has a STRICTLY LARGER PHI, remove the subset
			else if( candidate.subset.is_superset( MCs[i].subset ) )
			{
				if( candidate.min_phi > MCs[i].min_phi && !fequals(candidate.min_phi, MCs[i].min_phi) )
					MCs.erase( MCs.begin()+i);
				else
					removed_all=false;
			}
		}

		// if removed any matching supersets as well as any matching subsets
		// add the candidate to the MCs
		if( removed_all ) {
			MCs.push_back( candidate );
			ASSERT( MCs.size() <= (numunits*numunits) );
		}
	}

	
	// sort the t_phi_results
	vec_sort( MCs );

	// note that this CAN BE EMPTY.
	
	return MCs;
}


unsigned int t_consciousness::size_of_smallest_part( const t_partition& restrict P )
// takes a partition represented as list of bitmasks, returns the size of the smallest part
{
	// min_partsie starts as the size of the first part
	unsigned int min_partsize=P.size(0);
	
	// for each part...
	for( int i=1; i<P.size(); i++ )
	{
		unsigned int partsize = P.size(i);

		ASSERT( partsize < this->numunits );
		
		// the smallest value min_partsize can ever be is 1.
		// So if get a 1, return that.
		min_partsize = min(min_partsize, partsize);

		if( min_partsize == 1 )
			return 1;
	}
	
	return min_partsize;
}

long double t_consciousness::normalized_ei( const bitmask x1, const t_partition& restrict P )
// Calculates the normalized ei for a given partition
{
    assert( is_valid_mu1(x1, FULL_MASK) );
    
    double ei = this->ei( x1, P );
//	double ei = this->bracket_ei( P );
	
	if( fequals(ei,0.0) )
		return 0.0;
	
	double norm = this->normalization( P );
	long double z = ei/norm;
	
    //	assert( z > 0.0 );
    
	if( fequals(norm,1.0) )
		z = ei;
    
    //		if( myround(z) > ei ) {
    //			cout << "P=" << P->BAREstr() << endl;
    //			cout << "ei=" << ei << endl;
    //			cout << "ei/norm=" << z << endl;
    //		}
    
	if( FLAG__MIP_METHOD == "TONONI" ) {
		assert( 1.0 <= norm );
		assert( fabs(z) <= fabs(ei) || fequals(z,ei) );
	}
    
	return z;
}


long double t_consciousness::bracket_normalized_ei( const t_partition& restrict P )
// Calculates the normalized ei for a given partition
{
	double ei = this->bracket_ei( P );
	
	if( fequals(ei,0.0) )
		return 0.0;
	
	double norm = this->normalization( P );
	long double z = ei/norm;
	
//	assert( z > 0.0 );
		
	if( fequals(norm,1.0) )
		z = ei;

//		if( myround(z) > ei ) {
//			cout << "P=" << P->BAREstr() << endl;
//			cout << "ei=" << ei << endl;
//			cout << "ei/norm=" << z << endl;
//		}

	if( FLAG__MIP_METHOD == "TONONI" ) {
		assert( 1.0 <= norm );
		assert( fabs(z) <= fabs(ei) || fequals(z,ei) );
	}

	return z;
}

double t_consciousness::normalization( const t_partition& restrict P )
/* Calculates the normalization for a given partition */
{

	double z=0.0;
	
	// if no normalization, return 1
	if( FLAG__MIP_METHOD == "NONE" ) {
//		cout << "Doing none" << endl;
		return 1.0;
	}

	// if doing the Tononi normalization, do that now.
	else if( FLAG__MIP_METHOD == "TONONI" )
	{
//		cout << "doing tononi" << endl;
		// If this is the total partition, the normalization is the number of size of that whole partition
		// DO NOT USE this->numunits here.  Because sometimes we're getting the normalization for partitions of a *subset* of the network.
		// P->size(0) is always correct.  this->numunits is WRONG when calculating normalization of a subset of the network.
		if( P.size() == 1 ) {
			z = P.size(0);
			ASSERT( z <= this->numunits );
		}
		else {
			// Normalization is (K-1) * <size of smallest part>
			z = (double) (P.size()-1) * this->size_of_smallest_part( P );
		}

	}
	
	else {
		z = 1.0;
	}
	
	return z;

}

vector<t_partition> t_consciousness::MIPs( int x1, unsigned int subset )
// Goes through all partitions in generator and finds the MIP using ei(x1/P).
{
	// Have we already found the MIP for this x1 state?
	if( !subset )
		subset = this->FULL_MASK;

	// verify that the x1 and subset are valid.
	ASSERT( this->existing_states.count(x1) );
	ASSERT( subset <= this->FULL_MASK );
	assert( (x1|subset) == subset );
	
	//foreach line in the partitions file...	
    vector<t_partition> MIPs;
	double best_MIP_score=DBL_MAX * -1.0;
    
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// BEGIN the partitions generator
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int numnodes = numunits_in_subset( subset );
//	PartitionEnumerator pgen( numnodes );
	PARTITIONS_GENERATOR.reset( numnodes );
	
	// define a partition with the maximum number of parts
	t_partition P( this->max_num_parts, subset );
	
	// initialize node_map to be as large as it could ever be, numnodes
	// Now to create the array that maps from result -> incoming partition
	unsigned int nodemap[numnodes];
	if( numnodes < this->numunits )
		subset2array(subset, nodemap);
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// END partitions generator
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// now read the parts from GENERATOR...
	while( get_next_partition(&PARTITIONS_GENERATOR, &P, nodemap, numnodes) )
	{
		
		double this_MIP_score = MIP_score( x1, P );
        if( ! (-1*DBL_MAX < this_MIP_score) )
            cout << "this_MIP_score=" << this_MIP_score << endl;

		assert( -1*DBL_MAX < this_MIP_score );
		
//		if( debug ) {
//			cerr << "P: " << P.str() << "\tMIPscore=" << this_MIP_score << "\tBestMIPscore=" << best_MIP_score << endl;
//          if( this_MIP_score > best_MIP_score ) {
//              cerr << "Was better than Best!" << endl;
//            }
//        }
		//		cout << "subset=" << subset << "\t ei=" << ei << " \t norm=" << norm << endl;
  //      }
		
		
		// if this norm EI is EQUAL TO the current one, do...
		if( fequals(this_MIP_score, best_MIP_score) )
		{
			assert( !MIPs.empty() );
			
			if( MIPs.size() < MAXIMUM_NUMBER_OF_MIPS_TO_STORE )
				MIPs.push_back( P );
			else
			{
				for( vector<t_partition>::iterator it=MIPs.begin(); it != MIPs.end(); ++it )
				{
					// if P wins a tiebreaker, replace *it with P
					if( break_MIPscore_ties( *it, P ) )
					{
						// if the current P has MORE PARTS than a MIP, replace that MIP, and stop the loop
						MIPs.erase( it );
						MIPs.push_back( P );
						it = MIPs.end();
						break;
					}
				}
			}
		}
		
		// if this partition has a HIGHER MIP score, clear previous MIPs and set partition as the MIP.
		else if( this_MIP_score > best_MIP_score ) {
//            cerr << "Setting best MIP-so-far.  P:" << P.str() << endl;
			best_MIP_score = this_MIP_score;
			MIPs.clear();
			MIPs.push_back( P );			
		}
	}
	
	ASSERT( ! MIPs.empty() );	
//	vec_sort( min_partitions );
	
	return MIPs;
}


int t_consciousness::save_highestbracketphis_spectrum( string ofilename )
// Returns all partitions with the same minimum normalized EI.
{

	
	// TODO: go over each subset, calculate the <phi> over that subset.  Then print it with the MIP.
	
//	// get all of the highest_bracketphis
//	vector<t_phi_result> highest_bracketphis = this->highest_bracketphis( );
//	assert( ! highest_bracketphis.empty() );
//	t_subset MC_subset = (highest_bracketphis[0]).subset;
	
	double ei, norm, norm_ei;
	//foreach line in the partitions file...
//	t_partition P( max_num_parts  );

	ofstream ofile;
	ofile.open( ofilename.c_str() );
	if( ! ofile.is_open() ) {
		cerr << "* Error! -- Could not open filename " << ofilename << endl;
		return 1;
	}
	
	// Always output some basic network statistics
	ofile << "# statistics:" << "\tN=" << numunits << "\tmax_num_parts=" << max_num_parts << "\t#x1states=" << NUM_X1_STATES;
	ofile << "\tH(X1)=" << H_X1 << "\tH(X0|X1)=" << H_X0_GIVEN_X1 << "\tnormalization=" << FLAG__MIP_METHOD;
	
	// print comments, if we have any
	if( ! this->network_comments.empty() )
	{
		for( int i=0; i<this->network_comments.size(); i++ ) {
			ofile << "# comment:\t" << this->network_comments[i] << endl;
		}	
	}

	ofile << "# subset \t ei \t phi \t MIPs \t mip_score" << endl;
//	ofile << "# Main Complex: " << MC_subset.str() << endl;
	

	
	//foreach subset...calculate the MIP and the PHI
	//we start from the FULL_MASK because we can often skip steps that way
	for( unsigned int subset=FULL_MASK; subset>=1; subset-- )
	{

		t_phi_result r = this->bracketphi( subset );
		
		if( r.MIPs.empty() ) {
			assert( use_total_partition == false );
			assert( r.subset.numnodes() == 1 );
			continue;
		}
		

		// sort the MIPs
		for( int i=0; i<r.MIPs.size(); i++ )
			r.MIPs[i].sort();
		
		// print the phi for this subset
//		ofile << r.subset.str();
		ofile << r.subset.brep();		
		ofile << "\t" << r.ei;
		ofile << "\tbphi=" << r.min_phi;
		ofile << "\t" << partitions2str( r.MIPs );
		ofile << "\t" << r.mip_score << endl;
		
	}
	
	
	ofile << "# DONE!" << endl;
	ofile.flush();	//flush it
	ofile.close();	//close it
	
	return 0;
}

int t_consciousness::save_bracketphi_spectrum( string ofilename )
// Returns all partitions with the same minimum normalized EI.
{
	t_subset S( numunits, FULL_MASK );
	
	ofstream ofile;
	ofile.open( ofilename.c_str() );
	if( ! ofile.is_open() ) {
		cerr << "* Error! -- Could not open filename " << ofilename << endl;
		return 1;
	}
	
	// Always output some basic network statistics
	ofile << "# statistics:" << "\tN=" << this->numunits << "\tmax_num_parts=" << this->max_num_parts << "\t#x1states=" << this->NUM_X1_STATES;
	ofile << "\tnormalization=" << FLAG__MIP_METHOD;
	ofile << "\tH(X0|X1)=" << this->H_X0_GIVEN_X1 << "\tH(X1)=" << this->bracket_ei( this->FULL_MASK ) << endl;

//	ofile << "# partition \t |P| \t <ei(P)> \t norm \t MIP_score" << endl;
//	ofile << "# partition \t |P| \t synergy \t holism \t MIP_score" << endl;	

	
	
	// print comments, if we have any
	if( ! this->network_comments.empty() )
	{
		for( int i=0; i<this->network_comments.size(); i++ ) {
			ofile << "# comment:\t" << this->network_comments[i] << endl;
		}	
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// BEGIN the partitions generator
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int numnodes = S.numnodes();
	PARTITIONS_GENERATOR.reset( numnodes );	

	// define a partition with the maximum number of parts
	t_partition P( this->max_num_parts, S );
	
	// initialize node_map to be as large as it could ever be, numnodes
	// Now to create the array that maps from result -> incoming partition
	unsigned int nodemap[numnodes];
	if( numnodes < this->numunits )
		S.array( nodemap );

	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// END partitions generator
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// now read the parts from GENERATOR...
	while( get_next_partition(&PARTITIONS_GENERATOR, &P, nodemap, numnodes) )
	{	
		
//		if( P.size() != 2 )
//			continue;

//		if( P.size(0) >= 4 )
//			continue;
		
/*
		bool go_forward=true;
		
		for( int i=0; i<P.size(); i++ )
			if( P.size(i) != 3 )
				go_forward = false;

		if( go_forward == false )
			continue;
*/		
		P.sort();
//		cout << "calculating ei..." << endl;
//		cout.flush();
		
		
		
		double ei=bracket_ei(P);
		double score=MIP_score(P);

//		ofile << P.BAREstr() << "\t|P|=" << P.size() << "\tei=" << ei << "\tMIPscore=" << score;
		ofile << P.BAREstr() << "\t|P|=" << P.size();
		ofile << "\t<ei(P)>=" << ei;
//		ofile << "\tMIPscore=" << score;
        
        ofile << "\t ei(x1/P): \t ";
        
        // print the perstate phis as well.
		for( map<unsigned int,unsigned int>::iterator it=this->existing_states.begin(); it != this->existing_states.end(); it++ )
        {
			unsigned int x1 = (*it).first;
            double this_x1_ei = this->ei( x1, P );
            string bstring = binary_repr( x1, this->numunits );
			ofile << bstring << "=" << this_x1_ei << " ";
        }

		ofile << endl;
		ofile.flush();
	}
	
	//	cerr << "Finished all partitions!" << endl;
	//	cerr.flush();
	
	
	//	cout << "subset=" << subset << "\tFinished bracket_MIPs()" << endl;
	//	cout.flush();
	
	ofile << "# DONE!" << endl;
	ofile.flush();	//flush it
	ofile.close();	//close it

	return 0;
}

bool t_consciousness::get_next_partition( PartitionEnumerator* restrict pgen, t_partition* restrict P, const unsigned int* restrict nodemap, unsigned int numnodes )
// Takes the enumerator object, converts it into a subset, and slides it into partition P
// returns false when there are no more partitions
{
	
	const unsigned int* result = NULL;
	unsigned int numparts = 0;
	
	while( 1 ) {
		result = pgen->nextPartition();
		numparts = result[0];
		
		// if numparts is zero, we've explored all of the partitions
		if( numparts == 0 )
			return false;

		if( numparts==1 && use_total_partition==false )
			continue;

		if( numparts <= this->max_num_parts )
			break;
	}
	
	
	if( use_total_partition == false )
		ASSERT( 2 <= numparts );

	// if the partition paves the ENTIRE NETWORK (numnodes == this->numunits), then we don't need to 
	// do any remapping.
	if( numnodes == this->numunits ) {
		//		cout << "numnodes == this->numunits!" << endl;
		P->update( numparts, result+1 );
		return true;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// remap result -> incoming_P
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int* restrict incoming_P = (unsigned int*) calloc(numparts, sizeof(unsigned int));
	remapNodes( incoming_P, result, nodemap, numnodes );
	
	// copy incoming_P into P; free incoming_P
	P->update( numparts, incoming_P );
	free( incoming_P );
	
	
//	cout << "P:" << P->BAREstr() << endl;
//	cout << "\t P:" << P->BAREstr() << endl;
	
	return true;	
}


vector<t_partition> t_consciousness::bracket_MIPs( unsigned int subset, bool debug ) {
	t_subset S( numunits, subset );
	return bracket_MIPs( S, debug );
}


double t_consciousness::MIP_score( const bitmask x1, const t_partition& restrict P )
// returns a MIP score for each partition.  Higher values are more MIP-like.
{
	double z=0.0;
    
	if( FLAG__MIP_METHOD == "NONE" ) {
		// for a MIP Z should be SMALL
		z = this->ei( x1, P ) * -1;
        //		z = hurt * -1;
	}
	
	else if( FLAG__MIP_METHOD == "TONONI" ) {
		// for a MIP Z should be SMALL
		double hurt = this->normalized_ei( x1, P );
		z = hurt * -1;
	}

	else if( FLAG__MIP_METHOD == "ATOMIC" ) {
		// the MIP score becomes equal the number of parts in the partition.
		z = (double) P.size();
	}
    
    
	else { 
		cerr << "Error. Did not know which MIP method you speak of -- '" << FLAG__MIP_METHOD << "'" << endl;
		exit(-1);
	}
	if( fequals(z, 0.0 ) )
		z = 0.0;
	
    //	assert( 0.0 <= z );
	
    //	cerr << "Calculated a MIP score=" << z << endl;
    //	cerr.flush();
	
	return z;    
}
    
double t_consciousness::MIP_score( const t_partition& restrict P )
// returns a MIP score for each partition.  Higher values are more MIP-like.
{
	double z=0.0;

	if( FLAG__MIP_METHOD == "NONE" ) {
		// for a MIP Z should be SMALL
		z = this->bracket_ei( P ) * -1;
//		z = hurt * -1;
	}
    
	
	else if( FLAG__MIP_METHOD == "TONONI" ) {
		// for a MIP Z should be SMALL
		double hurt = this->bracket_normalized_ei( P );
		z = hurt * -1;
	}
	
	else if( FLAG__MIP_METHOD == "ATOMIC" ) {
		// the MIP score becomes equal the number of parts in the partition.
		z = (double) P.size();
	}
		   
	else if( FLAG__MIP_METHOD == "TONONI_PROD_M0" )
	{
		// for a MIP Z should be SMALL		
		double help = this->bracket_ei( P );
		unsigned int hurt = 1;
		
		// create the product of parts
		for( int i=0; i<P.size(); i++ )
			hurt *= P.size(i);
		
		z = (help / (double) hurt) * -1;
	}

	
	
	else if( FLAG__MIP_METHOD == "TONONI_PROD_M1" )
	{
		// for a MIP Z should be SMALL				
		double help = this->bracket_ei( P );
		double hurt = 1.0;
		double temp = DBL_NONE;
		
		// create the product of parts
		for( int i=0; i<P.size(); i++ ) {
			
			temp = H_M1( P[i] );
			// if temp != 0.0, don't multiply.
			
			if( ! fequals(temp, 0.0 ) )
				hurt *= temp;
		}

		z = (help/hurt) * -1;
	}
	

	else { 
		cerr << "Error. Did not know which MIP method you speak of -- '" << FLAG__MIP_METHOD << "'" << endl;
		exit(-1);
	}
	if( fequals(z, 0.0 ) )
		z = 0.0;
	
//	assert( 0.0 <= z );
	
//	cerr << "Calculated a MIP score=" << z << endl;
//	cerr.flush();
	
	return z;
}

vector<t_partition> t_consciousness::bracket_MIPs( const t_subset& restrict S, bool debug )
// Returns all partitions with the same minimum normalized EI.
{
	ASSERT( S > 0 );
//	cerr << "MIPs for " << S.str() << "..." << endl;
//	cerr.flush();

	unsigned int numnodes = S.numnodes();
	
	// special value
	double best_MIP_score=DBL_MAX * -1;

	//foreach line in the partitions file...
	vector<t_partition> MIPs;	
	MIPs.clear();

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// BEGIN the partitions generator
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	PARTITIONS_GENERATOR.reset( numnodes );
	

	// define a partition with the maximum number of parts
	t_partition P( max_num_parts, S );
	
	// initialize node_map to be as large as it could ever be, numnodes
	// Now to create the array that maps from result -> incoming partition
	unsigned int nodemap[numnodes];
	if( numnodes < numunits )
		S.array( nodemap );
	

	unsigned int count=0;
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// END partitions generator
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// now read the parts from GENERATOR...
	while( get_next_partition(&PARTITIONS_GENERATOR, &P, nodemap, numnodes) )	
	{		
		double this_MIP_score = MIP_score( P );
		assert( -1*DBL_MAX < this_MIP_score );
		
		if( debug )
			cerr << "P: " << P.str() << "\tMIPscore=" << this_MIP_score << endl;
		//		cout << "subset=" << subset << "\t ei=" << ei << " \t norm=" << norm << endl;
		

		// if this norm EI is EQUAL TO the current one, do...
		if( fequals(this_MIP_score, best_MIP_score) )
		{
			assert( !MIPs.empty() );
			
			if( MIPs.size() < MAXIMUM_NUMBER_OF_MIPS_TO_STORE )
				MIPs.push_back( P );
			else
			{
				for( vector<t_partition>::iterator it=MIPs.begin(); it != MIPs.end(); ++it )
				{
					// if P wins a tiebreaker, replace *it with P
					if( break_MIPscore_ties( *it, P ) )
					{
						// if the current P has MORE PARTS than a MIP, replace that MIP, and stop the loop
						MIPs.erase( it );
						MIPs.push_back( P );
						it = MIPs.end();
						break;
					}
				}
			}
		}
		
		// if this partition has a HIGHER MIP score, clear previous MIPs and set partition as the MIP.
		else if( this_MIP_score > best_MIP_score ) {
			best_MIP_score = this_MIP_score;
			MIPs.clear();
			MIPs.push_back( P );			
		}
		
//		if( debug )
//		count += 1;
	}


	//foreach MIP...
//	foreach( t_partition& MIP, MIPs ) {
//		cerr << "within:MIP=" << MIP.str() << "\t\tscore=" << MIP_score( MIP ) << endl;
//	}
	
//	cout << "bracketMIPs: best_score=" << best_MIP_score << endl;
		
	
	// sort the parts
	foreach( t_partition& MIP, MIPs )
		MIP.sort();
	
//	// sort the MIPs
	vec_sort( MIPs );
	
	if( debug ) {
		cout << "\t Tried " << count << " partitions for subset=" << S.binrep() << endl;
		cout << "\t Skipping Total Partition? " << PARTITIONS_GENERATOR.skip_total_partition << endl;	
		cout << "------------------" << endl;
		cout.flush();
	}
	

	//foreach MIP...
//	foreach( t_partition& MIP, MIPs ) {
//		cerr << "within2:MIP=" << MIP.str() << "\t\tscore=" << MIP_score( MIP ) << endl;
//	}
//	cout << "bracketMIPs: best_score=" << best_MIP_score << endl;
	

	if( use_total_partition )
		assert( ! MIPs.empty() );
	
	return MIPs;
}


t_phi_result t_consciousness::phi( const unsigned int x1, bitmask inp_subset )
/* Returns phi(x1) = ei( x1/P=MIP(x1) ) for state x1 searching over all partitions in partitions_filename */
{    
	//verify our inputs
	ASSERT( this->existing_states.count(x1) );
	ASSERT( inp_subset <= this->FULL_MASK );

	
	if( inp_subset==0 )
		inp_subset = this->FULL_MASK;

	const bitmask subset = inp_subset;

	// assert x1 is a subset of subset
	assert( x1 | (subset == subset) );
	
	
	t_phi_result z;
	
	// set the subset and x1
	z.subset = subset;
	z.x1 = x1;
	
	// set the MIPs
	z.MIPs = this->MIPs( x1, subset );
	
//	// set the constant normalized phi
//	z.phi_norm = this->ei( x1, z.MIPs[0], subset ) / this->normalization( z.MIPs[0] );
	
	z.ei = this->ei( x1, subset );
	
	// set z.phi to the minimum unnormalized effective information
	double ei, min_ei = DBL_MAX;
	
	(z.PHIs).clear();
	// go through each MIP, record the ei, and find the lowest unnormalized ei
	for( int i=0; i<z.MIPs.size(); i++ ) {
		ei = this->ei( x1, z.MIPs[i] );

		// append this ei to the PHIs
		(z.PHIs).push_back(ei);
		
		if( ei < min_ei )
			min_ei = ei;	
	}
	z.min_phi = min_ei;

	if( fequals(z.min_phi,z.ei) )
		z.min_phi = z.ei;

	
	// z.min_phi MAY BE GREATER than z.ei
	
	assert( 0.0 <= z.min_phi );
    
    return z;
}


double t_consciousness::ei( const unsigned int x1, const t_partition& P, const string force_version )
/* Returns ei(x1/P) for state x1 and partition P */
{
	const unsigned int subset = P.subset;
	
	assert( subset <= FULL_MASK );
	assert( (x1|subset) == subset );

    if( force_version != "UNITS" && force_version != "WIRES" )
        cout << "force_version=" << force_version << endl;
    
	assert( force_version == "UNITS" || force_version == "WIRES" );
	
	// is this the total partition?  If so, then ei(x1/P) = ei(x1)
	if( P.size() == 1 )
        return this->ei(x1, subset);
	
	//	cout << "trying masks=" << parts2str( *parts ) << endl;
	
	const double h_s0_given_s1 = H_S0_given_S1_equals_s1(x1, subset);
//	double ei = 1.0 / pow(2,h_s0_given_s1);
	double ei = 0.0;
	
	// the sum entropy of parts is the same regardless of the subset
	// for each part in parts...
	double entropy_of_parts=0.0;
	for(int i=0; i<P.size(); i++)
	{
		if( force_version == "WIRES" ) {
			entropy_of_parts += entropy_of_part_given_x1__WIRES( x1, P[i], subset );
		}
		else if( force_version == "UNITS" ) {
			entropy_of_parts += entropy_of_part_given_x1( x1, P[i], subset );
		}
		else if( FLAG__STATE_EI == "ANTI-CONDITIONED" || FLAG__STATE_EI == "WIRES" ) {
			entropy_of_parts += entropy_of_part_given_x1__WIRES( x1, P[i], subset );
		}
		
		else if( FLAG__STATE_EI == "NAKED" || FLAG__STATE_EI == "UNITS" ) {
			entropy_of_parts += entropy_of_part_given_x1( x1, P[i], subset );	
		}
		else {
			assert( 0 == 1 );
		}
			
	}
	
//	cout << "entropy_of_parts=" << entropy_of_parts << endl;
	ei = entropy_of_parts;
	ei -= h_s0_given_s1;

	assert( entropy_of_parts == entropy_of_parts );
	// assert that ei is not NaN
	assert( ei == ei );
	
	// remove any -0.0's	
	if( fequals(ei,0.0) )
		ei = 0.0;

	if( ei < 0.0 ) {
        cout << "state=" << x1 << endl;
		cout << "entropy_of_parts=" << entropy_of_parts << endl;
		cout << "H[S0|S1=s1]=" << h_s0_given_s1 << endl;
		cout << "ei=" << ei << endl;
        cout << endl << endl;
	}
	// check bounds on ei
	assert( 0.0 <= ei );

    return ei;
}

unsigned int t_consciousness::num_s1_states( const bitmask subset ) {
	t_subset S( numunits, subset );
	return num_s1_states( S );
}

unsigned int t_consciousness::num_s1_states( const t_subset S )
// returns the number of distinct s1 states in subset S
{
	vector<unsigned int> seen_s1_states;
	
	// for each unique x1 state...
	for( map<unsigned int,unsigned int>::iterator it=this->existing_states.begin(); it != this->existing_states.end(); it++ )
	{
		unsigned int x1 = (*it).first;
		unsigned int s1 = x1 & S.S;
		
		if( ! vec_in( s1, seen_s1_states ) )
			seen_s1_states.push_back( s1 );
	}

	return seen_s1_states.size();

}

unsigned int t_consciousness::numunits_in_subset( unsigned int mask, const bool allow_zero ) { return numunits_in_mask(mask, allow_zero); }

unsigned int t_consciousness::numunits_in_mask( unsigned int subset, const bool allow_zero )
// Returns the number of 1s in the subset, representing the number of units
{
	static unsigned int old_subset=0, old_z=0;

	// if this is the subset we just did, use the old value.
	if( subset == old_subset )
		return old_z;
	
    old_subset = subset;
	unsigned int z=0;	

	do	{
		if( (subset&1) )
			z += 1;
		
		subset >>= 1;
	}while( subset );
	

//	ASSERT( 0 < z || allow_zero );
	ASSERT( z <= numunits );

	old_z = z;
    
//    cout << "\t subset=" << subset << "\tz=" << z << endl;

	return z;
}

double t_consciousness::ei( const t_state& s1 ) {
    return this->ei( s1.value, s1.mask );
}

double t_consciousness::ei( const unsigned int s1, const bitmask subset )
/* Returns ei(x1, P=total partition) = numunits - H[X0|X1=x1]
 
 This function is called when ei(x1/P) is passed a total partition.
 
 If sunset is true, calculates H[S0] - H[S0|S1=s1].
 */
{
//	if( subset == 0 )
//		subset = FULL_MASK;
	assert( (s1|subset) == subset );
	assert( subset <= FULL_MASK );

	double ei = H_M0( subset ) - H_S0_given_S1_equals_s1(s1,subset);

	// remove -0.0's
	if( fequals(ei,0.0) )
		ei = 0.0;
	
	
	// check that ei is not NaN
	assert( ei == ei );
	
	// check bounds on ei
	assert( 0.0 <= ei );

	return ei;
}

double t_consciousness::average_phi( unsigned int subset )
/* This function computes the average phi by doing: 
 
 average_phi = SUM_x1{ p(x1) * phi(x1) }
 
 Calculating phi(x1) searching over all partitions.
 
 */
{
    // foreach x1...
    
    double sum = 0.0, sum_prob=0.0;
	unsigned int x1;
    double prob_this_x1=0.0;
    t_phi_result res;
	
    // Foreach output state x1...
    for( map<unsigned int,unsigned int>::iterator it=this->existing_states.begin(); it != this->existing_states.end(); it++ ) {
        x1 = (*it).first;
        prob_this_x1 = double((*it).second) / (double) this->numstates;
        res = this->phi( x1, subset );
		
		//        cout << "\tx1=" << x1 << "\tprob(X1=" << x1 << ")=" << prob_this_x1 << "\t phi(x1)=" << this_phi << "\t product=" << prob_this_x1 * this_phi << endl;
        sum_prob += prob_this_x1;
        sum += (prob_this_x1 * res.min_phi);
		//        cout << "sum=" << sum << endl;
	}
	
    // sanity check.
	//    cout << "sum_prob=" << sum_prob << endl;
    assert( sum_prob == 1.0 );
	
    return sum;
}

// alias to bracketphi() that uses t_subset
t_phi_result t_consciousness::bracketphi( unsigned int subset, bool require_atomic_partition )
{
	t_subset S( numunits, subset );
	return bracketphi( S, require_atomic_partition );
}


t_partition t_consciousness::atomic_partition( unsigned int subset )
// returns the atomic partition for a specified subset
{
	// if no subset is specified, we mean the full system.
	if( ! subset )
		subset = FULL_MASK;

//	cout << "FULL_MASK=" << FULL_MASK << endl;
	
	// define a partition with the maximum number of parts
	t_partition P( this->max_num_parts, FULL_MASK );

//	cout << "P.subset=" << P.subset << "\t\t subset=" << subset << endl;
	
	unsigned int subset_size = numunits_in_mask( subset );
	unsigned int* restrict nodes = new unsigned int[subset_size];
	for( int i=0; i<subset_size; i++ ) 
		nodes[i] = 1 << i;
		
//	cout << "FULL_MASK=" << FULL_MASK << endl;
//	cout << "P.correct_subset=" << P.correct_subset << "\t\tP.subset=" << P.subset << "\t\t subset=" << subset << endl;
	
/*	
	for( int i=0; i<subset_size; i++ )
		cout << " " << nodes[i];
	cout << endl;
*/	
	P.update( subset_size, nodes );

//	cout << "P.correct_subset=" << P.correct_subset << "\t\tP.subset=" << P.subset << "\t\t subset=" << subset << endl;
	
	delete [] nodes;

	return P;
}

t_phi_result t_consciousness::bracketphi( const t_subset& restrict S, bool require_atomic_partition )
// bracketphi is defined as the bracket_ei over the <MIP>.
{
	t_phi_result z;	
	z.subset = S;

//	cout << "z.subset=" << z.subset.S << "\tbinrep=" << z.subset.binrep() << "\tmax_numnodes=" << z.subset.max_numnodes << "\tmax_S=" << z.subset.max_S << endl; cout.flush();
	
	z.ei = bracket_ei( z.subset );
	
	//	cout << "\rGetting the MIPs...";
	//	cout.flush();
	// set the MIPs
	z.MIPs.clear();
	
	if( require_atomic_partition )
		z.MIPs.push_back( atomic_partition(z.subset.S) );
	else
		z.MIPs = bracket_MIPs( z.subset );

//	cout << "[3] bracketphi: " << "S=" << S.S << "\t\t" << S.binrep() << "(" << S.str() << ")" << endl; cout.flush();	
	
	// if there are no MIPs, require that use_total_partion is FALSE and there is only 1 node.
	if( z.MIPs.empty() ) {
		assert( use_total_partition == false );
		assert( z.subset.numnodes() == 1 );
		z.ei = 0.0;
		z.min_phi = 0.0;
//		z.phi_norm = 0.0;
		z.mip_score = 0.0;
		return z;
	}
	
	
	// if using total partition, we should NEVER get empty MIPs
	if( use_total_partition )
		assert( ! z.MIPs.empty() );
	
	// calculate the AVERAGE phi_norm
//	z.phi_norm = DBL_NONE;
	z.mip_score = 0.0;

//	cerr << "MIPs for " << S.str() << endl;
	// average the MIP score and normalized phi over all of the MIPs
	foreach( t_partition MIP, z.MIPs ) {
		double this_mipscore = this->MIP_score( MIP );
		z.mip_score += this_mipscore;
		
//		cerr << "MIP=" << MIP.BAREstr() << "\t\t score=" << this_mipscore << endl;
//		cerr.flush();
	}
	
	z.mip_score /= (double) z.MIPs.size();
	
//	cout << "[4] bracketphi: " << "S=" << S.S << "\t\t" << S.binrep() << "(" << S.str() << ")" << endl; cout.flush();		
	
	// sanity check that all have the same MIP score.
	foreach( t_partition MIP, z.MIPs ) {
		double this_mip_score = this->MIP_score( MIP );
		if( ! fequals( this_mip_score, z.mip_score ) ) { 
			cout << "this_mip_score=" << this_mip_score << endl;
			cout << "diff=" << this_mip_score - z.mip_score << endl;
		}
		assert( fequals(this_mip_score, z.mip_score) );
	}

	// set z.phi to the minimum unnormalized effective information
	double min_ei=DBL_MAX;
	
	z.PHIs.clear();
	// go through each MIP and find the lowest unnormalized ei
	foreach( t_partition MIP, z.MIPs ) {

//		cout << "[MIP] bracketphi: " << "S=" << S.S << "\t\t" << S.binrep() << "(" << S.str() << ") MIP=" << MIP.BAREstr() << endl; cout.flush();				
		double ei = bracket_ei( MIP );
//		double vsi = vsi( MIP );
		
//		cout << "[MIP-2] bracketphi: " << "S=" << S.S << "\t\t" << S.binrep() << "(" << S.str() << ") MIP=" << MIP.BAREstr() << endl; cout.flush();				
		// add this ei to our list of PHIs
		(z.PHIs).push_back(ei);

		// make the min_ei
		min_ei = min( min_ei, ei );
//		z.min_vsi = min( vsi, z.min_vsi );

	}
	z.min_phi = min( min_ei, z.ei );

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	// IT IS VERY VERY IMPORTANT WE DO NOT ROUND THE PHI_NORM.  THE BRACKET_EI IS DECIMAL_PRECISION and NORMALIZATION is DECIMAL_PRECISION
	// WE WILL GET ZERO FOR THE PHI_NORM SOMETIMES if we round that to decimal precision too!
	// z.phi_norm = XXXXXX /* myround( z.phi_norm ); */ XXXXXXXXX
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////		


	for( int i=0; i<z.MIPs.size(); i++ ) {	
		unsigned int mip_size = z.MIPs[i].size();

		// check that if there is a total partition, it is equal to the bracket_ei of the entire subset.
		if( FLAG__MIP_METHOD == "TONONI" && mip_size == 1 ) {
			assert( fequals(z.PHIs[i], this->bracket_ei(z.subset)) );
//			cout << "[check] bracketphi: " << "z.subset=" << z.subset.str() << "\t\t" << z.subset.binrep() << "(" << z.subset.str() << ") MIP=" << z.MIPs[i].BAREstr() << endl; cout.flush();
		}
	}
	
//	if( ! (z.min_phi <= this->bracket_ei(subset)) ) {
//		cerr << "subset=" << subset << "\t min_phi=" << z.min_phi << " \t bracket_ei=" << this->bracket_ei(subset) << endl; 
//	}
	
	ASSERT( z.min_phi <= z.ei || fequals(z.min_phi, z.ei) );
	
	return z;
}

double t_consciousness::bracket_ei( unsigned int subset )
// Returns \overline{ei}(P=total partition) = numunits - H[X0|X1] -- This function only exists for convenience so that humans can call it.
// if passed subset it returns |S| - H[S0|S1]
{
	ASSERT( subset > 0 );
	double z;	
	
	// is S the total partition?
	if( subset == FULL_MASK )
	{
		
		//		if( ! fequals( H_X0_GIVEN_X1, H_S0_given_S1(S) ) ) {
		//			cout << "H[X0|X1]=" << H_X0_GIVEN_X1 << "\t\t H[S0|S1]=" << H_S0_given_S1(S) << endl;
		//			cout << "H[X1]=" << this->H_X1 << endl;
		//		}
		
		assert( fequals( H_X0_GIVEN_X1, H_S0_given_S1(subset) ) );
		
		z = this->H_X1;
	}
	
	
    z = I_A0_B1( subset, subset );
    
	if( subset != FULL_MASK )
		cerr << "Computing <ei(subset)> when subset != FULL_MASK" << endl;
    
	// replace any -0.0's with 0.0
	if( z == -0.0 )
		z = fabs(z);
	
	
	return z;
}

// calculates bracket_ei for a subset 
double t_consciousness::bracket_ei( const t_subset& restrict S ) { return bracket_ei( S.S ); }


double t_consciousness::bracket_ei( const t_partition& restrict P )
// Returns \overline{ei}(P) for partition P
{
	assert( P.subset > 0 );
	
	// is this the total partition?  If so, then average ei(P) is simply bracket_ei()
	if( P.size() == 1 ) {
        assert( P.subset == FULL_MASK );
		return bracket_ei( P.subset );
    }

	assert( 2 <= P.size() && P.size() <= numunits );

	// keep these default values as they are.  It's important.
    double z = I_A0_B1( P.subset, P.subset );
    
    for( int i=0; i<P.size(); i++ ) {
        assertmsg( (P[i] | P.subset) == P.subset, "part was not a subset of the partition.subset" );
        z -= I_A0_B1( P[i], P[i] );
    }

	// replace any -0.0's with 0.0
	if( fequals(z,0.0) )
		z = 0.0;
	
	ASSERT( z == z );
	
    return z;
}


double t_consciousness::H_M0_GIVEN_M1( const unsigned int part_mask )
// returns H[M_0|M_1] for a given part specified by part_mask
// if anticond_mask is defined
{

	// if the cache_part_entropies is not NONE, then it's a valid cache
	if( FLAG__CACHE_M0_GIVEN_M1 && this->H_M0_GIVEN_M1_cache[part_mask] >= 0.0 )
		return this->H_M0_GIVEN_M1_cache[ part_mask ];

	// if doing the FULL_MASK, the answer is simply H[X0|X1]
	if( part_mask == FULL_MASK )
		return H_X0_GIVEN_X1;
	
	ASSERT( part_mask > 0 );
	ASSERT( part_mask < FULL_MASK );
	
	double z=0.0;
    z = H_M0_GIVEN_M1__ELEMENTS( part_mask );

	
	assert( z >= 0.0 );

	// remove any -0.0s
	z = fabs(z);
	
	if( FLAG__CACHE_M0_GIVEN_M1 ) {
		// Store the part entropy in the cache, if we're doing that.
		ASSERT( this->H_M0_GIVEN_M1_cache[part_mask] != this->H_M0_GIVEN_M1_cache[part_mask] );
		
		// set the part entropy
		this->H_M0_GIVEN_M1_cache[ part_mask ] = z;
	}
	
	return z;
}

void t_consciousness::load_transformations( t_cmdline_args &args, bool show )
{
	// clear the MIP caches
	this->network_comments.clear();

	FILE* f=fopen( args.networkfilename.c_str() ,"r+t");
	
	if( f == NULL ) {
		cerr << "Could not open filename '" << args.networkfilename << "'.  Exiting." << endl;
		perror( args.networkfilename.c_str() );
		exit(1);
	}
	
	fpos_t last_location;
	fgetpos(f, &last_location);	// get the location at the very beginning of the file
	
	// Go until we find a line that isn't blank and doesn't begin with a '#'
	char line[1024];
	fgets(line, sizeof(line), f );
	string stringline = string(line);
	
	//	cout << "posititon=" << last_location << endl;
	//	cout.flush();
	
	while( ! str_startswith(stringline, "type=" ) )
	{
		if( str_startswith(stringline, "comment=") || str_startswith(stringline, "comment:" ) )
		{
			
			// remove the text "comment:" from stringline
			stringline.erase(0,8);
			
			// remove any newline chars at the end
			stringline = str_rstrip(stringline, "\r\n" );
			stringline = str_strip( stringline, " ");
			stringline = str_strip( stringline, "\t");
			//			cout << "pushing back comment='" << stringline << "'" << endl;
			this->network_comments.push_back( stringline );			
		}
		
		fgetpos(f, &last_location);
		fgets(line, sizeof(line), f );
		stringline = string(line);
	}
	
	fsetpos(f, &last_location);
	
	// What type is this network?	
	char type[300];
	fscanf(f,"type=%s\n", type);
	this->networktype = tolower(string(type));
	
	assert( this->networktype == "neural" || this->networktype == "circuit" || networktype == "transition" );
	
	if(show)
		cout << "type=" << this->networktype << endl;
	
	//READIN NEURAL NETWORK
	if( this->networktype == "transition" )
	{
		fgets(line, sizeof(line), f );
		stringline = str_rstrip( string(line), "\r\n" );
		
		// Skip lines until find a line beginning with "nodes" -- but be sure to add any comments.
		while( ! str_istartswith(stringline, "nodes=" ) )
		{
			if( str_istartswith(stringline, "comment=") || str_istartswith(stringline, "comment:" ) )
			{
				// remove the text "comment:" from stringline
				stringline.erase(0,8);
				
				this->network_comments.push_back( stringline );			
			}

			fgets(line, sizeof(line), f );
			stringline = str_rstrip(string(line), "\r\n");			
		}

		string nodes_str = str_istrip( stringline, "nodes=" );
		this-> numunits = atoi( nodes_str.c_str() );

		if(show)
			cout << "numunits=" << this->numunits << endl;
		
	}
	else if( this->networktype == "neural" )
	{
		
		fgets(line, sizeof(line), f );
		stringline = str_rstrip( string(line), "\r\n" );
		
		// Skip lines until find a line beginning with "nodes" -- but be sure to add any comments.
		while( ! str_istartswith(stringline, "nodes=" ) )
		{
			if( str_istartswith(stringline, "comment=") || str_istartswith(stringline, "comment:" ) )
			{
				// remove the text "comment:" from stringline
				stringline.erase(0,8);
				
				this->network_comments.push_back( stringline );			
			}
			
			// if we're setting fire at threshold
			else if( str_istartswith(stringline, "fire_only_at_threshold=") )
			{
				string fire_only_at_threshold = tolower( str_istrip( stringline, "fire_only_at_threshold=") );
				
				if( fire_only_at_threshold == "true" || fire_only_at_threshold == "1" )
					this->neuron_fires_only_at_threshold = true;
				
				else if( fire_only_at_threshold == "false" || fire_only_at_threshold == "0" )
					this->neuron_fires_only_at_threshold = false;
				
				else {
					cerr << "Had an error setting FLAG__NEURON_FIRES_ONLY_AT_THRESHOLD=" << fire_only_at_threshold << endl;
					exit(-1);
				}
				
				assert( fire_only_at_threshold == "true" || fire_only_at_threshold == "1" || fire_only_at_threshold == "false" || fire_only_at_threshold == "0" );
				
			}
			
			fgets(line, sizeof(line), f );
			stringline = str_rstrip(string(line), "\r\n");
		}
		
		
		string nodes_str = str_istrip( stringline, "nodes=" );
		unsigned int numunits = atoi( nodes_str.c_str() );
		
		if(show)
			cout << "numunits=" << numunits << endl;
		
		assert( numunits >= 1 );
		this->node_thresholds.resize(numunits);
		
		
		double temp=0.0;
		fscanf(f,"thresholds=");
		if(show)
			cout << "thresholds\t=";
		
		for( int i=0; i<numunits; i++ ) {
			fscanf(f,"%lf",&temp);
			assert( temp >= 0.0 );
			this->node_thresholds[i] = temp;
			if(show)
				printf("%.3f ",this->node_thresholds[i]);
		}
		if(show)
			cout << endl;
		
		//Readin and print the connection matrix
		fscanf(f,"\n");
		this->weightmatrix.resize(numunits);	//define n rows
		if(show)
			cout << "cxnmatrix=" << endl;
		
		for( int i=0; i<numunits; i++ ) {
			(this->weightmatrix[i]).resize(numunits);	//define n columns
			
			for( int j=0; j<numunits; j++ ) {
				fscanf(f,"%lf",&this->weightmatrix[i][j]);
				
				if(show)
					printf("%+.3lf ", this->weightmatrix[i][j]);
				//				cout << this->weightmatrix[i][j] << " ";
			}
			fscanf(f,"\n");
			if(show) cout << endl;
		}
		
		// define numunits and numstates for cleaner code later on.
		this->numunits = numunits;		
	}
	
	//READIN CIRCUIT NETWORK
	else if( this->networktype == "circuit" )
	{
		t_transformation t;
		int I,A,B;
		char op;
		while(!(feof(f))){
			fscanf(f,"%i	=	%i	%c	%i\n",&I,&A,&op,&B);
			//			t.A=(unsigned char)A;
			//			t.B=(unsigned char)B;
			t.A=A;
			t.B=B;
			switch(op){
				case 'A':t.operation=0; break;
				case 'O':t.operation=1; break;
				case 'X':t.operation=2; break;
				case 'E':t.operation=3; break;
				case 'N':t.operation=4; break;
			}
			transformations.push_back(t);
		}
		//		show_rules();
		
		// define numunits and numstates for cleaner code later on.
		this->numunits = this->transformations.size();
		
	}
	fclose(f);
	
	assert( this->numunits <= 31 );
	this->numstates = 1 << this->numunits;	//numstates = 2**numunits
	
	
	if( this->neuron_fires_only_at_threshold )
		cout << "NEURON_FIRES_ONLY_AT_THRESHOLD=true" << endl;
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// if max_num_parts was defined, set it.  By default, set it to the number of nodes in the network (all parts)
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if( 1 <= args.max_num_parts && args.max_num_parts <= this->numunits ) 
		this->max_num_parts = args.max_num_parts;
	else
		this->max_num_parts = this->numunits;		//default
	
	args.max_num_parts = this->max_num_parts;

	cout << "numunits=" << this->numunits << "\t max_num_parts=" << this->max_num_parts << endl;
	
	assert(this->max_num_parts >= 1);
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Sanity checks
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	assert( 0 <= args.worker_index && args.worker_index <= args.numworkers-1 );
	assert( 0 <= args.numworkers );
	

}		


//makes a list of all X0 to X1
void t_consciousness::make_all_states( t_cmdline_args &args, bool show){
	
	unsigned int n=this->numunits, k=this->numstates;
//	if(show) cout << "Performing " << k << " different x0 states..." << endl;
	this->existing_states.clear();
		
	assert( this->states == NULL );
	this->states = new unsigned int[k];

	assert( this->numunits <= 31 );
	
	if( show )
		cout << "networktype = '" << this->networktype << "'" << endl;
	
	if( this->networktype == "transition" )
	{	
		
		// reopen the networkfilename, and move ahead until the transitions begin
		FILE* f=fopen( args.networkfilename.c_str() ,"r+t");
		assert( f != NULL );
		
		char line[1024];
		fgets(line, sizeof(line), f );
		string stringline = string(line);
		fpos_t last_location;
		

		// go through the file until find a line containing a " -> "
		while( stringline.find(" -> ") == string::npos )
		{

			fgetpos(f, &last_location);
			fgets(line, sizeof(line), f );
			stringline = string(line);
		}

		// set the location in the file to the first line with transitions
		fsetpos(f, &last_location);
		
		
		// set all entries in states to UINT_NONE
		for( int x0=0; x0<numstates; x0++ )
			states[x0] = UINT_NONE;
		
        unsigned int old_x0=0;
		// now read until end of file getting the transition table
		while(!(feof(f))){
			char chars_x0[numunits], chars_x1[numunits];

			// 1. Read in both as character-arrays.
			fscanf(f,"%s -> %s\n", chars_x0, chars_x1);

//			string str_x0 = string(chars_x0);
//			string str_x1 = string(chars_x1);
			
			// 2. If the lengths are EQUAL to NUMUNITS, treat as binary.			
			// 3. Otherwise, treat as unsigned integers.
			//unsigned int x0=UINT_NONE, x1=UINT_NONE;
			unsigned int x0=UINT_NONE, x1=UINT_NONE;
			
			( strlen(chars_x0) == numunits ) ? ( x0=bin2dec(chars_x0) ) : ( x0=atoi(chars_x0) );
			( strlen(chars_x1) == numunits ) ? ( x1=bin2dec(chars_x1) ) : ( x1=atoi(chars_x1) );

			// check bounds on x0 and x1
			assert( x0 < numstates );
			assert( x1 < numstates );
			
            assert( old_x0 <= x0 );
            
			// check that we haven't done this x0 before
			assert( states[x0] == UINT_NONE );
			
            // we set this to ensure all of the x0 states are IN ORDER.
            // (technically this isn't required, but lets be clean about it.)
            old_x0 = x0;
            
			this->states[x0] = x1;
			this->existing_states[x1]++;

		}

		// now assert that no UINT_NONE states exist.
		for( int x0=0; x0<numstates; x0++ )
        {
            if( states[x0] == UINT_NONE )
            {
                cerr << "Not output state for x0=" << x0 << " defined!" << endl;
                cerr << "x0=" << binary_repr(x0) << endl;
                cerr.flush();
            }
			assert( states[x0] != UINT_NONE );
        }

	}
	
	
	//Make the states for a circuit network
	else if( this->networktype == "circuit" )
	{
		int x1, p;
		
/*		We don't need the state_map (for now)		
		this->state_map.resize(k);
		for(i=0;i<k;i++)
		{
			this->state_map[i].resize(k);
			for(j=0;j<k;j++)
				this->state_map[i][j]=0;
		}
*/		
		for(int x0=0;x0<k;x0++)
		{
			x1=0;
			for(int j=0;j<n;j++)
			{
				p=0;
				switch(transformations[j].operation){
					case 0://AND
						p=((x0>>transformations[j].A)&1)&&((x0>>transformations[j].B)&1);
						break;
					case 1://OR
						p=((x0>>transformations[j].A)&1)||((x0>>transformations[j].B)&1);
						break;
					case 2://XOR
						p=((x0>>transformations[j].A)&1)^((x0>>transformations[j].B)&1);
						break;
					case 3://EQUAL
						p=(x0>>transformations[j].A)&1;
						break;
					case 4://NOT
						p=((x0>>transformations[j].A)&1);
						break;
					default:
						cerr << "Error: Did not understand specified transformation!";
						exit(-1);
				}
				x1 += (p<<j);
			}
			this->states[x0] = x1;
			this->existing_states[x1]++;
//			this->possible_x0s[x1].push_back( x0 );
			//			printf("%i	->	%i\n",i,ns);
//			this->state_map[i][ns]++;
		}
	}
	
	// Make the states for a neural network
	else if( this->networktype == "neural" )
	{
		
/*		We don't need the state_map (for now)
		this->state_map.resize(k);

		for(int x0=0;x0<k;x0++)
		{
			this->state_map[x0].resize(k);
			for(int x1=0;x1<k;x1++)
				this->state_map[x0][x1]=0;
		}
*/
		//foreach input state x0...
		for(unsigned int x0=0;x0<k;x0++)
		{
			//			cout << "running state x0=" << x0 << " \t k=" << k << endl;
			
			unsigned int x1=0;	// the x1 we get for state x0
			
			//foreach node of x1...
			for(int j=0; j<n; j++) {
				double m=0.0;				
				
				for(unsigned int i=0; i<n; i++) {
					unsigned int x0_bit = ( x0 >> i ) & 1;
					m += x0_bit * this->weightmatrix[i][j];
				}
				
				if( (m == this->node_thresholds[j]) || (m > this->node_thresholds[j] && this->neuron_fires_only_at_threshold==false) ) {
					x1 |= (1 << j);		//	x1 += pow(2,j);
										//	cout << "\t" << j << " exceeded threshold...  x1=" << x1 << endl;
				}
				
			}
			
			// Have our x1 state now.  Update stuff.
			this->states[x0] = x1;
			this->existing_states[x1]++;
//			this->state_map[x0][x1] = 1;			
//			this->possible_x0s[x1].push_back( x0 );
//			if(show)
//				printf("%i	->	%i\n",x0,x1);
		}
	}

	this->NUM_X1_STATES = this->existing_states.size();
	this->x1_states = (unsigned int*) calloc(NUM_X1_STATES, sizeof(unsigned int));
	this->num_x0s_to_this_x1 = (unsigned int*) calloc(NUM_X1_STATES, sizeof(unsigned int));	
	
	assert( this->x1_states != NULL );
	assert( this->num_x0s_to_this_x1 != NULL );

	// Now make H_X0_GIVEN_X1... and H[X0|X1=x1]
    this->H_X0_GIVEN_X1 = 0.0;
	this->H_X0_given_X1_equals_x1.clear();	//pre-computed H[X0|X1=x1]


	// Foreach output state x1...
	unsigned int x1_index=0;
    for( map<unsigned int,unsigned int>::iterator it=this->existing_states.begin(); it != this->existing_states.end(); it++ ) {
		const unsigned int x1=(*it).first;
		const unsigned int num_x0s=(*it).second;


		///////////////////////////////////////////////////////////////
		// 1. add this x1 to our list of x1 states
		// 2. add this num_x0s to our array of num_x0s states		
		this->x1_states[x1_index] = x1;
//		this->num_x0s_to_this_x1[x1_index] = num_x0s;
		x1_index += 1;
		///////////////////////////////////////////////////////////////
		
		// sanity check the possible_x0s data-structure.  The number of x0s should be equal to the num_x0s as reported by existing_states.
//		assert( this->possible_x0s[x1].size() == num_x0s );

		///////////////////////////////////////////////////////////////		
		// 3. make H[X0|X1=x1] and sanity check
		
		this->H_X0_given_X1_equals_x1[ x1 ] = log2(num_x0s);
		assert( 0 <= this->H_X0_given_X1_equals_x1[x1] );
		assert( this->H_X0_given_X1_equals_x1[x1] <= this->numunits );
		///////////////////////////////////////////////////////////////
		
		///////////////////////////////////////////////////////////////
		// 3. now do H[X0|X1]
		this->H_X0_GIVEN_X1 += num_x0s * log2( num_x0s );
		///////////////////////////////////////////////////////////////
	}
	
	// sort the x1_states in ASCENDING ORDER
	sort( x1_states, x1_states + NUM_X1_STATES );

	// put the num_x0s_to_this_x1 in sync with x1_states
	for( int x1_index=0; x1_index<NUM_X1_STATES; x1_index++ ) {
		unsigned int x1 = x1_states[x1_index];
		num_x0s_to_this_x1[x1_index] = existing_states[ x1 ];
	}

	
    this->H_X0_GIVEN_X1 /= double(this->numstates);	
//	this->H_X0_GIVEN_X1 = this->H_X0_GIVEN_X1;
	
	// Calculate H(X1) using...
	// I(X0:X1) = H(X1) - H(X1|X0) = H(X1) = N - H(X0|X1)
	this->H_X1 = this->numunits - this->H_X0_GIVEN_X1;

	assert( 0.0 <= this->H_X1 );
	assert( this->H_X1 <= this->numunits );

//	cout << "H(X1)=" << this->H_X1 << endl;
//	cout << "N-H(X0|X1)=" << this->numunits - this->H_X0_GIVEN_X1 << endl;	
		
	// make FULL_MASK
	this->FULL_MASK = (this->numstates-1);
	
	// clear some memory.  We don't need these anymore.
	this->node_thresholds.clear();		// for neural networks

	// assert that one of these is nonempty, but only 1.
	assert( (((int) this->weightmatrix.empty() + (int) this->transformations.empty()) == 1) || this->networktype == "transition" );

	// assert that the number of units isn't bigger than the number of bits in an unsigned int
	assert( this->numunits <= (sizeof(unsigned int)*8) );
	
	// initialize the cache(s)
	assert( H_M1_cache == NULL );
	assert( H_M0_GIVEN_M1_cache == NULL );
	assert( prob_s1_given_mu0__Vnodes == NULL );

	// this memory will be there for the entire time the program runs
	prob_s1_given_mu0__Vnodes = new unsigned int[numunits - 1];	
	
	
	if( FLAG__CACHE_M1 ) {
		H_M1_cache = new double[numstates];
		for( int i=0; i<numstates; i++ )
			H_M1_cache[i] = DBL_NONE;		
	}
	

	if( FLAG__CACHE_M0_GIVEN_M1 ) {
		H_M0_GIVEN_M1_cache = new double[numstates];
		for( int i=0; i<numstates; i++ )
			H_M0_GIVEN_M1_cache[i] = DBL_NONE;		
	}
	
}


void t_consciousness::show_rules(void)
//shows the rules loaded	
{
	
	assert( networktype == "circuit" );
	const char output_rule[] = {'A','O','X','E','N'};

	for(int i=0;i<this->numunits;i++)
		printf("rule %i: %i = %i %c %i\n",i,i,transformations[i].A,output_rule[transformations[i].operation],transformations[i].B);
}


inline unsigned int t_consciousness::Xstate2MUstate_exchg( const unsigned int* restrict XOUTnodes, const unsigned int* restrict MUDESTnodes, const unsigned int partsize, const unsigned int MUINmask, const unsigned int Xstate )
// converts an Xstate to a MUstate
// Has the requirement that Xstate0 and Xstate1 be PARTIAL BIT-MASKS ONLY CONTAINING 1-bits that are part of part.
// ex: Xstate0 &= partmask;
{
	//idea, get rid of the &1 by zeroing all of the bits not in the mask and then going in the mask from greatest -> least order?
	
	// For the Xvals with indices INSIDE of MU, just pass them along.
	unsigned int MUstate = Xstate & MUINmask;
	
	// For Xvals outside MU, get their value and shift by a MUDEST val.
	for( int i=0; i<partsize; i++ )
		MUstate |= ((Xstate >> XOUTnodes[i])&1) << MUDESTnodes[i];
	
	
	return MUstate;
}



inline void t_consciousness::Xstate2MUstate( const unsigned int* restrict part, const unsigned int partsize, const unsigned int Xstate0, unsigned int& MUstate0 )
// converts an Xstate to a MUstate
// Has the requirement that Xstate0 and Xstate1 be PARTIAL BIT-MASKS ONLY CONTAINING 1-bits that are part of part.
// ex: Xstate0 &= partmask;
{
	//idea, get rid of the &1 by zeroing all of the bits not in the mask and then going in the mask from greatest -> least order?
	
	
	// set the first bit
	MUstate0 =(( Xstate0 >>part[partsize-1])&1);
//	MUstate0 = Xstate0 >>part[partsize-1];
	

// for all future bits, first shift to the RIGHT, then add the bit on the end.
	for(int i=partsize-2; i>=0; i--)
		MUstate0 = (MUstate0<<1) | ((Xstate0 >>part[i])&1);
	
	
	// for all future bits, first shift to the RIGHT, then add the bit on the end.
//	for(int i=partsize-2; i>=0; i--) {
//		MUstate0 = (MUstate0<<1) | ((Xstate0 >>part[i])&1);
//		assert( part[i] < part[i+1] );
//	}
	
}



inline void t_consciousness::Xstate2MUstate( const unsigned int*  restrict part, const unsigned int partsize, const unsigned int Xstate0, unsigned int& MUstate0, const unsigned int Xstate1, unsigned int& MUstate1 )
// converts an Xstate to a MUstate
// converts two of them
// Has the requirement that Xstate0 and Xstate1 be PARTIAL BIT-MASKS ONLY CONTAINING 1-bits that are part of part.
// ex: Xstate0 &= partmask;
{
	// set the first bit
//	MUstate0=(( Xstate0 >> part[partsize-1])&1);
//	MUstate1=(( Xstate1 >> part[partsize-1])&1);
	
	MUstate0=Xstate0 >> (part[partsize-1])&1;
	MUstate1=Xstate1 >> (part[partsize-1])&1;
	
	
	
	// for all future bits, first shift to the RIGHT, then add the bit on the end.
	for(int i=partsize-2; i>=0; i--) {
		MUstate0 = (MUstate0<<1) | ((Xstate0>>part[i]) &1 );
		MUstate1 = (MUstate1<<1) | ((Xstate1>>part[i]) &1 );

//		MUstate0 = (MUstate0<<1) | ((Xstate0>>part[i]) );
//		MUstate1 = (MUstate1<<1) | ((Xstate1>>part[i]) );
		assert( part[i] < part[i+1] );		
	}
	
}



double t_consciousness::H_M0_GIVEN_M1__WIRES( unsigned int part_mask )
// Computes H[M_0|M_1] for a given part (M) of the network by PERTURBING THE WIRES
{
    // this function should work even if not FULL_MASK, but we haven't tested it in a while.
	assert( part_mask == FULL_MASK );

	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Define the stuff we need
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// create the nodes[] array
	unsigned int numnodes = numunits_in_subset( part_mask );
	unsigned int* restrict nodes = new unsigned int[numnodes];
	subset2array(part_mask, nodes);
	
	const unsigned int d = 1 << numnodes;		//d=2**part.size()	
	const double Dnumstates = (double) numstates;

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// now define a 2-dimensional array as a single contiguous block of memory
	// access the i,jth element via: myarray[i * ncolumns + j]
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	double* prob_mu1_given_mu0 = (double*) calloc(d * d, sizeof(double));
	assert( prob_mu1_given_mu0 != NULL );

	
	// temp matrix for calculating the intersection
	double* prob_mu0_mu1 = (double*) calloc(d * d, sizeof(double));
	assert( prob_mu0_mu1 != NULL );
	
	
	// define a block of memory to hold the number of times a mu1 occurs
	unsigned int* mu1_instances = (unsigned int*) calloc(d, sizeof(unsigned int));	
	assert( mu1_instances != NULL );
	
	ASSERT( sum(mu1_instances,d) == 0 );
	ASSERT( fequals(sum(prob_mu1_given_mu0,d*d), 0.0) );
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	

	//foreach x0...
	for(unsigned int temp_x0=0;temp_x0<=FULL_MASK;temp_x0++)
	{
		unsigned int temp_mu0, temp_mu1;
		this->Xstate2MUstate( nodes, numnodes, temp_x0, temp_mu0, this->states[temp_x0], temp_mu1 );

		////////////////////////////////////////////////////
		// Calculate p(mu1).  We will use this to calculate p(mu0 AND mu1) = p(mu0|mu1) * p(mu1)
		// this will score the #times each mu1 came up
		////////////////////////////////////////////////////
		mu1_instances[temp_mu1] += 1;
		////////////////////////////////////////////////////
		
		
		/////////////////////////////////////////////////////////
		// Now calculating p(mu1|mu0)
		/////////////////////////////////////////////////////////
		double this_prob_mu1_given_mu0 = 1.0;
		
		//foreach node_num.... (node_index)
//		cout << "Mask: " << binary_repr(part_mask) << "\t x0=" << binary_repr(temp_x0) << " \t Calculating p(mu1|mu0) over " << numnodes << " nodes." << endl;
		for( int node_index=0; node_index<numnodes; node_index++ )
		{
			this_prob_mu1_given_mu0 *= prob_node_matches_mu1( part_mask, node_index, temp_mu0, temp_mu1, nodes, numnodes );
			ASSERT( 0.0 < this_prob_mu1_given_mu0 );
			ASSERT( this_prob_mu1_given_mu0 <= 1.0 );
		}
		
		ASSERT( this_prob_mu1_given_mu0 == 1.0 || part_mask != FULL_MASK );

		prob_mu1_given_mu0[ (temp_mu1*d)+temp_mu0 ] += this_prob_mu1_given_mu0;    //*(1.0/d)*(1.0/prob_mu1);		
		prob_mu0_mu1[ (temp_mu1*d)+temp_mu0 ] += (this_prob_mu1_given_mu0 * (1.0/d));
	}

	// we're done with the nodes[] now.
	delete [] nodes;

/*	
	////////////////////////////////////////////////////////////////////
	// Creating the intersection matrix by normalizing by its sum
	////////////////////////////////////////////////////////////////////
	double bigsum = sum(prob_mu0_mu1, d*d);
	for( int i=0; i<d*d; i++ )
		prob_mu0_mu1[ i ] /= bigsum;
	
	////////////////////////////////////////////////////////////////////
	// Go from p(mu0,mu1) -> p(mu0|mu1) by 
	// dividing by p(mu1)
	////////////////////////////////////////////////////////////////////
	// make a pointer to the same location in memory.
	double* prob_mu0_given_mu1 = prob_mu1_given_mu0;
	prob_mu1_given_mu0 = NULL;		// null out the pointer

	for( int mu1=0; mu1<d; mu1++ ) {

//		if( mu1_instances[mu1] == 0 )
//			continue;
//		const double prob_mu1 = mu1_instances[mu1] / Dnumstates;
		
		const double prob_mu1 = sum( prob_mu0_mu1+(mu1*d), d );
		if( fequals( prob_mu1, 0.0 ) )
			continue;
		
		
		for( int mu0=0; mu0<d; mu0++ ) {
			prob_mu0_given_mu1[ (mu1*d)+mu0 ] = prob_mu0_mu1[(mu1*d)+mu0] / prob_mu1;
		}

	}

	cout << "sum=" << sum(prob_mu0_given_mu1,d*d) << "\t\t #s1states=" << num_s1_states(part_mask) << endl;
	ASSERT( fequals(sum(prob_mu0_given_mu1,d*d), num_s1_states(part_mask)) );
*/
	
	free( prob_mu0_mu1 );
	
	////////////////////////////////////////////////////////////////////
	// Go from p(mu1|mu0) -> p(mu0|mu1) by 
	// dividing by the sum down the columns of p(mu1|mu0)
	////////////////////////////////////////////////////////////////////
	// make a pointer to the same location in memory.
	double* prob_mu0_given_mu1 = prob_mu1_given_mu0;
	prob_mu1_given_mu0 = NULL;		// null out the pointer
	////////////////////////////////////////////////////////////////////
	
	for( int mu1=0; mu1<d; mu1++ ) {

		double sum1 = 0.0;
		for( int mu0=0; mu0<d; mu0++ )
			sum1 += prob_mu0_given_mu1[ (mu1*d)+mu0 ];
		
		// sum down the column
		// OPTIMIZATION:  THIS IS FASTER THAN A LOOP.
		const double col_sum = sum( prob_mu0_given_mu1+(mu1*d), d );
		assert( fequals( sum1, col_sum ) );

		
		if( fequals( col_sum, 0.0 ) ) {
			assert( col_sum == 0.0 );
			continue;
		}
		
		for( int mu0=0; mu0<d; mu0++ )
			prob_mu0_given_mu1[ (mu1*d)+mu0 ] /= col_sum;
	}

	cout << "sum=" << sum(prob_mu0_given_mu1,d*d) << "\t\t #s1states=" << num_s1_states(part_mask) << endl;
	ASSERT( fequals(sum(prob_mu0_given_mu1,d*d),num_s1_states(part_mask)) );

	
	////////////////////////////////////////////////////////////////////
	// Now compute the entropies
	////////////////////////////////////////////////////////////////////
	double z=0.0;
	double sum_prob1=0.0, sum_prob2=0.0;
	for(int mu1=0; mu1<d; mu1++)
	{
		if( mu1_instances[mu1] == 0 )
			continue;
		
		const double prob_mu1 = mu1_instances[mu1] / Dnumstates;
		sum_prob1 += prob_mu1;
		
		for( int mu0=0; mu0<d; mu0++ )
		{			
			if( fequals(prob_mu0_given_mu1[ (mu1*d)+mu0 ], 0.0) )
				continue;
						
			const double pr_mu0_given_mu1 = prob_mu0_given_mu1[ (mu1*d)+mu0 ];
			const double pr_mu0_mu1 = pr_mu0_given_mu1 * prob_mu1;
			
			// sanity check
			sum_prob2 += pr_mu0_mu1;
			
//			cout << "p(mu0,mu1)=" << pr_mu0_mu1 << "\t\t"  << "p(m0|mu1)=" << pr_mu0_given_mu1 << endl;
			z -= ( pr_mu0_mu1 ) * log2( pr_mu0_given_mu1 );
		}
		
	}
	
	// sanity check that the sum_prob over the joint distribution is 1.0
//	assert( fequals(sum(prob_mu0_mu1,d*d), 1.0) );
//	assert( fequals(sum(prob_mu1,d), 1.0) );	
	assert( fequals( sum_prob1, 1.0 ) );
	assert( fequals( sum_prob2, 1.0 ) );
	assert( prob_mu1_given_mu0 == NULL );	

	free(mu1_instances);
	free(prob_mu0_given_mu1);

	if( z < 0.0 ) {
		if( myceil(z) < 0.0 ) {
			cerr << "Warning: Had negative entropy in entropy_of_part() that didn't round to 0.0 (" << myceil(z) << ").  Maybe decimal precision is too high?" << endl;
			cerr.flush();
			assert( 0.0 <= myceil(z) );
		}
		z = 0.0;
	}

	// DO NOT DO ANY ROUNDING HERE!
	
	ASSERT( z == z );
	ASSERT( z >= 0.0 );	
	
	return z;
}


double t_consciousness::H_M0_GIVEN_M1__ELEMENTS( unsigned int part_mask )
// Computes H[M_0|M_1] for a given part (M) of the network PERTURBING THE STATES
{
	// if we've done this calculation before, use the cache (if we use that)
//	assert( FLAG__CALC_EI_PERTURBING == "ELEMENTS" || part_mask == FULL_MASK );	
//	assert( part_mask != FULL_MASK );

	unsigned int numnodes = numunits_in_subset( part_mask );
	unsigned int* restrict nodes = new unsigned int[numnodes];
	subset2array(part_mask, nodes);

	unsigned int d = 1 << numnodes;		//d=2**part.size()

	
	// define a block of memory to hold the number of times a mu1 occurs
	unsigned int* restrict mu1_instances = (unsigned int*) calloc(d, sizeof(unsigned int));	
	assert( mu1_instances != NULL );

	// now define a 2-dimensional array as a single contiguous block of memory
	// this requires us to create a big-ass 1-d array but we PRETEND it is a 2-d array
	// However, we must now perform subscript calculations manually, accessing the i,jth element via:
	// myarray[i * ncolumns + j]	
	unsigned int* restrict mu0mu1_instances = (unsigned int*) calloc(d * d, sizeof(unsigned int));
	assert( mu0mu1_instances != NULL );
	

	//////////////////////////////////////////////////////////////////////////////////////
	// ATTEMPTS HERE TO USE SPARSE MATRICES usually resulted in slower code             //
	// i.e: map<int,double> prob_mu1s; vector< map<int,double> > prob_mu0_given_mu1(d); //
	//////////////////////////////////////////////////////////////////////////////////////

	
	//foreach x0...
	unsigned int cur_M0_mask=0, cur_M1_mask=0;
	unsigned int old_M0_mask=UINT_MAX, old_M1_mask=UINT_MAX;
	unsigned int temp_mu0=0, temp_mu1=0;
		
	for(unsigned int temp_x0=0; temp_x0<this->numstates; temp_x0++ )
	{
		cur_M0_mask = temp_x0 & part_mask;
		cur_M1_mask = states[temp_x0] & part_mask;
		
//		int temp_x1=this->states[temp_x0];
		if( cur_M0_mask != old_M0_mask || cur_M1_mask != old_M1_mask ) {
//			Xstate2MUstate( nodes, numnodes, temp_x0, temp_mu0, this->states[temp_x0], temp_mu1 );
			Xstate2MUstate( nodes, numnodes, cur_M0_mask, temp_mu0, cur_M1_mask, temp_mu1 );			
			old_M0_mask = cur_M0_mask;
			old_M1_mask = cur_M1_mask;
		}

		//printf("%i	%i\n",temp_mu0,temp_mu1);
		mu1_instances[temp_mu1] += 1;

//		mu0mu1_instances[temp_mu1][temp_mu0] += 1;
		mu0mu1_instances[ (temp_mu1*d) + temp_mu0] += 1;
	}

	// don't need the nodes anymore.
	delete [] nodes;
	
	const double Dnumstates = (double) this->numstates;
	double z=0.0, sum_prob=0.0;
	
	for(int mu1=0; mu1<d; mu1++)
	{
		if( mu1_instances[mu1] == 0 )
			continue;
		
		const double Dmu1_instances = (double) mu1_instances[mu1];
		
		for( int mu0=0; mu0<d; mu0++ )
		{
			if( mu0mu1_instances[(mu1*d)+mu0] == 0 )
				continue;

			//Make each row a conditional probability by dividing by #mu1's/
//			//Now run over both to return the information using the probabilities of p(mu0,mu1) of p(mu0|mu1)
//			z -= ( mu0mu1_instances[mu1][mu0] / Dnumstates ) * log2( (mu0mu1_instances[mu1][mu0] / Dmu1_instances) );
			z -= ( mu0mu1_instances[(mu1*d)+mu0] / Dnumstates ) * log2( (mu0mu1_instances[(mu1*d)+mu0] / Dmu1_instances) );
			
			sum_prob += mu0mu1_instances[(mu1*d)+mu0] / Dnumstates;
		}
		
	}

	assert( fequals(sum_prob, 1.0) );
	
	// free the memory in both.
	free( mu0mu1_instances );
	free( mu1_instances );
	
	if( z < 0.0 ) {
		if( myceil(z) < 0.0 ) {
			cerr << "function: entropy_of_part__ELEMENTS" << endl;
			cerr << "Warning: Had negative entropy in entropy_of_part() that didn't round to 0.0 (" << myceil(z) << ").  Maybe decimal precision is too high?" << endl;
			assert( 0.0 <= myceil(z) );
		}
		z = 0.0;
	}


	// DO NOT DO ANY ROUNDING HERE!
	
	ASSERT( z >= 0.0 );
	
	return z;
}


double t_consciousness::prob_node_matches_mu1( unsigned int part_mask, unsigned int node_index, unsigned int mu0, unsigned int mu1, const unsigned int* restrict part_nodes, unsigned int numnodes )
// returns the probability of an arbitrary node in the network matching a passed state mu1
// part_mask = bitmask of the entire mask
// node_index = the index within the part for the node we care about
{
	assert( part_mask > 0 );
	assert( part_mask <= FULL_MASK );

	// ensure that mu0 and mu1 have no bits ON outside the partmask.
	assert( numunits_in_mask(part_mask) == numnodes );	
	assert( numunits_in_mask(mu0,true) <= numnodes );
	assert( numunits_in_mask(mu1,true) <= numnodes );
	assert( node_index <= numnodes );

	// the value of the node at timestep 1
	uint n = part_nodes[node_index];
	bool n1 = (mu1 >> node_index) & 1;
	
//	cout << "\t mu1=" << binary_repr(mu1) << "\t n1_" << node_index << " = " << n1 << endl;

	vector<int> x0states = this->states_to_attempt__WIRES( n, mu0, part_nodes, numnodes );
	
	// if there are no external wires to this part, then the probability of matching mu1 is 1.0
	if( x0states.empty() ) {
//		cout << "\t\t No states to try." << endl;
		return 1.0;
	}

//	cout << "# x0states:" << x0states.size() << endl;

	uint n1_instances = 0;
	for( int i=0; i<x0states.size(); i++ )
	{
		uint temp_x0 = x0states[i];
		uint temp_x1 = this->states[temp_x0];
		
//		cout << "\t\t" << "x0:" << binary_repr( temp_x0 ) << " -> " << binary_repr( temp_x1 ) << endl; 

//		int temp_mu1 = temp_x1 & part_mask;		
		bool temp_n1 = (temp_x1 >> n) & 1;

		if( temp_n1 == n1 )
			n1_instances += 1;		
	}
	
	ASSERT( 0 < n1_instances );
	
	double prob_node_matches_mu1 = n1_instances / (double) x0states.size();
	
	assert( 0.0 < prob_node_matches_mu1 );
	assert( prob_node_matches_mu1 <= 1.0 );

//	cout << "\t n1_instances / x0states.size() = " << n1_instances << "/" << x0states.size() << "=" << prob_node_matches_mu1 << endl;	
	
	return prob_node_matches_mu1;
}


double t_consciousness::entropy_of_part_given_x1__WIRES( const unsigned int s1, const bitmask partmask, const bitmask subset )
// Get the entropy of a part by perturbing the WIRES
{
	assertmsg( subset == FULL_MASK, "Don't support calculating p(s0|s1) for non-fullmask subsets yet" );
	assert( subset <= FULL_MASK );
	assert( (s1|subset) == subset);
	assert( (partmask|subset) == subset );
	
	const unsigned int partsize = numunits_in_mask(partmask);
	const unsigned int Ssize = numunits_in_mask(subset);
	
	unsigned int* restrict nodes = new unsigned int[partsize];
	subset2array(partmask, nodes);
	
	const unsigned int mu1 = s1 & partmask;
	const unsigned int d = 1 << partsize; // d = 2**part.size()
	
	// these could be converted to arrays of doubles
//	vector<double> prob_mu1_given_M0( d, 1.0 );			// d vector of 0.0
//	prob_mu1_given_M0.clear();
	
//	vector<double> prob_M0_given_mu1( d, 0.0 );			// d vector of 0.0
//	prob_M0_given_mu1.clear();

//	// how many x1 matches for each mu0 ?
//	vector<unsigned int> num_x1_matches(d, 0);			// d elements of 0

	// if there is no subset, set it to the complete mask.	
	unsigned int s1_instances=0;
//	unsigned int mu1_instances=0;
	
	
	/////////////////////////////////////////////////////////
	// Now calculate p(s1)
	/////////////////////////////////////////////////////////
	//foreach x0...
	for(unsigned int x0=0;x0<=FULL_MASK;x0++)
	{
		if( (states[x0]&subset) == s1 )
			s1_instances += 1;
	}
	
	const double prob_s0 = 1.0 / (double) (1 << Ssize);
	const double prob_mu0 = 1.0 / (double) (1 << partsize);		
	const double prob_s1 = (double) s1_instances / (double) numstates;	
	
	/////////////////////////////////////////////////////////
	// Now calculating p(mu1|M0)
	/////////////////////////////////////////////////////////

	//foreach mu0...
	double z=0.0;
	for( unsigned int x0=0; x0<=FULL_MASK; x0++ )
	{
		
		// If this x0 doesn't goto this s1, skip it.
		if( (states[x0]&subset) != s1 )
			continue;
		
		const double prob_s0_given_s1 = prob_s0 / prob_s1;
		
		assert( prob_s0_given_s1 > 0.0 );
		
		const unsigned int mu0 = x0 & partmask;
		
		const double this_prob_mu0_given_mu1 = prob_mu0_given_mu1__anticond( mu0, mu1, partmask, partsize );
//		const double this_prob_mu0_given_mu1 = prob_mu0_given_s1( mu0, partmask, mu1, partmask );
		
		// if it's 0.0, then we don't need to do this one.
//		if( this_prob_mu0_given_mu1 == 0.0 )
//			continue;

		const double logterm = -1.0 * log2( this_prob_mu0_given_mu1 );
//		cout << "\t p(s0|s1)=" << prob_s0_given_s1;
//		cout << "\t p(mu0|mu1)=" << this_prob_mu0_given_mu1 << "\t -log2( p(mu0|mu1) )=" << logterm << endl;
//		cout << "z=" << z << endl;
		
		z += (prob_s0_given_s1 * logterm);
	}
	
	// remove -0.0's
	if( fequals(z,0.0) )
		z = 0.0;
	
	assert( 0.0 <= z );
	
    // don't need this anymore!
    delete [] nodes;
    
//	cout << "s1=" << s1 << "\tz=" << z << endl;
	return z;	
}



double t_consciousness::entropy_of_part_given_x1( const unsigned int s1, const bitmask partmask, const bitmask subset )
// get the entropy of a part by perturbin the ELEMENTS
{
	// assert that both s1 and partmask are a subset of subset
	assert( (s1|subset) == subset);
	assert( (partmask|subset) == subset);	
	ASSERT( subset <= this->FULL_MASK );
	
	const unsigned int partsize = numunits_in_subset(partmask);
	unsigned int* restrict nodes = new unsigned int[partsize];	
//	unsigned int* restrict Mnodes = new unsigned int[M1size];

	mask2array( partmask, nodes );

	
//	cout << "\t\tentropy_of_part_given_x1(x1=" << x1 << "\t part=" << this->part2uint(part) << ")" << endl;
	
	const unsigned int d = 1 << partsize; // d = 2**part.size()
	vector<unsigned int> partition_vector_S1(d, 0);			// d elements of 0
	vector<unsigned int> partition_vector_M1(d, 0);			// d elements of 0
	
	// make the target mu1
	const unsigned int mu1 = s1 & partmask;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	unsigned int mu1_instances=0, s1_instances=0;
	
	//foreach x0...
	for(unsigned int x0=0;x0<=FULL_MASK;x0++)
	{
		if( (states[x0]&partmask) == mu1 )
		{
			// figure out what the mu0 is
			unsigned int mu0_index=0;
			Xstate2MUstate( nodes, partsize, x0, mu0_index );

			mu1_instances += 1;									// #times mu1 has appeared		
			partition_vector_M1[mu0_index] += 1;						// this mu0 leads to mu1
			
			// temp_x1 can only match x1 if temp_mu1 matches mu1
			if( (states[x0]&subset) == s1 ) {
				s1_instances += 1;
				partition_vector_S1[mu0_index] += 1;
			}
		}
	}

	delete [] nodes;
	
	assert( s1_instances <= mu1_instances );

	//convert to doubles so we don't have to typecast all the time.
	const double Dmu1_instances=(double) mu1_instances, Ds1_instances=(double) s1_instances;
	double z=0.0;
	
	for(unsigned int mu0=0;mu0<d;mu0++)
	{
		
		// if p(X1=x1) > 0.0, then p(M1=mu1) is guaranteed to be >0.0 so we don't have to check that: partition_vector_M1[mu0] != 0
		if( partition_vector_S1[mu0] == 0 )
			continue;
		
		// divide by the number of  s1_instances to get p(S0=i|S1=s1)
		// divide by the number of mu1_instances to get p(M0=i|M1=mu1)
		z -= (partition_vector_S1[mu0] / Ds1_instances) * log2( (partition_vector_M1[mu0] / Dmu1_instances) );
	}
	
	
	if( fequals(z,0.0) )
		z = 0.0;
	
	assert( 0.0 <= z );
	
	return z;
}

double t_consciousness::H_S0_given_S1( unsigned int subset )
// Returns H[S0|S1] by calculating H[M_0|M_1]
{
	double z = H_M0_GIVEN_M1( subset );
	
	if( z < 0.0 ) {
		assert( 0.0 <= myceil(z) );
		z = 0.0;
	}
	
	assert( z >= 0.0 );
	
	return z;	
}

// alias to the unsigned int version
double t_consciousness::H_S0_given_S1( const t_subset& restrict S ) { return H_S0_given_S1( S.S ); }


double t_consciousness::H_S0_given_S1_equals_s1( const unsigned int x1, const bitmask subset )
// Returns H[S0|S1=s1]
{

	// 0 <= s1 <= mask <= numstates
	ASSERT( 1 <= subset );
	ASSERT( subset <= this->FULL_MASK );

	assert( (x1|subset) == subset);
	
	const unsigned int mu1 = x1 & subset;
	unsigned int mu1_instances=0;
	
	const unsigned int d = 1 << this->numunits_in_subset(subset);		// d = 2**|S|
//	cerr << "doing d=2**N = 2**" << this->numunits_in_subset(subset) << " for subset=" << subset << "\t repr()=" << binary_repr(subset,this->numunits) << " and x1=" << x1 << " now." << endl;
//	cerr << "------------------------------------" << endl;

//	cerr << "\td=" << d << endl;	
//	cerr.flush();

	map<int, double> prob_mu0_given_mu1;
	prob_mu0_given_mu1.clear();
	ASSERT( prob_mu0_given_mu1.max_size() >= d );
	
	//foreach x0...
	for(unsigned int temp_x0=0;temp_x0<this->numstates;temp_x0++)
	{
		// the x1 of this x0
		unsigned int temp_x1=this->states[temp_x0];

		// make the temp_mu0 and temp_mu1
		int temp_mu0=temp_x0&subset, temp_mu1=temp_x1&subset;
	
		ASSERT( 0 <= temp_mu0 && temp_mu0 <= subset );
		ASSERT( 0 <= temp_mu1 && temp_mu1 <= subset );

		
//		cerr << "temp: t_x0=" << temp_x0 << "\t t_x1=" << temp_x1 << "\t t_mu0=" << temp_mu0 << "\t t_mu1=" << temp_mu1 << endl;
		
		if( temp_mu1 == mu1 )
		{
			prob_mu0_given_mu1[temp_mu0] += 1.0;					// this mu0 leads to mu1
			mu1_instances += 1;										// #times mu1 has appeared
		}
	}
	

	// this size should never be greater than #mu0 states
	ASSERT( prob_mu0_given_mu1.size() <= d );
	
//	cerr << "mu1_instances=" << mu1_instances << " \t #x0s->x1: " << this->existing_states[x1] << " \t numstates=" << this->numstates << endl;
	
	double z=0.0;
	for( map<int,double>::iterator it=prob_mu0_given_mu1.begin(); it != prob_mu0_given_mu1.end(); it++ )
	{
		int mu0=(*it).first;
		// we must divide by prob(mu1) = #mu1_instances to get a conditional probability
		prob_mu0_given_mu1[mu0] /= mu1_instances;
		z += prob_mu0_given_mu1[mu0] * log2(prob_mu0_given_mu1[mu0]);
	}

	
	z *= -1.0;

	// 0.0 <= z <= |S|
	ASSERT( 0.0 <= z );
	ASSERT( z <= this->numunits_in_subset(subset) );
	
	return z;
}


vector<t_phi_result> t_consciousness::highest_bracketphis()
// This function iterates over all possible subsets and returns the bracketphi, and subsets -> average_MIPs for subsets with bracketphi <= highest_found_bracketphi bracket
{
	
	// THIS INITIAL VALUE IS IMPORTANT.  DO NOT SET IT TO DBL_MIN.
	double highestphi=0.0;
	vector<t_phi_result> results;
	
	//foreach subset...calculate the MIP and the PHI
	//we start from the FULL_MASK because we can often skip steps that way
	for( unsigned int subset=FULL_MASK; subset>=1; subset-- )
	{
		//		cout << "Starting subset=" << subset << endl;
		//OPTIMIZATION: only attempt all partitions to find the MIPs if it's possible
		//for this subset to have higher or equal PHI to the current highestphi
		if( this->bracket_ei(subset) < highestphi )
			continue;
		
		t_phi_result result = this->bracketphi( subset );

		if( result.MIPs.empty() ) {
			assert( use_total_partition == false );
			assert( result.subset.numnodes() == 1 );
			continue;
		}
		
		// if this subset has the highest phi seen thus far, set it as the highest phi, MIPs, and subset
		if( result.min_phi > highestphi )	{
			results.clear();
			results.push_back( result );
			highestphi = result.min_phi;
			
			//			cout << "\tsubset=" << result.subset << " had a new highestphi!  phi=" << highestphi << "\tphi_norm=" << result.phi_norm << "\t#MIPs=" << result.MIPs.size() << endl;			
		}
		
		// if this PHI is EQUAL to the current highestphi, add this subset and its MIPs to the highestphis
		// only push back if we're BELOW the maximum number of complexes
		else if( fequals(result.min_phi, highestphi) && results.size() < MAXIMUM_NUMBER_OF_MAIN_COMPLEXES_TO_STORE ) {
			highestphi = result.min_phi = max(highestphi, result.min_phi);
			results.push_back( result );
		}
	}

	return results;
}



vector<t_phi_result> t_consciousness::highest_phis( unsigned int x1 )
/* This function iterates over all possible subsets and returns the phi, MIPs, and subsets that had the highest found phi for a given x1 state */
{
	double highestphi=-1.0;
	vector<t_phi_result> results;
	t_phi_result result;
	
	//foreach subset...calculate the MIP and the PHI
	// 0 <= subset <= (2**numstates)-2  -- this is because ei(subset=0) == ei(subset=numstates-1)
	// So we don't have to redo the same work, make it from [0...numstates-2]
	for( unsigned int subset=1; subset<=this->FULL_MASK; subset++ )
	{
		
		//OPTIMIZATION: only attempt all partitions to find the MIPs if it's possible
		//for this subset to have higher or equal PHI to the current highestphi
		// ei() requires that s1 be passed, not x1.
		if( this->ei(x1, subset) < highestphi )
			continue;
		
		result = this->phi( x1, subset );
		
		// if this subset has the highest phi seen thus far, set it as the highest phi, MIPs, and subset
		if( result.min_phi > highestphi )	{
			results.clear();
			results.push_back( result );
			highestphi = result.min_phi;
		}
		
		// if this PHI is EQUAL to the current highestphi, append the MIPs and subset
		else if( result.min_phi == highestphi ) {
			results.push_back( result );

			// stop if we have the maximum number of main complexes
			if( results.size() >= MAXIMUM_NUMBER_OF_MAIN_COMPLEXES_TO_STORE )
				break;			
		}
		
	}
	
	return results;
}


// just an alias to the bitmask version
double t_consciousness::prob_s1( const t_state& s1 ) { return prob_s1( s1.value, s1.mask ); }

double t_consciousness::prob_mu0( const t_state& mu0 )
// returns 1.0 / 2^(|M|)
{
    assert( is_valid_mu0(mu0) );
    
    double z= 1.0 / (1 << mu0.size());
    
    assert( 0.0 < z );
    assert( z < 1.0 );
    return z;
}

double t_consciousness::prob_mu0_given_s1( const t_state& mu0, const t_state& s1, const unsigned int MUsize )
// simply an alias to the bitmask version.
{
    return prob_mu0_given_s1( mu0.value, mu0.mask, s1.value, s1.mask, MUsize );
}


double t_consciousness::prob_a0_given_b0_c1( const unsigned int a0, const bitmask A0mask, const unsigned int b0, const bitmask B0mask, const unsigned int c1, const bitmask C1mask )
// Returns p(a0|b0,c1) = p(a0,b0,c1) / p(b0,c1)
{
    // assert our c1 is valid
    assert( is_valid_mu1( c1, C1mask ) );
    
    // assert A0mask and B0mask don't overlap.  They technically could overlap, but so far we've never wanted to do that.
    assert( (A0mask & B0mask) == 0 );

    
    double z = prob_a0_b1( a0|b0, A0mask|B0mask, c1, C1mask );
    
    if( z == 0.0 )
        return 0.0;

    z /= prob_a0_b1( b0, B0mask, c1, C1mask );
    
    
    assert( 0 <= z );
    assert( z <= 1.0 );
    
    return z;    
}

double t_consciousness::prob_a0_b1( const unsigned int a0, const bitmask A0mask, const unsigned int b1, const bitmask B1mask )
// Returns the probability of state a0 and b1 co-occurring.  ex: p(a0,b1) = p(a0) * p(b1|a0)
{
    assert( is_valid_mu1( b1, B1mask ) );
    
    const unsigned int Asize = numunits_in_mask(A0mask);
    // p(a0) = 1.0 / 2^{|A|}
    const double prob_a0 = 1.0 / (1 << Asize);
    
    double z = prob_a0;
    z *= prob_s1_given_mu0( b1, B1mask, a0, A0mask, Asize );

    
    assert( 0 <= z );
    assert( z <= 1.0 );
    
    return z;
}

int t_consciousness::save_to_file(string ofilename, bool perstate, bool bracketphi, bool highestbracketphi, bool bracketMCs, bool higheststatephi, bool highestbracketei, bool higheststateei )
// Dumps the states, averages, and whatever else to a file.
{
	cout.unsetf( ios::showpoint );
    cout.setf( ios::fixed );
    cout << setprecision(2);
    
//    cout.setpercision(4)
	
	// assert that at least one of the bools is true
	assert( perstate || bracketphi || highestbracketphi || higheststatephi || highestbracketei || higheststateei || bracketMCs );
	
	time_t t1 = time(NULL);

	// make the tempfilename...
	string t_filename = ofilename;
	
	this->output_filename = ofilename;

	// delete the character "/" from the string	
	while( t_filename.find("/") != string::npos )
		t_filename.erase( int(t_filename.find("/")), 1 );
	
	string tempfilename = DEFAULT_TEMP_DIRECTORY + string("/,") + string( t_filename );
	
	// output to a temporary file.  At the copy to ofilename
	ofstream tempfile;
	tempfile.open( tempfilename.c_str() );
	if( ! tempfile.is_open() ) {
		cerr << "* Error! -- Could not open filename " << tempfilename << endl;
		return -1;
	}
	
	// Always output some basic network statistics
	tempfile << "# statistics:" << "\tN=" << numunits << "\tmax_num_parts=" << max_num_parts << "\t#x1states=" << NUM_X1_STATES;
	tempfile << "\tH(X1)=" << H_X1 << "\tH(X0|X1)=" << H_X0_GIVEN_X1 << "\tnormalization=" << FLAG__MIP_METHOD;
	tempfile << "\tEI-MODIFIER=" << FLAG__STATE_EI << endl;
	
	// print comments, if we have any
	if( ! this->network_comments.empty() )
	{
		for( int i=0; i<this->network_comments.size(); i++ ) {
			tempfile << "comment:\t" << this->network_comments[i] << endl;
		}	
	}
	
	// Output the perstate and averaged phi, if we're doing that.
	if( perstate )
	{
		double average_phi=0.0, average_ei=0.0, average_IbS=0.0, average_eiweakestlink=0.0;
		
		
		for( map<unsigned int,unsigned int>::iterator it=this->existing_states.begin(); it != this->existing_states.end(); it++ )
        {
			unsigned int x1 = (*it).first;
			double prob_x1 = (*it).second / (double) this->numstates;
			t_phi_result res = this->phi( x1 );

			// for displaying the S.brep();
			t_subset S(numunits, x1);
			
			// this obviously doesn't have to be true, but we'd like to know anytime it isn't.
            if( res.MIPs.size() >= 2 )
                cerr << "Warning: Found more than one MIP for state x1=" << S.brep() << endl;
                
			average_phi += prob_x1 * res.min_phi;
			average_ei += prob_x1 * res.ei;

			
			tempfile << "state:" << "\tx1=" << x1;
			tempfile << "\tbstr=" << S.brep();
			tempfile << "\tp(x1)=" << prob_x1;
			
			tempfile << "\tei=" << res.ei;
			
            if( res.min_phi > res.ei )
                cerr << "** Warning: state x1=" << S.brep() << " has phi > ei.  This is a known problem in state-dependent ɸ.  ei=" << res.ei << " \t phi=" << res.min_phi << " **" << endl;
            
            tempfile << "\tmin_phi=" << res.min_phi;
            
			tempfile << "\t#MIPs=" << res.MIPs.size();
			tempfile << "\tPHIs=" << PHIs2str(res.PHIs);
            
            /////////////////////////////////////////////////////////////////////////////////////////////////////////
            // This is the code for handling the IbS
            /////////////////////////////////////////////////////////////////////////////////////////////////////////
            // ensure there is only one MIP.  If there's more than one MIP, then we have more code to write to handle printing them all out.
            if( res.MIPs.size() >= 2 ) {
                cerr << "Had more than one MIP!  Printing results only for the first MIP." << endl;
            }

            
            //create the t_state for this x1 state.
            t_state x1state( res.MIPs[0].subset, x1 );
            
            t_subset s_FULLMASK( numunits, FULL_MASK );
            /////////////////////////////////////////////////////////////////////////////////////////////////////////
            
			tempfile << "\tMIPs=" << partitions2str(res.MIPs);
			
//			tempfile << "\tPHI_wires=" << ei( x1, res.MIPs[0], "WIRES" );
//          tempfile << "\tPHI_units=" << ei( x1, res.MIPs[0], "UNITS" );

			tempfile << endl;
		}
		
		assertmsg( fequals(average_ei,H_X1), "state-averaged ei(x1) != H(X1)" );
		tempfile << "state-averaged:" << "\tei=" << average_ei << "\tphi=" << average_phi;
		tempfile << endl;
	}
	
    
	// Output the bracketphi, if we're doing that.
	if( bracketphi )
	{
		t_phi_result r = this->bracketphi( this->FULL_MASK );
		r.round();

		tempfile << "bracket:" << "\tei=" << r.ei << "\tMIPscore=" << r.mip_score << "\t#MIPs=" << r.MIPs.size() << "\tPHIs=" << PHIs2str(r.PHIs, r.ei) << "\tMIPs=" << partitions2str(r.MIPs);
        
        t_partition P = r.MIPs[0];
        

		tempfile << endl;
		
	}
	
	
	// if we're doing the atomic partition for bracketphi, do that now.
	if( SAVE_ATOMIC_PARTITION )
	{
		
		// true = requires the atomic partition
		t_phi_result r = this->bracketphi( FULL_MASK, true );
		r.round();
		tempfile << "atomic:" << "\tei=" << r.ei << "\tMIPscore=" << r.mip_score << "\t#MIPs=" << r.MIPs.size() << "\tMIPs=" << partitions2str(r.MIPs);

        // if there's more than one atomic partition thing something is deeply wrong.
        assert( r.MIPs.size() == 1 );
        
        tempfile << "\tPHIs=" << PHIs2str(r.PHIs);
        
        
		tempfile << endl;
	}
	
	// Output the subsets with the highest bracketphis, if we're doing that
	if( highestbracketphi )
	{
		vector<t_phi_result> results = highest_bracketphis();

		assert( ! results.empty() );
		
		// foreach subset with the highest minimum unnormalized phi...
		foreach( t_phi_result r, results )	{
			r.round();
			tempfile << "highestbracketphi:" << "\tsubset=" << r.subset.str() << "\tei=" << r.ei << "\tMIPscore=" << r.mip_score << "\t#MIPs=" << r.MIPs.size();
			tempfile << "\tPHIs=" << PHIs2str(r.PHIs) << "\tMIPs=" << partitions2str( r.MIPs );
			tempfile << "\t#s1states=" << num_s1_states( r.subset );
			tempfile << endl;
		}
	}
	
	if( bracketMCs )
	{
		vector<t_phi_result> MC_results = bracket_maincomplexes();
		
		foreach( t_phi_result r, MC_results )	{
			r.round();
			tempfile << "bracketmaincomplex:" << "\tsubset=" << r.subset.str() << "\tei=" << r.ei << "\tMIPscore=" << r.mip_score << "\t#MIPs=" << r.MIPs.size();
			tempfile << "\tPHIs=" << PHIs2str(r.PHIs) << "\tMIPs=" << partitions2str( r.MIPs );
			tempfile << "\t#s1states=" << num_s1_states( r.subset );

			// if we're computing the new_holism of the <MIP>s, do that now.

			
			
			
			tempfile << endl;			
		}
		
	}
		
	// output the subsets with the higest state phis, if we're doing that
	if( higheststatephi )
	{
		for( int i=0; i<NUM_X1_STATES; i++ ) {
			int x1 = this->x1_states[i];
			
			// each of these results will have the same phi_norm
			vector<t_phi_result> results = this->highest_phis( x1 );

			assert( ! results.empty() );
					
			for( int i=0; i<results.size(); i++ )
			{
				t_phi_result r = results[i];
				tempfile << "higheststatephi:" << "\tx1=" << r.x1 << "\tsubset=" << r.subset.str() << "\tei=" << r.ei << "\t#MIPs=" << r.MIPs.size() << "\tPHIs=" << PHIs2str(r.PHIs);
				tempfile << "\tMIPs=" << partitions2str(r.MIPs) << endl;
			}
		}
		
	}

	//output highest <ei>
	if( highestbracketei )
	{
		/*
		t_ei_result r = this->highest_bracket_ei();
		tempfile << "highest_bracket_ei:" << "\tei=" << r.ei << "\tsubsets=";
		for( int i=0; i<(r.subsets.size()-1); i++ ) {
			tempfile << this->subset2str( r.subsets[i] ) << ";";
		}
		// now print the very last subset
		tempfile << this->subset2str( r.subsets.back() ) << endl;
		 */
	}
	
	// output highest state ei
	if( higheststateei )
	{
		/*
		for( int i=0; i<this->x1_states.size(); i++ )
		{
			int x1 = this->x1_states[i];
			
			t_ei_result r = this->highest_ei( x1 );
			tempfile << "highest_state_ei:" << "\tx1=" << r.x1 << "\tei=" << r.ei << "\tsubsets=";
						
			for( int i=0; i<(r.subsets.size()-1); i++ ) {
				tempfile << this->subset2str( r.subsets[i] ) << ";";
			}
			tempfile << this->subset2str( r.subsets.back() ) << endl;
			
		}
		 */
		
	}
	
	
	// the lastline is always "Done!"
	tempfile << "DONE!" << endl;
	tempfile.flush();	//flush it
	tempfile.close();	//close it
	
	///////////////////////////////////////////////////////////////////////////////////////////////////
	// copy from temp_filename to ofilename
	///////////////////////////////////////////////////////////////////////////////////////////////////
	if( rename( tempfilename.c_str(), ofilename.c_str() ) ) {
		perror( NULL );
	}

	ifstream ofile;
	ofile.open( ofilename.c_str() );
	if( ! ofile.is_open() ) {
		cerr << "* Error! -- Could not open filename " << ofilename << endl;
		return -1;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////

	
	cerr << "\t Saved to '" << ofilename << "'";
//	cerr << "\t time: ~" << (((float)clock()-(float)c1))/(float)CLOCKS_PER_SEC << "s";
	cerr << "\t N=" << this->numunits;
	cerr << "\t time: ~" << difftime(time(NULL),t1) << "sec" << endl;	
	
	return 0;
}

vector<bitmask> t_consciousness::all_s1s( const bitmask Smask )
// returns a vector of all possible s1 states given the Smask.
{
	
	vector<bitmask> z;
	z.clear();
	
	//foreach x1 state...
	for( int x1_index=0; x1_index<NUM_X1_STATES; x1_index++ )
	{
		const unsigned int s1 = (x1_states[x1_index]) & Smask;

		// if this s1 IS NOT in z, add it to the list
		bool found_s1_in_z = false;
		
		foreach( unsigned int z_i, z ) {
			if( s1 == z_i ) {
				found_s1_in_z = true;
				break;
			}
		}
		
		if( found_s1_in_z == false )
			z.push_back( s1 );
		
	}
	
//	vec_sort( z );
	
	return z;
}




double t_consciousness::DKL_S0a1b1_and_S0a1_S0b1( const bitmask S0mask, const unsigned int a1, const bitmask A1mask, const unsigned int b1, const bitmask B1mask )
// calculates DKL[ p(X0,a1b1) || p(X0,a1)*p(X0,b1) ]
// which we arrived at using the DKL version of mutual information
// this expression is useful because it's non-negative!!!
{
	assert( S0mask && S0mask <= FULL_MASK );
	assert( is_valid_mu1(a1,A1mask) );
	assert( is_valid_mu1(b1,B1mask) );
	
	const unsigned int a1b1 = a1|b1;
	const bitmask A1B1mask = A1mask|B1mask;
	const unsigned int Ssize = numunits_in_mask( S0mask );
	
	// this doesn't have to be true, but it is for all the cases we're currently looking at.
	assert( S0mask == FULL_MASK );
	
	double summ=0.0;
	
	vector<double> indeps;
	
	//foreach s0 state...
	for( int possible_s0=0; possible_s0<=S0mask; possible_s0++ )
	{
		// if this possible_s0 has any bits ON outside of S0mask, skip it.
		if( (S0mask|possible_s0) > S0mask )
			continue;
		
		const unsigned int s0 = possible_s0;
		
		// the joint distribution
		double joint_distrib = prob_mu0_given_s1( s0, S0mask, a1b1, A1B1mask, Ssize );
		joint_distrib *= prob_s1( a1b1, A1B1mask );
		
		// if this is zero then we're going to get 0.0 * log( 0.0 ), thus we skip.
		if( joint_distrib == 0.0 )
			continue;

		const double prob_s0_given_a1 = prob_mu0_given_s1( s0, S0mask, a1, A1mask, Ssize );
		const double prob_s0_given_b1 = prob_mu0_given_s1( s0, S0mask, b1, B1mask, Ssize );
		
		double indep_distrib = prob_s0_given_a1 * prob_s0_given_b1;
		indep_distrib *= (prob_s1(a1,A1mask) * prob_s1(b1,B1mask));
		
		
		indeps.push_back( indep_distrib );
		
		// sanity check
		assert( 0.0 < indep_distrib );
		
		const double term = log2( joint_distrib / indep_distrib );
								 
		summ += joint_distrib * term;
	}
	
	assert( 0.0 <= summ );
	
	// TODO: check upper-bounds...
	
	return summ;
}




double t_consciousness::H_A0_given_b1( const bitmask A0mask, const unsigned int b1, const bitmask B1mask )
// computes H[ A0 | B1 = b1 ]
// = - \sum_{a0} p(a0|b1) * log p(a0|b1)
{
	assert( is_valid_mu1(b1,B1mask) );
	assert( A0mask && A0mask <= FULL_MASK );
	
	const unsigned int Asize = numunits_in_mask(A0mask);
	
	int count=0;
	double summ = 0.0;
	vector<double> probs;
	
//	cout << "A0mask=" << A0mask << "\t B1mask=" << B1mask << endl;
	
	for( unsigned int possible_a0=0; possible_a0<=A0mask; possible_a0++ )
	{
		// if there is anything ON outside of Amask, skip it
		if( (possible_a0 | A0mask) > A0mask )
			continue;
		
		const unsigned int a0 = possible_a0;
	
		const double this_prob_a0_given_b1 = prob_mu0_given_s1( a0, A0mask, b1, B1mask, Asize );
		
		if( this_prob_a0_given_b1 == 0.0 )
			continue;
		
//		cout << "\t p(a0=" << a0 << "|b1=" << b1 << ")=" << this_prob_a0_given_b1 << endl;
		
		probs.push_back( this_prob_a0_given_b1 );
		count += 1;
		summ -= this_prob_a0_given_b1 * log2( this_prob_a0_given_b1 );
	}
	
//	summ *= -1.0;

	assert( summ == summ );
	assert( 0.0 <= summ );
	
//	cout << "A[A0|B1=" << B1mask << "^b1=" << b1 << "]=" << summ << "\t count=" << count << endl;
//	foreach( double p, probs ) {
//		cout << p << " ";
//	} cout << endl;
	
//	cout << "===================================" << endl;
	
	return summ;
}


double t_consciousness::I_A0_B1_equals_b1( const bitmask Amask, const t_state& b1state )
// alias to I_A0_B1 without the t_state
{
    return I_A0_B1_equals_b1( Amask, b1state.mask, b1state.value );
}

double t_consciousness::I_A0_B1_equals_b1( const bitmask Amask, const bitmask Bmask, const unsigned int b1 )
// calculates the "specific-surprise" DKL, of I( A0 : B1 = b1 )
{
	// assert these are all greater than zero and <= 1
	assert( Amask && Amask <= FULL_MASK );
	assert( Bmask && Bmask <= FULL_MASK );
	assert( is_valid_mu1(b1,Bmask) );
	
	const unsigned int Asize=numunits_in_mask(Amask);
	
	// p(mu0) = 1.0 / 2^{|M|}
	const double prob_a0 = 1.0 / ((double) (1 << Asize));
	
	assert( (0.0 < prob_a0) && (prob_a0 < 1.0) );
	
	double summ=0.0;
	//foreach a0 state...
	for( uint x0=0; x0<=Amask; x0++ )
	{
		// if this x0 has bits ON outside of S0, skip it.
		if( (x0 | Amask) != Amask )
			continue;
		
		const unsigned int a0 = x0 & Amask;
		assert( a0 <= Amask );
		
		//		const double prob_a0_given_b1 = prob_mu0_given_s1( a0, Amask, b1, Bmask );
		const double prob_a0_given_b1 = prob_mu0_given_s1( a0, Amask, b1, Bmask, Asize );		
		
		// if this is 0.0, we don't need to do anymore.  Just skip the next one :)
		if( prob_a0_given_b1 == 0.0 )
			continue;
		
		//		const double this_prob_s1_given_mu0 = prob_s1_given_mu0( s1, Smask, mu0, Mmask );
		
		assert( 0.0 <= prob_a0_given_b1 );
		assert( prob_a0_given_b1 <= 1.0 );
		
		const double term = log2( prob_a0_given_b1 / prob_a0 );
		
		summ += prob_a0_given_b1 * term;
	}
	
	
	assert( 0.0 <= summ );
	
	return summ;	
}


double t_consciousness::perstate_ei( const bitmask S0mask, const unsigned int m1, const bitmask M1mask )
// computes I( S_0 : M_1 = m1 ) using both spec_info and spec_surprise
{
	// validate M1mask and state m1
	assert( S0mask && S0mask <= FULL_MASK );
	assert( is_valid_mu1(m1, M1mask) );
	
	double z = I_A0_B1_equals_b1( S0mask, M1mask, m1 );

	assert( 0.0 <= z );
    assert( z <= H_M0(S0mask) );

	return z;
}

double t_consciousness::perstate_ei( const unsigned int x1 )
// computes I( X_0 : X_1 = x1 ) using both spec_info and spec_surprise.
// This function is just for human convenience as shorthand for perstate_ei setting M1=FULL_MASK
{
	// assert this x1 is valid
	assert( is_valid_mu1( x1, FULL_MASK ) );
	
	return perstate_ei( FULL_MASK, x1, FULL_MASK );	
}




double t_consciousness::specific_info__specinfo__A0_and_b1( const bitmask Amask, const bitmask Bmask, const unsigned int b1 )
// computes the average specific information between part-state b1 and part Amask
// I( A0 : B1 = b1 ) = H[A0] - H[A0|B1=b1]
{
//	assert( FLAG__PIL_SPECIFIC_INFO == "SPEC_INFORMATION" );
	
//	cerr << "Running specific_info" << endl;
	
	// assert these are all greater than zero and <= 1
	assert( Amask && Amask <= FULL_MASK );
	assert( Bmask && Bmask <= FULL_MASK );
	assert( is_valid_mu1(b1,Bmask) );
	const unsigned int Asize = numunits_in_mask( Amask );	
	const double H_A0 = Asize;
	
//	// p(mu0) = 1.0 / 2^{|M|}
//	const double prob_a0 = 1.0 / ((double) (1 << Asize));
//	assert( (0.0 < prob_a0) && (prob_a0 < 1.0) );
	

	//foreach a0 state....
	double summ=0.0;
	for( uint x0=0; x0<=Amask; x0++ )
	{
		// if this x0 has bits ON outside of S0, skip it.
		if( (Amask|x0) > Amask )
			continue;
		
		const unsigned int a0 = x0 & Amask;
		assert( a0 <= Amask );
		
//		const double prob_a0_given_b1 = prob_mu0_given_s1( a0, Amask, b1, Bmask );
		const double prob_a0_given_b1 = prob_mu0_given_s1( a0, Amask, b1, Bmask, Asize );		

		// if this is 0.0, we don't need to do anymore.  Just skip the next one :)
		if( prob_a0_given_b1 == 0.0 )
			continue;
		
		summ -= prob_a0_given_b1 * log2( prob_a0_given_b1 );
	}
	
	double z = H_A0 - summ;
	
	
	if( fequals(z,0.0) )
		z = 0.0;

	assert( 0.0 <= z );
	assert( z <= H_A0 );
	
	return z;
}



bool t_consciousness::is_valid_mu0( const statemask mu0, const bitmask MUmask )
// returns TRUE if the mu0 state is a valid input state of part MUmask
// this function is here to to sanity-check internal calculations
{
	if( (MUmask && MUmask <= FULL_MASK) && ((MUmask|mu0) == MUmask) )
		return TRUE;
	
	return FALSE;
}

bool t_consciousness::is_valid_mu0( const t_state& restrict mu0 )
// returns TRUE if the mu0 state is a valid input state of part MUmask
// this function is here to to sanity-check internal calculations
{
    if( (mu0.mask && mu0.mask <= FULL_MASK) && ((mu0.mask|mu0.value) == mu0.mask) )
        return TRUE;
    
    return FALSE;
}

bool t_consciousness::is_valid_mu1( const statemask mu1, const bitmask MUmask )
// returns TRUE if the mu1 state is a valid output state of part MUmask
// this function is here to to sanity-check internal calculations
{
	// sanity check MUmask and mu1
    assert( MUmask <= FULL_MASK );
	assert( MUmask != 0 );
    assert( MUmask > 0 );

	assert( (MUmask|mu1) == MUmask );
	
	// Foreach output state x1...	
	for( int x1_index=0; x1_index<NUM_X1_STATES; x1_index++ )
	{
		const unsigned int x1 = x1_states[x1_index];
		
		// if this x1 matches s1, then add it's #x0s to s1
		if( (x1&MUmask) == mu1 )
			return TRUE;
	}
	
	return FALSE;
}

bool t_consciousness::is_valid_mu1( const t_state& restrict mu1 )
// returns TRUE if the mu1 state is a valid output state of part MUmask
// this function is here to to sanity-check internal calculations
{
	// sanity check MUmask and mu1
    assert( mu1.mask <= FULL_MASK );
	assert( mu1.mask != 0 );
    assert( mu1.mask > 0 );
    
	assert( (mu1.mask|mu1.value) == mu1.mask );
	
	// Foreach output state x1...	
	for( int x1_index=0; x1_index<NUM_X1_STATES; x1_index++ )
	{
		const unsigned int x1 = x1_states[x1_index];
		
		// if this x1 matches s1, then add it's #x0s to s1
		if( (x1&mu1.mask) == mu1.value )
			return TRUE;
	}
	
	return FALSE;
}

double t_consciousness::DKL_X1a1_from_X1b1( const unsigned int a1, const bitmask A1mask, const unsigned b1, const bitmask B1mask )
// calculates the Kullback-Leibiler diverence D_{KL}[ p(X0|a1) || p(X0|b1) ] 
// A1mask and B1mask *CAN* (and typically do) overlap
{
	// assert both masks are valid
	assert( A1mask && A1mask <= FULL_MASK );
	assert( B1mask && B1mask <= FULL_MASK );
	
	// assert a1,b1 is a subset of A1mask, B1mask
	assert( is_valid_mu1(a1,A1mask) );
	assert( is_valid_mu1(b1,B1mask) );
	
	double summ=0.0;
	
	// Foreach state x1...	
	for( int x1_index=0; x1_index<NUM_X1_STATES; x1_index++ )
	{
		const unsigned int x1 = x1_states[x1_index];
		
		const double prob_x1_given_a1 = prob_a1_given_mu1( x1, FULL_MASK, a1, A1mask );

		if( prob_x1_given_a1 == 0.0 )
			continue;
		
		const double prob_x1_given_b1 = prob_a1_given_mu1( x1, FULL_MASK, b1, B1mask );
		
		// if this is ever zero with a nonzero prob_x1_given_a1, then the DKL isn't defined.
		assert( 0.0 < prob_x1_given_b1 );

		const double term = log2( prob_x1_given_a1 / prob_x1_given_b1 );
		
		summ += prob_x1_given_a1 * term;		
	}
	
	return summ;
}

double t_consciousness::DKL_S0a1_from_S0b1( const bitmask S0mask, const unsigned int a1, const bitmask A1mask, const unsigned b1, const bitmask B1mask )
// calculates the Kullback-Leibiler diverence D_{KL}[ p(S0|a1) || p(S0|b1) ] 
// A1mask and B1mask don't have to be disjoint.  In fact typically aren't.
{
	// assert a1,b1 are valid output states
	assert( is_valid_mu1(a1,A1mask) );
	assert( is_valid_mu1(b1,B1mask) );
	assert( S0mask && S0mask <= FULL_MASK );
	
	const unsigned int Ssize = numunits_in_mask(S0mask);
	
	double summ=0.0;
	//foreach x0 state...
	for( unsigned int possible_s0=0; possible_s0<=S0mask; possible_s0++ )
	{
		// if this possible_s0 has any bits ON outside of S0mask, skip it
		if( (S0mask|possible_s0) > S0mask )
			continue;
		
		const unsigned int s0 = possible_s0;
		
		const double prob_s0_given_a1 = prob_mu0_given_s1( s0, S0mask, a1, A1mask, Ssize );
		
		// if this is zero then we can skip the next!
		if( prob_s0_given_a1 == 0.0 )
			continue;
		
		const double prob_s0_given_b1 = prob_mu0_given_s1( s0, S0mask, b1, B1mask, Ssize );
		
		// if this is ever zero with a nonzero prob_x0_given_a1, then the DKL isn't defined.
		assert( 0.0 < prob_s0_given_b1 );
		
		const double term = log2( prob_s0_given_a1 / prob_s0_given_b1 );
		
		summ += prob_s0_given_a1 * term;
	}
	
	assert( 0.0 <= summ );
	
	return summ;
}

double t_consciousness::DKL_X0a1_from_X0b1( const unsigned int a1, const bitmask A1mask, const unsigned b1, const bitmask B1mask )
// calculates the Kullback-Leibiler diverence D_{KL}[ p(X0|a1) || p(X0|b1) ] 
// A1mask and B1mask don't have to be disjoint.  In fact typically aren't.
{
	return DKL_S0a1_from_S0b1( FULL_MASK, a1, A1mask, b1, B1mask );
}


double t_consciousness::prob_a1_given_mu1( const unsigned int a1, const bitmask Amask, const unsigned int mu1, const bitmask MUmask )
// calculates p( a1 | mu1 ) with both parts at timestep 1
// note that A1mask and MUmask *CAN* overlap.
{
	// verify masks are subsets of the network
	assert( Amask <= FULL_MASK );
	assert( MUmask <= FULL_MASK );
	
	// assert that states are subsets of their parts
	assert( is_valid_mu1(a1,Amask) );
	assert( is_valid_mu1(mu1,MUmask) );

	unsigned int mu1_sightings=0, a1mu1_sightings=0;

	// Foreach output state x1...
	for( int x1_index=0; x1_index<NUM_X1_STATES; x1_index++ )
	{
		const unsigned int x1 = x1_states[x1_index];
		
		// if this x1 matches mu1, then add it's #mu1_sightings
		if( (x1&MUmask) == mu1 )
		{
//			mu1_sightings += (*it).second;
			mu1_sightings += num_x0s_to_this_x1[x1_index];

			// if this x1 matches mu1 AND a1, then add it's #a1mu1_sightings
			if( (x1&Amask) == a1 )
				a1mu1_sightings += num_x0s_to_this_x1[x1_index];
		}		
	}
	
	if( a1mu1_sightings == 0 )
		return 0.0;

	// p(a1|mu1) = p( a1,mu1 ) / p(mu1)
	double z = (double) a1mu1_sightings / (double) mu1_sightings;
	
	// this should be 0.0 < because we already checked for the 0.0 case
	assert( 0.0 < z );
	assert( z <= 1.0 );
	
	return z;
}

// alias to t_consciousness::binary_repr()
string t_consciousness::binrep( const unsigned int digit, unsigned int numplaces ) { return binary_repr( digit, numplaces ); }

string t_consciousness::binary_repr( const unsigned int digit, unsigned int numplaces )
{
	// if numplaces isn't defined, assume whole network
	if( numplaces == 0 )
		numplaces = this->numunits;
	
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

/*

void t_consciousness::MAPREDUCE_WORKER_MIPs( unsigned int worker_index, unsigned int number_of_workers, string generator_cmd )
// This function tries all x1 states, attempting every worker_index'th partition read from STDIN. 
{

	cerr << "Needs to be updated." << endl;
	exit(-1);
	
	// foreach x1 state...
	int SIZE=512;
	char buf[SIZE];	
	for( map<int,unsigned int>::iterator it=this->existing_states.begin(); it != this->existing_states.end(); it++ )	{
		int x1 = (*it).first;
		double prob_x1 = float((*it).second) / this->numstates;
		
		double min_ei=this->numunits, min_norm=1.0;
		double ei, norm=1.0;
		vector<vector<vector<int> > > min_partitions;
		
		vector<vector<int> > partition;
		partition.clear();	
		string line;
		unsigned int counter=0;
		
		FILE * my_pipe = this->get_partitions( 0, this->max_num_parts );
		
		// now read the parts from STDIN...
		while( fgets(buf, SIZE, my_pipe) != NULL )
		{
			counter++;
			line = string(buf);
			//			cout << "myline=" << line << endl;
			
			// Is this worker assigned to this partition?  If so, attempt it.
			if( counter % number_of_workers == worker_index ) 
			{
				
				//				cout << counter << " % " << number_of_workers << " = " << worker_index << "\ttrying parts=" << line << endl;				
				
				this->str2partition( line, &partition );
				ei = this->ei( x1, &partition );
				// if we are normalizing, do that now.						
				norm = this->normalization( &partition );


				// if this partition has a LOWER normal EI than the current MIP: clear previous MIP.  Set partition as the MIP.
				if( (ei/norm) < (min_ei/min_norm) ) {
					min_partitions.clear();
					min_partitions.push_back( partition );
					min_ei = ei;
					min_norm = norm;
					
					//					cout << "\tmin-partition:" << this->parts2str( partition ) << endl;
				}
				
				// if the current partition has THE SAME normalized EI as the current MIP
				else if( (ei/norm) == (min_ei/min_norm) ) {
					min_partitions.push_back( partition );
				}
			}
			//			else
			//				cout << counter << " % " << number_of_workers << " = " << counter % number_of_workers << "\t NOT trying parts=" << line << endl;							
		}
		
		pclose(my_pipe);
		
		
		// FINISHED FOR THIS x1 state -- PRINT IT.
		cout << "x1=" << x1 << "\tMIPs=";
		
		for( int i=0; i< min_partitions.size(); i++ ) {
			cout << this->partition2str( &min_partitions[i] );
			
			if( i < min_partitions.size()-1 )
				cout << ";";	//partition separator
		}
		
		cout << "|min_ei=" << min_ei;
		cout << "|normed_min_ei=" << min_ei/min_norm;
		cout << "|prob(x1)=" << prob_x1;
		cout << "|ei(x1)=" << this->ei(x1);
		cout << endl;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Compute the average EI min partition for the average effective information now.
	//////////////////////////////////////////////////////////////////////////////////////////////////
	
	double min_ei=this->numunits, min_norm=1.0;
	double ei, norm=1.0;
	vector<vector<vector<int> > > min_partitions;
	
	vector<vector<int> > partition;
	partition.clear();
	string line;
	unsigned int counter=0;
	
	FILE * my_pipe = this->get_partitions( 0, this->max_num_parts );
	
	// now read the parts from STDIN...
	while( fgets(buf, SIZE, my_pipe) != NULL )
	{
		counter++;
		line = string(buf);
		//			cout << "myline=" << line << endl;
		
		// Is this worker assigned to this partition?  If so, attempt it.
		if( counter % number_of_workers == worker_index ) 
		{
			
			//				cout << counter << " % " << number_of_workers << " = " << worker_index << "\ttrying parts=" << line << endl;				
			this->str2partition( line, &partition );
			ei = this->bracket_ei( &partition );
			// if we are normalizing, do that now.					
			norm = (double) this->normalization( &partition );
			
			// if this partition has a LOWER normal EI than the current MIP: clear previous MIP.  Set partition as the MIP.
			if( (ei/norm) < (min_ei/min_norm) ) {
				min_partitions.clear();
				min_partitions.push_back( partition );
				min_ei = ei;
				min_norm = norm;
			}
			//					cout << "\tmin-partition:" << this->parts2str( partition ) << endl;
			
			// if the current partition has THE SAME normalized EI as the current MIP
			else if( (ei/norm) == (min_ei/min_norm) ) {
				min_partitions.push_back( partition );
			}
		}
	}
	
	pclose(my_pipe);
	
	// FINISHED FOR THIS x1 state -- PRINT IT.
	cout << "x1=avr" << "\tMIPs=";
	
	for( int i=0; i< min_partitions.size(); i++ ) {
		cout << this->partition2str( &min_partitions[i] );
		
		if( i < min_partitions.size()-1 )
			cout << ";";	//partition separator
	}
	
	cout << "|min_ei=" << min_ei;
	cout << "|normed_min_ei=" << min_ei/min_norm;
	cout << "|prob(x1)=0.0";
	cout << "|ei(x1)=" << this->bracket_ei();
	cout << endl;
	
	
}	// end the worker function
*/


// returns the number of entries in the array
void subset2array( unsigned int subset, unsigned int* restrict z )
{

	unsigned int num_ON=0, power=0;
	
	
	do	{
		// this only works because we're starting at power=0 and num_ON=0
		if( subset & 1 ) {
			z[num_ON] = power;
			num_ON += 1;
		}
		
		subset >>= 1;
		power += 1;
	}while( subset );

}

inline bool fequals( const double a, const double b )
// if and a b are within PRECISION of each other, return true -- else false.
{

	return( fabs(a-b) <= PRECISION );
	
//	if( fabs(a-b) <= PRECISION ) {
//		ASSERT( fabs(b-a) <= PRECISION );
//		return true;
//	}
	
	
//	return false;
}


double assert_bounds( double lower, double value, double upper )
// check value against the specified upper and lower bounds.  Being mindful of float-point error at equality.
{
    // switch the bounds incase we need to.
    // this also takes care of floating-point error on the bounds
    lower = min(lower,upper);
    upper = max(lower,upper);
    
    // this should be true
    assert( lower <= upper );
    
    // address lower floating-point error
    if( fequals(lower,value) )
        value = max(lower,value);

    // address upper floating-point error    
    if( fequals(value,upper) )
        value = min(value,upper);
    
    // these should be true
    assert( lower <= value );
    assert( value <= upper );
    
    return value;
}


// floors a double to DECIMAL_PRECISION places
//long double myceil( long double x ) { return (1/pow( (long double) 10.0,DECIMAL_PRECISION-1)) * ceil( x * pow( (long double) 10.0,DECIMAL_PRECISION-1) ); }
long double myceil( double x ) { return (1/pow( (double) 10.0, (double) DECIMAL_PRECISION-1)) * ceil( x * pow( (double) 10.0, (double) DECIMAL_PRECISION-1) ); }

// ceils a double to DECIMAL_PRECISION places
//long double myfloor( long double x ) { return (1/pow( (long double) 10.0,DECIMAL_PRECISION-1)) * floor( x * pow( (long double) 10.0,DECIMAL_PRECISION-1) ); }
long double myfloor( double x ) { return (1/pow( (double) 10.0, (double) DECIMAL_PRECISION-1)) * floor( x * pow( (double) 10.0, (double) DECIMAL_PRECISION-1) ); }

// rounds a double to DECIMAL_PRECISION places
double myround( double x, double precision )
{
	if( precision==0.0 )
		precision = DECIMAL_PRECISION;
	
	return (1/pow( (double) 10.0,precision)) * round( x * pow( (double) 10.0,precision) );
//	return x;
}

string str_rstrip( string haystack, const string needles )
{
	// loop until we don't find a needle at the end
	while( 1 )
	{
		bool found_a_needle = false;
		
		for( int i=0; i<needles.length(); i++ ) {
			char temp_char = needles[i];
			
			if( temp_char == haystack[ haystack.length()-1 ] ) {
				haystack.erase( haystack.length()-1 );
				found_a_needle = true;
			}
			
		}
		
		
		if( found_a_needle==false )
			break;
		
	}
	
	return haystack;
}

string str_istrip( string haystack, string needle )
// returns the string haystack with string 'needle' removed from the beginning (if it exists at all)
{
	while( str_istartswith( haystack, needle ) ) {
		haystack.erase( 0, needle.length() );
	}
	
	return haystack;
}


bool str_istartswith( string haystack, string needle )
// returns TRUE if string haystrack starts with needle
// else returns FALSE
{
	for( int i=0; i<haystack.length(); i++ )
		haystack[i] = tolower(haystack[i]);
	
	for( int i=0; i<needle.length(); i++ )
		needle[i] = tolower(needle[i]);
	
	return str_startswith( haystack, needle );
}



void binary2partition( FILE* f, t_partition* P )
// reads from filestream f to return partition
{
	//it's IMPORTANT that numparts be initially set to ZERO
	unsigned int numparts=0;
	
	fread( &numparts, sizeof(unsigned int), 1, f );
	
	if( numparts == 0 ) {
		ASSERT( feof(f) == true );
		return;
	}
	
	ASSERT( 1 <= numparts );
	
	// 1. read the mask of each part into the masks array
	unsigned int masks[numparts];
	fread( &masks, sizeof(unsigned int), numparts, f );
	
	// new version using t_partition
	P->update( numparts, masks );
	
}


unsigned int mask2array( unsigned int inpmask, unsigned int* restrict z )
{
	//	unsigned int power = (unsigned int) floor( log2( subset ) );
	unsigned int power=0, num_ON=0;
	//	unsigned int temp_S;
	
	//	// if outside is TRUE, get the nodes that are OUTSIDE the subset
	//	if( outside == true )
	//		temp_S ^= max_S;
	
	do	{
		// this only works because we're starting at power=0 and num_ON=0
		if( inpmask&1 ) {
			z[num_ON] = power;
			num_ON += 1;
		}
		
		inpmask >>= 1;
		power += 1;
	}while( inpmask );
	
    
    return num_ON;
}

#endif
