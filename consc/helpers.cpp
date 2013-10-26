/////////////////////////////////////////////////////////////////////////////////////////////////////////
// helpers.cpp: CPP file for the handy helper functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef HELPERS_CPP
#define HELPERS_CPP

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <queue>
#include <math.h>
#include <algorithm>
#include <boost/foreach.hpp>
#include <sstream>
#include "t_partition.h"
#include "helpers.h"
#include <limits>


#ifndef foreach
	#define foreach BOOST_FOREACH
#endif

//extern const unsigned int MAXIMUM_NUMBER_OF_MAIN_COMPLEXES_TO_STORE;



int SIGN( const double a )
{
    if( 0.0 < a )
        return 1;
    if( a < 0.0 )
        return -1;
    else
        return 0.0;
}


double ZERO( const double a, const double b )
{
    const double adist = fabs(a);
    const double bdist = fabs(b);
    
    if( adist <= bdist )
        return a;
    if( bdist <= adist )
        return b;
    else
        assert( 0 == 1 );   // shouldn't ever get here.
}

double ZERO( vector<double> inps )
// return the element of inps that is closest to ZERO (can be positive or negative).
{
    double z = DBL_MAX;
    
    foreach( double inp, inps ) {
//        cout << "inp=" << inp << endl;
//        cout << "z=" << z << endl;
        const double this_dist_from_zero = fabs(inp);
        const double closest_so_far = fabs(z);
        
        if( this_dist_from_zero < closest_so_far )
            z = inp;
    }
    
//    cout << "* z =" << z << endl << endl;
    return z;
}

unsigned int fast_log2_v2( uint32_t inp )
// calculate the int( log2(inp) ) where inp is an 32-bit unsigned integer AS QUICKLY AS POSSIBLE
// !!IMPORTANT!!  THIS METHOD ASSUMES inp is one less than a power of TWO!!
// Algorithm taken from: http://www-graphics.stanford.edu/~seander/bithacks.html
{
	assertmsg( 0 == 1, "This one wasn't as fast as fast_log2()" );
	
	// input is ONE LESS than a POWER OF TWO.
	// but algorithm requires inp be a POWER OF TWO.
	// Thus we increment inp, and decrement the result, z.
	inp += 1;

	// fancy way to check that inp is a power of two
	assert( ((inp) % 2) == 0 );	
	assert( (inp && !(inp & (inp - 1))) == true);			
	
	static const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 
		0xFF00FF00, 0xFFFF0000};

	register unsigned int z = (inp & b[0]) != 0;
	
//	for(int i = 4; i > 0; i--) // unroll for speed...
//		z |= ((inp & b[i]) != 0) << i;
	
	z |= ((inp & b[4]) != 0 ) << 4;
	z |= ((inp & b[3]) != 0 ) << 3;
	z |= ((inp & b[2]) != 0 ) << 2;
	z |= ((inp & b[1]) != 0 ) << 1;
		

	return (z-1);
}



inline unsigned int fast_log2( const uint32_t inp )
// calculate the int( log2(inp) ) where inp is an 32-bit unsigned integer AS QUICKLY AS POSSIBLE
// !!IMPORTANT!!  THIS METHOD ASSUMES inp is one less than a power of TWO!!
// This algorithm is the fastest one thus far
// Algorithm taken from: http://www-graphics.stanford.edu/~seander/bithacks.html
{
	static const int MultiplyDeBruijnBitPosition[32] = 
	{
		0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
		8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
	};	

//	const uint32_t inp_plus_one = inp + 1;		
//	// fancy way of checking that inp_plus_one is a power of 2
//	assert( ((inp_plus_one) % 2) == 0 );
//	assert( (inp_plus_one && !(inp_plus_one & (inp_plus_one - 1))) == true);		

	return MultiplyDeBruijnBitPosition[(uint32_t)(inp * 0x07C4ACDDU) >> 27];	
}

/*
void reverse_array(unsigned int* inver_a, const unsigned int size)
{	
	for(int i=0; i<(size/2); i++)
	{
		const unsigned int temp=inver_a[i];
		inver_a[i]=inver_a[size-i-1];
		inver_a[size-i-1]=temp;
	}
}
*/

double simpler_double( double x, double y )
// for two equal doubles, returns the double that is "simpler"
// the double that is closest to a whole number
{
    // assert that these two floats ARE NOT the same.
    // assert that the two floats ARE the same within 1 decimal point

    assert( fabs(x-y) <= 0.001 );
    
    // get the distance from the nearest whole number
    double x_dist = fabs( x - round(x) );
    double y_dist = fabs( y - round(y) );
    
    if( x_dist <= y_dist )
        // x_dist is smaller
        return x;
    else 
        // y_dist is smaller
        return y;
    
}

////////////////////////////////////////////////////////////////////////////////////////////////
//	functions for dealing with combinations.
//      * printc()              -   prints our a combination stored as an array
//      * next_combinat()       -   stores the next combination into comb[]
////////////////////////////////////////////////////////////////////////////////////////////////


// Prints out a combination like {1, 2}
void printc(int comb[], int k) {
	printf("{");
	for(int i = 0; i < k; ++i)
		printf("%d, ", comb[i] );
	printf("\b\b}\n");
}



////////////////////////////////////////////////////////////////////////////////////////////////
//	next_comb(int comb[], int k, int n)
//		Generates the next combination of n elements as k after comb
//
//	comb => the previous combination ( use (0, 1, 2, ..., k) for first)
//	k => the size of the subsets to generate
//	n => the size of the original set
//
//	Returns: 1 if a valid combination was found
//		0, otherwise
////////////////////////////////////////////////////////////////////////////////////////////////



int next_combination( unsigned int comb[], const unsigned int k, const unsigned int n) {
	int i = k - 1;
	++comb[i];
    //	while ((i >= 0) && (comb[i] >= n - k + 1 + i)) {
	while ((i > 0) && (comb[i] >= n - k + 1 + i)) {
		--i;
		++comb[i];
	}
    
	if (comb[0] > n - k) // Combination (n-k, n-k+1, ..., n) reached
		return 0; // No more combinations can be generated
    
	// comb now looks like (..., x, n, n, n, ..., n).
	// Turn it into (..., x, x + 1, x + 2, ...) 
	for (i = i + 1; i < k; ++i)
		comb[i] = comb[i - 1] + 1;
    
	return 1;
}



LatticeAtoms lattice_nodes_below( const t_partition& restrict P, const bool P_timestep )
// this function returns all of the nodes below the LatticeNode "{P_1, ... P_K}" for a given P
{
	///////////////////////////////////////////////////////////////////////////
	// if this isn't true then there will will be NO lattice nodes beneath it
	assert( 2 <= P.size() );

	assert( P_timestep == 0 );
	///////////////////////////////////////////////////////////////////////////
	
	
	// we're going to have exactly P.size() of these.
	LatticeAtoms z;
	z.reserve( P.size() );
	
	// creates |P| nodes, each containing a SINGLE atom
	for( int i=0; i<P.size(); i++ )
	{
		// atom_mask contains everything except for the i'th part
		const unsigned int atom_mask = P.subset & (~P[i]);
		
		// make the atom
		LatticeAtom new_atom( atom_mask, P_timestep );

		// make the node, and push the atom onto it
//		LatticeNode new_node;
//		new_node.push_back( new_atom );		
		
		z.push_back( new_atom );
	}
	
	assert( z.size() == P.size() );	
/*	
	// if there is an exclusion mask, add the mask of itself to every node
	if( exclusion_mask )
	{
		// make an atom for the whole partition
		LatticeAtom partition_atom( P.subset, Ptimestep );
		
		
		foreach( LatticeNode this_node, z )
		{
			// if this isn't true something is amiss
			assert( this_node.size() == 1 );
			
			if( exclusion_timestep )
				this_node[0].X1mask |= exclusion_mask;
			else
				this_node[0].X0mask |= exclusion_mask;
			
			
			this_node.push_back( partition_atom );
			
			// if this isn't true something is wrong
			assert( this_node.size() == 2 );
		}
		
	}
*/
	
	assert( z.size() == P.size() );
	return z;
}

unsigned int bin2dec(char* bin)
// converts a string of binary characters into an unsigned integer
{
	int b, n;
	unsigned int z=0; 
	
	const unsigned int len = strlen(bin) - 1;
	for(int k = 0; k<=len; k++) 
	{
		b = 1;
		n = (bin[k] - '0'); // char to numeric value
		if ((n > 1) || (n < 0)) 
		{
			puts("\n\n ERROR! BINARY has only 1 and 0!\n");
			return (0);
		}
		b = b<<(len-k);
		// sum it up
		z += n * b;
		//printf("%d*%d + ",n,b);  // uncomment to show the way this works
	}
	
	return z;
}

unsigned int NchooseK( const unsigned int n, const unsigned int k )
// compute N choose K really fast
{
	assert( k <= n );
	
	vector<unsigned int> b( n+1 );
	
	b[0] = 1;
	for( unsigned int i=0; i<=n; i++ )
	{
		b[i] = 1;
		for( unsigned int j=1-1U; j>0; j-- )
			b[j] += b[ j-1U ];
	}
	
	return b[k];	
}

string PHIs2str( const vector<double> PHIs, double ei )
// converts a vector<double> to a string for printing purposes
// if ei is not zero, then print the normalized value in parens.
{
    assert( ei >= 0.0 );
    
    
	string z = "";
	char buf[64];
	
	for( int i=0; i<PHIs.size(); i++ )
	{
		sprintf(buf, "%.4f", PHIs[i] );
		z.append( buf );
        
        // if ei > 0.0, then print the normalized value.
        if( ei > 0.0 ) {
            sprintf(buf, " (%.4f)", PHIs[i] / ei );
            z.append( buf );
        }
		
		if( i<PHIs.size()-1 )
			z.append(";");
	}
	
	return z;	
}

string str_strip( string haystack, string needle )
// returns the string haystack with string 'needle' removed from the beginning (if it exists at all)
{
	while( str_startswith( haystack, needle ) )
		haystack.erase( 0, needle.length() );
	
	return haystack;
}


bool str_startswith( string haystack, string needle )
// returns TRUE if string haystrack starts with needle
// else returns FALSE
{
	int nlength = needle.length();
	int hlength = haystack.length();
	if( hlength < nlength )
		return false;
	
	if( haystack.substr(0,nlength) == needle )
		return true;
	else
		return false;
}

string tolower( string z )
{	
	for( int i=0; i<z.length(); i++ ) {
		z[i] = tolower(z[i]);
	}	
	return z;
}

string toupper( string z )
{	
	for( int i=0; i<z.length(); i++ ) {
		z[i] = toupper(z[i]);
	}	
	return z;
}


string binary_repr( unsigned int digit, int numplaces )
{	
	string z="";
	for( int i=0; i<numplaces; i++ )
	{
		int bit = ((digit >> i)&1);
		
		if( bit )
			z.append( "1" );
		else
			z.append( "0" );
	}
	
	return z;
}


string partitions2str( const vector<t_partition> Ps )
// convert a list of partitions into a string
{
	string z="";
	
	for( int i=0; i<Ps.size(); i++ )
	{
		z.append( Ps[i].BAREstr() );
		
		if( i < (Ps.size()-1) )
			z.append( string(";") );
		
	}
	
	return z;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


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
