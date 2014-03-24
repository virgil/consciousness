/////////////////////////////////////////////////////////////////////////////////////////////////////////
// helpers.h : A bunch of handy helper functions.  Some of them are nessecary for t_consc to work
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma once
#ifndef HELPERS_H
#define HELPERS_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <queue>
#include <math.h>
#include "t_partition.h"
#include <algorithm>
#include <boost/foreach.hpp>
#include <sstream>
#include <limits>

using namespace std;

#ifndef foreach
	#define foreach BOOST_FOREACH
#endif

#ifndef reverse_foreach
	#define reverse_foreach BOOST_REVERSE_FOREACH
#endif

#define restrict __restrict__

#ifndef assertmsg
	#define assertmsg(exp, msg) if((exp)==false) { cerr << endl << "* ASSERT ERROR: '" << msg << "'" << endl; cout.flush(); cerr.flush(); assert(exp); }
#endif

typedef unsigned int uint;
typedef unsigned long ulong;

// improves code clarity and makes easier to
// change to unsigned long's should we need to.
typedef unsigned int statemask;
typedef unsigned int partmask;
typedef unsigned int bitmask;


// define 'TRUE/True' and 'FALSE/False' (C++ style!)
const bool FALSE=false;
const bool TRUE=true;

// default values for unsigned int's and int's 
const unsigned int UINT_NONE = numeric_limits<unsigned int>::max(); // UINT_MAX;
const int INT_NONE = numeric_limits<int>::max(); //INT_MIN;
const double DBL_NONE = numeric_limits<double>::quiet_NaN();
const float FLT_NONE = numeric_limits<float>::quiet_NaN();
const double DBL_NaN = numeric_limits<double>::quiet_NaN();
const float FLT_NaN = numeric_limits<float>::quiet_NaN();

const double DBL_MAX = numeric_limits<double>::max();
const double DBL_MIN = numeric_limits<double>::min();


// returns the sign of the passed double
int SIGN( const double a );

// returns which double is "simpler"
double simpler_double( double x1, double x2 );

// calculate the int( log2(inp) ) where inp is an 32-bit unsigned integer AS QUICKLY AS POSSIBLE
//unsigned int fast_log2_v2( uint32_t inp );
unsigned int fast_log2( const uint32_t inp );

// reverse the elements in an array of length 'len'
void reverse_array(unsigned int* inver_a, const unsigned int len);


// converts a string of binary characters into an unsigned integer
unsigned int bin2dec(char* bin);

// compute N choose K
unsigned int NchooseK( const unsigned int n, const unsigned int k );

double ZERO( vector<double> inps );
double ZERO( const double a, const double b );

void printc(int comb[], int k);
int next_combination(unsigned int comb[], const unsigned int k, const unsigned int n);


string partitions2str( const vector< t_partition > partitions );
string PHIs2str( const vector< double > PHIs, double ei=0.0 );

string str_istrip( string haystack, string needle );
bool str_startswith( string haystack, string needle );
bool str_istartswith( string haystack, string needle );

string tolower( string z );
string toupper( string z );


string str_rstrip( string haystack, const string needles="\r\n" );
string str_strip( string haystack, string needle );
string str_istrip( string haystack, string needle );
string binary_repr( unsigned int number, int padding=1 );

///////////////////////////////////////////////////////////////////////
// FUNCTION PROTOTYPES
template <class T> void vec_sort( vector<T>& inp );
template <class T> T sum( vector<T> inp );
template <class T> T sum( vector<vector<T> > inp );
template <class T> T sum( const T* inp, unsigned int n );
///////////////////////////////////////////////////////////////////////


template <class T> string vec_str( const vector<T>& restrict inp );
/*
// return the vector as a string
{
	string z = "[ ";

	stringstream ss (stringstream::in | stringstream::out);
	
	foreach( T ele, inp ) {
//		z.append( string(ele) );

		ss << ele;
		z.append( ss.str() );
		
		z.append(" ");
	}

//	z = str_rstrip( z, ", " );

	//	// remove the final two characters
//	z.erase( z.end()-2, z.end() );
//	z.erase( z.end() );
	
//	z.append(" ]");
	z.append("]");

	return z;
}
*/

template <class T> void vec_sort( vector<T>& restrict inp )
// sort the passed vector
{
	sort( inp.begin(), inp.end() );	
}


template <class T> void vec_remove( vector<T>& restrict haystack, const T& restrict needle )
// remove a needle from the vector haystack
{
	for( int i=0; i<haystack.size(); i++ ) {
		if( haystack[i] == needle )
			haystack.erase( i );
	}
	
}

template <class T> void vec_extend( vector<T>& restrict v1, const vector<T>& restrict v2 )
// remove multiple needles from the vector haystack
{
	v1.reserve( v1.size() + v2.size() );
	
	foreach( T ele, v2 )
        v1.push_back( ele );
	
}



template <class T> void vec_remove( vector<T>& haystack, vector<T>& needles )
// remove multiple needles from the vector haystack
{
	// sort
	vec_sort( needles );
	vec_sort( haystack );
	
	int count=0;
	foreach( T needle, needles ) {
		
		
		for( ; count<haystack.size(); count++ ) {
			if( haystack[count] == needle )
				haystack.remove( count );
			
			// goto the next needle if the haystack has exceeded this needle value
			else if( haystack[count] > needle )
				break;
		}
		
	}
	
}


template <class T> void vec_unique( vector<T>& restrict inp )
// uniqueify the values in the vector
{
	// you MUST DO THIS SORT HERE -- unique() only removes CONSECUTIVE SAME VALUES
	vec_sort( inp );
	inp.resize( unique(inp.begin(), inp.end()) - inp.begin() );		//uniquify	
}

template <class T> T max( const vector<T>& restrict v )
// return the highest value in the vector
{
	return *max_element(v.begin(), v.end() );	
}

template <class T> T min( const vector<T>& restrict v )
// return the highest value in the vector
{
	return *min_element(v.begin(), v.end() );
	
}


template <class T> bool vec_in( const T& restrict needle, const vector<T>& restrict haystack )
// returns true if the needle is in the vector
{
	for( int i=0; i<haystack.size(); i++ ) {
		if( haystack[i] == needle )
			return true;
	}
	
	return false;
}


// alias to the other way
template <class T> bool vec_in( const vector<T>& restrict haystack, const T& needle ) { return vec_in( needle, haystack); }

template <class T> T sum( const T* inp, unsigned int n )
{
	T z = 0.0;
	
	for( int i=0; i<n; i++ )
		z += inp[i];
	
	return z;
}

template <class T> T sum( vector<T> inp )
{	
	T z = 0.0;
	
	for( int i=0; i<inp.size(); i++ )
		z += inp[i];
	
	return z;
}

template <class T> T sum( vector<vector<T> > inp )
{
	T z = 0.0;
	
	for( int i=0; i<inp.size(); i++ ) {
		for( int j=0; j<inp[i].size(); j++ )
			z += inp[i][j];
	}
	
	return z;
}


#endif