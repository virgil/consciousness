#ifndef T_SUBSET_CPP
#define T_SUBSET_CPP

#include <assert.h>
#include <string>
#include <iostream>
#include <sstream>
#include <list>
#include <boost/foreach.hpp>
#include <stdio.h>

#include "t_subset.h"

#ifndef DO_CAPITALIZED_ASSERTS
	#define DO_CAPITALIZED_ASSERTS 0
#endif

#if DO_CAPITALIZED_ASSERTS
	#define ASSERT(exp) assert(exp)
	#define ASSERTMSG(exp, msg) assertmsg(exp,msg)
#else
	#define ASSERT(exp)				//defined as nothing
	#define ASSERTMSG(exp, msg)		//defined as nothing
#endif

#ifndef foreach
	#define foreach BOOST_FOREACH
#endif

#define restrict __restrict__

using namespace std;

const unsigned int MAX_NUMUNITS = 31;
const unsigned int DEFAULT_SUBSET = 0;

void t_subset::update_num_nodes()
{
	this->num_nodes = 0;
	
	unsigned int tempS = this->S;
	
	do	{
		if( tempS&1 )
			this->num_nodes += 1;
		
		tempS >>= 1;
	}while( tempS );
	
}

t_subset::t_subset( unsigned int inp_numunits, unsigned int inp_S )
// the constructor
{
	if( inp_numunits > MAX_NUMUNITS && inp_S == DEFAULT_SUBSET ) {
		inp_S = inp_numunits;
		inp_numunits = MAX_NUMUNITS;
	}
	
//	cout << "running constructor..." << endl;
	this->max_numnodes = inp_numunits;
	this->max_S = (1 << max_numnodes)-1;
	this->S = inp_S;

	this->num_nodes = 0;
	this->update_num_nodes();
	
	assert( num_nodes <= max_numnodes );
	
//	cout << "inp_numunits=" << inp_numunits << endl;
//	cout << "INIT: S=" << S << "\tmax_S=" << max_S << endl;
//	cout.flush();
	
	assert( 1 <= max_numnodes );
	assert( max_numnodes <= MAX_NUMUNITS );
	assert( S <= max_S );
}

t_subset::t_subset()
// constructor for no arguments
{
	this->max_numnodes = MAX_NUMUNITS;
	this->max_S = (1 << max_numnodes)-1;
	this->S = DEFAULT_SUBSET;
	update_num_nodes();
}

void t_subset::array( unsigned int* restrict z, bool outside ) const
{
	//	unsigned int power = (unsigned int) floor( log2( subset ) );
	unsigned int power=0, num_ON=0;
	unsigned int temp_S=S;

	// if outside is TRUE, get the nodes that are OUTSIDE the subset
	if( outside == true )
		temp_S ^= max_S;
	
	do	{
		// this only works because we're starting at power=0 and num_ON=0
		if( temp_S&1 ) {
			z[num_ON] = power;
			num_ON += 1;
		}
		
		temp_S >>= 1;
		power += 1;
	}while( temp_S );

}

vector<unsigned int> t_subset::nodes(bool outside) const
// returns all of the bits that are turned ON in the subset
{
	// clear the old part
	vector<unsigned int> z;
	z.clear();
	
	//	unsigned int power = (unsigned int) floor( log2( subset ) );
	unsigned int power=0, temp_S=S;

	// if outside is TRUE, get the nodes that are OUTSIDE the subset
	if( outside == true )
			temp_S ^= max_S;
	
	do	{
		// this only works because we're starting at i=0
		if( (temp_S&1) )
			z.push_back(power);
		
		temp_S >>= 1;
		power += 1;
	}while( temp_S );
		
	return z;
}

vector<unsigned int> t_subset::masks() const
// return an vector of bitmasks --- one bitmask for each node.
// this is implemented by takes the nodes and turning the node into a bitmask.
{
    vector<unsigned int> z;
    z.clear();
    
    z = this->nodes();
    
    // convert each node to a bitmask
    for( int i=0; i<z.size(); i++ ) {
        z[i] = 1 << z[i];
        
        // z[i] should never exceed the full_mask of the subset
        assert( z[i] <= this->S );
    }
        
    return z;
    
}


// returns whether bit 'i' of S is ON
bool t_subset::operator[] (unsigned int i) const { return this->at( i ); }
bool t_subset::in_subset( unsigned int i ) const { return this->at( i ); }

bool t_subset::at( unsigned int node_index ) const
// returns the value of S at node_index
{
//	cout << "node_index=" << node_index << "\tmax_mask=" << max_S << endl;	
	assert( ((1<<node_index)-1) <= max_S );
	return (bool) ((S >> node_index)&1);
}

string t_subset::str() const
// returns the nodes of S, sorted.
{
	string z = "";
	char buf[32];

	foreach( unsigned int node, this->nodes() ) {
		sprintf(buf, "%i ", node );		
		z.append( buf );
	}
	
	// now remove the final space
	
	// check that the final space is there
	assert( z[ z.size()-1 ] == ' ' );
	
	z.resize( z.size()-1 );
	
	// remove the final space
	assert( z[ z.size()-1 ] != ' ' );					
	
	return z;
}

// alais to binrep()
string t_subset::brep() const { return binrep(); }

string t_subset::binrep() const
// returns a string version of S
{
	string z = "";
	z.reserve( max_numnodes );
	
//	int length = this->size();
	
	for( int i=0; i<max_numnodes; i++ )
	{
		if( this->at(i) == true )
			z.insert(0, string("1"));
		else
			z.insert(0, string("0"));

	}
	
	return z;
}

// sets the node_index in S to ON
void t_subset::set_on( unsigned int node_index ) { S |= (1<<node_index); }

// sets the node_index in S to OFF
void t_subset::set_off( unsigned int node_index ) { S &= (~(1<<node_index)); }


unsigned int t_subset::numnodes() const
// returns the number fo bits that are ON
{
	
	return this->num_nodes;
/*	
	unsigned int z=0, tempS = this->S;
	
	do	{
		if( tempS&1 )
			z += 1;
		
		tempS >>= 1;
	}while( tempS );
	
	assert( z <= max_numnodes );
	
	return z;
*/	
}

// alias to numnodes()
unsigned int t_subset::size() const { return this->num_nodes; }

//return a vector of all supersets of the current S
vector<t_subset> t_subset::super_sets() const
{
	vector<unsigned int> outside_nodes = nodes(true);

	vector<t_subset> z;
	z.reserve( outside_nodes.size() );

	t_subset temp = *this;
	
	for( int i=0; i<outside_nodes.size(); i++ ) {
		assert( at(outside_nodes[i]) == false );

		t_subset temp( this->max_numnodes, S );
        temp.set_on( outside_nodes[i] );
		
		z.push_back( temp );
		
//		cout << "Pushing: " << temp.str() << endl;
	}
	
	assert( z.size() == outside_nodes.size() );
	
	return z;
}

vector<t_subset> t_subset::sub_sets() const
{
	vector<unsigned int> inside_nodes = nodes();
	
	vector<t_subset> z;
	z.reserve( inside_nodes.size() );
	
	t_subset temp = *this;
	
	for( int i=0; i<inside_nodes.size(); i++ ) {
		assert( at(inside_nodes[i]) == true );
		
		t_subset temp( this->max_numnodes, S );
        temp.set_off( inside_nodes[i] );
		
		z.push_back( temp );
		
//		cout << "Pushing: " << temp.str() << endl;
	}

	assert( z.size() == inside_nodes.size() );	
	
	return z;
}


//various equality operators with unsigned int
t_subset& t_subset::operator=(const unsigned int& inp_S) throw()
// set two t_partitions equal to another
{
//	cout << "running operator=" << endl;
	assert( inp_S <= max_S );
	
	this->S = inp_S;
	this->update_num_nodes();
	
	return *this;
}

t_subset& t_subset::operator=(const t_subset& restrict inp) throw()
// set two t_partitions equal to another
{
	// Check for self-assignment!
    if(this == &inp)		// Same object?
		return *this;       // Yes, so skip assignment, and just return *this.
	
	//	cout << "running operator=" << endl;
	this->S = inp.S;
	this->max_S = inp.max_S;
	this->max_numnodes = inp.max_numnodes;
	this->num_nodes = inp.num_nodes;
	
	assert( num_nodes <= max_numnodes );

	
	assert( 1 <= max_numnodes );
	assert( max_numnodes <= MAX_NUMUNITS );
	assert( S <= max_S );
	
	return *this;
}


//equality operators with int
t_subset& t_subset::operator=(const int& inp_S) throw()
// set two t_partitions equal to another
{
//	cout << "running operator=" << endl;
	
	assert( inp_S <= max_S );
	this->S = (unsigned int) inp_S;
	this->update_num_nodes();
	
	return *this;
}

bool t_subset::is_subset( const t_subset& restrict inp ) const
// returns TRUE if this subset is a IMPROPER NONSTRICT subset of inp.
{
    if( inp.S == this->S )
        return true;
    
    // if this->S is *bigger* (>=) than inp, return false.
	if( (inp.S | this->S) == this->S )
		return false;
	
	return true;
}

bool t_subset::is_superset( const t_subset& restrict inp ) const
// returns TRUE if this subset is an IMPROPER NONSTRICT superset of inp.
{
    if( inp.S == this->S )
        return true;
    
    // if this->S is *bigger* (>=) than inp, return true.
	if( (inp.S | this->S) == this->S )
		return true;

	return false;
}



// check if two subsets are equal.
bool t_subset::operator==(const t_subset &other) const { return (this->S == other.S); }
bool t_subset::operator!=(const t_subset &other) const { return !(*this == other ); }

bool t_subset::operator<(const t_subset &other) const { return (this->S < other.S); }
bool t_subset::operator<=(const t_subset &other) const { return (this->S <= other.S); }
bool t_subset::operator>=(const t_subset &other) const { return (this->S >= other.S); }
bool t_subset::operator>(const t_subset &other) const { return (this->S > other.S); }


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Stuff for comparing subsets and ints
////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool t_subset::operator==(const unsigned int &inp_int) const { return (this->S == inp_int ); }
bool t_subset::operator!=(const unsigned int &inp_int) const { return !(*this == inp_int ); }

bool t_subset::operator==(const int &inp_int) const { return (this->S == (unsigned int) inp_int ); }
bool t_subset::operator!=(const int &inp_int) const { return !(*this == inp_int ); }

bool t_subset::operator<(const unsigned int &inp_int) const { return (S < inp_int); }
bool t_subset::operator<(const int &inp_int) const { return ( S < (unsigned int) inp_int ); }
bool t_subset::operator<=(const unsigned int &inp_int) const { return (S <= inp_int); }
bool t_subset::operator<=(const int &inp_int) const { return ( S <= (unsigned int) inp_int ); }

bool t_subset::operator>(const unsigned int &inp_int) const { return (S > inp_int); }
bool t_subset::operator>(const int &inp_int) const { return ( S > (unsigned int) inp_int ); }
bool t_subset::operator>=(const unsigned int &inp_int) const { return (S >= inp_int); }
bool t_subset::operator>=(const int &inp_int) const { return ( S >= (unsigned int) inp_int ); }


////////////////////////////////////////////////////////////////////////////////////////////////////////////

//operators to tinker with ints
t_subset& t_subset::operator|=(const unsigned int &inp_int) throw() {	
	this->S |= inp_int;
	assert( this->S <= max_S );	
	this->update_num_nodes();
	
	return *this;
}

// returns x & S
unsigned int t_subset::operator&(const unsigned int &inp_int) throw() { return inp_int & S; }

t_subset& t_subset::operator+=(const unsigned int &inp_int) throw() { return operator|=( inp_int ); }

t_subset& t_subset::operator-=(const unsigned int &inp_int) throw() { S &= (~inp_int); return *this; }

t_subset& t_subset::operator++() throw() {
// this function makes ++t_subset work
	this->S++;
	this->update_num_nodes();

	assert( this->S <= max_S );
	return *this;
}

t_subset& t_subset::operator--() throw() {
	// this function makes ++t_subset work
	this->S--;
	this->update_num_nodes();	
	assert( this->S <= max_S );
	return *this;
}

	
// this function makes t_subset++ work
void t_subset::operator++(int) throw() { operator++(); }
void t_subset::operator--(int) throw() { operator--(); }
	
	

/*
int main()
{ 
	unsigned int test=19;
	t_subset S(6, test);
	
	cout << "subset=" << test << endl;
	cout << "string=" << S.str() << endl;
	cout << "---------------------------" << endl;
	
	cout << "Setting S = 5" << endl;
	S = 5;
	cout << "subset=" << test << endl;
	cout << "string=" << S.str() << endl;
	cout << "---------------------------" << endl;
	
	cout << "S += 2" << endl;
	S += 2;
	cout << "subset=" << test << endl;
	cout << "string=" << S.str() << endl;
	cout << "---------------------------" << endl;
	
	cout << "S++" << endl;
	S++;
	cout << "subset=" << test << endl;
	cout << "string=" << S.brep() << endl;
	cout << "---------------------------" << endl;
	cout << "T = 99" << endl;
	t_subset T = 99;
	
	cout << "T=" << T.brep() << endl;
	S = T = 100;
	cout << "S=" << S.brep() << endl;
	cout << "T=" << T.brep() << endl;
	cout << "---------------------------" << endl;
	
	t_subset W(8, 100 );
	cout << "W=\t\t" << W.brep() << " \t=\t" << W.str() << endl;		
	vector<t_subset> superset = W.super_sets();
	for( int i=0; i<superset.size(); i++ ) { 
		cout << "super[" << i << "]=\t" << superset[i].str() << endl;
	}

	
	cout << "W=\t\t" << W.brep() << " \t=\t" << W.str() << endl;		
	vector<t_subset> subset = W.sub_sets();
	for( int i=0; i<subset.size(); i++ ) { 
		cout << "sub[" << i << "]=\t\t" << subset[i].str() << endl;
	}
	
	cout << "================================" << endl;
	cout << "Adding node 3" << endl;
	W |= (1<<3);
	cout << "W=\t\t" << W.brep() << " \t=\t" << W.str() << endl;			
}

*/


#endif