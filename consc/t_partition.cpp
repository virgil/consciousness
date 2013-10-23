#ifndef T_PARTITION_CPP
#define T_PARTITION_CPP

#include <algorithm>
#include <assert.h>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "t_partition.h"
#include "helpers.h"
#include "t_consc.h"
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

//DIE( msg, __FILE__, __LINE__ )

using namespace std;


//int mask_compare(const void *a, const void *b) { return (*(int*)a - *(int*)b); }
//int mask_compare(const void *a, const void *b) { return (*(unsigned int*)b - *(unsigned int*)a); }
int mask_compare(const void *a, const void *b)
// this function compares two masks 'a' and 'b'. It returns which one should come first.
// this function is ONLY to be called by the qsort() function.
{
	const unsigned int s1=numunits_in_mask(*(unsigned int*)a), s2=numunits_in_mask(*(unsigned int*)b);

	if( s1 == s2 )
		return (*(unsigned int*)a) - (*(unsigned int*)b);

	return s2-s1;
}

void t_partition::sort()
// sorts the masks by the number of 1s in each mask
{
	
	qsort( masks, numparts, sizeof(unsigned int), mask_compare );
	
//	std::sort( this->masks, this->masks + numparts );
/*	
	assertmsg( 0 == 1, "this function corrupts things!!!!" );
	
	// foreach part...
	for( int i=0; i<(numparts-1); i++ )
	{
		// if the next mask is BIGGER than the previous, switch them.
		unsigned int s1=numunits_in_mask(masks[i]), s2=numunits_in_mask(masks[i+1]);
		
		// Put the MANY NODES part first.
		// But if two parts have EQUAL number of nodes, the one with the lower indices goes first
		if( s1 < s2 || (s1 == s2 && masks[i+1] < masks[i]) )		
		{
			unsigned int temp = masks[i];
			masks[i] = masks[i+1];
			masks[i+1] = temp;
			i = -1;
		}
		
	}
*/
}

void t_partition::clear()
// delete everything from the class
{
	// set the masks to 0
	memset( this->masks, 0, this->numparts * sizeof(unsigned int) );
	
	this->numparts = 0;
	this->subset = 0;

}

void t_partition::resize( unsigned int new_maxnumparts, unsigned int new_correct_subset )
// delete everything from the class
{
	// clear everything
	this->clear();
	delete [] this->masks;

	if( new_correct_subset > 0 ) {
		//		std::cout << "subset=" << subset << std::endl;
		//		std::cout.flush();
		new_maxnumparts = min(new_maxnumparts, numunits_in_mask(new_correct_subset) );
	}
	
	assert( new_maxnumparts > 0 );
	this->maxnumparts = new_maxnumparts;
	
	this->masks = new unsigned int[this->maxnumparts];

	
	// initialize the arrays to zero
	memset( this->masks, 0, this->maxnumparts * sizeof(unsigned int) );
	
	this->numparts = 0;
	this->subset = 0;
	this->correct_subset = new_correct_subset;
}

t_partition::t_partition( const vector<vector<int> >& restrict parts, unsigned int inp_correct_subset )
// the constructor using the vector-vector-ints
{

	assert( parts.size() >= 1 );

	if( inp_correct_subset ) {
		assert( parts.size() <= numunits_in_mask(inp_correct_subset) );
	}

	this->maxnumparts = parts.size();
	this->numparts = 0;
	this->subset = 0;
	this->correct_subset = inp_correct_subset;

	this->masks = new unsigned int[this->maxnumparts];
	this->update( parts );

}

t_partition::t_partition( const t_partition& restrict ref )
// make a copy
{
//	cout << "running copy constructor..." << endl;
//	cout.flush();

	this->maxnumparts = ref.size();
	this->numparts = ref.size();
	this->subset = ref.subset;
	this->correct_subset = ref.subset;
	
	this->masks = new unsigned int[maxnumparts];

	this->update( numparts, ref.masks );
/*	
	// make new array
//	unsigned int* temp = new unsigned int[numparts];
//	unsigned int temp[numparts];

//	cout << "numparts=" << numparts << endl;
//	cout.flush();

//	cout << "incoming partition: " << ref.MASKstr() << endl;
	
	// copy the inputmasks over the masks
//	memcpy( temp, ref.masks, sizeof(unsigned int) * numparts );

	for( int i=0; i<numparts; i++ )  {
		temp[i] = ref.at(i);
		assert( temp[i] == ref.at(i) );
		cout << "temp[" << i << "]=" << temp[i] << endl;
	}
*/
}

t_partition::t_partition( unsigned int inp_maxnumparts, unsigned int inp_correct_subset )
// the constructor using max_num_parts
{

	// if subset is defined, then the maximum number of parts
	// if the number of bits set to 1 in the subset
	if( inp_correct_subset > 0 ) {
		maxnumparts = min(inp_maxnumparts, numunits_in_mask(inp_correct_subset) );
	}
	
	assert( inp_maxnumparts > 0 );
	
	this->masks = new unsigned int[inp_maxnumparts];
	
	// initialize the arrays to zero
//	memset( this->masks, 0, inp_maxnumparts * sizeof(unsigned int) );

	this->maxnumparts = inp_maxnumparts;	
	this->numparts = 0;
	this->subset = 0;
	this->correct_subset = inp_correct_subset;
	
}
t_partition::t_partition( unsigned int inp_maxnumparts, const t_subset& restrict inp_correct_subset )
{
	// if subset is defined, then the maximum number of parts
	// if the number of bits set to 1 in the subset
	if( inp_correct_subset > 0 ) {
		maxnumparts = min(inp_maxnumparts, inp_correct_subset.numnodes() );
	}
	
	assert( inp_maxnumparts > 0 );
	
	this->masks = new unsigned int[inp_maxnumparts];
	
	// initialize the arrays to zero
	//	memset( this->masks, 0, inp_maxnumparts * sizeof(unsigned int) );
	
	this->maxnumparts = inp_maxnumparts;
	this->numparts = 0;
	this->subset = 0;
	this->correct_subset = inp_correct_subset.S;
}

// returns the number of parts
unsigned int t_partition::size() const { return this->numparts; }

// return the size of a given part
unsigned int t_partition::size(unsigned int part_index) const
{
	ASSERT( 0 <= part_index );
	ASSERT( part_index < this->numparts );
//	cout << "MASKs: " << this->MASKstr() << "\t part_index=" << part_index << endl;
	assert( this->masks[part_index] );
	
	return numunits_in_mask( this->masks[part_index] );

}

bool t_partition::empty() const
// returns TRUE if numparts is zero
{
	if( this->numparts == 0 )
		return true;
	else
		return false;
}


unsigned int t_partition::numnodes() const
// returns the number of nodes in the partition -- the number of bits in the subset
{
	return numunits_in_mask( this->subset );
}


bool t_partition::update( const vector<vector<int> >& restrict parts )
// If we get an update via a vector<vector<int> >s
// merely convert to an array pass the result to the real update()
{
	unsigned int numparts = parts.size();
	unsigned int temp_masks[numparts];

	for( int i=0; i<numparts; i++ ) {
		temp_masks[i] = this->part2mask( parts[i] );
	}

	return this->update( numparts, temp_masks );
}

bool t_partition::update( const vector<unsigned int>& restrict inputmasks )
// If we get an update via a vector<unsigned int>s
// merely convert to an array pass the result to the real update()
{
	unsigned int numparts = inputmasks.size();
	unsigned int temp_masks[numparts];

	for( int i=0; i<numparts; i++ ) {
		temp_masks[i] = inputmasks[i];
	}
	return this->update( numparts, temp_masks );
}

bool t_partition::update( unsigned int inp_numparts, const unsigned long* restrict inputmasks )
{
	
	assert( inp_numparts <= this->maxnumparts );
	
	unsigned int temp[inp_numparts];

	for( int i=0; i<inp_numparts; i++ ) {
		temp[i] = (unsigned int) inputmasks[i];
	}
	
	return this->update( inp_numparts, temp );
	
}

bool t_partition::update( unsigned int inp_numparts, const unsigned int* restrict inputmasks )
// updates numparts, the masks, the sizes, the subsets, the numunits
// returns TRUE if update was successful, otherwise returns false
{
	if( inp_numparts == 0 ) {
		cerr << "inp_numparts == " << inp_numparts << endl;
	}
	ASSERT( 1 <= inp_numparts );
	ASSERT( inp_numparts <= maxnumparts );
	
	// update numparts
	this->numparts = inp_numparts;
	this->subset = 0;

	// copy the inputmasks over the masks
	memcpy( this->masks, inputmasks, sizeof(unsigned int) * this->numparts );

	// assert that this->masks and inputmasks are identical
	assert( (memcmp(masks, inputmasks, sizeof(unsigned int) * numparts)) == 0 );
	
	// check the masks, create the subset
	for( int i=0; i<numparts; i++ )
	{
		// zero parts not allowed
		assert( masks[i] != 0 );

		// there should be NO overlap among parts
		assert( (subset & masks[i]) == 0 );

//		// update subset
		this->subset |= masks[i];
	}

	
//	for( int i=0; i<this->maxnumparts; i++ )
//		cout << "mask[" << i << "] = " << binary_repr(this->masks[i], maxnumparts+1) << endl;

	
//	cout << "[-1] MASKS: " << this->MASKstr() << endl;

	if( this->correct_subset == 0 )
		this->correct_subset = this->subset;
	
	assertmsg( this->subset == this->correct_subset, "the merged inputmasks did not match the correct_subset" );
	
//	cout << "input: " << binary_repr(subset,maxnumparts+1) << " \t\t correct:  " << binary_repr(correct_subset,maxnumparts+1) << endl;

	return true;
}

// return the first element in masks
unsigned int t_partition::front() const { return masks[0]; }

// return the last element in masks
unsigned int t_partition::back() const {
	ASSERT( this->masks[ numparts-1 ] != 0 );

	return masks[ this->numparts-1 ];
}

vector<int> t_partition::nodes( int part_index ) const
// returns mask[i] converted to a part
{
	assert( 0 <= part_index );
	assert( part_index <= this->numparts );
	
	// clear the old part
	vector<int> z;
	
	unsigned int mask = masks[part_index];
	unsigned int bit_index=0;
	
	do	{
		// this only works because we're starting at bit_index=0
		if( (mask&1) )
			z.push_back(bit_index);
		
		mask >>= 1;
		bit_index += 1;
	}while( mask );
	
	
	return z;
	
}

// returns a t_subset of the partition
t_subset t_partition::get_subset() const {
	t_subset z( numnodes(), subset );
	
	
	return z;
}


vector< vector<int> > t_partition::nodes() const
// returns the bit_masks as parts
{
	
	// make a vector of size numparts
	vector<vector<int> > z;
	z.reserve( this->numparts );

	for( int i=0; i<this->numparts; i++ ) {
		z.push_back( this->nodes(i) );
	}
	
	return z;
}

unsigned int t_partition::part2mask(const vector<int>& restrict part) const
// returns the mask of a part
{
	unsigned int z=0;
	
	for( int i=0; i<part.size(); i++ )
		z |= 1 << part[i];	// z += 2**part[i]

	return z;	
}

vector<unsigned int> t_partition::vec_masks() const
// returns the masks as a vector<unsigned int>
{
	// make a vector of size numparts
	vector<unsigned int> z;
	z.reserve( numparts );
	
	for( int i=0; i<numparts; i++ ) {
		z.push_back( masks[i] );
	}
	
	return z;
}

t_partition& t_partition::operator=(const t_partition& rhs) throw()
// set two t_partitions equal to another
{
//	cout << "running operator=" << endl;

	
    // Check for self-assignment!
    if(this == &rhs)		// Same object?
		return *this;       // Yes, so skip assignment, and just return *this.
	
	// if the incoming numparts is GREATER than our container
	// delete our container and make it bigger.
	if( rhs.numparts > maxnumparts ) {
		delete [] this->masks;
		this->masks = new unsigned int[rhs.numparts];
		this->maxnumparts = rhs.numparts;
	}
	
	this->correct_subset = rhs.correct_subset;
	
	// make sure there's enough space!
	assert( rhs.numparts <= maxnumparts );

	this->update( rhs.numparts, rhs.masks );
	
//	ASSERT( (memcmp(this->masks, rhs.masks, sizeof(unsigned int) * this->numparts)) == 0 );	
	return *this;
}


bool t_partition::operator==(const t_partition &rhs) const
// check if two partitions are equal
{
	if( correct_subset != rhs.correct_subset || subset != rhs.subset || numparts != rhs.numparts )
		return false;
	
	// check that every part is the same
	for( int i=0; i<numparts; i++ ) {
		if( masks[i] != rhs[i] )
			return false;
	}
	
	return true;
}

// is this partition LESS THAN the rhs?
bool t_partition::operator<(const t_partition &rhs) const
{
	// if this has MORE parts than rhs, this comes first
	if( numparts > rhs.numparts )
		return true;

	// if rhs has more parts, rhs comes first.
	else if( numparts < rhs.numparts )
		return false;
	
	// if they have the SAME number of parts, the part with the HIGHER product-of-partsizes comes first.
	return ( product_of_partsizes() > rhs.product_of_partsizes() );		
}



bool t_partition::operator!=(const t_partition &rhs) const { return !(*this == rhs ); }

unsigned int& t_partition::operator[] (int part_index) const
// returns mask[part_index]
{
	ASSERT( 0 <= part_index );
	ASSERT( part_index <= numparts );

	return this->masks[part_index];
}

unsigned int t_partition::at(int part_index) const
// returns mask[part_index]
{
	ASSERT( 0 <= part_index );
	ASSERT( part_index <= numparts );
	
	return this->masks[part_index];
}

// alias to BAREstr()
string t_partition::str() const { return this->BAREstr(); }

string t_partition::BAREstr() const
// return a BARE representation of the partition
{
	char temp[24];	
	string z="";

	for( int part_index=0; part_index<numparts; part_index++ )
	{
		vector<int> part_nodes = this->nodes(part_index);
		
		for( int j=0; j<part_nodes.size()-1; j++ )
		{
			snprintf(temp, sizeof(temp), "%d ", part_nodes[j]);
			z.append( temp );		// the element
		}
		
		snprintf(temp, sizeof(temp), "%d", part_nodes.back() );
		z.append( temp );			// the final element in a part
		
		
		if( part_index<numparts-1 )
			z.append("/");			// the part separator
	}

	return z;

	
}

unsigned int t_partition::product_of_partsizes() const
// returns the product of all of the numunits multipled together.
{
	unsigned int z=1.0;
	for( int i=0; i<numparts; i++ )
		z *= numunits_in_mask( masks[i] );

	assertmsg( 1 <= z, "Got a product-of-partsizes less than 1. That is impossible." );
	
	return z;	
}

string t_partition::MASKstr() const
// return a MASK representation of the partition
{
	assertmsg( 1<=numparts, "numparts was 0.  Cleared?" );
	string z="";
	
	for( int i=0; i<(this->numparts-1); i++ ) {
		stringstream ss;
		unsigned int test = this->masks[i];
		ss << test;
		z.append( ss.str() );
		z.append( string("/") );
		ss.flush();		
	}

	stringstream ss;
	ss << this->masks[ numparts-1 ];
	z.append( ss.str() );
	ss.flush();	
	
	return z;
	
}

t_partition::~t_partition( )
// deletes the arrays
{
	delete [] this->masks;
//	assert( this->masks == NULL );
}

unsigned int numunits_in_mask( unsigned int mask )
// returns the number of bits in the mask are 1
{
	unsigned int z=0;	
	
	do	{
		if( mask&1 )
			z += 1;
		
		mask >>= 1;
	}while( mask );
	
	assert( z != 0 );
	
	return z;
}

#endif