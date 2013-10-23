#pragma once
#ifndef T_PARTITION_H
#define T_PARTITION_H

#include <vector>
#include <assert.h>
#include "t_subset.h"
// WARNING: do NOT include helpers.h

/////////////////////////////////////////////////////////////////////////////////////
// CONFIGURABLE PARAMETERS
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


#define restrict __restrict__
using namespace std;

unsigned int numunits_in_mask( const unsigned int mask );

class t_partition {	
public:

	//////////////////////////////////////////////////////
	// Variables
	//////////////////////////////////////////////////////	
	
	// sum of the sizes -- equivalent to the number of units in the subset
	unsigned int subset;
	unsigned int correct_subset;

	unsigned int* masks;	
	
	
	//////////////////////////////////////////////////////
	// Functions
	//////////////////////////////////////////////////////	

	// 'correct subset' is a sanity check that is checked after an update()
	t_partition( unsigned int maxnumparts, unsigned int correct_subset=0 );
	t_partition( unsigned int maxnumparts, const t_subset& restrict correct_subset );	
	t_partition( const vector<vector<int> >& restrict parts, unsigned int correct_subset=0 );	
	t_partition( const t_partition& restrict ref );
	
//	t_partition( unsigned int numparts, unsigned int* inputmasks );
	~t_partition();
	
	
	// returns a t_subset of the partition
	t_subset get_subset() const;
	
	// updates numparts, the masks, the sizes, the subsets, the numunits
	bool update( unsigned int numparts, const unsigned int* restrict inputmasks );
	bool update( unsigned int numparts, const unsigned long* restrict inputmasks );
	
	bool update( const vector<unsigned int>& restrict inputmasks );	
	bool update( const vector<vector<int> >& restrict parts );
	bool update( const t_partition& restrict P );

	// returns the bitmasks as a vector<int>
	vector<unsigned int> vec_masks() const;
	
	// returns the the nodes of the entire partition
	vector< vector<int> > nodes() const;

	// returns the nodes of a given part
	vector<int> nodes( int part_index ) const;
	
	void clear();			// delete everything	
	bool empty() const;		// TRUE if has anything in it
	
	// the inline functions
	////////////////////////////////////////////
	unsigned int front() const;
	unsigned int back() const;
	unsigned int& operator[] (int part_index) const;		// returns mask[part_index]


	unsigned int  at(int part_index) const;			// returns mask[part_index]		

	unsigned int size(void) const;							// returns the number of parts
	unsigned int size(unsigned int part_index) const;		// returns size of a part

	// BARE string representation
	string BAREstr() const;
	string str() const;	
	string MASKstr() const;


	// returns the number of nodes in the partition -- the number of bits in the subset
	unsigned int numnodes() const;	
	t_partition& operator=(const t_partition& rhs) throw();

	bool operator==(const t_partition &rhs) const;
	bool operator!=(const t_partition &rhs) const;
	bool operator<(const t_partition &rhs) const;
		
	void resize( unsigned int newsize, unsigned int new_correct_subset=0 );

	// sorts the masks by the number of units
	void sort();

	// returns the size of all parts multipled together
	unsigned int product_of_partsizes() const;
	
private:
//	unsigned int correct_subset;
	unsigned int numparts;
	unsigned int maxnumparts;
		
	unsigned int part2mask( const vector<int>& restrict ) const;

	
};

#endif


