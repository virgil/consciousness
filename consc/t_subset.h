#pragma once
#ifndef T_SUBSET_H
#define T_SUBSET_H

#include <vector>
#include <list>

typedef unsigned int bitmask;

#define restrict __restrict__
/////////////////////////////////////////////////////////////////////////////////////
// CONFIGURABLE PARAMETERS
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

using namespace std;

class t_subset {	
public:

	// Constructors
	t_subset( unsigned int numunits, unsigned int subset=0 );
	t_subset();

	
	
	//////////////////////////////////////////////////////
	// Variables
	//////////////////////////////////////////////////////	
	
	// sum of the sizes -- equivalent to the number of units in the subset that 

	// returns the the nodes of the entire partition
	// the inline functions
	////////////////////////////////////////////
	bool operator[] (unsigned int node_index) const;		// returns mask[part_index]
	bool at(unsigned int node_index) const;		
	bool in_subset( unsigned int node_index ) const;

	// sets a node in S to ON or OFF
	void set_on( unsigned int node_index );
	void set_off( unsigned int node_index );
	
	bool is_subset( const t_subset& restrict inp ) const;
	bool is_superset( const t_subset& restrict inp ) const;

	/////////////////////////////////////////////////////////////////////////////////
	
	unsigned int size() const;
	unsigned int numnodes() const;

	////////////////////////////////////////////////////////////////////
	// various equality operators to work with ints
	////////////////////////////////////////////////////////////////////
	t_subset& operator=(const unsigned int& inp_S) throw();
	t_subset& operator=(const int& inp_S) throw();
	t_subset& operator=(const t_subset& S) throw();
	t_subset& operator+=(const unsigned int &inp_int) throw();
	t_subset& operator-=(const unsigned int &inp_int) throw();
	t_subset& operator|=(const unsigned int &inp_int) throw();


	unsigned int operator&(const unsigned int &inp_int) throw();
	
	
	t_subset& operator+(const int &inp_int) throw();
	t_subset& operator-(const int &inp_int) throw();	
	
	bool operator==(const unsigned int &inp_int) const;
	bool operator==(const int &inp_int) const;
	bool operator!=(const unsigned int &inp_int) const;
	bool operator!=(const int &inp_int) const;

	bool operator<(const unsigned int &inp_int) const;
	bool operator<(const int &inp_int) const;
	bool operator<=(const unsigned int &inp_int) const;
	bool operator<=(const int &inp_int) const;	
	bool operator>(const unsigned int &inp_int) const;
	bool operator>(const int &inp_int) const;
	bool operator>=(const unsigned int &inp_int) const;
	bool operator>=(const int &inp_int) const;
	
	////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////
	// various equality operators to work with itself
	////////////////////////////////////////////////////////////////////	
	bool operator==(const t_subset &other) const;
	bool operator!=(const t_subset &other) const;

	t_subset& operator++() throw();		// ++t_subset
	void operator++(int) throw();		// t_subset++

	t_subset& operator--() throw();		// --t_subset
	void operator--(int) throw();		// t_subset--

	bool operator<(const t_subset &other) const;
	bool operator<=(const t_subset &other) const;
	bool operator>=(const t_subset &other) const;
	bool operator>(const t_subset &other) const;	
	////////////////////////////////////////////////////////////////////	
	
	
	void update_num_nodes();

	vector<unsigned int> masks() const;	
	vector<unsigned int> nodes( bool outside=false ) const;
	void array( unsigned int* array, bool outside=false ) const;

	vector<t_subset> super_sets() const;
	vector<t_subset> sub_sets() const;
	
	
	// return the subset as a string
	string str() const;
	string binrep() const;	
	string brep() const;		// alias to binrep()

	unsigned int S;
	unsigned int max_numnodes;
	unsigned int max_S;
	unsigned int num_nodes;
	
		
private:
	// numunits is the highest allowed value of S.
//	unsigned int max_numnodes;
//	unsigned int max_S;
};



#endif