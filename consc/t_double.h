#ifndef T_SUBSET_H
#define T_SUBSET_H

#include <vector>
#include <list>
#include <assert.h>
/////////////////////////////////////////////////////////////////////////////////////
// CONFIGURABLE PARAMETERS
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

using namespace std;

// custom class t_double that is a double that overrides ==, >=, <=, !=
class t_double {	
public:

	bool operator==(const unsigned int &inp_int) const;
	bool operator==(const int &inp_int) const;
	
	bool operator!=(const unsigned int &inp_int) const;
	bool operator!=(const int &inp_int) const;

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
	
	vector<unsigned int> nodes( bool outside=false ) const;

	vector<t_subset> super_sets() const;
	vector<t_subset> sub_sets() const;
	
	
	// return the subset as a string
	string str() const;
	string binrep() const;	
	string brep() const;		// alias to binstr()

	unsigned int S;
	
		
private:
	// numunits is the highest allowed value of S.
	unsigned int max_numnodes;
	unsigned int max_S;
};



#endif