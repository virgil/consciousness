/////////////////////////////////////////////////////////////////////////////////////
// This defines a class 'state' that represents the state of the network
/////////////////////////////////////////////////////////////////////////////////////

// #pragma once
#ifndef T_STATE_H
#define T_STATE_H

#include <string>
#include <boost/cstdint.hpp>
#include "helpers.h"

//typedef unsigned int bitmask;
//#include <vector>
//#include <assert.h>
//#include "t_subset.h"

/////////////////////////////////////////////////////////////////////////////////////
// CONFIGURABLE PARAMETERS
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


#define restrict __restrict__
using namespace std;

class t_state
{

    public:
        bitmask mask;
        unsigned int value;

        t_state( const bitmask inp_mask, const unsigned int inp_val );
        // ~t_state();


        void update( const bitmask inp_mask, const unsigned int inp_value );

        // return the string representation of this state
        string lstr() const;
        string str(const unsigned int FULL_MASK) const;
        string MASKstr( const bitmask FULL_MASK=0 ) const;

        void clear();			// delete everything
        bool empty() const;		// TRUE if has anything in it

    //	bool operator[] (const unsigned int node_index) const;		// returns mask[part_index]
    //	bool at(const unsigned int node_index) const;

        // return the number of bits ON in this_mask
        unsigned int size() const;

        const bool& operator[] (const unsigned int bit_index) const;	// returns the bool of bit_index. 0=RIGHT SIDE.
        const bool at(const unsigned int bit_index) const;			// returns the bool of bit_index. 0=RIGHT SIDE.


    private:
        // return the number of bits that are ON in inp
        unsigned int numbits_on( const bitmask inp ) const;

        // returns an array of the indices for which the bits of inp are ON
        // returns the number of bits that were ON
        void bits_on( bitmask inp, unsigned int* restrict z, const unsigned int z_length ) const;


};


#endif
