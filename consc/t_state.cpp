/////////////////////////////////////////////////////////////////////////////////////
// This defines the functions for class 'state' that represents the state of the network
/////////////////////////////////////////////////////////////////////////////////////

#ifndef T_STATE_CPP
#define T_STATE_CPP

#include <vector>
#include <assert.h>
#include <string>
#include "t_state.h"


// defines a base37 representation for a digit
const char base37[] = { '0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z' };

void t_state::clear()
// delete everything from the class
{
    this->mask = 0;
    this->value = 0;
    
    assert( mask == 0 );
    assert( value == 0 );
}


bool t_state::empty() const
// is this class empty?
{
    // assert that value is a subset of mask
    assert( (mask|value) == mask);
    
    if( mask == 0 )
        return true;
    
    return false;
}

const bool& t_state::operator[] (const unsigned int bit_index) const
// returns the bool of bit_index. 0=RIGHT SIDE.
{
    // we have to do this because the operator[] always returns a reference.
    static bool z;
    
    // the bit_index should never be less than the size
//    assert( bit_index < this->size() );

    z = this->at(bit_index);
    
    return z;
}

const bool t_state::at(const unsigned int bit_index) const
// returns the bool of bit_index. 0=RIGHT SIDE.
{
    // assert bit_index is <= 32
    assert( bit_index <= 32 );
    
    // assert bit_index is a subset of mask
    assert( ( (1 << bit_index)|mask) == mask );
    
    // get the i'th bit of value
    return (const bool) ((value >> bit_index)&1);
}

t_state::t_state( const bitmask inp_mask, const unsigned int inp_value )
// create a new state from inp_mask and inp_value
{
    this->update( inp_mask, inp_value );
}


void t_state::update( const bitmask inp_mask, const unsigned int inp_value )
// create a new state from inp_mask and inp_value
{
    // assert that the inp_mask exists.
    assert( inp_mask );    
    
    
    // assert that inp_value is a subset of inp_mask
    assert( (inp_value|inp_mask) == inp_mask );
    
    
    this->mask = inp_mask;
    this->value = inp_value;
}


string t_state::lstr() const
// returns a LONG string (lstr) representation of value and mask
// will return MMMMM:VVVVVV
// where M = mask_index from this->mask and V = {0,1} from this->value.
{
    const unsigned int MASKsize = numbits_on( this->mask );
    unsigned int* restrict MASKindices = new unsigned int[MASKsize];

    // fill the MASKindices
    bits_on( this->mask, MASKindices, MASKsize );
    

	stringstream z (stringstream::in | stringstream::out);
	
    for( int i=(MASKsize-1); i>=0; i-- ) {
        assert( MASKindices[i] <= 36 );
        z << base37[ MASKindices[i] ];
    }
    
    z << ":";
    
//    z.append( ss.str() );
    
    
    //Now for each node of this->mask determine the node's value
    for( int i=(MASKsize-1); i>=0; i-- )
    {
        // assert this MASKindex is infact part of the mask
        const bitmask node_mask = (1 << MASKindices[i]);
        const unsigned int node_value = (this->value & node_mask);
        const unsigned int node_value_shifted = node_value >> MASKindices[i];

        // assert that node_mask is a subset of mask
        assert( (node_mask|mask) == mask );        
        assert( node_value_shifted == 0 || node_value_shifted == 1 );
        
        z << node_value_shifted;
    }
    
//    z.append( ss.str() );
    

    delete [] MASKindices;
    
    // the length of string should be 2*MASKsize+1
    assert( z.str().length() == (2*MASKsize)+1);
    
//    cout << "z='" << z.str() << "'" << endl;
    
	return z.str();
}


string t_state::str(const unsigned int FULL_MASK) const
// returns the string representation with the value and mask together.
// will return _v_vv__v
// where '_' is an entry where FULL_MASK is 1 but this->mask is zero.
// and 'v' is a {0,1} from this->value where this->mask is one.
{
    // assert FULL_MASK is bigger than this->mask
    assert( (this->mask|FULL_MASK)==FULL_MASK );
    
    const unsigned int FMsize = numbits_on( FULL_MASK );
    
    // make a character array of size FMsize
    vector<char> chars;
    
    // set each entry of chars to '_'
    for( int i=0; i<FMsize; i++ )
        chars.push_back( '_' );

    
    assert( chars.size() == FMsize );
    
    // now determine the MASKindices, and set those entries of chars
    
    const unsigned int MASKsize = numbits_on( this->mask );
    unsigned int* restrict MASKindices = new unsigned int[MASKsize];
    
    // fill the MASKindices
    bits_on( this->mask, MASKindices, MASKsize );
    
    
    for( int i=0; i<MASKsize; i++ ) {
        unsigned int node_index = MASKindices[i];
        bitmask node_mask = (1 << MASKindices[i]);
        unsigned int node_value = value & node_mask;
        unsigned int node_value_shifted = node_value >> node_index;
        
        // node_value_shifted should be 0 or 1
        assert( node_value_shifted == 0 || node_value_shifted==1 );
        
//        cout << "\t setting node_index=" << node_index << " to val=" << node_value_shifted << endl;
        
        assert( node_index < FMsize );
        if( node_value_shifted == 0 )
            chars[ node_index ] = '0';
        else
            chars[ node_index ] = '1';
//        chars[ node_index ] = node_value_shifted;
    }
    
    delete [] MASKindices;
    
	stringstream z (stringstream::in | stringstream::out);

    
    // now make a string
    // remember to ***REVERSE*** the chars
    for( int i=chars.size()-1; i>=0; i-- ) {
        z << chars[i];
    }

//    assert( chars.size() == z.length() );

	return z.str();
}


unsigned int t_state::size() const
// returns the number of bits ON in mask
// We do the static oldsize,oldmask for efficiency purposes.
{
    static bitmask oldmask=0;
    static unsigned int oldsize=0;
    
    // if we've seen this oldmask before, just return the oldsize
    if( mask == oldmask )
        return oldsize;
    
    // mask didn't match the oldmask, so set it and return the 
    oldmask = mask;
    oldsize = numbits_on( mask );
    
    return oldsize;
}


//////////////////////////////////////////////////////////
// PRIVATE FUNCTIONS
//////////////////////////////////////////////////////////
unsigned int t_state::numbits_on( bitmask inp ) const
// Returns the number of 1s in the inp
{

	unsigned int z=0;	
    
	do	{
        
		if( inp&1 )
			z += 1;
        
		inp >>= 1;
        
	}while( inp );
	
    
	return z;
}



void t_state::bits_on( const bitmask input, unsigned int* restrict z, const unsigned int z_length ) const
// returns the indices of the bits are the ON in inp
{
	unsigned int num_ON=0, power=0, inp=input;
	
	do	{
		// this only works because we're starting at power=0 and num_ON=0
		if( inp & 1 ) {
			z[num_ON] = power;
			num_ON += 1;
		}
		
		inp >>= 1;
		power += 1;
	}while( inp );
    
    assert( num_ON == z_length );
    
//    return num_ON;
}


#endif
