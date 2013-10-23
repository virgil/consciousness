#ifndef T_LATTICEATOM_CPP
#define T_LATTICEATOM_CPP

#include <vector>
#include <list>
#include "t_latticeatom.h"

using namespace std;

t_LatticeAtom::t_LatticeAtom( const unsigned int bitmask, const bool timestep )
// if a bitmask and a timestep are passed, set them.
{
	
	if( timestep ) {
		this->X0mask = 0;		
		this->X1mask = bitmask;
	}
	else {
		this->X0mask = bitmask;
		this->X1mask = 0;
	}	
}


t_LatticeAtom::t_LatticeAtom()
// definition for if nothing is passed.
{
	this->X0mask = 0;
	this->X1mask = 0;
}


#endif