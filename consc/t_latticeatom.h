#pragma once
#ifndef T_LATTICEATOM_H
#define T_LATTICEATOM_H

#include <vector>
#include <list>


#define restrict __restrict__
/////////////////////////////////////////////////////////////////////////////////////
// CONFIGURABLE PARAMETERS
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

using namespace std;

// this defines a single entry/slice (ex: {12} or {123})  of a lattice node
class t_LatticeAtom {

	public:
		unsigned int X0mask;
		unsigned int X1mask;

	t_LatticeAtom( const unsigned int bitmask, const bool timestep );
	t_LatticeAtom();
};


typedef t_LatticeAtom LatticeAtom;
typedef vector<t_LatticeAtom> LatticeAtoms;
typedef vector<t_LatticeAtom> LatticeNode;
typedef vector<LatticeNode> LatticeNodes;


#endif