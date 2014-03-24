///////////////////////////////////////////////////////////////////////////////////////////
// These are functions that are part of the class t_consciousness
// but they are rarely used anymore.  These rarely-used functions are moved
// to more easily work with the more recent code.
// by Virgil Griffith -- virgil@caltech.edu
//////////////////////////////////////////////////////////////////////////////////////////
#ifndef T_CONSC_PART_TWO_CPP
#define T_CONSC_PART_TWO_CPP
//////////////////////////////////////////////////////////////////////////////////////////
#include "t_consc.h"
//#include "helpers.h"
//#include "helpers.cpp"
#include "t_consc.cpp"
//////////////////////////////////////////////////////////////////////////////////////////

//#ifndef ASSERT
//	#define ASSERT(exp) //defined as nothing
//#endif

extern const unsigned int MAXIMUM_NUMBER_OF_MIPS_TO_STORE;
extern const unsigned int MAXIMUM_NUMBER_OF_MAIN_COMPLEXES_TO_STORE;

extern void subset2array( unsigned int subset, unsigned int* restrict z );
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

inline unsigned int t_consciousness::env( const bitmask part, bitmask subset )
// returns the ENVIRONMENT of all nodes in part
{
    if( subset == 0 )
        subset = FULL_MASK;
    
    // assert that FULL_MASK fully encompasses the subset
    assert( (subset && FULL_MASK) | (subset == FULL_MASK) );

    // assert that subset fully encompasses the part.
    assert( (part && subset) | (part == subset) );
               
    // get all of the bits in subset are NOT in the part
	const bitmask z = subset & (~part);
    const bitmask z2 = subset ^ part;
    
    assert( z == z2 );
    
//    cout << "FULL_MASK \t=" << binary_repr( FULL_MASK, numunits ) << endl;
//    cout << "subset \t=" << binary_repr( subset, numunits ) << endl;
//    cout << "part \t=" << binary_repr( part, numunits ) << endl;    
//    cout << "z    \t=" << binary_repr( z, numunits ) << endl;
//    cout.flush();
    
    // assert that subset fully encompasses z
    assert( subset| (z == subset) );

	return z;
}




bool t_consciousness::can_calc_by_perturbing_elements( const bitmask part_mask )
// returns FALSE if at least one outside node projects to multiple (>=2) nodes in the part
// else returns TRUE
// If returns TRUE then entropy_of_part__WIRES() == entropy_of_part__ELEMENTS()
{
	// 0. Calculate all of the part_nodes
	unsigned int numnodes = numunits_in_subset( part_mask );
	unsigned int* nodes = new unsigned int[numnodes];
	subset2array(part_mask, nodes);
	//	vector<int> nodes = mask2nodes( part_mask );
	//	unsigned int numnodes = nodes.size();
	
	
	cout << "part=" << binary_repr(part_mask) << endl;
	
	for( int i=0; i<numnodes; i++ ) {
		vector<int> wires = external_nodes_with_wires_into_node( nodes[i], nodes, numnodes );
		
		cout << "\t node: " << nodes[i] << " <- "; 
		for( int j=0; j<wires.size(); j++) {
			cout << wires[j] << " ";
		}
		cout << "\t # incoming wires=" << wires.size() << endl;
	}
	
	
	// Make a vector indicating the nodes that have outside connections to this part
	//	vector<bool> outside_cxns( numunits );
	unsigned int outside_cxns = 0;
	
	if( networktype == "neural" )
	{
		for( int from_node=0; from_node<numunits; from_node++ )
		{
			// skip all from_nodes within the part
			if( (part_mask>>from_node)&1 )
				continue;
			
			for( int i=0; i<numnodes; i++ ) {
				int to_node = nodes[i];
				
				// this ensures that outside_cxns never contains a node
				// within the part
				
				if( weightmatrix[from_node][to_node] != 0.0 && ((part_mask >> to_node)&1) )
				{
					
					
					//					if( outside_cxns[from_node] == true ) {
					if( (outside_cxns >> from_node) & 1 ) {
						cout << "from_node=" << from_node << "\t bitset=" << binary_repr(outside_cxns) << endl;
						cout << "WILL NOT BE SAME!" << endl;
						return false;
					}
					else {
						//						outside_cxns[from_node] = true;
						outside_cxns |= (1 << from_node);
						cout << "from_node=" << from_node << "\t bitset=" << binary_repr(outside_cxns) << endl;						
					}
				}
			}
			
		}
	}
	
	cout << "WILL BE SAME! -- TRUE" << endl;	
	return true;
	
	/*
	 // If this is a circuit, use A and B
	 else {
	 //		cerr << endl << "getting transformations[" << node_index << "]." << endl;
	 if( this->transformations[node_index].A >= 0 )
	 z.push_back( this->transformations[node_index].A );
	 if( this->transformations[node_index].A >= 0 && this->transformations[node_index].B != this->transformations[node_index].A )
	 z.push_back( this->transformations[node_index].B );
	 }
	 
	 
	 */
	
}



vector<int> t_consciousness::external_nodes_with_wires_into_node( unsigned int node_index, const unsigned int* restrict part, unsigned int numnodes, unsigned int subset )
// returns the the nodes in the network that have connections to node_index yet ARE NOT within the part
{
	
	// OPTIMIZATION -- REMOVE ANY NODES THAT ARE NOT IN THE SUBSET
	
	// 1. Get all nodes that connect to node node_index.
	vector<int> z;
	z.clear();
	
	if( this->networktype == "neural" )
	{
		for( int i=0; i<this->numunits; i++ )
		{
			if( this->weightmatrix[i][node_index] != 0.0 )
				z.push_back( i );
		}
	}
	
	// If this is a circuit, use A and B
	else {
		//		cerr << endl << "getting transformations[" << node_index << "]." << endl;
		if( this->transformations[node_index].A >= 0 )
			z.push_back( this->transformations[node_index].A );
		if( this->transformations[node_index].A >= 0 && this->transformations[node_index].B != this->transformations[node_index].A )
			z.push_back( this->transformations[node_index].B );
	}
	
	ASSERT( this->networktype == "circuit" || this->networktype == "neural" );
	
	// sort the vector z
	sort( z.begin(), z.end() );
	
	//remove all elements of z that are in part
	for( int i=0; i<z.size(); i++ )
	{
		for(int j=0; j<numnodes;j++ )
		{
			if( z[i] == part[j] ) {
				z.erase( z.begin() + i );
				i=-1;
				break;
			}
		}
	}
	
	return z;
}


vector<int> t_consciousness::states_to_attempt__WIRES( unsigned int node_index, int mu0, const unsigned int* restrict part_nodes, unsigned int numnodes )
// returns a vector of the x0 states to run for node_index in part, with partstate mu0
{
	vector<int> z;
	z.clear();
	
	// Now to rotate the external nodes that this node is connected to.
	vector<int> wire_nodes = this->external_nodes_with_wires_into_node( node_index, part_nodes, numnodes );
	
	// if there were no wire_nodes, we are now done.
	if( wire_nodes.empty() )
		return z;
	
	int num_wirestates = (1 << wire_nodes.size() );
	z.reserve(num_wirestates);
	
	// make the partstate
	int proto_x0=0;	
	for( int i=0; i<numnodes; i++ )
	{
		// take the i'th bit of mu0 and shift it by the node_index[i]
		proto_x0 |= ( ((mu0 >> i)&1) << part_nodes[i] );			//partstate += ( ((mu0 >> i)&1) <<part_nodes[i]);
	}
	
	
	for( int wirestate=0; wirestate<num_wirestates; wirestate++ )
	{
		
		int x0state = proto_x0;
		
		// add the bits of the wire to the x0state, and append it.
		for( int i=0; i<wire_nodes.size(); i++ )
		{
			//take the i'th bit of wirestate and shift it by the wire_ondes[i]
			x0state |= ( ((wirestate >> i)&1) << wire_nodes[i]);	//x0state += ( ((wirestate >> i)&1) <<wire_nodes[i]);
		}
		
		z.push_back( x0state );
	}
	
	if( z.empty() )
		assert( num_wirestates == 0 );
	
	return z;
}



#endif