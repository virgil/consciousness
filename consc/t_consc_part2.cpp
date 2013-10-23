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

extern const bool FLAG__CACHE_M1_GIVEN_ENV;
extern const bool FLAG__CORRELATION_CASCADE_FORCE_POSITIVE_WHOLE;
extern const bool FLAG__MMI_IS_INTERACTION_INFORMATION;
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



double t_consciousness::H_M1_GIVEN_ENV( unsigned int M1 )
// Calculates H[M_1|env(M1)] = H[M_1|e0,e1]
{
	//	if( M1 == 0 )
	//		return 0.0;
	
	if( M1 == FULL_MASK )
		return H_X1;
	
	
	assert( 0 < M1 );
	assert( M1 < FULL_MASK );
	
	unsigned int ENV = env( M1 );
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////
	// Do we have this entropy in the cache?  If so, return that.
	////////////////////////////////////////////////////////////////////////////////////////////////
	// if the cache_part_entropies is not NONE, then it's a valid cache
	if( FLAG__CACHE_M1_GIVEN_ENV && H_M1_GIVEN_ENV_cache[M1] >= 0.0 ) {
		//		cerr << "returning cached entropy of H[M1=" << binrep(M1) << "|env(M1)=" << H_M1_GIVEN_ENV_cache[M1] << endl;
		return H_M1_GIVEN_ENV_cache[M1];
	}
	////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	const unsigned int M1size=numunits_in_mask(M1);
	const unsigned int ENVsize=numunits - M1size;
	
	const unsigned int M1d = 1 << M1size;		//d=2**part.size()
	const unsigned int ENVd = 1 << ENVsize ;	//d=2**part.size()	
	
	unsigned int* restrict env0env1mu1_instances = (unsigned int*) calloc(ENVd*ENVd*M1d, sizeof(unsigned int));
	assert( env0env1mu1_instances != NULL );
	
	
	unsigned int* restrict M1nodes = new unsigned int[M1size];
	unsigned int* restrict ENVnodes = new unsigned int[ENVsize];
	
	mask2array( M1, M1nodes );
	mask2array( ENV, ENVnodes );
	
	//////////////////////////////////////////////////////////////////////////////////////
	// ATTEMPTS HERE TO USE SPARSE MATRICES usually resulted in slower code             //
	// i.e: map<int,double> prob_mu1s; vector< map<int,double> > prob_mu0_given_mu1(d); //
	//////////////////////////////////////////////////////////////////////////////////////	
	//foreach x0...
	//	unsigned int mu0=0, mu1=0;
	//	unsigned int old_X0OUT_mask=UINT_MAX, old_X1OUT_mask=UINT_MAX;
	unsigned int old_S0_mask = UINT_MAX, old_M1_mask=UINT_MAX;
	
	for(unsigned int x0=0; x0<numstates; x0++ )
	{
		const unsigned int x1 = states[x0];		
		unsigned int mu1=0, env0=0, env1=0;
		
		// define: env0, env1, mu1
		//		Xstate2MUstate( ENVnodes, ENVsize, x0, env0 );
		//		Xstate2MUstate( ENVnodes, ENVsize, x1, env1 );
		Xstate2MUstate( ENVnodes, ENVsize, x0, env0, x1, env1 );
		Xstate2MUstate( M1nodes, M1size, x1, mu1 );
		
		assert( mu1 < M1d );
		assert( env0 < ENVd );
		assert( env1 < ENVd );
		
		
		env0env1mu1_instances[ env0*(ENVd*M1d) + env1*(M1d) + mu1 ] += 1;
	}
	
	// don't need the nodes anymore.
	delete [] ENVnodes;
	delete [] M1nodes;
	
	
	const double Dnumstates = (double) this->numstates;
	//	double H_M1_ENV=0.0;
	double z=0.0;
	double sum_prob=0.0;
	
	//	cout << "total=" << Dnumstates << " numerators: ";
	// create H( M1,S0 )
	for( int env0=0; env0<ENVd; env0++ )
	{
		for( int env1=0; env1<ENVd; env1++ )
		{
			const unsigned int start = env0*(ENVd*M1d) + env1*(M1d);
			unsigned int env0env1_instances = 0;
			for( int i=0; i<M1d; i++ ) 
				env0env1_instances += env0env1mu1_instances[ start + i ];
			
			const double prob_env0env1 = env0env1_instances / Dnumstates;
			
			for(int mu1=0; mu1<M1d; mu1++)
			{
				//			cout << mu0mu1_instances[(mu0*M1d)+mu1] << " ";
				
				// TODO: Go by 1 instead of recalculating the new position every time
				if( env0env1mu1_instances[ env0*(ENVd*M1d) + env1*(M1d) + mu1 ] == 0 )
					continue;
				
				const double prob_env0env1mu1 = env0env1mu1_instances[ env0*(ENVd*M1d) + env1*(M1d) + mu1 ] / Dnumstates;
				assert( prob_env0env1mu1 <= prob_env0env1 );
				
				// H( M1 | ENV ) = p(mu1,env0,env1) * -log2( p(mu1|env0,env1) )
				z -= prob_env0env1mu1 * log2( prob_env0env1mu1 / prob_env0env1 );
				sum_prob += prob_env0env1mu1;
			}
		}
	}
	
	
	// free the memory
	free( env0env1mu1_instances );
	
	//	cerr << "sum_prob=" << sum_prob << endl;
	assert( fequals(sum_prob, 1.0) );
	
	
	if( z < 0.0 ) {
		if( myceil(z) < 0.0 ) {
			cerr << "function: entropy_of_part__ELEMENTS" << endl;
			cerr << "Warning: Had negative entropy in entropy_of_part() that didn't round to 0.0 (" << myceil(z) << ").  Maybe decimal precision is too high?" << endl;
			assert( 0.0 <= myceil(z) );
		}
		z = 0.0;
	}
	
	
	// remove -0.0s
	z = fabs(z);
	
	
	if( FLAG__CACHE_M1_GIVEN_ENV ) {
		// Store the part entropy in the cache, if we're doing that.
		ASSERT( this->H_M1_GIVEN_ENV_cache[M1] != this->H_M1_GIVEN_ENV_cache[M1] );
		
		// set the part entropy
		this->H_M1_GIVEN_ENV_cache[ M1 ] = z;
	}
	
	
	
	// DO NOT DO ANY ROUNDING HERE!	
	ASSERT( z >= 0.0 );
	
	//	cerr << "- H[M1|env]=" << z << endl;
	//	cerr.flush();
	
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


t_ei_result t_consciousness::highest_bracket_ei()
// this functions returns the subsets with the highest bracket_ei
{
	double highest_ei=-1.0;
	t_ei_result z;
	
	//we start from the FULL_MASK because we can often skip steps that way
	for( unsigned int subset=this->FULL_MASK; subset>=1; subset-- )
	{
		//		cout << "Starting subset=" << subset << endl;
		//OPTIMIZATION: only attempt all partitions to find the MIPs if it's possible
		//for this subset to have higher or equal PHI to the current highestphi
		double this_ei = this->bracket_ei(subset);
		
		if( this_ei > highest_ei ) {
			z.subsets.clear();
			z.ei = this_ei;
			z.subsets.push_back(subset);
		}
		else if( this_ei == highest_ei ) {
			z.subsets.push_back(subset);
			
			// break if we have the maximum number of main complexes
			if( z.subsets.size() >= MAXIMUM_NUMBER_OF_MAIN_COMPLEXES_TO_STORE )
				break;
		}
		//		cout << "\tsubset=" << subset << "\t ei=" << result.ei << "\t phi=" << result.min_phi << endl;
		
		
	}
	
	return z;	
}



t_ei_result t_consciousness::highest_ei( int x1 )
// this functions returns the subsets with the highest bracket_ei
{
	t_ei_result z;
	z.x1 = x1;
	z.ei = -1.0;
	
	
	for( unsigned int subset=1; subset<=this->FULL_MASK; subset++ )
	{
		//		cout << "Starting subset=" << subset << endl;
		//OPTIMIZATION: only attempt all partitions to find the MIPs if it's possible
		//for this subset to have higher or equal PHI to the current highestphi
		double this_ei = this->ei(x1, subset);
		
		if( this_ei > z.ei ) {
			z.subsets.clear();
			z.ei = this_ei;
			z.subsets.push_back(subset);
		}
		else if( this_ei == z.ei ) {
			z.subsets.push_back(subset);
			
			// stop if we have the maximum number of main complexes
			if( z.subsets.size() >= MAXIMUM_NUMBER_OF_MAIN_COMPLEXES_TO_STORE )
				break;
		}
		//		cout << "\tsubset=" << subset << "\t ei=" << result.ei << "\t phi=" << result.min_phi << endl;
		
		
	}
	
	
	return z;	
}



double t_consciousness::I_A0_B1_GIVEN_ENV( const uint A0, const uint B1 )
// Calculates the information directly transferred from A0 -> B1
// Calculates I[A0:B1|env(A,B)] = 
// H[ A_0, e_1^{AB} ] - H[ e_1^{AB} ] - |A| + H[A_1|e_0^A e_1^A]
{
	if( A0 == B1 )
		return H_M1_GIVEN_ENV( A0 );
	
	// Although technically there's no reason why A and B can't overlap, 
	// In this code if we have an overlapping A and B we're probbly doing something
	// deeply wrong.
	assert( (A0 & B1) == 0);
	
	const uint ENV = env( A0 | B1 );
	
	if( ENV == 0 )
		return I_A0_B1( A0, B1 );
	
	
	double z = 0.0;
	
	// Calculate joint entropy of H[ A_0, e_1^{AB} ]
	
	z -= H_M1( ENV );
	z -= numunits_in_mask(A0 );
	z += H_M1_GIVEN_ENV( A0 );
	
	return z;	
}

#endif