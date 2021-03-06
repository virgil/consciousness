/*
 *  t_consc.h
 *  
 *
 *  Created by Virgil Griffith (Caltech) on 2/15/09.
 *  Core algorithms derived originally from Arend Hintze (Keck Graduate Institute).
 *  Partition Enumerator by Nicolas Chaumont (Keck Graduate Institute).
 *
 */

//GCC supports binary literals in C programs using the 0b prefix, pretty much like you would use 0x for hexadecimal literals. In the following example, we initialize an integer using a binary literal for its value: 
//int integer=0b10111001;

#ifndef T_CONSC_H
#define T_CONSC_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <math.h>
#include <time.h>
#include <algorithm>
#include "t_partition.h"
#include "t_subset.h"
#include "t_state.h"
#include "PartitionEnumerator.h"
#include <boost/foreach.hpp>
#include "helpers.h"


using namespace std;


struct t_transformation {
    int A, B, C, D, E;
    unsigned char operation;
    t_transformation(): A(-1), B(-1), C(-1), D(-1), E(-1), operation(' ') {}
};


// -----------------------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------------------

// returns the number of entries in the array
void subset2array( const bitmask subset, unsigned int* restrict z );
unsigned int mask2array( const bitmask inpmask, unsigned int* restrict z );


//void  binary2partition( FILE* f, t_partition* P );
double myround( double x, double precision = 0.0 );
long double myfloor( long double x );
long double myfloor( double x );
long double myceil( long double x );
long double myceil( double x );
unsigned int bin2dec(char* bin);

bool fequals( const double a, const double b );
double assert_bounds( double lower, double value, double upper );

struct t_cmdline_args {
    string networkfilename;
    unsigned int max_num_parts;
    unsigned int worker_index;
    unsigned int numworkers;
    string pgenerator_cmd;                  // shell command that outputs the partitions    
    string ofilename;                       // the output filename
    string bracketphi_logfile;              // output log file for the distribution of bracketphi
    string highestbracketphis_logfile;      // output log file for the distribution of highestbracketphis
    t_cmdline_args(): networkfilename(""), max_num_parts(UINT_NONE), worker_index(0), numworkers(1), pgenerator_cmd(""), ofilename(""), bracketphi_logfile(""), highestbracketphis_logfile("") {}
};


// what the psi( x1 ) function returns
class t_psi_result {
    
public:
    
    t_psi_result(): lowerbound(DBL_NONE), upperbound(DBL_NONE), ei(DBL_NONE), x1(0,0) {}
    
    double lowerbound;
    double upperbound;
    double ei;
    t_state x1;
};


class t_phi_result {

public:
    double min_phi;         // the minimum unnormalized phi.
//  long double phi_norm;
    double ei;              // the raw effective information ei(x1) or bracket_ei()
    vector< t_partition > MIPs;
    vector< double > PHIs;    // the ei(x/P) over each MIP
    long double mip_score;
    
    double atomic_IbD;
    
    // Variables that are sometimes defined.
    t_subset subset;    // the subset, if it's defined.  Default=0
    int x1;             // the x1 state, if it's defined. Default=-1;

    t_phi_result(): min_phi(DBL_NONE), mip_score(DBL_NONE), ei(DBL_NONE), x1(UINT_NONE), atomic_IbD(DBL_NONE) {}
    
    bool operator==(const t_phi_result &rhs) const
    {
    
        if( subset != rhs.subset || x1 != rhs.x1 )
            return false;
        
        if( ! fequals(min_phi, rhs.min_phi) )
            return false;

        if( ! fequals(atomic_IbD, rhs.atomic_IbD) )
            return false;
        
        
        for( int i=0; i<PHIs.size(); i++ )
        {
            if( MIPs[i] != rhs.MIPs[i] || ! fequals( PHIs[i], rhs.PHIs[i] ) ) 
                return false;
        }
        
        return true;
    }


    bool operator<(const t_phi_result &rhs) const
    {
        // sort by the min_phi
        if( min_phi > rhs.min_phi ||
           // sort by the number of units in the subset
           (fequals(min_phi, rhs.min_phi) && subset.size() > rhs.subset.size() ) ||
           // sort by ei
           (fequals(min_phi, rhs.min_phi) && subset.size() == rhs.subset.size() && (ei > rhs.ei) ) ||
           // sort by integer value of subset
           (fequals(min_phi, rhs.min_phi) && subset.size() == rhs.subset.size() && fequals(ei,rhs.ei) && subset < rhs.subset )
           )
            return true;

        else
            return false;


    }

        
    bool operator!=(const t_phi_result &rhs) const { return !(*this == rhs ); }
    
    
    void round() { 
        min_phi = myround(min_phi);
//      phi_norm = myround(phi_norm);
        ei = myround(ei);
        atomic_IbD = myround(atomic_IbD);
    }
};



class t_consciousness{
public:
    
    
    // for the PSI measure
    double ei( const t_state& restrict s1 );
    t_psi_result psi( const t_state x1 );
    double psi_lowerbound( const t_state& restrict x1 );
    double psi_upperbound( const t_state& restrict x1 );
    
    t_psi_result bracketpsi();
//    double bracketpsi();
    double bracketpsi_lowerbound();
    double bracketpsi_upperbound();
    


    // returns the perstate total ei.  This should be the UPPER BOUND of the perstate total_holism and total_integration.
	double total_perstate_ei( const t_partition& P, const unsigned int x1 );
	double perstate_ei( const bitmask S0mask, const unsigned int m1, const bitmask M1mask );	// for a part of a subset
	double perstate_ei( const unsigned int x1 );												// for the whole system

	vector<bitmask> all_s1s( const bitmask Smask );
	
    
    double prob_m0( const t_state& m0 );
    
    double prob_s1_given_m0( const t_state& restrict s1, const t_state& restrict m0 );
    double prob_s1_given_m0( const register unsigned int s1, const register bitmask Smask, const unsigned int mu0, const bitmask MUmask, const unsigned int MUsize );
    double prob_s1_given_m0__slower( const unsigned int s1, const bitmask Smask, const unsigned int mu0, const bitmask MUmask, const unsigned int MUsize );
    
    double prob_m0_given_s1( const t_state& m0, const t_state& s1 );

    
    double prob_m1_given_s1( const t_state& restrict m1, const t_state& restrict s1 );

    
    // These functions are for calculating new_holism over a subset
    double prob_s1( const unsigned int s1, const unsigned int partmask );
    double prob_s1( const t_state& s1 );
    
    double prob_m0_s1( const t_state& restrict m0, const t_state& restrict s1 );
    double prob_m0_s1( const unsigned int m0, const bitmask M0mask, const unsigned int s1, const bitmask S1mask );
    

    double prob_m0_given_s0_r1( const t_state& m0state, const t_state& s0state, const t_state& r1state );

    ////////////////////////////////////////////////////////////	
	// current scratch space
    ////////////////////////////////////////////////////////////	
	double H_A0_given_b1( const bitmask A0mask, const unsigned int b1, const bitmask B1mask );
	double H_A0_given_B0_c1( const bitmask A0mask, const bitmask B0mask, const t_state& restrict c1state );
	
	double I_A0_B1_equals_b1( const bitmask Amask, const bitmask Bmask, const unsigned int b1 );
	double I_A0_B1_equals_b1( const bitmask Amask, const t_state& b1 );

	
    ////////////////////////////////////////////////////////////
    // For network type 'neural'
    ////////////////////////////////////////////////////////////
    vector< double > node_thresholds;
    vector<vector< double> > weightmatrix;
	
    ////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////
    // For a cache of H[M_0|M_1]
    ////////////////////////////////////////////////////////////

    // cache of each possible H[M_0|M_1]
    double* restrict H_M0_GIVEN_M1_cache;
    
    
    ////////////////////////////////////////////////////////////
    // For network type 'circuit'
    ////////////////////////////////////////////////////////////
    vector< t_transformation > transformations;
    ////////////////////////////////////////////////////////////
    vector< string > network_comments;    // stores any comments in the .txt file
    string networktype;                 // can be either 'neural' or 'circuit'

    
    
    
    // Whether or not to use the total partition
    bool use_total_partition;
    
    ////////////////////////////////////////////////////////////
    // For saving the spectrum
    ////////////////////////////////////////////////////////////    
    int save_bracketphi_spectrum( string ofilename );
    int save_highestbracketphis_spectrum( string ofilename );
    
    bool neuron_fires_only_at_threshold;
    unsigned int max_num_parts; // default: this->numunits;
    unsigned int FULL_MASK;     // bit-mask for the entire system.  (= this->numstates-1)
    ////////////////////////////////////////////////////////////

    
    // Index is the x0, value is the x1.
    // Stores the x1 for a given x0.
    unsigned int* restrict states;

    // vector< int > x1_states;                      // vector of all valid x1 states
    map< int, double > H_X0_given_X1_equals_x1;    // pre-computed H[X0|X1=x1]

    //  H(X_1) -- the entropy of the X1 states
    double H_X1;

    // key: x1, value: #x0 states that lead to that x1
    map< unsigned int, unsigned int > existing_states;
    
    unsigned int numunits;
    unsigned int numstates;
    double H_X0_GIVEN_X1;
    unsigned int NUM_X1_STATES;

	// these two arrays are the replacements for the map of existing_states
	unsigned int* restrict x1_states;			// array of all valid x1 states	
	unsigned int* restrict num_x0s_to_this_x1;	// the number of x0s that went to this x1 state
	
    void load_transformations( t_cmdline_args &args, bool show = false);
    void make_all_states( t_cmdline_args &args, bool show = false);
    void show_rules(void);

    double entropy_of_part_given_x1( const unsigned int x1, const bitmask part_mask, unsigned int subset ); 
    double entropy_of_part_given_x1__WIRES( const unsigned int x1, const bitmask part_mask, const bitmask subset );
    
    double H_S0_given_S1_equals_s1( const unsigned x1, const bitmask subset );

    bool can_calc_by_perturbing_elements( unsigned int part_mask );

    
    // returns the atomic partition of a particular subset
    t_partition atomic_partition( unsigned int subset = 0 );
    
    double bracket_ei( const vector<vector<int> >& restrict, unsigned int subset );
    double bracket_ei( const t_partition& restrict P );
    double bracket_ei( unsigned int subset );    // calculates bracket_ei for total partition
    double bracket_ei( const t_subset& restrict S );    // calculates bracket_ei for total partition    

    /////////////////////////////////////////////////////////////////////////////
    // Functions for computing INFORMATION FLOW via Ay, Polani (2006)
    /////////////////////////////////////////////////////////////////////////////
    
    // TODO: remove these functions
    double I_A0_ARROW_B1_IMPOSING_S0( const bitmask Amask, const bitmask Bmask, const bitmask Smask );
    double I_A0_ARROW_B1_IMPOSING_s0( const bitmask Amask, const bitmask Bmask, const bitmask Smask, const unsigned int s0 );
    double prob_b1_IMPOSING_a0s0( const bitmask Bmask, const unsigned int b1,
                                  const bitmask ASmask, const unsigned int a0s0 );
    double prob_n1_GIVEN_a0s0( const unsigned int n1, const unsigned int n_index, 
                               const unsigned int a0s0, const bitmask ASmask, const unsigned int ASsize );
    double I_A0_ARROW_B1( const bitmask Amask, const bitmask Bmask );
	

    /////////////////////////////////////////////////////////////////////////////
    

    
    // computes the TOTAL CORRELATION among parts (assuming timestep=1) anticonditioned on the anticond_mask.
    
    double average_phi( const bitmask subset );
    // Bracketphi is defined as the average ei over the average MIP.    
    
    vector<t_partition> bracket_MIPs( unsigned int subset, bool debug = false );
    vector<t_partition> bracket_MIPs( const t_subset& restrict S, bool debug = false );

    double ei( const unsigned int s1, const t_partition&, const string force_version="UNITS" );      // calculates ei(x1/P) and ei(s1/P)  
    double ei( const unsigned int s1, const bitmask subset );   // calculates ei(s1)
    vector< t_partition > MIPs( int, unsigned int subset );
    
    t_phi_result bracketphi( unsigned int subset, bool require_atomic_partition = false );
    t_phi_result bracketphi( const t_subset& restrict S, bool require_atomic_partition = false ); 
    t_phi_result phi( const unsigned int x1, const bitmask subset = 0 );
    
    // number of distinct s1 states in subset S
    unsigned int num_s1_states( const bitmask S );
    unsigned int num_s1_states( const t_subset S );
    
    double normalization( const t_partition& restrict P );

    inline unsigned int size_of_smallest_part( const t_partition& restrict P );

    
    
    //misc
    int save_to_file(string ofilename, bool perstate, bool bracketphi, bool highestbracketphi, bool bracketMCs, bool higheststatephi, bool highestbracketei, bool higheststateei );

    string subset2str( const unsigned int subset );
    
    vector< t_phi_result > highest_bracketphis();
    vector< t_phi_result > highest_phis( const unsigned int x1 );

    long double normalized_ei( const bitmask x1, const t_partition& restrict P );
    long double bracket_normalized_ei( const t_partition& restrict P );
    
    vector< int > external_nodes_with_wires_into_node( unsigned int node_index, const unsigned int* restrict part, unsigned int numnodes, unsigned int subset = 0 );
    
    vector< t_phi_result > bracket_maincomplexes();

    
    string binary_repr( const unsigned int number, unsigned int numplaces=0 );
    string binrep( const unsigned int number, unsigned int numplaces=0 );   
    
    
    double H_M0( const t_subset& restrict S );
    double H_M0( const bitmask S );
    double H_M1( const bitmask S );

    
    double I_A0_B1( const bitmask A0, const bitmask B1 );
	double H_A0_B1( const bitmask A0, const bitmask B1 );
	
    double I_A0_B1_GIVEN_C0( const bitmask A0, const bitmask B1, const bitmask C0 );
    double I_A0_B1_GIVEN_C1( const bitmask A0, const bitmask B1, const bitmask C1 );
    double I_A0_B1_GIVEN_C0_D1( const bitmask A0, const bitmask B1, const bitmask C0, const bitmask D1 );
	
	
    double MIP_score( const bitmask x1, const t_partition& restrict P );
    double MIP_score( const t_partition& restrict P );
    
    bool break_MIPscore_ties( const t_partition& restrict P0, const t_partition& restrict P1 ); 
    
    inline unsigned int cachekey( unsigned int S0, unsigned int M1 );
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // OBSOLETE
    // FILE* get_partitions( unsigned int subset=0, unsigned int max_num_parts=0 );
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    t_consciousness( void );
    ~t_consciousness( void );

    
    
private:
	
	// this is just here so we can delete [] it and set it to NULL at the end of the program
	unsigned int* prob_s1_given_mu0__Vnodes;

    inline bitmask env( const bitmask incoming_mask, bitmask subset=0 );
    
	bool is_valid_m0( const unsigned int mu0, const bitmask MUmask );
	bool is_valid_m1( const unsigned int mu1, const bitmask MUmask );
	bool is_valid_m1( const t_state& restrict mu1 );
	bool is_valid_m0( const t_state& restrict mu0 );
	
    ////////////////////////////////////////////////////////////
    // For perturbing wires
    ////////////////////////////////////////////////////////////    
    double prob_n1_given_s0__anticond( const unsigned int n1, const unsigned int n_index, const unsigned int s0,
									  const bitmask S0mask, const unsigned int anticond_size );

    double prob_mu1_given_s0__anticond( const unsigned int mu1, const unsigned int s0, const bitmask M1mask, 
									   const bitmask S0mask, const unsigned int anticond_mask, 
									   const unsigned int M1size, const unsigned int S0size, 
									   const unsigned int anticond_size );

    
    double prob_mu1_given_mu0__anticond( const unsigned int mu1, const unsigned int mu0, const unsigned int partmask, const unsigned int partsize );
    double prob_mu0_given_mu1__anticond( const unsigned int mu0, unsigned int mu1, unsigned int partmask, unsigned int partsize );
    


    double prob_node_matches_mu1( const unsigned int part_mask, const unsigned int node_index, const unsigned int mu0, 
								 const unsigned int mu1, const unsigned int* restrict nodes, const unsigned int numnodes );

    vector< int > states_to_attempt__WIRES( unsigned int node_index, int mu0, 
										 const unsigned int* restrict part_nodes, unsigned int numnodes );   

    // For getting the partitions.
    bool get_next_partition( PartitionEnumerator* restrict pgen, t_partition* restrict P, const unsigned int* restrict nodemap, unsigned int numnodes );    
    
    inline void Xstate2MUstate( const unsigned int* restrict part, const unsigned int partsize, 
							   const unsigned int Xstate0, unsigned int& MUstate0 );
	
    inline void Xstate2MUstate( const unsigned int* restrict part, const unsigned int partsize, 
							   const unsigned int Xstate0, unsigned int& MUstate0, unsigned int Xstate1, unsigned int& MUstate1 );
	

    double H_M0_GIVEN_s1( const bitmask Mmask, const t_state& restrict s1 );
    double H_M0_GIVEN_S1( const bitmask M0, const bitmask S1 );
    double H_M0_GIVEN_M1( const bitmask M );

    
    double H_M0_GIVEN_M1__ELEMENTS( const bitmask part_mask );
    double H_M0_GIVEN_M1__WIRES( const bitmask part_mask );

    
    double H_M1_GIVEN_S0( const bitmask M, const bitmask S );
    double H_M1_GIVEN_M0( const bitmask partmask );


    // subsets
    inline unsigned int numunits_in_subset( const unsigned int, bool allow_zero = false );
    inline unsigned int numunits_in_mask( const unsigned int, bool allow_zero = false);
    
    string output_filename;

};

extern t_consciousness consciousness;



void sort_phi_results( vector< t_phi_result >& restrict results );
void sort_partitions( vector< t_partition >& restrict partitions );

#endif

