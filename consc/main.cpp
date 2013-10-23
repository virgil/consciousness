#ifndef MAIN_CPP
#define MAIN_CPP

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>

#define NDEBUG
#include <assert.h>

#include <math.h>
#include <iomanip>

#include "t_consc.h"
#include "t_state.h"

using namespace::std;

//////////////////////////////////////////////////////////////////////////////////////////
// CONFIGURABLE VARIABLES
//////////////////////////////////////////////////////////////////////////////////////////

bool SAVETOFILE__DO_BRACKETPHI = true;
bool SAVETOFILE__DO_BRACKET_MAINCOMPLEXES = false;
bool SAVETOFILE__DO_HIGHEST_BRACKETPHI = false;
bool SAVETOFILE__DO_HIGHEST_BRACKET_EI = false;

bool SAVETOFILE__DO_PERSTATE_PHI = true;
bool SAVETOFILE__DO_HIGHEST_PERSTATE_PHI = false;
bool SAVETOFILE__DO_HIGHEST_PERSTATE_EI = false;

bool VERBOSE_DEBUG_MODE = false;

//////////////////////////////////////////////////////////////////////////////////////////
// FUNCTION PROTOTYPES
//////////////////////////////////////////////////////////////////////////////////////////
void print_usage();


t_cmdline_args read_cmdline_args( int argc, char* argv[], bool show=false );
//////////////////////////////////////////////////////////////////////////////////////////

void print_usage() 
{
    cout << "================================================================" << endl;
    cout << "HOW TO RUN:" << endl;
	cout << "./consciousness [opts] path/to/networkspec.txt [outputfilename] " << endl << endl;

    cout << " \t first arg: network.txt" << endl;
	cout << " \t second arg: outputfilename (string)  - filename to send the result to.  default: <firstarg>.out" << endl;
		
	cout << "Options: " << endl;
    cout << " \t -k max_num_parts (int)       - Maximum Number of Parts. 1<=K<=#nodes. default: K=#nodes" << endl;
	cout << " \t -s file.log                  - log distribution of <\\phi> of whole network" << endl;
	cout << " \t -S file.log                  - log spectrum of the subsets" << endl;
	cout << " \t                                if multiple MCs becomes file-1.logm file-2.log, etc." << endl;
	cout << " \t -v                           - turn ON verbose mode for debugging. " << endl;	
	cout << " \t [+]code                      - turn ON savetofile for this" << endl;
	cout << " \t [-]code                      - turn OFF savetofile for this" << endl;
	cout << " \t codes: b=<ϕ>  \t B=<Φ> (main complex)" << endl;
    cout << "================================================================" << endl;    
    exit(0);
}


int main(int argc, char* argv[])
{
	
	
	
	t_cmdline_args args = read_cmdline_args( argc, argv, VERBOSE_DEBUG_MODE );
	cerr << "- system: " << args.networkfilename << endl;
	cerr.flush();
	
	
	if( VERBOSE_DEBUG_MODE ) {
		cout << "VERBOSE MODE is ON" << endl;
		cout.flush();
		cerr.flush();
	}

	t_consciousness consciousness;

	
	// Don't mess with this.  This allows the partitions generator to be re-used.	
	consciousness.load_transformations( args, VERBOSE_DEBUG_MODE );
	
	if( VERBOSE_DEBUG_MODE ) {
		cout << "networkfilename=" << args.networkfilename << endl;
		cout << "max_num_parts=" << args.max_num_parts << endl;	
		cout << "ofilename=" << args.ofilename << endl;
		cout.flush();
		cerr.flush();
	}

	
	// run the input states.
	consciousness.make_all_states( args, VERBOSE_DEBUG_MODE );

	cout.flush();
	
	if( VERBOSE_DEBUG_MODE )
	{		
		cout << "x0->x1 transitions:" << endl;
		for( int x0=0; x0<consciousness.numstates; x0++ )
		{
			
//			t_subset temp_x0( consciousness.numunits, x0 );
//			t_subset temp_x1( consciousness.numunits, consciousness.states[x0] );
			
            t_state temp_x0( consciousness.FULL_MASK, x0 );
            t_state temp_x1( consciousness.FULL_MASK, consciousness.states[x0] );            
            
            cout << "\t x0: " << temp_x0.str( consciousness.FULL_MASK ) << " -> " << temp_x1.str( consciousness.FULL_MASK ) << endl;
            
//			cout << "\t x0: " << binary_repr(x0,consciousness.numunits) << " -> " << binary_repr(consciousness.states[x0], consciousness.numunits) << endl;
//			cout << "\t" << x0 << " -> " << consciousness.states[x0] << endl;			
		}
	}
	
	cout.flush();	

	// saving distribution of bracketphi?
	if( ! args.bracketphi_logfile.empty() )
	{
		cerr << "- Saving distribution of bracketphi to: " << args.bracketphi_logfile << " ...";
		cerr.flush();
		consciousness.save_bracketphi_spectrum( args.bracketphi_logfile );
		cerr << "done." << endl;
	}

	// saving distribution of bracketphi?
	if( ! args.highestbracketphis_logfile.empty() )
	{
		cerr << "- Saving distribution of main complexes to: " << args.highestbracketphis_logfile << " ...";		
		cerr.flush();
		consciousness.save_highestbracketphis_spectrum( args.highestbracketphis_logfile );	
		cerr << "done." << endl;		
	}


	if( SAVETOFILE__DO_BRACKETPHI || SAVETOFILE__DO_HIGHEST_BRACKETPHI || SAVETOFILE__DO_BRACKET_MAINCOMPLEXES || SAVETOFILE__DO_PERSTATE_PHI || SAVETOFILE__DO_HIGHEST_PERSTATE_PHI || SAVETOFILE__DO_HIGHEST_BRACKET_EI || SAVETOFILE__DO_HIGHEST_PERSTATE_EI )
	{
		int success = consciousness.save_to_file( args.ofilename, 
			SAVETOFILE__DO_PERSTATE_PHI, SAVETOFILE__DO_BRACKETPHI, 
			SAVETOFILE__DO_HIGHEST_BRACKETPHI, SAVETOFILE__DO_BRACKET_MAINCOMPLEXES, SAVETOFILE__DO_HIGHEST_PERSTATE_PHI,
			SAVETOFILE__DO_HIGHEST_BRACKET_EI, SAVETOFILE__DO_HIGHEST_PERSTATE_EI );
		
		 if( success < 0 ) {
			cerr << "* Error!  Could not write to '" << args.ofilename << "'" << endl;
		}
	}

	
	return 0;
}

	
t_cmdline_args read_cmdline_args( int argc, char* argv[], bool show )
/*
 This function returns the arguments as strings including the:
 1) networkfilename
 2) maximum number of parts
 */
{

	// make a vector of string of args.  This will be much easier to work with than directly working with char* argv[]
	vector<string> args;
	vector<string> params;
	args.reserve( argc );

	for( int i=1; i<argc; i++ ) {
		args.push_back( string(argv[i]) );

	}

	if( args.empty() )
		print_usage();

	if( args.empty() ) {
		string default_network_filename = "/consciousness/networks/p-theta2.txt";
		cout << "Using network '" << default_network_filename << "'..." << endl;
		cout.flush();
		args.push_back( default_network_filename );
	}

	
	// initialize cmdline_args with empty values
	t_cmdline_args z;


	///////////////////////////////////////////////////////////////////////////////////////////
	// PROCESS ALL OPTIONS
	///////////////////////////////////////////////////////////////////////////////////////////	
	
	// go through each arg in args until none begin with a dash
	while( args.size() )
	{

		string opt=args[0], value="";
		cout.flush();
		
		if( args.size() >= 2 )
			value = args[1];


		if( opt[0] == '+' )
		{
			//turn ON <\phi>
			if( opt[1] == 'b' ) {
				SAVETOFILE__DO_BRACKETPHI = true;
				assert( SAVETOFILE__DO_BRACKETPHI );
				cout << "* PROCESSING <phi> (whole system)" << endl;
				// only delete the first part because there is no 2nd argument for -b.
				args.erase( args.begin() );
			}
			
			//turn ON <\Phi>
			else if( opt[1] == 'B' ) {				
				SAVETOFILE__DO_BRACKET_MAINCOMPLEXES = true;
				assert( SAVETOFILE__DO_BRACKET_MAINCOMPLEXES );
				cout << "* PROCESSING <Phi> (all main complexes)" << endl;
				
				// only delete the first part because there is no 2nd argument for -B.
				args.erase( args.begin() );
			}			
		}
		
		else if( opt[0] == '-' )
		{
			cout << "found option: " << opt << endl;
						
			// doing -k?
			if( opt[1] == 'k' ) {
				z.max_num_parts = atoi( value.c_str() );
//				cout << " \t setting z.max_num_parts = " << z.max_num_parts << endl;
				assert( z.max_num_parts > 0 );
				args.erase( args.begin(), args.begin()+2 );				
			}

			else if( opt[1] == 's' ) {
				z.bracketphi_logfile = value;				
				assert( z.bracketphi_logfile != "" );
				args.erase( args.begin(), args.begin()+2 );
				
				cout << "z.bracketphi_logfile=" << z.bracketphi_logfile << endl;
				cout.flush();
			}

			else if( opt[1] == 'S' ) {
				z.highestbracketphis_logfile = value;				
				assert( z.highestbracketphis_logfile != "" );
				args.erase( args.begin(), args.begin()+2 );				
			}
			
			// turn on debug mode
			else if( opt[1] == 'v' ) {
				VERBOSE_DEBUG_MODE = true;
				// only delete the first part because there is no 2nd argument for -v.
				args.erase( args.begin() );
			}
			
			//turn OFF <\phi>
			else if( opt[1] == 'b' ) {
				SAVETOFILE__DO_BRACKETPHI = false;
				assert( ! SAVETOFILE__DO_BRACKETPHI );
				cout << "* SKIPPING <phi> (whole system)" << endl;
				
				// only delete the first part because there is no 2nd argument for -b.
				args.erase( args.begin() );
			}
			
			//turn OFF <\Phi>
			else if( opt[1] == 'B' ) {
				SAVETOFILE__DO_BRACKET_MAINCOMPLEXES = false;
				assert( ! SAVETOFILE__DO_BRACKET_MAINCOMPLEXES );
				cout << "* SKIPPING <Phi> (all main complexes)" << endl;	
				// only delete the first part because there is no 2nd argument for -B.
				args.erase( args.begin() );

			}
			
			else {
				cerr << "* Error: unrecognized option: " << opt << endl;
				print_usage();
			}
			
			// delete the first two entries in args

		}
		
		else {
			params.push_back( args[0] );
			args.erase( args.begin() );
		}
		
	}

	if( show )
		cout << "# params=" << params.size() << endl;
	
	///////////////////////////////////////////////////////////////////////////////////////////
	// process the input filename and the .out filename
	///////////////////////////////////////////////////////////////////////////////////////////	
	if( params.size() == 0 ) {
		cerr << "* Error: Must specify an input network filename" << endl;
		print_usage();
	}
	else if( params.size() == 1 || params.size() == 2 )
	{
		// do the network filename
		z.networkfilename = params[0];
		params.erase( params.begin() );

		// Set the ofilename if it's there
		if( params.size() == 1 ) {
			z.ofilename = params[0];
			params.erase( params.begin() );
		}

	}
	else {
		cerr << "* Error: You specified too many arguments.  Maximum is two." << endl;
		print_usage();
	}

	
	if( z.ofilename.empty() ) {
		z.ofilename = z.networkfilename + string(".out");
	}
	
	
	///////////////////////////////////////////////////////////////////////////////////////////
	// CHECK THAT FILES EXIST
	///////////////////////////////////////////////////////////////////////////////////////////		
	FILE* f = fopen( z.networkfilename.c_str() , "r");
	if( f == NULL ) {
		cerr << "* Error: network filename '" << z.networkfilename << "' doesn't exist." << endl;
		exit(-1);
	}
	fclose(f);
	
	// the args should be empty at this point
	assert( args.empty() );
	assert( ! z.networkfilename.empty() );
	assert( ! z.ofilename.empty() );

	
	return z;
}


#endif