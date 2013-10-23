// MR-mapper.cpp -- The Mapper part of the MapReduce.
//

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
//#include <cstdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <math.h>
#include <time.h>

#define NORMALIZE

#include "t_consc.h"

using namespace std;

	
void print_usage() 
{
    cout << "================================================================" << endl;
    cout << "HOW TO RUN:" << endl;
    cout << "first arg: network.txt" << endl;
    cout << "second arg: directory (biparts|allparts|bitriparts)" << endl;
    cout << "================================================================" << endl;    
    exit(-1);
}


void StringExplode(string str, string separator, vector<string>* results) {
	//    int found;
	//    found = str.find_first_of(separator);
    string::size_type found = str.find( separator );
	
    while(found != string::npos){
        if(found > 0){
            results->push_back(str.substr(0,found+separator.length()));
        }
		//        str = str.substr(found+1);
        str = str.substr(found+separator.length() );
        found = str.find( separator );
		//        found = str.find_first_of(separator);
    }
    if(str.length() > 0){
        results->push_back(str);
    }
}

void ParsePartsString( string line, vector<vector<int> >* results )
/* This function takes a JSON partsstring, and converts it to a vector<vector<int> >'s. */
{
	vector<string> pieces;
    pieces.clear();
    (*results).clear();
    
    
    StringExplode( line, "], ", &pieces );
    for( int i=0; i<pieces.size(); i++ )
    {
        string part = pieces[i];
		//    	cout << "part='" << part << "'" << endl;
		
        //remove ']'s
        string::size_type found = part.find( "]" );
        while(found != string::npos) {
            part.erase(found, 1);
            found = part.find( "]" );
        }
		
		//				cout << "\tpart=" << part << endl;                
        //remove '['s
        found = part.find( "[" );
        while(found != string::npos) {
            part.erase(found, 1);
            found = part.find( "[" );
        }
		
        //if there's a ", " at the end, remove it.
        if( part[part.size()-2] == ',' && part[part.size()-1] == ' ' ) {
            part.erase(part.size()-2, 2);
        }
		
        // a single entry of results
        vector<int> result;
        result.clear();
        
        found = part.find( ", " );
        while( found != string::npos) {
			//            cout << "\t\tpush_back()" << (part.substr(0,found)) << endl;
            result.push_back( atoi( (part.substr(0,found)).c_str()) );
            part = part.substr(found+2);
            found = part.find( ", " );
        }
        
        result.push_back( atoi(part.c_str()) );
		//        cout << "\t\tpush_back()" << part << endl;
		
        (*results).push_back( result );
        
    }
	
}


void print_statistics()
{

	cout << "===Statistics=========================================================" << endl;    	

	cout << "States:" << endl;
	for( int i=0; i<consciousness.states.size(); i++ ) {
		cout << "\tstates[" << i << "]=" << consciousness.states[i] << endl;
	}

	cout << "state_map: " << endl;
	consciousness.show_vec_vec_int(consciousness.state_map);

	cout << "existing_states" << endl;
	for( map<int,int>::iterator it=consciousness.existing_states.begin(); it != consciousness.existing_states.end(); it++ ) {
		int x1 = (*it).first;
		int numx0s = (*it).second;
		cout << x1 << " -> " << numx0s << endl;
	}
	cout << "===END=Statistics=========================================================" << endl;    		
}

void print_averages( )
{
    cout << "===AVERAGES=========================================================" << endl;    
	vector<vector <int> > average_MIP = consciousness.average_MIP( partitions_filename );
	string MIPstr = consciousness.partition2str( average_MIP );
//    cout << "H[X0|X1]=" << consciousness.H_X0_GIVEN_X1 << endl;
    cout << "average ei       \t=" << consciousness.average_ei() << endl;
	cout << "average MIP      \t=" << MIPstr << endl;
	cout << "\t average_ei( average_MIP )=" << consciousness.average_ei( average_MIP ) << endl;

	vector<vector<vector <int> > > average_MIPs = consciousness.average_MIPs( partitions_filename );
	if( average_MIPs.size() >= 2 )	{
		for( int i=0; i<average_MIPs.size(); i++ )
			cout << "\tavrMIPs[" << i << "]=" << consciousness.partition2str( average_MIPs[i] ) << "\taverage_ei(MIPs[" << i << "])=" << consciousness.average_ei( average_MIPs[i] ) << endl;
	}
	
    cout << "average PHI (raw)\t=" << consciousness.average_phi( partitions_filename ) << endl;
    cout << "===END=AVERAGES=====================================================" << endl;
}

t_cmdline_args read_args_from_stdin() {
    t_cmdline_args z;
	z.worker_index = 0;
	z.numworkers = 0;
	z.max_num_parts = 0;
	
	// are we reading things from STDIN ?
	if( ! feof(stdin) )	{
		char temp[1024];
		char temp3[400];
		char temp4[1024];
		fgets( temp, 1024, stdin );
		
		int n = sscanf(temp, "%d %d %s %d %s", &z.worker_index, &z.numworkers, &temp3, &z.max_num_parts, &temp4 );
		assert( n == 5 );
		assert( 0 <= z.worker_index && z.worker_index <= z.numworkers );

		z.networkfilename = string(temp3);
		z.pgenerator_cmd = string(temp4);

		// check that networkfilename and pgenerator_cmd are valid files...
		FILE * tempfile;
		tempfile = fopen( z.networkfilename.c_str(), "r" );
		if(tempfile == NULL) {
			perror ("Error opening networkfilename.");
			exit(1);
			}
		fclose(tempfile);

//		cerr << "Opening " << z.pgenerator_cmd << endl;		
		tempfile = fopen( z.pgenerator_cmd.c_str(), "r" );
		if(tempfile == NULL) {
			perror ("Error opening generator.");
			exit(1);
			}
		fclose(tempfile);
	}

	else
        print_usage();

	return z;

}

int main(int argc, char* argv[])
{
	
    // Load the network, get the filenames.
	t_cmdline_args args;

    args = read_args_from_stdin();
	
	consciousness.load_transformations(args, false);
	
	
	cerr << "networkfilename=" << args.networkfilename << endl;
	cerr << "max_num_parts=" << args.max_num_parts << endl;	
	cerr << "pgenerator_cmd=" << args.pgenerator_cmd << endl;
	cerr << "pfilename=" << args.pfilename << endl;	
	cerr << "worker_index=" << args.worker_index << endl;
	cerr << "numworkers=" << args.numworkers << endl;	

	//    cout << "1. parts filename   \t = '" << partitions_filename << "'" << endl;    
    
	//    cout << "2. parts filename   \t = '" << partitions_filename << "'" << endl;    
	
    // run the input states.
    consciousness.make_all_states(false);
    partitions_filename = (char*) args.pfilename.c_str();
		
//	print_averages();
//	print_statistics();
//	print_perstate();
//	print_statistics();

//	print_average_singleline( argv );
	consciousness.MAPREDUCE_WORKER_MIPs( args.worker_index, args.numworkers, args.pgenerator_cmd );
	
//	do_worker( args.worker_index, args.numworkers );
	

/*    
	clock_t t1,t2;
	t1=clock();    
	t2=clock();
	printf("time: %f\n",((float)t2-(float)t1)/(float)CLOCKS_PER_SEC);	
*/
	
	return 0;
}
