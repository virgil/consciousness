#include <iostream>
#include <assert.h>
#include "boost/foreach.hpp"
#include <vector>

#ifndef foreach
	#define foreach BOOST_FOREACH
#endif

#define restrict __restrict__


using namespace std;
typedef unsigned int uint;

//unsigned int x2mu_smallex();
unsigned int FULL_MASK = (1 << 9)-1;



vector<unsigned int> to_vector(unsigned int S)
// returns all of the bits that are turned ON in the subset
{
	assert( S <= FULL_MASK );
	// clear the old part
	vector<unsigned int> z;
	z.clear();
	
	//	unsigned int power = (unsigned int) floor( log2( subset ) );
	unsigned int power=0;
	
	do	{
		// this only works because we're starting at i=0
		if( (S&1) )
			z.push_back(power);
		
		S >>= 1;
		power += 1;
	}while( S );
	
	return z;
}


unsigned int numunits( unsigned int x ) {
	unsigned int z = 0;

	do	{
		if( x&1 )
			z += 1;
		
		x >>= 1;
	}while( x );

	return z;
}

string str( unsigned int x, unsigned int custom_full_mask=0 )
// returns a string version of S, passing an optional full_mask
{
	
	if( custom_full_mask == 0 )
		custom_full_mask = numunits(FULL_MASK);

	string z = "";
	z.reserve( numunits(x) );
	
	//	int length = this->size();
	
	for( int i=0; i<custom_full_mask; i++ )
	{
		if( ((x>>i)&1) == true )
			z.insert(0, string("1"));
		else
			z.insert(0, string("0"));
		
	}
	
	return z;
}



uint x2mu_naive( unsigned int Xval, unsigned int Xmask, uint MUmask, vector<unsigned int> Xnodes, vector<unsigned int> MUnodes )
{
	
	// remove the bits of X not in the mask	
	assert( Xval <= Xmask );

	unsigned int Xsize = Xnodes.size();	
	unsigned int MUval = Xval >> Xnodes[Xsize-1];
	
	
	// for all future bits, first shift to the RIGHT, then add the bit on the end.
	for(int i=Xsize-2; i>=0; i--) {
//		MUval = (MUval<<1) | ((Xval >>Xnodes[i])&1);
		MUval <<= 1;
		MUval |= (Xval >>Xnodes[i])&1;
	}	
	
	return MUval;
}


uint x2mu_exchg( uint Xval, uint Xmask, uint MUmask, uint XOUTmask, uint MUINmask, vector<uint> XOUTnodes, vector<uint> MUnodes )
{

	assert( XOUTnodes.size() == MUnodes.size() );

	// For Xvals inside MU, pass them along.
	unsigned int MUval = Xval & MUINmask;
	
//	cout << "\t inside MUval=" << str(MUval) << endl;
	
	// For Xvals outside MU, get their value and shift by a MUIN val.
	for( int i=0; i<XOUTnodes.size(); i++ )
		MUval |= ((Xval >> XOUTnodes[i])&1) << MUnodes[i];
		
	
	return MUval;
}

int main()
{
	int inp[4] = { 0, 1, 4, 8 };

	unsigned int Xmask = 0;

	for( int i=0; i<4; i++ )
		Xmask |= (1<<inp[i]);

	unsigned int Xmask_size = numunits( Xmask );
	vector<unsigned int> Xnodes = to_vector( Xmask );

	//////////////////////////////////////////////////////////////////////////////	
	// set the MU stuff
	//////////////////////////////////////////////////////////////////////////////
	unsigned int MUmask = (1 << Xmask_size) - 1;
	vector<unsigned int> MUnodes = to_vector(MUmask);
	unsigned int MUsize = numunits(MUmask);		
	const unsigned int MUd = 1 << MUsize;		//d=2**part.size()
	//////////////////////////////////////////////////////////////////////////////	

	//////////////////////////////////////////////////////////////////////////////	
	// set the Xoutside stuff
	//////////////////////////////////////////////////////////////////////////////
	unsigned int XOUTmask = Xmask & ~MUmask;
	vector<uint> XOUTnodes = to_vector( XOUTmask );
	unsigned int XOUTsize = numunits(XOUTmask);

	unsigned int MUINmask = Xmask & MUmask;
	unsigned int MUdest = MUmask & ~MUINmask;
	vector<uint> MUdest_nodes = to_vector( MUmask & ~MUINmask );
//	vector<uint> MUINnodes = to_vector( MUINmask );
	//////////////////////////////////////////////////////////////////////////////	
	
	
	cout << "FULLmask	\t=" << str(FULL_MASK) << endl;	
	cout << "Xmask		\t=" << str(Xmask) << endl;
	cout << "MUmask		\t=" << str(MUmask) << endl;
	cout << "XOUTmask	\t=" << str(XOUTmask) << endl;
	cout << "MUINmask	\t=" << str(MUINmask) << endl;
	cout << "MUdest		\t=" << str(MUdest) << endl;
	
	for( int Xval=0; Xval<=FULL_MASK; Xval++ )
	{
		uint Xfiltered = Xval & Xmask;
		uint MUval_naive = x2mu_naive( Xfiltered, Xmask, MUmask, Xnodes, MUnodes );
		

		uint MUval_exchg = x2mu_exchg( Xfiltered, Xmask, MUmask, XOUTmask, MUINmask, XOUTnodes, MUdest_nodes );		
		
		
		cout << "\t X=" << str(Xval) << " -> " << str(Xfiltered) << " --> " << str(MUval_naive);
		cout << " --> " << str(MUval_exchg) << endl << endl;
		
		assert( MUval_naive == MUval_exchg );
	}

	return 0;
}