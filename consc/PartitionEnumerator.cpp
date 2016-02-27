#ifndef PARTITION_ENUMERATOR_CPP
#define PARTITION_ENUMERATOR_CPP

//#define NDEBUG
#include <assert.h>
#include <iostream>
#include "PartitionEnumerator.h"
#include <stdio.h>
#include <string.h>


using namespace std;

PartitionEnumerator::PartitionEnumerator(unsigned short initial_network_size, bool inp_skip_total_partition )
{
	outputted_last_partition = false;
	network_size = initial_network_size;

	assert( network_size <= 32 );


	partition			= new unsigned int[network_size + 1];
	partition[0]		= 1;
	memset( partition+1, 0, network_size * sizeof(unsigned int) );
	assert( partition[network_size] == 0 );


	partition_sizeList	= new unsigned int[network_size];
	partition_sizeList[0] = 1;


	element_list		= new unsigned int[network_size];
	memset( element_list, 0, network_size * sizeof(unsigned int) );


	node_posList		= new unsigned int[network_size];
	for(unsigned int pos=0; pos<network_size; ++pos) node_posList[pos] = 1;


	partition_increaseList	= new unsigned int[network_size];
	memset( partition_increaseList, 0, network_size * sizeof(unsigned int) );
	partition_increaseList[0] = 1;

	element_countList	= new unsigned int[network_size];
	untouched_nodes = 0;

	if( inp_skip_total_partition ) {
//		cout << "SKIPPING TOTAL PARTITION" << endl;
		this->skip_total_partition = true;

		const unsigned int* temp = this->nextPartition();

		// assert that this is the total partition
		assert( temp[0] == 1 );

		// again, assert that this is the total partition
		assert( temp[1] == ((1<<network_size) - 1) );

	}

}


PartitionEnumerator::~PartitionEnumerator()
{
	if(partition)				delete[] partition;
	if(partition_sizeList)		delete[] partition_sizeList;
	if(element_list)			delete[] element_list;
	if(node_posList)			delete[] node_posList;
	if(partition_increaseList)	delete[] partition_increaseList;
	if(element_countList)		delete[] element_countList;
}


unsigned int PartitionEnumerator::size()
{
	unsigned int partition_count = 0;

	switch(network_size)
	{
	case 1:
		partition_count = 1;
		break;
	case 2:
		partition_count = 2;
		break;
	case 3:
		partition_count = 5;
		break;
	case 4:
		partition_count = 15;
		break;
	case 5:
		partition_count = 52;
		break;
	case 6:
		partition_count = 203;
		break;
	case 7:
		partition_count = 877;
		break;
	case 8:
		partition_count = 4140;
		break;
	case 9:
		partition_count = 21147;
		break;
	case 10:
		partition_count = 115975;
		break;
	case 11:
		partition_count = 678570;
		break;
	case 12:
		partition_count = 4213597;
		break;
	case 13:
		partition_count = 27644437;
		break;
	case 14:
		partition_count = 190899322;
		break;
	case 15:
		partition_count = 1382958545;
		break;
	default:
		std::cout << "Number of different partitions is too big and does not fit in an unsigned int." << std::endl;
	};


	return partition_count;
}


void PartitionEnumerator::reset(unsigned short new_network_size, unsigned short maxnumparts )
{
	this->resize( new_network_size );
}


void PartitionEnumerator::resize(unsigned short new_network_size)
{

	this->outputted_last_partition = false;
	this->untouched_nodes = 0;

	// if the sizes are THE SAME, simply re-initialize everything so nextPartition() starts over
	// This saves us from re-creating all of the memory.
	if( new_network_size == this->network_size )
	{
		//		cout << "Same-size resize!()  N=" << network_size << endl;
		//////////////////////////////////////////////////////////////////////////////////////////
		assert( partition != NULL );
		memset( partition+1, 0, (network_size) * sizeof(unsigned int) );
		partition[0]		= 1;

		//////////////////////////////////////////////////////////////////////////////////////////
		assert( partition_sizeList != NULL );
		partition_sizeList[0] = 1;

		//////////////////////////////////////////////////////////////////////////////////////////
		assert( element_list != NULL );
		memset( element_list, 0, network_size * sizeof(unsigned int) );

		//////////////////////////////////////////////////////////////////////////////////////////
		assert(node_posList != NULL );
		for(unsigned int i=0; i<network_size; ++i)
			node_posList[i] = 1;

		//////////////////////////////////////////////////////////////////////////////////////////
		assert( partition_increaseList != NULL );
		memset( partition_increaseList, 0, network_size * sizeof(unsigned int) );
		partition_increaseList[0] = 1;


		//////////////////////////////////////////////////////////////////////////////////////////
		assert( element_countList != NULL );
		//////////////////////////////////////////////////////////////////////////////////////////

		return;
	}

	this->network_size = new_network_size;

	if(partition)
	{
		delete[] partition;
		partition = NULL;
	}
	partition = new unsigned int[network_size+1];
	partition[0]		= 1;
	for(unsigned int pos=1; pos<network_size+1; ++pos) partition[pos] = 0;


	if(partition_sizeList)
	{
		delete[] partition_sizeList;
		partition_sizeList = NULL;
	}
	partition_sizeList = new unsigned int[network_size];
	partition_sizeList[0] = 1;


	if(element_list)
	{
		delete[] element_list;
		element_list = NULL;
	}
	element_list = new unsigned int[network_size];
	for(unsigned int p=0; p<network_size; ++p) element_list[p] = 0;


	if(node_posList)
	{
		delete[] node_posList;
		node_posList = NULL;
	}
	node_posList = new unsigned int[network_size];
	for(unsigned int pos=0; pos<network_size; ++pos) node_posList[pos] = 1;


	if(partition_increaseList)
	{
		delete[] partition_increaseList;
		partition_increaseList = NULL;
	}
	partition_increaseList = new unsigned int[network_size];

	partition_increaseList[0] = 1;
	for(unsigned int pos=1; pos<network_size; ++pos) partition_increaseList[pos] = 0;


	if(element_countList)
	{
		delete[] element_countList;
		element_countList = NULL;
	}
	element_countList = new unsigned int[network_size];


	if( this->skip_total_partition ) {
		const unsigned int* temp = this->nextPartition();

		// assert that this is the total partition
		assert( temp[0] == 1 );

		// again, assert that this is the total partition
		assert( temp[1] == ((1<<network_size) - 1) );
	}

}



const unsigned int* PartitionEnumerator::firstPartition()
{
	// Initialize the partition:
	partition[0] = 1;						// Partition size is 1
	partition[1] = (1 << network_size) - 1;	// One single element with all the nodes in it

	// Initialize the partition increase history
	partition_increaseList[1] = 1;			// Increased by 1 the first time
	memset(&partition_increaseList[0] + 2, 0, sizeof(unsigned int) * network_size+1);	// No subsequent increase

	// Initialize the partition size after addition
	partition_sizeList[0] = 1;

	// Initialize the node addition history (all additions at pos zero)
	memset(&node_posList[0], 0, sizeof(unsigned int) * network_size);

	// Only 1 node to modify next time
	untouched_nodes = network_size - 1;

	element_list[0]++;

	return partition;
}


const unsigned int* PartitionEnumerator::currentPartition()
{
	return partition;
}


const unsigned int* PartitionEnumerator::nextPartition()
{
	if( this->outputted_last_partition )
	{
		// the <= network_size IS IMPORTANT because the 0'th entry is the number of entries
		assert( partition[0] == this->network_size || partition[0] == 0 );

		// set the partition to all zeros.
		memset( partition, 0, sizeof(unsigned int) * (network_size+1) );

//		for( int i=0; i<=network_size; i++ )
//			assert( partition[i] == 0 );


		return partition;
	}


	// Revert the outdated modifications made to the partition
	for(unsigned int d=untouched_nodes; d<network_size; ++d)
	{
		temp = network_size - 1 - d + untouched_nodes;
		partition[node_posList[temp]]	&= ~(1 << temp);	// Remove the current node
		partition[0] -= partition_increaseList[temp];	// Update partition size
	}

	// All modifications are relevant to the current partition. We now can add nodes to it
	for(unsigned int d=untouched_nodes; d<network_size; ++d)
	{
		if(d)	element_countList[d] = partition_sizeList[d-1] + 1;
		else	element_countList[d] = 1;

		// Should we add an element to this partition?
		if(element_list[d] == element_countList[d]-1)
		{
			if(d)	partition_sizeList[d] = partition_sizeList[d-1] + 1;
			else	partition_sizeList[d] = 1;

			// Add new node to the new element
			// Equivalent to partition.push_back(1 << d)
			partition[0]++;
			partition[partition[0]] = 1 << d;

			// Track the changes in the partition so that we can revert them later
			partition_increaseList[d] = 1;	// The partition grew by 1
			node_posList[d] = partition[0];	// Record where the node was added
			untouched_nodes++;				// We have a new valid hierarchical level
		}
		// Just need to add the new node to the d-th partition element
		else
		{
			if(d)	partition_sizeList[d] = partition_sizeList[d-1];
			else	partition_sizeList[d] = 1;

			// Add the new node to an existing element
			partition[element_list[d]+1] |= (1 << d);

			// Track the changes in the partition so that we can revert them later
			partition_increaseList[d] = 0;			// The partition hasn't grown
			node_posList[d] = element_list[d]+1;	// Record where the node was added
			untouched_nodes++;						// We have a new valid hierarchical level
		}


		// We reached a leaf, increment the number of elements we've added
		if(d == network_size-1)
		{
			element_list[d]++;
			untouched_nodes--;
		}

		// If the current depth is saturated, then find the next branch in the enumeration tree
		// that is not saturated
		if(element_list[d] == element_countList[d])
		{
			// -----------------------------------------------------------------------------------------------
			// The current depth is saturated, i.e. no more partition can be enumerated from this parent depth.
			// We have to 'roll back' in the enumeration tree to the next node that has partitions left
			// to enumerate.
			// -----------------------------------------------------------------------------------------------
			unsigned int offset = 0;
			// Go up the hierarchy stepwise
			while((element_list[d-offset] == element_countList[d-offset]) && (d != offset))
			{
				// Clear changes made to the current node
				element_list[d-offset]--;	// Reset the element list locus to a legal value
				untouched_nodes--;			// This node should be reassigned

				// Update the next level to allow cascading updates if necessary
				element_list[d-offset-1]++;

				offset++;
			}

			// If path went up to the root, we are at the last partition. Correct the first node so that
			// the element_list remains valid if the user wants this partition again
			if(d == offset)
			{
				element_list[0] = 0;
			}
			// If we didn't reach the root, then we found a node from which to enumerate. We just have to
			// reset those entries in element_list that correspond to the depths below that node.
			else
			{
				for(unsigned int depth_to_reset=d-offset+1; depth_to_reset<network_size; ++depth_to_reset)
				{
					element_list[depth_to_reset] = 0;
				}
			}
		}
	}


	if( partition[0] == this->network_size )
		this->outputted_last_partition = true;

	return partition;
}


// -----------------------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------------------
void remapNodes(unsigned int* z, const unsigned int* __restrict__ inp, const unsigned int* __restrict__ node_map, const unsigned int map_size)
{

	for(unsigned int part=0; part<inp[0]; ++part)	{

		// OPTIMIZATION: node will always be >= part, so we start at part.
		// instead of starting at node=0
		for(unsigned int node=part; node<map_size; ++node)	{

			if(inp[part+1] & (1 << node)) {
				z[part] |= (1 << node_map[node]);
//				blah = max(part,node);
//				assert( part <= node );
			}
		}
	}
}


#endif
