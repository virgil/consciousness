#ifndef _PARTITION_ENUMERATOR_H_
#define _PARTITION_ENUMERATOR_H_

//#include <vector>

class PartitionEnumerator
{
private:
	unsigned short	network_size;				// Number of elements you want to partition
	unsigned short	untouched_nodes;			// Number of nodes that don't need processing
		
	unsigned int* __restrict__ partition;					// Partition, which is a list of elements. First entry is the partition size
	unsigned int* __restrict__ partition_sizeList;			// Partition size before adding a new node
	unsigned int* __restrict__ element_list;				// List of elements that were chosen to make up the current partition
	unsigned int* __restrict__ element_countList;			// Number of elements at a given depth
	unsigned int* __restrict__ node_posList;				// Used to remove the node from the right partition element
	unsigned int* __restrict__ partition_increaseList;		// Used to recover the partition size (result[0])
	

	// used by nextPartition to know when it's reach the end.
	bool outputted_last_partition;
	
	int temp;
	unsigned short max_num_parts;
	
public:
	PartitionEnumerator(unsigned short initial_network_size, bool inp_skip_total_partition=false );	
	virtual ~PartitionEnumerator();

	unsigned int size();
	void reset(unsigned short new_network_size, unsigned short maxparts=0 );
	
	void resize(unsigned short new_network_size);	
	
	// ---------------------------------------------------------------------------------------------------------------------------
	// The following functions return a pointer to a valid partition. The format is as follows:
	//	- 1st entry:			# of elements (E) in the partition.
	//	- 'E' next entries:		Encoding of the 'E' partition elements (node 'n' is in a partition iff n-th bit is set to 'true').  
	//	- Remaining entries:	Unused garbage space.
	// ---------------------------------------------------------------------------------------------------------------------------
	const unsigned int*	firstPartition();
	const unsigned int*	currentPartition();
	const unsigned int*	nextPartition();

	bool skip_total_partition;
	
};


////////////////////////////////////////////////////////////////////////////////////
// helper functions
////////////////////////////////////////////////////////////////////////////////////
void remapNodes(unsigned int* z, const unsigned int* __restrict__ inp, const unsigned int* __restrict__ node_map, unsigned int map_size);

#endif

