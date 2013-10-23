#!/usr/bin/python
############################################################
# This program generates the todo/ directory for performing a Hadoop MapReduce job.
# Arguments:
#       networkfilename
#       num_workers
#       maxnumparts
############################################################



NETWORK_FILENAME = '/consciousness/circuits/circuit.txt'
PARTITIONS_GENERATOR = '/consciousness/gen_partitions'
NUM_WORKERS = 10

# Set this to a high value to always attempt all partitions
MAXNUMPARTS = 99999

OUTPUT_DIRECTORY = '/consciousness/mapreduce/todo'

from os import mkdir
from os import system
import os.path

import sys

# Make the directory if it doesn't exist.
if os.path.isdir( OUTPUT_DIRECTORY ):
    system('rm -rf %s' % OUTPUT_DIRECTORY )
else:
    mkdir( OUTPUT_DIRECTORY )
    
for worker_index in range(NUM_WORKERS):
    filename = 'job-%s-of-%s' % (worker_index, NUM_WORKERS)
    output_filename = "%s/%s" % ( OUTPUT_DIRECTORY, filename )
    
    f = open( output_filename, 'w' )
    string_to_write = "%(worker_index)s %(NUM_WORKERS)s %(NETWORK_FILENAME)s %(MAXNUMPARTS)s  %(PARTITIONS_GENERATOR)s\n" % locals()    
    f.write( string_to_write )
    
    f.close()
    print "Finished making %s" % output_filename

print "Done!"
