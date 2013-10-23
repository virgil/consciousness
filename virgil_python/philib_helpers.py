#!/usr/bin/python

import numpy
from math import sqrt
from numpy import array
import os.path

def INPUT_savetofile( ofilename, cxnmatrix, thresholds, comments=None, remarks=None ):
	'''returns True is the save was successful.	 comments = comments:
	remarks = comments that begin with a #'''
	
	dirname = os.path.dirname( ofilename )
	if not os.path.isdir( dirname ):
		os.makedirs( dirname )

	###################################################################################################
	# Process cxnmatrix and thresholds
	###################################################################################################
	matrix_string = '\n'.join([ ' '.join([ "%+.3f" % ele for ele in row ]) for row in cxnmatrix ])
	threshold_string = ' '.join([ "%+.3f" % x for x in thresholds])

	numnodes = len(thresholds)
	assert len(thresholds) > 0, "* Error: # thresholds == 0"
	assert len(cxnmatrix) == len(thresholds), "* Error: matrix not same length as thresholds"

	for row in cxnmatrix:
		assert len(cxnmatrix) == len(row), "#rows=%s but #cols=%s" % ( len(cxnmatrix), len(row) )

	###################################################################################################
	# Process comments
	###################################################################################################
	comments_string = None
	if type(comments) is str:
		if not commments.startswith('comment: '):
			comments_string = 'comment: '
			comment_string += comments
	elif type(comments) is list:
		comment_string = '\n'.join([ 'comment: %s' % x for x in comments ])
	elif type(comments) is dict:
			comments_string = '\n'.join([ 'comment: %s=%s' % (arg.lower(),value) for arg, value in comments.iteritems() ] )

	###################################################################################################
	# Process remarks
	###################################################################################################

	remarks_string = ''
	if remarks:
		assert type(remarks) is list, "remarks must be a list!!"
		remarks = [ x for x in remarks if x.startswith('#') ]
		remarks.extend( [ '# %s' % x for x in remarks if not x.startswith('#') ] )		
		remarks_string = '\n'.join( remarks )

		# prepend the remarks to the comments
		comments_string = remarks_string + '\n' + comments_string
	
	comments_string = comments_string.replace('\n\n','\n')

	output_string = \
'''%(comments_string)s

type=neural

nodes=%(numnodes)s
thresholds=%(threshold_string)s

%(matrix_string)s
''' % locals()

	#print output_string
	f = open(ofilename, 'w')
	f.write( output_string )
	f.close()
	
	return os.path.isfile( ofilename )

	
	
def std_err( values ):
	'''returns the standard errors of the passed values'''
	if len(values) <= 1:
		return numpy.inf
	
	return numpy.std(values) / sqrt(len(values)-1) 

def is_contiguous( cxnmatrix ):
	'''This function returns True if the cxnmatrix is contiguous.  False otherwise.'''
	
	cxnmatrix = array( cxnmatrix )
	
	assert cxnmatrix.shape[0] == cxnmatrix.shape[1], "cxnmatrix must be square"
	
	N = cxnmatrix.shape[0]
	# null case of the matrix only containing a single element
	if N == 1:
		return True

	assert( N >= 2 )
	unvisited_nodes, unseen_nodes = range(N), range( N )

	# now recursively drill-down
	return drill( [0], cxnmatrix, unvisited_nodes, unseen_nodes )


def subset_distance( s1, s2 ):
	'''returns the distance between [0.0,1.0] between subsets s1 and s2'''

#	print "\t s1=%s" % s1
#	print "\t s2=%s" % s2

	s1, s2 = set(s1), set(s2)
	both = s1.union(s2)
	
#	print "\t set(s1)=%s" % s1
#	print "\t set(s2)=%s" % s2

	return (len(s2-s1) + len(s1-s2)) / float(len(both))
	
def subsets_distance( ss1, ss2, norm=False ):
	'''returns the distance between two lists of subsets ss1 and ss2.
	if norm=True returns between [0.0,1.0], else is between [0.0, max(len(ss1),len(ss2))]'''

	if not ss1 and not ss2:
		return None

	# swap if ss2 is bigger
	if len(ss2) > len(ss1):
		ss1, ss2 = ss2, ss1

	if not ss2:
		return len(ss1)

#	min_distances = []
	# Sorting might help when have multiple s2's with the same min_distance
	# If nothing else, sorting never hurts
	ss1.sort()
	ss2.sort()
	
	distance = 0.0
	assert distance >= 0.0, "* Error -- invalid distance"
	
	for s1 in ss1:
#		print "ss1=%s" % ss1
#		print "s1=%s" % s1
#		print "ss2=%s" % ss2
		
		min_dist = min( [ subset_distance(s1,s2) for s2 in ss2 ] )
#		min_distances.append( min_dist )
		distance += min_dist

	assert 0 <= distance <= len(ss1), "Error -- got an invalid distance"

	if norm:
		distance /= float(len(ss1))
	
	return distance

	
def drill( current_locations, cxnmatrix, unvisited_nodes, unseen_nodes ):
	'''Checks to see if the network is contiguous yet'''

	# add the current locations to the processed_nodes
	unvisited_nodes = [ x for x in unvisited_nodes if x not in current_locations ]

	if unseen_nodes == []:
		return True

	
	# Get all nodes current_locations connects to
	connected_nodes = []	
	origins, dests = numpy.where(cxnmatrix != 0)	
	for (origin_node, dest_node) in zip(origins, dests):

		if origin_node in current_locations:
			connected_nodes.append( dest_node )
		if dest_node in current_locations:
			connected_nodes.append( origin_node )

	#uniquify the list
	connected_nodes = list(set(connected_nodes))

	# remove all connected_nodes from unseen_nodes
	for connected_node in connected_nodes:
		if connected_node in unseen_nodes:	  
			unseen_nodes.remove(connected_node)

	nodes_to_visit = [ n for n in connected_nodes if n in unvisited_nodes ]

	# if there are unseen_nodes and we have no where else to go, the network is disconnected.  Return FALSE
	if not nodes_to_visit and unseen_nodes != []:
		return False

	return drill( nodes_to_visit, cxnmatrix, unvisited_nodes[:], unseen_nodes[:] )


if __name__ == '__main__':
	'''Testing the is_contiguous function...'''

	As = []
	As.append( array([[1,0],[0,0]]) )
	As.append( array([[1,0,1],[1,0,0],[0,0,1]]) )
	As.append( array([[0,1,0],[0,0,0],[0,0,1]]) )
	As.append( array([[0,0,1,0],[0,0,1,0],[0,0,0,0],[0,0,0,1]]) )
	As.append( array([[0,1,1,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]]) )
	print "----------------------------------------"		

	for A in As:
		print A
		origins, dests = numpy.where(A != 0)
		for o, d in zip(origins,dests):
			print "%s -> %s" % (o,d)

		print "is_contiguous()=%s" % is_contiguous(A)
		print "----------------------------------------"