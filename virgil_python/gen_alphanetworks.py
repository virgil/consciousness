#!/usr/bin/python
############################################################
# This program generates random neural networks.
# USAGE:  ./gen_randnetworks.py output_dir [#to_make] [minN-maxN] [avr_#connections]
############################################################

import philib_helpers
import os, sys, getopt, string
from random import randint, uniform, choice, sample
from numpy import zeros, ones, arange
from TerminalController import TerminalController

######################################################################################################
# These must be specified as command-line parameters

# List of the allowed network sizes, thresholds, and fan-in
N = 12
Ks = [2,3,4]
ALPHA_inc, THETA_inc = 0.05, 0.5
ALPHAs = list( arange(0.0,1.0+ALPHA_inc,ALPHA_inc) )	
THETAs = list( arange(1.0,7.0+THETA_inc,THETA_inc) )
INSIDE_WEIGHT = 1.0

# if this is true, use random 10-letter a-z0-9.txt filenames
OUTPUT_DIRECTORY = None
PRINT_OUTPUT_FILENAMES = True
######################################################################################################
term = TerminalController()

def group( node_index, groups ):
	'''returns the group index of the node'''
	
	group_index = [ group_index for group_index, group in enumerate(groups) if node_index in group ]

	assert len(group_index) == 1, "group index != 1"

	group_index = group_index[0]
	
	return group_index

'''
def parsearg( opt_string ):

	z = []
	if '-' in opt_string:
		pieces = map(int,opt_string.split('-'))
		assert len(pieces) == 2
		pieces.sort()
		z = range(pieces[0], pieces[1]+1)

	elif '...' in opt_string:
		pieces = map(int,opt_string.split('...'))
		assert len(pieces) == 2
		pieces.sort()
		z = range(pieces[0], pieces[1]+1)

	elif ',' in opt_string:
		z = map(float,opt_string.split(','))
		z.sort()
		
	elif opt_string.isdigit():
		z = [float(opt_string)]

	else:
		print "error! don't know argument '%s'" % opt_string
		sys.exit(1)

	return z
'''


def print_usage():
	print term.NORMAL
	print '''\
############################################################
# This program generates random neural networks.
# USAGE:  ./gen_alphanetworks.py [output_dir]

	'''
	
	sys.exit(0)

def main():	
	global OUTPUT_DIRECTORY
	
	if OUTPUT_DIRECTORY is None:
		if len(sys.argv) == 2:
			OUTPUT_DIRECTORY = sys.argv[1].rstrip('/')
		else:
			print_usage()
	
	if not os.path.isdir( OUTPUT_DIRECTORY ):
		print "- creating directory", OUTPUT_DIRECTORY
		os.makedirs( OUTPUT_DIRECTORY )
	
	for K in Ks:
		GROUPS = []
		# define the K groups
		if K == 2:
			GROUPS = [ [0,1,2,3,4,5], [6,7,8,9,10,11] ]
		elif K == 3:
			GROUPS = [ [0,1,2,3], [4,5,6,7], [8,9,10,11] ]
		elif K == 4:
			GROUPS = [ [0,1,2], [3,4,5], [6,7,8], [9,10,11] ]
		else:
			print "K must be 2, 3 or 4"
			sys.exit(1)

		total_made, total_num_to_make = 0, len(ALPHAs) * len(THETAs)
	   	filename_index = 1

		print "total_to_make=%s" % total_num_to_make
		print "groups=%s" % GROUPS

		for ALPHA in ALPHAs:
			for THETA in THETAs:
				comments = { 'ALPHA': ALPHA, 'THETA': THETA, 'K': K }			
				thresholds = [THETA] * N

				print "Writing theta=%(THETA)s \t alpha=%(ALPHA)s" % locals()
				# assume EVERYTHING is in the INSIDE group			
				cxnmatrix = ones( (N, N) ) * INSIDE_WEIGHT

				for from_node in range(N):
					for to_node in range(N):
						# if the FROM and TO node are in DIFFERENT groups, set weight to ALPHA
						if group(from_node, GROUPS) != group(to_node, GROUPS):
							cxnmatrix[from_node][to_node] = ALPHA

				if OUTPUT_DIRECTORY:
					ofilename = OUTPUT_DIRECTORY.rstrip('/')
				ofilename += '/k%s/a%s_t%s.txt' % ( K, ALPHA, THETA )

				# don't regenerate any files we've already made
				if os.path.isfile( ofilename ) or os.path.isfile( ofilename.replace('.txt','.out') ):
					continue
				
				if philib_helpers.INPUT_savetofile( ofilename, cxnmatrix, thresholds, comments ):
					total_made += 1 
					print "\r%%Done=%s%%" % ( int(round((float(total_made) / total_num_to_make) * 100.0)) ),
					if PRINT_OUTPUT_FILENAMES:
						print " \t -> \t '%s'" % ofilename

				else:
					print "Had error writing to file %s" % ofilename

				sys.stdout.flush()

	print "\n"
	print "Generated %s random networks." % total_made
	print "Finished writing to directory '%s'" % OUTPUT_DIRECTORY

if __name__ == '__main__':
	main()

