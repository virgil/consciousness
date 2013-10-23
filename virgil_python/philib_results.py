#!/usr/bin/python
#########################################################################################################
# This class reads in a .out file and returns a happy data-structure representing its contents.
#########################################################################################################
from pprint import pprint
import os, os.path, string
from numpy import array, flatnonzero
from TerminalController import TerminalController
term = TerminalController()
import sys

MAXIMUM_FILESIZE_TO_PARSE = 2101832
DELETE_BIG_FILES = True

def PHI_Inputs( dirs, print_errors=True, get_done=False ):

#	 print "I am here!!!! 1"
#	 sys.stdout.flush()


	if type(dirs) is str and os.path.isfile(dirs) and not os.path.isdir(dirs):
		full_filename = dirs

		print "doing filename=%s" % full_filename
				
		if not full_filename.endswith('.txt') and not get_done:
			yield (None, None)

		try:
			result = PHI_Input( full_filename )
		except AssertionError, e: 
			if print_errors:
				print "error=", e
				print "* Assert Error: '%(full_filename)s'.	 Skipping." % locals()
		except ValueError:
			if print_errors:
				print 'error=', e
				print "* Value Error: '%(full_filename)s'.	Skipping." % locals()
		else:
			yield (result, full_filename)
	
	else:
#		 print "I am here!!!!"
		if type(dirs) is str:
			dirs = [dirs]

		# What to do if passed a list of directories...
		dirs = [ d.rstrip('/ ') for d in dirs if os.path.isdir(d) ]
#		 print "dirs=%s" % dirs
		
		for d in dirs:
 #			 print "\t d=%s" % d
			# get the .txt's in the directory, and sort the list
			if get_done:
			    filenames = [ x for x in os.listdir(d) if x.endswith('.txt') or x.endswith('.done') ]
			else:
			    filenames = [ x for x in os.listdir(d) if x.endswith('.txt') ]

			filenames.sort()
			
			for filename in filenames:

				full_filename = d + '/' + filename
			
				try:
					result = PHI_Input( full_filename )
				except AssertionError, e: 
					if print_errors:
						print "error=", e
						print "* Assert Error: '%(full_filename)s'.	 Skipping." % locals()
				except ValueError,e :
					if print_errors:
						print "error=", e
						print "* Value Error: '%(full_filename)s'.	Skipping." % locals()
				except UnknownTag:
					if print_errors:
						print "error=", e						
						print "* Unknown Tag: %s.  Delete?" % full_filename
#						os.remove( full_filename )

				else:
					yield (result, full_filename)



def PHI_Results( dirs, print_errors=True ):
	global MAXIMUM_FILESIZE_TO_PARSE, DELETE_BIG_FILES
#	 print "dirs=%s" % dirs
	
	# What to do if this is a filename...	 
	if type(dirs) is str and os.path.isfile(dirs) and not os.path.isdir(dirs):
		full_filename = dirs

		print "doing filename=%s" % full_filename
				
		if not full_filename.endswith('.out'):
			yield (None, None)

		try:
			result = PHI_Output( full_filename )
			
		except FileNotDone:
#			if print_errors:
			print "* File '%(full_filename)s' not done.	 Deleting it." % locals()			
#			os.remove( full_filename )

		except AssertionError, e: 
			if print_errors:
				print 'error=', e
				print "* Assert Error: '%(full_filename)s'.	 Skipping." % locals()
		except ValueError, e:
			if print_errors:
				print 'error=', e
				print "* Value Error: '%(full_filename)s'.	Skipping." % locals()
		except UnknownTag, e:
			if print_errors:
				print 'error=', e
				print "* Unknown Tag: %s." % full_filename
#				os.remove( full_filename )

		else:
			
			yield (result, full_filename)
	
	else:		 
		if type(dirs) is str:
			dirs = [dirs]

		# What to do if passed a list of directories...
		dirs = [ d.rstrip('/ ') for d in dirs if os.path.isdir(d) ]
#		 print "dirs=%s" % dirs
		
		for d in dirs:
 #			 print "\t d=%s" % d
			# get the .out's in the directory, and sort the list
			filenames = [ x for x in os.listdir(d) if x.endswith('.out') ]
			filenames.sort()			
#			print "filenames=%s" % filenames

			for filename in filenames:
				full_filename = d + '/' + filename

				if os.path.getsize(full_filename) > MAXIMUM_FILESIZE_TO_PARSE:
					print "Skipping file=%s for being too big" % full_filename
					if DELETE_BIG_FILES:
						print "Deleting %s..." % full_filename
						os.remove( full_filename )
					sys.stdout.flush()
					continue

	#			print "full_filename=%s" % full_filename
			
				try:
#					 print "processing %s..." % full_filename
#					 sys.stdout.flush()
					result = PHI_Output( full_filename )
				except FileNotDone:
#					if print_errors:

					print "* File '%(full_filename)s' not done.	 Delete it?" % locals()
#					os.remove( full_filename )
				except AssertionError, e: 
					if print_errors:
						print 'error=', e
						print "* Assert Error: '%(full_filename)s'.	 Skipping." % locals()
				except ValueError, e:
					if print_errors:
						print 'error=', e
						print "* Value Error: '%(full_filename)s'.	Skipping." % locals()
				except CorruptTag, e:
					if print_errors:
						print 'error=', e						
						print "* Unknown Tag: %s.  Delete?" % full_filename
#						os.remove( full_filename )

				except UnknownTag, e:
					if print_errors:
						print 'error=', e						
						print "* Unknown Tag: %s.  Delete?" % full_filename
#						os.remove( full_filename )
						
				else:
					yield (result, full_filename)

			

def strip_equals( x ):
	'''strip anything before the first equals sign of variable x'''
	if type(x) is str and '=' in x:
		x = x[x.index('=')+1:]		
	return x

class FileNotDone(Exception):
	def __init__(self, value):
		self.parameter = value
	def __str__(self):
		return repr(self.parameter)

class UnknownTag(Exception):
	def __init__(self, value):
		self.parameter = value
	def __str__(self):
		return repr(self.parameter)

class CorruptTag(Exception):
	def __init__(self, value):
		self.parameter = value
	def __str__(self):
		return repr(self.parameter)

class FileExists(Exception):
	def __init__(self, value):
		self.parameter = value
	def __str__(self):
		return repr(self.parameter)


class PHI_Input:
	'''stores all of the data from a .txt file'''

	def __str__(self):
		'''return the string representation of PHI_Input'''
		z = """\
\t type=%s \t numunits=%s
\t thresholds=%s \t indegrees=%s
\t cxnmatrix:
%s
""" % ( self.type, self.numunits, self.thresholds, self.indegrees, self.cxnmatrix )
		
		if self.rawcomments:
			z += "\t rawcomments=%s" % (self.rawcomments)
		if self.comments:
			z += "\n\t comments=%s" % (self.comments)
		
		return z
			
	def __init__(self, inputfilename):

		if not os.path.isfile( inputfilename ):
			raise ValueError, "inputfilename='%s' is not a valid file" % inputfilename
			
		self.rawcomments, self.thresholds, self.indegrees, self.cxnmatrix = [], [], [], []
		self.numunits, self.type, self.filename = None, None, None
		self.comments = {}

		tempdigits = string.digits + '+-'
		
		f = open( inputfilename )
		for line in f:
			
			line = line.rstrip()
			
			if not line:
				continue

			if line.startswith('comment:'):
				# add the raw comment
				line = line[len('comment:')+1:]
				self.rawcomments.append( line )

				# try to parse the comment
				if '=' in line:
					pieces = line.split('=')
					key, value = pieces[0], pieces[1]

					assert len(pieces) == 2, "pieces in comment was >2 !"
					assert key not in self.comments, "key='%s' was already in comments"

					try:
						value = float(value)
					except:
						value = pieces[1]

					self.comments[key] = value
				
							
			elif line.startswith('type='):
				self.type = strip_equals(line)
				
			elif line.startswith('nodes='):
				self.numunits = int(strip_equals(line))
				
			elif line.startswith('thresholds='):
				self.thresholds = map(float, [ x for x in strip_equals(line).split(' ') ] )
				
			elif line[0] in tempdigits and ( len(line.split(' ')) == self.numunits ) :
				row = map(float,line.split(' '))
				self.cxnmatrix.append( row )
			
		f.close()
		
		self.cxnmatrix = array(self.cxnmatrix)
		
#		print ".T=%s" % self.cxnmatrix.T
		self.indegrees = [ len(flatnonzero( col )) for col in (self.cxnmatrix.T) ] 

		self.cleanup();


	def write( self, outputfilename=None, overwrite=False ):
		''' write the current network to outputfilename'''
		
		f = None
		
		# Delete the outputfilename if it exists
		if not outputfilename is None and os.path.isfile( outputfilename ):
			if overwrite:
				os.path.remove( outputfilename )
			else:
				print "* Error: filename='%s' already exists!  Can pass parameter overwrite=True to overwrite" % outputfilename
				raise FileExists
		
		if outputfilename is not None:
			f = open( outputfilename, 'w' )
			assert f, "could not open filename %s!" % outputfilename		

		# write the comments
		if self.rawcomments:
			for rawcomment in self.rawcomments:
				# how to handle the thetas...
				if 'thetas=' in rawcomment:
					f.write("comment: thetas=%s\n" % ' '.join( map(str,list(set(self.thresholds))) ) )
				else:
					f.write( "comment: %s\n" % rawcomment )
		
		# write the type
		assert self.type in ['neural','circuit'], "type='%s' is not neural or circuit!" % self.type
		type_string = self.type
		
		# write the nodes
		nodes_string = self.numunits
		
		# write the thresholds
		threshold_string = ' '.join([ "%+.3f" % x for x in self.thresholds])
#		 f.write("thresholds=%s\n" % threshold_string )
	
		matrix_string = '\n'.join([ ' '.join([ "%+.3f" % ele for ele in row ]) for row in self.cxnmatrix ])	   
		
		f.write( '''
type=%(type_string)s
nodes=%(nodes_string)s
thresholds=%(threshold_string)s
		
%(matrix_string)s
''' % locals() )
		

		f.close()

		return True
		
	def cleanup(self):
		'''do some asserts, set blanks to None'''

		assert self.type is None or self.type.upper() == 'NEURAL' or self.type.upper() == 'CIRCUIT', "invalid type.	type='%s'" % self.type
		assert self.numunits > 0, "numunits must be >0"



		assert len(self.indegrees) == self.numunits, "len(indegrees list) != numunits"
		assert len(self.thresholds) == self.numunits, "#thresholds != numunits"
		assert min(self.thresholds) >= 0, "Thresholds '%s' had a value <0" % self.thresholds
		
		if not self.comments:
			self.comments = None
		if not self.rawcomments:
			self.rawcomments = None

		assert self.cxnmatrix.shape[0] == self.cxnmatrix.shape[1], "cxnmatrix isn't square"
		assert self.cxnmatrix.shape[0] == self.numunits, "cxnmatrix isn't the size of numunits"		   



class PHI_Output:
	'''stores all of the data from a .out file'''
	def __str__(self):
		z = '''
stats: %s
comments: %s
bracketphi: %s
highestbracketphis: %s 
highest_bracket_ei: %s
<MCs>: %s\
''' % ( self.stats, self.comments, self.bracketphi, self.highestbracketphis, self.highest_bracket_ei, self.bracket_maincomplexes )

		if self.x1states:
			z += '\n x1states:', self.x1states
		if self.AVERAGED:
			z += '\n AVERAGED: %s' % str(self.AVERAGED)

		if self.highest_state_eis:
			z += '\n also: highest_state_eis'
		if self.x1s:
			z += '\n also: x1s[]'
		if self.highestx1s:
			z += '\n also: highestx1s'

		return z


	def __init__(self, inputfilename ):
		# strip baggage from the beginning/end
		inputfilename = inputfilename.strip()
		
		# require that the inputfilename exists and is done.
		self.AssertDone( inputfilename )
#		print "inputfilename='%s'" % inputfilename


		entries = []

		# initialize everything to None
		self.stats, self.x1states, self.AVERAGED = None, [], None
		self.x1s, self.highestx1s = {}, {}
		self.bracketphi, self.highestbracketphis = None, None
		self.highest_state_eis, self.highest_bracket_ei = None, None
		self.comments = None
		self.bracket_maincomplexes = None

		f = open( inputfilename )
#		print term.RED + "ENTRIES=" + term.NORMAL,
		for line in f:
			# remove anything with # and spaces
			line = line.lstrip('# ')
			line = line.rstrip()
			lineparts = line.split('\t')
			entry = lineparts.pop(0).upper()

			if not entry:
				continue
			
			entries.append( entry )

			
			# convert the lineparts to a dictionary
			contents = {}
			for linepart in lineparts:
				key, arg = linepart.split('=')
				contents[key] = arg			
			
#			print term.WHITE + "%s" % entry.rstrip(':') + term.NORMAL,
#			print "contents=%s" % contents
			
			
			if entry == 'STATISTICS:':
				assert self.stats is None
				for tag in 'N','max_num_parts','#x1states':
					assert tag in contents, "didnt find required tag='%s'" % tag
				
				self.stats = statistics( contents, inputfilename )

			elif entry == 'STATE:':
				for tag in 'x1','p(x1)','ei','#MIPs','PHIs','MIPs':
					assert tag in contents, "didnt find required tag='%s'" % tag

				x1 = int(contents['x1'])
				assert x1 not in self.x1s, "x1 was not in self.x1s!"

				self.x1s[x1] = PHI_result( contents, self.stats )
#				print self.x1s[x1]
#				raw_input('...')

			# Content on AVERAGED: goes into a dictionary called 'AVERAGED'
			elif entry == 'STATE-AVERAGED:':
				assert self.AVERAGED is None
				for tag in 'ei','phi':
					assert tag in contents, "didnt find required tag='%s'" % tag
				
				self.AVERAGED = {}
				for key, arg in contents.iteritems():
					self.AVERAGED[key] = float(arg)
			elif entry == 'COMMENT:':
				if self.comments is None:
					self.comments = {}
				for key, value in contents.iteritems():
					assert key not in self.comments, "key=%s was already in self.comments!" % key
					try:
						value = float(value)
					except:
						value = str(value)

					finally:
						self.comments[key] = value
				
#				 print self.comments
				
			elif entry == 'HIGHESTSTATEPHI:':
				for tag in 'x1','subset','ei','#MIPs','PHIs','MIPs':
					assert tag in contents, "didnt find required tag='%s'" % tag

				x1 = int(contents['x1'])
				# if this x1 is new (99% of the time), initialize the list to []
				self.highestx1s.setdefault(x1,[])				
				self.highestx1s[x1].append( PHI_result(contents, self.stats) )
#				pprint( self.highestx1s )
								
			elif entry == 'BRACKET:':
				assert self.bracketphi is None, "More than 1 BRACKETPHI:!"
				for tag in 'ei','#MIPs','PHIs','MIPs':
					assert tag in contents, "didnt find required tag='%s'" % tag
									
				self.bracketphi = PHI_result( contents, self.stats )
				
			elif entry == 'HIGHESTBRACKETPHI:':
				for tag in 'subset','ei','#MIPs','PHIs','MIPs':
					assert tag in contents, "didnt find required tag='%s'" % tag

				if self.highestbracketphis is None:
					self.highestbracketphis = []
				
				self.highestbracketphis.append( PHI_result(contents, self.stats) )
				
				# the main complex must have a phi >= 0
				if min(self.highestbracketphis[-1].PHIs) <= 0.0:
					del self.highestbracketphis[-1]
								
#				 print "res=%s" % self.highestbracketphis[-1]
				
			elif entry == 'HIGHEST_BRACKET_EI:':
				assert self.highest_bracket_ei is None, "Found two highest <ei>'s! "
				for tag in 'ei', 'subsets':
					assert tag in contents, "didnt find required tag='%s'" % tag				
				
				self.highest_bracket_ei = EI_result( contents )
				

			elif entry == 'HIGHEST_STATE_EI:':
				for tag in 'x1', 'ei', 'subsets':
					assert tag in contents, "didnt find required tag='%s'" % tag
				
				# If this is the first highest_state_ei, set it {}
				if self.highest_state_eis is None:
					self.highest_state_eis = {}
				
				r = EI_result( contents )
				x1 = r.x1
				assert x1 not in self.highest_state_eis, "state x1=%s was already in highest state eis" % x1
				self.highest_state_eis[x1] = r
			
			elif entry == 'BRACKETMAINCOMPLEX:':
				for tag in 'subset','ei','#MIPs','PHIs','MIPs','#s1states':
					assert tag in contents, "didnt find required tag='%s'" % tag

				if self.bracket_maincomplexes is None:
					self.bracket_maincomplexes = []
				
				r = PHI_result(contents, self.stats)
				if min( r.PHIs ) <= 0.0:
					print "* Error: Had a zero phi for a main complex.  Exiting"
					sys.exit(1)
				
				self.bracket_maincomplexes.append( r )
				
			elif entry == 'DONE!':
				break
			
			else:
#				 print "entry='%s'" % entry
#				 print "filename='%s'" % inputfilename
				raise UnknownTag, "Did not know command '%s'" % entry
				


#		entries_hist = [ (x, entries.count(x)) for x in set(entries) if x != 'DONE!' ]
#		entries_hist.sort()

#		print ""
		
#		for entry, num in entries_hist:
#			entry = str(entry).ljust(30)
#			print "%s" % entry + "%s".rjust(5) % num

		# if this is over 20, delete it.
#		if len(self.highestbracketphis) >= 20:
#			print "filename=%s had over 20 maincomplexes.  Deleting." % ( inputfilename )
#			os.remove( inputfilename )
#			raise FileNotDone, "File wasn't done!"
				
		f.close()
		
		
		
	def AssertDone(self, inputfilename ):
		'''verifies that the final line of the inputfilename is "DONE!" '''

		if not os.path.isfile( inputfilename ):
			raise ValueError, "inputfilename '%s' is not a valid file" % inputfilename
		
		f = open( inputfilename )
		lastline = (f.readlines()[-1]).strip().upper()		
		f.close()
		
		if lastline != 'DONE!':
			raise FileNotDone, "'%s' was not done." % inputfilename


class statistics:
	''' stores the first line of the .out file'''
	def __init__(self, inputs, filename):
		
#		pprint( inputs )		
		try:
			self.numunits = int( inputs['N'] )
			self.max_num_parts = int( inputs['max_num_parts'] )
			self.num_x1_states = int( inputs['#x1states'] )
			self.H_X1 = float( inputs['H(X1)'] )
		except:
			raise ValueError, "Could not convert arguments to integers"

		# do the normalization
		try:
			self.normalization = str( inputs['normalization'] ).upper()
		except KeyError:
			self.normalization = 'TONONI'
		
		try:
			self.EI = str( inputs['EI'] ).upper()
		except KeyError:
			self.EI = None
		
		assert self.H_X1 >= 0.0, "H(X1) must be >=0.0"
		assert self.numunits > 0, "numunits must be >0"
		assert self.max_num_parts > 0, "max_num_parts must be >0"
		assert self.num_x1_states > 0, "num_x1_states must be >0"
		
		assert self.max_num_parts <= self.numunits, "Error.	 This is impossible. #maxparts must be <= #units "
		
		assert self.normalization in ['TONONI','NONE','ATOMIC','ADAMI','UNITY_SCORE_SUM','SYNERGY_PLUS_ANTERGY','HOLISM','UNITY_SCORE_PRODUCT','VIRGIL_INFO_RATIO','EDLUND_INFO_RATIO','TONONI_PROD_M0','KOCH_M0','KOCH_M1'], "Normalization '%s' wasn't recognized" % self.normalization
		assert self.EI in ['ANTI-CONDITIONED','CONDITIONED','NAKED','NEW_HOLISM','NEW_SYNERGY','TONONI_ANTI-CONDITIONED'], "EI '%s' was not recognized." % self.perturb

		self.filename = strip_equals(filename)

		# aliases
		self.N = self.numunits
		self.norm = self.normalization

	def __str__(self):
		return "file=%s \t N=%s \t K=%s \t norm=%s" % ( self.filename, self.N, self.max_num_parts, self.norm )


class EI_result:
	'''Data structure to hold a x1state, subsets, and ei.'''
	def __str__(self):
		if self.x1 is None:
			return "ei=%s \t subsets=%s" % (self.ei, self.subsets )
		else:
			return "x1=%s ei=%s \t subsets=%s" % (self.x1, self.ei, self.subsets )

	def __init__( self, inputs ):
		assert type(inputs) is dict
		self.ei, self.x1, self.subsets = None, None, None
		
		# do ei
		self.ei = float( inputs['ei'] )
		assert self.ei >= 0.0, "* ei_result: Invalid ei"

		# do subsets
		self.subsets = [ map(int,x.split(' ')) for x in (inputs['subsets']).split(';') ]

		# sort each entry
		for i in range(len(self.subsets)):
			self.subsets[i].sort()
			
		# sort the whole list
		self.subsets.sort()
		
		# do x1, if we're doing that
		if 'x1' in inputs:
			self.x1 = int(inputs['x1'])
			assert self.x1 >= 0, "* ei_result: Invalid x1 state"

class PHI_result:
	'''Data-structure to store the result of a phi for a particular state or bracketphi'''
	def __str__(self):
		z = """\
		ei=%s \t subset=%s
		PHIs=%s \t MIPs=%s
		""" % (self.ei, self.subset, self.PHIs, self.MIPs )
		
		if self.num_s1states:
			z += "\t num_s1states=%s" % self.num_s1states

		if self.x1 and self.prob:
			z += "\t x1=%s \t prob(x1)=%s" % (self.x1, self.prob)
		
		return z
		
		
	def __init__( self, inputs, stats ):
		assert type(inputs) is dict

		self.ei = float( inputs['ei'] )
#		self.phi_norm = float( inputs['phi_norm'] )
		self.PHIs = map(float, inputs['PHIs'].split(';') )
		partitions = [ x.split('/') for x in inputs['MIPs'].split(';') ]
		
		for i, partition in enumerate(partitions[:]):
			for j, part in enumerate(partition):
				if ' ' in part:
					part = map(int, part.split(' ') )
				else:
					part = [int(part)]
				
				part.sort()				
				partitions[i][j] = part
		
		self.MIPs = partitions
		################################################################################################
		# derive the subset from the MIPs
		################################################################################################		
		self.subset = []
		for part in self.MIPs[0]:
			self.subset.extend( part )
		self.subset.sort()

		################################################################################################
		# Get the HOLISMs if they were defined
		################################################################################################		
		self.holisms = None
		if 'AGGREGATE_HOLISMs' in inputs:
			self.holisms = [ float(x) for x in inputs['AGGREGATE_HOLISMs'].split(';') ]
			assert len(self.holisms)
			
			
		self.total_holisms = None
		if 'AGGREGATE_HOLISMs' in inputs:
			self.total_holisms = [ float(x) for x in inputs['AGGREGATE_HOLISMs'].split(';') ]
			assert len(self.total_holisms), 'total_holisms was empty!'

		################################################################################################
		# Get the SYNERGYs and ANTERGYs if they were defined
		################################################################################################		
		self.integrations, self.antergys = None, None
		if 'AGGREGATE_COOPs' in inputs:
			self.integrations = [ float(x) for x in inputs['AGGREGATE_COOPs'].split(';') ]
		if 'ANTERGYs' in inputs:
			self.antergys = [ float(x) for x in inputs['ANTERGYs'].split(';') ]

		self.total_integrations = None
		if 'AGGREGATE_TOTAL_COOPs' in inputs:
			self.total_integrations = [ float(x) for x in inputs['AGGREGATE_TOTAL_COOPs'].split(';') ]
				
		################################################################################################
		# do the x1 and prob_x1 if those were specified
		################################################################################################		
		self.x1, self.prob = None, None
		if 'x1' in inputs:
			self.x1 = float( inputs['x1'] )
		if 'p(x1)' in inputs:
			self.prob = float( inputs['p(x1)'] )
			
		################################################################################################
		# do the #s1states if that's present
		################################################################################################		
		self.num_s1states = None
		if '#s1states' in inputs:
			self.num_s1states = int( inputs['#s1states'] )

		################################################################################################
		# sanity checks
		################################################################################################
		assert len(self.PHIs) == len(self.MIPs)
		
#		for phi in self.PHIs:
#			print "0.0 <= %s <= %s <= %s" % (self.phi_norm, phi, self.ei)
#			assert 0.0 <= self.phi_norm <= phi+.003 or stats.norm == 'ADAMI', "0.0 </= phi_norm </= phi"
#			assert phi-.003 <= self.ei, "phi </= self.ei"
		
		# if x1 and prob are specified, check their bounds
		if self.x1 is not None:
			assert 0 <= self.x1
			# Do not require that prob_x1 exists.  Because on HIGHEST_STATE_PHI there will be no prob_x1
		if self.prob is not None:
			assert 0.0 < self.prob <= 1.0
			assert self.x1 is not None, "prob was specified but not x1"
		
		# if there's an x1 and no subset, make sure the probability exists
		if 'x1' in inputs and 'subset' not in inputs:
			assert 'p(x1)' in inputs
		
		# if subset was specified, check it against the derived subset
		if 'subset' in inputs:
			# the if x takes care of any ''
			subset_check = [ int(x) for x in inputs['subset'].split(' ') if x ]
			subset_check.sort()
			assert self.subset == subset_check, "subset doesn't equal subset_check!"
		################################################################################################				
		

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print "You must specify a .out filename."
		sys.exit(1)

	filenames = list(set([ x for x in sys.argv[1:] if os.path.isfile(x) ]))
	print "filenames=%s" % filenames
	
	units_ext = lambda r: r.stats.numunits
	
	for filename in filenames:
		output = PHI_Output( filename )
		print "numunits=%s" % units_ext( output )
	