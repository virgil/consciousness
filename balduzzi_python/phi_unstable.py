#!/usr/bin/python
# phi.py
# functions for computing integrated information

# The heart of the library is JEI, which computes EI. It is vectorized, so probably difficult to read.

# OPTIMIZATION.
# The obvious function to optimize is PHI. As it stands, it calls JEI every time, and recomputes the
# same basic probabilities over and over againpy. This could be short-circuited by pulling some of the
# work out of JEI into PHI. This would require keeping of track of ORDERING informationpy.

# import fancy stuff
import TerminalController
import numpy as npy
import simplejson as json

# import from normal libraries...
import sys, logging, math, os.path

# now import things into the main namespace
from pprint import pprint, PrettyPrinter
from numpy import array, matrix, ndarray, zeros, log2

# COLORS:     BLACK = BLUE = GREEN = CYAN = RED = MAGENTA = YELLOW = WHITE
term = TerminalController.TerminalController()
BELL_biglist = []
        
'''
logging.basicConfig(level=logging.ERROR,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename='phi_v.log',
                    filemode='w')
'''
pp = PrettyPrinter( indent=4, depth=3, width=10 )


class phi:
    """This class defines all of the needed functions for calculating the PHI of a binary system."""
    
    # This dictionary applies the function to input x.
    # Ex: NODE_FUNCTIONS['and2']([1 2 0 4 5]) = [0 1 0 1 1]
    NODE_FUNCTIONS = {
        # the [0] returns an array instead of a matrix
        'life': lambda x: array( array(x==3,npy.int8)+array(x==13,npy.int8)+array(x==14,npy.int8) > 0 ),
        'parity': lambda x: array(mod(x,2)==1,npy.int8),
        'and2': lambda x: array(x>=2,npy.int8),        
        'threshold2': lambda x: array(x>=2,npy.int8),
        'alwaysfire': lambda x: 1
    }

    # default plots label.  Valid values are 'decimal','binary','hex'
    plotslabel = 'binary'
    
    def __init__(self, N=None, cxnmatrix=None, nodefunc=None):

        ##################################################
        # Are the inputs valid?
        ##################################################
        if N is None and cxnmatrix is not None and self._is_valid_cxnmatrix(cxnmatrix):
            print "N=%s" % N
            N = len(cxnmatrix)

        if N is not None and cxnmatrix is not None and self._is_valid_cxnmatrix(cxnmatrix) and N != len(cxnmatrix):
            raise ValueError, "Passed N=%s is not the correct dimension of the cxnmatrix" % N
            
        if type(N) is not int or N <= 0:
            raise ValueError, 'Number of units must be a positive integer.  You specified=%s' % N

        if cxnmatrix is not None and not self._is_valid_cxnmatrix(cxnmatrix):
            raise ValueError, 'You specified an invalid connection matrix.'

        if nodefunc is not None and nodefunc not in self.NODE_FUNCTIONS and type(nodefunc) is not type(lambda x: x):
#            print "nodefunc=%s" % nodefunc
#            print type(nodefunc)
#            raw_input('...')
            raise ValueError, 'You specified an invalid node functionpy.'
        ##################################################

        # set number of units, number of states
        self.numunits, self.numstates = N, long(2**N)
        
        # no matrix, no nodefunc
        if cxnmatrix is None and nodefunc is None:
            # generate random causal table.
            raise ValueError, "not ready yet!"
            pass

        # no matrix, specified valid nodefunc
        elif cxnmatrix is None and nodefunc is not None and nodefunc in self.NODE_FUNCTIONS:
            # generate random cxngraph, make causal table.
            raise ValueError, "not ready yet!"
            pass

        # specified both, both are valid.
        elif cxnmatrix is not None and nodefunc is not None and self._is_valid_cxnmatrix(cxnmatrix) and (nodefunc in self.NODE_FUNCTIONS or type(nodefunc) is type(lambda x: x)):
            # generate the causal table from the cxnmatrix and nodefunc.

            # If the nodefunc is a lambda function, just use that.
            if type(nodefunc) is type(lambda x: x):
                self.APPLY_GATE = nodefunc
            else:
                self.APPLY_GATE = self.NODE_FUNCTIONS[nodefunc]
            cxnmatrix = self._cm2fastcm(cxnmatrix)
            self.x0_to_x1 = {}


            # we keep a decimal count of the state so x0_to_x1 can use it
            state = 0
            # foreach state x0 of the network...
            for binstate in self.enumerate_binary_states(self.numunits):
#                print "state=%s binstate=%s" % (state, binstate )
                self.x0_to_x1[state] = self._nextstate(binstate, cxnmatrix, nodefunc)
                state += 1
            
#            logging.debug( "x0_to_x1=%s" % self.x0_to_x1 )

        else:
            raise ValueError, "You messed up soemthing.  Did you specify a node function?"

        ##################################################
        

        # 1. Create the list of all valid output states.
        
        self.VALID_OUTPUTS = list(set([ self._binstate2dec(x) for x in self.x0_to_x1.values() ]))
        self.VALID_OUTPUTS.sort()
#        logging.debug( "VALID_OUTPUT_STATES=%s" % self.VALID_OUTPUTS )
        
        # 0. Make the stats for each output state
        self.x1_STATS = {}
        for x1 in self.VALID_OUTPUTS:
            self.x1_STATS[x1] = {}
            self.x1_STATS[x1]['possible_x0s'] = []
        
        for x0, x1state in self.x0_to_x1.iteritems():
            x1 = self._binstate2dec(x1state)
            self.x1_STATS[x1]['possible_x0s'].append( x0 )                 
        

        # make the main H_X0_GIVEN_x1 first, then we'll dool out the entries to VOs.        
        H_X0_GIVEN_x1 = self._make_H_X0_given_x1( self.x0_to_x1 )

        
        # the prob_this_x1 is determined by summing the probability of each x0 state the leads to it.
        for x1 in self.VALID_OUTPUTS:
            self.x1_STATS[x1]['prob_this_x1'] = 1.0/self.numstates * len( self.x1_STATS[x1]['possible_x0s'] )
            self.x1_STATS[x1]['H_X0_given_x1'] = H_X0_GIVEN_x1[x1]

#        logging.debug( self.x1_STATS )


        # now time to make the STATE TRANSITION MATRIX
        self.tmatrix = npy.zeros( (self.numstates,self.numstates) )
        # foreach x1 state...
        for x1, stats in self.x1_STATS.iteritems():
            for x0 in stats['possible_x0s']:
                self.tmatrix[x0][x1] = 1
        
        # store the cxnmatrix so we can calculate the number of external wires going into a node inside of a part
        self.cxnmatrix = cxnmatrix
        
        
    def external_wires_into_nodes_of_part(self, part ):
        """returns a 2-d array of the number of external wires going into each node of part        
Ex: cxnmatrices.fiver
>>> myphi = phi(N=5, cxnmatrix=array([[ True, False, True, False, True], [ True, True, True, True, True], [ True, False, True, False, True], [False, True, False, True, False], [False, True, False, True, False]]), nodefunc='and2')
>>> myphi.external_wires_into_nodes_of_part( array([1, 2]) )
[[3, 4], [0]]
>>> myphi.external_wires_into_nodes_of_part( array([0, 3, 4]) )
[[1, 2], [1], [1, 2]]
        """
        ######## CHECKS
#        if not part:
        if not part.size:
            print "thispart=%s" % part
            raise ValueError, "part cannot be blank."

        # check to see if part is not sorted.
        prev_ele = None
        for ele in part:
            if ele < prev_ele and ele is not None:
                raise ValueError, "Your part list should be sorted first."
            prev_ele = ele
        ######## END CHECKS
                
        
        # origin elements = everything outside of the part = rows
        source_nodes = [ x for x in xrange(self.numunits) if x not in part ]

#        logging.debug("source_nodes=%s" % source_nodes)
#        logging.debug("receiving_nodes=%s" % part )

        # receiving elements = the part = columns
#        externalwires = array(self.cxnmatrix[source_nodes][:,part].T)
        externalwires = self.cxnmatrix[source_nodes][:,part].T

        # OPTIMIZATION -- use the 'npy.where' instead of this for loop to make this faster.
        indices_external_wires_going_into_partnodes = []
        # for each row in the externalwires
        for row in externalwires:
            indices_going_into_partnode = [ source_nodes[i] for i, weight in enumerate(row) if weight ]
            indices_external_wires_going_into_partnodes.append( indices_going_into_partnode )

#        logging.debug( "external wire indices=%s" % indices_external_wires_going_into_partnodes )
        
        return indices_external_wires_going_into_partnodes



    def _make_H_X0_given_x1( self, x0_to_x1 ):
        x1s = [ self._binstate2dec(x) for x in x0_to_x1.values() ]
#        logging.debug("_make_H_X0_given_x1: x1s=%s" % x1s )
        
        H_X0_GIVEN_x1={}
        for x1 in self.VALID_OUTPUTS:     # for each valid output state...
#            print "trying x1=%s" % x1
            # H(X0|x1) = log2( number of times x1 was an output )
            H_X0_GIVEN_x1[x1] = log2( x1s.count(x1) )
        
#        logging.debug( "H_X0_GIVEN_x1=%s" % H_X0_GIVEN_x1 )

        return H_X0_GIVEN_x1
        

    def _make_prob_x0x1( self, x0_to_x1 ):
        """ returns a joint matrix with x0 on columns and all valid x1s on the rows.  Each entry is p(X0=x0 and X1=x1) . """
        
        prob_x0x1 = npy.zeros( (len(self.VALID_OUTPUTS),self.numstates) )
#        prob_x0x1 = npy.zeros( (self.numstates,self.numstates) )
        
        # all x1s with repeats
        x1s = [ self._binstate2dec(x) for x in x0_to_x1.values() ]
        
        # for each input/output pair...
        for inputstate, outputstate in x0_to_x1.iteritems():
            # convert the outputstate to it's decimal value 
            print "\t Doing input=%s output=%s" % (inputstate, outputstate )
#            logging.debug( "\t Doing input=%s output=%s" % (inputstate, outputstate ) )
            outputstate = self._binstate2dec(outputstate)
            outputstate_index = self.VALID_OUTPUTS.index(outputstate)
#            outputstate_index = outputstate

            # set the p(X0=x0 and X1=x1).
            # p(X0=x0 and X1=x1) = p(X1=x1|X0=x0) * p(X0=x0)
            #                    = p(X1=x1|X0=x0) * 1.0/2^N         // p(X0=x0) is constant.
            #                    = p(X1=x1|X0=x0) * 1.0/numstates   // p(X0=x0) is constant.            
            #                    = 1 * (1.0/2^N)                    // p(X1=x1|X0=x0) is 1.0 because we're only doing the correct input/output pairs.
            # Finally, recall that in this code the input_index = inputstate.
            print "\t -setting prob_x0x1[%s][%s] = %s" % ( outputstate_index, inputstate, 1.0 / self.numstates )
            prob_x0x1[outputstate_index][inputstate] = 1.0 / self.numstates
            
        pp.pprint( prob_x0x1 )
            
        return prob_x0x1

        
    def _make_prob_x0_given_x1( self, x0_to_x1 ):
        """ returns a matrix storing the conditional probability of each x0 given x1."""
        prob_x0_given_x1 = npy.zeros( (len(self.VALID_OUTPUTS),self.numstates) )
        x1s = [ self._binstate2dec(x) for x in x0_to_x1.values() ]
        
        # for each input/output pair...
        for inputstate, outputstate in x0_to_x1.iteritems():
            # convert the outputstate to it's decimal value 
#            logging.debug( "\t Doing input=%s output=%s" % (inputstate, outputstate ) )
            outputstate = self._binstate2dec(outputstate)
            outputstate_index = self.VALID_OUTPUTS.index(outputstate)

            # set the p(X0=x0|X1=x1).
            # p(X0=x0|X1=x1)    = p(X1=x1|X0=x0) * p(X0=x0) / p(X1=x1) 
            #                   = 1 * (1.0/2^N) / p(X1=x1) = (1.0/2^N) / SUM_x0 { p(X0=x0) * p(X1=x1|X0=x0) }
            #                   = (1.0/2^N) / SUM_x0 { (1.0/2^N) * 1 } = (1.0/2^N) / SUM_x0 { (1.0/2^N) }
            #                   = (1.0/2^N) / [ (#instances that lead to x1)*(1.0/2^N) ]
            #                   = 1.0 / (#instances that lead to x1)                      
            # Finally, recall that in this code the input_index = inputstate.
#            logging.debug( "\t -setting prob_x0_given_x1[%s][%s] = %s" % ( outputstate_index, inputstate, 1.0 / x1s.count( outputstate ) ) )
            prob_x0_given_x1[outputstate_index][inputstate] = 1.0 / x1s.count( outputstate )
            
        return prob_x0_given_x1
                    
    def _nextstate( self, current_state, cxn=None, nodefunc=None ):
        """Runs the update rule for a given starting state current_state. """
        
        if nodefunc is not None and type(nodefunc) is not type(lambda x: x):
            apply_gate = self.NODE_FUNCTIONS[nodefunc]
        else:
            apply_gate = self.APPLY_GATE
        if cxn is None:
            cxn = self.cxnmatrix

        assert len(current_state) == cxn.shape[0] == cxn.shape[1], "current_state not same length as cxnmatrix"
        assert type(current_state) is npy.ndarray, "current_state must be a numpy array"
#        assert type(cxn) is npy.core.defmatrix.matrix, "cxnmatrix must be a numpy matrix"

        # APPLY_GATE returns an array.  npy.dot multiplies the current_state by the connection matrix.        
        nextstate = apply_gate( npy.dot( current_state, cxn ) )
        assert type(nextstate) is npy.ndarray, "nextstate must be a numpy array.  type(nextstate)=%s" % (type(nextstate))

#        print "%s -> %s; " % (self._binstate2dec(current_state), self._binstate2dec(nextstate) ),

#        logging.debug("current_state=%s input=%s  nextstate=%s" % (current_state, ii, nextstate) )
        return nextstate
        
        
    def _dec2binstate(self, dec, padding=None ):
        """returns the binary of integer n"""
        if padding is None:
            padding = self.numunits

        # briefly tried the algorithm:
#        return array([(dec >> y) & 1 for y in range(padding-1, -1, -1)], dtype=npy.int8)
        # but it was slower than the one below
        return array( list( npy.binary_repr( dec, width=padding ) ), dtype=npy.int8 )


    
    def _binstate2dec(self, binstate):
        """returns the decimal value of a binary state."""

        # another algorithm...
        return long( ''.join(map(str,binstate)), 2 )
        
    def _is_valid_cxnmatrix( self, cxnmatrix ):
        """ returns True or False whether a passed matrix is a valid connection matrix."""
        
        # Try converting to an honest matrix and see if it errors.
        try:
            npy.matrix( cxnmatrix )
        except:
            return False
        
        # is the matrix square?
        if cxnmatrix.shape[0] != cxnmatrix.shape[1]:
            return False

#        logging.debug("cxnmatrix=%s" % cxnmatrix)
        
        return True

    def _cm2fastcm( self, cm ):
        """Converts a valid connection matrix to a fast connection matrix of booleans"""
        assert self._is_valid_cxnmatrix( cm ), "_cm2fastcm received an invalid connection matrix."
        
#        return array(cm, dtype=bool)
#        return array(cm, dtype=npy.int8)
        return array(cm, dtype=npy.float)
        # this must be a matrix because we're going to do alot of matrix mutliplcation, etc.
#        return npy.matrix(cm, dtype=bool)


    def dec2label( self, num, label=plotslabel, pad_to=None ):
        """Converts a decimal number to a label for plotting.  Has an option 'pad_to' variable that specifies the padding."""

#        # if no label is specified, use the default label.
#        if label is None:
#            label = self.plotslabel
            
        # what did to use for all of the bases.
        basestring='0123456789abcdefghijklmnopqrstuvwxyz'
        prepend = ''
    
        if label == 'decimal':
            return num
        elif label == 'binary':
            n, prepend = 2, ''
        elif label == 'hex' or label == 'hexadecimal':
            n, prepend = 16, 'x'
        elif type(label) is int and 2 <= label <= 36:
            n, prepend = label, label + '*'
        else:
            raise ValueError, "Passed label '%s' is not a valid label." % label
    
        labelstring=''
        current=num
    
        while current != 0:
           remainder=current%n       
           remainder_string = basestring[remainder]

           labelstring=remainder_string+labelstring
           current=current/n

        # if num is a zero, and the labelstring isn't set, set it to '0'
        if num == 0 and labelstring == '':
            labelstring = '0'

        # do the padding, if we're doing that.  But don't do any padding for decimals.  That's silly.
        if pad_to and label != 'decimal':
            padding = math.log( pad_to, n )
            labelstring = labelstring.rjust( int(padding), '0')
    
        # prepend any code we wanted to add.
        labelstring = prepend + labelstring
    
    
        return labelstring    



    def _is_valid_parts( self, parts ):
        """Returns True or False for whether a given list of parts is a valid bipartition of the sytem."""

        alleles = npy.concatenate( parts )
#        print "alleles=%s" % alleles

        # If any of the elements are not integers, return False
        if [ x for x in alleles if x != int(x) ]:
            raise TypeError, "Parts must be composed solely of ints."
            return False
        
        # The two parts must be disjoint, and have no duplicates, and pave the entire system
        if len(alleles) != len(set(alleles)) != self.numunits:
            raise ValueError, "Parts are not disjoint or have duplicates."
            return False

        # Passed all tests.  Is valid set of parts.
        return True

    def _make_prob_mu0_given_mu1( self, part, mu0, mu1, perturb='STATES' ):
        """ This algorithm is faster than summing over all states of the probability matrix p(M0=mu0,M1=mu1).  Instead, this algorithm:
        
        1. Calculates p(mu1|M0) for all mu0 in M0
        2. Calculates p(mu1 AND M0 ) by multiplying 1.0 / 2^{ numelements in M0 } = 1.0 / 2^{ |part| }
        3. Calculates p(mu1) by summing all of #2
        4. Calculates p(mu1 AND mu0 ) by selecting an element from #2
        
        User can specify to perturb either 'wires' or 'states'
        """
        numstates_in_part = 2**len(part)

        # 1. get the p(mu1|M0)
        prob_mu1_given_M0 = self._prob_mu1_given_M0( part, mu1, perturb=perturb )
        assert numstates_in_part == len(prob_mu1_given_M0), "num states in part should be the same as numelements in prob_mu1_given_M0 !"
                
        # p( mu1 AND M0 ) = p(M0=i) * p(mu1|M0=i) = 1.0/numstates_in_part * p(mu1|M0=i) = p(mu1|M0=i) / numstates_in_part
        # 2. p(M1=mu1) = Sum over intersection matrix prob_mu1_M0        
        prob_mu1_M0 = prob_mu1_given_M0 / float(numstates_in_part)


        # 3. p(mu0|mu1) = p(M0=mu0 AND M1=mu1) / p(M1=mu1)
        prob_mu0_given_mu1 = prob_mu1_M0[mu0] / npy.sum( prob_mu1_M0 )
        return prob_mu0_given_mu1        
                

    def _prob_mu1_given_M0( self, part, mu1, perturb ):
        """This function simply calls either the wires or the states version. """
        
        perturb = perturb.lower()
        
        if perturb == 'wires':
            return self._prob_mu1_given_M0_wires( part, mu1 )
        else:
            return self._prob_mu1_given_M0_states( part, mu1 )

    def _prob_mu1_given_M0_states( self, part, mu1 ):
        """ This returns p(mu1|M0) by perturbing the WIRES.
        This algorithm is faster than summing over all states of the probability matrix p(M0=mu0,M1=mu1).  Instead, this algorithm:        
            Calculates p(mu1|M0) for all mu0 in M0 via:
            1. Make p(X0 and mu1) by selecting over matrix p(X0 and X1) for all states x1 matching mu1...
            2. Make p(M0 and mu1) by selecting over matrix p(X0 and mu1) for all mu0states...
            3. p(mu1|M0) = p(M0 and mu1) / p(M0)
        """
        num_partnodes, numstates_in_part = len(part), 2**len(part)
        mu1state = self._dec2binstate( mu1, padding=num_partnodes )

#        print term.WHITE + "Calculating p(M0|mu1) for part=%s  mu1state=%s" % ( str(part), mu1state ) + term.NORMAL

        # discover all indices from 0...2^N with binaryform matching the mu1state
        x1_indices = [ index for index, binstate in enumerate(self.enumerate_binary_states(self.numunits)) if npy.all(binstate[part] == mu1state) ]
#        print term.YELLOW + "matching indices: %s" % x1_indices + term.NORMAL

        # p(X0 and X1) = (self.tmatrix * 1.0/self.numstates)
        prob_X0mu1 = (self.tmatrix * 1.0/self.numstates)[:,x1_indices]
#        print "p(X0mu1)[%s]..." % str(x1_indices)
#        pprint( prob_X0mu1 )
#        raw_input('--')

        # Now, create p(M0 and mu1) by iterating over all mu0 states, and summing over p(X0 and mu1)
        prob_M0mu1 = npy.zeros( (numstates_in_part,1) ).flatten()
        
        for mu0, mu0state in enumerate(self.enumerate_binary_states(num_partnodes) ):
#            print "mu0=%s" % mu0
            x0_indices = [ index for index, binstate in enumerate(self.enumerate_binary_states(self.numunits)) if npy.all(binstate[part] == mu0state) ]
#            print "x0_indices=%s" % x0_indices
#            print "p(M0 and mu1)..."
#            pprint( prob_M0mu1 )
            prob_M0mu1[mu0] = npy.sum( prob_X0mu1[x0_indices] )
        
#        print "\t p(M0 and mu1)=%s" % prob_M0mu1
#        raw_input('--')
        
        # Now make p(mu1|M0) via p(mu1 and M0) / p(M0)        
        return (prob_M0mu1 * 1.0/numstates_in_part)
                
    def _prob_mu1_given_M0_wires( self, part, mu1 ):
        """ This returns p(mu1|M0) by perturbing the WIRES.
        This algorithm is faster than summing over all states of the probability matrix p(M0=mu0,M1=mu1).  Instead, this algorithm:        
            Calculates p(mu1|M0) for all mu0 in M0 via:
            1. Iterating over all nodes of M...
            2. Iterating over all of the states of mu0...
            3. Calculates how often the updated node matches mu1.
        """
        
        num_partnodes, numstates_in_part = len(part), 2**len(part)
        
        # get the index of the external wires into each node of the part
        wire_indices_ofeach_partnode = self.external_wires_into_nodes_of_part( part )

        # set the mu1state.
        mu1state = self._dec2binstate( mu1, padding=num_partnodes )

        # define the vector p(mu1|M0).        
        # this vector will be filled in and returned.
        # We will gradually fill in this by multiplying the probabilities of each node matching mu1.
        # Because we are getting the product, the vector is initialized with all 1.0s
        prob_mu1_given_M0_vector = npy.ones( numstates_in_part )
        
        # for each node of mu1...
        # Calculate p(node=mu1[i]|mu0)
        for node_partindex, mu1_for_partnode in enumerate(mu1state):
            # we have to do to reset the shape of the x0states from the previous node to make the shape for the current node.
            # this must be zeros, not npy.empty.  We leave all non-mu0-related states at zero.
            x0states = npy.zeros( (1,self.numunits) ) 

            # get all of the incoming connections into this node
            node_sysindex = part[node_partindex]
            cxnmatrix_node = self.cxnmatrix[:,node_sysindex]
            
            wire_indices = wire_indices_ofeach_partnode[node_partindex]
#            print "wire_indices=%s" % wire_indices


            num_outside_wires = len( wire_indices )
            num_wirestates = 2**num_outside_wires
            # First, expand x0states to accomodate all of the different wirestates.
            # We convert x0states from a single vector to a big-matrix of states -- one row for each wirestate.
            # Now get all possibilities of the wirestates, and set the wire_indices equal to all possibilities of the wirestates.
            if num_outside_wires:
                # Yes, 2**num_outside_wires is correct.  Because npy.repeat( x, N, axis=0 ) returns N rows.
                x0states = npy.repeat( x0states, num_wirestates, axis=0 )
                x0states[:,wire_indices] = self.binary_states_matrix( num_outside_wires )
                                                            

            # p( node = mu1[i] | mu0 )  = 1.0/2^{num_outside_wires} * Sum_over_all_states_of_outside_wires { p(node=mu1[i]| mu0,wirestate) }
            #                           = 1.0/2^{num_outside_wires} * < number_times_node=mu1[i] >
            #                           = < number_times_node=mu1[i] > / 2^{num_outside_wires}
            # for each possible mu0state...
            for mu0, mu0state in enumerate(self.enumerate_binary_states(num_partnodes) ):
            
                # OPTIMIZATION: If the probability of this state matching mu0 is already zero,
                # then we know the product is zero.  So don't bother trying this state for any of the other nodes.
                if prob_mu1_given_M0_vector[mu0] == 0:
                    continue

                # set all of the part indices in x0states to the mu0state.
                x0states[:,part] = mu0state
                        
                # now multiply the x0states by the vector of incoming cxns into the node, and see how many matched the mu1state for this node
                # we then get the probability of this node matching the mu1state by dividing by number of wirestates
                # Each node is treated as independent from the other nodes.  Thus to get the probability
                # of mu1 we multiply the probability that each node is correct.
                # The probability of p(mu1|mu0) = product{ p(nodes|mu0) }
#                prob_node_matches_mu1_given_mu0 = float(num_matches_to_mu1) / num_wirestates
#                num_matches_to_mu1 = npy.sum(self.APPLY_GATE( npy.dot(x0states, cxnmatrix_node) ) == mu1_for_partnode)
#                prob_mu1_given_M0_vector[mu0] *= prob_node_matches_mu1_given_mu0

                # Below is the one line version of the above.
                prob_mu1_given_M0_vector[mu0] *= npy.sum(self.APPLY_GATE( npy.dot(x0states, cxnmatrix_node) ) == mu1_for_partnode) / float(num_wirestates)

        # Now normalize so that we have a probability distribution, and done!
        prob_mu1_given_M0_vector /= npy.sum(prob_mu1_given_M0_vector)
        
        # replace all instances of NaN with 0.0
        prob_mu1_given_M0_vector[npy.isnan(prob_mu1_given_M0_vector)] = 0.0
        
        # Have now made p(mu1|M0)            
        return prob_mu1_given_M0_vector

    def _weave_mu0_and_wirestate_into_x0(self, part, mu0state, wire_indices, wirestate):
        """Returns the system state (x0) given the state of a part, and the state of the external wires.  If there are any leftover states they are filled in with zeros."""
        
        # create an empty array
        x0state = array( [0] * self.numunits, npy.int8 )

        # fill in the part units with mu0
        for node_partindex, node_sysindex in enumerate(part):
#            print term.YELLOW + "\t\t\t tweave: setting x0[%s] = mu[%s] = %s" % (node_sysindex,node_partindex,mu0state[node_partindex]) + term.NORMAL
            x0state[node_sysindex] = mu0state[node_partindex]

        # fill in the externalwire units with wirestate
        if wire_indices is not None and wirestate is not None:
            for externalwireindex, sysindex in enumerate(wire_indices):
                #            print term.RED + "\t\t\t weave: setting x0[%s] = wire[%s] = %s" % (sysindex,externalwireindex,wirestate[externalwireindex]) + term.NORMAL        
                x0state[sysindex] = wirestate[externalwireindex]

        # and the rest are all already zeros.
#        print "x0state=%s" % x0state
        return x0state

        phi = self.ei( x1, MIP, perturb=perturb )
                
    def ei( self, x1, parts=None, cx=None, perturb='STATES' ):
        """This function returns the EFFECTIVE INFORMATION that is generated when a system enters state x1.
        
        Alternatively, if passed parts, returned the effective information generated above and beyond the information created
        by the parts.        
        """
        #################################################################
        # 1. Do some checks and conversions on the inputs.
        #################################################################
        x1, perturb = long(x1), perturb.upper()

        # User passed an invalid x1 state.  Return None.
        if x1 not in self.VALID_OUTPUTS:
            return None

        if parts is not None and not self._is_valid_parts(parts):
            print "_is_valid_parts(parts) = %s" % self._is_valid_parts(parts)
            sys.stdout.flush()
            raise ValueError, "Parts specifiction was not valid.  parts=%(parts)s" % locals()
        if cx is not None:
            raise ValueError, "Are not ready to do things across Complexes."

        if not perturb in ('WIRES','STATES'):
            raise ValueError, "You passed an invalid perturb.  Valid perturbs are 'WIRES' and 'STATES'."

        #################################################################
        # 2. Are we calculating ei() for the entire system (or across the total partition)?  If so, use the simple equation.
        # ei( X0 -> x1 )    = H[ X0|X1=x1,pmax(X0) ] - H(X0|X1=x1)
        #                   = numunits - H[X0|X1=x]
        #################################################################        
        if parts is None or max(map(len,parts)) == self.numunits:
#            print term.YELLOW + "\tDoing ei(%s/%s) as whole..." % (x1, parts) + term.NORMAL
            return (self.numunits - self.x1_STATS[x1]['H_X0_given_x1'])
        
        #################################################################
        # 3.  If we're calculating the EI generated BEYOND that of the parts, then we have two parts to do it.
        #     Here we decide whether to perturb the WIRES or the STATES
        #################################################################

        # parts will ALWAYS be a tuple of numpy arrays.
        parts = tuple(map(array,parts))
        
        if perturb == 'WIRES':
            return self._ei_wires( x1, parts, cx )
        else:
            return self._ei_states( x1, parts, cx )


    def _ei_wires( self, x1, parts, cx=None ):
        """This function returns the effective information (EI) that is generated when a system enters state x1.
            x1 = calculates the EI of enterting state x1
            parts = If specified, calculates the EI *beyond* the EI generated by the parts.
        """

        assert len(parts) == 2, "Number of parts must be 2."

        # sort the parts
        part1, part2 = parts
        part1.sort()
        part2.sort()
        parts = (part1,part2)
        # 2. Computing ei() across a particular partitionpy.
        # MASTER EQUATION:  Sum_x0_in_X0{ p(X0=x0|X1=x1) * -log( PROD{ p(M0=mu0|M1=partstate) } ) } - H[X0|X1=x1]
        #                   assuming all input states are at uniformly probability, the probability of an input state given the output is just 1.0/(# inputs that lead to x1)
        #                   = 1./<#x0s that lead to x1> * Sum_x0_in_X0{ -log( PROD{ p(M0=mu0state|M1=mu1state) } ) } - H[X0|X1=x1]
                
        x1state = self._dec2binstate(x1)
        # PART I: Sum_x0_in_X0{ ...
#        Sum_x0_in_X0 = 0
        productofparts = 1.0
        # For each x0 that leads to this x1...        
        for x0 in self.x1_STATS[x1]['possible_x0s']:
                        
            # Now to compute the product of the part probabilities.
            x0state = self._dec2binstate(x0)
            
            # PART II: PROD{ p(M0=mu0|M1=partstate)
            for part in parts:

                mu0state, mu1state = x0state[part], x1state[part]
                mu0, mu1 = self._binstate2dec(mu0state), self._binstate2dec(mu1state)

                prob_mu0_given_mu1 = self._make_prob_mu0_given_mu1( part, mu0, mu1 )

                # PART III: PROD{ p(M0=mu0|M1=partstate) }
                productofparts *= prob_mu0_given_mu1
            
            # PART IV: Sum_x0_in_X0{ p(X0=x0|X1=x1) * -log( PROD{ p(M0=mu0|M1=partstate) } ) }
#            Sum_x0_in_X0 += -log2(productofparts)
        try:
            Sum_x0_in_X0 = -log2(productofparts)
        except OverflowError:
            print "HAD Error calculating -log2(productofparts)!  productofparts=%s" % productofparts

        # PART VI: - H[X0|X1=x1]
        ei_across_partition = (2**(-self.x1_STATS[x1]['H_X0_given_x1']) * Sum_x0_in_X0) - self.x1_STATS[x1]['H_X0_given_x1']

#        print "  v: H[X0|X1=x]=%s" % self.x1_STATS[x1]['H_X0_given_x1']
#        print "  v: ei(%s/P)=ei(%s/%s)=%s - %s = %s" % ( x1, self._dec2binstate(x1), parts, Sum_x0_in_X0, self.x1_STATS[x1]['H_X0_given_x1'], ei_across_partition )

        assert ei_across_partition >= 0, "Eeep!  Had invalid ei(x1,P)!  ei(%s/P)=ei(%s/%s)=%s - %s = %s" % ( x1, self._dec2binstate(x1), parts, Sum_x0_in_X0, self.x1_STATS[x1]['H_X0_given_x1'], ei_across_partition )

        return ei_across_partition

    def _ei_states( self, x1, parts, cx=None ):
        """This function returns the EI generated *beyond that* the EI that is generated by the parts.
        
        This calculation is done by perturbing the STATES.        
        """

        ################################################################################################################        
        # 2. Computing ei(x,P)
        # MASTER EQUATION:  Sum_x0_in_X0{ p(X0=x0|X1=x1) * -log( PROD{ p(M0=mu0|M1=partstate) } ) } - H[X0|X1=x1]
        #                   assuming all input states are at uniformly probability, the probability of an input state given the output is just 1.0/(# inputs that lead to x1)
        #                   = 1./<#x0s that lead to x1> * Sum_x0_in_X0{ -log( PROD{ p(M0=mu0|M1=partstate) } ) } - H[X0|X1=x1]

        # We will *= this for each part.
        
        x1state = self._dec2binstate(x1)
        productofparts = 1.0
#       print 'x0s that lead to x1: %s' % [ (int(x0), list(self._dec2binstate(int(x0)))) for x0 in self.x1_STATS[x1]['possible_x0s'] ]        
        for x0 in self.x1_STATS[x1]['possible_x0s']:
#            print "possible_x0=%s" % x0

            # Now to compute the product of the part probabilities.
            x0state = self._dec2binstate(x0)
            
            # prob( mu0 | mu1 ) = p(mu0 and mu1) / p(mu1)
            # Calculate #matches to mu1, and calculate #matches to x0
            # Calculate PROD{ p(M0=mu0|M1=mu1)


            for part in parts:
                # get the mu0 and mu1

                instances_M1_equals_mu1, instances_M1_equals_mu1_AND_M0_equals_mu0 = 0, 0
                mu0state, mu1state = x0state[part], x1state[part]

                for tempx1 in self.VALID_OUTPUTS:

                    # found a match for this mu1!
                    if npy.all(self._state_of_part( tempx1, part ) == mu1state):
                        instances_M1_equals_mu1 += len(self.x1_STATS[tempx1]['possible_x0s'])
                    
                        # now for the #instances mu0 matched
                        instances_M1_equals_mu1_AND_M0_equals_mu0 += len( [ tempx0 for tempx0 in self.x1_STATS[tempx1]['possible_x0s'] if npy.all(self._state_of_part( tempx0, part ) == mu0state) ] )

                assert instances_M1_equals_mu1_AND_M0_equals_mu0 <= instances_M1_equals_mu1, "#mu0mu1 matches > #mu1 matches!  This is impossible!"

                prob_mu0_given_mu1 = float(instances_M1_equals_mu1_AND_M0_equals_mu0) / float(instances_M1_equals_mu1)

                productofparts *= prob_mu0_given_mu1
#                print "\t\t prob(mu0|mu1)=%s \t productofparts=%s" % ( prob_mu0_given_mu1, productofparts )

            # PART IV: Sum_x0_in_X0{ p(X0=x0|X1=x1) * -log( PROD{ p(M0=mu0|M1=partstate) } ) }
        try:
            Sum_x0_in_X0 = -log2(productofparts)
        except OverflowError:
            print "Error computing ei(x1=%s/P=%s) // x0s: %s" % (x1, str(parts), len(self.x1_STATS[x1]['possible_x0s']) )
            raise ValueError, "HAD Error calculating -log2(productofparts)!  productofparts=%s" % productofparts

        # PART VI: - H[X0|X1=x1]
#        print "parts=%s \t sum_x0_in_X0=%s" % (str(parts), Sum_x0_in_X0)
#        Sum_x0_in_X0 = ((1.0/len(self.x1_STATS[x1]['possible_x0s'])) * logterm)
        ei_across_partition = (2**(-self.x1_STATS[x1]['H_X0_given_x1']) * Sum_x0_in_X0) - self.x1_STATS[x1]['H_X0_given_x1']


#        print "  v: H[X0|X1=x]=%s" % self.x1_STATS[x1]['H_X0_given_x1']
#        print "  v: ei(%s/P)=ei(%s/%s)=%s*%s - %s = %s" % ( x1, self._dec2binstate(x1), parts, 1.0/len(self.x1_STATS[x1]['possible_x0s']), Sum_x0_in_X0, self.x1_STATS[x1]['H_X0_given_x1'], ei_across_partition )
        assert ei_across_partition >= 0, "Eeep!  Had invalid ei(x1,P)!  ei(%s/P)=ei(%s/%s)=%s - %s = %s" % ( x1, self._dec2binstate(x1), parts, Sum_x0_in_X0, self.x1_STATS[x1]['H_X0_given_x1'], ei_across_partition )

        return ei_across_partition

    def Minimum_Information_Bipartition(self, x1, totalpart=True, perturb='STATES' ):
        """ Returns the minimum information bipartition as a sequence of (partA, partB).
        Where partA, partB are lists.
        
        totalpart=True specifies to also check the total partitionpy.
        """
        
        # min_ei and min_norm are set at their maximum possible values.
        # These magic numbers are better than setting them to None because it doesn't have to do a is None check on min_ei
        min_ei, min_norm, min_parts = float(self.numunits), 1.0, None
        
        for parts in self.all_possible_bipartitions( totalpart=totalpart ):
#            print "trying parts=%s" % str(parts)
#            ei = self.ei( x1, parts=parts )
            ei = self.ei( x1, parts=parts, perturb=perturb )

            # Normalization = (K-1) * Hmax( size of smallest part ) 
            #               = (K-1) * (size of smallest part)
            #               = (size of smallest part)            // for K=2
            # In the case of the total partition, Normalization = (size of total partition)
            # the if x takes care of the case of the total partition (when there will be an empty part, [], which evaluates false
            norm = min([ len(part) for part in parts if part.size ])

            # if the normalized_ei is lower, use that.  If the normalized ei's are equal, use the one with the smaller unnormalized ei.
            if (ei/norm < min_ei/min_norm) or (ei/norm == min_ei/min_norm and ei < min_ei):
                min_ei, min_norm, min_parts = ei, norm, parts

                # OPTIMIZATION.  We know that the lowest ei(x,P) can be is zero.  Thus, we get a value of ei(x1,P) == 0, stop searching.
                if min_ei == 0.0:
                    break
                                
        # we want to return python lists, not numpy arrays.
        return tuple(map(list,min_parts))
        
    def Minimum_Information_Partition(self, x1, perturb='STATES' ):
        """ Returns the minimum information bipartition as a sequence of (partA, partB, partC...).
        Where each part is a list.
        
        totalpart=True specifies to also check the total partitionpy.
        """
        
        # min_ei and min_norm are set at their maximum possible values.
        # These magic numbers are better than setting them to None because it doesn't have to do a is None check on min_ei
        min_ei, min_norm, min_parts = float(self.numunits), 1.0, None

        # BELL_all_partitions assigns ALL POSSIBLE PARTITIONS to BELL_biglist which we can then iterate over.
        global BELL_biglist
        BELL_biglist = []
        BELL_all_partitions( [], range(self.numunits) )
        print "Attempting %s partitions..." % len(BELL_biglist)

        
        for parts in BELL_biglist:
            parts = map( array, parts )
#            print "trying parts=%s" % str(parts)
            ei = self.ei( x1, parts=parts, perturb=perturb )


            # Normalization = (K-1) * Hmax( size of smallest part ) 
            #               = (K-1) * (size of smallest part)
            # In the case of the total partition, Normalization = (size of total partition)
            # take care of the totalpartition case.
            if parts is None or max(map(len,parts)) == self.numunits:
                norm = self.numunits
            else:
                # take care of the other cases partitionpy.
                norm = (len(parts)-1) * min(map(len,parts))

            # if the normalized_ei is lower, use that.  If the normalized ei's are equal, use the one with the smaller unnormalized ei.
            if (ei/norm < min_ei/min_norm) or (ei/norm == min_ei/min_norm and ei < min_ei):
                min_ei, min_norm, min_parts = ei, norm, parts
                print term.YELLOW + "\tei(%s/P=%s) \t =%s \t / norm=%s \t =%s" % (x1, str(map(list,parts)), ei, norm, ei/norm) + term.NORMAL

                # OPTIMIZATION.  We know that the lowest ei(x,P) can be is zero.  Thus, we get a value of ei(x1,P) == 0, stop searching.
                if min_ei == 0.0:
                    break

        # we want to return the min_parts as a list, not a numpy array
        del BELL_biglist
        
        return tuple(map(list, min_parts))

    def binary_states_matrix( self, elements ):
        """ This function does the same thing as enumerate_binary_states,
        except that it returns the result as a single matrix rather than one at a time.
        Each binary number is returned on it's own row.
        It fills in the binary matrix "visually" using the Balduzzi method, a pretty darn fast method.
        """
        mat = npy.empty( (2**elements, elements), npy.int8 )
        # loop from [n-1, 0]
        # for each column in the matrix...
        for colindex in xrange(elements):
            length = elements - colindex - 1
            mat[:,colindex] = npy.tile( npy.concatenate( (npy.zeros(2**length), npy.ones(2**length)) ), 2**(colindex) )
        return mat
        
        # Alternative method for calculating the binary_states_matrix from Virgil.
        # Method is slightly slower than Balduzzi's "visual fill"
        '''
        endat = 2**elements
        return array( map(list, map(npy.binary_repr,range(startat,endat),[elements]*(endat-startat)) ), dtype=npy.int8)        
        '''

    def enumerate_binary_states( self, elements, startat=0 ):
        """ This function enumerates all binary values of size npy.  Returns a lists of binary values 0...2^n-1."""
        
        if not int(elements) >= 1:
            raise ValueError, "elements must be >= 1."

        elements=int(elements)
        
        if startat >= 2**elements:
            raise ValueError, "startat must be < 2**elements"


        pos=startat
        endat=(2**elements)-1
        while pos <= endat:
            binstate = array( list( npy.binary_repr(pos,width=elements) ), dtype=npy.int8 )
#            print "yeilding %s" % binstate
            yield binstate
            pos += 1
        '''
        # start with binary string for where we start at.
        state = self._dec2binstate(startat, padding=elements )
        yield state
        
        while 1:
            # get the last zero in state.  Call it's index 'finalzero'
            try:
                finalzero = npy.max( npy.where(state==0) )
            except ValueError:
                # if no zeros were found that means it's composed of all ones.
                # That means we're done.  break.
                break

            # set the index at finalzero to a 1, everything following it should be zero.
            state[finalzero] = 1
            state[finalzero+1:] = npy.zeros(elements-finalzero-1, dtype=npy.int8)
            yield state
        '''
        
    
    def all_possible_bipartitions(self, totalpart=True):
        """ Returns a list of all possible bipartitions of the units."""

        # iterate across all binary digits from [1, 2^(self.numunits-1)-1]
        # the number of 1s in each digit will be partB's elements.
        
        for binaryarray in self.enumerate_binary_states( self.numunits-1, startat=1 ):

            # partB is the number of 1s in each digit
            # partA is everything in the system not in partB
            partB = npy.where(binaryarray==1)[0]
#            print "B=%s" % partB
#            print "type(B)=%s" % type(partB)
#            partA = [ x for x in range(self.numunits) if x not in partB ]
            partA = npy.append( npy.where(binaryarray==0), self.numunits-1 )

            # if partA is smaller than partB, flip them
            if len(partA) < len(partB):
                partA, partB = partB, partA

            # return this partA, partB
            yield (partA, partB)
            
        # if we're doing the total partition, return that one now.
        if totalpart:
            yield ( npy.arange(self.numunits), )

    def _state_of_part( self, binstate, part ):
        """ This function converts binstate to a binary array and returns the slice representing the state of the part."""
        # if the user passes a statearray instead of a decstate, use that state instead.
        if type(binstate) is long or type(binstate) is int:
            binstate = self._dec2binstate( binstate )

        return binstate[part]

    def _state2mask( self, binstate, part ):
        """ This function converts binstate to a binary array and returns the slice representing the state of the part."""
        # if the user passes a statearray instead of a decstate, use that state instead.
        if type(binstate) is long or type(binstate) is int:
            binstate = self._dec2binstate( binstate )
        
#        part=list(part)
#        print "part=%s" % part
        # make a list with ones for the bits in part
        bitsIcareabout = npy.zeros( len(binstate), npy.int8 )
        for i, ele in enumerate(bitsIcareabout):
            if i in part:
                bitsIcareabout[i] = 1
#        print "\tbitsIcareabout=%s" % bitsIcareabout
        
        entrieswith1 = npy.zeros( len(binstate), npy.int8 )
        for i, ele in enumerate(binstate):
            if i in part and ele == 1:
                entrieswith1[i] = 1
#        print "\tentrieswith1=%s" % entrieswith1
        
        mask, entrieswith1 = ''.join(map(str,bitsIcareabout)), ''.join(map(str,entrieswith1))

        return mask, entrieswith1

    def average_ei( self, parts, perturb='STATES' ):
        """Computes the average effective information over all valid output states, \overline{ei}.  Or, computes \overline{ei}(P) if passed parts."""
#        print "%s: running parts=%s" % (perturb, str(parts))
        if not self._is_valid_parts(parts):
            raise ValueError, "average_ei: You passed invalid parts=%s" % parts

        H_X0_GIVEN_X1 = sum( [ self.x1_STATS[x1]['prob_this_x1'] * self.x1_STATS[x1]['H_X0_given_x1'] for x1 in self.VALID_OUTPUTS ] )
        
        # if passed no arguments, or passed the total partition, compute the average ei of the total system (total partition)
        # ei(x) = N - H[X0|X1=x] 
        # average_ei    = Sum_x1_in_X1 { p(X1=x1) (N - H[X0|X1=x]) } = N - Sum_x1_in_X1 { p(X1=x1) * H[X0|X1=x] } 
        #               = N - H[X0|X1]
        if parts is None or max( [ len(part) for part in parts if part.size ] ) == self.numunits:
            return float(self.numunits) - H_X0_GIVEN_X1

        # Computing the average_ei for a given partition P is:
        # Sum_i{ H[M^i_0|M^i_1] } - H[X0|X1]
        # Sum_i{ H[M^i_0|M^i_1] }   = p(M_0 AND M_1) * -log( p(M_0|M_1) ) 
        #                           = p(M_0 AND M_1) * -log( p(M_0 AND M_1) / p(M_1) )
        # average_ei(P) = Sum_i{ p(M_0 AND M_1) * -log( p(M_0 AND M_1) / p(M_1) ) } - H[X0|X1]
        sumparts = sum( [ self.M0_GIVEN_M1(part, perturb=perturb) for part in parts ] )
        average_ei = sumparts - H_X0_GIVEN_X1

#        print term.CYAN + "\t %s - %s" % (sumparts, H_X0_GIVEN_X1) + term.NORMAL
                            
        assert 0.0 <= average_ei <= self.numunits, "average_ei: Error!  average_ei must be between [0,%s].  However, it is %s!" % ( self.numunits, average_ei )
        return average_ei

    def M0_GIVEN_M1( self, part, perturb='STATES' ):
        """Calculates conditional entropy for a given part M, H[M_0|M_1], averaged over all states mu0 and mu1.
            H[M_0|M_1]              = SUM_mu0{ SUM_mu1{ p(M_0 AND M_1) * -log( p(M_0|M_1) ) } }
                                    = SUM_mu0{ SUM_mu1{ p(M_0 AND M_1) * -log( p(M_0 AND M_1) / p(M_1) ) } }
                                    
           As is calculates H[M0|M1] by perturbing the WIRES.
        """

        ###########################################################
        # 0. Do checking of inputs
        ###########################################################
        if type(part) is list:
            part = array(part)

#        print term.WHITE + "\t MO_GIVEN_M1: Calculating H[M0|M1] for part=%s" % (str(part)) + term.NORMAL
        # 1. Make   p(M_0 AND M_1)
        #           p(M_0 AND M_1) = p(M1|M0) * p(M0)
        prob_M1_given_M0 = npy.zeros( (2**len(part), 2**len(part)) )
        
        # foreach possible mu1, make p(mu1|M0)
        # each mu1 is a different column
        for mu1 in xrange(2**len(part)):
            # !!! This part determines whether we're perturbing the STATES or WIRES
            prob_M1_given_M0[:,mu1] = self._prob_mu1_given_M0( part, mu1, perturb=perturb )

        # renormalize p(M1|M0) to make each row sum to 1.0...
        prob_M1_given_M0 = (prob_M1_given_M0.T * 1.0/prob_M1_given_M0.sum(axis=1)).T
        
#        print "p(M1|M0)..."
#        pprint( prob_M1_given_M0 )

                
        prob_M0M1 = prob_M1_given_M0 * 1.0/2**len(part)
        assert npy.sum(prob_M0M1) == 1.0, "M0_GIVEN_M1: Error -- Sum over all p(M0 AND M1) should be 1.0.  But it is %s" % npy.sum(prob_M0M1)

#        print "p(M1 and M0)..."        
#        pprint( prob_M0M1 )

        # Don't p(M1|M0) anymore!  Clear it.
        del prob_M1_given_M0
        
        # 2. Make p(M_0|M_1) = p(M_0 AND M_1) * 1.0/p(M_1)
        # Making p(M_1)...
        prob_M1 = prob_M0M1.sum(axis=0)
        assert len(prob_M1) == 2**len(part), "M0_GIVEN_M1: Error -- length of prob_M1 should be 2**%s.  Shape is %s" % str(prob_M1.shape)

        # now take each column of prob_M0M1 and divide by that p(mu1).  Convert NaNs to zero.
        prob_M0_given_M1 = prob_M0M1 * 1.0/prob_M1
        prob_M0_given_M1[npy.isnan(prob_M0_given_M1)] = 0.0
        
#        print "p(M0|M1)..."
#        pprint( prob_M0_given_M1 )

        # 3. H[M0|M1] = SUM{ p(M_0 and M1) * -log( p(M0|M1) ) }
        assert prob_M0M1.shape == prob_M0_given_M1.shape, "Error -- p(M0 and M1) and p(M0|M1) should have same shape.  But their shapes are: %s and %s!" % (str(prob_M0M1.shape), str(prob_M0_given_M1.shape))

        # logterm = -log[ p(M0|M1) ].  Convert infs to zero.
        logterm = -log2(prob_M0_given_M1)
        logterm[npy.isinf(logterm)] = 0.0

#        pprint( logterm )
        
        H_M0_given_M1 = npy.sum( prob_M0M1 * logterm )

        # catch any roundoff error...
        if H_M0_given_M1 < 0.0:
            # if the H is really negative, then Error.  But if -0.001 < H < 0.0 then just set it to zero.
            assert H_M0_given_M1 >= -0.001, "Error -- H[M0|M1] must be >= 0.  However, it is: %s" % H_M0_given_M1
            H_M0_given_M1 = 0.0
        
#        print term.RED + "\t%s=%s" % (part, H_M0_given_M1) + term.NORMAL        
        return H_M0_given_M1
        
    
    def pull_partitions_from_file( self, partitionsfilename ):
        """This function reads partitionsfilename line by line yields the JSON encoded partitions."""
        
        if not os.path.exists( partitionsfilename ):
            return ValueError, "filename '%s' does not exist." % partitionsfilename
        
        partitionsfile = open( partitionsfilename, 'r' )
        
        for line in partitionsfile:
            parts = array(json.loads(line))
            print "yielding parts=%s" % parts
            raw_input('--')
            yield parts 
        
#        partitionsfile.close()
        
    def average_MIP( self, partitionsfilename=None, perturb='STATES' ):
        """Returns the MIP of ei(P).
            partitionsfilename = If None, checks all bipartitions.
                                 If a filename is passed, instead checks all partitions in the file.
        """
        min_ei, min_norm, min_parts = float(self.numunits), 1, None

        if partitionsfilename:
            print "reading partitions from file='%s'" % partitionsfilename
            raw_input('--')
            partitions_iterator = self.pull_partitions_from_file( partitionsfilename )
        else:
            partitions_iterator = self.all_possible_bipartitions
            
        for parts in partitions_iterator():
#            print "trying parts=%s..." % str(parts)

            # average_ei = average_ei for this part
            average_ei = self.average_ei( parts, perturb=perturb )

            # Normalization = (K-1) * Hmax( size of smallest part ) 
            #               = (K-1) * (size of smallest part)
            #               = (size of smallest part)            // for K=2
            # In the case of the total partition, Normalization = (size of total partition)
            # the if x takes care of the case of the total partition (when there will be an empty part, [], which evaluates false
            if parts is None or max(map(len,parts)) == self.numunits:
                norm = self.numunits
            else:
                # take care of the other cases partitionpy.
                norm = (len(parts)-1) * min([ len(part) for part in parts if part.size ])

            # if the normalized_ei is lower, use that.  If the normalized ei's are equal, use the one with the smaller unnormalized ei.
            if (average_ei/norm < min_ei/min_norm) or (average_ei/norm == min_ei/min_norm and average_ei < min_ei):
                min_ei, min_norm, min_parts = average_ei, norm, parts
                print term.RED + "\t parts=%s \t\t quotient: %s \t min_ei=%s \t min_norm=%s" % (str(parts), min_ei/min_norm, min_ei, min_norm) + term.NORMAL
            
                # OPTIMIZATION.  We know that the lowest ei(P) can be is zero.  Thus, we get a value of ei(P) == 0, stop searching.
                if min_ei == 0.0:
                    break

        # return the min_parts as a LIST, not an array
        return map(list,min_parts)
        
        
    def PHI( self, x1, perturb='STATES', allpartitions=False ):
        """ Returns the PHI, MIP of a particular state decstate. 
        
        totalpart = Try the totalpartition?
        perturb = Are we going to calculate the ei via perturbing the states or the wires?
        allpartitions = Try all Bell's Number of partitions.  If False, only attempt bipartitions.
        
        """

        if type(x1) is npy.ndarray or type(x1) is list:
            x1 = self._binstate2dec( x1 )

        # Make the Prob( node = x1[i] | X0 ) matrix
        
        # phi is defined as the ei across the minimum information (bi)partition
        if allpartitions:
            MIP = self.Minimum_Information_Partition( x1, perturb=perturb )
        else:
            MIP = self.Minimum_Information_Bipartition( x1, perturb=perturb )

#        print "MIP=%s  \t type(MIP)=%s" % (MIP, type(MIP))
        phi = self.ei( x1, MIP, perturb=perturb )
            
        return phi, MIP
        
def _test():
    """ Run all of our unit tests. """
    import doctest
    doctest.testmod()


def main():
    print "- Running standalone"    
    import cxnmatrices
#    import phi_balduzzi as Bphi
#    import timeit
    
#    tm = cxnmatricesvFig2b
    tm = cxnmatrices.and4
#    tm = cxnmatrices.testand
#    tm = cxnmatrices.fiver
#    tm = cxnmatrices.g2_nd
#    tm=cxnmatrices.checker6
#    tm=cxnmatrices.checker7
#    tm=cxnmatrices.c1
#    tm=cxnmatrices.g3_1d

    numunits=len(tm)
    print "tm=%s" % tm
    
    Vphi = phi(N=numunits, cxnmatrix=tm, nodefunc='and2')

    print "================================================================================"
    # for each partition

    for x0, x1 in Vphi.x0_to_x1.iteritems():
        x1 = Vphi._binstate2dec( x1 )
        print "%(x0)s: \t %(x1)s" % locals()
    pp.pprint( Vphi.x0_to_x1)

    # try every possible output state.
#    for x1 in xrange(2**numunits):
    '''
    for y in Vphi.VALID_OUTPUTS:
        for x in xrange(Vphi.numstates):
            print "RUNNING MAP REDUCE FOR x=%s y=%s" % ( x, y)
            Vphi.MapReduce( x, y )
    '''

    for x1 in Vphi.VALID_OUTPUTS:
        x1state = Vphi._dec2binstate( x1 )
        print "x1state=%s" % x1state
        
#        bphi, partB, partA = Bphi.phi(tm,range(numunits),x1state,T=-1,iself=True)
#        bphi = round(bphi, 5)
#        print "balduzzi (wires): \tMIP=%s \t PHI(x)=%s" % ( (partA,partB), bphi )
#        sys.stdout.flush()
        
        # get the phi(x1) and MIP        
        vphi, vMIP = Vphi.PHI(x1, totalpart=True, perturb='STATES')
        vphi = round(vphi,5)
        print "virgil (wires): \tMIP=%s \t PHI(x)=%s" % ( vMIP, vphi )
        vphi, vMIP = Vphi.PHI(x1, totalpart=True, perturb='STATES')
        vphi = round(vphi,5)
        print "virgil (states): \tMIP=%s \t PHI(x)=%s" % ( vMIP, vphi )
#        print "virgil ei(x)=%s \t prob=%s" % (Vphi.ei(x1), Vphi.x1_STATS[x1]['prob_this_x1'] )
        sys.stdout.flush()
        
        print ""
        sys.stdout.flush()
#    '''        
#    print "average_ei=%s" % Vphi.average_ei()


from copy import deepcopy
from pprint import pprint

def BELL_add_entry( myinput, c, w ):
    '''appends 'w' to list with index 'c' ''' 
    if len(myinput) <= c:
        myinput.append( [] )

    (myinput[c]).append( w )

    return myinput


def BELL_all_partitions( result, cue ):
    '''this function returns the bell numbers.'''

    if cue:
        j = cue.pop(0)
        for i in range( len(result)+1 ):
             BELL_all_partitions( BELL_add_entry(deepcopy(result),i,j), deepcopy(cue) )
    else:
#        print "SHOULD YIELD THIS: result=%s" % result
        global BELL_biglist
        #    print "\t\t f: result=%s \t cue=%s" % (result, cue)
        BELL_biglist.append( result )

#        yield result

            
if __name__ == '__main__':
    PROFILE_CODE = False

    if PROFILE_CODE:
        import profile, pstats
        p = profile.Profile()
        p.run("main()")

        s = pstats.Stats(p)
        s.sort_stats("time", "name").print_stats()

    else:
        # run the tests        
        _test()
#        main()
