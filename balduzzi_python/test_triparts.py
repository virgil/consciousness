#!/usr/bin/python
################################################################################################################################################
# generic modules
import numpy as N
import matplotlib.pylab as p
from numpy import array, matrix

# PHI-specific modules
import phi as Virgilphi
import phiplots, cxnmatrices, sys, TerminalController
from pprint import pprint
from math import ceil

term = TerminalController.TerminalController()

################################################################################################################################################
# This code calculates the PHI of a perceptron
# 1. Plot PHI of 0 and 1 and average phi with BIAS=-1.0 with increasning number of nodes.
################################################################################################################################################
# What is the bias of the output node?

NUM_INPUT_NODES = [2,3,4,5]
#NUM_INPUT_NODES = [5]
#NUM_INPUT_NODES = [2]
#NUM_INPUT_NODES = [1,2,3,4,5]

# If True, plot the PHI using the Virgil code.  If False, plot the PHI using the Balduzzi code
USE_VIRGIL_PHI = True
label='binary'

################################################################################################################################################

def output_nodefunction( suminput ):
    return array(suminput >= 2, dtype=N.int8)

################################################################################################################################################    
# PART 1 -- Get the values of PHIs and EIs with increasing number of nodes.
# If SET_BIAS_TO_NUMINPUT_UNITS is True, then  OUTPUTNODE_BIAS is OVERWRITTEN.  Else is kept constant
################################################################################################################################################    

print term.YELLOW + "\nRunning tripartition test..." + term.NORMAL

    # set the filename_root and cxnmatrix for this number of input nodes

fn_root='/PHIplots/test_triparts'
cxnmatrix = cxnmatrices.disc3ands
#    fn_root='backcxns/ag-N=%s' % num_input_nodes    
#    cxnmatrix = cxnmatrices.andgates_with_backcns[num_input_nodes]
numunits=len(cxnmatrix)
print "cxnmatrix=%s" % cxnmatrix

################################################################################################################################################    
Vphi = Virgilphi.phi(N=numunits, cxnmatrix=cxnmatrix, nodefunc=output_nodefunction)
################################################################################################################################################    

print term.WHITE + "Valid Outputs=%s" % Vphi.VALID_OUTPUTS + term.NORMAL

MIPs, eis, phis = [], [], []
    

for x1 in Vphi.VALID_OUTPUTS:
    x1state = Vphi._dec2binstate( x1 )
    print "x1=%s \t x1state=%s" % (x1, x1state)
    print term.CYAN + "possible x0s: %s" % Vphi.x1_STATS[x1]['possible_x0s'] + term.NORMAL
    
    possible_x0states = map( Vphi._dec2binstate, Vphi.x1_STATS[x1]['possible_x0s'] )
    print "number of x0s=%s" % len(Vphi.x1_STATS[x1]['possible_x0s'])
#        pprint( possible_x0states )
    
    print "\t virgil ei(x)=%s \t prob=%s" % (Vphi.ei(x1), Vphi.x1_STATS[x1]['prob_this_x1'] )
    print "---------------------------"
    sys.stdout.flush()

    if USE_VIRGIL_PHI:
        vphi, vMIP = Vphi.PHI(x1, totalpart=True, perturb='states', allpartitions=True )
        vphi, vei = round(vphi,5), round(Vphi.ei(x1),5)
        print "\t virgil (states): \tMIP=%s \t PHI(x)=%s" % ( vMIP, vphi )
        sys.stdout.flush()

        phi, MIP, ei = vphi, vMIP, vei

    phis[x1] = phi
    eis[x1] = ei

    MIPs.append( MIP )
    
phiplots.plotmatrix( cxnmatrix, title='cxn matix of '+fn_root, filename=fn_root+'M.png' )
phiplots.plot_cmatrix( cxnmatrix, fn_root+'G.png', caption='cxn-graph of '+ fn_root, prog='dot' )
phiplots.plot_tmatrix( Vphi.tmatrix, fn_root+'T.png', caption='t-graph of ' + fn_root )

phiplots.states_and_hist( eis, fn_root+'EI.png', show_invalidstates=False, label=label, plotting_eis=True, title='$ei(X_1)$ \t $E[ei(X_1)]=%s$' % expected_value_of_ei )
phiplots.states_and_hist( phis, fn_root+'PHIs.png', show_invalidstates=False, label=label, title='$\phi(X_1)$ \t $E[\phi(X_1)]=%s$' % expected_value_of_phi )
phiplots.plot_MIPs_for_allstates( MIPs, title=fn_root+' parts', filename=fn_root+'P.png', show_invalidstates=False, label=label)

