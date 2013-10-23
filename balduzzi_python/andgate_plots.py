#!/usr/bin/python
################################################################################################################################################
# generic modules
import numpy as N
import matplotlib.pylab as p
from numpy import array, matrix
# PHI-specific modules
import phi_balduzzi as Bphi
import phi as Virgilphi
import cxnmatrices, sys, TerminalController
#import phiplots
from pprint import pprint
from math import ceil


term = TerminalController.TerminalController()

################################################################################################################################################
# This code calculates the PHI of a perceptron
# 1. Plot PHI of 0 and 1 and average phi with BIAS=-1.0 with increasning number of nodes.
################################################################################################################################################
# What is the bias of the output node?

#NUM_INPUT_NODES = [2,3,4,5]
NUM_INPUT_NODES = [2,4]
#NUM_INPUT_NODES = [5]
#NUM_INPUT_NODES = [2]
#NUM_INPUT_NODES = [1,2,3,4,5]

# If True, plot the PHI using the Virgil code.  If False, plot the PHI using the Balduzzi code
USE_VIRGIL_PHI = False
label='binary'

# Size of the points
PROFILE_CODE = False
MARKER_SIZE = 5.0
################################################################################################################################################

def output_nodefunction( suminput ):
    return array(suminput >= 2, dtype=N.int8)

#Bphi.WHICH = output_nodefunction
################################################################################################################################################


################################################################################################################################################
# Initialize the lists that we are going to fill and plot
#zero_phis, one_phis, expected_phis = [], [], []
#zero_eis, one_eis, expected_eis = [], [], []
#MIPs = []

################################################################################################################################################    
# PART 1 -- Get the values of PHIs and EIs with increasing number of nodes.
# If SET_BIAS_TO_NUMINPUT_UNITS is True, then  OUTPUTNODE_BIAS is OVERWRITTEN.  Else is kept constant
################################################################################################################################################    

expected_phis, expected_eis = [], []

for num_input_nodes in NUM_INPUT_NODES:
    print term.YELLOW + "\nRuning andgate with %s input nodes..." % num_input_nodes + term.NORMAL

    # set the filename_root and cxnmatrix for this number of input nodes

    fn_root='/PHIplots/andgateplots2/ag-N=%s' % num_input_nodes    
    cxnmatrix = cxnmatrices.andgates[num_input_nodes]
#    fn_root='backcxns/ag-N=%s' % num_input_nodes    
#    cxnmatrix = cxnmatrices.andgates_with_backcns[num_input_nodes]
    numunits=len(cxnmatrix)
    print "cxnmatrix=%s" % cxnmatrix

    # Initialize the statistics we want to compute
    phis, eis, MIPs = [None]*2**numunits, [None]*2**numunits, []
    expected_value_of_phi, expected_value_of_ei = 0.0, 0.0
    
    ################################################################################################################################################    
    Vphi = Virgilphi.phi(N=numunits, cxnmatrix=cxnmatrix, nodefunc=output_nodefunction)
    ################################################################################################################################################    

    print term.WHITE + "Valid Outputs=%s" % Vphi.VALID_OUTPUTS + term.NORMAL

    
    print "===AVERAGES=========================================="
    print "average_ei=\t%s" % Vphi.average_ei( array([array(range(Vphi.numunits))]) )
    averageMIP = Vphi.average_MIP()
    print "average_MIP=\t%s \t average_ei(average_MIP)=%s" % (str(averageMIP), Vphi.average_ei(averageMIP))
    print "average PHI (raw)=\t%s" % Vphi.average_phi()
    print "==END=AVERAGES======================================="

    print "===PERSTATE=========================================="
    for x1 in Vphi.VALID_OUTPUTS:
        x1state = Vphi._dec2binstate( x1 )
        print "x1=%s \t x1state=%s" % (x1, x1state)
#        print term.CYAN + "possible x0s: %s" % Vphi.x1_STATS[x1]['possible_x0s'] + term.NORMAL
        
        possible_x0states = map( Vphi._dec2binstate, Vphi.x1_STATS[x1]['possible_x0s'] )
        print "number of x0s=%s" % len(Vphi.x1_STATS[x1]['possible_x0s'])
#        pprint( possible_x0states )
        
        bphi, partB, partA = Bphi.phi(cxnmatrix,range(numunits),x1state,T=-1,iself=True)
#        bphi, bMIP, bei = (round(bphi, 5), (partB,partA), round(Bphi.tcut( cxnmatrix, range(numunits), x1state ),5) )
#        print "\t virgil ei(x)=%s \t prob=%s" % (Vphi.ei(x1), Vphi.x1_STATS[x1]['prob_this_x1'] )
        print "\t balduzzi (wires): \tMIP=%s \t PHI(x)=%s" % ( bMIP, bphi )
        sys.stdout.flush()

        if USE_VIRGIL_PHI:

            # get the phi(x1) and MIP using wires
#            vphi, vMIP = Vphi.PHI(x1, perturb='wires')
#            vphi = round(vphi,5)
#            print "\t virgil (wires): \tMIP=%s \t PHI(x)=%s" % ( vMIP, vphi )
            vphi, vMIP = Vphi.PHI(x1, perturb='states', allpartitions=False )
            vphi, vei = round(vphi,5), round(Vphi.ei(x1),5)
            print "\t virgil (states): \tMIP=%s \t PHI(x)=%s" % ( vMIP, vphi )
            sys.stdout.flush()

#            assert vphi == bphi, "vphi != bphi!  vphi=%s  bphi=%s" % (vphi, bphi)
#            assert vei == bei, "vei != bei!  vei=%s  bei=%s" % (vei, bei)
            phi, MIP, ei = vphi, vMIP, vei

        # Using the Balduzzi PHI
        else:
            phi, MIP, ei = bphi, bMIP, bei

        phis[x1] = phi
        eis[x1] = ei

        # What is this state's contribution to the expected value of phi
        expected_value_of_phi += Vphi.x1_STATS[x1]['prob_this_x1'] * phi
        expected_value_of_ei += Vphi.x1_STATS[x1]['prob_this_x1'] * ei

        MIPs.append( (x1,MIP) )
        print "---------------------------"

    print "==END=PERSTATE======================================="
    # Calculate the average_MIP
#    average_MIP__wires = Vphi.average_MIP(perturb='wires')    
#    print term.WHITE + "\t virgil (wires):  \t average MIP=%s" % str(average_MIP__wires) + term.NORMAL
#    average_MIP__states = Vphi.average_MIP(perturb='states', allpartitions=False)
#    print term.WHITE + "\t virgil (states): \t average MIP=%s" % str(average_MIP__states) + term.NORMAL
        
    expected_value_of_ei, expected_value_of_phi = round(expected_value_of_ei,4), round(expected_value_of_phi,4)

    expected_eis.append( (num_input_nodes,expected_value_of_ei) )
    expected_phis.append( (num_input_nodes,expected_value_of_phi) )

    # switch the labels if we're plotting a lot of nodes
    if num_input_nodes >= 5:
        label = 'decimal'
    else:
        label = 'binary'
    
#    phiplots.plotmatrix( cxnmatrix, title='cxn matix of '+fn_root, filename=fn_root+'M.png' )
#    phiplots.plot_cmatrix( cxnmatrix, fn_root+'G.png', caption='cxn-graph of '+ fn_root, prog='dot' )
#    phiplots.plot_tmatrix( Vphi.tmatrix, fn_root+'T.png', caption='t-graph of ' + fn_root )

#    phiplots.states_and_hist( eis, fn_root+'EI.png', show_invalidstates=False, label=label, plotting_eis=True, title='$ei(X_1)$ \t $E[ei(X_1)]=%s$' % expected_value_of_ei )
#    phiplots.states_and_hist( phis, fn_root+'PHIs.png', show_invalidstates=False, label=label, title='$\phi(X_1)$ \t $E[\phi(X_1)]=%s$' % expected_value_of_phi )
#    phiplots.plot_MIPs_for_allstates( MIPs, title=fn_root+' parts', filename=fn_root+'P.png', label=label, averageMIP=('ei(P)',average_MIP__states) )

####################################################################################
# Now to plot the expected phi and ei
####################################################################################
'''
print "zero_phis: "
pprint( zero_phis )
print "one_phis: "
pprint( one_phis )
print "expected_phis:"
pprint( expected_phis )
####################################################################################
'''
p.clf()            

p.xlabel('Number of input units')
p.ylabel('$\phi(X_1)$ and $ei(X_1)$')
#        p.xlim(0, maxlength)
#X1s, Y1s = range(len(zero_phis)), zero_phis
#    X2s, Y2s = range(len(one_phis)), one_phis        
X0s, Y0s = [ x for x, y in expected_phis ], [ y for x, y in expected_phis ]
X1s, Y1s = [ x for x, y in expected_eis ], [ y for x, y in expected_eis ]

#line1 = p.plot( X1s, Y1s )
#line2 = p.plot( X2s, Y2s)
#line0 = p.plot( X0s, Y0s, color='red', marker='o', label='$\phi(0)$', markersize=MARKER_SIZE )
#line1 = p.plot( X1s, Y1s, color='green', marker='o', label='$\phi(1)$', markersize=MARKER_SIZE )
line0 = p.plot( X0s, Y0s, color='red', marker='s', label='$E[\phi(X_1)]$', linewidth=2.5 )
line1 = p.plot( X1s, Y1s, color='green', marker='s', label='$E[ei(X_1)]$', linewidth=2.5 )

#p.legend((line0,line1), ('0 state', '1 state'), shadow = True )        
#p.legend(loc='best', shadow=True)
p.legend( loc='upper left', shadow=True)

p.title('Average phi and ei of andgates')

# Set the Xticks
xticks = list(set(X0s + X1s))
xticks.sort()
p.xticks( xticks )
if Y0s + Y1s:
    p.axis(ymin=-0.1, ymax=ceil(max(Y0s+Y1s))+0.5)


p.savefig( 'andgate-averages.png', dpi=300 )
p.clf()

