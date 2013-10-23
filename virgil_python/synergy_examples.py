#!/usr/bin/epd-python
# -*- coding: utf-8 -*-
# #!/usr/bin/python

from pprint import pprint
from numpy import log2
import numpy
from copy import deepcopy
from Xjointstate import Xjointstate
import sys
import TerminalController as TC
term = TC.TerminalController()
import scipy.optimize

#import Pstar


import tables
from tables import SUITE_N2

#SUITE_N2 = [ tables.AND, tables.X1 ]
#SUITE_N2 = [ tables.X1 ]

#SUITE_N2 = [ ]

PRECISION = 5
eps = 10**(-PRECISION)
epps = 10 * eps

#TABLE = tables.AND_skew2
#TABLE = tables.INDEPENDENCE_TWO_ab
#TABLE = tables.CHRIS_ELLISON
#TABLE = tables.SNAGGER
#TABLE = tables.SNAGGER3
#TABLE = tables.AND_NONDET

#TABLE = tables.INDEPENDENCE_TWO_ab
#TABLE = tables.SNAGGER3

#TABLE = tables.MAUER03

###############################
# N=2
###############################
#TABLE = tables.INDEPENDENCE_TWO
#TABLE = tables.INDEPENDENCE_TWO_ab
#TABLE = tables.XOR
#TABLE = tables.REDUNDANCY_TWO
#TABLE = tables.ANY_DOUBLET
#TABLE = tables.X0
#TABLE = tables.X1
#TABLE = tables.OR
#TABLE = tables.AND_skew
#TABLE = tables.AND_skew2
#TABLE = tables.AND_skew_notuniform
#TABLE = tables.AND_X0           # no Markov Model passed all 5 tests ):
#TABLE = tables.AND_X1
#TABLE = tables.AND_OR
#TABLE = tables.AND_XOR
#TABLE = tables.XOR_X1
#TABLE = tables.SYN_REDUN_INDEP_TWO
#TABLE = tables.SYN_EQUALS_REDUN_TWO
#TABLE = tables.XORx1x2rdn
#TABLE = tables.LATHAM_FIG4
#TABLE = tables.AND_mod
#TABLE = tables.AND_DUP
#TABLE = tables.XOR_z_x_xor_y
#TABLE = tables.TriAND
#TABLE = tables.XOR_DUPLICATE

#TABLE = tables.XOR_TEST 
#TABLE = tables.ONE_TO_ONE
#TABLE = tables.ONE_TO_MANY
#TABLE = tables.TRIPLE_ADD
#TABLE = tables.DOUBLE_XOR

#TABLE = tables.ADD_CIRCUIT

#TABLE = tables.XOR_UNIQUE

# Two nice holism examples
## #TABLE = tables.HOLISM_AB_or_C
## #TABLE = tables.HOLISM_AB_xor_C

#TABLE = tables.AND
#TABLE = tables.ANY_DOUBLET
#TABLE = tables.PARITY_REDUN_REDUN
TABLE = tables.NEG_SYN_TEST2






def p_x_y_z1z2( xstate, yvalue, z1state, z2state ):
    '''computes p(x,y,z1,z2) = p(x,y) * p(z1|xy) * p(z2|xy)'''

    #assert get_numXs() == 2
    #check_xs( xstate.indices, xstate.states )
    #check_y_value( yvalue )
    check_xs( z1state.indices, z1state.values )
    check_xs( z2state.indices, z2state.values )

    returnval = Prob_Xs_y( xstate.indices, xstate.values, yvalue )
    returnval *= Prob_xs_given_condxs_y( z1state.indices, z1state.values, xstate.indices, xstate.values, yvalue )
    returnval *= Prob_xs_given_condxs_y( z2state.indices, z2state.values, xstate.indices, xstate.values, yvalue )

    assert 0.0 <= returnval <= 1.0

    return returnval



def Iintersect_lowerbound_v1():
    '''This uses the state-dependent WholeMinusSum extracting any negative values to get a lowerbound on the average redundancy.'''

    assert get_numXs() == 2
    X1, X2, X1X2 = [0], [1], [0,1]

    z = 0.0

    for y in Y_states():
        prob_y = Prob_Y_equals_y(y)

        whole = specific_info_with_Y_equals_y(X1X2, y )
        X1term = specific_info_with_Y_equals_y(X1, y )
        X2term = specific_info_with_Y_equals_y(X2, y )

        diff = whole - X1term - X2term
        print "\t y=%s \t diff=%s \t z=%s" % (y,diff, z)
        # if there is REDUNDANCY for this state, add it to z
        if diff < 0.0:
            z += prob_y * abs(diff)

    # get the average WholeMinusSum, and check that it is closer to zero than z
    average_WMS = mutual_info_with_Y(X1X2) - mutual_info_with_Y(X1) - mutual_info_with_Y(X2)

    if average_WMS < 0.0:
        assert abs(average_WMS) <= (z+eps), "average_WMS=%s \t z=%s" % (average_WMS,z)

    assert 0.0 <= z

    assert z <= mutual_info_with_Y(X1)+eps
    assert z <= mutual_info_with_Y(X2)+eps
    assert z <= mutual_info_between_Xs(X1,X2)+eps, "z=%s I(X1:X2)=%s" % (z,mutual_info_between_Xs(X1,X2))


    return z






def reduce_Xstates( c1, c2, incoming_table ):
    '''returns the REDUCED version of incoming_table collapsing states 'x1' and 'x2'.'''

    assert c1.indices == c2.indices, "Indices should be the same"
    assert c1.states != c2.states, "But the states should be different"

    # the collapse index should be a single X_index
    collapse_index = c1.indices[0]
    c1state, c2state = c1.states[0], c2.states[0]


    z = { }

    for entry_xs, entry_y in incoming_table.iteritems():

        # split the Xs by a space
        new_entry_xs = entry_xs.split()
        new_entry_y = entry_y
        ########## PART ONE ##########
        # if the collapse index isn't either of our c1 or c2 states, add it to the output and goto the next.
        if new_entry_xs[collapse_index] not in [c1state, c2state]:
            z[entry_xs] = entry_y
            continue


        ########## PART TWO ##########
        # assert the new_entry_xs is c1state or c2state
        assert new_entry_xs[collapse_index] in [c1state, c2state]

        # set the state the COLLAPSED (joined) version
        new_entry_xs[collapse_index] = '%s%s' % (c1state,c2state)
        assert new_entry_xs[collapse_index] == str(c1state + c2state)

        # and rejoin the new_entry_xs
        new_entry_xs = ' '.join(new_entry_xs)

        # if this new_entry_xs isn't in z, add it along with the y, and goto the next.
        if new_entry_xs not in z:
            z[new_entry_xs] = new_entry_y
            continue

        ########## PART THREE ##########
        # assert that new_entry_xs is already in the output table
        assert new_entry_xs in z

        old_y = z[new_entry_xs]

        # if either y is a str, make it a list
        if type(new_entry_y) is str:
            new_entry_y = [new_entry_y]

        if type(old_y) is str:
            old_y = [old_y]

        # set the output table to the concatenation of new_entry_y and old_y
        z[new_entry_xs] = old_y + new_entry_y

    return z

def Prob_caret_x1_x2_y( x1, x2, y_value ):
    '''computes p( x1 ^ x2, y ).   p( x1 ^ x2, y ) = p(x1,x2, y) if and only if p(x1) = p(x2) = p(x1,x2).  Otherwise it equals zero.'''

    check_xs( x1.indices, x1.states )
    check_xs( x2.indices, x2.states )
    check_y_value( y_value )

    x1x2 = Xjointstate( x1.indices+x2.indices, x1.states+x2.states )
    z = Prob_caret_x1_x2(x1, x2)

    if z:
        return Prob_Xs_y(x1x2.indices, x1x2.states, y_value )

    return 0.0

def Prob_caret_x1_x2( x1, x2 ):
    '''computes p( x1 ^ x2 ).   p( x1 ^ x2 ) = p(x1,x2) if and only if p(x1) = p(x2) = p(x1,x2).  Otherwise it equals zero.'''

    check_xs( x1.indices, x1.states )
    check_xs( x2.indices, x2.states )

    prob_x1 = Prob_Xs( x1.indices, x1.states )
    prob_x2 = Prob_Xs( x2.indices, x2.states )
    prob_x1x2 = Prob_Xs( x1.indices+x2.indices, x1.states+x2.states )

    if prob_x1 == prob_x2 == prob_x1x2:
        return prob_x1x2

    return 0.0




def PRINT_SYNERGY_MEASURES( Xs ):
    '''prints out a bunch of synergy measures about Xs about target Y '''
    check_Xs( Xs )
    assert Xs == allXs

    print "===================================================================="
    print "Measures \t\t bits \t\t  normalized"
    print "===================================================================="
    print "WholeMinusSum(Xs:Y) \t\t %s \t\t\t %s" % ( Chechik(Xs), Chechik(Xs,norm=True) )
    print "WMMS(Xs:Y) \t\t %s \t\t\t %s" % ( Synmin(Xs), Synmin(Xs,norm=True) )
    print "I.I.(Xs;Y) \t\t %s \t\t\t %s" % ( Interaction_Information(Xs), Interaction_Information(Xs,norm=True) )
    print "I.I.v2(Xs;Y) \t\t %s" % ( Interaction_Information_v2(Xs) )
#    print "I.I.v3(Xs;Y) \t\t %s" % ( Interaction_Information_v3(Xs) )
#    print "I.I.v4(Xs;Y) \t\t %s" % ( Interaction_Information_v4(Xs) )
    print "--------------------------------------------------------------------"
    print "DeltaI(Xs;Y) \t\t %s \t\t\t %s" % ( Delta_I(Xs), None )
    print "S_max(Xs:Y) \t\t %s \t\t\t %s" % ( Imax_Synergy(Xs), Imax_Synergy(Xs,norm=True) )
    print "H_max(Xs:Y) \t\t %s \t\t\t %s" % ( Holism(Xs), Holism(Xs,norm=True) )
    print "CoR(Xs:Y) \t\t %s" % ( Cost_of_Removal(Xs) )
    print "===================================================================="


    print_redundancy_atoms( allXs )

    #    gacs_korner = common_random_variable( [0], [1], True )
    #    print "Gács-Körner=%s" % gacs_korner


    print "Imin=%s" % round(Imin_coalition_Y( [[0],[1]]),PRECISION)
    print "Ismw=%s" % round(Ismw_Y(), PRECISION)
    print "I(X0:X1)=%s" % str(round(H_Xs([0]) + H_Xs([1]) - H_Xs([0,1]),PRECISION))





def run_gce_tests2():
    '''run all of the GCE tests on the given markov model -- with the ORDER of models/suites FLIPPED'''
    global TABLE
    old_TABLE = TABLE

    print "---------GCE--TESTS-----------------------------"

    # for each system in the n=2 suite...
    for system in SUITE_N2:
        print term.YELLOW + "\rSYSTEM:\t",
        pprint( system )
        print term.NORMAL,

        TABLE = system

        gce1 = generalized_conditional_entropy_test1()
        print "------------------------------------------------"
        gce2 = generalized_conditional_entropy_test2()
        print "------------------------------------------------"
        gce3 = generalized_conditional_entropy_test3_union()
        #       gce3 = generalized_conditional_entropy_test3_intersect()
        #       print "------------------------------------------------"
        gce4 = generalized_conditional_entropy_test4_union()
        #       gce4 = generalized_conditional_entropy_test4_intersect()
        print ''

        if not gce1 or not gce2 or not gce3 or not gce4:
            passed_all_tests = False
            break

        #            raw_input('...continue...')

    # set it back to None
    TABLE = old_TABLE


def generalized_conditional_entropy_test1():
    '''A test to see if H(Y|{X1,X1}=H(Y|X1)'''

    global MARKOV_MODEL

    if get_numXs() != 2:
        return

    cond1, cond2 = H_Y_given_Xs([0]), H_Y_given_Xs([1])
    #cond11, cond112 = H_Y_given_Xs([0,0]), H_Y_given_Xs([0,0,1])
    #print "cond11=%s \t cond112=%s" % (cond11,cond112)

    #cond_1, cond_2 = H_Y_given_multiXs( [[0]] ), H_Y_given_multiXs( [[1]] )
    cond_1_1 = H_Y_given_multiXs( [[0],[0]] )
    cond_2_2 = H_Y_given_multiXs( [[1],[1]] )

    passed_test = True

    print "gce1: H(Y)=%.4f" % H_Y()
    if abs(cond1-cond_1_1) > epps:
        print term.RED,
        passed_test = False
    #    print "\rgce1: H(Y|X1)=%.4f \t H(Y|{X1,X1})=%.4f" % ( cond1, cond_1_1 ) + term.NORMAL
    print "\rgce1: H(Y|X1)=%.4f \t H(Y|X1 U X1)=%.4f" % ( cond1, cond_1_1 ) + term.NORMAL

    if abs(cond2-cond_2_2) > epps:
        print term.RED,
        passed_test = False
    #    print "\rgce1: H(Y|X2)=%.4f \t H(Y|{X2,X2})=%.4f" % ( cond2, cond_2_2 ) + term.NORMAL
    print "\rgce1: H(Y|X2)=%.4f \t H(Y|X2 U X2)=%.4f" % ( cond2, cond_2_2 ) + term.NORMAL

    ############################################################################
    # ensure that p(y|x1) == p(y|x1,x1) == p(y|x1 ∪ x1)  \forall y, x1 states
    # ensure that p(y|x2) == p(y|x2,x2) == p(y|x2 ∪ x2)  \forall y, x2 states
    ############################################################################

    # for every state of y...
    for y_state in Y_states():
        p_y = Prob_Y_equals_y( y_state )

        # for every x1 state...
        for x1 in joint_states([0], True):
            x1x1 = Xjointstate( x1.indices + x1.indices, x1.states + x1.states )
            t1 = Prob_y_given_xs( y_state, x1.indices, x1.states )
            t1_1 = pstar( MARKOV_MODEL, y=y_state, GIVENx1=x1, GIVENx2=x1 )
            t11 = Prob_y_given_xs( y_state, x1x1.indices, x1x1.states )

            assert (t1-t11) <= eps, "Probabilities weren't all equal!"

            if abs(t1-t1_1) > epps:
                print term.RED,
                print u"\ry=%s x1=%s \t p(y)=%.4f \t p(y|x1x1)=%.4f \t p(y|x1 U x1)=%.4f \t p(y|x1)=%.4f" % ( y_state, x1.states[0], p_y, t11, t1_1, t1 )
                sys.stdout.flush()
                print term.NORMAL,
                if passed_test:
                    raw_input('....failed for this but passed_test = TRUE')


        for x2 in joint_states([1], True):
            x2x2 = Xjointstate( x2.indices + x2.indices, x2.states + x2.states )
            t2 = Prob_y_given_xs( y_state, x2.indices, x2.states )
            t2_2 = pstar( MARKOV_MODEL, y=y_state, GIVENx1=x2, GIVENx2=x2 )
            t22 = Prob_y_given_xs( y_state, x2x2.indices, x2x2.states )

            assert abs(t2 - t22) <= eps, "Probabilities weren't all equal!"

            if abs(t2-t2_2) > epps:
                print term.RED,
                print u"\ry=%s x2=%s \t p(y)=%.4f \t p(y|x2x2)=%.4f \t p(y|x2 U x2)=%.4f \t p(y|x2)=%.4f" % ( y_state, x2.states[0], p_y, t22, t2_2, t2 )
                sys.stdout.flush()
                print term.NORMAL,
                if passed_test:
                    raw_input('....failed for this but passed_test = TRUE')



                #    assert abs(cond1-cond_1_1) <= epps
                #    assert abs(cond2-cond_2_2) <= epps

    return passed_test


def generalized_conditional_entropy_test2():
    '''Test GCE symmetry--that H(Y|{X1,X2}) == H(Y|{X2,X1}'''

    cond_1_2 = H_Y_given_multiXs([[0],[1]])
    cond_2_1 = H_Y_given_multiXs([[1],[0]])

    passed_test = True

    if abs(cond_1_2-cond_2_1) > epps:
        print term.RED,
        passed_test = False
    #    print "\rgce2: H(Y|{X1,X2})=%.4f == H(Y|{X2,X1})=%.4f" % ( cond_1_2, cond_2_1 ) + term.NORMAL
    print "\rgce2: H(Y|X1 U X2)=%.4f == H(Y|X2 U X1)=%.4f" % ( cond_1_2, cond_2_1 ) + term.NORMAL

    #    assert (cond_1_2 - cond_2_1) <= epps

    return passed_test


def generalized_conditional_entropy_test3_intersect():
    '''Test the bounds of GCE---that \max_i H(Y|X_i) <= H(Y|{X1,X2}) <=  H(Y)'''

    upper = H_Y()
    cond_1_2 = H_Y_given_multiXs([[0],[1]])
    max_cond1_cond2 = max( H_Y_given_Xs([0]), H_Y_given_Xs([1]) )

    passed_test = True

    if not ((max_cond1_cond2-epps) <= cond_1_2):
        print term.RED,
        passed_test = False
    print "\rgce3: max=%.4f <= H(Y|{X1,X2})=%.4f" % ( max_cond1_cond2, cond_1_2 ) + term.NORMAL
    #    print "\rgce3: max=%.4f <= H(Y|X1 U X2)=%.4f" % ( max_cond1_cond2, cond_1_2 ) + term.NORMAL

    if not (cond_1_2 <= (upper+epps)):
        print term.RED,
        passed_test = False
    print "\rgce3: H(Y|{X1,X2})=%.4f <= H(Y)=%.4f" % ( cond_1_2, upper ) + term.NORMAL
    #    print "\rgce3: H(Y|X1 U X2)=%.4f <= H(Y)=%.4f" % ( cond_1_2, upper ) + term.NORMAL


    #    assert 0.0 <= cond12 <= (cond_1_2+epps)
    #    assert cond_1_2 <= (min_cond1_cond2+epps)

    return passed_test

def generalized_conditional_entropy_test3_union():
    '''Test the bounds of GCE---that H(Y|X1,X2) <= H(Y|{X1,X2}) <=  min[ H(Y|X1), H(Y|X2) ]'''

    cond12 = H_Y_given_Xs([0,1])
    cond_1_2 = H_Y_given_multiXs([[0],[1]])
    min_cond1_cond2 = min( H_Y_given_Xs([0]), H_Y_given_Xs([1]) )

    passed_test = True

    if not ((cond12-epps) <= cond_1_2):
        print term.RED,
        passed_test = False
    #    print "\rgce3: H(Y|X12)=%.4f <= H(Y|{X1,X2})=%.4f" % ( cond12, cond_1_2 ) + term.NORMAL
    print "\rgce3: H(Y|X12)=%.4f <= H(Y|X1 U X2)=%.4f" % ( cond12, cond_1_2 ) + term.NORMAL

    if not (cond_1_2 <= (min_cond1_cond2+epps)):
        print term.RED,
        passed_test = False
    #    print "\rgce3: H(Y|{X1,X2})=%.4f <= min=%.4f" % ( cond_1_2, min_cond1_cond2 ) + term.NORMAL
    print "\rgce3: H(Y|X1 U X2)=%.4f <= min[H(Y|X1),H(Y|X2)]=%.4f" % ( cond_1_2, min_cond1_cond2 ) + term.NORMAL


    #    assert 0.0 <= cond12 <= (cond_1_2+epps)
    #    assert cond_1_2 <= (min_cond1_cond2+epps)

    return passed_test

def generalized_conditional_entropy_test4_intersect():
    '''Test the redundancy ({1}{2}) gotten when using this GCE against the known bounds for redundancy/synergy.'''

    # if this isn't true we can't use these bounds
    assert H_Y_given_Xs([0,1]) == 0.0

    # get the difference
    mi1, mi2 = mutual_info_with_Y([0]), mutual_info_with_Y([1])
    D = mutual_info_with_Y([0,1]) - mi1 - mi2

    cond_1_2 = H_Y_given_multiXs([[0],[1]])
    gce_rdn = H_Y() - cond_1_2
    gce_syn = D + gce_rdn

    min_rdn = max(0.0,-D)

    max_rdn = min(mi1,mi2)

    min_syn = max(0.0,D)
    max_syn = D + min(mi1,mi2)

    passed_test = True

    print "\rgce4: D=%.4f \t min(mi1,mi2)=%.4f" % (D, min(mi1,mi2))

    # if gce_rdn is not within bounds, print in RED
    if not ( (min_rdn-epps) <= gce_rdn <= (max_rdn+epps) ):
        print term.RED,
        passed_test = False
    print "\rgce4: min_rdn=%.4f <= %.4f <= %.4f" % (min_rdn, gce_rdn, max_rdn) + term.NORMAL

    # if gce_syn is not within bounds, print in RED
    if not ( (min_syn-epps) <= gce_syn <= (max_syn+epps) ):
        print term.RED,
        passed_test = False
    print "\rgce4: min_syn=%.4f <= %.4f <= %.4f" % (min_syn, gce_syn, max_syn) + term.NORMAL

    return passed_test


def generalized_conditional_entropy_test4_union():
    '''Test the synergy ({12}) gotten when using this GCE against the known bounds for synergy.'''

    # if this isn't true we can't use these bounds
    assert H_Y_given_Xs([0,1]) == 0.0

    # get the difference
    mi1, mi2 = mutual_info_with_Y([0]), mutual_info_with_Y([1])
    D = mutual_info_with_Y([0,1]) - mi1 - mi2

    cond_1_2 = H_Y_given_multiXs([[0],[1]])
    gce_syn = cond_1_2
    gce_rdn = gce_syn - D

    min_rdn = max(0.0,-D)

    max_rdn = min(mi1,mi2)

    min_syn = max(0.0,D)
    max_syn = D + min(mi1,mi2)

    passed_test = True

    print "\rgce4: D=%.4f \t min(mi1,mi2)=%.4f" % (D, min(mi1,mi2))

    # if gce_rdn is not within bounds, print in RED
    if not ( (min_rdn-epps) <= gce_rdn <= (max_rdn+epps) ):
        print term.RED,
        passed_test = False
    print "\rgce4: min_rdn=%.4f <= %.4f <= %.4f" % (min_rdn, gce_rdn, max_rdn) + term.NORMAL

    # if gce_syn is not within bounds, print in RED
    if not ( (min_syn-epps) <= gce_syn <= (max_syn+epps) ):
        print term.RED,
        passed_test = False
    print "\rgce4: min_syn=%.4f <= %.4f <= %.4f" % (min_syn, gce_syn, max_syn) + term.NORMAL

    return passed_test


def H_Y_given_Xs( Xs ):
    '''returns H[Y | Xs ]'''

    check_Xs( Xs )

    X_states = joint_states( Xs, True )

    z = 0.0
    for x_state in X_states:
        prob_x = Prob_Xs( x_state.indices, x_state.values )

        summ = 0.0
        for y_state in Y_states():

            prob_y_given_x = Prob_y_given_xs( y_state, x_state.indices, x_state.values )

            if prob_y_given_x == 0.0:
                continue

            summ += prob_y_given_x * log2( 1.0 / prob_y_given_x )

        z += prob_x * summ


    assert 0.0 <= z
    assert z <= (H_Y() + epps)

    # Now evaluate H(Y|Xs) = \sum_{y} p(y) * \sum_{x1} p(x1|y) \log[ 1 / p(y|x1) ]
    z2 = 0.0
    for y_state in Y_states():
        prob_y = Prob_Y_equals_y( y_state )
        if prob_y:
            summ = sum([ Prob_Xs_given_y( x1.indices, x1.values, y_state ) * log2(1.0/Prob_y_given_xs( y_state, x1.indices, x1.values )) for x1 in X_states if Prob_Xs_given_y( x1.indices, x1.values, y_state ) ])
            z2 += prob_y * summ

    assert 0.0 <= z2 <= (H_Y() + epps)
    assert abs(z-z2) <= epps



    # Now evaluate H(Y|Xs) = \sum_{x1,y} p(x1,y) \log[ 1 / p(y|x1) ]
    z3 = 0.0
    for x1 in X_states:
        for y in Y_states():
            prob_x1_y = Prob_Xs_y(x1.indices, x1.values, y)

            if prob_x1_y:
                z3 += prob_x1_y * log2(1.0/Prob_y_given_xs( y, x1.indices, x1.values ))

            #    print "z=%s \t z2=%s \t z3=%s" % (z,z2,z3)
    assert 0.0 <= z3 <= (H_Y() + epps)
    assert abs(z-z3) <= epps

    return z

def conditional_independence_test1():
    '''This determines whether p(y|x1,x2) = p(y|x1) * p(y|x2) for all states where p(y,x1,x2) > 0.'''
    if get_numXs() != 2:
        return

    all_X1s = joint_states( [0], True )
    all_X2s = joint_states( [1], True )

    all_matched = True
    # foreach state x1, x2, y...
    for x1 in all_X1s:
        for x2 in all_X2s:
            for y_state in Y_states():

                x1x2 = Xjointstate( x1.indices + x2.indices, x1.states + x2.states )
                p_yx1x2 = Prob_Xs_y( x1x2.indices, x1x2.states, y_state )

                if p_yx1x2 == 0.0:
                    continue

                p_y_given_x1x2 = Prob_y_given_xs( y_state, x1x2.indices, x1x2.states )

                # This expression is from Latham (2005) Equation (1) and (1.5)
                top = Prob_Y_equals_y(y_state) * Prob_y_given_xs(y_state,x1.indices,x1.states) * Prob_y_given_xs(y_state,x2.indices,x2.states)
                bottom = Prob_Xs_ind( x1x2.indices, x1x2.states )
                pind_y_given_x1_x2 = top / bottom

                if abs(p_y_given_x1x2 - pind_y_given_x1_x2) > eps:
                #                    print term.YELLOW + "\rTerms y=%s x1=%s x2=%s did not match. \t " % (y_state, x1.states[0], x2.states[0]),
                #                    print "p(y|x1,x2)=%s \t p_{ind}(y|x1,x2)=%s" % (p_y_given_x1x2, pind_y_given_x1_x2)
                #                    print term.NORMAL,
                    all_matched = False
                #                else:
                #                    print term.NORMAL + "Terms y=%s x1=%s x2=%s matched!" % (y_state, x1.states[0], x2.states[0])
                #                    print "p(y|x1,x2)=%s \t p_{ind}(y|x1,x2)=%s" % (p_y_given_x1x2, pind_y_given_x1_x2)
                #                    print term.NORMAL,

    if all_matched:
        print term.NORMAL + "\rp_ind test1: All p(y|x1,x2) == p_ind(y|x1,x2)!"
    else:
        print term.YELLOW + "\rp_ind test1: All p(y|x1,x2) != p_ind(y|x1,x2)!" + term.NORMAL

def conditional_independence_test2():
    '''This determines whether p(x1,x2|y) = p_{ind}(x1,x2|y) for all states where p(y,x1,x2) > 0.'''
    if get_numXs() != 2:
        return

    all_X1s = joint_states( [0], True )
    all_X2s = joint_states( [1], True )

    all_matched = True
    # foreach state x1, x2, y...
    for x1 in all_X1s:
        for x2 in all_X2s:
            for y_state in Y_states():

                x1x2 = Xjointstate( x1.indices + x2.indices, x1.states + x2.states )
                p_yx1x2 = Prob_Xs_y( x1x2.indices, x1x2.states, y_state )

                if p_yx1x2 == 0.0:
                    continue

                p_x1x2_given_y = Prob_Xs_given_y( x1x2.indices, x1x2.states, y_state )

                # This expression is from Latham (2005) Equation (1)
                pind_x1x2_given_y = Prob_Xs_given_y(x1.indices, x1.states, y_state) * Prob_Xs_given_y(x2.indices, x2.states, y_state)

                if abs(p_x1x2_given_y - pind_x1x2_given_y) > eps:
                #                    print term.YELLOW + "Terms y=%s x1=%s x2=%s did not match." % (y_state, x1.states[0], x2.states[0])
                #                    print "p(x1,x2|y)=%s \t p_{ind}(x1,x2|y)=%s" % (p_x1x2_given_y, pind_x1x2_given_y)
                #                    print term.NORMAL,
                    all_matched = False
                #                else:
                #                    print term.NORMAL + "Terms y=%s x1=%s x2=%s matched!" % (y_state, x1.states[0], x2.states[0])
                #                    print "p(y|x1,x2)=%s \t p_{ind}(y|x1,x2)=%s" % (p_y_given_x1x2, pind_y_given_x1_x2)
                #                    print term.NORMAL,


    if all_matched:
        print term.NORMAL + "\rp_ind test2: All p(x1,x2|y) == p_ind(x1,x2|y)!"
    else:
        print term.YELLOW + "\rp_ind test2: All p(x1,x2|y) != p_ind(x1,x2|y)!" + term.NORMAL

    raw_input('...')

def mutual_info_with_Y_multilowercond( X, Cs ):
    return H_Y() - H_Y_given_multiXs( X + Cs[0])



def pstar( mm, x1=None, x2=None, y=None, X1=None, X2=None, GIVENx1=None, GIVENx2=None, GIVENy=None ): #, GIVENx1=None, GIVENx2=None, GIVENy=None):
    '''Computes p*( ... | ... ) for any passed arguments.  If X1 or X2 specified, it means we average over that'''

    # one of these must be true
    assert x1 or x2 or y

    assert mm is not None

    # assert any x1, x2, y is valid
    if x1:
        assert X1 is None and GIVENx1 is None
        check_xs( x1.indices, x1.states )

    if x2:
        assert X2 is None and GIVENx2 is None
        check_xs( x2.indices, x2.states )

    if y:
        assert GIVENy is None
        check_y_value( y )

    ############ assert the GIVEN inputs are valid
    if GIVENx1:
        assert x1 is None and X1 is None
        check_xs( GIVENx1.indices, GIVENx1.states )
        # this is for programming convenience, we'll compute the joint p*, and then divide by the GIVEN p*
        x1 = GIVENx1

    if GIVENx2:
        assert x2 is None and X2 is None
        check_xs( GIVENx2.indices, GIVENx2.states )
        # this is merely for programming convenience
        x2 = GIVENx2

    if GIVENy:
        assert y is None
        check_y_value( GIVENy )
        y = GIVENy

    if X1:
        assert x1 is None and GIVENx1 is None
        check_Xs( X1 )

    if X2:
        assert x2 is None and GIVENx2 is None
        check_Xs( X2 )



    z = None
    if x1 and x2 and y:
        z = pdollar_a_y_cs( mm, x1, y, [x2] )

    elif x1 and x2 and y is None:
        # sum over all y
        z = sum([ pstar(mm, x1=x1,x2=x2,y=temp_y) for temp_y in Y_states() ])

    elif x1 and y and x2 is None and X2:
        # sum over all x2
        z = sum([ pstar(mm, x1=x1,x2=temp_x2,y=y) for temp_x2 in joint_states(X2,True) ])

    elif x2 and y and x1 is None and X1:
        # sum over all x2
        z = sum([ pstar(mm, x1=temp_x1,x2=x2,y=y) for temp_x1 in joint_states(X1,True) ])

    elif x1 and (x2 is y is None) and X2:
        # sum over all x2,y
        z = sum([ pstar(mm, x1=x1,x2=temp_x2,y=temp_y) for temp_x2 in joint_states(X2,True) for temp_y in Y_states() ])

    elif x2 and (x1 is y is None) and X1:
        # sum over all x1,y
        z = sum([ pstar(mm, x1=temp_x1,x2=x2,y=temp_y) for temp_x1 in joint_states(X1,True) for temp_y in Y_states() ])

    elif y and (x1 is x2 is None) and X1 and X2:
        # sum over all x1,x2
        z = sum([ pstar(mm, x1=temp_x1,x2=temp_x2,y=y) for temp_x1 in joint_states(X1,True) for temp_x2 in joint_states(X2,True) ])
    else:
        assert 0 == 1, "Should never get here."

    assert 0.0 <= z <= 1.0

    if z == 0.0:
        return 0.0

    # if there was anything GIVEN, divide by p* of the GIVENs.
    if GIVENx1 or GIVENx2 or GIVENy:

        if GIVENx1 is None and X1 is None:
            X1 = x1.indices
        if GIVENx2 is None and X2 is None:
            X2 = x2.indices

        z /= pstar(mm, x1=GIVENx1, x2=GIVENx2, y=GIVENy, X1=X1, X2=X2)

    assert 0.0 <= z <= 1.0

    return z



def PI_tests_n2():


    if get_numXs() != 2:
        return

    #    run_gce_tests2()

    int_12, int_1, int_2 = mutual_info_with_Y( [0,1] ), mutual_info_with_Y( [0] ), mutual_info_with_Y( [1] )
    #int_1_2 = mutual_info_with_Y([0]) - mutual_info_with_Y_multilowercond( [0], [[1]])
    #int_2_1 = mutual_info_with_Y([1]) - mutual_info_with_Y_multilowercond( [1], [[0]])

    #new_int_1_2 = I_intersect_with_Y( [[0],[1]] )
    #I_intersect_pdollar = round(I_intersect_with_Y_dollar( [[0],[1]] ),PRECISION-1)




    #assert round(int_1_2,PRECISION-1) == round(int_2_1,PRECISION-1), "Didn't come out equal!!!"
    """
    pi_1_2, pi_2_1 = int_1_2, int_2_1
    pi_1 = int_1 - pi_1_2
    pi_2 = int_2 - pi_1_2
    pi_12 = int_12 - pi_1 - pi_2 - pi_1_2

    pi_sum = pi_12 + pi_1 + pi_2 + pi_1_2

    new_int_1 = I_intersect_with_Y_dollar( [ [0] ] )
    new_int_2 = I_intersect_with_Y_dollar( [ [1] ] )
    new_int_1_1 = I_intersect_with_Y_dollar( [ [0],[0] ] )
    new_int_2_2 = I_intersect_with_Y_dollar( [ [1],[1] ] )
    print "I(X1:Y)=%s \t I(X2:Y)=%s" % ( mutual_info_with_Y([0]), mutual_info_with_Y([1]) )
    print "I_∩({1}:Y)=%s \t I_∩({2}:Y)=%s" % (new_int_1, new_int_2)
    print "I_∩({1}{1}:Y)=%s \t I_∩({2}{2}:Y)=%s" % (new_int_1_1, new_int_2_2)
    print "I_∩({1}{2}:Y)=%s \t I_∩({2}{1}:Y)=%s" % (int_1_2, int_2_1)

    """

    #MI = mutual_info_with_Y( [0] )
    #MI_X1_Y_lowercond_X2 = mutual_info_with_Y_lowercondXs( [0], [1] )
    #MI_X1_Y_cond_X2 = mutual_info_with_Y_condXs( [0], [1] )

    #print "I(X0:Y)=%s" % MI

    ###########################################
    # Now we're going to check the bounds for multicond.
    ###########################################
    #cond_1, cond_2 = mutual_info_with_Y_multilowercond( [0,1], [[0]]), mutual_info_with_Y_multilowercond( [0,1], [[1]])
    #multicond = mutual_info_with_Y_multilowercond( [0,1], [[0],[1]])
    #jointcond = mutual_info_with_Y_multilowercond( [0,1], [[0,1]] )

    #print "I(X_12:Y|X1)=%s \t I(X_12:Y|X2)=%s" % (cond_1,cond_2)
    #print "I(X_12:Y|X1,X2)=%s \t (multicond)" % multicond
    #print "I(X_12:Y|X_12)=%s \t (jointcond)" % jointcond

    #assert jointcond <= multicond
    #assert multicond <= min([cond_1,cond_2])

    #print "\t\t ** %s <= %s <= %s **" % (jointcond, multicond, min([cond_1,cond_2]))
    ###########################################




    #    print u"I(X_12:Y↓X1)=%s" % mutual_info_with_Y_multilowercond( [0,1], [[0]] )
    #    print u"I(X_12:Y↓X2)=%s" % mutual_info_with_Y_multilowercond( [0,1], [[1]] )
    #    print u"I(X_12:Y↓X1 or X2)=%s" % mutual_info_with_Y_multilowercond( [0,1], [[0],[1]] )

    """
    #print term.WHITE,
    print "\r\t\t {12}=%s" % pi_12
    print "\t {1}=%s \t {2}=%s" % ( pi_1, pi_2 )
    print "\t\t {1}{2}=%s" % (pi_1_2)

    print "\nI(Xs:Y)=%s \t sum_pis=%s" % ( int_12, pi_sum )
    print "hardcoded I_∩({1}{2}:Y)=%s" % new_int_1_2
    print "p$ I_∩({1}{2}:Y)=%s" % I_intersect_pdollar
    #print term.NORMAL
    """

    ###########################################
    # Synergy bounds for each state Y
    ###########################################
    print "-----------------------------------"
    for y_state in Y_states():
        p1, p2, both = specific_info_with_Y_equals_y( [0], y_state ), specific_info_with_Y_equals_y( [1], y_state ), specific_info_with_Y_equals_y( [0,1], y_state )
        min_p1p2 = min([p1,p2])
        D = both - p1 - p2

        min_syn, max_syn = max([0,D]), D+min_p1p2
        min_rdn, max_rdn = abs(max([-D,0])), min_p1p2
        print "Y=%s: \t " % y_state,
        print "p1p2=%s  \t p1=%s \t p2=%s" % (both, p1, p2) ,
        print "\t %s <= synergy <= %s" % (min_syn, max_syn),
        print "\t %s <= rdn <= %s" % (min_rdn, max_rdn)

    p1, p2, both = mutual_info_with_Y([0]), mutual_info_with_Y([1]), mutual_info_with_Y([0,1])
    min_p1p2 = min([p1,p2])
    D = both - p1 - p2
    min_syn, max_syn = max([0,D]), min([both, D+min_p1p2])
    min_rdn, max_rdn = p1 + p2 + min_syn - both, p1 + p2 + max_syn - both
    print "Average:",
    print "p1p2=%s  \t p1=%s \t p2=%s" % (both, p1, p2) ,
    print "\t %s <= synergy <= %s" % (min_syn, max_syn),
    print "\t %s <= rdn <= %s" % (min_rdn, max_rdn)

def I_intersect_with_Y_dollar( Xs ):

    for X_i in Xs:
        check_Xs( X_i )

    #    print Xs
    z = mutual_info_with_Y( Xs[0] )

    if len(Xs) == 1:
        return z

    elif len(Xs) == 2:
        z -= mutual_info_with_Y_multilowercond( Xs[0], [Xs[1]] )

    elif len(Xs) == 3:
        z -= mutual_info_with_Y_multilowercond( Xs[0], [Xs[1]] )
        z -= mutual_info_with_Y_multilowercond( Xs[0], [Xs[2]] )
        z += mutual_info_with_Y_multilowercond( Xs[0], [Xs[1],Xs[2]] )

    else:
        assert 0 == 1, "Shouldn't have gotten here."

    #    assert 0.0 <= z
    return z



def I_intersect_with_Y( Xs ):
    '''returns I_intersect({Xs_1}{Xs_2}...{Xs_k}:Y) =
        \sum_y p(y) \sum_{xs \in Xs} \prod_{i=1}^k p(x_i|y) log [ \sum_{y`} p(y`) \prod_{i=1}^k p(x_i|y`) / (\prod_{i=1}^k p(x_i) ]
    '''

    for X_i in Xs:
        check_Xs( X_i )

    Xs_jointstates = joint_states_of_Cs( Xs )

    z = 0.0

    for y_state in Y_states():
        this_prob_y = Prob_Y_equals_y( y_state )

        assert 0.0 < this_prob_y <= 1.0

        # for each jointstate of the Xs...
        summ = 0.0
        for xs_jointstate in Xs_jointstates:

            term1 = 1.0
            for x_i in xs_jointstate:
                term1 *= Prob_Xs_given_y( x_i.indices, x_i.states, y_state )

            if term1 == 0.0:
                continue

            # the bottom term of the log
            bottom = 1.0
            for x_i in xs_jointstate:
                bottom *= Prob_Xs( x_i.indices, x_i.states )

            # the top term of the log
            top = 0.0
            for yprime_state in Y_states():
                this_prob_yprime_state = Prob_Y_equals_y( yprime_state )

                top_prod = 1.0
                for x_i in xs_jointstate:
                    top_prod *= Prob_Xs_given_y( x_i.indices, x_i.states, yprime_state )

                top += this_prob_yprime_state * top_prod

            logterm = log2( top / bottom )

            summ += term1 * logterm

        z += this_prob_y * summ


    assert 0.0 <= z
    min_MI = min([ mutual_info_with_Y(X_i) for X_i in Xs ])

    if round(z,PRECISION) > round(min_MI,PRECISION):
        print "z=%s  min_MI=%s" % (z, min_MI)

    #assert round(z,PRECISION) <= round(min_MI,PRECISION)


    return round(z, PRECISION + 1)


def mutual_info_dollar_with_Y( Xs ):
    '''I^$( Xs : Y ) = DKL[ p$(x,y) || p$(x) p$(y) ]'''

    check_Xs( Xs )

    # C is everything outside of Xs
    #    C = [ x_index for x_index in get_Xindices() ]
    #    check_Xs( C )


    z = 0.0
    # define the C as everything outside of X
    X_states = joint_states( Xs, True )

    for x in X_states:
        for y_state in Y_states():

            prob_xy = pdollar_a_y_noC( x, y_state )

            if prob_xy == 0.0:
                continue

            # calculate p$(a) and p$(y) by summing over all of the Y and X states.
            prob_x, prob_y = 0.0, 0.0
            for temp_y in Y_states():
                prob_x += pdollar_a_y_noC( x, temp_y )
            for temp_x in X_states:
                prob_y += pdollar_a_y_noC( temp_x, y_state )

            top = prob_xy
            bottom = prob_x * prob_y

            z += prob_xy * log2( top / bottom )

    #print "z=%s" % z
    assert 0.0 <= z
    return round(z,PRECISION)


def test_rdn_measures_with_Y( X0, X1 ):

    check_Xs( X0 )
    check_Xs( X1 )
    both = X0 + X1
    check_Xs(both)

    # foreach y_state...
    for y_state in Y_states():
        this_minI = min( specific_info_with_Y_equals_y(X0,y_state), specific_info_with_Y_equals_y(X1,y_state) )
        rdn_Nihat = specific_info_with_Y_equals_y(X0, y_state) - specific_info_with_y_lowercondXs( X0, y_state, X1 )


        this_specific_cond = specific_info_with_y_given_condXs( X0, y_state, X1 )

        print "\t y='%s' \t minI=%s \t rdn_Nihat=%s \t rdn_Virgil=%s \t I(X0:y|X1)=%s" % ( y_state, this_minI, rdn_Nihat, rdn_Virgil, this_specific_cond)

    min_I = min( mutual_info_with_Y( X0 ), mutual_info_with_Y( X1 ) )
    rdn_Nihat = mutual_info_with_Y( X0 ) - mutual_info_with_Y_lowercondXs( X0, X1 )
    rdn_Nihat2 = mutual_info_with_Y( X1 ) - mutual_info_with_Y_lowercondXs( X1, X0 )


    print "AVR min_I=%s" % min_I
    print "AVR rdn_Nihat=%s \t rdn_Nihat2=%s" % (rdn_Nihat, rdn_Nihat2)

    assert round(rdn_Nihat,PRECISION-1) == round(rdn_Nihat2,PRECISION-1)

def pdollar_a_y_multi_cs( a, y_state, cs ):
    '''computes the pdollar for multiple c singletons'''

    # check our inputs
    assert len(cs) >= 2
    check_xs( a.indices, a.states )
    check_y_value( y_state )
    for c in cs:
        check_xs( c.indices, c.states )

        # for now only doing k=2.
    #    assert len(cs) == 2

    #    z = Prob_Xs( a.indices, a.states )
    #    z *= Prob_y_given_xs( y_state, a.indices, a.states )
    #    for c in cs:
    #        z *= Prob_Xs_given_y( c.indices, c.states, y_state )
    #    assert 0.0 <= z <= 1.0

    # z3 = p(y) p(a|y) \prod_{c in cs} p(c|y)
    z3 = Prob_Y_equals_y( y_state )
    z3 *= Prob_Xs_given_y( a.indices, a.states, y_state )
    for c in cs:
        z3 *= Prob_Xs_given_y( c.indices, c.states, y_state )


    assert 0.0 <= z3 <= 1.0

    return z3


def pdollar_a_y_noC( a, y_state ):
    '''computes p$(a,y) given that C is None'''

    global MARKOV_MODEL
    check_xs( a.indices, a.states )
    check_y_value( y_state )

    assert type(MARKOV_MODEL) is tuple and len(MARKOV_MODEL) == 3, "MARKOV_MODEL wasn't valid for this"
    indeg_A, indeg_Y, indeg_C = MARKOV_MODEL

    assert indeg_A in [0,1,2,3,4]
    assert indeg_Y in [0,1,2,3,4]

    z = None
    if indeg_A in [0,2]:
        z = Prob_Xs( a.indices, a.states )
    elif indeg_A in [1,3,4]:
        z = Prob_Xs_given_y( a.indices, a.states, y_state )
    else:
        assert 0 == 1, "should be impossible"

    if indeg_Y in [0,2]:
        z *= Prob_Y_equals_y( y_state )
    elif indeg_Y in [1,3,4]:
        z *= Prob_y_given_xs( y_state, a.indices, a.states )
    else:
        assert 0 == 1, "should be impossible"

    assert 0.0 <= z <= 1.0

    return z

def pdollar_x1_x2_y( x1, x2, y_state ):
    '''comutes p$(a,y,c_1 ... c_k) = p(y_state) * p(x1|y) * p(x2|y)'''

    assert mm is not None

    check_xs( a.indices, a.states )
    check_y_value( y_state )

    for c in cs:
        check_xs( c.indices, c.states )

    if len(cs) >= 2:
        #print "- conditioning on %s Cs!" % len(cs)
        return pdollar_a_y_multi_cs( a, y_state, cs )

    z = None
    c = cs[0]

    # make a state of the joint random variable of ac
    ac = Xjointstate( a.indices + c.indices, a.states + c.states )


    assert len(cs) == 1

    if not (0.0 <= z <= 1.0):

        if mm == 'AY_PRODUCT':
            top = Prob_y_given_xs( y_state, a.indices, a.states ) * Prob_y_given_xs( y_state, c.indices, c.states )
            bottom = sum([ Prob_y_given_xs( yprime, a.indices, a.states ) * Prob_y_given_xs( yprime, c.indices, c.states ) for yprime in Y_states() ])
            z = Prob_Xs( a.indices, a.states ) * Prob_Xs( c.indices, c.states ) * (top/bottom)

        elif mm == 'AY_SUM':
            top = Prob_y_given_xs( y_state, a.indices, a.states ) + Prob_y_given_xs( y_state, c.indices, c.states )
            bottom = sum([ Prob_y_given_xs( yprime, a.indices, a.states ) + Prob_y_given_xs( yprime, c.indices, c.states ) for yprime in Y_states() ])
            z = Prob_Xs( a.indices, a.states ) * Prob_Xs( c.indices, c.states ) * (top/bottom)

        elif mm == 'MASKELL2008':
            # uses equation (4) from Maskell (2008), "A Bayesian approach to fusing uncertain, imprecise, and conflicting information"
            # p( y | x1 U x2 ) = p(y)*p(x1|y)*p(x2|y) / p(x1,x2)

            top = Prob_Y_equals_y( y_state )
            top *= Prob_Xs_given_y(a.indices, a.states, y_state)
            top *= Prob_Xs_given_y(c.indices, c.states, y_state)
            bottom = Prob_Xs(ac.indices, ac.states)

            #            top = Prob_Xs_y(a.indices,a.states, y_state)
            #            top *= Prob_Xs_y(c.indices,c.states, y_state)
            #            bottom = Prob_Xs( ac.indices, ac.states )


            if top == 0.0:
                z = 0.0

            elif top and bottom == 0.0:
                z = 0.0
                print "top=%s \t bottom=%s" % (top,bottom)
                raw_input('...')
            else:
                z = top / bottom

            if z > 1.0:
                print "[high] z=%s \t top=%s \t bottom=%s" % (z, top, bottom)
                raw_input('...')

                #    elif MARKOV_MODEL == 'X->Y->Z':
                #        z = Prob_Xs( a.indices, a.states ) * Prob_y_given_xs( y_state, a.indices, a.states ) * Prob_Xs_given_y( c.indices, c.states, y_state )
                #    elif MARKOV_MODEL == 'Z->Y->X':
                #        z = Prob_Xs( c.indices, c.states ) * Prob_y_given_xs( y_state, c.indices, c.states ) * Prob_Xs_given_y( a.indices, a.states, y_state )
                #    elif MARKOV_MODEL == 'SPOKE_OUTWARDS':
                #        z = Prob_Y_equals_y( y_state ) * Prob_Xs_given_y( a.indices, a.states, y_state ) * Prob_Xs_given_y( c.indices, c.states, y_state )
                #    elif MARKOV_MODEL == 'SPOKE_INWARDS':
                #        z = Prob_Xs( a.indices, a.states ) * Prob_Xs( c.indices, c.states ) * Prob_y_given_xs( y_state, ac.indices, ac.states )
                #    elif MARKOV_MODEL == 'SPOKE_UNDIRECTED':
                #        z = Prob_Xs_given_y( a.indices, a.states, y_state ) * Prob_y_given_xs( y_state, ac.indices, ac.states ) * Prob_Xs_given_y( c.indices, c.states, y_state )
                #    elif MARKOV_MODEL == 'SPOKE4':
                #        z = Prob_Xs_given_y( a.indices, a.states, y_state ) * Prob_y_given_xs( y_state, a.indices, a.states ) * Prob_Xs_given_y( c.indices, c.states, y_state )
                #    elif MARKOV_MODEL == 'SPOKE5':
                #        z = Prob_Xs_given_y( a.indices, a.states, y_state ) * Prob_y_given_xs( y_state, ac.indices, acb.states ) * Prob_Xs( c.indices, c.states )
                #    elif MARKOV_MODEL == 'TRIANGLE':
                #        z = Prob_xs_given_condxs( a.indices, a.states, c.indices, c.states ) * Prob_y_given_xs( y_state, ac.indices, ac.states ) * Prob_xs_given_condxs( c.indices, c.states, a.indices, a.states )
                #    elif MARKOV_MODEL == 'TRIANGLE2':
                #        z = Prob_xs_given_condxs( a.indices, a.states, c.indices, c.states ) * Prob_y_given_xs( y_state, ac.indices, ac.states ) * Prob_Xs( c.indices, c.states )
                #    elif MARKOV_MODEL == 'TRIANGLE3':
                #        z = Prob_Xs( a.indices, a.states ) * Prob_y_given_xs( y_state, ac.indices, ac.states ) * Prob_xs_given_condxs( c.indices, c.states, a.indices, a.states )
                #    elif MARKOV_MODEL == 'SPOKE6':
                #        z = Prob_Xs_given_y( a.indices, a.states, y_state ) * Prob_y_given_xs( y_state, c.indices, c.states ) * Prob_Xs_given_y( c.indices, c.states, y_state )
                #    elif MARKOV_MODEL == 'SPOKE7':
                #        z = Prob_Xs( a.indices, a.states ) * Prob_y_given_xs( y_state, ac.indices, ac.states ) * Prob_Xs_given_y( c.indices, c.states, y_state )
        elif mm == 'NULL':
            z = Prob_Xs( a.indices, a.states ) * Prob_y_given_xs( y_state, a.indices, a.states ) * Prob_xs_given_condxs_y( c.indices, c.states, a.indices, a.states, y_state )
        elif mm == 'LATHAM':
            z = Prob_Y_equals_y(y_state) * Prob_Xs_given_y(a.indices,a.states,y_state) * Prob_Xs_given_y(c.indices,c.states,y_state)

        elif type(mm) is tuple and len(mm) == 3:
            indeg_A, indeg_Y, indeg_C = mm

            assert indeg_A in [0,1,2,3,4,5]
            assert indeg_Y in [0,1,2,3,4,5]
            assert indeg_C in [0,1,2,3,4,5]

            z1, z2, z3 = None, None, None

            if indeg_A == 0:
                z1 = Prob_Xs( a.indices, a.states )
            elif indeg_A == 1:
                z1 = Prob_Xs_given_y( a.indices, a.states, y_state )
            elif indeg_A == 2:
                z1 = Prob_xs_given_condxs( a.indices, a.states, c.indices, c.states )
            elif indeg_A == 3:
                z1 = Prob_xs_given_condxs_y( a.indices, a.states, c.indices, c.states, y_state )
            elif indeg_A == 4:
                z1 = Prob_Xs_given_y( a.indices, a.states, y_state ) * Prob_xs_given_condxs( a.indices, a.states, c.indices, c.states )
            elif indeg_A == 5:
                z1 = Prob_Xs_given_y( a.indices, a.states, y_state ) * Prob_xs_given_condxs( a.indices, a.states, c.indices, c.states )
                if z1:
                    z1 /= sum([ Prob_Xs_given_y( aprime.indices, aprime.states, y_state ) * Prob_xs_given_condxs( aprime.indices, aprime.states, c.indices, c.states ) for aprime in joint_states(a.indices,True) ])
            else:
                assert 0 == 1, "should be impossible"

            if indeg_Y == 0:
                z2 = Prob_Y_equals_y( y_state )
            elif indeg_Y == 1:
                z2 = Prob_y_given_xs( y_state, a.indices, a.states )
            elif indeg_Y == 2:
                z2 = Prob_y_given_xs( y_state, c.indices, c.states )
            elif indeg_Y == 3:
                z2 = Prob_y_given_xs( y_state, ac.indices, ac.states )
            elif indeg_Y == 4:
                z2 = Prob_y_given_xs( y_state, a.indices, a.states ) * Prob_y_given_xs( y_state, c.indices, c.states )
            elif indeg_Y == 5:
                z2 = Prob_y_given_xs( y_state, a.indices, a.states ) * Prob_y_given_xs( y_state, c.indices, c.states )
                if z2:
                    z2 /= sum([ Prob_y_given_xs( yprime, a.indices, a.states ) * Prob_y_given_xs( yprime, c.indices, c.states ) for yprime in Y_states() ])
            else:
                assert 0 == 1, "should be impossible"

            if indeg_C == 0:
                z3 = Prob_Xs( c.indices, c.states )
            elif indeg_C == 1:
                z3 = Prob_xs_given_condxs( c.indices, c.states, a.indices, a.states )
            elif indeg_C == 2:
                z3 = Prob_Xs_given_y( c.indices, c.states, y_state )
            elif indeg_C == 3:
                z3 = Prob_xs_given_condxs_y( c.indices, c.states, a.indices, a.states, y_state )
            elif indeg_C == 4:
                z3 = Prob_xs_given_condxs( c.indices, c.states, a.indices, a.states ) * Prob_Xs_given_y( c.indices, c.states, y_state )
            elif indeg_C == 5:
                z3 = Prob_xs_given_condxs( c.indices, c.states, a.indices, a.states ) * Prob_Xs_given_y( c.indices, c.states, y_state )
                if z3:
                    z3 /= sum([ Prob_xs_given_condxs( cprime.indices, cprime.states, a.indices, a.states ) * Prob_Xs_given_y( cprime.indices, cprime.states, y_state ) for cprime in joint_states(c.indices,True) ])
            else:
                assert 0 == 1, "should be impossible"

            z = z1*z2*z3

        else:
            assert 0 == 1, "Didn't know the MARKOV_MODEL='%s'" % mm

    assert 0.0 <= z <= 1.0

    return z

def pdollar_a_given_cs( a, cs ):
    '''computes p$( a | c_1 ... c_k ) = p$( a, c_1 ... c_k ) / p$( c_1 ... c_k )'''

    check_xs( a.indices, a.states )
    for c in cs: check_xs( c.indices, c.states )

    top = pdollar_a_cs( a, cs )
    bottom = pdollar_cs( cs, a.indices )

    z = top / bottom

    assert 0.0 <= z <= 1.0

    return z







def joint_states_of_Cs( Cs ):
    '''returns a big list of sequences of Xjointstates.
    '''
    for C in Cs: check_Xs( C )

    assert 1 <= len(Cs) <= 4

    z = None
    if len(Cs) == 1:
        z = [ (c_1,) for c_1 in joint_states(Cs[0],True) ]
    elif len(Cs) == 2:
        z = [ (c_1,c_2) for c_1 in joint_states(Cs[0],useclass=True) for c_2 in joint_states(Cs[1],useclass=True) ]
    elif len(Cs) == 3:
        z = [ (c_1,c_2,c_3) for c_1 in joint_states(Cs[0],useclass=True) for c_2 in joint_states(Cs[1],useclass=True) for c_3 in joint_states(Cs[2],useclass=True) ]
    elif len(Cs) == 4:
        z = [ (c_1,c_2,c_3,c_4) for c_1 in joint_states(Cs[0],useclass=True) for c_2 in joint_states(Cs[1],useclass=True) for c_3 in joint_states(Cs[2],useclass=True) for c_4 in joint_states(Cs[3],useclass=True) ]
    else:
        assert 0 == 1, "Should be impossible to get here."

    return z

def mutual_info_with_Y_condXs( Xs, condXs ):
    '''returns the conditional mutual information I(Xs:Y|condXs) = H(Xs,condXs) + H(Y,condXs) - H(Xs,condXs,Y) - H(condXs)'''

    check_Xs( Xs )
    check_Xs( condXs )

    Xs_and_condXs = Xs + condXs
    check_Xs( Xs_and_condXs )

    z = H_Xs( Xs_and_condXs )
    z += H_Xs_Y( condXs )
    z -= H_Xs_Y( Xs_and_condXs )
    z -= H_Xs( condXs )

    z2 = 0.0
    for y_state in Y_states():
        this_prob = Prob_Y_equals_y( y_state )
        z2 += this_prob * specific_info_with_y_given_condXs( Xs, y_state, condXs )


    z3 = 0.0
    all_condXstates = joint_states( condXs )
    for condXstates in all_condXstates:
        this_prob = Prob_Xs( condXs, condXstates )
        z3 += this_prob * mutual_info_with_Y_conditional_statedep( Xs, condXs, condXstates )


    assert 0.0 <= round(z,PRECISION)
    assert 0.0 <= round(z2,PRECISION)
    assert 0.0 <= round(z3,PRECISION)

    assert abs(z-z2) <= eps and abs(z2-z3) <= eps

    #print "z=%s \t z2=%s \t z3=%s" % (round(z,PRECISION), round(z2,PRECISION), round(z3,PRECISION))

    return z

def mutual_info_with_Y_lowercond_statedep( Xindices, condXindices, condXstates ):
    '''Computes I( Xs: Y \downarrow condXs=condxs ) = DKL[ p*(x,y|condx) || p*(x|condx) p*(y|condx) ]'''

    check_Xs( Xindices )
    check_xs( condXindices, condXstates )

    all_Xstates = joint_states( Xindices )

    z = 0.0
    # foreach Xstate and Ystate...
    for Xstates in all_Xstates:
        for y_state in Y_states():
            term1 = pstar_xs_y_given_condxs( Xindices, Xstates, y_state, condXindices, condXstates )

            if term1 == 0.0:
                continue

            term2 = pstar_xs_given_condxs( Xindices, Xstates, condXindices, condXstates )
            term2 *= pstar_y_given_condxs( y_state, condXindices, condXstates, Xindices )

            z += term1 * log2( term1 / term2 )

    assert 0.0 <= z

    # check z's upperbound, z <= I(X:Y) and z <= I(X:Y|Z=z)
    #    bound1 = mutual_info_with_Y_conditional_statedep( Xindices, condXindices, condXstates )
    #    bound2 = mutual_info_with_Y( Xindices )

    #    print "z=%s \t bound1=%s \t bound2=%s" % (z,bound1,bound2)
    #    assert z <= bound1
    #    assert z <= bound2

    return z


def mutual_info_with_Y_conditional_statedep( Xindices, condXindices, condXstates ):
    '''Computes I( Xs: Y | condXs=condxs ) = DKL[ p(x,y|condx) || p(x|condx)*p(y|condx) ]'''

    check_Xs( Xindices )
    check_xs( condXindices, condXstates )

    all_Xstates = joint_states( Xindices )

    z = 0.0
    # foreach Xstate and Ystate...
    for Xstates in all_Xstates:
        for y_state in Y_states():
            term1 = Prob_xs_y_given_condxs( Xindices, Xstates, y_state, condXindices, condXstates )

            if term1 == 0.0:
                continue

            term2 = Prob_xs_given_condxs( Xindices, Xstates, condXindices, condXstates )
            term2 *= Prob_y_given_xs( y_state, condXindices, condXstates )

            z += term1 * log2( term1 / term2 )

    assert 0.0 <= z

    return z

def Prob_xs_y_given_condxs( Xindices, Xstates, y_state, condXindices, condXstates ):
    '''Compuates p(x,y|condx) = p(x,y,condx) / p(condx) '''

    check_xs( Xindices, Xstates )
    check_xs( condXindices, condXstates )
    check_y_value( y_state )

    both_indices = Xindices + condXindices
    both_states = Xstates + condXstates

    top = Prob_Xs_y( both_indices, both_states, y_state )
    bottom = Prob_Xs( condXindices, condXstates )

    z = top / bottom

    assert 0.0 <= z <= 1.0

    return z

def specific_info_with_y_given_condXs( Xs, y_state, condXs ):
    '''computesI( Xs : Y=y | condXs ) = DKL[ p*(xs, condxs| Y=y ) || p(xs | condxs ) * p(condxs|Y=y) ]'''

    check_Xs( Xs )
    check_Xs( condXs )

    Xs_and_condXs = Xs + condXs
    check_Xs( Xs_and_condXs )

    all_Xstates = joint_states( Xs )
    all_condXstates = joint_states( condXs )

    z = 0.0
    # foreach Xstate and condXstate...
    for Xstates in all_Xstates:
        for condXstates in all_condXstates:
            term1 = Prob_Xs_given_y( Xs_and_condXs, Xstates + condXstates, y_state )

            if term1 == 0.0:
                continue

            #term2 = Prob_xs_given_condxs_y( Xs, Xstates, condXs, condXstates, y_state )
            term2 = Prob_xs_given_condxs( Xs, Xstates, condXs, condXstates )
            term2 *= Prob_Xs_given_y( condXs, condXstates, y_state )

            z += term1 * log2( term1 / term2 )


    assert 0.0 <= z
    #print "regularcond: Xs=%s y=%s condXs=%s" % ( Xs, y_state, condXs )
    return z

def Prob_xs_given_condxs_y( Xindices, Xstates, condXindices, condXstates, y_state ):
    '''calculates p(x|condx,y) = p(x,y,condx) / p(y,condx)'''

    check_xs( Xindices, Xstates )
    check_xs( condXindices, condXstates )
    check_y_value( y_state )

    both_indices = Xindices + condXindices
    both_states = Xstates + condXstates

    top = Prob_Xs_y( both_indices, both_states, y_state )
    bottom = Prob_Xs_y( condXindices, condXstates, y_state )

    if top == 0.0:
        return 0.0

    z = top / bottom

    assert 0.0 <= z <= 1.0
    return z

def Prob_xs_given_condxs( Xs, Xstates, condXs, condXstates ):
    '''returns p(xs|condxs) = p(xs,condxs) / p(condxs)'''

    check_xs( Xs, Xstates )
    check_xs( condXs, condXstates )

    both = Xs+condXs
    both_states = Xstates+condXstates

    top = Prob_Xs( both, both_states )
    bottom = Prob_Xs( condXs, condXstates )

    z = top / bottom

    assert 0.0 <= z <= 1.0

    return z


def D( Xs, ystate ):
    '''computes the D of whole - the parts for this state Y=y'''
    check_Xs( Xs )
    check_y_value( ystate )

    z = specific_info_with_Y_equals_y( Xs, ystate )

    part_infos = [ specific_info_with_Y_equals_y([X_index],ystate) for X_index in Xs ]

    z -= sum(part_infos)

    return z

def min_max_synergy_Y( Xs ):
    '''returns the AVERAGED minimum and maximum synergy'''

    check_Xs( Xs )

    avr_min, avr_max = 0.0, 0.0

    for ystate in Y_states():
        this_prob = Prob_Y_equals_y( ystate )

        this_min, this_max = min_max_synergy_y( Xs, ystate )

        if this_min is None and this_max is None:
            return [None,None]

        avr_min += this_prob * this_min
        avr_max += this_prob * this_max


    if avr_min is None and avr_max is None:
        return [None,None]

    # 0 <= avr_min <= avr_max <= I(Xs:Y)
    assert 0.0 <= avr_min
    assert avr_min <= avr_max
    assert avr_max <= mutual_info_with_Y( Xs )+eps

    #avr_min, avr_max = round(avr_min,PRECISION), round(avr_max,PRECISION)

    return (avr_min, avr_max)

def min_max_synergy_y( Xs, ystate ):
    '''returns the minimum and maximum synergy for a given Y=ystate'''

    check_Xs( Xs )
    check_y_value( ystate )

    if len(Xs) != 2:
        return (None,None)

    # define the DIFFERENCE, D
    whole, part1, part2 = specific_info_with_Y_equals_y( Xs, ystate ), specific_info_with_Y_equals_y( [Xs[0]], ystate ), specific_info_with_Y_equals_y( [Xs[1]], ystate )
    D = whole - part1 - part2

    z_min = max( [0.0, D] )
    z_max = min( [whole, D + min([part1,part2]) ] )

    z_min, z_max = round(z_min,PRECISION+1), round(z_max,PRECISION+1)

    #    if test_this_synergy_value is not None:

    return (z_min, z_max)



def Delta_I( Xs ):
    '''Computes the Niremberg- \Delta-I measure.  A measure of ``correlation importance''.
    It is defined as:

    TC( X_1 ; ... ; X_n | Y ) - DKL[ P(X_{1...n}|| Pr_{ind}(X_{1...n}) ]
    '''

    check_Xs( Xs )

    ########################################################################################
    z = total_correlation_given_Y( Xs )

    allstates = joint_states( Xs )

    summ = 0.0
    for Xstate in allstates:
        this_prob = Prob_Xs( Xs, Xstate )

        this_prob_ind = Prob_Xs_ind( Xs, Xstate )

        assert 0.0 < this_prob <= 1.0

        summ += this_prob * log2( this_prob / this_prob_ind )

    #   print "summ=%s" % summ

    # z = TC( X_1 ; ... ; X_n | Y ) - DKL[ P(X_{1...n}|| Pr_{ind}(X_{1...in}) ]
    z = total_correlation_given_Y( Xs ) - summ

    ########################################################################################
    assert 0.0 <= z

    if len(Xs) != 2:
        return round(z,PRECISION)
        ########################################################################################

    # This is the typical way of computing \Delta I for n=2 via the DKL[ p(y|x1,x2) || p*(y|x1,x2) ]
    # we compute it this way as a sanity check for the MARKOV_MODEL
    assert len(Xs) == 2
    X1, X2 = [Xs[0]], [Xs[1]]
    x1s, x2s = joint_states(X1,True), joint_states(X2,True)

    return round(z,PRECISION)


def Prob_Xs_ind( Xs, Xstates ):
    '''Calculates the Prob_ind of Prob( Xs = Xstate ) = sum_{y in Y} Prob(Y = y) \prod_i Pr( x_i | y )'''

    check_Xs( Xs )

    z = 0.0

    for y_value in Y_states():
        prob_y = Prob_Y_equals_y( y_value )

        product = 1.0

        for this_X, this_Xstate in zip(Xs,Xstates):
            product *=  Prob_Xs_given_y( [this_X], [this_Xstate], y_value )

        z += prob_y * product

    assert 0.0 <= z <= 1.0

    return z

def Interaction_Information_v4( Xs, norm=False):
    '''Computes the interaction information among {X_1, ... , X_n ,Y}.  Note this is the BELL FORMULATION, not McGill.'''

    #    print "Doing interaction information.  Xs=%s" % Xs

    # We will generate ALL SUBSETS of len(Xs)+1.
    # We define the index of len(Xs) to be Y
    size_of_set = len(Xs)

    summ = 0.0

    for i, S in enumerate(genallsubsets( size_of_set )):
    #        print "%s) S=%s" % (i, S)
        Scopy = deepcopy(S)
        term1 = (-1.0)**(size_of_set + 1 - len(Scopy))
        term2 = H_Y_given_Xs( Scopy )

        summ += term1 * term2

    # include the case there S={}
    summ += (-1)**(size_of_set + 1) * H_Y()

    # Just a manual check to not have -0.0
    if summ == -0.0:
        summ = 0.0

    if norm:
        return None

    return summ

def Interaction_Information_v3( Xs, norm=False):
    '''Computes the interaction information among {X_1, ... , X_n ,Y}.  Note this is the BELL FORMULATION, not McGill.'''

    #    print "Doing interaction information.  Xs=%s" % Xs

    # We will generate ALL SUBSETS of len(Xs)+1.
    # We define the index of len(Xs) to be Y
    index_of_Y = len(Xs)

    numXs = len(Xs)
    size_of_full_set = len(Xs) + 1

    summ = 0.0

    for i, S in enumerate(genallsubsets( size_of_full_set )):
    #        print "%s) S=%s" % (i, S)
        Scopy = deepcopy(S)
        term1 = (-1.0)**(numXs - len(Scopy))


        term2 = None

        if index_of_Y in Scopy:
            Scopy.remove( index_of_Y )
            if Scopy:
                term2 = H_Xs_Y( Scopy )
            else:
                term2 = H_Y()

        else:
            term2 = H_Xs( Scopy )

        assert index_of_Y not in Scopy

        summ += term1 * term2


    # Just a manual check to not have -0.0
    if summ == -0.0:
        summ = 0.0

    if norm:
        return None

    return summ

def Interaction_Information_v2( Xs ):
    '''Computes the Interaction Information among {X_1, ... X_n, Y} using the \sum_{S \subset {X_1, ... , X_n} I(S:Y) .  This should be equivalent to the Bell expression'''

    # We will generate ALL SUBSETS of len(Xs)+1.
    # We define the index of len(Xs) to be Y
    size_of_set = len(Xs)

    summ = 0.0
    #summ = mutual_info_with_Y(Xs)

    for i, S in enumerate(genallsubsets( size_of_set )):
        #print "%s) S=%s" % (i, S)
        Scopy = deepcopy(S)
        term1 = (-1.0)**(size_of_set - len(Scopy) )

        assert len(Scopy) >= 1, "ran the null set, S={}"

        term2 = mutual_info_with_Y(Scopy)

        summ += term1 * term2


    # Just a manual check to not have -0.0
    if summ == -0.0:
        summ = 0.0


    return summ

def Interaction_Information( Xs, norm=False):
    '''Computes the interaction information among {X_1, ... , X_n ,Y}.  Note this is the BELL FORMULATION, not McGill.'''

    #    print "Doing interaction information.  Xs=%s" % Xs

    # We will generate ALL SUBSETS of len(Xs)+1.
    # We define the index of len(Xs) to be Y
    index_of_Y = len(Xs)

    size_of_full_set = len(Xs) + 1

    summ = 0.0

    for i, S in enumerate(genallsubsets( size_of_full_set )):
    #        print "%s) S=%s" % (i, S)
        Scopy = deepcopy(S)
        term1 = (-1.0)**(size_of_full_set + 1 - len(Scopy))


        term2 = None

        if index_of_Y in Scopy:
            Scopy.remove( index_of_Y )
            if Scopy:
                term2 = H_Xs_Y( Scopy )
            else:
                term2 = H_Y()

        else:
            term2 = H_Xs( Scopy )

        assert index_of_Y not in Scopy

        summ += term1 * term2


    # Just a manual check to not have -0.0
    if summ == -0.0:
        summ = 0.0

    if norm:
        return None

    return summ

def Synmin( Xs, norm=False ):
    '''Computes the Varadan Synergy, or SYN-MIN, SYN-MIN(Xs:Y)'''

    check_Xs( Xs )
    z = mutual_info_with_Y( Xs )

    if norm:
        return None

    #print "All partitions of size len(Xs)=%s" % len(Xs)

    max_over_partitions = None
    # The partition that reached max_over_partitions
    max_partition = None

    for i, P in enumerate(partitions( Xs )):

        # skip the total partition
        if len(P) == 1:
            continue

        info_in_this_partition = sum( [ mutual_info_with_Y(part) for part in P ] )

        #        print "\t %s) P=%s \t\t info_in_this=%s \t\t max_so_far=%s \t max_P=%s" % (i, P, info_in_this_partition, max_over_partitions, max_partition)

        if max_over_partitions is None:
            max_over_partitions = info_in_this_partition
            max_partition = deepcopy(P)

        if info_in_this_partition >= max_over_partitions:
            assert max_over_partitions is not None
            max_partition = deepcopy(P)

        max_over_partitions = max(max_over_partitions, info_in_this_partition)

    #    print "max_partition=%s" % (max_partition)
    z -= max_over_partitions


    #    print "z=%s  \t Checkik=%s" % (z, Chechik( Xs ) )
    assert z <= Chechik( Xs )+epps

    return z


def Chechik( Xs, norm=False ):
    '''Computes the Chechik Synergy'''
    check_Xs(Xs)

    z1 = mutual_info_with_Y( Xs )

    for index in Xs:
        z1 -= mutual_info_with_Y( [index] )

    #    z1 = round(z1,PRECISION)

    z2 = total_correlation_given_Y( Xs ) - total_correlation( Xs )
    #    z2 = round(z2,PRECISION)

    assert abs(z1-z2) <= PRECISION, "z1=%s \t z2=%s" % (z1,z2)

    if norm:
    #        z1 /= mutual_info_with_Y( Xs )
        z1 = None
    #        assert -1.0 <= z1 <= 1.0
    #        z1 = round(z1,PRECISION)

    return z1

def print_redundancy_atoms( Xs ):
    '''prints some redundancy atoms that interest us'''
    check_Xs( Xs )
    if len(Xs) != 3:
        return

    summ=0.0

    for y_value in Y_states():
        this_prob = Prob_Y_equals_y( y_value )
        assert 0.0 < this_prob

        this_min = min([specific_info_with_Y_equals_y( [0,1], y_value ), specific_info_with_Y_equals_y( [2], y_value ) ])

        summ += this_prob * this_min

    print "* {0,12}={1,02}={2,01}=%s" % summ
										  

def Imax_Synergy( Xs, norm=False ):
    '''returns the Imax-synergy among the Xs about target Y.  S_1(Xs:Y) = S(Xs:Y)'''

    check_Xs( Xs )

    summ = 0.0
    for y_value in Y_states():
        this_prob = Prob_Y_equals_y( y_value )
        assert 0.0 < this_prob

        this_infos = [ specific_info_with_Y_equals_y( [Xindex], y_value ) for Xindex in Xs ]
        this_maxbelow = max(this_infos)

        summ += this_prob * this_maxbelow

    z = mutual_info_with_Y( Xs ) - summ

    # if we're normalizing do that now
    if norm:
        z /= mutual_info_with_Y( Xs )
        z = round(z,PRECISION)
        assert 0.0 <= z <= 1.0
        return z


    if z < 0.0 and round(z,PRECISION) == 0.0:
        z = 0.0

    assert 0.0 <= z <= mutual_info_with_Y( Xs ), "Imax_Synergy=%s but I(Xs:Y)=%s" % (z, mutual_info_with_Y(Xs))

    z = round(z,PRECISION)

    return z

def Cost_of_Removal( Xs ):
    '''returns the "cost of removal" holism measure.'''

    check_Xs( Xs )

    max_below = 0.0
    for i, Xindex in enumerate(Xs):
        this_below_Xs = deepcopy(Xs)
        temp = this_below_Xs.pop(i)
        assert temp == Xindex
        max_below = max( max_below, mutual_info_with_Y(this_below_Xs) )

    z = mutual_info_with_Y( Xs ) - max_below

    if z < 0.0 and round(z2,PRECISION) == 0.0:
        z = 0.0

    assert 0.0 <= z

    z = round(z,PRECISION)

    return z

def Holism( Xs, norm=False ):
    '''returns the virgil-holism among the Xs about target Y.  HOLISM( Xs: Y )'''

    check_Xs( Xs )

    #    z1 = mutual_info_with_Y( Xs )
    #
    #    max_below = 0.0
    #    for i, Xindex in enumerate(Xs):
    #        this_below_Xs = deepcopy(Xs)
    #        temp = this_below_Xs.pop(i)
    #        assert temp == Xindex
    #        max_below = max( max_below, mutual_info_with_Y(this_below_Xs) )

    #    z1 -= max_below

    #    if norm:
    #        z1 /= mutual_info_with_Y( Xs )
    #        z1 = round(z1, PRECISION)
    #        assert 0.0 <= z1 <= 1.0
    #        return z1

    ##########################################################
    ## NOW TRY THE STATE-DEPENDENT WAY
    ##########################################################

    summ = 0.0
    for y_value in Y_states():
        this_prob = Prob_Y_equals_y( y_value )
        assert 0.0 < this_prob

        this_infos = []
        for i, Xindex in enumerate(Xs):
            this_below_Xs = deepcopy(Xs)
            temp = this_below_Xs.pop(i)
            assert temp == Xindex
            this_infos.append( specific_info_with_Y_equals_y( this_below_Xs, y_value ) )

        this_maxbelow = max(this_infos)

        summ += this_prob * this_maxbelow

    z2 = mutual_info_with_Y( Xs ) - summ

    if norm:
        z2 /= mutual_info_with_Y( Xs )
        z2 = round(z2, PRECISION)
        assert 0.0 <= z2 <= 1.0
        return z2

    if z2 >= Imax_Synergy(Xs)+eps:
        print "z=%s \t S(Xs:Y)=%s" % ( z2, Imax_Synergy(Xs) )

    if z2 < 0.0 and round(z2,PRECISION) == 0.0:
        z2 = 0.0

    assert 0.0 <= z2 <= Imax_Synergy( Xs )+eps, "z2=%s but Imax_Synergy=%s" % (z2, Imax_Synergy(Xs))

    z2 = round(z2,PRECISION)

    return z2

def check_collection( collection ):
    '''returns TRUE if the collection is valid'''

    for i, coalition in enumerate(collection):
        assert type(coalition) is list, "Coalition %s was not a list" % i
        assert check_Xs( coalition )

    return True

def check_y_value( y_value ):
    '''check that the Y_value is valid'''

    assert type(y_value) is str
    assert y_value in Y_states()

    return True

def check_Conds( Conds ):
    '''A ``Cond'' is a joint-predictor (Xs).  So a list of Conds is a list of joint-predictors.  So we check each joint predictor. '''

    assert type(Conds) is list or type(Conds) is tuple

    # make sure each Cond is a valid joint-predictor
    for Cond in Conds:
        check_Xs( Cond )

    return True

def check_Xs( Xs ):
    '''do the various checks to ensure this list of X indices is valid'''

    assert type(Xs) is not int

    if type(Xs) is list:
        Xs = tuple(Xs)

    #    if type(Xs) is int:
    #        Xs = [Xs]

    #    print "Xindices=%s" % str(Xs)
    #    assert type(Xs) is list or type(Xs) is tuple
    assert type(Xs) is tuple

    # uniquify the list
    #print "Xs=%s" % Xs
    Xs = list(set(Xs))
    assert 1 <= len(Xs) <= 4

    for X in Xs:
        assert X in get_Xindices()

    return True


def get_reducible_Xstates( inputTABLE ):
    '''returns the states that can be reduced per Gács-Körner..  If none can be reduced, returns (None, None)

    Two states 'x' and 'w' can be REDUCED if and only if p(x) = p(w) = p(x,w)
    '''

    global TABLE
    oldTABLE = TABLE

    TABLE = inputTABLE

    assert get_numXs() == 2, "get_numXs=%s" % get_numXs()

    Xstates = joint_states([0], True )
    Wstates = joint_states([1], True )

    z = []
    for Xstate in Xstates:
        assert len(z) == 0, "z=%s" % z
        for Wstate in Wstates:

            # if the conditional probability is 0.0 or 1.0, then skip it.
            if Prob_xs_given_condxs( Wstate.indices, Wstate.states, Xstate.indices, Xstate.states ) in [0.0, 1.0]:
                continue

            #print "-- Discovered that p(x2=%s | x1=%s) = %s" % ( Wstate.states, Xstate.states, Prob_xs_given_condxs( Wstate.indices, Wstate.states, Xstate.indices, Xstate.states ) )
            z.append( Wstate )

            if len(z) == 2:
                TABLE = oldTABLE
                return z

        assert len(z) == 0, "z=%s" % z

    # now do the same thing the other way...
    for Wstate in Wstates:
        assert len(z) == 0
        for Xstate in Xstates:

            # if the conditional probability is 0.0 or 1.0, then skip it.
            if Prob_xs_given_condxs( Xstate.indices, Xstate.states, Wstate.indices, Wstate.states ) in [0.0, 1.0]:
                continue

            #            print "-- Discovered that p(x1=%s | x2=%s) = %s" % ( Xstate.states, Wstate.states, Prob_xs_given_condxs( Xstate.indices, Xstate.states, Wstate.indices, Wstate.states ) )
            z.append( Xstate )

            if len(z) == 2:
                TABLE = oldTABLE
                return z

        assert len(z) == 0


    TABLE = oldTABLE

    return (None, None)


def check_xs( Xindices, Xstates, require_positive_prob=True ):
    '''do the various checks to ensure this list of X indices is valid'''

    assert type(Xindices) is tuple or type(Xindices) is list
    assert type(Xstates) is tuple or type(Xstates) is list

    #    print "Xindices=%s" % str(Xindices)

    # first check the indices
    check_Xs( Xindices )

    if type(Xstates) is str:
        Xstates = (Xstates,)
    elif type(Xstates) is list:
        Xstates = tuple(Xstates)

    #    assert type(Xindices) is tuple

    # assert both are the same length
    assert len(Xindices) == len(Xstates)

    # assert that each state exists as a singleton
    for x_index, x_state in zip(Xindices,Xstates):
        assert x_state in X_values(x_index)


    if require_positive_prob:
        # assert that the JOINT STATE of them all exists
        # now check the state
        allstates = joint_states( Xindices )

        if Xstates not in allstates:
            print "Xindices=%s \t Xstates=%s" % (Xindices, Xstates)
            print "allstates=%s" % allstates

        assert Xstates in allstates, "check_xs(): Xstate=%s wasn't found in the specified Xindices" % Xstates

    return True


def Imin_coalition_Y( collection ):
    '''calculates the *average* minimum specific-information over each coalition within the collection, I_min( coalition: Y )'''

    check_collection( collection )

    z = 0.0

    for y_value in Y_states():
        this_prob = Prob_Y_equals_y( y_value )
        this_Imin = Imin_coalition_Y_equals_y( collection, y_value )
        z += this_prob * this_Imin

    assert 0.0 <= z
    return z

def Imin_coalition_Y_equals_y( collection, y_value ):
    '''calculates the minimum specific-information over each coalition within the collection, min I(coalition:Y=y)'''

    check_y_value( y_value )
    check_collection( collection )


    z = min( [ specific_info_with_Y_equals_y(coalition, y_value) for coalition in collection ] )

    #z = z.round(PRECISION)

    assert 0.0 <= z
    return z



def specific_info_with_Y_equals_y( Xs, y_value ):
    '''computes the specific information I(Xs:Y=y) using the Kullback-Liebler Divergence'''

    check_Xs( Xs )
    check_y_value( y_value )

    #print "In I(Xs:Y=y)!",

    z = 0.0
    Xstates = joint_states(Xs)

    for Xstate in Xstates:

        term1 = Prob_Xs_given_y(Xs,Xstate,y_value)

        if term1 == 0.0:
            continue

        top, bottom = Prob_Xs_given_y(Xs,Xstate,y_value), Prob_Xs(Xs,Xstate) # Prob_Y_equals_y(y_value)

        term2 = log2( top / bottom )

        z += term1 * term2

    z = round(z, PRECISION+1)
    #z = round(z,PRECISION)

    #print "z=%s" % z
    assert 0.0 <= z

    return z

def mutual_info_between_Xs( firstXs, secondXs ):
    '''Calculates I(first_Xs, second_Xs) = H(first) + H(second) - H(first,second)'''

    check_Xs( firstXs )
    check_Xs( secondXs )
    check_Xs( firstXs+secondXs )

    z = H_Xs(firstXs) + H_Xs(secondXs) - H_Xs(firstXs+secondXs)

    assert 0.0 <= z

    assert z <= min([H_Xs(firstXs),H_Xs(secondXs)])+eps

    return z

def mutual_info_with_Y( Xs ):
    '''Calculates I(Xs:Y) = H(Xs) + H(Y) - H(Xs,Y)'''
    check_Xs( Xs )

    z1 = H_Xs(Xs) + H_Y() - H_Xs_Y(Xs)

    # if z1 is just ever-so-barely negative, make it zero.
    if z1 < 0 and abs(z1) < epps:
        z1 = 0

    if min( H_Xs(Xs), H_Y() )+eps < z1 or z1 < 0:
        print "min(H_Xs(Xs),H_Y())=%s < z1=%s" % ( min( H_Xs(Xs), H_Y() ), z1 )


    assert 0.0 <= z1 <= min( H_Xs(Xs), H_Y() )+eps

    # Now calculate this a second way
    z2 = 0.0
    for y_value in Y_states():
        z2 += Prob_Y_equals_y( y_value ) * specific_info_with_Y_equals_y( Xs, y_value )
    z2 = round(z2,PRECISION)

    # Now calculate this a THIRD way
    z3 = H_Xs(Xs) - H_Xs_given_Y( Xs )
    z3 = round(z3,PRECISION)
    assert abs(z1-z3) <= 2*eps, "z1=%s \t z3=%s The two methods for calculating mutual information didn't come out equal" % (z1,z3)

    #    print "z1=%s   \t   z2=%s    \t   z3=%s" % (z1,z2,z3)

    return z1

def H_Xs_given_Y( Xs ):
    '''returns the conditional entropy H(Xs|Y) = H(Xs,Y) - H(Y)
    '''
    check_Xs(Xs)

    z = H_Xs_Y(Xs) - H_Y()

    if z < 0.0 or z > H_Xs(Xs)+eps:
        print "z=%s \t H_Xs()=%s" % (z, H_Xs(Xs) )

    assert 0.0 <= z <= H_Xs(Xs)+eps, "Conditional entropy H(Xs|Y) didn't obey the bounds"

    return z


def H_Xs_given_Y_equals_y( Xs, ystate ):
    '''returns the conditional entropy H(Xs|Y=y) = \sum_x p(x|y) log[ 1/p(x|y) ]'''
    check_Xs(Xs)
    check_y_value( ystate )

    z = 0.0
    Xstates = joint_states(Xs)

    for Xstate in Xstates:
        term1 = Prob_Xs_given_y(Xs, Xstate, ystate)
        if term1 == 0.0:
            continue

        term2 = log2( 1.0 / term1 )
        z += term1 * term2


    assert 0.0 <= z

    return z


def total_correlation_given_Y_equals_y( Xs, ystate ):
    '''returns the CONDITIONAL TOTAL CORRELATION, TC(Xs|Y) = [ sum_i H(X_i|Y) ] - H(Xs|Y)
    '''

    check_Xs(Xs)
    check_y_value( ystate )

    z = sum([ H_Xs_given_Y_equals_y([index],ystate) for index in Xs])

    z -= H_Xs_given_Y_equals_y( Xs, ystate )

    assert 0.0 <= z

    z = round(z,PRECISION+1)

    return z


def total_correlation_given_Y( Xs ):
    '''returns the CONDITIONAL TOTAL CORRELATION, TC(Xs|Y) = [ sum_i H(X_i|Y) ] - H(Xs|Y)
    '''

    check_Xs(Xs)
    z = sum([ H_Xs_given_Y([index]) for index in Xs])

    z -= H_Xs_given_Y(Xs)

    assert 0.0 <= z

    return z

def total_correlation(Xs):
    '''returns the TOTAL CORRELATION among the X terms'''

    check_Xs(Xs)

    z = sum([ H_Xs([index]) for index in Xs])
    z -= H_Xs( Xs )

    if z < 0.0 and abs(z) <= eps:
        z = 0

    assert 0.0 <= z

    return z

def joint_states( Xs, useclass=False ):
    '''returns the joint states of the passed X indices'''

    check_Xs(Xs)

    z = None
    if len(Xs) == 1:
        z = [ (x_1,) for x_1 in X_values(Xs[0]) if Prob_Xs(Xs,[x_1]) ]
    elif len(Xs) == 2:
        z = [ (x_1,x_2) for x_1 in X_values(Xs[0]) for x_2 in X_values(Xs[1]) if Prob_Xs(Xs,[x_1,x_2]) ]
    elif len(Xs) == 3:
        z = [ (x_1,x_2,x_3) for x_1 in X_values(Xs[0]) for x_2 in X_values(Xs[1]) for x_3 in X_values(Xs[2]) if Prob_Xs(Xs,[x_1,x_2,x_3]) ]
    elif len(Xs) == 4:
        z = [ (x_1,x_2,x_3,x_4) for x_1 in X_values(Xs[0]) for x_2 in X_values(Xs[1]) for x_3 in X_values(Xs[2]) for x_4 in X_values(Xs[3]) if Prob_Xs(Xs,[x_1,x_2,x_3,x_4]) ]
    else:
        assert 0 == 1, "Should be impossible to get here."


    if useclass == False:
        return z


    # if the user has passed useclass=True, then we convert each entry of Z to an xjointclass.
    z2 = []

    for state in z:
        to_append = Xjointstate( Xs, state )
        z2.append( to_append )

    return z2


def joint_states_with_Y( Xs, useclass=False ):
    '''returns the joint states of the passed X indices with Y '''

    check_Xs(Xs)


    z = None
    if len(Xs) == 1:
        z = [ ((x_1,),y) for x_1 in X_values(Xs[0]) for y in Y_states() if Prob_Xs_y(Xs,(x_1,), y) ]
    elif len(Xs) == 2:
        z = [ ((x_1,x_2),y) for x_1 in X_values(Xs[0]) for x_2 in X_values(Xs[1]) for y in Y_states() if Prob_Xs_y(Xs,(x_1,x_2),y) ]
    elif len(Xs) == 3:
        z = [ ((x_1,x_2,x_3),y) for x_1 in X_values(Xs[0]) for x_2 in X_values(Xs[1]) for x_3 in X_values(Xs[2]) for y in Y_states() if Prob_Xs_y(Xs,(x_1,x_2,x_3),y) ]
    elif len(Xs) == 4:
        z = [ ((x_1,x_2,x_3,x_4),y) for x_1 in X_values(Xs[0]) for x_2 in X_values(Xs[1]) for x_3 in X_values(Xs[2]) for x_4 in X_values(Xs[3]) for y in Y_states() if Prob_Xs_y(Xs,(x_1,x_2,x_3,x_4),y) ]
    else:
        assert 0 == 1, "Should be impossible to get here."

        #for i, (Xstate, Ystate) in enumerate(z):
        #print "%s) Xstate=%s \t Ystate=%s \t prob=%s" % ( i, Xstate, Ystate, Prob_Xs_y(Xs,Xstate,Ystate) )

    if useclass == False:
        return z

    # if the user has passed useclass=True, then we convert each entry of Z to an Xjointstate.
    z2 = [ Xjointstate(Xs,state) for state in z ]
    return z2


def H_Xs_Y( Xs ):
    '''calculate the joint entropy of the set of Xs with Y.  H(Xs,Y)'''

    ## Check input
    check_Xs(Xs)
    #################################

    allstates = joint_states_with_Y( Xs )

    allprobs = [ Prob_Xs_y(Xs,Xstates,Ystate) for Xstates,Ystate in allstates ]

    # remove all probabilities that are 0.0
    allprobs = [ p for p in allprobs if p ]

    assert (sum(allprobs)-1.0) <= eps

    z = sum([ p * log2(1.0/p) for p in allprobs ])

    #    z = round(z,PRECISION)

    assert 0.0 <= z

    return z


def H_Xs( Xs ):
    '''calculate the entropy of the set of Xs'''

    ## Check input
    check_Xs(Xs)
    #################################

    allstates = joint_states( Xs )

    allprobs = [ Prob_Xs(Xs,this_state) for this_state in allstates ]

    # remove all probabilities that are 0.0
    allprobs = [ p for p in allprobs if p ]

    #    summ = sum(allprobs)
    #    if allprobs != 1.0:
    #        print "Probs=%s" % allprobs
    #        print "summ=%s" % summ

    assert (sum(allprobs)-1.0) < eps

    z = sum([ p * log2(1.0/p) for p in allprobs ])

    #    z = round(z,PRECISION)
    assert 0.0 <= z

    return z

def H_Y():
    '''calculate the entropy of the Ys'''

    ProbYs = [ Prob_Y_equals_y(x) for x in Y_states() ]

    summ = sum(ProbYs)
    #    if summ != 1.0:
    #        print "p(Y)=%s" % ProbYs
    #        print "summ=%s" % summ

    assert (sum(ProbYs)-1.0) <= eps

    #    z = 0.0
    z = sum([ prob_y * log2(1.0/prob_y) for prob_y in ProbYs])
    #    for prob_y in ProbYs:
    #        z += (prob_y * log2( 1.0 / prob_y ))
    #    z = round(z,PRECISION)
    assert 0.0 <= z

    return z

def Prob_Xs_y( Xs, Xvalues, y_value ):
    "Compute the joint probability of all of these Xs occurring along with Y=y"
    # sanity check our inputs

    #    print "* INCOMING: Xs=%s X_values=%s  \t Y=%s" % ( Xs, Xvalues, y_value )

    check_xs( Xs, Xvalues, False )
    check_y_value( y_value )
    #    if type(Xs) is int and type(Xvalues) is str:
    #        Xs, Xvalues = [Xs], [Xvalues]

    #    for X, Xval in zip(Xs,Xvalues):
    #        assert 0 <= X < get_numXs()
    #        assert Xval in X_values(X)

    #    check_Xs(Xs)

    global TABLE
    num_instances = 0

    #print "\t num_instances=%s" % num_instances,
    #print "Xs=%s Xvalues=%s y=%s" % (Xs, Xvalues, y_value)

    z = 0.0

    uniform_weight = 1.0 / float(len(TABLE))

    #    print "uniform_weight=%s" % uniform_weight

    # for each input string in TABLE...
    for tableX, tableY in TABLE.iteritems():

        tableXs = tableX.split(' ')
        assert type(tableY) in [str, list]

        if type(tableY) is str:
            tableY = [tableY]

        # do any y_values match?  If not, skip.
        if not tableY.count(y_value):
            continue

        # Do all of the X-values match?
        match=True
        for X, Xval in zip(Xs,Xvalues):
            if tableXs[X] != Xval:
                match = False
                break

        if not match:
            continue

        #        print "Xs=%s \t Xvalues=%s \t y_value=%s" % (Xs, Xvalues, y_value)
        #        print "tableX=%s \t tableY=%s" % (tableXs, tableY)

        contrib = (tableY.count(y_value) / float(len(tableY)))
        #        print "contrib=%s" % contrib
        # increment z by the amount of y_values that match, weighted by the uniform weight
        z += uniform_weight * contrib

    assert 0.0 <= z <= 1.0, "Pr( Xs = xs ) wasn't between (0,1]."

    assert z <= Prob_Xs( Xs, Xvalues ), "p(Xs)=%s  \t p(Xs,y)=%s" % ( Prob_Xs(Xs,Xvalues), z )
    assert z <= Prob_Y_equals_y( y_value )

    return z

def Prob_Xs( Xs, Xvalues ):
    "Compute the joint probability of all of these Xs occurring"
    # sanity check our inputs

    check_Xs( Xs )
    #    print Xs
    #    print Xvalues
    assert len(Xs) == len(Xvalues)

    for X, Xval in zip(Xs,Xvalues):
        assert 0 <= X < get_numXs()
        assert Xval in X_values(X)


    global TABLE
    num_instances = 0

    #    print "Xs=%s  Xvalues=%s" % ( Xs, Xvalues)

    #    print "TABLE=%s" % TABLE.keys()
    # for each input string in TABLE...
    for keystring in TABLE.keys():
        values = keystring.split(' ')

        #        print "to_match values=%s" % str(values)

        match=True
        # If that string matches, num_instances += 1
        for X, Xval in zip(Xs,Xvalues):
            if values[X] != Xval:
                match = False

        if match:
            num_instances += 1
        #            print "Matched!  num_instances=%s" % num_instances
        #            raw_input('....')

    z = float(num_instances) / float(len(TABLE))

    assert 0.0 <= z <= 1.0, "Pr( Xs = xs ) wasn't between [0,1]."
    return z

def Prob_y_given_xs( y_value, Xs, Xvalues ):
    '''returns p(Y=y|Xs=xs) = p(Xs=xs,Y=y) / p(Xs=xs)'''

    check_y_value( y_value )
    check_xs( Xs, Xvalues, False )

    top = Prob_Xs_y( Xs, Xvalues, y_value )

    if top == 0.0:
        return 0.0

    bottom = Prob_Xs( Xs, Xvalues )

    #    assert top <= bottom

    z = top / bottom

    assert 0.0 <= z <= 1.0, "Not a valid probability! Eep! y_value=%s \t Xs=%s \t Xvalues=%s \n top=%s \t bottom=%s \t z=%s" % (y_value, Xs, Xvalues, top, bottom, z)

    return z


def Prob_Xs_given_y( Xs, Xvalues, y_value ):
    '''returns p(Xs=xs|Y=y) = p(Xs=xs,Y=y) / p(Y=y)'''

    check_xs( Xs, Xvalues, False )
    check_y_value( y_value )

    top = Prob_Xs_y(Xs, Xvalues, y_value )

    if top == 0.0:
        return 0.0

    bottom = Prob_Y_equals_y(y_value)

    z = top / bottom

    assert top <= bottom

    assert 0.0 <= z <= 1.0, "Not a valid probability! Eep! Xs=%s \t Xvalues=%s \t y_value=%s \n top=%s \t bottom=%s \t z=%s" % (Xs, Xvalues, y_value, top,bottom,z)

    return z

def genallsubsets(n):   #generator version.
    S = [[i] for i in range(n)] #initial.
    for s in S:
        yield s
    N = S
    for k in range(n):
        if N:
            S = N
            N = []
            for e in S:
                for i in range(e[-1]+1,n):
                    f = e[:]
                    f.append(i)
                    yield f
                    N.append(f)

def partitions(set_):
    if not set_:
        yield []
        return

    for i in xrange(2**len(set_)/2):
        parts = [set(), set()]
        for item in set_:
            parts[i&1].add(item)
            i >>= 1
        for b in partitions(parts[1]):
            yield [ list(x) for x in [parts[0]]+b ]


def Y_states():
    '''Returns the sorted list of Y values.'''
    global TABLE

    z = []

    # add all str values to z
    z = [ x for x in TABLE.values() if type(x) is str ]

    # make a list of the entries that aren't strings...
    non_strs = [ x for x in TABLE.values() if not (type(x) is str) ]

    # add each list of states to z.
    for entry in non_strs:
        z.extend( entry )

    # remove duplicate from it, sort it, return it.
    return sorted(list(set(z)))

def Prob_Y_equals_y( y_value ):
    assert y_value in Y_states()
    global TABLE

    # For each input state in the table, determine the amount of probability this y_value gets
    num_entries = float(len(TABLE))

    z = 0.0
    for entry in TABLE.values():

        # if the entry is a string, it gets the full contribution of 1.0 / num_entries
        if type(entry) is str:
            if entry == y_value:
                z += 1.0 / num_entries
                continue

        # if the entry is a LIST, it gets 1.0/num_entries multiplied by it's proportion among the outcomes
        elif type(entry) is list:
            z += 1.0 / num_entries * float(entry.count(y_value)) / float(len(entry))
            continue

        else:
            assert 0 == 1, "should never get here"


    assert 0.0 < z <= 1.0, "Pr(Y = y ) wasn't between (0,1]."
    return z



def X_values( X_index ):
    '''Returns the sorted list of X_i values.'''
    assert 0 <= X_index < get_numXs()
    global TABLE

    z = [ key.split()[X_index] for key in TABLE ]
    #    z = []
    #    for key in TABLE:
    #        values = key.split()
    #        z.append(values[X_index])

    #    return sorted(list(set(z)))
    # !!!!!!!
    #return list(set(z))
    return sorted(list(set(z)))



def get_numXs():
    global TABLE

    numX = None
    for key in TABLE:
        Xs = key.split(' ')
        if numX is None:
            numX = len(Xs)
        assert numX == len(Xs)

    return numX

def get_Xindices():
    return range( get_numXs() )

if __name__ == '__main__':
#    pprint(TABLE)
    numXs = get_numXs()
    allXs = get_Xindices()

    # if n=2, run the tests for it.
    #    PI_tests_n2()
    #    PI_tests_n3()


    #    print numXs
    print term.MAGENTA,
    print "\rPossible X-states:"
    for index in allXs:
        this_values = X_values(index)

        print "\tX_%s = {%s} " % (index, str(this_values)[1:-1] )

        this_probs = [ Prob_Xs([index], [x]) for x in this_values ]

        for val, prob in zip(this_values,this_probs):
            print "\t   Prob(X_%s = %s ) = %s" % ( index, val, prob )

    print term.NORMAL,
    print "\r=================================================="

    print term.YELLOW,
    print "\rPossible Y-states: {%s}" % ( str(Y_states())[1:-1] )

    Ystates = Y_states()
    Prob_Ys = [ Prob_Y_equals_y(y) for y in Ystates ]

    #    print "allXs=%s" % allXs

    for val, prob in zip(Ystates, Prob_Ys):
        # get the lowest and highest possible synergy for this state Y=y
        this_Isyn_min, this_Isyn_max = min_max_synergy_y( allXs, val )
        #        this_Isyn = Delta_I_y(allXs, val)

        print "\t y=%s \t Prob(Y=%s) = %s" % (val, val, prob)

        print "\t I(X1,X2:Y=%s) = %s" % ( val, specific_info_with_Y_equals_y(allXs, val) )
        print "\t I(X1:Y=%s)=%s" % ( val, specific_info_with_Y_equals_y( [0], val ) )
        print "\t I(X2:Y=%s)=%s" % ( val, specific_info_with_Y_equals_y( [1], val ) )
        print "\t I(X1:X2|Y=%s)=%s" % ( val, total_correlation_given_Y_equals_y([0,1], val) )
        print ''

        print "\t\t D=%s" % D( allXs, val )
        #        print u"\t\t ΔI(X;Y=%s) = %s" % ( val, this_Isyn )

        #        if not (this_Isyn_min <= this_Isyn <= this_Isyn_max ):
        #            print term.RED,
        #        print u"\t\t min=%s ≤ %s ≤ max=%s" % (this_Isyn_min, this_Isyn, this_Isyn_max),
        #        if not (this_Isyn_min <= this_Isyn <= this_Isyn_max ):
        #            print term.YELLOW,

        print '\n'
    print term.NORMAL,

    # Now print the AVERAGE minimum and maximum synergy
    Isyn_min, Isyn_max = min_max_synergy_Y( allXs )
    Isyn = Delta_I( allXs )

    if not (Isyn_min <= Isyn <= Isyn_max ):
        print term.RED,
    print u"\rAverage ΔI: \t min=%s ≤ %s ≤ max=%s" % (Isyn_min, Isyn, Isyn_max),
    if not (Isyn_min <= Isyn <= Isyn_max ):
        print term.NORMAL,

    print ""
    print "=================================================="
    print "allXs=%s" % allXs

    print "H(Y)=%s" % H_Y()
    print "H(Xs)=%s" % H_Xs(allXs),
    print "\t TC(Xs)=%s \t TC(Xs|Y)=%s" % (total_correlation(allXs), total_correlation_given_Y(allXs))

    for i, (Xstate, Ystate) in enumerate(joint_states_with_Y( allXs )):
        print "%s) Xstate=%s \t Ystate=%s \t prob=%s" % ( i, Xstate, Ystate, Prob_Xs_y(allXs,Xstate,Ystate) )

    print "H(Xs,Y)=%s" % H_Xs_Y( allXs )
    print "I(Xs:Y)=%s" % mutual_info_with_Y( allXs )
    print "H(Y|Xs)=%s" % H_Y_given_Xs( allXs )
    for X_index in range(numXs):
        print "\t H(X%s)=%s \t I(X%s:Y)=%s" % ( X_index, H_Xs([X_index]), X_index, mutual_info_with_Y([X_index]) )



    print "I(X0:X1|Y) = %s" % total_correlation_given_Y( [0,1] )
    print "I(X0:X1) = %s" % float(H_Xs([0]) + H_Xs([1]) - H_Xs([0,1]))

    MI_X1_Y_cond_X2 = mutual_info_with_Y_condXs( [0], [1] )
    MI_X2_Y_cond_X1 = mutual_info_with_Y_condXs( [1], [0] )

    #print "I(X0:Y)=%s" % MI
    #print u"I(X0:Y↓X1)=%s" % MI_X1_Y_lowercond_X2
    print "I(X0:Y|X1)=%s" % MI_X1_Y_cond_X2
    print "I(X1:Y|X0)=%s" % MI_X2_Y_cond_X1


    PRINT_SYNERGY_MEASURES( allXs )



#    print "I(X0:Y \downarrow X1 )=%s" % I_X_Y_downarrow_Z( [0], [1] )

print "Finished!  Yay!"