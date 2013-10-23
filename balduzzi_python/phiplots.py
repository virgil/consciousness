#!/usr/bin/python
# Examples

# !!! Virgil
# generic modules
import matplotlib.pyplot as pyplot
import matplotlib.pylab as pylab
import matplotlib.colors as pycolors
import numpy as npy
import sys, math, logging
import pygraphviz as pgv
from math import ceil
from pprint import pprint
from os import getpid
from scipy.stats import variation
from scipy import mean
from numpy import array, matrix

# define the log2 function
log2 = lambda x: math.log(x, 2)

################################################
# Configuration
logging.basicConfig(level=logging.WARNING,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename='/tmp/%s.log' % getpid(),
                    filemode='w')

################################################
# NOTE: call PHI.SP and PHI.SW before doing anything

def dec2label( num, label, pad_to=None ):
    """Converts a decimal number to a label for plotting.  Has an option 'pad_to' variable that specifies the padding."""

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
        n, prepend = label, label + 'b'
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




def MIPlist_2_MIPdict( MIPs, show_invalidstates=True ):
    """Returns the MIPs in 'dictionary form' of label -> array .  Also returns the maximum number of parts"""

    numstates = len( MIPs )
    
    # join all of the arrays into into a single array, then count the number of element in it.
    # determine the number of units
    for row in MIPs:
        if len(row):
            numunits = len(npy.concatenate(row))
            break

    print "numunits=%s \t numstates=%s" % ( numunits, numstates )


    # initial max number of parts is zero.  Will be increased
    maxnum_nonemptyparts = 0
    # all units are assumed to not be part of any part ( -1 )
    MIPmatrix = npy.ones( (numstates, numunits) ) * -1

#    pprint(MIPmatrix)
    
    # for each state...
    for stateindex, MIP in enumerate(MIPs):
#        print "stateindex=%s" % stateindex
        numparts = len(MIP)        
        num_nonemptyparts = len([x for x in MIP if x != [] ])
        
        if num_nonemptyparts > maxnum_nonemptyparts:
            maxnum_nonemptyparts = num_nonemptyparts
        
        # for each unit...
        for unit in range(numunits):
#            print "\tunit=%s" % unit

            # look across all of the parts until we find the part that is unit is in
            for partnum, part in enumerate(MIP):
                if unit in part:
#                    print "\t\tpartnum=%s\tpart=%s" % (partnum, part)
                    MIPmatrix[stateindex][unit] = partnum
#                    break
#                else:
#                    print "\t\t%s not in part %s: %s" % ( unit, partnum, part )
    
    # !!! HACK, convert the matrix to a dictionary and return that.
    MIPdict={}
    
    for index, row in enumerate(MIPmatrix):
#        print "index=%s  row=%s" % (index, row )
 #       print "row=%s  sum(row)=%s" % (row, sum(row))
        # if this row is an unreachable state...
        if sum(row) >= 0.0:
            MIPdict[index] = row
    
        # add an invalid_state if we're doing that.
        elif show_invalidstates and sum(row) < 0:
            MIPdict[index] = row
            
    
    
    return MIPdict, maxnum_nonemptyparts


def HD( v1, v2 ):
    """Calculates the hamming distance between two vectors"""

    if len(v1) != len(v2):
        raise ValueError, "vectors v1 and v2 must be of the same length."

    hamdist = abs( len(v1) - len(v2) )
    
    for i in range( min( len(v1),len(v2)) ):
        if v1[i] != v2[i]:
            hamdist += 1

    return hamdist

def plot_MIPs_for_allstates( states_and_MIPs, title=None, filename=None, label='hex', averageMIP=None ):
    """prints a matrix plot of the MIP for each state"""

    MY_COLORMAP_COLORS = ['r','g','b','y','m','c','w']
    numunits = len(npy.concatenate(states_and_MIPs[0][1])) # get the length of the first array 
    maxnum_nonemptyparts = max([ len(MIP) for x1, MIP in states_and_MIPs] )
    print "numunits=%s" % numunits
    
    # if we have averageMIP, append it to rest so it goes through the same processing.
    if averageMIP:
        states_and_MIPs.append( averageMIP )
    
    # convert the tuple MIP to a matrix-suitable row
    new_states_and_MIPs = []
    for state, MIP in states_and_MIPs:
        newMIP = npy.ones( (numunits,1) ).flatten() * -1  # initially assume all parts are -1s
        for partnumber, indices in enumerate(MIP):
            newMIP[indices] = partnumber
        new_states_and_MIPs.append( (state,newMIP) )

#    print "newMIPs=" ,
#    pprint( new_states_and_MIPs )

    # Make a list of (ylabel, MIP) from MIPdict
    ylabel_and_MIPs = new_states_and_MIPs
    

    # sort the rows of the MIPmatrix by haming distance.
    # 1. Make new datastructure of [  (hammingdistance,row) ]
    # 2. Sort by hamingdistance
    # 3. put all the rows back into the MIPmatrix

    # if there is an averageMIP, put that one top.  Otherwise, just use the first one.
    if averageMIP:
        firstrow = ylabel_and_MIPs[-1][1]
    else:
        firstrow = ylabel_and_MIPs[0][1]
    hamdist_and_ylabel_rows = []

    # invert the firstrow if we need to.
    if firstrow[0] == 0 and max(firstrow)==1 and min(firstrow)==0:
        firstrow = abs(firstrow-1)
    
    for ylabel, row in ylabel_and_MIPs:

        # Make the first node always the same color.  (Assuming they are only bipartitions)
        if row[0] == 0 and max(row)==1 and min(row)==0:
#            print "oldrow=%s" % row
#            print "negrow=%s" % abs(row-1)
            row = abs(row-1)
        
        hamdist_and_ylabel_rows.append( [HD(firstrow, row), ylabel, list(row)] )
        
#        print "\tfirstrow=%s \t row=%s   HD=%s" % (firstrow, row, HD(firstrow,row) )

#    pprint( hamdist_and_ylabel_rows )
    hamdist_and_ylabel_rows.sort()
    # Force the average MIP to be on top, if we have one
#    pprint( hamdist_and_ylabel_rows )

    # overwrite the old ylabel_and_MIPs
    ylabel_and_MIPs = [ (ylabel,array(row)) for hamdist, ylabel, row in hamdist_and_ylabel_rows ]

    # if we do have an average MIP, FORCE it to be on top.
    if averageMIP:
        converted_averageMIP = [ (ylabel, row) for ylabel, row in ylabel_and_MIPs if ylabel == averageMIP[0] ]
        assert len(converted_averageMIP) == 1, "should have only had 1 entry with the average ylabel!"
        converted_averageMIP = converted_averageMIP[0]       # there will only be 1 item in this list
        ylabel_and_MIPs.remove( converted_averageMIP )    # remove the old one in whatever location it's in
        ylabel_and_MIPs.insert(0, converted_averageMIP )  # and re-insert it at position zero.

    # create the MIPmatrix from the ylabel_and_MIPs
    MIPmatrix = None

    for ylabel, row in ylabel_and_MIPs:
        # this is a weird construction with append being a function instead of numpy instead of property of array, but whatever.

        if MIPmatrix is None:
            MIPmatrix = row
        else:
            MIPmatrix = npy.vstack( [MIPmatrix, row] )


#        pprint( MIPmatrix )
#        MIPmatrix = npy.append( MIPmatrix, row )


#    pprint( MIPdict )
#    pprint( MIPmatrix )

    # set the size of the figure -- default size is (8,6)
    figuresize = ( max(8,numunits/2.0), max(6,len(MIPmatrix)/3.0) )
    fig = pyplot.figure(figsize=figuresize )

    my_cmap = pycolors.ListedColormap( MY_COLORMAP_COLORS, N=maxnum_nonemptyparts)
    
    pyplot.pcolor( MIPmatrix[ ::-1,:], cmap=my_cmap, shading='faceted' )
#    pyplot.matshow( MIPmatrix, cmap=my_cmap, alpha=1.0, origin='upper' )

    # if we have a title, print it.
    if title:
        pyplot.title( title + ' r vs g')

    # plot the labels
    pyplot.xlabel('Unit')
    pyplot.ylabel('Network State')


    # set the x and y ticks (the units)
    xtick_labels = [ int(x) for x in pyplot.xticks()[0] if x % 1 == 0 and x >= 1]  # labels [1...n]
    # we're zero based.  So pop off the last one.
    xticks = [ x-0.5 for x in xtick_labels ] # positions [ 0.5 ... n+0.5 ]
    pyplot.xticks( xticks, xtick_labels )


    # set the limits of the matrix plot
    pyplot.xlim( (0,MIPmatrix.shape[1]) )
    pyplot.ylim( (0,MIPmatrix.shape[0]) )
#    print "\tY-limits=%s" % str(pyplot.ylim())
    # if we are not showing invalid_states, we must put a label for every entry of MIPdict.
    # If we are showing the invalid_states, just leave as the default


    # we do -1 because the graph starts at 0.
    highest_state = len(ylabel_and_MIPs) - 1
    
    # add the positions of the yticks
    # every state must have a position.  So do x+.5
    yticks_positions = [ x+.5 for x in range(len(ylabel_and_MIPs)) ] # positions [ 0.5 ... n+0.5 ]            

    # Now make our ytick labels
    ytick_labels = []
    for ylabel, row in ylabel_and_MIPs:
        try:
            ytick_label = dec2label(int(ylabel),label,pad_to=2**numunits)
        except:
            ytick_label = ylabel
        ytick_labels.append( ytick_label )
    
    # always center the ytick_positions to be at
    ytick_labels.reverse()
    print "ytick_positions=%s" % yticks_positions    
    print "ytick_labels=%s" % str(ytick_labels)

    pyplot.yticks( yticks_positions, ytick_labels )
    

    if filename:
        print "saving to %s" % filename
        pyplot.savefig(filename)
    else:
        pyplot.show()


def _phis_2_philist_with_labels( phis, show_invalidstates, label ):
    """This function returns an array of phiarr[position]=(label, value) from a list of phis"""
    
    philist = []
    
    for index, value in enumerate(phis):

        # if we're doing binary labels, convert the index to it's binary form.
        indexlabel = dec2label(index, label, pad_to=len(phis) )

        # if value is valid, add it.
        if value is not None:
            philist.append( (indexlabel, value) )

        # if value is invalid and show_invalidstates is on, append it.
        elif value is None and show_invalidstates:
            philist.append( (indexlabel,value) )
    
    return philist

def states_and_hist( phis, filename=None, label='decimal', show_invalidstates=True, title=None, plotting_eis=False ):
    """Make a barchart of PHI for each state as well as a PHI histogram.
    
    If plotting_eis is True, will adjust the labels, etc. to say ei(x_1) instead of $phi$
    
    """
    
#    print "Plotting states and histogram for phis=%s" % phis
    
    phis = _phis_2_philist_with_labels( phis, show_invalidstates, label )
    pprint( phis )

    numstates_to_show = len(phis)
    
    fig = pyplot.figure()
    
    # define the top subplot
    ax = fig.add_subplot(2,1,1)

    # plot the phis as a bar-chart
    # create arrays for the position and value of the reachable_values
    reachable_values,xposition_of_reachable_values = [], []
    for position, (label, val) in enumerate(phis):        
        if val is not None:
            xposition_of_reachable_values.append( position )
            reachable_values.append(val)

    xtick_positions = range(numstates_to_show)

#    pprint( xposition_of_reachable_values )
#    pprint( reachable_values )

    ax.bar( xposition_of_reachable_values, reachable_values, align='center' )


    # set the xtick labels
    xtick_labels = [ label for label, value in phis ]    
#    print "xtick_positions=%s\tlabels=%s" % ( xtick_positions, xtick_labels )
    pyplot.xticks( xtick_positions, xtick_labels )



    # calculate the height of the "special values (the red Xs and blue 0s)"
    # height of a special is 5% higher than the y_minimum
    height_of_specials = abs(pyplot.ylim()[0] + 0.05 * abs((pyplot.ylim()[1] - pyplot.ylim()[0])))

    # plot the zero-phi states.
    states_with_zerophi = [ index for index, (label, val) in enumerate(phis) if val == 0.0 ]
#    print "States with zerophi=%s height=%s" % (states_with_zerophi, height_of_specials)


    for state_with_zerophi in states_with_zerophi:
#        print "\tplotting at x=%s, height=%s" % ( state_with_zerophi, height_of_specials)
        ax.text( state_with_zerophi, height_of_specials, '0', family='monospace', size='large', color='b', horizontalalignment='center' )
#    raw_input('showed zerophi.')

    # print the invalidstates, if we're doing that.
    if show_invalidstates:
        # plot the unreachable phi states
        unreachable_states = [ i for i, (label, val) in enumerate(phis) if val is None]
        ax.plot( unreachable_states, [height_of_specials] * len(unreachable_states), marker='x', color='r', markersize=10.0, linewidth=0.0 )
        logging.debug("Unreachable states=%s" % unreachable_states)

    if title:
        ax.set_title( title )
    elif plotting_eis:
#        ax.set_title('$ei(X_1)$ \t $\mu$=%s  $\sigma^2$=%g' % ( mean(reachable_values), round(variation(reachable_values),3) ) )        
        ax.set_title('$ei(X_1)$')
    else:
#        ax.set_title('$\phi(X_1)$ \t $\mu$=%s  $\sigma^2$=%g' % ( mean(reachable_values), round(variation(reachable_values),3) ) )
        ax.set_title('$\phi(X_1)$')
   
   # if we're showing binary labels, the labels should be going vertically, not horizontal (for space efficiency reasons)

#        pylab.set(pylab.gca(), 'xticks', [1,2,3,4])
#        pylab.set( pylab.gca() )
#        mylabels = fig.set( fig.gca('xticklabels') )
#        pylab.set( mylabels, 'rotation', 'vertical')


    # add the xtick for the final state if it's not already there.
#    xtick_locs = pyplot.xticks()[0]
#    if (numstates-1) not in xtick_locs:
#        xtick_locs = npy.append( xtick_locs, numstates-1 )
#        logging.debug("Changed xtick_locations.  new xtick_locations=%s" % ( pyplot.xticks()[0] ) )


    # apply the xtick_locs and xtick_labels
#    pyplot.xticks( xtick_locs )

    # set the bounds of the plot
    pyplot.xlim( -0.5, (numstates_to_show-1) + 0.5 )
    pyplot.ylim(ymin=0)
    
    # HACK: if the highest Y is very low, make it a little taller.  This makes the display for all zeros look alright.
    if pyplot.ylim()[-1] < 10*height_of_specials:
        pyplot.ylim(ymax=10*height_of_specials)
        

    # set the labels
    ax.set_xlabel('$x_1$')
    if plotting_eis:
        ax.set_ylabel('$ei(x_1)$')
    else:
        ax.set_ylabel('$\phi(x_1)$')
            
    #######################################################################################
    # PART II -- plot the phi histogram
    #######################################################################################
    
    logging.debug("Ploting histogram of phis=%s" % ( reachable_values ) )
    ax2 = fig.add_subplot(2,1,2)

    # guess the number of bins that we need.
    n, guessedbins, patches = ax2.hist( reachable_values, normed=False, align='center',visible=False, figure=None )
#    print "bins=%s" % guessedbins

    # we require a minimum binwidth of .05.  If the binwidth was going to be greater than that
    # require that binwidth by setting minnumbins.
    minnumbins = (max(reachable_values) - min(reachable_values)) * 20
    numbins = max( minnumbins, len(guessedbins) )

    # plot the histogram with the proper number of bins.
    ax2.hist( reachable_values, bins=numbins, normed=False, align='center', facecolor='red' )

    # set the x/y ticks.
    # set the yticks to only have integers.
    ax2.set_yticks( [ x for x in ax2.get_yticks() if x % 1 == 0 ] )
    
    # set the x/y labels.
    # lowest possible x-value is zero.  So don't only need enough space to show the bar for 0.0
    if ax2.get_xlim()[0] < -0.05:
        ax2.set_xlim(xmin=-0.05)

    # set our labels.
    ax2.set_xlabel('$\phi$ Histogram')
    ax2.set_ylabel('# occurrences')
    
    if filename:
        pyplot.savefig(filename)
    else:
        pyplot.show()





def plotmatrix( matrix, title=None, filename=None, cmap=pyplot.cm.gray ):
    """Plot the display of a binary matrix"""

    # create a new figure
    fig = pyplot.figure()

    pyplot.delaxes()
    
#############################################################################################
# About colors:
# Default is that WHITE is connected.  But we want it so that BLACK is connected.
#############################################################################################
#    pyplot.pcolor( matrix[ ::-1,:], cmap=cmap, shading='faceted', alpha=0.8 )
#    pyplot.pcolor( abs(matrix[ ::-1,:]-1), cmap=cmap, shading='faceted', alpha=0.8 )    
    pyplot.pcolor( -matrix[ ::-1,:], cmap=cmap, shading='faceted', alpha=0.85 )    
#############################################################################################

#    pyplot.matshow( matrix, cmap=cmap, alpha=1.0, origin='upper' )

    # set the x and y ticks
    xtick_labels = [ int(x) for x in pyplot.xticks()[0] if x % 1 == 0  and x >= 1 ]  # labels [1...n]
    # we're zero based.  So pop off the last one.

    xticks = [ x-.5 for x in xtick_labels ] # positions [ 0.5 ... n+0.5 ]
    pyplot.xticks( xticks, xtick_labels )
    
    # Now for the yticks
    ytick_labels = [ int(x) for x in pyplot.yticks()[0] if x % 1 == 0 and x >= 1]  # labels [1...n]

    # we're zero based.  So pop off the last one.

    yticks = [ x-.5 for x in ytick_labels ] # positions [ 0.5 ... n+0.5 ]

    # flip the ytick_labels just like we flipped the matrix
    ytick_labels.reverse()
    pyplot.yticks( yticks, ytick_labels )



    # if we have a title, print it.
    if title:
        pyplot.title( title + '   black=connected')

    # plot the labels
    pyplot.xlabel('Destination Node')
    pyplot.ylabel('Origin Node')
    
    if filename:
        pyplot.savefig(filename)
    else:
        pyplot.show() 

def plot_cmatrix( m, filename, prog='neato', labelbase=10, caption='', **kwargs ):
    '''Prints a connection matrix as a directed graph.  Nonzero entries are connections.
        m = the matrix to be graphed
        filename = the output filename
        prog = what layout algorithm to plot with
        labelbase = the base of the index of each node.  Can be any integer value >=1
        caption = caption for the graph
    '''
    
    # check that the matrix is SQUARE
    if m.shape[0] != m.shape[1]:
        raise ValueError, "plot_as_graph -- passed matrix is not square."

    # Setting strict=False allows nodes to have self-connections
    # Also, set any graph properties in kwargs.
    A=pgv.AGraph(directed=True, strict=False, rankdir='BT', overlap='scale')
    A.graph_attr.update( kwargs )


    # Add each edge.
    fromnodes, tonodes = npy.nonzero(m)    
    for from_index, to_index in zip(fromnodes, tonodes):
        A.add_edge( from_index, to_index )
    
    # Now make the node labels
    if labelbase != 10:
        total_label_length = int(ceil(max(set( npy.concatenate( (fromnodes,tonodes) ) )) / int(labelbase)) +1)
    else:
        total_label_length = 0
    for node_index in set( npy.concatenate( (fromnodes,tonodes) ) ):
        node = A.get_node( node_index)
        label_length = len(npy.base_repr( node_index, base=labelbase ))
        node.attr['label'] = npy.base_repr( node_index+1, base=labelbase, padding=total_label_length - label_length )

    # Add the caption
    A.graph_attr['label'] = "%s (base %s)" % (caption, labelbase)
    
    # write it
    A.draw( filename, prog=prog )    
        
def plot_tmatrix( m, filename, selfcxn_shape='octagon', prog='neato', labelbase=2, caption='', noinputs_shape='triangle', **kwargs ):
    '''prints a state-transition matrix as a directed graph.  Nonzero entries are transitions.
    
        m = the matrix to be graphed
        filename = the output filename
        changeshapes = if True, use the noinputs_shape and self-cxn-shape
        prog = what layout algorithm to plot with
        labelbase = the base of the index of each node.  Can be any integer value >=1
        caption = caption for the graph
        noinputs_shape = Changes the shape of nodes with no incoming cxns to noinputs_shapeas
    '''
    
    # check that the matrix is SQUARE
    if m.shape[0] != m.shape[1]:
        raise ValueError, "plot_as_graph -- passed matrix is not square."

    # These must all be STRINGS or True/False values
    A=pgv.AGraph( directed=True, strict=False, overlap='scale', dpi='175' )
    A.graph_attr.update( kwargs )
    
    fromnodes, tonodes = npy.nonzero(m)

    # Add each edge
    for from_index, to_index in zip(fromnodes, tonodes):
        A.add_edge( from_index, to_index )
    
        # check to make sure the self-connections are present.
        if from_index == to_index:
            assert A.has_edge( from_index, to_index ), "Plot did not make the edge for a self-connection between %s and %s!" % (from_index, to_index)
            # Make self-cxns octagons
            if selfcxn_shape:
                (A.get_node( to_index )).attr['shape'] = selfcxn_shape

    
    # Now make the node labels
    if labelbase != 10:
        total_label_length = int(ceil(max(set( npy.concatenate( (fromnodes,tonodes) ) )) / int(labelbase)) +1)
    else:
        total_label_length = 0
    for node_index in set( npy.concatenate( (fromnodes,tonodes) ) ):
        node = A.get_node( node_index)
        label_length = len(npy.base_repr( node_index, base=labelbase ))
        node.attr['label'] = npy.base_repr( node_index, base=labelbase, padding=total_label_length - label_length )

#        print "in_degree=%s" % A.in_degree( node )
        # Change the shape of no-input nodes, if we're doing that.
        if noinputs_shape and not A.in_degree(node) and node.attr['shape'] != selfcxn_shape:
            node.attr['shape'] = noinputs_shape

    # Add the caption
    A.graph_attr['label'] = "%s (base %s)" % (caption, labelbase)

    # write it
    A.write( filename + '.neato' )
    
    print "Saving to '%s'..." % filename
    A.draw( filename, prog=prog )
