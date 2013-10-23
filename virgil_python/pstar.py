#!/usr/bin/python

class Pstar():
	'''this class computes the probability of different events using modified probability distributions'''
	
	def __init__(self, TABLE, mm=None):
		self.TABLE = TABLE
		self.mm = None
	
	def base_p(self, x1, x2, y):
		'''computes the base probability p(x1,x1,y) using the markov model mm.
		We define the basep(x1,x2,y) = p(x1)*p(x2)*p(y|x1)*p(y|x2)		
		OR We define the basep(x1,x2,y) = p(y)*p(x1|y)*p(x2|y)
		'''
		
		z = 1.0
		# using the definition basep(x1,x2,y) = p(x1)*p(x2)*p(y|x1)*p(y|x2)		
		z *= Prob_xs( x1 )
		z *= Prob_xs( x2 )
		z *= Prob_y_given_xs(y,x1)
		z *= Prob_y_given_xs(y,x2)


		# using the definition basep(x1,x2,y) = p(y)*p(x1|y)*p(x2|y)
#		z *= Prob_y(y)
#		z *= Prob_xs_given_y(y,x1)
#		z *= Prob_xs_given_y(y,x2)


	def Prob_y(self, y_value ):
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

	def Prob_y_given_xs( self, xs, y_value ):
		'''returns p(X=x|Y=y) = p(Xs,Y=y) / p(Y=y)'''
		 top = Prob_xs_y( xs, y_value )

	    if top == 0.0:
	        return 0.0

	    bottom = Prob_y( y_value )
	    z = top / bottom

	    assert 0.0 <= z <= 1.0, "Not a valid probability! Eep! y_value=%s \t Xs=%s \t Xvalues=%s \n top=%s \t bottom=%s \t z=%s" % (y_value, xs.indices, xs.values, top, bottom, z)

	    return z
	    

	def Prob_y_given_xs( self, y_value, xs ):
	    '''returns p(Y=y|Xs=xs) = p(Xs=xs,Y=y) / p(Xs=xs)'''

	    top = Prob_xs_y( xs, y_value )

	    if top == 0.0:
	        return 0.0

	    bottom = Prob_xs( xs )

	    #    assert top <= bottom

	    z = top / bottom

	    assert 0.0 <= z <= 1.0, "Not a valid probability! Eep! y_value=%s \t Xs=%s \t Xvalues=%s \n top=%s \t bottom=%s \t z=%s" % (y_value, xs.indices, xs.values, top, bottom, z)

	    return z



	def Prob_xs( self, xs ):
	    "Compute the joint probability of all of these Xs occurring"
	    # sanity check our inputs

	    num_instances = 0

	    #    print "Xs=%s  Xvalues=%s" % ( Xs, Xvalues)

	    #    print "TABLE=%s" % TABLE.keys()
	    # for each input string in TABLE...
	    for keystring in TABLE.keys():
	        values = keystring.split(' ')

	        #        print "to_match values=%s" % str(values)

	        match=True
	        # If that string matches, num_instances += 1
	        for X, Xval in zip(xs.indices,xs.values):
	            if values[X] != Xval:
	                match = False

	        if match:
	            num_instances += 1
	        #            print "Matched!  num_instances=%s" % num_instances
	        #            raw_input('....')

	    z = float(num_instances) / float(len(TABLE))

	    assert 0.0 <= z <= 1.0, "Pr( Xs = xs ) wasn't between [0,1]."
	    return z
	
	def Prob_xs_y( self, xs, y_value ):
	    "Compute the joint probability of all of these Xs occurring along with Y=y"

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
