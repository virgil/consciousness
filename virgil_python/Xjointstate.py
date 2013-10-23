#!/usr/bin/python

class Xjointstate(object):
    def __init__(self, Xindices, Xvalues ):
        '''define data-structures'''

        assert len(Xindices) == len(Xvalues)

        assert sorted(Xindices) == Xindices, "Xindices must be sorted"

        # TODO: SORT THE LIST OF Xindices.  Then sort Xstates the same way.
        self.indices = tuple(Xindices)
        self.values = tuple(Xvalues)

    def __eq__(self, other):
        '''test equality of two Xjointstates'''

#        print "running it.!"
#        raw_input('...')

        # if other is not an Xjointstate, it's false.
        if not isinstance(other, Xjointstate):
            return False

        if len(self) != len(other):
            return False

        for i_self, i_other in zip(self.indices, other.indices):
            if i_self != i_other:
                return False

        for v_self, v_other in zip(self.states, other.values):
            if v_self != v_other:
                return False

        return True
    
    def str(self):
        '''returns the string version of this Xstate'''
        z = 'indices: %s / values: %s' % (self.indices, self.values)
        return z

    def __len__(self):
        '''returns the length of this Xjointstate'''

#        print "running __len__()!"
#        raw_input('...')

        return len(self.indices)
