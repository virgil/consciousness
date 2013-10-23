#!/usr/bin/python
###############################################################
# 
###############################################################


import sys
from pprint import pprint

LOWEST = {}

while 1:
    line = sys.stdin.readline()
    
    if not line:
        break
    
    line = line.rstrip()
    
    x1, rest = line.split('\t')    
    restparts = rest.split('|')
#    print "restparts=%s" % restparts
    
    x1 = x1.split('=')[-1]
    ei = float( ([ x.split('=')[-1] for x in restparts if x.startswith('min_ei=') ])[0] )
    ei_over_norm = float(restparts[-1])
#    print "ei=%s" % ei


    if not x1 in LOWEST or ei_over_norm < LOWEST[x1][1]:
        LOWEST[x1] = (ei, ei_over_norm)
        

for x1, (ei, ei_over_norm) in LOWEST.iteritems():
    print "%(x1)s  ->  %(ei)s" % locals()
