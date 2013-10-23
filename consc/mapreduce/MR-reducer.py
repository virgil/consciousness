#!/usr/bin/python
###############################################################
# This takes the minimum ei of each x1 state.
# Format of the STDIN input...:
'''
x1=0	MIPs=[[0, 1, 2, 6], [3, 4, 5]];[[0, 1, 2], [3, 4, 5, 6]]|min_ei=0.333333|normed_min_ei=0.111111|prob(x1)=0.421875|ei(x1)=1.24511
x1=4	MIPs=[[0, 1, 2, 6], [3, 4, 5]]|min_ei=0.333333|normed_min_ei=0.111111|prob(x1)=0.140625|ei(x1)=2.83007
x1=32	MIPs=[[0, 1, 2], [3, 4, 5, 6]]|min_ei=0.333333|normed_min_ei=0.111111|prob(x1)=0.140625|ei(x1)=2.83007
x1=36	MIPs=[[0, 1, 2, 6], [3, 4, 5]];[[0, 1, 2], [3, 4, 5, 6]]|min_ei=0.333333|normed_min_ei=0.111111|prob(x1)=0.046875|ei(x1)=4.41504
x1=64	MIPs=[[0, 1, 2, 5, 6], [3, 4]];[[0, 1, 2, 5, 6], [3], [4]];[[0, 1], [2, 3, 4, 5, 6]];[[0, 1], [2, 5, 6], [3, 4]];[[0, 3], [1, 4], [2, 5, 6]];[[0, 4], [1, 3], [2, 5, 6]];[[0], [1], [2, 3, 4, 5, 6]];[[0], [1], [2, 5, 6], [3], [4]]|min_ei=0.415037|normed_min_ei=0.207519|prob(x1)=0.140625|ei(x1)=2.83007
x1=68	MIPs=[[0, 1, 2, 5, 6], [3, 4]];[[0, 1, 2, 5, 6], [3], [4]]|min_ei=0.415037|normed_min_ei=0.207519|prob(x1)=0.046875|ei(x1)=4.41504
x1=96	MIPs=[[0, 1], [2, 3, 4, 5, 6]];[[0], [1], [2, 3, 4, 5, 6]]|min_ei=0.415037|normed_min_ei=0.207519|prob(x1)=0.046875|ei(x1)=4.41504
x1=100	MIPs=[[0, 1, 2, 6], [3, 4, 5]];[[0, 1, 2], [3, 4, 5, 6]]|min_ei=1|normed_min_ei=0.333333|prob(x1)=0.015625|ei(x1)=6
'''
###############################################################


import sys
from pprint import pprint

LOWEST = {}

while 1:
    line = (sys.stdin.readline()).rstrip()
    
    if not line:
        break

    # ignore all comments
    if line.startswith('#'):
        continue
    
    x1, rest = line.split('\t')    
    x1 = x1.split('=')[-1]
    restparts = rest.split('|')
    
    MIPs_str = ([x.split('=')[-1] for x in restparts if x.startswith('MIPs=') ])[0]
    ei = float( ([x.split('=')[-1] for x in restparts if x.startswith('ei(x1)=') ])[0] )
    min_ei = float( ([x.split('=')[-1] for x in restparts if x.startswith('min_ei=') ])[0] )
    normed_min_ei = float( ([x.split('=')[-1] for x in restparts if x.startswith('normed_min_ei=') ])[0] )
    prob_x1 = float( ([x.split('=')[-1] for x in restparts if x.startswith('prob(x1)=') ])[0] )
    
    MIPs = tuple(set(MIPs_str.split(';')))
    
    
    # if we found a new lowest normalized ei, replace.
    if not x1 in LOWEST or normed_min_ei < LOWEST[x1]['normed_min_ei']:
        LOWEST[x1] = { 'MIPs': MIPs, 'ei': ei, 'min_ei': min_ei, 'normed_min_ei': normed_min_ei, 'prob_x1': prob_x1 }

    # if two workers found partitions with the SAME normalied_min_ei, take the union of their partitions
    elif normed_min_ei == LOWEST[x1]['normed_min_ei']:
        LOWEST[x1]['MIPs'] = tuple( set(list(LOWEST[x1]['MIPs'])+list(MIPs)) )
        
    # Now do some sanity checks...
    assert prob_x1 == LOWEST[x1]['prob_x1'], "p(x1) from each worker should be the same!"
    assert ei == LOWEST[x1]['ei'], "ei(x1) from each worker should be the same!"    
    

for x1, parts in LOWEST.iteritems():
    print "x1=%s \t-> " % x1 ,
    
    print "phi=%s\t" % parts['min_ei'],
    print "MIPs=%s\t" % str(parts['MIPs']),
    print "prob(x1)=%s\t" % parts['prob_x1'],
    print "ei(x1)=%s\t" % parts['ei'], 
    print "phi_norm=%s" % parts['normed_min_ei']


# now print the average PHI
average_phi, sum_prob = 0.0, 0.0
for x1, parts in LOWEST.iteritems():
    sum_prob += parts['prob_x1']
    average_phi += parts['prob_x1'] * parts['min_ei']
    

assert 0.99999 <= sum_prob <= 1.00001, "sum_probability should be 1.0 but it is %s !" % sum_prob

print "Average PHI=%s" % round(average_phi,4)