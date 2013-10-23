#!/usr/bin/python
# phi.py
# functions for computing integrated information
# David Balduzzi 5/2008

# The heart of the library is JEI, which computes EI. It is vectorized, so probably difficult to read.

# OPTIMIZATION.
# The obvious function to optimize is PHI. As it stands, it calls JEI every time, and recomputes the
# same basic probabilities over and over again. This could be short-circuited by pulling some of the
# work out of JEI into PHI. This would require keeping of track of ORDERING information.

from numpy import *
from probstat import *
# PROBSTAT produces fast Combinations from a list.
import random

global WHICH # choose rule applied by all elts
global PERM # allow complicated boolean rules, specific to each elt

# set some nice default values incase they aren't already defined.
try:
    PERM
except:
    PERM = 'sksk'

try:
    WHICH
except:
    WHICH = lambda suminput: array(suminput >= 5.0, dtype=int)
#    WHICH = 'and2'
    print "Setting WHICH=%s" % WHICH    


#from pprint import pprint



##########################################################################################
def BinaryNumberGenerator(n):
    # create binary numbers, 0...2^n-1
    pg=ones((2**n,n),int)
    for i in range(0,n):
        pg[range(0,2**i)*(2**(n-i-1)) + repeat(range(0,2**(n-i-1)),2**i)*(2**(i+1)),n-i-1]=range(0,1)*(2**(n-1))
    return pg
##########################################################################################
# random stuff
def ranseq(n,p=.5):
    # create random sequence of {0,1}'s, with p(1)=p
    seq=ones(n,int)
    for i in range(0,n):
        if p<random.random():
            seq[i]=0
    return seq

def smooth_hopf(adj,r=.8):
    # creates random connection matrix, with sparse entries, multiples of 1/8
    nadj=adj
    for i in range(0,8):
        for k in range(0,8):
            p=random.random()
            if p>r:
                nadj[i,k]=float(round((random.random()*2-1)*8))/8
    for i in range(0,8):
        nadj[i,i]=0
    overs=len(nadj[nadj!=0])-36
    if overs>0:
        for i in range(8):
            for k in range(8):
                p=random.random()
                if p<(float(2*overs)/64):
                    nadj[i,k]=0
    return nadj

def overlap(s1,s2):
    # overlap between Hopfield network states
    s1=s1*2-1; s2=s2*2-1;
    return float(dot(s1,s2))/len(s1)


# low level snippets
#
# SW and SP: set WHICH and PERM
# PINGS: output list of perturbations of n elts
# APPLY_GATE: apply rule chosen by WHICH to poss_states
# CINDICES: keeps track of indices when marginalizing
def sw(wval):
    # set value of WHICH inside module, for use in apply_gate
    # options: 'life', 'parity', 'and1', 'and2', 'and3', 'and4'
    global WHICH; WHICH=wval

def sp(wval):
    # set to 'permute' or whatever for use in d_p0s
    # if set to permute, ADJ variable fed to jei should look like output of perm_adj
    global PERM; PERM=wval

def perm_adj(permutn):
    # creates permutation matrix
    # meaning, each row implicitly represents input state for system, and explicitly represents output state
    n=int(log2(len(permutn)))
    adj=zeros((2**n,n),int)
    vecs=pings(n)
    for i in range(2**n):
        adj[i,:]=vecs[permutn[i],:]
    return adj

def entropy(vect):
    # normalizes (n x m) array of numbers and computes entropy of rows
    nvect=transpose(array(vect,float))/sum(transpose(vect),axis=0)  # normalize rows
    nvect=nvect+array(equal(nvect,0),float)     # convert 0s to 1s so doesn't blow up
    nvect=log10(nvect)*nvect
    return -sum(nvect,axis=0)/log10(2)

def relent(v1,v2,norml=False):
    # relative entropy
    if norml:
        v1=array(v1,float)/sum(v1); v2=array(v2,float)/sum(v2)
    return -dot(v1[v1>0],log10(v2[v1>0]))/log10(2) - entropy(v1)

"""Should be rewritten with a yield statement to give each binary number"""
def pings(n):
    # create binary numbers, 0...2^n-1
    pg=ones((2**n,n),int)
    for i in range(0,n):
        pg[range(0,2**i)*(2**(n-i-1)) + repeat(range(0,2**(n-i-1)),2**i)*(2**(i+1)),n-i-1]=range(0,1)*(2**(n-1))
    return pg

def complement(list,sub):
    # calculates L\S
    n=len(list); s=len(sub)
    comp=zeros(n-s,int); cnt=0
    for i in range(0,n):
        yes=True
        for k in range(0,s):
            if list[i]==sub[k]:
                yes=False
        if yes:
            comp[cnt]=list[i]
            cnt+=1
    return comp

def apply_gate(old):
    # applies element rules to pre-processed values
    global WHICH
#    print "WHICH=%s" % WHICH
    
    if WHICH=='life':
        temp=array(old==3,int)+array(old==13,int)+array(old==14,int)
        return array(temp>0,int)
    elif WHICH=='parity':
        return array(mod(old,2)==1,int)
    elif WHICH=='and1':
        return array(old>0,int)
    elif WHICH=='and2':
        return array(old>1,int)
    elif WHICH=='and3':
        return array(old>2,int)
    elif WHICH=='and4':
        return array(old>3,int)
    # if WHICH is a FUNCTION, run it on the input.
    elif type(WHICH) is type(lambda x: x):
        return array(WHICH(old),int)

def cindices(gaps,ents,reps,shift):
    # generate indices summed over when averaging out elements
    basvec=arange(ents)*shift
    stepfn=arange(reps)*gaps
    return tile(basvec,(reps,1))+transpose(tile(stepfn,(ents,1)))

def cindices2(ents,reps,shift):
    # generate indices summed over when averaging out elements
    basvec=arange(ents)*shift
    stepfn=arange(reps)
    step2=repeat(arange(int(reps/shift))*shift,shift)
    return tile(basvec,(reps,1))+transpose(tile(stepfn+step2,(ents,1)))

def marginalize(pdist,nodes):
    # [probably, marginalize is slow, this could be fixed by permuting]
    # marginalize onto nodes
    n=int(log2(len(pdist))); other=sort(n-1-complement(range(n),nodes))
    temp=pdist; cnt=0
    for i in other:
        cind=cindices2(2,2**(n-cnt-1),2**(i-cnt))
        temp=sum(transpose(temp.take(cind)),axis=0)
        cnt+=1
    return temp

def create_ff(drule):
    # take rules for D layer, and create connection matrix for feedforward network of S-->D layers
    M,n=drule.shape; m=int(log2(M))
    adj=zeros((2**(m+n),m+n),int)
    adj[:,m:m+n]=reshape(tile(drule,(1,2**n)),(2**(m+n),n))
    return adj

def makeadj(xi):
    # creates Hopfield connection matrix; accepts xi's of {0,1}s
    xi=xi*2-1
    sh=xi.shape
    if len(sh)>1:
        n=sh[1]; p=sh[0]
    else:
        p=1; n=sh[0]
    ADJ=zeros((n,n))
    for i in range(0,n):
        for k in range(0,n):
            for l in range(0,p):
                ADJ[i,k]=ADJ[i,k]+xi[l,i]*float(xi[l,k])/n;
    return ADJ


# game of life converters
#
# convert grd to ADJ and curr_state
def grd2line(grd):
    # converts grid to a line
    m,n=grd.shape
    return squeeze(reshape(grd,(1,m*n)))

def line2grd(ll,m,n):
    # converts line to mxn grid
    return reshape(ll,(n,m))

def grd2lpart(agr,m,n):
    # converts subset of a grid to a line subset
    a=agr.shape[0]; ss=zeros(a,int)
    for i in range(0,a):
        ss[i]=agr[i,0]*n+agr[i,1]
    return ss

def create_life(n,m,de=True):
    # creates connectivity matrix for n*m game of life grid
    grd=zeros((n*m,n*m),int)
    counter=0
    if de:
        # dead edges
        for i in range(0,n):
            for k in range(0,m):
                for l1 in range(-1,2):
                    for l2 in range (-1,2):
                        pos=(i+l1)*m+k+l2
                        if ((k+l2)>=0) and ((i+l1)>=0) and ((k+l2)<m) and ((i+l1)<n):
                            if (pos>=0) and (pos<(n*m)):
                                grd[pos,i*m+k]=1
    else:
        # torus edges
        for i in range(0,n):
            for k in range(0,m):
                for l1 in range(-1,2):
                    for l2 in range (-1,2):
                        pos=(i+l1)*m+k+l2
                        if ((k+l2)>=0) and ((i+l1)>=0) and ((k+l2)<m) and ((i+l1)<n):
                            if (pos>=0) and (pos<(n*m)):
                                grd[pos,i*m+k]=1
                        else:
                            if (k+l2)<0:
                                pos=pos+m
                            if (i+l1)<0:
                                pos=pos+n*m
                            if (k+l2)>=m:
                                pos=pos-m
                            if (i+l1)>=n:
                                pos=pos-n*m
                            if (pos>=0) and (pos<(n*m)):
                                grd[pos,i*m+k]=1
    return grd + array(eye(n*m),int)*10


# ei and phi
#
# mainc: computes main complex, calls phi
# phi: computes phi, calls jei (and sometimes tcut)
# tcut: EI over self-cut
# jei: EI
# d_p0s: computes probabilities
def d_p0s(ADJ,cell,curr_state):
    # computes probabilities for deterministic elements
    n=len(curr_state)
    poss_states=zeros((2**n,n),int)

#   print "cell (A/B)="
#   pprint(cell)

    poss_states[:,cell]=pings(n)
    if PERM=='permute':
        conv=2**(n-1-arange(n))
        new_states=ADJ[dot(poss_states,conv)]
    else:
        new_states=apply_gate(poss_states*matrix(ADJ))


        
#   print "poss_states="
#   pprint(poss_states)
#   print "new_states="
#   pprint(new_states)
    
    p0=array(new_states[:,cell]==curr_state,float)
#   print "p0="
#   pprint(p0)

    p0 = p0/sum(p0,axis=0)
#   print "p0="
#   pprint(p0)
        
    return p0

def jei(ADJ,cells,cx,curr_state,T=-1):
    # relent computation of effective information
    # ADJ: nxn connection matrix, unless PERM='permute'
    # cells: list of arrays, describing the partition
    # cx: list of elements in complex under consideration (should be concate of cells), can be length<njei
    # curr_state: state of whole system, not just cx
    # T: if \neq -1, then treat system as Hopfield net with T, else, deterministic
    parts=len(cells); c=len(cx); n=len(curr_state);
    ordering=concatenate(cells); a=zeros(parts,int);
    if T==-1:
        # deterministic network
        if c<n:
            C=complement(arange(n),ordering)
            tot=concatenate((ordering,C))
            temp=d_p0s(ADJ,tot,curr_state[tot])
        else:
            temp=d_p0s(ADJ,ordering,curr_state[ordering])
    else:
        # Hopfield network
        curr_state=curr_state*2-1
        poss_state=zeros((2**n,n),int)
        if c<n:
            C=complement(arange(n),ordering)
            poss_state[:,concatenate((ordering,C))]=pings(n)*2-1
        else:
            poss_state[:,ordering]=pings(n)*2-1
        temp=array(.5*(1+tanh(poss_state*matrix(ADJ)/T)))
        temp=temp[:,ordering]*curr_state[ordering]+array(curr_state[ordering]==-1,float)

#    print "ordering=%s" % ordering
#    print "TEMP=%s" % temp
    # marginalize over perturbations outside cx
    cind=cindices(2**(n-c),2**(n-c),2**c,1)
    p0=array(sum(temp[cind,:],axis=1),float)
    joint_x=prod(p0[:,range(len(ordering))],axis=1)
    x=sum(joint_x)
    if x==0 or x!=x:
        return -1
    else:
        # a posteriori of whole cx
        joint_x=joint_x/x
    # compute marginals
    
    header=0; p_a0=[]; add=zeros((parts))
    for k in range(0,parts):
        a[k]=len(cells[k])
        cind=cindices(1,2**header,2**(c-header),2**(c-header))
        ta=sum(p0[cind,:],axis=1)
        
        cind=cindices(2**(c-a[k]-header),2**(c-a[k]-header),2**a[k],1)
        header += a[k]
        temp_a=sum(ta[cind,:],axis=1)
        temp_a=prod(temp_a[:,header-a[k]:header],axis=1)
        
        # a posteriori of part k
        p_a0.append(temp_a/sum(temp_a))
    
    # create product distribution
    joint_aa=p_a0[0]; tot=a[0]
    for k in range(1,parts):
        tot += a[k]
        joint_aa=reshape(outer(p_a0[k],joint_aa),2**tot,1)
    joint_aa=joint_aa/sum(joint_aa)
    
#    print "  b: H[X0|X1=x]=%s" % ( entropy(joint_x) )
#    print "  b: joint_x=%s" % str(joint_x)
#    print "  b: mu0s=%s" % -log2(joint_aa[joint_aa>0])
#    print "  b: Sum_x0_in_X0=%s" % ( dot(joint_x[joint_aa>0],-log2(joint_aa[joint_aa>0])) )
#    print "  b: joint_aa=%s" % str(joint_aa)
#    print "  b: joint_x[aa>0]=%s" % str(joint_x[joint_aa>0])
    return -entropy(joint_x) - dot(joint_x[joint_aa>0],log2(joint_aa[joint_aa>0]))
#    return -entropy(joint_x) - dot(joint_x[joint_aa>0],log10(joint_aa[joint_aa>0]))/log10(2)


def tcut(ADJ,cx,curr_state,T=-1):
    # as for jei, but computes information generated by self-cut
    c=len(cx); n=len(curr_state); ordering=cx;
    if T==-1:
        # deterministic network
        if c<n:
            C=complement(arange(n),ordering)
            tot=concatenate((ordering,C))
            temp=d_p0s(ADJ,tot,curr_state[tot])
        else:
            temp=d_p0s(ADJ,ordering,curr_state[ordering])
    else:
        # Hopfield network
        curr_state=curr_state*2-1
        poss_state=zeros((2**n,n),int)
        if c<n:
            C=complement(arange(n),ordering)
            poss_state[:,concatenate((ordering,C))]=pings(n)*2-1
        else:
            poss_state[:,ordering]=pings(n)*2-1
        temp=array(.5*(1+tanh(poss_state*matrix(ADJ)/T)))
        temp=temp[:,ordering]*curr_state[ordering]+array(curr_state[ordering]==-1,float)
        
    # marginalize over perturbations outside cx
    cind=cindices(2**(n-c),2**(n-c),2**c,1)
    p0=array(sum(temp[cind,:],axis=1),float)
    joint_x=prod(p0[:,range(len(ordering))],axis=1)
    x=sum(joint_x)
    if x==0:
        return -1
    else:
        joint_x=joint_x/x
    return c-entropy(joint_x)

# phi(ADJ,lj,curr_state,T)

def phi(ADJ,cx,curr_state,T=-1,iself=False):
    """
    This function returns (PHI,partA,partB) for state curr_state over units cx, using cxnmatrix ADJ.
    
Ex: samples.vFig2b
>>> phi( array([ [0,1,1], [1,0,1], [1,1,0] ]), range(3), array([0,0,0]), iself=False )
(0.6699250014423126, [0], [1, 2])

>>> phi( array([ [0,1,1], [1,0,1], [1,1,0] ]), range(3), array([0,0,0]), iself=True )
(1.0, [0, 1, 2], [])

Ex: samples.testand
>>> phi( array([[0, 0, 1], [0, 0, 1], [0, 0, 0]]), range(3), array([0,0,0]), iself=False )
(0.33333333333333393, [0], [1, 2])
>>> phi( array([[0, 0, 1], [0, 0, 1], [0, 0, 0]]), range(3), array([0,0,0]), iself=True )
(0.41503749927884392, [0, 1, 2], [])
>>> phi( array([[0, 0, 1], [0, 0, 1], [0, 0, 0]]), range(3), array([1,1,1]), iself=False )
(None, [], [])

Ex: samples.c1
>>> phi( array([[0, 1, 0, 0, 0, 0, 0, 1],[1, 0, 1, 0, 0, 0, 0, 0],[0, 1, 0, 1, 0, 0, 0, 0],[0, 0, 1, 0, 1, 0, 0, 0],[0, 0, 0, 1, 0, 1, 0, 0],[0, 0, 0, 0, 1, 0, 1, 0],[0, 0, 0, 0, 0, 1, 0, 1],[1, 0, 0, 0, 0, 0, 1, 0]]), range(8), array([1, 0, 1, 0, 0, 1, 0, 1]), iself=True) 
(0.0, [0, 7], [1, 2, 3, 4, 5, 6])

>>> phi( array([[0, 1, 0, 0, 0, 0, 0, 1],[1, 0, 1, 0, 0, 0, 0, 0],[0, 1, 0, 1, 0, 0, 0, 0],[0, 0, 1, 0, 1, 0, 0, 0],[0, 0, 0, 1, 0, 1, 0, 0],[0, 0, 0, 0, 1, 0, 1, 0],[0, 0, 0, 0, 0, 1, 0, 1],[1, 0, 0, 0, 0, 0, 1, 0]]), range(8), array([1, 0, 1, 0, 1, 0, 0, 1]), iself=True) 
(None, [], [])
"""
    # compute phi for bipartitions
    # iself: do we care about the self-cut
    # use probstat here to compute Combinations
    nphi=1000; curr_norm=1; A=0; B=0
    
    c=len(cx); top=int(c/2); cxa=array(cx,int)
    
    for i in range(1,top+1):
#        print "cx=%s" % cx
#        print "i=%s" % i
        items=Combination(cx,i)     
        if (2*i)==c:
            items=items[0:int(len(items)/2)]
        for lj in items:
            compl=complement(cx,lj)
            tphi=jei(ADJ,[array(lj,int),array(compl,int)],cxa,curr_state,T)
            small=min(i,c-i)
            if (nphi/curr_norm) > (tphi/small):
                curr_norm=small; nphi=tphi; A=lj; B=compl
    if iself:
        tphi=tcut(ADJ,cx,curr_state,T)
        if (nphi/curr_norm) > (tphi/c):
            nphi=tphi; curr_norm=c; A=cx; B=[]
            
    if nphi == -1:
        nphi, A, B = None, [], []

            
    return (nphi, list(A), list(B))

def mainc(ADJ,cx,curr_state,T=-1):
    # find complex contained in cx with highest phi
    nphi=0; n=len(cx); A=0; B=0
    for i in range(2,n+1):
        items=Combination(cx,i)
        for lj in items:
            tphi,tA,tB=phi(ADJ,lj,curr_state,T)
            if tphi>nphi:
                nphi=tphi
                A, B = tA, tB
    return nphi,A,B


# game of life wrappers
#
# uses converters to take GRD, turn it into ADJ + CURR_STATE; then output EI or PHI
def l_tcut(grd,de=True):
    # Wrapper for computing EI on self-cut for game of life.
    # input is a mxn grid, with whatever configuration
    # Note: check WHICH for apply_gate!
    # de: are edges dead or wrap-around (torus)?
    m,n=grd.shape;
    adj=create_life(m,n,de);
    curr_state=grd2line(grd);
    cx=range(m*n)
    
    return tcut(adj,cx,curr_state)

def l_jei(grd,gcells,de=True):
    # Wrapper for computing EI on game of life.
    # gcells: list of arrays. Array entries should be grid-coords FLATTENED onto LINE
    # Note: check apply_gate!
    m,n=grd.shape;
    adj=create_life(m,n,de);
    curr_state=grd2line(grd);
    cx=concatenate(gcells)
    return jei(adj,gcells,cx,curr_state)

def l_phi(grd,cx,de=True):
    # compute phi on subset of a grid
    m,n=grd.shape;
    adj=create_life(m,n,de);
    curr_state=grd2line(grd);
    return phi(adj,cx,curr_state)



def _test():
    """ Run all of our unit tests. """
    import doctest
    doctest.testmod()

if __name__ == "__main__":

    _test()


