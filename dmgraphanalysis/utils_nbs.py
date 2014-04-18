import numpy as np
import networkx as netwx

import time


import scipy.stats as stat

import itertools as it

import os


#from define_variables import conf_interval_binom,t_test_thresh

#from define_variables import conf_interval_binom_fdr,t_test_thresh_fdr

############################################################################### NBS ##############################################################
################################################################### original function by Zalesky #################################################
##################################################################################################################################################

def ttest2(X,Y):
    """ Compute the two-sided t-statistic of X,Y
"""
    t = np.mean(X) - np.mean(Y)
    n1 = len(X) * 1.
    n2 = len(Y) * 1.
    s = np.sqrt( ( (n1-1) * np.var(X,ddof=1) + (n2-1)*np.var(Y,ddof=1) ) / (n1+n2-2.) )
    t = t / (s*np.sqrt(1/n1+1/n2))
    return t

#def compute_nbs(X, Y, t_test_thresh, K = 1000, TAIL = 'both'):
    #""" Computes the network-based statistic (NBS) as described in [1].
#Performs the NBS for populations X and Y for a
#T-statistic threshold of t_test_thresh. The third dimension of X and Y
#references a particular member of the populations. The first two
#dimensions reference the connectivity value of a particular edge
#comprising the connectivity matrix. For example, X[i,j,k] stores the
#connectivity value corresponding to the edge between i and j for the
#kth memeber of the population. PVAL is a vector of corrected p-values
#for each component identified. If at least one of the p-values is
#less than 0.05, then the omnibus null hypothesis can be rejected at
#5% significance. The null hypothesis is that the value of
#connectivity at each edge comes from distributions of equal mean
#between the two populations.
#Parameters
#----------
#X, Y : ndarray
#THRES : float
#K : integer, default = 1000, optional.
#Enables specification of the number of
#permutations to be generated to estimate the empirical null
#distribution of maximal component size.

#TAIL : {'equal', 'left', 'right'}, optional
#Enables specification of the type
#of alternative hypothesis to test. If TAIL:
#'equal' - alternative hypothesis is means are not equal (default)
#'left' - mean of population X < mean of population Y
#'right' - mean of population X > mean of population Y
#Returns
#-------
#PVAL : ndarray
#p-values for each component
#ADJ : ndarray
#Returns an adjacency matrix identifying the edges comprising each component.
#Edges corresponding to the first p-value stored in the vector PVAL are assigned
#the value 1 in the adjacency matrix ADJ, edges corresponding to the second
#p-value are assigned the value 2, etc.
#NULL : ndarray
#Returns a vector of K samples
#from the the null distribution of maximal component size.

#ALGORITHM DESCRIPTION
#The NBS is a nonparametric statistical test used to isolate the
#components of an N x N undirected connectivity matrix that differ
#significantly between two distinct populations. Each element of the
#connectivity matrix stores a connectivity value and each member of
#the two populations possesses a distinct connectivity matrix. A
#component of a connectivity matrix is defined as a set of
#interconnected edges.

#The NBS is essentially a procedure to control the family-wise error
#rate, in the weak sense, when the null hypothesis is tested
#independently at each of the N(N-1)/2 edges comprising the
#connectivity matrix. The NBS can provide greater statistical power
#than conventional procedures for controlling the family-wise error
#rate, such as the false discovery rate, if the set of edges at which
#the null hypothesis is rejected constitues a large component or
#components.
#The NBS comprises fours steps:
#1. Perfrom a two-sample T-test at each edge indepedently to test the
#hypothesis that the value of connectivity between the two
#populations come from distributions with equal means.
#2. Threshold the T-statistic available at each edge to form a set of
#suprathreshold edges.
#3. Identify any components in the adjacency matrix defined by the set
#of suprathreshold edges. These are referred to as observed
#components. Compute the size of each observed component
#identified; that is, the number of edges it comprises.
#4. Repeat K times steps 1-3, each time randomly permuting members of
#the two populations and storing the size of the largest component
#identified for each permuation. This yields an empirical estimate
#of the null distribution of maximal component size. A corrected
#p-value for each observed component is then calculated using this
#null distribution.

#[1] Zalesky A, Fornito A, Bullmore ET (2010) Network-based statistic:
#Identifying differences in brain networks. NeuroImage.
#10.1016/j.neuroimage.2010.06.041

#Written by: Andrew Zalesky, azalesky@unimelb.edu.au
#Rewritten for Python: Stephan Gerhard, connectome@unidesign.ch

#"""

    ## check input matrices
    #Ix,Jx,nx = X.shape
    #Iy,Jy,ny = Y.shape
    
    #assert Ix == Iy
    #assert Jx == Jy
    #assert Ix == Jx
    #assert Iy == Jy
    
    ## number of nodes
    #N = Ix

    ## Only consider elements above upper diagonal due to symmetry
    #ind_mask = ( np.triu(np.ones( (N,N) ),1) == 1 )
    #ind_i, ind_j = np.nonzero( np.triu(np.ones( (N,N) ),1) )
    
    ## Number of edges
    #M = N * (N - 1) / 2
    
    ## Look up table
    #ind2ij = np.zeros( (M,2) , dtype = np.int16)
    #ind2ij[:,0] = ind_i
    #ind2ij[:,1] = ind_j
    
    ## Vectorize connectivity matrices
    ## Not necessary, but may speed up indexing
    ## Uses more memory since cmat temporarily replicates X
    #cmat = np.zeros( (M, nx) )
    #pmat = np.zeros( (M, ny) )
    #for i in range(nx):
        #cmat[:,i] = X[ind2ij[:,0], ind2ij[:,1],i].ravel()
    #for i in range(ny):
        #pmat[:,i] = Y[ind2ij[:,0], ind2ij[:,1],i].ravel()
    
    ## Perform T-test at each edge
    #t_stat = np.zeros( M )
    #for i in range(M):
        ## compute ttest2, assume independent random samples
        #t = ttest2(cmat[i,:], pmat[i,:])
    
        #t_stat[i] = t

    #if TAIL == 'both':
        #t_stat = np.abs( t_stat )
    #elif TAIL == 'left':
        #t_stat = -t_stat
    #elif TAIL == 'right':
        #pass
    #else:
        #raise('Tail option not recognized')
 
    ## Threshold
    #ind_t = np.where( t_stat > t_test_thresh )
    
    ## Suprathreshold adjacency matrix
    #ADJ = np.zeros( (N,N) )
    #reledg = ind2ij[ind_t[0]] # relevant edges
    #ADJ[ reledg[:,0], reledg[:,1] ] = 1 # this yields a binary matrix, selecting the edges that are above threshold
    #ADJ = ADJ + ADJ.T

    ## Find network components
    #G = netwx.from_numpy_matrix(ADJ)
    ## Return connected components as subgraphs.
    #comp_list = netwx.connected_component_subgraphs(G)

    ## store the number of edges for each subgraph component
    #nr_edges_per_component = np.zeros( len(comp_list) )
    #nr_edges_per_component_bigenough = []
    
    #for idx, componentG in enumerate(comp_list):
        #nr_edges_per_component[idx] = componentG.number_of_edges()
        
        ## if number of edges bigger than zero
        #if nr_edges_per_component[idx] > 0:
            #nr_edges_per_component_bigenough.append(nr_edges_per_component[idx])
            
            ## store the number of edges of the component as value in the adjacency matrix
            #for ed in componentG.edges():
                #ADJ[ed[0],ed[1]] = ADJ[ed[1],ed[0]] = idx + 1
                ## if we would like to store the number of edges per component
                ## ADJ[ed[0],ed[1]] = ADJ[ed[1],ed[0]] = nr_edges_per_component[idx]
    
    ## renaming
    #sz_links = nr_edges_per_component_bigenough
    
    #if len(sz_links) > 0:
        #max_sz = np.max(sz_links)
    #else:
        #max_sz = 0

    #if False:
        ## additionally, store all the components in the matrix with the value of their number of edges
        #all_components = np.zeros( (N,N) )
        #for idx, componentG in enumerate(comp_list):
            
            #tmp_max = netwx.to_numpy_matrix( componentG , nodelist = range(N) )
            ## set nonzero to number of edges
            #tmp_max[tmp_max!=0.0] = componentG.number_of_edges()
            #all_components[:,:] = all_components[:,:] + tmp_max
    
    #print "Max component size is: %s" % max_sz
        
    ## Empirically estimate null distribution of maximum component size by
    ## generating K independent permutations.
    #print "====================================================="
    #print "Estimating null distribution with permutation testing"
    #print "====================================================="
    
    #hit=0.0
    #NULL = np.zeros( (K, 1) )
    ## stack matrices for permutation
    #d_stacked = np.hstack( (cmat, pmat) )

    #for k in range(K):
        
        #t1 = time.time()
        
        ## Randomize
        #indperm = np.random.permutation( nx+ny )
        #d = d_stacked[:, indperm].copy()

        ##################
        
        ## Perform T-test at each edge
        #t_stat_perm = np.zeros( M )
        #for i in range(M):
            ## assume independent random samples
            #t = ttest2(d[i,:nx], d[i,nx:nx+ny])
        
            #t_stat_perm[i] = t
        
        #if TAIL == 'both':
            #t_stat_perm = np.abs( t_stat_perm )
        #elif TAIL == 'left':
            #t_stat_perm = -t_stat_perm
        #elif TAIL == 'right':
            #pass
        #else:
            #raise('Tail option not recognized')
     
        ## Threshold
        #ind_t = np.where( t_stat_perm > t_test_thresh )
        
        ## Suprathreshold adjacency matrix
        #adj_perm = np.zeros( (N,N) )
        #reledg = ind2ij[ind_t[0]] # relevant edges
        #adj_perm[ reledg[:,0], reledg[:,1] ] = 1 # this yields a binary matrix, selecting the edges that are above threshold
        #adj_perm = adj_perm + adj_perm.T
        
        ## Find network components
        #G = netwx.from_numpy_matrix(adj_perm)
        ## Return connected components as subgraphs.
        #comp_list = netwx.connected_component_subgraphs(G)
        
        ## store the number of edges for each subgraph component
        #nr_edges_per_component = np.zeros( len(comp_list) )
        #for idx, componentG in enumerate(comp_list):
            #nr_edges_per_component[idx] = componentG.number_of_edges()
        
        ## more then one node (= at least one edge)
        #nr_edges_per_component_bigenough = nr_edges_per_component[nr_edges_per_component>0]
        
        ## renaming
        #sz_links_perm = nr_edges_per_component_bigenough
        
        #if len(sz_links_perm) > 0:
            #sz_links_perm_max = np.max(sz_links_perm)
        #else:
            #sz_links_perm_max = 0
    
        #NULL[k] = sz_links_perm_max

        ## if the component size of this random permutation is bigger than
        ## the component size of the group difference computed above, this is a hit
        #if NULL[k] >= max_sz:
            #hit = hit + 1
            
        
        #t2 = time.time()
        
        #print "Perm %d of %d (took %d s). Perm max is: %d. Observed max is: %d. P-val estimate is: %0.3f" % ((k+1),K,int(t2-t1),NULL[k],max_sz,hit/(k+1))

    ## Calculate p-values for each component
    #PVAL = np.zeros( len(sz_links) )
    #for i in range( len(sz_links) ):
        #PVAL[i] = len( NULL[NULL >= sz_links[i]] ) * 1.0 / K
        
    #return (PVAL, ADJ, NULL, sz_links)

    ################################################################ Modified from Zalesky (integrate some intermediate save) ###################################################################################
    ################################################### ttest

def find_biggest_component_size_ttest(X,Y,t_test_thresh):

    # number of nodes
    N = X.shape[0]
   
    # Perform binomial test at each edge
    ADJ = np.zeros((N,N),dtype = 'int')
    
    for i,j in it.combinations(range(N), 2):
        
        t_val = ttest2(X[i,j,:],Y[i,j,:])
        
        if t_val > t_test_thresh:
                ADJ[i,j] = ADJ[j,i] = 1
        
    

    # Find network components
    G = netwx.from_numpy_matrix(ADJ)
    # Return connected components as subgraphs.
    comp_list = netwx.connected_component_subgraphs(G)

    # store the number of edges for each subgraph component
    nr_edges_per_component = np.zeros( len(comp_list) )
    nr_edges_per_component_bigenough = []
    
    #ADJ_signif = np.zeros( (N,N) )
    
    for idx, componentG in enumerate(comp_list):
        nr_edges_per_component[idx] = componentG.number_of_edges()
        
        # if number of edges bigger than zero
        if nr_edges_per_component[idx] > 0:
            nr_edges_per_component_bigenough.append(nr_edges_per_component[idx])
            
    
            # store the number of edges of the component as value in the adjacency matrix
            for ed in componentG.edges():
                ADJ[ed[0],ed[1]] = ADJ[ed[1],ed[0]] = idx + 1
                # if we would like to store the number of edges per component
                # ADJ_signif[ed[0],ed[1]] = ADJ_signif[ed[1],ed[0]] = nr_edges_per_component[idx]
    
        #else:
            #print "Warning, no edges have for component " + str(idx)
    
    # renaming
    sz_links = nr_edges_per_component_bigenough
    
    if len(sz_links) > 0:
        max_sz = np.max(sz_links)
    else:
        max_sz = 0

    return max_sz,sz_links,ADJ
        
        
def compute_nbs_ttest(d_stacked, nx, ny , t_test_thresh, K = 1000):
    """ Computes the network-based statistic (NBS) as described in [1].
Performs the NBS for populations X and Y for a
T-statistic threshold of t_test_thresh. The third dimension of X and Y
references a particular member of the populations. The first two
dimensions reference the connectivity value of a particular edge
comprising the connectivity matrix. For example, X[i,j,k] stores the
connectivity value corresponding to the edge between i and j for the
kth memeber of the population. PVAL is a vector of corrected p-values
for each component identified. If at least one of the p-values is
less than 0.05, then the omnibus null hypothesis can be rejected at
5% significance. The null hypothesis is that the value of
connectivity at each edge comes from distributions of equal mean
between the two populations.
Parameters
----------
X, Y : ndarray
conf_interval_binom : float (between 0. and 1.0), default = 0.95, optional.
percentage of confidence for test of difference
K : integer, default = 1000, optional.
Enables specification of the number of
permutations to be generated to estimate the empirical null
distribution of maximal component size.

Returns
-------
PVAL : ndarray
p-values for each component
ADJ : ndarray
Returns an adjacency matrix identifying the edges comprising each component.
Edges corresponding to the first p-value stored in the vector PVAL are assigned
the value 1 in the adjacency matrix ADJ, edges corresponding to the second
p-value are assigned the value 2, etc.

Modified from the original so that only components surviving a given threshold are exported

NULL : ndarray
Returns a vector of K samples
from the the null distribution of maximal component size.

ALGORITHM DESCRIPTION
The NBS is a nonparametric statistical test used to isolate the
components of an N x N undirected connectivity matrix that differ
significantly between two distinct populations. Each element of the
connectivity matrix stores a connectivity value and each member of
the two populations possesses a distinct connectivity matrix. A
component of a connectivity matrix is defined as a set of
interconnected edges.

The NBS is essentially a procedure to control the family-wise error
rate, in the weak sense, when the null hypothesis is tested
independently at each of the N(N-1)/2 edges comprising the
connectivity matrix. The NBS can provide greater statistical power
than conventional procedures for controlling the family-wise error
rate, such as the false discovery rate, if the set of edges at which
the null hypothesis is rejected constitues a large component or
components.
The NBS comprises fours steps:
1. Perfrom a two-sample T-test at each edge indepedently to test the
hypothesis that the value of connectivity between the two
populations come from distributions with equal means.
2. Threshold the T-statistic available at each edge to form a set of
suprathreshold edges.
3. Identify any components in the adjacency matrix defined by the set
of suprathreshold edges. These are referred to as observed
components. Compute the size of each observed component
identified; that is, the number of edges it comprises.
4. Repeat K times steps 1-3, each time randomly permuting members of
the two populations and storing the size of the largest component
identified for each permuation. This yields an empirical estimate
of the null distribution of maximal component size. A corrected
p-value for each observed component is then calculated using this
null distribution.

[1] Zalesky A, Fornito A, Bullmore ET (2010) Network-based statistic:
Identifying differences in brain networks. NeuroImage.
10.1016/j.neuroimage.2010.06.041

Written by: Andrew Zalesky, azalesky@unimelb.edu.au
Rewritten for Python: Stephan Gerhard, connectome@unidesign.ch

"""

    print d_stacked.shape
    
    assert d_stacked.shape[2] == nx + ny
    
    t1 = time.time()
    
    max_sz,sz_links,nbs_adj_matrix = find_biggest_component_size_ttest(d_stacked[:,:,:nx],d_stacked[:,:,nx:nx+ny],t_test_thresh)

    
    
    print 'save nbs sz_links file'
    nbs_stats_file  = os.path.abspath('nbs_sz_links_'+ str(t_test_thresh) +'.txt')
    np.savetxt(nbs_stats_file,sz_links,fmt = '%f')
    
    
    
    print 'save nbs adj mat file' 
    nbs_adj_mat_file= os.path.abspath('nbs_adj_matrix_t_test_thresh_'+ str(t_test_thresh) +'.npy')
    np.save(nbs_adj_mat_file,nbs_adj_matrix)
    
    nbs_adj_text_file= os.path.abspath('nbs_adj_matrix_t_test_thresh_'+ str(t_test_thresh) +'.txt')
    np.savetxt(nbs_adj_text_file,nbs_adj_matrix,fmt = '%d')
    
    
    t2 = time.time()
    
    #if False:
        ## additionally, store all the components in the matrix with the value of their number of edges
        #all_components = np.zeros( (N,N) )
        #for idx, componentG in enumerate(comp_list):
            
            #tmp_max = netwx.to_numpy_matrix( componentG , nodelist = range(N) )
            ## set nonzero to number of edges
            #tmp_max[tmp_max!=0.0] = componentG.number_of_edges()
            #all_components[:,:] = all_components[:,:] + tmp_max
    
    print "Max component size is: %s, time = %d" % (max_sz,int(t2-t1))
        
    # Empirically estimate null distribution of maximum component size by
    # generating K independent permutations.
    print "====================================================="
    print "Estimating null distribution with permutation testing"
    print "====================================================="
    
    hit=0.0
    NULL = np.zeros( (K, 1) )
    # stack matrices for permutation
   
    
    for k in range(K):
        
        t1 = time.time()
        
        # Randomize
        indperm = np.random.permutation( nx+ny )
        
        
        d_perm = d_stacked[:,:, indperm]

        #################
        
        max_sz_perm,sz_links_perm,nbs_adj_matrix = find_biggest_component_size_ttest(d_perm[:,:,:nx],d_perm[:,:,nx:nx+ny],t_test_thresh)
        
        NULL[k] = max_sz_perm

        # if the component size of this random permutation is bigger than
        # the component size of the group difference computed above, this is a hit
        if NULL[k] >= max_sz:
            hit = hit + 1
            
        
        t2 = time.time()
        
        print "Perm %d of %d (took %d s). Perm max is: %d. Observed max is: %d. P-val estimate is: %0.3f" % ((k+1),K,int(t2-t1),NULL[k],max_sz,hit/(k+1))

    # Calculate p-values for each component
    PVAL = np.zeros( len(sz_links) )
    for i in range( len(sz_links) ):
        PVAL[i] = len( NULL[NULL >= sz_links[i]] ) * 1.0 / K
        
        
    
    ## Suprathreshold adjacency matrix
    #ADJ_signif = np.zeros( (N,N) )
    
    #for idx, componentG in enumerate(comp_list):
        #nr_edges_per_component[idx] = componentG.number_of_edges()
        
        ## if number of edges bigger than zero
        #if PVAL[idx] > 0.1:
            
            ## store the number of edges of the component as value in the adjacency matrix
            #for ed in componentG.edges():
                #ADJ_signif[ed[0],ed[1]] = ADJ_signif[ed[1],ed[0]] = idx + 1
                ## if we would like to store the number of edges per component
                ## ADJ_signif[ed[0],ed[1]] = ADJ_signif[ed[1],ed[0]] = nr_edges_per_component[idx]
    
        
    return (PVAL, NULL, nbs_adj_mat_file)
        
    #return (PVAL, ADJ, ADJ_signif, NULL, sz_links)

        
    ########################################################################### Binomial distribution (special coclass)
    
def binon_CI_test(X,Y,conf_interval_binom):
    """ Compute binomial comparaison
"""
    nX = len(X) * 1.
    nY = len(Y) * 1.
    
    pX = np.sum(X == 1)/nX

    pY = np.sum(Y == 1)/nY

    #print pX,pY,np.absolute(pX-pY) 

    SE = np.sqrt(pX * (1-pX)/nX + pY * (1-pY)/nY)

    #print SE
   
    norm =  stat.norm.ppf(1-conf_interval_binom/2) * SE
    
    #if (np.absolute(pX-pY) > norm) == True:
        #print pX,pY,np.absolute(pX-pY),norm
   
    
    return np.absolute(pX-pY) > norm

    ### adapted from Zalesky for binomial distrib ############
def find_biggest_component_size_binom(X,Y,conf_interval_binom):

    # number of nodes
    N = X.shape[0]
   
    # Perform binomial test at each edge
    ADJ = np.zeros((N,N),dtype = 'int')
    
    for i,j in it.combinations(range(N), 2):
        
        ADJ[i,j] = ADJ[j,i] = binon_CI_test(X[i,j,:],Y[i,j,:],conf_interval_binom)
        
    

    # Find network components
    G = netwx.from_numpy_matrix(ADJ)
    # Return connected components as subgraphs.
    comp_list = netwx.connected_component_subgraphs(G)

    # store the number of edges for each subgraph component
    nr_edges_per_component = np.zeros( len(comp_list) )
    nr_edges_per_component_bigenough = []
    
    #ADJ_signif = np.zeros( (N,N) )
    
    for idx, componentG in enumerate(comp_list):
        nr_edges_per_component[idx] = componentG.number_of_edges()
        
        # if number of edges bigger than zero
        if nr_edges_per_component[idx] > 0:
            nr_edges_per_component_bigenough.append(nr_edges_per_component[idx])
            
    
            # store the number of edges of the component as value in the adjacency matrix
            for ed in componentG.edges():
                ADJ[ed[0],ed[1]] = ADJ[ed[1],ed[0]] = idx + 1
                # if we would like to store the number of edges per component
                # ADJ_signif[ed[0],ed[1]] = ADJ_signif[ed[1],ed[0]] = nr_edges_per_component[idx]
    
        #else:
            #print "Warning, no edges have for component " + str(idx)
    
    # renaming
    sz_links = nr_edges_per_component_bigenough
    
    if len(sz_links) > 0:
        max_sz = np.max(sz_links)
    else:
        max_sz = 0

    return max_sz,sz_links
        
        
def compute_nbs_binom(d_stacked, nx, ny , conf_interval_binom, K = 1000):
    """ Computes the network-based statistic (NBS) as described in [1].
Performs the NBS for populations X and Y for a
T-statistic threshold of t_test_thresh. The third dimension of X and Y
references a particular member of the populations. The first two
dimensions reference the connectivity value of a particular edge
comprising the connectivity matrix. For example, X[i,j,k] stores the
connectivity value corresponding to the edge between i and j for the
kth memeber of the population. PVAL is a vector of corrected p-values
for each component identified. If at least one of the p-values is
less than 0.05, then the omnibus null hypothesis can be rejected at
5% significance. The null hypothesis is that the value of
connectivity at each edge comes from distributions of equal mean
between the two populations.
Parameters
----------
X, Y : ndarray
conf_interval_binom : float (between 0. and 1.0), default = 0.95, optional.
percentage of confidence for test of difference
K : integer, default = 1000, optional.
Enables specification of the number of
permutations to be generated to estimate the empirical null
distribution of maximal component size.

Returns
-------
PVAL : ndarray
p-values for each component
ADJ : ndarray
Returns an adjacency matrix identifying the edges comprising each component.
Edges corresponding to the first p-value stored in the vector PVAL are assigned
the value 1 in the adjacency matrix ADJ, edges corresponding to the second
p-value are assigned the value 2, etc.

Modified from the original so that only components surviving a given threshold are exported

NULL : ndarray
Returns a vector of K samples
from the the null distribution of maximal component size.

ALGORITHM DESCRIPTION
The NBS is a nonparametric statistical test used to isolate the
components of an N x N undirected connectivity matrix that differ
significantly between two distinct populations. Each element of the
connectivity matrix stores a connectivity value and each member of
the two populations possesses a distinct connectivity matrix. A
component of a connectivity matrix is defined as a set of
interconnected edges.

The NBS is essentially a procedure to control the family-wise error
rate, in the weak sense, when the null hypothesis is tested
independently at each of the N(N-1)/2 edges comprising the
connectivity matrix. The NBS can provide greater statistical power
than conventional procedures for controlling the family-wise error
rate, such as the false discovery rate, if the set of edges at which
the null hypothesis is rejected constitues a large component or
components.
The NBS comprises fours steps:
1. Perfrom a two-sample T-test at each edge indepedently to test the
hypothesis that the value of connectivity between the two
populations come from distributions with equal means.
2. Threshold the T-statistic available at each edge to form a set of
suprathreshold edges.
3. Identify any components in the adjacency matrix defined by the set
of suprathreshold edges. These are referred to as observed
components. Compute the size of each observed component
identified; that is, the number of edges it comprises.
4. Repeat K times steps 1-3, each time randomly permuting members of
the two populations and storing the size of the largest component
identified for each permuation. This yields an empirical estimate
of the null distribution of maximal component size. A corrected
p-value for each observed component is then calculated using this
null distribution.

[1] Zalesky A, Fornito A, Bullmore ET (2010) Network-based statistic:
Identifying differences in brain networks. NeuroImage.
10.1016/j.neuroimage.2010.06.041

Written by: Andrew Zalesky, azalesky@unimelb.edu.au
Rewritten for Python: Stephan Gerhard, connectome@unidesign.ch

"""

    print d_stacked.shape
    
    assert d_stacked.shape[2] == nx + ny
    
    t1 = time.time()
    
    max_sz,sz_links,nbs_adj_mat_file = find_biggest_component_size_binom_save(d_stacked[:,:,:nx],d_stacked[:,:,nx:nx+ny],conf_interval_binom)

    print 'save nbs adj mat file' 
    nbs_adj_mat_file= os.path.abspath('nbs_adj_matrix_'+ str(conf_interval_binom) +'.npy')
    np.save(nbs_adj_mat_file,ADJ)
    
    
    print 'save nbs sz_links file'
    nbs_stats_file  = os.path.abspath('nbs_sz_links_'+ str(conf_interval_binom) +'.txt')
    np.savetxt(nbs_stats_file,sz_links,fmt = '%f')
    
    t2 = time.time()
    
    #if False:
        ## additionally, store all the components in the matrix with the value of their number of edges
        #all_components = np.zeros( (N,N) )
        #for idx, componentG in enumerate(comp_list):
            
            #tmp_max = netwx.to_numpy_matrix( componentG , nodelist = range(N) )
            ## set nonzero to number of edges
            #tmp_max[tmp_max!=0.0] = componentG.number_of_edges()
            #all_components[:,:] = all_components[:,:] + tmp_max
    
    print "Max component size is: %s, time = %d" % (max_sz,int(t2-t1))
        
    # Empirically estimate null distribution of maximum component size by
    # generating K independent permutations.
    print "====================================================="
    print "Estimating null distribution with permutation testing"
    print "====================================================="
    
    hit=0.0
    NULL = np.zeros( (K, 1) )
    # stack matrices for permutation
   
    
    for k in range(K):
        
        t1 = time.time()
        
        # Randomize
        indperm = np.random.permutation( nx+ny )
        
        
        d_perm = d_stacked[:,:, indperm]

        #################
        
        max_sz_perm,sz_links_perm = find_biggest_component_size_binom(d_perm[:,:,:nx],d_perm[:,:,nx:nx+ny],conf_interval_binom)
        
        NULL[k] = max_sz_perm

        # if the component size of this random permutation is bigger than
        # the component size of the group difference computed above, this is a hit
        if NULL[k] >= max_sz:
            hit = hit + 1
            
        
        t2 = time.time()
        
        print "Perm %d of %d (took %d s). Perm max is: %d. Observed max is: %d. P-val estimate is: %0.3f" % ((k+1),K,int(t2-t1),NULL[k],max_sz,hit/(k+1))

    # Calculate p-values for each component
    PVAL = np.zeros( len(sz_links) )
    for i in range( len(sz_links) ):
        PVAL[i] = len( NULL[NULL >= sz_links[i]] ) * 1.0 / K
        
        
    
    ## Suprathreshold adjacency matrix
    #ADJ_signif = np.zeros( (N,N) )
    
    #for idx, componentG in enumerate(comp_list):
        #nr_edges_per_component[idx] = componentG.number_of_edges()
        
        ## if number of edges bigger than zero
        #if PVAL[idx] > 0.1:
            
            ## store the number of edges of the component as value in the adjacency matrix
            #for ed in componentG.edges():
                #ADJ_signif[ed[0],ed[1]] = ADJ_signif[ed[1],ed[0]] = idx + 1
                ## if we would like to store the number of edges per component
                ## ADJ_signif[ed[0],ed[1]] = ADJ_signif[ed[1],ed[0]] = nr_edges_per_component[idx]
    
        
    return (PVAL, NULL, nbs_adj_mat_file)
        
    #return (PVAL, ADJ, ADJ_signif, NULL, sz_links)
