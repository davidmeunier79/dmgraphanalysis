# -*- coding: utf-8 -*-

import scipy.stats as stat
import numpy as np
import itertools as it

def info_CI(X,Y):
    """ Compute binomial comparaison
"""
    nX = len(X) * 1.
    nY = len(Y) * 1.
    
    pX = np.sum(X == 1)/nX

    pY = np.sum(Y == 1)/nY

    #print pX,pY,np.absolute(pX-pY) 

    SE = np.sqrt(pX * (1-pX)/nX + pY * (1-pY)/nY)

    #print SE
    
    #if (np.absolute(pX-pY) > norm) == True:
        #print pX,pY,np.absolute(pX-pY),norm
   
    
    return np.absolute(pX-pY),SE,np.sign(pX-pY)

    ###################################################################################### pairwise/nodewise stats ###########################################################################################################
    
    
def return_signif_bin_vect(p_values,fdr_alpha = 0.05):

    print p_values
    
    
    N = p_values.shape[0]
    
    
    
    order =  p_values.argsort()
    
    init_order = range(N)
    
    
    print order
    
    sorted_p_values= p_values[order]
    
    print sorted_p_values
    
    ### by default, code = 1 (cor at 0.05)
    sorted_order = np.ones(shape = N)
    
    ################ uncor #############################
    
    sorted_order[p_values > fdr_alpha] = 0
    
    ################ fdr ###############################
    seq = np.arange(N,0,-1)
    
    seq_fdr_p_values = fdr_alpha/seq
    
    print seq_fdr_p_values
    
    sorted_order[sorted_p_values < seq_fdr_p_values] = 2
    
    ################# bonferroni #######################
    
    sorted_order[sorted_p_values < fdr_alpha/N] = 4
    
    return sorted_order
    
def return_signif_pval_vect(sort_np_list_diff,t_test_thresh_fdr):

    n = sort_np_list_diff.shape[0]
    
    ############### uncor #############################
    
    sort_np_list_diff[sort_np_list_diff[:,1] > t_test_thresh_fdr,2] = 0
    
    ############### fdr ###############################
    seq = np.arange(n,0,-1)
    
    #print seq
    
    seq_fdr_p_values = t_test_thresh_fdr/seq
    
    sort_np_list_diff[sort_np_list_diff[:,1] < seq_fdr_p_values,2] = sort_np_list_diff[sort_np_list_diff[:,1] < seq_fdr_p_values,2] * 2
    
    ################ bonferroni #######################
    
    sort_np_list_diff[sort_np_list_diff[:,1] < t_test_thresh_fdr/n,2] = sort_np_list_diff[sort_np_list_diff[:,1] < t_test_thresh_fdr/n,2] * 2
    
    return sort_np_list_diff
    
    
def return_signif_pval_mat(sort_np_list_diff,t_test_thresh_fdr):

    n = sort_np_list_diff.shape[0]
    
    ############### uncor #############################
    
    sort_np_list_diff[sort_np_list_diff[:,2] > t_test_thresh_fdr,3] = 0
    
    ############### fdr ###############################
    seq = np.arange(n,0,-1)
    
    #print seq
    
    seq_fdr_p_values = t_test_thresh_fdr/seq
    
    sort_np_list_diff[sort_np_list_diff[:,2] < seq_fdr_p_values,3] = sort_np_list_diff[sort_np_list_diff[:,2] < seq_fdr_p_values,3] * 2
    
    ################ bonferroni #######################
    
    sort_np_list_diff[sort_np_list_diff[:,2] < t_test_thresh_fdr/n,3] = sort_np_list_diff[sort_np_list_diff[:,2] < t_test_thresh_fdr/n,3] * 2
    
    return sort_np_list_diff
    
    
    
    
    
def return_signif_binom_mat(sort_np_list_diff,conf_interval_binom_fdr):

    n = sort_np_list_diff.shape[0]
    
    ############ uncor
    
    uncor_Z_val = stat.norm.ppf(1-conf_interval_binom_fdr/2)
    
    uncor_norm = uncor_Z_val * sort_np_list_diff[:,3]
    
    sort_np_list_diff[sort_np_list_diff[:,2] < uncor_norm,4] = 0
    
    ############ fdr
    
    seq = np.arange(n,0,-1)
    
    print seq
    
    seq_fdr_p_values = conf_interval_binom_fdr/seq
    
    print seq_fdr_p_values
    #seq_fdr_p_values
    
    seq_Z_val = stat.norm.ppf(1-seq_fdr_p_values/2)
    
    print 
    
    
    print seq_Z_val
    
    seq_norm =  seq_Z_val * sort_np_list_diff[:,3]
    
    print seq_norm
    
    print np.sum(np.array(sort_np_list_diff[:,2] > seq_norm,dtype = 'int'))
    
    sort_np_list_diff[sort_np_list_diff[:,2] > seq_norm,4] = sort_np_list_diff[sort_np_list_diff[:,2] > seq_norm,4] * 2
    
    ################# bonferroni
    
    bon_Z_val = stat.norm.ppf(1-conf_interval_binom_fdr/(2*n))
    
    bon_norm = bon_Z_val * sort_np_list_diff[:,3]
    
    sort_np_list_diff[sort_np_list_diff[:,2] > bon_norm,4] = sort_np_list_diff[sort_np_list_diff[:,2] > bon_norm,4] * 2
    
    return sort_np_list_diff
    
    
    
    #signif_signed_vect = np.zeros((N),dtype = 'int')
    
    #signif_signed_vect[np.array(signif_fdr[:,0],dtype = int)] = np.array(signif_fdr[:,2],dtype = int)
    
    #print signif_signed_vect
    
    #return signif_signed_vect

    
    
    
#def return_signif_fdr_values_pval_vect(sort_np_list_diff):

    #n = sort_np_list_diff.shape[0]
    
    #seq = np.arange(n,0,-1)
    
    #print seq
    
    #seq_fdr_p_values = t_test_thresh_fdr/seq
    
    #print np.sum(np.array(sort_np_list_diff[:,1] < seq_fdr_p_values,dtype = 'int'))
    
    #print sort_np_list_diff[sort_np_list_diff[:,1] < seq_fdr_p_values,:]
    
    #return sort_np_list_diff[sort_np_list_diff[:,1] < seq_fdr_p_values,:]
    
#def return_signif_fdr_values_binom_mat(sort_np_list_diff):

    #n = sort_np_list_diff.shape[0]
    
    #seq = np.arange(n,0,-1)
    
    #print seq
    
    #seq_fdr_p_values = conf_interval_binom_fdr/seq
    
    #print seq_fdr_p_values
    ##seq_fdr_p_values
    
    #seq_Z_val = stat.norm.ppf(1-seq_fdr_p_values/2)
    
    #print 
    
    
    #print seq_Z_val
    
    #seq_norm =  Z_val * sort_np_list_diff[:,3]
    
    #print seq_norm
    
    #print np.sum(np.array(sort_np_list_diff[:,2] > seq_norm,dtype = 'int'))
    
    #print sort_np_list_diff[sort_np_list_diff[:,2] > seq_norm,:]
    
    #return sort_np_list_diff[sort_np_list_diff[:,2] > seq_norm,:]
    
    
def compute_pairwise_binom(X,Y,conf_interval_binom):

    # number of nodes
    N = X.shape[0]
   
    # Perform binomial test at each edge
    ADJ = np.zeros((N,N),dtype = 'int')
    
    for i,j in it.combinations(range(N), 2):
        
        ADJ[i,j] = ADJ[j,i] = binom_CI_test(X[i,j,:],Y[i,j,:],conf_interval_binom)
        
    return ADJ

    
def compute_pairwise_ttest_fdr(X,Y,t_test_thresh_fdr):
    
    # number of nodes
    N = X.shape[0]
   
    list_diff = []
    
    for i,j in it.combinations(range(N), 2):
        
        #t_stat_zalewski = ttest2(X[i,j,:],Y[i,j,:])
        
        t_stat,p_val = stat.ttest_ind(X[i,j,:],Y[i,j,:])
        
        print t_stat,p_val
        
        list_diff.append([i,j,p_val,np.sign(np.mean(X[i,j,:])-np.mean(Y[i,j,:]))])
        
    #print list_diff
        
    np_list_diff = np.array(list_diff)
   
    print np_list_diff
    
    print np_list_diff
    
    order =  np_list_diff[:,2].argsort()
    
    print order
    
    sort_np_list_diff = np_list_diff[order[::-1]]
    
    print sort_np_list_diff.shape[0]
    
    signif_fdr = return_signif_pval_mat(sort_np_list_diff,t_test_thresh_fdr)
    
    #print signif_fdr
    
    
    print np.sum(signif_fdr[:,3] == 0.0),np.sum(signif_fdr[:,3] == 1.0),np.sum(signif_fdr[:,3] == 2.0),np.sum(signif_fdr[:,3] == 4.0),np.sum(signif_fdr[:,3] == -1.0),np.sum(signif_fdr[:,3] == -2.0),np.sum(signif_fdr[:,3] == -4.0)
    
    signif_signed_adj_mat = np.zeros((N,N),dtype = 'int')
    
        
    signif_i = np.array(signif_fdr[:,0],dtype = int)
    signif_j = np.array(signif_fdr[:,1],dtype = int)
    
    signif_sign = np.array(signif_fdr[:,3],dtype = int)
    
    print signif_i,signif_j
    
    print signif_signed_adj_mat[signif_i,signif_j] 
    
    #print signif_sign
    
    
    
    signif_signed_adj_mat[signif_i,signif_j] = signif_signed_adj_mat[signif_j,signif_i] = signif_sign
    
    print signif_signed_adj_mat
    
    return signif_signed_adj_mat

    
def compute_pairwise_binom_fdr(X,Y,conf_interval_binom_fdr):

    # number of nodes
    N = X.shape[0]
   
    # Perform binomial test at each edge
    
    
    list_diff = []
    
    for i,j in it.combinations(range(N), 2):
        
        abs_diff,SE,sign_diff = info_CI(X[i,j,:],Y[i,j,:])
         
        list_diff.append([i,j,abs_diff,SE,sign_diff])
        
    print list_diff
    
    np_list_diff = np.array(list_diff)
    
    
    print np_list_diff
    
    order =  np_list_diff[:,2].argsort()
    
    print order
    
    sort_np_list_diff = np_list_diff[order[::-1]]
    
    
    print sort_np_list_diff
    
    print sort_np_list_diff.shape[0]
    
    signif_fdr = return_signif_binom_mat(sort_np_list_diff,conf_interval_binom_fdr)
    
    
    print np.sum(signif_fdr[:,4] == 0.0),np.sum(signif_fdr[:,4] == 1.0),np.sum(signif_fdr[:,4] == 2.0),np.sum(signif_fdr[:,4] == 4.0),np.sum(signif_fdr[:,4] == -1.0),np.sum(signif_fdr[:,4] == -2.0),np.sum(signif_fdr[:,4] == -4.0)
    
        
    signif_signed_adj_mat = np.zeros((N,N),dtype = 'int')
    
    #signif_indexes = np.array(signif_fdr[:,0:2],dtype = int)
    
    signif_i = np.array(signif_fdr[:,0],dtype = int)
    signif_j = np.array(signif_fdr[:,1],dtype = int)
    
    signif_sign = np.array(signif_fdr[:,4],dtype = int)
    
    print signif_i,signif_j
    
    print signif_signed_adj_mat[signif_i,signif_j] 
    
    #print signif_sign
    
    
    
    signif_signed_adj_mat[signif_i,signif_j] = signif_signed_adj_mat[signif_j,signif_i] = signif_sign
    
    
    print signif_signed_adj_mat
    #signif_pos_mat_indexes = signif_fdr[signif_fdr[:,4] == 1,0:2]
    
    #print signif_pos_mat_indexes
    
    return signif_signed_adj_mat

#def compute_nodewise_t_values(X,Y):

    ## number of nodes
    #N = X.shape[0]
   
    ## Perform binomial test at each edge
    #t_val_vect = np.zeros((N),dtype = 'float')
    
    #for i in range(N):
        
        #t_val_vect[i] = ttest2(X[i,:],Y[i,:])
        
    #return t_val_vect

    
def compute_nodewise_t_values_fdr(X,Y,t_test_thresh_fdr):

    # number of nodes
    N = X.shape[0]
   
    list_diff = []
    
    for i in range(N):
        
        #t_stat_zalewski = ttest2(X[i,:],Y[i,:])
        
        t_stat,p_val = stat.ttest_ind(X[i,:],Y[i,:])
        
        print t_stat,p_val
        
        list_diff.append([i,p_val,np.sign(np.mean(X[i,:])-np.mean(Y[i,:]))])
        
    #print list_diff
        
    np_list_diff = np.array(list_diff)
   
    print np_list_diff
    
    sort_np_list_diff = np_list_diff[np_list_diff[:,1].argsort()]
    
    print sort_np_list_diff.shape[0]
    
    signif_fdr = return_signif_pval_vect(sort_np_list_diff,t_test_thresh_fdr)
    
    #print signif_fdr
    
    
    print np.sum(signif_fdr[:,2] == 0.0),np.sum(signif_fdr[:,2] == 1.0),np.sum(signif_fdr[:,2] == 2.0),np.sum(signif_fdr[:,2] == 4.0),np.sum(signif_fdr[:,2] == -1.0),np.sum(signif_fdr[:,2] == -2.0),np.sum(signif_fdr[:,2] == -4.0)
            
    signif_signed_vect = np.zeros((N),dtype = 'int')
    
    signif_signed_vect[np.array(signif_fdr[:,0],dtype = int)] = np.array(signif_fdr[:,2],dtype = int)
    
    print signif_signed_vect
    
    return signif_signed_vect

    
    
    
    
#def compute_pairwise_binom_adj_mat(d_stacked, nx, ny):

    #print d_stacked.shape
    
    #assert d_stacked.shape[2] == nx + ny
    
    #t1 = time.time()
    
    #binom_adj_mat = compute_pairwise_binom(d_stacked[:,:,:nx],d_stacked[:,:,nx:nx+ny])

    #t2 = time.time()
    
    #print "computation took %f" %(t2-t1)
    
    #return binom_adj_mat
    
    
#def compute_pairwise_binom_fdr_sign_adj_mat(d_stacked, nx, ny):

    #print d_stacked.shape
    
    #assert d_stacked.shape[2] == nx + ny
    
    #t1 = time.time()
    
    #signif_pos_adj_mat,signif_neg_adj_mat = compute_pairwise_binom_fdr(d_stacked[:,:,:nx],d_stacked[:,:,nx:nx+ny])

    
    #t2 = time.time()
    
    #print "computation took %f" %(t2-t1)
    
    #return  signif_pos_adj_mat,signif_neg_adj_mat
    
    
    
def compute_nodewise_t_test_vect(d_stacked, nx, ny):

    print d_stacked.shape
    
    assert d_stacked.shape[1] == nx + ny
    
    t1 = time.time()
    
    t_val_vect = compute_nodewise_t_values(d_stacked[:,:nx],d_stacked[:,nx:nx+ny])

    t2 = time.time()
    
    print "computation took %f" %(t2-t1)
    
    return t_val_vect
    
######################## correl ######################################

def compute_pairwise_correl_fdr(X,behav_score,correl_thresh_fdr):


    from scipy.stats.stats import pearsonr

    # number of nodes
    N = X.shape[0]
   
    list_diff = []
    
    for i,j in it.combinations(range(N), 2):
        
        #t_stat_zalewski = ttest2(X[i,j,:],Y[i,j,:])
        
        r_stat,p_val = pearsonr(X[i,j,:],behav_score)
        
        #print i,j,p_val,r_stat
        
        list_diff.append([i,j,p_val,np.sign(r_stat)])
        
    #print list_diff
        
    np_list_diff = np.array(list_diff)
   
    #print np_list_diff
    
    order =  np_list_diff[:,2].argsort()
    
    #print order
    
    sort_np_list_diff = np_list_diff[order[::-1]]
    
    #print sort_np_list_diff.shape[0]
    
    signif_fdr = return_signif_pval_mat(sort_np_list_diff,correl_thresh_fdr)
    
    #print signif_fdr
    
    
    print np.sum(signif_fdr[:,3] == 0.0),np.sum(signif_fdr[:,3] == 1.0),np.sum(signif_fdr[:,3] == 2.0),np.sum(signif_fdr[:,3] == 4.0),np.sum(signif_fdr[:,3] == -1.0),np.sum(signif_fdr[:,3] == -2.0),np.sum(signif_fdr[:,3] == -4.0)
    
    signif_signed_adj_mat = np.zeros((N,N),dtype = 'int')
    
        
    signif_i = np.array(signif_fdr[:,0],dtype = int)
    signif_j = np.array(signif_fdr[:,1],dtype = int)
    
    signif_sign = np.array(signif_fdr[:,3],dtype = int)
    
    print signif_i,signif_j
    
    print signif_signed_adj_mat[signif_i,signif_j] 
    
    #print signif_sign
    
    
    
    signif_signed_adj_mat[signif_i,signif_j] = signif_signed_adj_mat[signif_j,signif_i] = signif_sign
    
    print signif_signed_adj_mat
    
    return signif_signed_adj_mat