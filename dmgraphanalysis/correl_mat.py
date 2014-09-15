# -*- coding: utf-8 -*-

   
import nibabel as nib
import nipy as nip



### script Thirion
import os 

import numpy as np

#from nibabel import load, save

#import nipy.labs.spatial_models.mroi as mroi
#from nipy.labs.spatial_models.discrete_domain import grid_domain_from_image
#import nipy.labs.spatial_models.hroi as hroi

#import nipy.labs.statistical_mapping as stat_map

import itertools as iter
    
import scipy.spatial.distance as dist
    
########################################### extract time series ################################################

def return_ts(orig_ts,peak_position):

    #print orig_ts.shape
    
    if check_dimensions(peak_position,orig_ts.shape[:-1]):
    
        peak_ts = np.array(orig_ts[peak_position[0],peak_position[1],peak_position[2],:],dtype = 'int64')
        
    else:
        print "Warning, position %d %d %d out of indexes, returning 0 ts " %(peak_position[0],peak_position[1],peak_position[2])
    
        peak_ts = np.zeros(orig_ts.shape[3],dtype = 'int64')
    
    #print peak_ts
    
    return peak_ts
    
def return_mean_ts(orig_ts,peak_position,neighbourhood = 1):

    neigh_range = range(-neighbourhood,neighbourhood+1)
    
    list_neigh_ts = []
    list_neigh_coords = []
    
    for relative_coord in iter.product(neigh_range, repeat=3):

        neigh_coord = peak_position + relative_coord
        
        list_neigh_coords.append(neigh_coord)
        list_neigh_ts.append(return_ts(orig_ts,neigh_coord))
        
        
    neigh_coords = np.array(list_neigh_coords,dtype = 'int16')
    
    neigh_ts = np.array(list_neigh_ts,dtype = 'f')
    
    mean_ts = np.mean(neigh_ts,axis = 0)
    
    return mean_ts,neigh_coords

def compute_mean_ts_from_labelled_mask(file_4D,indexed_rois_file,coord_rois_file,min_BOLD_intensity = 50):
    
    
    import os
    import numpy as np
    import nibabel as nib
    
    from dmgraphanalysis.utils_plot import plot_signals
    
    
    ## loading ROI coordinates
    coord_rois = np.loadtxt(coord_rois_file)
    
    #print "coord_rois: " 
    #print coord_rois.shape
    
    ## loading ROI indexed mask
    indexed_rois_img = nib.load(indexed_rois_file)
    
    indexed_mask_rois_data = indexed_rois_img.get_data()
    
    #print "indexed_mask_rois_data: "
    #print indexed_mask_rois_data.shape
    
    ### loading time series
    orig_ts = nib.load(file_4D).get_data()
    
    #print "orig_ts shape:"
    #print orig_ts.shape
        
    ### extrating ts by averaging the time series of all voxel with the same index
    sequence_roi_index = np.array(np.unique(indexed_mask_rois_data),dtype = int)
    
    if sequence_roi_index[0] == -1.0:
        sequence_roi_index = sequence_roi_index[1:]
    
    #print "sequence_roi_index:"
    #print sequence_roi_index
    
    if sequence_roi_index.shape[0] != coord_rois.shape[0]:
        print "Warning, indexes in template_ROI are incompatible with ROI coords"
        return
    
    mean_masked_ts = []
    subj_coord_rois = []
    
    for roi_index in sequence_roi_index:
    #for roi_index in sequence_roi_index[0:1]:
        
        #print "Roi index " + str(roi_index)
        
        index_roi_x,index_roi_y,index_roi_z = np.where(indexed_mask_rois_data == roi_index)
        
        #print index_roi_x,index_roi_y,index_roi_z
        
        all_voxel_roi_ts = orig_ts[index_roi_x,index_roi_y,index_roi_z,:]
        
        #print "all_voxel_roi_ts shape: "
        #print all_voxel_roi_ts.shape
        
        mean_roi_ts = np.mean(all_voxel_roi_ts,axis = 0)
        
        ### testing if mean_roi_ts if higher than minimal BOLD intensity
        
        if sum(mean_roi_ts > min_BOLD_intensity) == mean_roi_ts.shape[0]:
            
            #print "Roi selected: " + str(roi_index)
            #print "coord_rois shape: "
            #print coord_rois.shape
            
            subj_coord_rois.append(coord_rois[roi_index,])
            mean_masked_ts.append(mean_roi_ts)
        else:
            print "ROI " + str(roi_index) + " was not selected"
            
        
    mean_masked_ts = np.array(mean_masked_ts,dtype = 'f')
    subj_coord_rois = np.array(subj_coord_rois,dtype = 'float')
    
    print mean_masked_ts.shape
        
    ### saving time series
    mean_masked_ts_file = os.path.abspath("mean_masked_ts.txt")
    np.savetxt(mean_masked_ts_file,mean_masked_ts,fmt = '%.3f')
    
    ### saving subject ROIs
    subj_coord_rois_file = os.path.abspath("subj_coord_rois.txt")
    np.savetxt(subj_coord_rois_file,subj_coord_rois,fmt = '%.3f')
    
    
    print "plotting mean_masked_ts"
    
    plot_mean_masked_ts_file = os.path.abspath('mean_masked_ts.eps')    
    
    plot_signals(plot_mean_masked_ts_file,mean_masked_ts)
    
    
    return mean_masked_ts_file,subj_coord_rois_file
    
def mean_select_ts_with_mask(file_4D,mask_file,suffix):

    import os
    
    import nibabel as nib
    import numpy as np 

    from dmgraphanalysis.utils_cor import mean_select_mask_data
    from dmgraphanalysis.utils_plot import plot_signals
    
    print 'in select_ts_with_mask'
    
    print "loading img data " + file_4D

    ### Reading 4D volume file to extract time series
    img = nib.load(file_4D)
    img_data = img.get_data()
    
    print img_data.shape

    
    print "loading mask data " + mask_file

    ### Reading 4D volume file to extract time series
    mask_data = nib.load(mask_file).get_data()
    
    print mask_data.shape

    print "mean_select_mask_data"
    
    ### Retaining only time series who are within the mask + non_zero
    mean_masked_ts = mean_select_mask_data(img_data,mask_data)
    
    print "saving mean_masked_ts"
    mean_masked_ts_file = os.path.abspath('mean_' + suffix + '_ts.txt')    
    np.savetxt(mean_masked_ts_file,mean_masked_ts,fmt = '%.3f')
    
    
    print "plotting mean_masked_ts"
    
    plot_mean_masked_ts_file = os.path.abspath('mean_' + suffix + '_ts.eps')    
    
    plot_signals(plot_mean_masked_ts_file,mean_masked_ts)
    
    return mean_masked_ts_file
    
############################ functions used in workflow ##############################################################

def regress_mvt_param(masked_ts_file,rp_file):

    import os
    
    import nibabel as nib
    import numpy as np
    
    #from scipy import stats

    from dmgraphanalysis.utils_cor import regress_movement_parameters
    
    print "load masked_ts_file"
    
    data_mask_matrix = np.loadtxt(masked_ts_file)
    
    print data_mask_matrix.shape

    ### Reading movement parameters 
    print "load rp parameters"
    
    print rp_file
    
    rp = np.genfromtxt(rp_file)
    #rp = np.loadtxt(rp_file,dtype = np.float)
    
    #print rp.shape
    
    #### remove first values (first scans are discarded)
    #rp = rp[nb_scans_to_remove:,]
    
    print rp.shape
    
    print "regress rp parameters"
    
    ### regression movement parameters and computing z-score on the residuals
    resid_data_matrix = regress_movement_parameters(data_mask_matrix,rp)
    print resid_data_matrix.shape

    print "saving resid_ts"
    
    resid_ts_file = os.path.abspath('resid_ts.npy')
    np.save(resid_ts_file,resid_data_matrix )

    #return resid_data_matrix,resid_ts_file
    return resid_ts_file

def regress_covariates(masked_ts_file,rp_file,mean_wm_ts_file,mean_csf_ts_file):

    import rpy,os
    
    import nibabel as nib
    import numpy as np
    
    #from dmgraphanalysis.utils_cor import regress_movement_wm_csf_parameters
    from dmgraphanalysis.utils_cor import regress_movement_wm_csf_parameters_and_filter
    from dmgraphanalysis.utils_plot import plot_signals,plot_sep_signals
    
    print "load mean_csf_ts_file" + str(mean_csf_ts_file)
    
    mean_csf_ts = np.loadtxt(mean_csf_ts_file)
    
    print mean_csf_ts
    
    print "load mean_wm_ts_file"
    
    mean_wm_ts = np.loadtxt(mean_wm_ts_file)
    
    print "load masked_ts_file"
    
    data_mask_matrix = np.loadtxt(masked_ts_file)
    
    print data_mask_matrix.shape

    print "load rp parameters"
    
    print rp_file
    
    rp = np.genfromtxt(rp_file)
    #rp = np.loadtxt(rp_file,dtype = np.float)
    
    print rp.shape
    
    
    
    ### was used in ua_avec_image1
    #print "remove rp first values (discarded scans)"
    
    #### remove first values (first scans are discarded)
    #rp = rp[nb_scans_to_remove:,]
    
    #print rp.shape
    
    print "regress rp parameters"
    
    ### regression movement parameters and computing z-score on the residuals
    #resid_data_matrix = regress_movement_wm_csf_parameters(data_mask_matrix,rp,mean_wm_ts,mean_csf_ts)
    resid_data_matrix,resid_filt_data_matrix,z_score_data_matrix = regress_movement_wm_csf_parameters_and_filter(data_mask_matrix,rp,mean_wm_ts,mean_csf_ts)
    
    print resid_data_matrix.shape

    print "saving resid_ts"
    
    resid_ts_file = os.path.abspath('resid_ts.npy')
    np.save(resid_ts_file,z_score_data_matrix )

    print "plotting resid_ts"
    
    plot_resid_ts_file = os.path.abspath('resid_ts.eps')
    
    plot_sep_signals(plot_resid_ts_file,z_score_data_matrix)
    
    
    print "plotting diff filtered and non filtered data"
    
    plot_diff_filt_ts_file = os.path.abspath('diff_filt_ts.eps')
    
    plot_signals(plot_diff_filt_ts_file,np.array(resid_filt_data_matrix - resid_data_matrix,dtype = 'float'))
    
    #plot_sep_signals(plot_resid_ts_file,resid_data_matrix)
    
    #return resid_data_matrix,resid_ts_file
    return resid_ts_file
    
#### find regressor corresponding to condition and run_index ####
def return_regressor_SPM(spm_mat_file,regressor_name,run_index):

    import scipy.io
    import numpy as np
    import os

    #print spm_mat_file
    
    ##Reading spm.mat for regressors extraction:
    d = scipy.io.loadmat(spm_mat_file)
    
    #print d
    
    
    ##Choosing the column according to the regressor name
    #_,col = np.where(d['SPM']['xX'][0][0]['name'][0][0] == u'Sn(1) ' + regressor_name)
    
    cond_name = u'Sn(' + str(run_index) + ') ' + regressor_name+'*bf(1)'
    
    print cond_name
    
    _,col = np.where(d['SPM']['xX'][0][0]['name'][0][0] == cond_name)
    
    #print col
    
    ## reformating matrix (1,len) in vector (len)
    regressor_vect = d['SPM']['xX'][0][0]['X'][0][0][:,col].reshape(-1)
    

    #print regressor_vect
    
    regressor_vect[regressor_vect < 0] = 0
    
    print "Saving extract_cond"
    regressor_file = os.path.abspath('extract_cond.txt')

    np.savetxt(regressor_file,regressor_vect)

    return regressor_file

def merge_several_runs_coords(resid_ts_files,regressor_files,coord_rois_files):
    
    import numpy as np
    import os
    
    if len(resid_ts_files) != len(regressor_files):
        
        print "Warning, time series and regressors have different length (!= number of runs)"
        return 0
        
    if len(resid_ts_files) != len(coord_rois_files):
        
        print "Warning, time series and number of coordinates have different length (!= number of runs)"
        return 0
          
    ### concatenate time series
    for i,resid_ts_file in enumerate(resid_ts_files):
        
        resid_data_matrix = np.load(resid_ts_file)
        
        print resid_data_matrix.shape
        
        ## loading ROI coordinates
        coord_rois = np.loadtxt(coord_rois_files[i])
        
        print coord_rois.shape
        
        if i == 0:
            resid_data_matrix_all_runs = np.empty((resid_data_matrix.shape[0],0),dtype = resid_data_matrix.dtype)
            
            coord_rois_all_runs = np.array(coord_rois,dtype = 'float')
            
            
        if coord_rois_all_runs.shape[0] != coord_rois.shape[0]:
            
            print "ROIs do not match for all different sessions "
            
            print os.getcwd()
            
            print "Warning, not implemented yet.... "
            
            ### pris de http://stackoverflow.com/questions/8317022/get-intersecting-rows-across-two-2d-numpy-arrays
            ### à tester....
            ### finir également la partie avec resid_data_matrix_all_runs, en supprimant les colonnes qui ne sont pas communes à tous les runs...
            
            0/0
            
            A = coord_rois_all_runs
            B = coord_rois
            
            nrows, ncols = A.shape
            dtype={'names':['f{}'.format(i) for i in range(ncols)],
                'formats':ncols * [A.dtype]}

            C = np.intersect1d(A.view(dtype), B.view(dtype))

            # This last bit is optional if you're okay with "C" being a structured array...
            C = C.view(A.dtype).reshape(-1, ncols)

            coord_rois_all_runs = C

        resid_data_matrix_all_runs = np.concatenate((resid_data_matrix_all_runs,resid_data_matrix),axis = 1)
        
        print resid_data_matrix_all_runs.shape
        
    ### save times series for all runs
    resid_ts_all_runs_file = os.path.abspath('resid_ts_all_runs.npy')
     
    np.save(resid_ts_all_runs_file,resid_data_matrix_all_runs)
     
    ### save coords in common for all runs
    coord_rois_all_runs_file = os.path.abspath('coord_rois_all_runs.txt')
     
    np.savetxt(coord_rois_all_runs_file,coord_rois_all_runs, fmt = '%2.3f')
    
    
    
    
    regress_data_vector_all_runs = np.empty(shape = (0), dtype = float)
    
     ### Sum regressors
    for i,regress_file in enumerate(regressor_files):
        
        regress_data_vector = np.loadtxt(regress_file)
        
        #print regress_data_vector.shape
        
        if regress_data_vector.shape[0] != 0:
            
            if regress_data_vector_all_runs.shape[0] == 0:
                
                regress_data_vector_all_runs = regress_data_vector
            else:
                regress_data_vector_all_runs = regress_data_vector_all_runs + regress_data_vector
            
        print np.sum(regress_data_vector_all_runs != 0.0)
    
    regress_data_vector_all_runs_file = os.path.abspath('regress_data_vector_all_runs.txt')

    np.savetxt(regress_data_vector_all_runs_file,regress_data_vector_all_runs)

    return resid_ts_all_runs_file,regress_data_vector_all_runs_file,coord_rois_all_runs_file
    
def merge_several_runs(resid_ts_files,regressor_files):
    
    import numpy as np
    import os
    
    if len(resid_ts_files) != len(regressor_files):
        
        print "Warning, time series and regressors have different length (!= number of runs)"
        return 0
        
    ### concatenate time series
    for i,resid_ts_file in enumerate(resid_ts_files):
        
        resid_data_matrix = np.load(resid_ts_file)
        
        print resid_data_matrix.shape
        
        if i == 0:
            print resid_data_matrix.shape[0]
            
            print resid_data_matrix.dtype
            
            resid_data_matrix_all_runs = np.empty((resid_data_matrix.shape[0],0),dtype = resid_data_matrix.dtype)
            
            print resid_data_matrix_all_runs
            
            print resid_data_matrix_all_runs.shape
        
        resid_data_matrix_all_runs = np.concatenate((resid_data_matrix_all_runs,resid_data_matrix),axis = 1)
        
        print resid_data_matrix_all_runs.shape
        
    resid_ts_all_runs_file = os.path.abspath('resid_ts_all_runs.npy')
     
    np.save(resid_ts_all_runs_file,resid_data_matrix_all_runs)
     
     
     
     
     ### Sum regressors
    for i,regress_file in enumerate(regressor_files):
        
        regress_data_vector = np.loadtxt(regress_file)
        
        #print regress_data_vector.shape
        
        if i == 0:
            regress_data_vector_all_runs = regress_data_vector
        else:
            regress_data_vector_all_runs = regress_data_vector_all_runs + regress_data_vector
        
        print np.sum(regress_data_vector_all_runs != 0.0)
    
    regress_data_vector_all_runs_file = os.path.abspath('regress_data_vector_all_runs.txt')

    np.savetxt(regress_data_vector_all_runs_file,regress_data_vector_all_runs)

    return resid_ts_all_runs_file,regress_data_vector_all_runs_file
    
################################## Weighted Correlation matrix computation #########################################
def compute_var_correlation_matrix(resid_ts_file,regressor_file):
    
    import rpy,os
    import nibabel as nib
    import numpy as np
    
    from dmgraphanalysis.utils_cor import return_var_cor_mat
    
    
    print 'load regressor_vect'
    
    regressor_vect = np.loadtxt(regressor_file)
    
    
    print 'load resid data'
    
    resid_data_matrix = np.load(resid_ts_file)
    
    print "compute return_Z_cor_mat"
    
    print resid_data_matrix.shape
    print np.transpose(resid_data_matrix).shape
    
    cor_mat,sderr_cor_mat,pval_cor_mat  = return_var_cor_mat(np.transpose(resid_data_matrix),regressor_vect)
    
    print cor_mat.shape
    
    print sderr_cor_mat.shape
    
    print pval_cor_mat.shape
    
    ### 
    print "saving cor_mat as npy"
    
    cor_mat_file = os.path.abspath('cor_mat.npy')
    
    np.save(cor_mat_file,cor_mat)
    
    ### 
    print "saving sderr_cor_mat as npy"
    
    sderr_cor_mat_file = os.path.abspath('sderr_cor_mat.npy')
    
    np.save(sderr_cor_mat_file,sderr_cor_mat)
        
    ### 
    print "saving pval_cor_mat as npy"
    
    pval_cor_mat_file = os.path.abspath('pval_cor_mat.npy')
    
    np.save(pval_cor_mat_file,pval_cor_mat)
    
    return cor_mat_file,sderr_cor_mat_file,pval_cor_mat_file
    
    
def compute_conf_correlation_matrix(resid_ts_file,regressor_file,conf_interval_prob):
    
    import rpy,os
    import nibabel as nib
    import numpy as np
    
    from dmgraphanalysis.utils_cor import return_conf_cor_mat
    
    
    print 'load regressor_vect'
    
    regressor_vect = np.loadtxt(regressor_file)
    
    
    print 'load resid data'
    
    resid_data_matrix = np.load(resid_ts_file)
    
    print "compute return_Z_cor_mat"
    
    print resid_data_matrix.shape
    print np.transpose(resid_data_matrix).shape
    
    cor_mat,Z_cor_mat,conf_cor_mat = return_conf_cor_mat(np.transpose(resid_data_matrix),regressor_vect,conf_interval_prob)
    
    print cor_mat.shape
    
    print Z_cor_mat.shape
    
    print conf_cor_mat.shape
    
    ### 
    print "saving cor_mat as npy"
    
    cor_mat_file = os.path.abspath('cor_mat.npy')
    
    np.save(cor_mat_file,cor_mat)
    
    print "saving conf_cor_mat as npy"
    
    conf_cor_mat_file = os.path.abspath('conf_cor_mat.npy')
    
    np.save(conf_cor_mat_file,conf_cor_mat)
    
    print "saving Z_cor_mat as npy"
    
    Z_cor_mat_file = os.path.abspath('Z_cor_mat.npy')
    
    np.save(Z_cor_mat_file,Z_cor_mat)
    
    return cor_mat_file,Z_cor_mat_file,conf_cor_mat_file

def compute_Z_correlation_matrix(resid_ts_file,regressor_file):
    
    import rpy,os
    import nibabel as nib
    import numpy as np
    
    from dmgraphanalysis.utils_cor import return_Z_cor_mat
    
    
    print 'load regressor_vect'
    
    regressor_vect = np.loadtxt(regressor_file)
    
    
    print 'load resid data'
    
    resid_data_matrix = np.load(resid_ts_file)
    
    print "compute return_Z_cor_mat"
    
    print resid_data_matrix.shape
    print np.transpose(resid_data_matrix).shape
    
    Z_cor_mat = return_Z_cor_mat(np.transpose(resid_data_matrix),regressor_vect)
    
    print Z_cor_mat.shape
    
    ### 
    print "saving Z_cor_mat as npy"
    
    Z_cor_mat_file = os.path.abspath('Z_cor_mat.npy')
    
    np.save(Z_cor_mat_file,Z_cor_mat)
    
    return Z_cor_mat_file

    
def plot_hist_var_cor_mat(cor_mat_file,sderr_cor_mat_file,pval_cor_mat_file):

    import os
    import numpy as np
    
    from dmgraphanalysis.utils_plot import plot_hist,plot_cormat
    
    import nibabel as nib
    
    from nipype.utils.filemanip import split_filename as split_f
    
    ########### cor_mat
    
    cor_mat = np.load(cor_mat_file)
    
    #### heatmap 
    
    print 'plotting cor_mat heatmap'
    
    plot_heatmap_cor_mat_file =  os.path.abspath('heatmap_cor_mat.eps')
    
    plot_cormat(plot_heatmap_cor_mat_file,cor_mat,list_labels = [])
    
    #### histogram 
    
    print 'plotting cor_mat histogram'
    
    plot_hist_cor_mat_file = os.path.abspath('hist_cor_mat.eps')
    
    plot_hist(plot_hist_cor_mat_file,cor_mat,nb_bins = 100)
    
    
    ########## sderr_cor_mat 
    
    sderr_cor_mat = np.load(sderr_cor_mat_file)
    
    #### heatmap 
    
    print 'plotting sderr_cor_mat heatmap'
    
    plot_heatmap_sderr_cor_mat_file =  os.path.abspath('heatmap_sderr_cor_mat.eps')
    
    plot_cormat(plot_heatmap_sderr_cor_mat_file,sderr_cor_mat,list_labels = [])
    
    #### histogram 
    
    print 'plotting sderr_cor_mat histogram'
    
    plot_hist_sderr_cor_mat_file = os.path.abspath('hist_sderr_cor_mat.eps')
    
    plot_hist(plot_hist_sderr_cor_mat_file,sderr_cor_mat,nb_bins = 100)
    
    ############# pval cor_mat
    
    pval_cor_mat = np.load(pval_cor_mat_file)
    
    #### heatmap 
    
    print 'plotting pval_cor_mat heatmap'
    
    plot_heatmap_pval_cor_mat_file =  os.path.abspath('heatmap_pval_cor_mat.eps')
    
    plot_cormat(plot_heatmap_pval_cor_mat_file,pval_cor_mat,list_labels = [])
    
    #### histogram 
    
    print 'plotting pval_cor_mat histogram'
    
    plot_hist_pval_cor_mat_file = os.path.abspath('hist_pval_cor_mat.eps')
    
    plot_hist(plot_hist_pval_cor_mat_file,pval_cor_mat,nb_bins = 100)
    
    return plot_hist_cor_mat_file,plot_heatmap_cor_mat_file,plot_hist_sderr_cor_mat_file,plot_heatmap_sderr_cor_mat_file,plot_hist_pval_cor_mat_file,plot_heatmap_pval_cor_mat_file
    
def plot_hist_conf_cor_mat(cor_mat_file,Z_cor_mat_file,conf_cor_mat_file):

    import os
    import numpy as np
    
    import nibabel as nib
    
    from nipype.utils.filemanip import split_filename as split_f
    
    from dmgraphanalysis.utils_plot import plot_hist,plot_cormat
    
    ############ cor_mat
    
    cor_mat = np.load(cor_mat_file)
    
    #### heatmap 
    
    print 'plotting cor_mat heatmap'
    
    plot_heatmap_cor_mat_file =  os.path.abspath('heatmap_cor_mat.eps')
    
    plot_cormat(plot_heatmap_cor_mat_file,cor_mat,list_labels = [])
    
    #### histogram 
    
    print 'plotting cor_mat histogram'
    
    plot_hist_cor_mat_file = os.path.abspath('hist_cor_mat.eps')
    
    plot_hist(plot_hist_cor_mat_file,cor_mat,nb_bins = 100)
    
    ############ Z_cor_mat
    
    Z_cor_mat = np.load(Z_cor_mat_file)
    
    #### heatmap 
    
    print 'plotting Z_cor_mat heatmap'
    
    plot_heatmap_Z_cor_mat_file =  os.path.abspath('heatmap_Z_cor_mat.eps')
    
    plot_cormat(plot_heatmap_Z_cor_mat_file,Z_cor_mat,list_labels = [])
    
    #### histogram 
    
    print 'plotting Z_cor_mat histogram'
    
    plot_hist_Z_cor_mat_file = os.path.abspath('hist_Z_cor_mat.eps')
    
    plot_hist(plot_hist_Z_cor_mat_file,Z_cor_mat,nb_bins = 100)
    
    ############ conf_cor_mat
    
    conf_cor_mat = np.load(conf_cor_mat_file)
    
    #### heatmap 
    
    print 'plotting conf_cor_mat heatmap'
    
    plot_heatmap_conf_cor_mat_file =  os.path.abspath('heatmap_conf_cor_mat.eps')
    
    plot_cormat(plot_heatmap_conf_cor_mat_file,conf_cor_mat,list_labels = [])
    
    #### histogram 
    
    print 'plotting conf_cor_mat histogram'
    
    plot_hist_conf_cor_mat_file = os.path.abspath('hist_conf_cor_mat.eps')

    plot_hist(plot_hist_conf_cor_mat_file,conf_cor_mat,nb_bins = 100)
    
    return plot_hist_cor_mat_file,plot_heatmap_cor_mat_file,plot_hist_Z_cor_mat_file,plot_heatmap_Z_cor_mat_file,plot_hist_conf_cor_mat_file,plot_heatmap_conf_cor_mat_file
    
def plot_hist_Z_cor_mat(Z_cor_mat_file):

    import os
    import numpy as np
    
    import nibabel as nib
    
    from nipype.utils.filemanip import split_filename as split_f
    
    from dmgraphanalysis.utils_plot import plot_hist,plot_cormat
    
    ########## Z_cor_mat
    
    Z_cor_mat = np.load(Z_cor_mat_file)
    
    #### heatmap 
    
    print 'plotting Z_cor_mat heatmap'
    
    plot_heatmap_Z_cor_mat_file =  os.path.abspath('heatmap_Z_cor_mat.eps')
    
    plot_cormat(plot_heatmap_Z_cor_mat_file,Z_cor_mat,list_labels = [])
    
    #### histogram 
    
    print 'plotting Z_cor_mat histogram'
    
    plot_hist_Z_cor_mat_file = os.path.abspath('hist_Z_cor_mat.eps')
    
    plot_hist(plot_hist_Z_cor_mat_file,Z_cor_mat,nb_bins = 100)
    
    return plot_hist_Z_cor_mat_file,plot_heatmap_Z_cor_mat_file
    
##### spm_mask and all infos from contrast index
def extract_signif_contrast_mask(spm_contrast_index):

    import os
    from define_variables import nipype_analyses_path,peak_activation_mask_analysis_name,ROI_mask_prefix
    
    #### indexed_mask
    indexed_mask_rois_file =  os.path.join(nipype_analyses_path,peak_activation_mask_analysis_name, "indexed_mask-" + ROI_mask_prefix + "_spm_contrast" + str(spm_contrast_index) + ".nii")
        
    #### saving ROI coords as textfile
    ### ijk coords
    coord_rois_file =  os.path.join(nipype_analyses_path,peak_activation_mask_analysis_name, "coords-" + ROI_mask_prefix + "_spm_contrast" + str(spm_contrast_index) + ".txt")

    ### coords in MNI space
    MNI_coord_rois_file =  os.path.join(nipype_analyses_path,peak_activation_mask_analysis_name, "coords-MNI-" + ROI_mask_prefix + "_spm_contrast" + str(spm_contrast_index) + ".txt")

    #### saving ROI coords as textfile
    label_rois_file =  os.path.join(nipype_analyses_path,peak_activation_mask_analysis_name, "labels-" + ROI_mask_prefix + "_spm_contrast" + str(spm_contrast_index) + ".txt")
    #label_rois_file =  os.path.join(nipype_analyses_path,peak_activation_mask_analysis_name, "labels-" + ROI_mask_prefix + "_jane.txt")
        
    #### all info in a text file
    info_rois_file  =  os.path.join(nipype_analyses_path,peak_activation_mask_analysis_name, "info-" + ROI_mask_prefix + "_spm_contrast" + str(spm_contrast_index) + ".txt")

    return indexed_mask_rois_file,coord_rois_file,MNI_coord_rois_file,label_rois_file,info_rois_file
