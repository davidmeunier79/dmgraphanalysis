# -*- coding: utf-8 -*-

from dmgraphanalysis.plot_igraph import *

############################################## Compute thresholded Z_list ###########################################

def compute_signif_conf_Z_list(cor_mat_file,conf_cor_mat_file,coords_file):       
        
    import rpy,os
    import nibabel as nib
    import numpy as np
    
    from dmgraphanalysis.utils_cor import export_List_net_from_list,export_Louvain_net_from_list
    from dmgraphanalysis.utils_cor import return_signif_conf_net_list
    from dmgraphanalysis.utils_plot import plot_cormat
    
    print "loading cor_mat_file"
    
    cor_mat = np.load(cor_mat_file)
    
    print "loading conf_cor_mat_file"
    
    conf_cor_mat = np.load(conf_cor_mat_file)
    
    print 'load coords'
    
    coords = np.array(np.loadtxt(coords_file),dtype = int)
    
    print "computing net_list by thresholding conf_cor_mat based on distance and net_threshold"
    
    net_list,binary_signif_matrix = return_signif_conf_net_list(cor_mat,conf_cor_mat)
    
    print binary_signif_matrix.shape
    
    print "saving binary_signif_matrix"
    
    binary_signif_matrix_file = os.path.abspath('binary_signif_matrix.npy')
    
    np.save(binary_signif_matrix_file,binary_signif_matrix)
    
    print "plotting binary_signif_matrix"
    
    plot_binary_signif_matrix_file = os.path.abspath('binary_signif_matrix.eps')
    
    plot_cormat(plot_binary_signif_matrix_file,binary_signif_matrix,list_labels = [])
    
    ## Z correl_mat as list of edges
    
    print "saving net_list as list of edges"
    
    net_List_file = os.path.abspath('net_List_signif_conf.txt')
    
    export_List_net_from_list(net_List_file,net_list)
    
    ### Z correl_mat as Louvain format
    
    print "saving net_list as Louvain format"
    
    net_Louvain_file = os.path.abspath('net_Louvain_signif_conf.txt')
    
    export_Louvain_net_from_list(net_Louvain_file,net_list,coords)
    
    #net_List_file = ''
    #net_Louvain_file = ''
    
    return net_List_file, net_Louvain_file
    
    
    
#def compute_dist_thr_conf_list(cor_mat_file,conf_cor_mat_file,coords_file):       
        
    #import rpy,os
    #import nibabel as nib
    #import numpy as np
    
    #from define_variables import cor_conf_thr,min_dist_between_voxels
    #from dmgraphanalysis.utils_cor import export_List_net_from_list,export_Louvain_net_from_list
    
    #from dmgraphanalysis.utils_cor import compute_dist_matrix,return_dist_thr_conf_net_list
    
    #print "loading cor_mat_file"
    
    #cor_mat = np.load(cor_mat_file)
    
    #print "loading conf_cor_mat_file"
    
    #conf_cor_mat = np.load(conf_cor_mat_file)
    
    #print 'load coords'
    
    #coords = np.array(np.loadtxt(coords_file),dtype = int)
    
    #print 'compute dist matrix'
    
    #dist_mat = compute_dist_matrix(coords)
    
    #print dist_mat.shape
    
        
    #print "computing net_list by thresholding conf_cor_mat based on distance and net_threshold"
    
    #net_list = return_dist_thr_conf_net_list(cor_mat,conf_cor_mat,dist_mat,cor_conf_thr)
    
    ### Z correl_mat as list of edges
    
    #print "saving dist matrix"
    
    #dist_mat_file = os.path.abspath('dist_matrix.npy')
    
    #np.save(dist_mat_file,dist_mat)
    
    #print "saving net_list as list of edges"
    
    #net_List_file = os.path.abspath('net_List_thr_' + str(cor_conf_thr) + '_dist_' + str(min_dist_between_voxels) +'.txt')
    
    #export_List_net_from_list(net_List_file,net_list)
    
    #### Z correl_mat as Louvain format
    
    #print "saving net_list as Louvain format"
    
    #net_Louvain_file = os.path.abspath('net_Louvain_thr_' + str(cor_conf_thr)  + '_dist_' + str(min_dist_between_voxels) +'.txt')
    
    #export_Louvain_net_from_list(net_Louvain_file,net_list,coords)
    
    ##net_List_file = ''
    ##net_Louvain_file = ''
    
    #return net_List_file, net_Louvain_file, dist_mat_file
    
#def compute_dist_thr_Z_list(conf_cor_mat_file,coords_file):       
        
    #import rpy,os
    #import nibabel as nib
    #import numpy as np
    
    #from dmgraphanalysis.utils_cor import export_List_net_from_list,export_Louvain_net_from_list
    
    #from dmgraphanalysis.utils_cor import compute_dist_matrix,return_dist_thr_net_list
    
    #print "loading conf_cor_mat_file"
    
    #conf_cor_mat = np.load(conf_cor_mat_file)
    
    #print 'load coords'
    
    #coords = np.array(np.loadtxt(coords_file),dtype = int)
    
    #print 'compute dist matrix'
    
    #dist_mat = compute_dist_matrix(coords)
    
    #print dist_mat.shape
    
        
    #print "computing Z_list by thresholding conf_cor_mat based on distance and Z_threshold"
    
    #Z_list = return_dist_thr_net_list(conf_cor_mat,dist_mat,Z_thr,min_dist_between_voxels)
    
    ### Z correl_mat as list of edges
    
    #print "saving dist matrix"
    
    #dist_mat_file = os.path.abspath('dist_matrix.npy')
    
    #np.save(dist_mat_file,dist_mat)
    
    #print "saving Z_list as list of edges"
    
    #net_List_file = os.path.abspath('net_List_thr_' + str(Z_thr) + '_dist_' + str(min_dist_between_voxels) +'.txt')
    
    #export_List_net_from_list(net_List_file,Z_list)
    
    #### Z correl_mat as Louvain format
    
    #print "saving Z_list as Louvain format"
    
    #net_Louvain_file = os.path.abspath('net_Louvain_thr_' + str(Z_thr)  + '_dist_' + str(min_dist_between_voxels) +'.txt')
    
    #export_Louvain_net_from_list(net_Louvain_file,Z_list,coords)
    
    ##net_List_file = ''
    ##net_Louvain_file = ''
    
    #return net_List_file, net_Louvain_file, dist_mat_file
    
    
#def compute_maxdist_thr_Z_list(conf_cor_mat_file,coords_file):       
        
    #import rpy,os
    #import nibabel as nib
    #import numpy as np
    
    #from define_variables import Z_thr,max_dist_between_voxels
    #from dmgraphanalysis.utils_cor import export_List_net_from_list,export_Louvain_net_from_list
    
    #from dmgraphanalysis.utils_cor import return_maxdist_thr_net_list,compute_dist_matrix
    
    #print "loading conf_cor_mat_file"
    
    #conf_cor_mat = np.load(conf_cor_mat_file)
    
    #print 'load coords'
    
    #coords = np.array(np.loadtxt(coords_file),dtype = int)
    
    #print 'compute dist matrix'
    
    #dist_mat = compute_dist_matrix(coords)
    
    #print dist_mat.shape
    
    #print "computing Z_list by thresholding conf_cor_mat based on distance and Z_threshold"
    
    #Z_list = return_maxdist_thr_net_list(conf_cor_mat,dist_mat)
    
    ### Z correl_mat as list of edges
    
    #print "saving dist matrix"
    
    #dist_mat_file = os.path.abspath('dist_matrix.npy')
    
    #np.save(dist_mat_file,dist_mat)
    
    
    #print "saving Z_list as list of edges"
    
    #net_List_file = os.path.abspath('net_List_maxdist_' + str(max_dist_between_voxels) +'.txt')
    
    #export_List_net_from_list(net_List_file,Z_list)
    
    #### Z correl_mat as Louvain format
    
    #print "saving Z_list as Louvain format"
    
    #net_Louvain_file = os.path.abspath('net_Louvain_maxdist_' + str(max_dist_between_voxels) +'.txt')
    
    #export_Louvain_net_from_list(net_Louvain_file,Z_list,coords)
    
    ##net_List_file = ''
    ##net_Louvain_file = ''
    
    #return net_List_file, net_Louvain_file, dist_mat_file
    
#def compute_maxdist_Z_list(conf_cor_mat_file,coords_file):       
        
    #import rpy,os
    #import nibabel as nib
    #import numpy as np
    
    #from define_variables import Z_thr,max_dist_between_voxels
    #from dmgraphanalysis.utils_cor import export_List_net_from_list,export_Louvain_net_from_list
    
    #from dmgraphanalysis.utils_cor import return_maxdist_net_list,compute_dist_matrix
    
    #print "loading conf_cor_mat_file"
    
    #conf_cor_mat = np.load(conf_cor_mat_file)
    
    #print 'load coords'
    
    #coords = np.array(np.loadtxt(coords_file),dtype = int)
    
    #print 'compute dist matrix'
    
    #dist_mat = compute_dist_matrix(coords)
    
    #print dist_mat.shape
    
    #print "saving dist matrix"
    
    #dist_mat_file = os.path.abspath('dist_matrix.npy')
    
    #np.save(dist_mat_file,dist_mat)
    
    
    #print "computing Z_list by thresholding conf_cor_mat based on distance and Z_threshold"
    
    #Z_list = return_maxdist_net_list(conf_cor_mat,dist_mat)
    
    ### Z correl_mat as list of edges
    
    #print "saving dist matrix"
    
    #dist_mat_file = os.path.abspath('dist_matrix.npy')
    
    #np.save(dist_mat_file,dist_mat)
    
    
    #print "saving Z_list as list of edges"
    
    #net_List_file = os.path.abspath('net_List_maxdist_' + str(max_dist_between_voxels) +'.txt')
    
    #export_List_net_from_list(net_List_file,Z_list)
    
    #### Z correl_mat as Louvain format
    
    #print "saving Z_list as Louvain format"
    
    #net_Louvain_file = os.path.abspath('net_Louvain_maxdist_' + str(max_dist_between_voxels) +'.txt')
    
    #export_Louvain_net_from_list(net_Louvain_file,Z_list,coords)
    
    ##net_List_file = ''
    ##net_Louvain_file = ''
    
    #return net_List_file, net_Louvain_file, dist_mat_file
    
#def compute_dist_Z_list(conf_cor_mat_file,coords_file):       
        
    #import rpy,os
    #import nibabel as nib
    #import numpy as np
    
    #from define_variables import min_dist_between_voxels
    #from dmgraphanalysis.utils_cor import export_List_net_from_list,export_Louvain_net_from_list
    
    #from dmgraphanalysis.utils_cor import compute_dist_matrix,return_dist_net_list
    
    #print "loading conf_cor_mat_file"
    
    #conf_cor_mat = np.load(conf_cor_mat_file)
    
    #print 'load coords'
    
    #coords = np.array(np.loadtxt(coords_file),dtype = int)
    
    #print 'compute dist matrix'
    
    #dist_mat = compute_dist_matrix(coords)
    
    #print dist_mat.shape
    
    #print "computing Z_list by thresholding conf_cor_mat based on distance and Z_threshold"
    
    #Z_list = return_dist_net_list(conf_cor_mat,dist_mat)
    
    ### Z correl_mat as list of edges
    
    #print "saving dist matrix"
    
    #dist_mat_file = os.path.abspath('dist_matrix.npy')
    
    #np.save(dist_mat_file,dist_mat)
    
    
    #print "saving Z_list as list of edges"
    
    #net_List_file = os.path.abspath('net_List_dist_' + str(min_dist_between_voxels) +'.txt')
    
    #export_List_net_from_list(net_List_file,Z_list)
    
    #### Z correl_mat as Louvain format
    
    #print "saving Z_list as Louvain format"
    
    #net_Louvain_file = os.path.abspath('net_Louvain_dist_' + str(min_dist_between_voxels) +'.txt')
    
    #export_Louvain_net_from_list(net_Louvain_file,Z_list,coords)
    
    ##net_List_file = ''
    ##net_Louvain_file = ''
    
    #return net_List_file, net_Louvain_file, dist_mat_file
    
#def compute_thr_Z_list(conf_cor_mat_file,coords_file):       
        
    #import rpy,os
    #import nibabel as nib
    #import numpy as np
    
    #from define_variables import Z_thr
    #from dmgraphanalysis.utils_cor import export_List_net_from_list,export_Louvain_net_from_list
    
    #from dmgraphanalysis.utils_cor import return_thr_net_list
    
    #print "loading conf_cor_mat_file"
    
    #conf_cor_mat = np.load(conf_cor_mat_file)
    
    #print 'load coords'
    
    #coords = np.array(np.loadtxt(coords_file),dtype = int)
    
    ### compute Z_list 
    
    #print "computing Z_list by thresholding conf_cor_mat"
    
    #Z_list = return_thr_net_list(conf_cor_mat,Z_thr)
    
    ### Z correl_mat as list of edges
    
    #print "saving Z_list as list of edges"
    
    #net_List_file = os.path.abspath('Z_List_thr_' + str(Z_thr) +'.txt')
    
    #export_List_net_from_list(net_List_file,Z_list)
    
    #### Z correl_mat as Louvain format
    
    #print "saving Z_list as Louvain format"
    
    #net_Louvain_file = os.path.abspath('Z_Louvain_thr_' + str(Z_thr) +'.txt')
    
    #export_Louvain_net_from_list(net_Louvain_file,Z_list,coords)
    
    #return net_List_file, net_Louvain_file
    
    
    
################################################ Louvain method  ###########################################

def convert_list_Louvain(Louvain_list_file,louvain_bin_path):
    import os
    
    from nipype.utils.filemanip import split_filename as split_f
    
    #path, fname, ext = '','',''
    path, fname, ext = split_f(Louvain_list_file)
    
    Louvain_bin_file = os.path.abspath(fname + '.bin')
    
    Louvain_node_file = os.path.abspath(fname  + '.node.txt')
    
    Louvain_conf_file = os.path.abspath(fname  + '.conf')
    
    cmd = os.path.join(louvain_bin_path,'slicer') + ' -i ' + Louvain_list_file + ' -o ' +  Louvain_bin_file + ' -n ' + Louvain_node_file + ' -c ' + Louvain_conf_file + ' -u' 
    
    print "executing command " + cmd
    
    os.system(cmd)
    
    return Louvain_bin_file,Louvain_node_file,Louvain_conf_file
    
def community_list_Louvain(Louvain_bin_file,Louvain_conf_file,louvain_bin_path):
    
    import os
    from nipype.utils.filemanip import split_filename as split_f
    
    #path, fname, ext = '','',''
    path, fname, ext = split_f(Louvain_bin_file)
    
    Louvain_mod_file = os.path.abspath(fname + '.mod')
    
    cmd = os.path.join(louvain_bin_path,'community') + ' ' + Louvain_bin_file + ' ' +  Louvain_conf_file + ' > ' + Louvain_mod_file
    
    print "executing command " + cmd
    
    os.system(cmd)
    
    return Louvain_mod_file
    
    
def export_mod_mask_file(Louvain_mod_file,Louvain_node_file,coords_file,mask_file):

    import numpy as np
    import nibabel as nib
    import os
    
    from dmgraphanalysis.utils_cor import return_mod_mask_corres,read_Louvain_corres_nodes,read_mod_file

    print 'Loading node_corres'
    
    node_corres = read_Louvain_corres_nodes(Louvain_node_file)
    
    print node_corres
    print node_corres.shape
    
    print 'Loading coords'
    
    coords = np.array(np.loadtxt(coords_file),dtype = int)
    
    print coords.shape
    
    print 'Loading mask parameters'
    
    mask = nib.load(mask_file)
    
    data_mask_shape = mask.get_data().shape
    
    mask_header = mask.get_header().copy()
    
    mask_affine = np.copy(mask.get_affine())
    
    print "Loading community belonging file" + Louvain_mod_file

    community_vect = read_mod_file(Louvain_mod_file)
    
    #print community_vect
    print community_vect.shape
    
    print "transforming to nii file"
    mod_mask_data = return_mod_mask_corres(community_vect,node_corres,coords,data_mask_shape)

    #print mod_mask_data
    print mod_mask_data.shape


    print "saving npy file"

    mod_mask_file = os.path.abspath("mod_mask_data.npy")

    np.save(mod_mask_file ,np.array(mod_mask_data,dtype = int))
    
    print "saving nii file"

    mod_mask_file = os.path.abspath("mod_mask_data.nii")

    nib.save(nib.Nifti1Image(np.array(mod_mask_data,dtype = int),mask_affine,mask_header),mod_mask_file)

    print "returning"

    return mod_mask_file
    

################################# radatools  ########################################################

def prep_radatools(List_net_file,radatools_prep_path):

    import os
    
    from nipype.utils.filemanip import split_filename as split_f
    
    #path, fname, ext = '','',''
    path, fname, ext = split_f(List_net_file)
    
    Pajek_net_file = os.path.abspath(fname + '.net')
    
    cmd = os.path.join(radatools_prep_path,'List_To_Net.exe') + ' ' + List_net_file + ' ' +  Pajek_net_file + ' U' 
    
    print "executing command " + cmd
    
    os.system(cmd)
    
    return Pajek_net_file
    
    
    
def community_radatools(Pajek_net_file,optim_seq,radatools_comm_path):
    
    import os
    
    from nipype.utils.filemanip import split_filename as split_f
    
    #path, fname, ext = '','',''
    path, fname, ext = split_f(Pajek_net_file)
    
    rada_lol_file = os.path.abspath(fname + '.lol')
    rada_log_file = os.path.abspath(fname + '.log')
    
    #cmd = os.path.join(radatools_comm_path,'Communities_Detection.exe') + ' v WS f 1 ' + Pajek_net_file + ' ' + rada_lol_file + ' > ' + rada_log_file
    cmd = os.path.join(radatools_comm_path,'Communities_Detection.exe') + ' v ' + optim_seq + ' ' + Pajek_net_file + ' ' + rada_lol_file + ' > ' + rada_log_file
    
    
    print "executing command " + cmd
    
    os.system(cmd)
    
    return rada_lol_file,rada_log_file
    
    
def export_lol_mask_file(rada_lol_file,Pajek_net_file,coords_file,mask_file):

    import numpy as np
    
    import nibabel as nib
    import os
    from dmgraphanalysis.utils_cor import return_mod_mask_corres,read_lol_file,read_Pajek_corres_nodes

    print 'Loading Pajek_net_file for reading node_corres'
    
    node_corres = read_Pajek_corres_nodes(Pajek_net_file)
    
    print node_corres.shape
    
    print 'Loading coords'
    
    coords = np.array(np.loadtxt(coords_file),dtype = int)
    
    print coords.shape
    
    print 'Loading mask parameters'
    
    mask = nib.load(mask_file)
    
    data_mask_shape = mask.get_data().shape
    
    mask_header = mask.get_header().copy()
    
    mask_affine = np.copy(mask.get_affine())
    
    print "Loading community belonging file" + rada_lol_file

    community_vect = read_lol_file(rada_lol_file)
    
    print community_vect
    
    print "transforming to nii file"
    lol_mask_data = return_mod_mask_corres(community_vect,node_corres,coords,data_mask_shape)

    #print lol_mask_data
    print lol_mask_data.shape


    #print "saving npy file"

    #mod_mask_file = os.path.abspath("mod_mask_data.npy")

    #np.save(mod_mask_file ,np.array(mod_mask_data,dtype = int))
    
    print "saving nii file"

    lol_mask_file = os.path.abspath("lol_mask_data.nii")

    nib.save(nib.Nifti1Image(np.array(lol_mask_data,dtype = int),mask_affine,mask_header),lol_mask_file)

    #print "returning"

    #lol_mask_file = ""
    
    return lol_mask_file
    

############################# computation on modules 
    
    
def compute_mod_average_ts_louvain(ts_mat_file,coords_file,Louvain_mod_file,Louvain_node_file):

    import os
    import numpy as np

    from dmgraphanalysis.utils_cor import compute_average_ts_by_module_corres,read_mod_file,read_Louvain_corres_nodes

    print 'load coords'
    
    coords = np.loadtxt(coords_file)
    
    print 'load time series'
    
    ts_mat = np.load(ts_mat_file)
    
    print 'load Louvain community'
    
    community_vect = read_mod_file(Louvain_mod_file)
    
    print 'load node corres Louvain node file'
    
    node_corres = read_Louvain_corres_nodes(Louvain_node_file)
    
    print 'compute_average_ts_by_module'
    
    mod_average_ts,mod_average_coords = compute_average_ts_by_module_corres(ts_mat,coords,community_vect,node_corres)

    #print mod_average_ts
    #print mod_average_coords


    print mod_average_ts.shape
    print mod_average_coords.shape

    print "saving mod average time series"
    mod_average_ts_file = os.path.abspath('mod_average_ts_mat.npy')

    np.save(mod_average_ts_file,mod_average_ts)


    print "saving mod average coordinates"
    mod_average_coords_file = os.path.abspath('mod_average_coords.txt')

    np.savetxt(mod_average_coords_file,mod_average_coords)

    return mod_average_ts_file,mod_average_coords_file

def compute_mod_average_ts_rada(ts_mat_file,coords_file,rada_lol_file,Pajek_net_file):

    import os
    import numpy as np

    from dmgraphanalysis.utils_cor import compute_average_ts_by_module_corres,read_lol_file,read_Pajek_corres_nodes

    print 'load coords'
    
    coords = np.loadtxt(coords_file)
    
    print 'load time series'
    
    ts_mat = np.load(ts_mat_file)
    
    print 'load Louvain community'
    
    community_vect = read_lol_file(rada_lol_file)
    
    print 'load node corres from Pajek file'
    
    node_corres = read_Pajek_corres_nodes(Pajek_net_file)
    
    
    print 'compute_average_ts_by_module_corres'
    
    mod_average_ts,mod_average_coords = compute_average_ts_by_module_corres(ts_mat,coords,community_vect,node_corres)

    #print mod_average_ts
    #print mod_average_coords


    print mod_average_ts.shape
    print mod_average_coords.shape

    print "saving mod average time series"
    mod_average_ts_file = os.path.abspath('mod_average_ts_mat.npy')

    np.save(mod_average_ts_file,mod_average_ts)


    print "saving mod average coordinates"
    mod_average_coords_file = os.path.abspath('mod_average_coords.txt')

    np.savetxt(mod_average_coords_file,mod_average_coords)

    
    return mod_average_ts_file,mod_average_coords_file

def compute_mod_cor_mat(mod_average_ts_file,regressor_file):

    import os
    import numpy as np

    from dmgraphanalysis.utils_cor import compute_weighted_cor_mat_non_zeros

    print 'load regressor_vect'
    
    regressor_vect = np.loadtxt(regressor_file)
    
    print 'load mod_average_ts_mat'
    
    mod_average_ts = np.load(mod_average_ts_file)
    
    print 'compute_weighted_cor_mat_non_zeros'
    
    mod_cor_mat,mod_Z_cor_mat = compute_weighted_cor_mat_non_zeros(np.transpose(mod_average_ts),regressor_vect)

    print mod_cor_mat
    print mod_Z_cor_mat

    print "saving mod cor mat"
    mod_cor_mat_file = os.path.abspath('mod_cor_mat.txt')

    np.savetxt(mod_cor_mat_file,mod_cor_mat,fmt = '%2.2f')

    print "saving mod Z cor mat"
    mod_Z_cor_mat_file = os.path.abspath('mod_Z_cor_mat.txt')

    np.savetxt(mod_Z_cor_mat_file,mod_Z_cor_mat,fmt = '%2.2f')

    
    return mod_cor_mat_file,mod_Z_cor_mat_file

    
################################################################## plotting ###############################################
def plot_dist_matrix(dist_mat_file):

    import os
    import numpy as np
    
    import nibabel as nib
    
    from nipype.utils.filemanip import split_filename as split_f
    from dmgraphanalysis.utils_plot import plot_cormat,plot_hist
    
    ######### dist_mat
    
    dist_mat = np.load(dist_mat_file)
    
    #### heatmap 
    
    print 'plotting distance matrix heatmap'
    
    plot_heatmap_dist_mat_file =  os.path.abspath('heatmap_distance_matrix.eps')
    
    plot_cormat(plot_heatmap_dist_mat_file,dist_mat,list_labels = [])
    
    #### histogram 
    
    print 'plotting distance matrix histogram'
     
    plot_hist_dist_mat_file = os.path.abspath('hist_distance_matrix.eps')
    
    plot_hist(plot_hist_dist_mat_file,dist_mat,nb_bins = 100)
    
    return plot_hist_dist_mat_file,plot_heatmap_dist_mat_file
    
        
########################### plot igraph #####################################""

    
def plot_igraph_modules_conf_cor_mat_louvain(Louvain_mod_file,Louvain_node_file,coords_file,net_Louvain_file,gm_mask_coords_file):

    import numpy as np
    import nibabel as nib
    import os
    import csv
    
    from dmgraphanalysis.utils_cor import return_mod_mask_corres,read_Louvain_corres_nodes,read_mod_file,read_Louvain_net_file
    from dmgraphanalysis.plot_igraph import plot_3D_igraph_modules_Z_list

    print 'Loading node_corres'
    
    node_corres = read_Louvain_corres_nodes(Louvain_node_file)
    
    print node_corres
    print node_corres.shape
    
    print 'Loading coords'
    
    #with open(coords_file, 'Ur') as f:
        #coords_list = list(tuple(map(float,rec))[0:2] for rec in csv.reader(f, delimiter=' '))
    
    coords = np.array(np.loadtxt(coords_file),dtype = int)
    
    print coords.shape
    
    print 'Loading gm mask coords'
    
    gm_mask_coords = np.array(np.loadtxt(gm_mask_coords_file),dtype = 'int64')
    
    print gm_mask_coords.shape
    
    print "Loading community belonging file" + Louvain_mod_file

    community_vect = read_mod_file(Louvain_mod_file)
    
    #print community_vect
    print community_vect.shape
    
    print "loading net_Louvain_file as list"
    
    Z_list = read_Louvain_net_file(net_Louvain_file)
    
    #print Z_list
    
    print 'extracting node coords'
    
    node_coords = coords[node_corres,:]
    
    print node_coords.shape
    
    print "plotting conf_cor_mat_modules_file with igraph"
    
    Z_list_all_modules_file = plot_3D_igraph_modules_Z_list(community_vect,node_coords,Z_list,gm_mask_coords)
    
    #Z_list_all_modules_file = ''
    #conf_cor_mat_big_modules_file = ''
    
    return Z_list_all_modules_file
    
    
def plot_igraph_modules_conf_cor_mat_rada(rada_lol_file,Pajek_net_file,coords_file,net_List_file,gm_mask_coords_file):

    import numpy as np
    import nibabel as nib
    import os
    import csv
        
    from dmgraphanalysis.utils_cor import return_mod_mask_corres,read_lol_file,read_Pajek_corres_nodes,read_List_net_file
    
    from dmgraphanalysis.plot_igraph import plot_3D_igraph_modules_Z_list

    print 'Loading node_corres'
    
    node_corres = read_Pajek_corres_nodes(Pajek_net_file)
    
    print node_corres
    print node_corres.shape
    
    print 'Loading coords'
    
    
    #with open(coords_file, 'Ur') as f:
        #coords_list = list(tuple(map(float,rec))[0:2] for rec in csv.reader(f, delimiter=' '))
    
    coords = np.array(np.loadtxt(coords_file),dtype = 'int64')
    
    print coords.shape
    
    print 'Loading gm mask coords'
    
    gm_mask_coords = np.array(np.loadtxt(gm_mask_coords_file),dtype = 'int64')
    
    print gm_mask_coords.shape
    
    
    print "Loading community belonging file" + rada_lol_file

    community_vect = read_lol_file(rada_lol_file)
    
    #print community_vect
    print community_vect.shape
    
    print "loading net_List_net as list"
    
    Z_list = read_List_net_file(net_List_file)
    
    print Z_list
    
    print 'extracting node coords'
    
    node_coords = coords[node_corres,:]
    
    print node_coords.shape
    
    print "plotting conf_cor_mat_modules_file with igraph"
    
    Z_list_all_modules_file = plot_3D_igraph_modules_Z_list(community_vect,node_coords,Z_list,gm_mask_coords)
    
    
    return Z_list_all_modules_file
    
    
def plot_igraph_matrix(mod_cor_mat_file,mod_average_coords_file):

    import os
    #import igraph as ig
    import numpy as np
    
    from dmgraphanalysis.plot_igraph import plot_3D_igraph_weighted_signed_matrix
    
    print 'loading module (node) coordinates'
    
    #mod_average_coords = np.loadtxt(mod_average_coords_file)
    

    print 'load coords'
    
    mod_average_coords = np.loadtxt(mod_average_coords_file)
    
    
    #with open(mod_average_coords_file, 'Ur') as f:
        #mod_average_coords_list = list(tuple(map(float,rec))[0:2] for rec in csv.reader(f, delimiter=' '))
    
    #print mod_average_coords
    
    print "loading mod cor mat"
    
    mod_cor_mat = np.loadtxt(mod_cor_mat_file)
    
    i_graph_file = plot_3D_igraph_weighted_signed_matrix(mod_cor_mat,mod_average_coords)
    
    return i_graph_file
    