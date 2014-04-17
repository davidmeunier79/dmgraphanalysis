
from dmgraphanalysis.utils_cor import *

    ################### common to radatools and louvain
def plot_signif_nbs_adj_mat(nbs_adj_mat_file,gm_mask_coords_file,gm_mask_file):

    import os
    import numpy as np
    import nibabel as nib
    import csv
    
    from dmgraphanalysis.plot_igraph import plot_igraph_3D_int_label_mat
    
    from nipype.utils.filemanip import split_filename as split_f
    
    from dmgraphanalysis.utils_plot import plot_cormat
    
    from dmgraphanalysis.utils_cor import return_img
    
    print 'load adj matrix'
    
    signif_adj_matrix = np.load(nbs_adj_mat_file)
    
    print signif_adj_matrix.shape
    
    print 'load gm mask'
    
    #with open(gm_mask_coords_file, 'Ur') as f:
        #gm_mask_coords_list = list(tuple(map(float,rec))[0:2] for rec in csv.reader(f, delimiter=' '))
    
    gm_mask_coords = np.array(np.loadtxt(gm_mask_coords_file),dtype = 'int64')
    
    print gm_mask_coords.shape
    
    #print gm_mask_coords
    
    print 'plotting igraph 3D'
    
    ######## igraph 3D
    plot_3D_nbs_adj_mat_file = os.path.abspath('plot_igraph_3D_signif_adj_mat.eps')
        
    plot_igraph_3D_int_label_mat(signif_adj_matrix,gm_mask_coords,plot_3D_nbs_adj_mat_file)
    
    ####### plot heat map
       
    #### heatmap
    print 'plotting signif_adj_matrix heatmap'
    
    plot_heatmap_nbs_adj_mat_file =  os.path.abspath('heatmap_signif_adj_mat.eps')
    
    plot_cormat(plot_heatmap_nbs_adj_mat_file,signif_adj_matrix,list_labels = [])
    
    #### significant coclassification degree in MNI space
    
    
    print 'format degree mask'
    
    signif_degree_img = return_img(gm_mask_file,gm_mask_coords,np.sum(signif_adj_matrix,axis = 0))
    
    signif_degree_img.update_header()
    
    print 'saving degree mask'
    
    signif_degree_img_file = os.path.abspath('signif_degree.nii')
    
    nib.save(signif_degree_img,signif_degree_img_file)
    
    return plot_3D_nbs_adj_mat_file,plot_heatmap_nbs_adj_mat_file,signif_degree_img_file

def plot_reorder_nbs_adj_matrix(nbs_adj_mat_file,node_order_vect_file,gm_mask_coords_file,gm_mask_file):
    
    import os
    import igraph as ig
    import numpy as np
    
    import csv
    
    print 'load adj matrix'
    
    signif_adj_matrix = np.load(nbs_adj_mat_file)
    
    print signif_adj_matrix.shape
    
    print 'load node_order_vect'
    
    node_order_vect = np.load(node_order_vect_file)
    
    print node_order_vect
    
    print "reorder signif_adj_matrix"
    
    
    signif_adj_matrix= signif_adj_matrix[node_order_vect,: ]
    
    reorder_signif_adj_matrix  = signif_adj_matrix[:, node_order_vect]
    
    
    print 'load gm mask'
    
    with open(gm_mask_coords_file, 'Ur') as f:
        gm_mask_coords_list = list(tuple(map(float,rec))[0:2] for rec in csv.reader(f, delimiter=' '))
    
    ####### plot heat map
    
    #### heatmap
    print 'plotting reorder_signif_adj_mat heatmap'
    
    plot_heatmap_reorder_signif_adj_mat_matrix_file =  os.path.abspath('heatmap_' + fname + '.eps')
    
    plot_cormat(plot_heatmap_reorder_signif_adj_mat_matrix_file,reorder_signif_adj_mat_matrix,list_labels = [])
    
    return heatmap_reorder_nbs_adj_mat_file
    

    ################################################# Z_list modular
    
def norm_coclass_to_net_list_thr(coclass_matrix_file,threshold):

    import os
    
    import numpy as np
    
    from dmgraphanalysis.utils_cor import return_int_net_list,export_List_net_from_list
    

    coclass_mat = np.load(coclass_matrix_file)
    
    #print coclass_mat.shape
    
    print "computing list by thresholding coclass_mat"
    
    coclass_mat [coclass_mat < threshold] = 0
    
    net_list = return_int_net_list(coclass_mat)
    
    #net_list = return_thr_int_net_list(coclass_mat,0)
    
    #print net_list
    
    ## Z correl_mat as list of edges
    
    print "saving net_list as list of edges"
    
    net_list_file = os.path.abspath('net_list_norm_coclass_thr_half.txt')
    
    export_List_net_from_list(net_list_file,net_list)
    
    print export_List_net_from_list
    
    return net_list_file
    
    
def coclass_to_net_list_thr(coclass_matrix_file):

    import os
    
    import numpy as np
    
    from dmgraphanalysis.utils_cor import return_thr_int_net_list,export_List_net_from_list
    

    coclass_mat = np.load(coclass_matrix_file)
    
    #print coclass_mat.shape
    
    print "computing list by thresholding coclass_mat"
    
    net_list = return_thr_int_net_list(coclass_mat,np.amax(coclass_mat)/2)
    
    #net_list = return_thr_int_net_list(coclass_mat,0)
    
    #print net_list
    
    ## Z correl_mat as list of edges
    
    print "saving net_list as list of edges"
    
    net_list_file = os.path.abspath('net_list_coclass_thr_half.txt')
    
    export_List_net_from_list(net_list_file,net_list)
    
    print export_List_net_from_list
    
    return net_list_file
    
def export_lol_mask_coclass_file(rada_lol_file,Pajek_net_file,gm_coords_file,mask_file):

    import numpy as np
    
    import nibabel as nib
    import os
    from dmgraphanalysis.utils_cor import return_mod_mask_corres,read_lol_file,read_Pajek_corres_nodes

    print 'Loading Pajek_net_file for reading node_corres'
    
    node_corres = read_Pajek_corres_nodes(Pajek_net_file)
    
    print node_corres.shape
    
    print 'Loading coords'
    
    coords = np.array(np.loadtxt(gm_coords_file),dtype = int)
    
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

    lol_mask_file = os.path.abspath("lol_mask_coclass_data.nii")

    nib.save(nib.Nifti1Image(np.array(lol_mask_data,dtype = int),mask_affine,mask_header),lol_mask_file)

    #print "returning"

    #lol_mask_file = ""
    
    return lol_mask_file
    
    
    
def plot_igraph_modules_coclass_rada(rada_lol_file,Pajek_net_file,net_list_file,gm_mask_coords_file,labels_file):
    
    import numpy as np
    import nibabel as nib
    import os
    import csv
        
    from dmgraphanalysis.utils_cor import return_mod_mask_corres,read_lol_file,read_Pajek_corres_nodes,read_List_net_file
    
    from dmgraphanalysis.plot_igraph import plot_3D_igraph_modules_net_list

    print 'loading labels'
    
    labels = [line.strip() for line in open(labels_file)]
    
    
    print 'Loading node_corres'
    
    node_corres = read_Pajek_corres_nodes(Pajek_net_file)
    
    print np.min(node_corres),np.max(node_corres)
    print node_corres.shape
    
    print 'Loading gm mask coords'
    
    gm_mask_coords = np.array(np.loadtxt(gm_mask_coords_file),dtype = 'int64')
    
    print gm_mask_coords.shape
    
    
    print "Loading community belonging file" + rada_lol_file

    community_vect = read_lol_file(rada_lol_file)
    
    #print community_vect
    print community_vect.shape
    
    print "loading net_list_net as list"
    
    net_list = read_List_net_file(net_list_file)
    
    #print net_list
    
    print 'extracting node coords'
    
    node_coords = gm_mask_coords[node_corres,:]
    
    print node_coords
    
    
    print "plotting Z_cor_mat_modules_file with igraph"
    
    coclass_all_modules_file = plot_3D_igraph_modules_net_list(community_vect,node_coords,net_list,gm_mask_coords,labels)
    
    
    return coclass_all_modules_file
    
    
    
    ################################################# reorder
    
def reorder_hclust_matrix(coclass_matrix_file,method_hie = 'ward'):
    
    #import matplotlib.pyplot as plt
    import os
    import sys
    sys.setrecursionlimit(2000)
    
    import numpy as np
    
    
    from dmgraphanalysis.utils_cor import return_hierachical_order
    
    coclass_mat = np.load(coclass_matrix_file)
    
    #print coclass_mat.shape
    
    reorder_coclass_matrix,node_order_vect = return_hierachical_order(coclass_mat,method_hie)
    
    #print node_order_vect
    
    reorder_coclass_matrix_file = os.path.abspath('reorder_coclass_matrix.npy')
    
    np.save(reorder_coclass_matrix_file,reorder_coclass_matrix)
    
    node_order_vect_file = os.path.abspath('node_order_vect.npy')
    
    np.save(node_order_vect_file,node_order_vect)
    
    node_order_vect_file = os.path.abspath('node_order_vect.txt')
    
    np.savetxt(node_order_vect_file,node_order_vect,fmt = '%d')
    
    return reorder_coclass_matrix_file,node_order_vect_file
    
def reorder_hclust_matrix_labels(coclass_matrix_file,labels_file,info_file,method_hie):
    
    #import matplotlib.pyplot as plt
    import os
    import sys
    sys.setrecursionlimit(2000)
    
    import numpy as np
    
    
    from dmgraphanalysis.utils_cor import return_hierachical_order
     
    
    
    print 'loading labels'
    labels = [line.strip() for line in open(labels_file)]
    
    np_labels = np.array(labels,dtype = 'string')
    
    #print np_labels
    
    
    print 'loading info'
    
    info = [line.strip() for line in open(info_file)]
    
    np_info = np.array(info,dtype = 'string')
    
    print np_info
    
    
    coclass_mat = np.load(coclass_matrix_file)
    
    #print coclass_mat.shape
    
    reorder_coclass_matrix,node_order_vect = return_hierachical_order(coclass_mat,method_hie)
    
    print node_order_vect
    
    #print node_order_vect
    np_reordered_labels = np_labels[node_order_vect]
    
    np_reordered_info = np_info[node_order_vect]
    
    
    #print 'reordered labels'
    #print np_reordered_labels
    
    reorder_coclass_matrix_file = os.path.abspath('reorder_coclass_matrix.npy')
    
    np.save(reorder_coclass_matrix_file,reorder_coclass_matrix)
    
    node_order_vect_file = os.path.abspath('node_order_vect.npy')
    
    np.save(node_order_vect_file,node_order_vect)
    
    node_order_vect_file = os.path.abspath('node_order_vect.txt')
    
    np.savetxt(node_order_vect_file,node_order_vect,fmt = '%d')
    
    reordered_labels_file = os.path.abspath('reordered_labels.txt')
    
    np.savetxt(reordered_labels_file,np_reordered_labels, fmt = '%s')
    
    
    reordered_info_file = os.path.abspath('reordered_info.txt')
    
    np.savetxt(reordered_info_file,np_reordered_info, fmt = '%s')
    
    
    
    return reorder_coclass_matrix_file,node_order_vect_file,reordered_labels_file,reordered_info_file
  
  
def reorder_coclass_matrix_labels_coords_with_force_order(coclass_matrix_file,labels_file,coords_file,node_order_vect_file):

    import numpy as np
    import os
    
    from nipype.utils.filemanip import split_filename as split_f
    
    from dmgraphanalysis.utils_plot import plot_ranged_cormat
    #from dmgraphanalysis.utils_plot import plot_cormat
    
    print 'loading node order'
    
    node_order_vect = np.array(np.loadtxt(node_order_vect_file),dtype = 'int')
    
    print 'loading labels'
    labels = [line.strip() for line in open(labels_file)]
    
    np_labels = np.array(labels,dtype = 'string')
    
    print 'reordering labels'
    #print node_order_vect
    np_reordered_labels = np_labels[node_order_vect]
    
    list_reordered_labels = np_reordered_labels.tolist()
    
    #print np_labels
    
    print 'loading coords'
    coords = np.loadtxt(coords_file)
    
    print coords
    
    print 'reordering coords'
    
    reordered_coords = coords[node_order_vect,:]
    
    print reordered_coords
    
    print 'loading coclass'
    coclass_mat = np.load(coclass_matrix_file)
    
    print coclass_mat.shape
    
    print 'reordering coclass'
    
    mat = coclass_mat[node_order_vect,: ]
    
    reordered_coclass_mat = mat[:, node_order_vect]
    
    ### saving
    print "saving reordered coclass"
    
    path,fname,ext = split_f(coclass_matrix_file)
    
    reordered_coclass_matrix_file =  os.path.abspath('reordered_' + fname + '.npy')
    
    np.save(reordered_coclass_matrix_file,reordered_coclass_mat)
    
    print "saving reordered labels"
    
    path,fname,ext = split_f(labels_file)
    
    reordered_labels_matrix_file =  os.path.abspath('reordered_' + fname + '.txt')
    
    np.savetxt(reordered_labels_matrix_file,list_reordered_labels, fmt = '%s')
    
    print "saving reordered coords"
    
    path,fname,ext = split_f(coords_file)
    
    reordered_coords_file =  os.path.abspath('reordered_' + fname + '.txt')
    
    np.savetxt(reordered_coords_file,reordered_coords, fmt = '%d')
    
    return reordered_coclass_matrix_file,reordered_labels_matrix_file,reordered_coords_file
    
    ################################ plot coclass matrix ############################
  
  
def plot_igraph_coclass_matrix_labels(coclass_matrix_file,gm_mask_coords_file,threshold,labels_file):

    import numpy as np
    import os
    import pylab as pl
    
    from dmgraphanalysis.plot_igraph import plot_igraph_3D_int_mat_labels
    
    from nipype.utils.filemanip import split_filename as split_f
    
    print 'loading labels'
    
    labels = [line.strip() for line in open(labels_file)]
    
    
    print 'loading coclass_matrix'
    coclass_matrix = np.load(coclass_matrix_file)
    
    path,fname,ext = split_f(coclass_matrix_file)
    
    
    print 'loading gm mask corres'
    
    gm_mask_coords = np.loadtxt(gm_mask_coords_file)
    
    print gm_mask_coords.shape
        
        
    print 'plotting igraph'
    
    coclass_matrix[coclass_matrix < threshold] = 0
    
    plot_igraph_3D_coclass_matrix_file = os.path.abspath('plot_igraph_3D_coclass_matrix.eps')
    
    plot_igraph_3D_int_mat_labels(coclass_matrix,gm_mask_coords,plot_igraph_3D_coclass_matrix_file,labels = labels)
    
    return plot_igraph_3D_coclass_matrix_file    

def plot_igraph_coclass_matrix(coclass_matrix_file,gm_mask_coords_file,threshold):

    import numpy as np
    import os
    import pylab as pl
    
    from dmgraphanalysis.plot_igraph import plot_igraph_3D_int_mat
    
    from nipype.utils.filemanip import split_filename as split_f
    
    print 'loading coclass_matrix'
    coclass_matrix = np.load(coclass_matrix_file)
    
    path,fname,ext = split_f(coclass_matrix_file)
    
    
    print 'loading gm mask corres'
    
    gm_mask_coords = np.loadtxt(gm_mask_coords_file)
    
    print gm_mask_coords.shape
        
        
    print 'plotting igraph'
    
    coclass_matrix[coclass_matrix < threshold] = 0
    
    plot_igraph_3D_coclass_matrix_file = os.path.abspath('plot_igraph_3D_coclass_matrix.eps')
    
    plot_igraph_3D_int_mat(coclass_matrix,gm_mask_coords,plot_igraph_3D_coclass_matrix_file)
    
    return plot_igraph_3D_coclass_matrix_file    

    ############################## plot coclass
    
def plot_coclass_matrix(coclass_matrix_file):

    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import pylab as pl
    
    from nipype.utils.filemanip import split_filename as split_f
    
    print 'loading sum_coclass_matrix'
    
    sum_coclass_matrix = np.load(coclass_matrix_file)
    
    path,fname,ext = split_f(coclass_matrix_file)
    
    print 'plotting heatmap'
    
    plot_heatmap_coclass_matrix_file =  os.path.abspath('heatmap_' + fname + '.eps')
    
    #fig1 = figure.Figure()
    fig1 = plt.figure()
    ax = fig1.add_subplot(1,1,1)
    im = ax.matshow(sum_coclass_matrix)
    im.set_cmap('spectral')
    fig1.colorbar(im)
    
    fig1.savefig(plot_heatmap_coclass_matrix_file)
    
    plt.close(fig1)
    #fig1.close()
    del fig1
    
    
    ############# histogram 
    
    print 'plotting distance matrix histogram'
     
    #plt.figure()
    
    #plt.figure.Figure()
    
    
    plot_hist_coclass_matrix_file = os.path.abspath('hist_coclass_matrix.eps')
    
    #fig2 = figure.Figure()
    fig2 = plt.figure()
    ax = fig2.add_subplot(1,1,1)
    y, x = np.histogram(sum_coclass_matrix, bins = 100)
    ax.plot(x[:-1],y)
    #ax.bar(x[:-1],y, width = y[1]-y[0])
    fig2.savefig(plot_hist_coclass_matrix_file)
    
    plt.close(fig2)
    #fig2.close()
    del fig2
    
    return plot_hist_coclass_matrix_file,plot_heatmap_coclass_matrix_file
    
    ################### plot coclass with labels
    
def plot_coclass_matrix_labels(coclass_matrix_file,labels_file):

    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import pylab as pl
    
    from nipype.utils.filemanip import split_filename as split_f
    
    from dmgraphanalysis.utils_plot import plot_hist,plot_cormat
    
    print 'loading labels'
    labels = [line.strip() for line in open(labels_file)]
    
    #print labels
    
    print 'loading coclass_matrix'
    
    coclass_matrix = np.load(coclass_matrix_file)
    
    path,fname,ext = split_f(coclass_matrix_file)
    
    #### heatmap
    
    print 'plotting coclass matrix heatmap'
    
    plot_heatmap_coclass_matrix_file =  os.path.abspath('heatmap_' + fname + '.eps')
    
    plot_cormat(plot_heatmap_coclass_matrix_file,coclass_matrix,labels)
    
    #### histogram 
    
    print 'plotting coclass matrix histogram'
     
    plot_hist_coclass_matrix_file = os.path.abspath('hist_coclass_matrix.eps')
    
    plot_hist(plot_hist_coclass_matrix_file,coclass_matrix)
    
    return plot_hist_coclass_matrix_file,plot_heatmap_coclass_matrix_file
    
    ########## reorder and plot coclass + labels 
    
def plot_order_coclass_matrix_labels(coclass_matrix_file,node_order_vect_file,labels_file):

    import numpy as np
    import os
    
    from nipype.utils.filemanip import split_filename as split_f
    
    from dmgraphanalysis.utils_plot import plot_cormat
    
    print 'loading labels'
    labels = [line.strip() for line in open(labels_file)]
    
    np_labels = np.array(labels,dtype = 'string')
    
    #print np_labels
    
    #print coclass_mat.shape
    
    print 'loading coclass'
    coclass_mat = np.load(coclass_matrix_file)
    
    print 'loading node order'
    
    node_order_vect = np.array(np.loadtxt(node_order_vect_file),dtype = 'int')
    
    print 'reordering labels'
    #print node_order_vect
    np_reordered_labels = np_labels[node_order_vect]
    
    list_reordered_labels = np_reordered_labels.tolist()
    
    print 'reordering coclass'
    
    mat = coclass_mat[node_order_vect,: ]
    
    reordered_coclass_mat = mat[:, node_order_vect]
    
    path,fname,ext = split_f(coclass_matrix_file)
    
    print 'plotting reorder heatmap'
    
    plot_reordered_coclass_matrix_file =  os.path.abspath('reordered_heatmap_' + fname + '.eps')
    
    plot_cormat(plot_reordered_coclass_matrix_file,reordered_coclass_mat,list_reordered_labels)
    
    return plot_reordered_coclass_matrix_file
    
def plot_coclass_matrix_labels_range(coclass_matrix_file,labels_file,list_value_range):

    import numpy as np
    import os
    
    from nipype.utils.filemanip import split_filename as split_f
    
    from dmgraphanalysis.utils_plot import plot_ranged_cormat
    #from dmgraphanalysis.utils_plot import plot_cormat
    
    print 'loading labels'
    labels = [line.strip() for line in open(labels_file)]
    
    np_labels = np.array(labels,dtype = 'string')
    
    #print np_labels
    
    #print coclass_mat.shape
    
    print 'loading coclass'
    coclass_mat = np.load(coclass_matrix_file)
    
    
    print 'plotting heatmap'
    
    path,fname,ext = split_f(coclass_matrix_file)
    
    plot_coclass_matrix_file =  os.path.abspath('heatmap_' + fname + '.eps')
    
    plot_ranged_cormat(plot_coclass_matrix_file,coclass_mat,labels,fix_full_range = list_value_range)
    
    return plot_coclass_matrix_file
    
    ########### reorder and plot coclass + labels + coords
    
#def plot_order_coclass_matrix_labels_coords(coclass_matrix_file,node_order_vect_file,labels_file,coords_file):

    #import numpy as np
    #import os
    
    #from nipype.utils.filemanip import split_filename as split_f
    
    #from dmgraphanalysis.utils_plot import plot_cormat
    
    #print 'loading labels'
    #labels = [line.strip() for line in open(labels_file)]
    
    #np_labels = np.array(labels,dtype = 'string')
    
    ##print np_labels
    
    #print 'loading coords'
    #coords = np.readtxt(coords_file)
    
    #print coords
    
    #reordered_coords = coords[node_order_vect,:]
    
    #print reordered_coords
    
    #0/0
    
    
    
    #print 'loading coclass'
    #coclass_mat = np.load(coclass_matrix_file)
    
    ##print coclass_mat.shape
    
    #print 'loading node order'
    
    #node_order_vect = np.array(np.loadtxt(node_order_vect_file),dtype = 'int')
    
    #print 'reordering labels'
    ##print node_order_vect
    #np_reordered_labels = np_labels[node_order_vect]
    
    #list_reordered_labels = np_reordered_labels.tolist()
    
    #print 'reordering coclass'
    
    #mat = coclass_mat[node_order_vect,: ]
    
    #reordered_coclass_mat = mat[:, node_order_vect]
    
    #path,fname,ext = split_f(coclass_matrix_file)
    
    #print 'plotting reorder heatmap'
    
    #plot_reordered_coclass_matrix_file =  os.path.abspath('reordered_heatmap_' + fname + '.eps')
    
    #plot_cormat(plot_reordered_coclass_matrix_file,reordered_coclass_mat,list_reordered_labels)
    
    #return plot_reordered_coclass_matrix_file
    
################################################################################### radatools only #######################################################################################


def compute_coclass_by_coclassmod(mod_files,coords_files,node_corres_files,gm_mask_coords_file,lol_coclassmod_file,node_corres_coclassmod_file):
    
    import numpy as np
    import os

    #import nibabel as nib
    
    from dmgraphanalysis.utils_cor import return_coclass_mat_list_by_module,return_coclass_mat
    from dmgraphanalysis.utils_cor import read_Pajek_corres_nodes,read_lol_file
    
    print 'loading gm mask corres'
    
    gm_mask_coords = np.loadtxt(gm_mask_coords_file)
    
    print gm_mask_coords.shape
        
        
    print 'loading coclassmod_file'
    
    community_vect = read_lol_file(lol_coclassmod_file)
    print "community_vect coclassmod:"
    print community_vect.shape
    
    node_corres_vect = read_Pajek_corres_nodes(node_corres_coclassmod_file)
    print "node_corres_vect coclassmod:"
    print node_corres_vect.shape
    
    
    print "coords_subj coclassmod:"
    corres_coords = gm_mask_coords[node_corres_vect,:]
        
    coclass_mat_list,possible_edge_mat = return_coclass_mat_list_by_module(community_vect,corres_coords,gm_mask_coords)
    
    for coclass_mat_mod in coclass_mat_list:
    
        np.fill_diagonal(coclass_mat_mod,0)
        
        
    np.fill_diagonal(possible_edge_mat,1)
            
    print len(coclass_mat_list)
    
    if len(mod_files) != len(coords_files) or len(mod_files) != len(node_corres_files):
        print "warning, length of mod_files, coords_files and node_corres_files are imcompatible {} {} {}".format(len(mod_files),len(coords_files),len(node_corres_files))
    
    
    density_coclass_by_mod_by_subj = np.zeros((len(mod_files),len(coclass_mat_list)),dtype = 'float')
        
    
    for index_file in range(len(mod_files)):
    #for index_file in range(1):
            
        print mod_files[index_file]
        
        if os.path.exists(mod_files[index_file]) and os.path.exists(node_corres_files[index_file]) and os.path.exists(coords_files[index_file]):
        
            community_vect = read_lol_file(mod_files[index_file])
            print "community_vect:"
            print community_vect.shape
            
            node_corres_vect = read_Pajek_corres_nodes(node_corres_files[index_file])
            print "node_corres_vect:"
            print node_corres_vect.shape
            
            
            coords = np.loadtxt(coords_files[index_file])
            print "coords_subj:"
            print coords.shape
            
            
            corres_coords = coords[node_corres_vect,:]
            print "corres_coords:"
            print corres_coords.shape
            
            
            coclass_mat,possible_edge_mat = return_coclass_mat(community_vect,corres_coords,gm_mask_coords)
            
            np.fill_diagonal(coclass_mat,0)
            
            np.fill_diagonal(possible_edge_mat,1)
            
            for mod_index,coclass_mat_mod in enumerate(coclass_mat_list):
            
                print coclass_mat_mod.shape
                
                correl_density_mod = np.mean(coclass_mat[coclass_mat_mod == 1],axis = 0)
                correl_density_full = np.mean(coclass_mat[possible_edge_mat == 1],axis = 0)
                
                #print correl_density_mod
                print correl_density_mod,correl_density_full
                
                density_coclass_by_mod_by_subj[index_file,mod_index] = correl_density_mod/correl_density_full
            
    print density_coclass_by_mod_by_subj
    
    print 'saving density_coclass_by_mod_by_subj'
    
    density_coclass_by_mod_by_subj_file = os.path.abspath('density_coclass_by_mod_by_subj.txt')
    
    np.savetxt(density_coclass_by_mod_by_subj_file,density_coclass_by_mod_by_subj,fmt = "%1.3f")
    
    return density_coclass_by_mod_by_subj_file
    
def prepare_nbs_stats_rada(mod_files,coords_files,node_corres_files,gm_mask_coords_file):
    
    import numpy as np
    import os

    #import nibabel as nib
    
    from dmgraphanalysis.utils_cor import return_coclass_mat,return_hierachical_order
    from dmgraphanalysis.utils_cor import read_Pajek_corres_nodes,read_lol_file
    
    print 'loading gm mask corres'
    
    gm_mask_coords = np.loadtxt(gm_mask_coords_file)
    
    print gm_mask_coords.shape
        
    #### read matrix from the first group
    #print Z_cor_mat_files
    
    sum_coclass_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = int)
    sum_possible_edge_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = int)
    
    #print sum_coclass_matrix.shape
    
            
    group_coclass_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0],len(mod_files)),dtype = float)
    
    print group_coclass_matrix.shape
    
    if len(mod_files) != len(coords_files) or len(mod_files) != len(node_corres_files):
        print "warning, length of mod_files, coords_files and node_corres_files are imcompatible {} {} {}".format(len(mod_files),len(coords_files),len(node_corres_files))
    
    for index_file in range(len(mod_files)):
    #for index_file in range(1):
            
        print mod_files[index_file]
        
        if os.path.exists(mod_files[index_file]) and os.path.exists(node_corres_files[index_file]) and os.path.exists(coords_files[index_file]):
        
            community_vect = read_lol_file(mod_files[index_file])
            print "community_vect:"
            print community_vect.shape
            
            node_corres_vect = read_Pajek_corres_nodes(node_corres_files[index_file])
            print "node_corres_vect:"
            print node_corres_vect.shape
            
            
            coords = np.loadtxt(coords_files[index_file])
            print "coords_subj:"
            print coords.shape
            
            
            corres_coords = coords[node_corres_vect,:]
            print "corres_coords:"
            print corres_coords.shape
            
            
            coclass_mat,possible_edge_mat = return_coclass_mat(community_vect,corres_coords,gm_mask_coords)
            
            np.fill_diagonal(coclass_mat,0)
            
            np.fill_diagonal(possible_edge_mat,1)
            
            sum_coclass_matrix += coclass_mat
            
            sum_possible_edge_matrix += possible_edge_mat
            
            group_coclass_matrix[:,:,index_file] = coclass_mat
            
            
        else:
            print "Warning, one or more files between " + mod_files[index_file] + ',' + node_corres_files[index_file] + ', ' + coords_files[index_file] + " do not exists"
        
        
    group_coclass_matrix_file= os.path.abspath('group_coclass_matrix.npy')
    
    np.save(group_coclass_matrix_file,group_coclass_matrix)
    
        
    print 'saving coclass matrix'
    
    sum_coclass_matrix_file = os.path.abspath('sum_coclass_matrix.npy')
    
    np.save(sum_coclass_matrix_file,sum_coclass_matrix)
    
    print 'saving possible_edge matrix'
    
    sum_possible_edge_matrix_file = os.path.abspath('sum_possible_edge_matrix.npy')
    
    np.save(sum_possible_edge_matrix_file,sum_possible_edge_matrix)
    
    
    #### save norm_coclass_matrix
    print 
    
    norm_coclass_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = int)
    
    print np.where(np.array(sum_possible_edge_matrix == 0))
        
    norm_coclass_matrix = np.divide(np.array(sum_coclass_matrix,dtype = float),np.array(sum_possible_edge_matrix,dtype = float)) * 100
    
    
    #0/0
    
    print 'saving norm coclass matrix'
    
    norm_coclass_matrix_file =  os.path.abspath('norm_coclass_matrix.npy')
    
    np.save(norm_coclass_matrix_file,norm_coclass_matrix)
    
    
    return group_coclass_matrix_file,sum_coclass_matrix_file,sum_possible_edge_matrix_file,norm_coclass_matrix_file
        
########################################################### only function to be use in run_mean_correl.py #################################################

def prepare_mean_correlation_matrices(cor_mat_files,coords_files,gm_mask_coords_file):
    
    import numpy as np
    import os

    #import nibabel as nib
    
    from dmgraphanalysis.utils_cor import return_corres_correl_mat
    #from utils_cor import read_Pajek_corres_nodes,read_lol_file
    
    print 'loading gm mask corres'
    
    gm_mask_coords = np.loadtxt(gm_mask_coords_file)
    
    print gm_mask_coords.shape
        
    #### read matrix from the first group
    #print Z_cor_mat_files
    
    sum_cor_mat_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = float)
    print sum_cor_mat_matrix.shape
    
    sum_possible_edge_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = int)
    print sum_possible_edge_matrix.shape
    
            
    group_cor_mat_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0],len(cor_mat_files)),dtype = float)
    print group_cor_mat_matrix.shape
    
    if len(cor_mat_files) != len(coords_files):
        print "warning, length of cor_mat_files, coords_files are imcompatible {} {} {}".format(len(cor_mat_files),len(coords_files))
    
    for index_file in range(len(cor_mat_files)):
        
        print cor_mat_files[index_file]
        
        if os.path.exists(cor_mat_files[index_file]) and os.path.exists(coords_files[index_file]):
        
            Z_cor_mat = np.load(cor_mat_files[index_file])
            print Z_cor_mat.shape
            
            
            coords = np.loadtxt(coords_files[index_file])
            #print coords.shape
            
            
            
            corres_cor_mat,possible_edge_mat = return_corres_correl_mat(Z_cor_mat,coords,gm_mask_coords)
            
            
            np.fill_diagonal(corres_cor_mat,0)
            
            np.fill_diagonal(possible_edge_mat,1)
            
            sum_cor_mat_matrix += corres_cor_mat
            
            sum_possible_edge_matrix += possible_edge_mat
            
            
            group_cor_mat_matrix[:,:,index_file] = corres_cor_mat
            
            
        else:
            print "Warning, one or more files between " + cor_mat_files[index_file] + ', ' + coords_files[index_file] + " do not exists"
        
        
    group_cor_mat_matrix_file= os.path.abspath('group_cor_mat_matrix.npy')
    
    np.save(group_cor_mat_matrix_file,group_cor_mat_matrix)
    
    
    print 'saving sum cor_mat matrix'
    
    sum_cor_mat_matrix_file = os.path.abspath('sum_cor_mat_matrix.npy')
    
    np.save(sum_cor_mat_matrix_file,sum_cor_mat_matrix)
    
    print 'saving sum_possible_edge matrix'
    
    sum_possible_edge_matrix_file = os.path.abspath('sum_possible_edge_matrix.npy')
    
    np.save(sum_possible_edge_matrix_file,sum_possible_edge_matrix)
    
    print 'saving avg_cor_mat_matrix'
    
    avg_cor_mat_matrix_file  = os.path.abspath('avg_cor_mat_matrix.npy')
    
    avg_cor_mat_matrix = np.zeros((gm_mask_coords.shape[0],gm_mask_coords.shape[0]),dtype = int)
    
    if (np.where(np.array(sum_possible_edge_matrix == 0)) != 0):
    
            avg_cor_mat_matrix = np.divide(np.array(sum_cor_mat_matrix,dtype = float),np.array(sum_possible_edge_matrix,dtype = float))
            
            np.save(avg_cor_mat_matrix_file,avg_cor_mat_matrix)
    
    else:
            print "!!!!!!!!!!!!!!!!!!!!!!Breaking!!!!!!!!!!!!!!!!!!!!!!!!, found 0 elements in sum_cor_mat_matrix"
            return
            
    return group_cor_mat_matrix_file,sum_cor_mat_matrix_file,sum_possible_edge_matrix_file,avg_cor_mat_matrix_file
        
        