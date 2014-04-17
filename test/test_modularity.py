import dmgraphanalysis as ga

def test_plot_igraph_matrix():

    mod_cor_mat_file = "/media/speeddata/Data/Nipype-Episodic/nipype_analyses/Graph_analysis_Model9_Only3events_radatools_signif_conf_correl-real_optim100/_cond_Odor_anticipation_subject_num_S02/mod_cor_mat_rada/mod_cor_mat.txt"
    
    mod_average_coords = "/media/speeddata/Data/Nipype-Episodic/nipype_analyses/Graph_analysis_Model9_Only3events_radatools_signif_conf_correl-real_optim100/_cond_Odor_anticipation_subject_num_S02/mod_ts_rada/mod_average_coords.txt"
    
    ga.modularity.plot_igraph_matrix(mod_cor_mat_file,mod_average_coords)

if __name__ =='__main__':
    
    test_plot_igraph_matrix()