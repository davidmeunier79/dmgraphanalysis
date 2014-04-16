import dmgraphanalysis as ga

def test_plot_hist_Z_cor_mat():

    Z_cor_mat_file = "/media/speeddata/Data/Nipype-Episodic/nipype_analyses/Correl_analyses-filtered-ROI_peaks-Model9_Only3events_WWW_What-cor_mat_analysis_name/_cond_Odor_anticipation_Hit-What_subject_num_S02/compute_conf_cor_mat/Z_cor_mat.npy"
    
    ga.correl_mat.plot_hist_Z_cor_mat(Z_cor_mat_file)

if __name__ =='__main__':
    
    test_plot_hist_Z_cor_mat()