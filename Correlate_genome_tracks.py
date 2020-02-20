###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - to calculate and plot correlation matrix of a set of genome tracks (wig).
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import scipy
import scipy.cluster.hierarchy as sch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#Path to folder with wig files.
PWD_wig="C:\Sutor\Science\E_coli_DRIP-Seq\WIG\\"
PWD_wig_2="C:\Sutor\Science\TopoIV-Topo-Seq\WIG\\"

#Input: Continuous data (e.g. EcTopoI fold enrichment) (WIG).
#Dictionary of replicas 
#'Track name' : 'Path to wig file'
Dict_of_tracks={'DRIP-Seq CTD-/Rif- 1+'   : PWD_wig + "TopA-SPA_DRIP-Seq_1_S1_edt_forward_depth.wig",
                'DRIP-Seq CTD-/Rif- 1-'   : PWD_wig + "TopA-SPA_DRIP-Seq_1_S1_edt_reverse_depth.wig",
                'DRIP-Seq CTD-/Rif- 2+'   : PWD_wig + "TopA-SPA_DRIP-Seq_2_S3_edt_forward_depth.wig",
                'DRIP-Seq CTD-/Rif- 2-'   : PWD_wig + "TopA-SPA_DRIP-Seq_2_S3_edt_reverse_depth.wig",
                'DRIP-Seq CTD-/Rif- 3+'   : PWD_wig + "TopA-SPA_DRIP-Seq_3_S5_edt_forward_depth.wig",
                'DRIP-Seq CTD-/Rif- 3-'   : PWD_wig + "TopA-SPA_DRIP-Seq_3_S5_edt_reverse_depth.wig",
                'DRIP-Seq CTD-/Rif+ 1+'   : PWD_wig + "TopA-SPA_DRIP-Seq_1rif_S7_edt_forward_depth.wig",
                'DRIP-Seq CTD-/Rif+ 1-'   : PWD_wig + "TopA-SPA_DRIP-Seq_1rif_S7_edt_reverse_depth.wig",
                'DRIP-Seq CTD-/Rif+ 2+'   : PWD_wig + "TopA-SPA_DRIP-Seq_2rif_S9_edt_forward_depth.wig",
                'DRIP-Seq CTD-/Rif+ 2-'   : PWD_wig + "TopA-SPA_DRIP-Seq_2rif_S9_edt_reverse_depth.wig",
                'DRIP-Seq CTD-/Rif+ 3+'   : PWD_wig + "TopA-SPA_DRIP-Seq_3rif_S11_edt_forward_depth.wig",
                'DRIP-Seq CTD-/Rif+ 3-'   : PWD_wig + "TopA-SPA_DRIP-Seq_3rif_S11_edt_reverse_depth.wig",  
                'DRIP-Seq CTD+/Rif- 1+'   : PWD_wig + "DRIP_CTD_1_S13_edt_forward_depth.wig",
                'DRIP-Seq CTD+/Rif- 1-'   : PWD_wig + "DRIP_CTD_1_S13_edt_reverse_depth.wig",
                'DRIP-Seq CTD+/Rif- 2+'   : PWD_wig + "DRIP_CTD_2_S15_edt_forward_depth.wig",
                'DRIP-Seq CTD+/Rif- 2-'   : PWD_wig + "DRIP_CTD_2_S15_edt_reverse_depth.wig",
                'DRIP-Seq CTD+/Rif- 3+'   : PWD_wig + "DRIP_CTD_3_S17_edt_forward_depth.wig",
                'DRIP-Seq CTD+/Rif- 3-'   : PWD_wig + "DRIP_CTD_3_S17_edt_reverse_depth.wig", 
                'DRIP-Seq CTD+/Rif+ 1+'   : PWD_wig + "DRIP_CTD_r1_S19_edt_forward_depth.wig",
                'DRIP-Seq CTD+/Rif+ 1-'   : PWD_wig + "DRIP_CTD_r1_S19_edt_reverse_depth.wig",
                'DRIP-Seq CTD+/Rif+ 2+'   : PWD_wig + "DRIP_CTD_r2_S20_edt_forward_depth.wig",
                'DRIP-Seq CTD+/Rif+ 2-'   : PWD_wig + "DRIP_CTD_r2_S20_edt_reverse_depth.wig",
                'DRIP-Seq CTD+/Rif+ 3+'   : PWD_wig + "DRIP_CTD_r3_S22_edt_forward_depth.wig",
                'DRIP-Seq CTD+/Rif+ 3-'   : PWD_wig + "DRIP_CTD_r3_S22_edt_reverse_depth.wig",                 
                'TopoIV Topo-Seq Cfx+/IP+ 1+'   : PWD_wig_2 + "ParCplusCfxplusIP1_S130_edt_forward_depth.wig",
                'TopoIV Topo-Seq Cfx+/IP+ 1-'   : PWD_wig_2 + "ParCplusCfxplusIP1_S130_edt_reverse_depth.wig",
                'TopoIV Topo-Seq Cfx+/IP+ 2+'   : PWD_wig_2 + "ParCplusCfxplusIP2_S134_edt_forward_depth.wig",
                'TopoIV Topo-Seq Cfx+/IP+ 2-'   : PWD_wig_2 + "ParCplusCfxplusIP2_S134_edt_reverse_depth.wig",                                                                
                'TopoIV Topo-Seq Cfx-/IP+ 3+'   : PWD_wig_2 + "ParC-CfxplusIP2_S135_edt_forward_depth.wig",
                'TopoIV Topo-Seq Cfx-/IP+ 3-'   : PWD_wig_2 + "ParC-CfxplusIP2_S135_edt_reverse_depth.wig",                 
                }

#Path to folder with output files.
PWD_out="C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Strands_anti_correlation\\"


#######
##Parses WIG file.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    print(len(NE_values))
    return NE_values

#######
##Read, combine data into dataframe.
#######

def read_combine(Dict_of_tracks):
    #Contains data of all replicas in separate arrays.
    dict_of_replicas={}
    samples_names_array=[]
    i=0
    for replica_name, replica_path in Dict_of_tracks.items():
        i+=1
        print('Progress: ' + str(i) + '/' + str(len(Dict_of_tracks)))
        samples_names_array.append(replica_name)
        dict_of_replicas[replica_name]=wig_parsing(replica_path)
    
    Replicas_dataframe=pd.DataFrame(dict_of_replicas, columns=samples_names_array)
    print(Replicas_dataframe.head(10))    
    return Replicas_dataframe

#########
##Mask regions (deletions and multiplicated genes).
#########

def mask_dataframe(Replicas_dataframe):
    print('Dataframe size before masking: ' + str(Replicas_dataframe.shape))
    regions_to_mask=[[274500, 372148], [793800, 807500], [1199000, 1214000], [656600, 657860], 
                     [934135, 934635], [453535,454140], [365000,366805], [2795440,2796115], 
                     [2064250,2065200], [3582520,3583940], [4104940,4105675], 
                     [3424400, 3429656], [3690771, 3695995], [2725847, 2730935], [4212776, 4217871], [3466047, 3471144], [223771, 229004], [3597167,3602272]] #Deletions and duplicated regions in E. coli w3110 strain that correspond to DY330 strain.
    
    #Remove items falling into delelted or masked regions.
    maska=(Replicas_dataframe.index>=0)
    len_of_masked=0
    for j in range(len(regions_to_mask)):
        len_of_masked+=(regions_to_mask[j][1]-regions_to_mask[j][0])
        maska=maska & (~((Replicas_dataframe.index<regions_to_mask[j][1]) & (Replicas_dataframe.index>regions_to_mask[j][0])))
    Replicas_dataframe_masked=Replicas_dataframe[maska]  
    print('Length of masked regions: ' + str(len_of_masked))
    print('Dataframe size after masking: ' + str(Replicas_dataframe_masked.shape))
    return Replicas_dataframe_masked

#########
##Compute correlation matrix and draw heatmaps.
#########

#Plot diagonal correlation matrix.
def make_correlation_matrix_plot(df, cor_method, title, outpath_folder, file_name):
    fig=plt.figure(figsize=(10,10), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    #Create correlation matrix and heatmap.
    df_cor_matrix=df.corr(method=cor_method)
    df_cor_matrix.to_csv(outpath_folder+file_name+'.csv', sep='\t', header=True, index=True)
    color_ax=ax1.imshow(df_cor_matrix, interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1, aspect="equal")
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    #Label ticks.
    labels=list(df)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    ax1.set_ylim(sorted(ax1.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    #Create text annotation for heatmap pixels.
    #for i in range(len(labels)):
    #    for j in range(len(labels)):
    #        text = ax1.text(i, j, round(df_cor_matrix[labels[i]][labels[j]], 2), ha="center", va="center", color="black")    
    #Add colorbar.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    #axins1=inset_axes(ax1, width="5%",  height="50%",  loc='upper right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax1.transAxes, borderpad=0)    #From here: https://matplotlib.org/3.1.1/gallery/axes_grid1/demo_colorbar_with_inset_locator.html 
    fig.colorbar(color_ax, ticks=[-1.00, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00], shrink=0.5, panchor=(1.0, 0.5))
    plt.tight_layout()
    plt.savefig(outpath_folder+file_name+'.png', dpi=400, figsize=(10, 10))
    plt.show()
    plt.close()
    return df_cor_matrix


#########
##Plot correlation matrix.
#########

def correlation_matrix_plot(df_cor_matrix, cor_method, title, outpath_folder, file_name):
    fig=plt.figure(figsize=(10,10), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    #Create correlation matrix and heatmap.
    color_ax=ax1.imshow(df_cor_matrix, interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1, aspect="equal")
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    #Label ticks.
    labels=list(df_cor_matrix)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    ax1.set_ylim(sorted(ax1.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    #Create text annotation for heatmap pixels.
    #for i in range(len(labels)):
    #    for j in range(len(labels)):
    #        text = ax1.text(i, j, round(df_cor_matrix[labels[i]][labels[j]], 2), ha="center", va="center", color="black")    
    #Add colorbar.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    #axins1=inset_axes(ax1, width="5%",  height="50%",  loc='upper right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax1.transAxes, borderpad=0)    #From here: https://matplotlib.org/3.1.1/gallery/axes_grid1/demo_colorbar_with_inset_locator.html 
    fig.colorbar(color_ax, ticks=[-1.00, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00], shrink=0.5, panchor=(1.0, 0.5))
    plt.tight_layout()
    plt.savefig(outpath_folder+file_name+'.png', dpi=400, figsize=(10, 10))
    plt.show()
    plt.close()
    return


#######
#Identify clusters in a corralation matrix (hierarchy clustering).
#Code stolen from https://github.com/TheLoneNut/CorrelationMatrixClustering/blob/master/CorrelationMatrixClustering.ipynb
#######

def Clustering(clust_matrix, outpath_folder, file_name):
    X = clust_matrix.values
    d = sch.distance.pdist(X)   # vector of pairwise distances
    L = sch.linkage(d, method='complete')
    ind = sch.fcluster(L, 0.5*d.max(), 'distance')
    columns = [clust_matrix.columns.tolist()[i] for i in list((np.argsort(ind)))]
    clust_matrix = clust_matrix.reindex(columns, axis=1)
    clust_matrix = clust_matrix.reindex(columns, axis=0)
    clust_matrix.to_csv(outpath_folder+file_name+'.csv', sep='\t', header=True, index=True)
    return clust_matrix

def read_wig_correlate_plot(Dict_of_tracks, corr_type, PWD_out):
    Replicas_dataframe=read_combine(Dict_of_tracks)
    Replicas_dataframe_masked=mask_dataframe(Replicas_dataframe)
    Correlation_matrix=make_correlation_matrix_plot(Replicas_dataframe_masked, corr_type, 'Correlation of samples', PWD_out, "All_cond_DRIP_and_ChIP_strands_correlation_matrix_masked_plus_rRNA")
    Correlation_matrix_clusterized=Clustering(Correlation_matrix, PWD_out, "All_availiable_tracks_correlation_matrix_clusterized_masked_plus_rRNA")
    correlation_matrix_plot(Correlation_matrix_clusterized, corr_type, 'Correlation of samples clusterized', PWD_out, "All_cond_DRIP_and_ChIP_strands_correlation_matrix_clusterized_masked_plus_rRNA")    
    return

read_wig_correlate_plot(Dict_of_tracks, 'pearson', PWD_out)


def read_matrix_plot(PWD_out, corr_type):
    Correlation_matrix=pd.read_csv(PWD_out+'All_availiable_tracks_correlation_matrix_masked.csv', sep='\t', header=0, index_col=0)
    Datasets_order=['CsiR Aquino', 'Nac Aquino', 'NtrC Aquino', 'OmpR Aquino', 'Fur Beauchene', 'BolA Dressaire', 'NsrR Mehta', 'FNR Myers', 'Lrp Kroner',                    
                    'EcTopoI score', 'CRP Singh', 'RpoS Peano', 'HNS Kahramanoglou', 'EcTopoI CTD-/Rif+', 'EcTopoI CTD+/Rif-', 'EcTopoI CTD+/Rif+', 'EcTopoI CTD-/Rif-',
                    'Cra Kim', 'ArcA Park', 'GadE Seo', 'GadW Seo', 'OxyR Seq', 'RpoS Seo', 'SoxR Seo', 'SoxS Seo', 'GadX Seo', 'RpoD Myers', 'Dps Antipov',         
                    'RpoB Borukhov', 'RpoB Kahramanoglou', 'Fis Kahramanoglou', 'MatP Nolivos', 'MukB Nolivos', 'RNA-Seq Sutormin', 'TopoIV Sayyed',       
                    'TopoIV Sutormin', 'Gyrase Sutormin', 'Gyrase Rif Sutormin', 'GC']    
    Correlation_matrix=Correlation_matrix[Datasets_order]
    Correlation_matrix=Correlation_matrix.reindex(Datasets_order)
    correlation_matrix_plot(Correlation_matrix, corr_type, 'Correlation of samples clusterized', PWD_out, "All_availiable_tracks_correlation_matrix_masked_refined_order")    
    return

#read_matrix_plot(PWD_out, 'pearson')


"""
Dict_of_tracks_1={'ParC_1' : "F:\Sayyed_TopoIV_data\Cov_depth\ParC_1_FE.wig",
                  'ParE_1' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_1_FE.wig",
                  'ParE_2' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_2_FE.wig",
                  'ParE_G1' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_G1_FE.wig",
                  'ParE_S20min' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_S20min_FE.wig",
                  'ParE_S40min' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_S40min_FE.wig",  
                  'ParE_G2' : "F:\Sayyed_TopoIV_data\Cov_depth\ParE_G2_FE.wig", 
                  'NorflIP_ParC' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParC_1_FE.wig", 
                  'NorflIP_ParE' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParE_1_FE.wig",
                  'NorflIP_ParE' : "F:\Sayyed_TopoIV_data\Cov_depth\\NorflIP_ParE_2_FE.wig",
                  'Cfx_10mkM_2' : "F:\TopoIV_Topo-Seq\Fold_enrichment\Cfx_10mkM_2_FE.wig",
                  'S83L_Cfx_10mkM' : "F:\TopoIV_Topo-Seq\Fold_enrichment\S83L_Cfx_10mkM_FE.wig",
                  'S83L_Cfx_100mkM' : "F:\TopoIV_Topo-Seq\Fold_enrichment\S83L_Cfx_100mkM_FE.wig"}
Dict_of_tracks={'TopA_CTD-Rif-_1' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_minus_1_FE.wig",
                'TopA_CTD-Rif-_2' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_minus_2_FE.wig",
                'TopA_CTD-Rif-_3' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_minus_3_FE.wig",
                'TopA_CTD-Rif+_1' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_plus_1_FE.wig",
                'TopA_CTD-Rif+_2' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_plus_2_FE.wig",
                'TopA_CTD-Rif+_3' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_minus_Rif_plus_3_FE.wig",
                'TopA_CTD+Rif-_1' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_minus_1_FE.wig",
                'TopA_CTD+Rif-_2' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_minus_2_FE.wig",
                'TopA_CTD+Rif-_3' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_minus_3_FE.wig",
                'TopA_CTD+Rif+_2' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_plus_2_FE.wig",
                'TopA_CTD+Rif+_3' : "C:\Sutor\Science\TopoI-ChIP-Seq\Fold_enrichment\TopA_ChIP_CTD_plus_Rif_plus_3_FE.wig",}
                
Dict_of_tracks={'CTD-Rif- 1' : PWD_wig + "TopA_ChIP_CTD_minus_Rif_minus_1_FE.wig",
                'CTD-Rif- 2' : PWD_wig + "TopA_ChIP_CTD_minus_Rif_minus_2_FE.wig",
                'CTD-Rif- 3' : PWD_wig + "TopA_ChIP_CTD_minus_Rif_minus_3_FE.wig",
                'CTD-Rif+ 1' : PWD_wig + "TopA_ChIP_CTD_minus_Rif_plus_1_FE.wig",
                'CTD-Rif+ 2' : PWD_wig + "TopA_ChIP_CTD_minus_Rif_plus_2_FE.wig",
                'CTD-Rif+ 3' : PWD_wig + "TopA_ChIP_CTD_minus_Rif_plus_3_FE.wig",
                'CTD+Rif- 1' : PWD_wig + "TopA_ChIP_CTD_plus_Rif_minus_1_FE.wig",
                'CTD+Rif- 2' : PWD_wig + "TopA_ChIP_CTD_plus_Rif_minus_2_FE.wig",
                'CTD+Rif- 3' : PWD_wig + "TopA_ChIP_CTD_plus_Rif_minus_3_FE.wig",
                'CTD+Rif+ 1' : PWD_wig + "TopA_ChIP_CTD_plus_Rif_plus_1_FE.wig",
                'CTD+Rif+ 2' : PWD_wig + "TopA_ChIP_CTD_plus_Rif_plus_2_FE.wig",}   
                
                
Dict_of_tracks={'HNS Kahramanoglou'   : PWD_wig + "Kahramanoglou_HNS_IP_ME.wig",
                'CsiR Aquino'         : PWD_wig + "Aquino_CsiR_FE_av.wig",
                'Nac Aquino'          : PWD_wig + "Aquino_Nac_FE_av.wig",
                'NtrC Aquino'         : PWD_wig + "Aquino_NtrC_FE_Rep1.wig",
                'OmpR Aquino'         : PWD_wig + "Aquino_OmpR_FE_av.wig",
                'Fur Beauchene'       : PWD_wig + "Beauchene_Fur_FE_av.wig",
                'BolA Dressaire'      : PWD_wig + "Dressaire_BolA_FE.wig",
                'Cra Kim'             : PWD_wig + "Kim_Cra_FE_av.wig",
                'NsrR Mehta'          : PWD_wig + "Mehta_NsrR_FE_av.wig",
                'FNR Myers'           : PWD_wig + "Myers_FNR_anaerobic_FE_av.wig",
                'Lrp Kroner'          : PWD_wig + "Kroner_Lrp_delta_Lrp_FE_av.wig",
                'ArcA Park'           : PWD_wig + "Park_anaerobic_ArcA_FE_av.wig",
                'GadE Seo'            : PWD_wig + "Seo_GadE_FE_av.wig",
                'GadW Seo'            : PWD_wig + "Seo_GadW_FE_av.wig",
                'OxyR Seq'            : PWD_wig + "Seo_OxyR_FE_av.wig",
                'RpoS Seo'            : PWD_wig + "Seo_RpoS_FE_Rep1.wig",
                'SoxR Seo'            : PWD_wig + "Seo_SoxR_FE_av.wig",
                'SoxS Seo'            : PWD_wig + "Seo_SoxS_FE_av.wig",
                'GadX Seo'            : PWD_wig + "Seo_GadX_FE_av.wig",
                'EcTopoI score'       : "C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Motif_scanning\ChIP-Munk\Rep12_thr_0.001\EcTopoI_motif_w3110_scanned_both.wig", 
                'EcTopoI CTD-/Rif-'   : PWD_wig + "Sutormin_TopA_ChIP_CTD_minus_Rif_minus_FE_av_123.wig",
                'EcTopoI CTD-/Rif+'   : PWD_wig + "Sutormin_TopA_ChIP_CTD_minus_Rif_plus_FE_av_123.wig",
                'EcTopoI CTD+/Rif-'   : PWD_wig + "Sutormin_TopA_ChIP_CTD_plus_Rif_minus_FE_av.wig",
                'EcTopoI CTD+/Rif+'   : PWD_wig + "Sutormin_TopA_ChIP_CTD_plus_Rif_plus_FE_av.wig",                
                'CRP Singh'           : PWD_wig + "Singh_CRP_ME_FE_av.wig",
                'RpoS Peano'          : PWD_wig + "Peano_RpoS_FE_av.wig",                
                'RpoD Myers'          : PWD_wig + "Myers_RpoD_FE_av.wig",
                'Dps Antipov'         : PWD_wig + "Antipov_Dps_FE_Rep1.wig",
                'RpoB Borukhov'       : PWD_wig + "Borukhov_RpoB_Pol_Sofi_LB_FE.wig",
                'RpoB Kahramanoglou'  : PWD_wig + "Kahramanoglou_RpoB_IP_ME.wig",                
                'Fis Kahramanoglou'   : PWD_wig + "Kahramanoglou_Fis_IP_ME.wig",
                'MatP Nolivos'        : PWD_wig + "Nolivos_MatP_IP_av.wig",
                'MukB Nolivos'        : PWD_wig + "Nolivos_MukB_IP_av.wig",
                'RNA-Seq Sutormin'    : PWD_wig + "Sutormin_RNA_Seq_Exponential_av.wig",
                'TopoIV Sayyed'       : PWD_wig + "Sayyed_TopoIV_ParC_1_FE.wig",
                'TopoIV Sutormin'     : PWD_wig + "Sutormin_TopoIV_Cfx_FE_av.wig",                
                'Gyrase Sutormin'     : PWD_wig + "Sutormin_Gyrase_Cfx_10mkM_FE_av.wig",
                'Gyrase Rif Sutormin' : PWD_wig + "Sutormin_Gyrase_RifCfx_122mkM_10mkM_FE_av.wig",
                'GC'                  : "C:\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\E_coli_w3110_Mu_GC_133bp.wig",
                }                
"""