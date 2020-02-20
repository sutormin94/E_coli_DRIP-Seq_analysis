###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes wig files for IP and mock control,
#Performes normalization on the number of reads,
#Calculates fold enrichment of IP over the mock control.
#Calculates correlation between tracks.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.cluster.hierarchy as sch
from scipy import stats
from Bio import SeqIO
from Bio.SeqUtils import GC as GC_count
import matplotlib.patheffects as PathEffects
from matplotlib_venn import venn2, venn3, venn3_circles
from matplotlib import cm as cm
import collections
from collections import OrderedDict
import pandas as pd
from pandas import DataFrame

#Dict of wig files to be read with the number of reads indicated.
Wig_input_dict={'Rep1_mock' : ['F:\TopoI_ChIP-Seq\Ec_TopoI_data\Cov_depth_nodup\DSu_67_S54_ns_fm_ps_nodup.wig', 5849528],
                'Rep1_IP' : ['F:\TopoI_ChIP-Seq\Ec_TopoI_data\Cov_depth_nodup\DSu_69_S56_ns_fm_ps_nodup.wig', 7185055],
                'Rep2_mock' : ['F:\TopoI_ChIP-Seq\Ec_TopoI_data\Cov_depth_nodup\DSu_68_S55_ns_fm_ps_nodup.wig', 5998210],
                'Rep2_IP' : ['F:\TopoI_ChIP-Seq\Ec_TopoI_data\Cov_depth_nodup\DSu_70_S57_ns_fm_ps_nodup.wig', 7369847],
                'Rif_Rep1_mock' : ['F:\TopoI_ChIP-Seq\Ec_TopoI_data\Cov_depth_nodup\DSu_63_S50_ns_fm_ps_nodup.wig', 9677259],
                'Rif_Rep1_IP' : ['F:\TopoI_ChIP-Seq\Ec_TopoI_data\Cov_depth_nodup\DSu_65_S52_ns_fm_ps_nodup.wig', 5802478],
                'Rif_Rep2_mock' : ['F:\TopoI_ChIP-Seq\Ec_TopoI_data\Cov_depth_nodup\DSu_64_S51_ns_fm_ps_nodup.wig', 9000904],
                'Rif_Rep2_IP' : ['F:\TopoI_ChIP-Seq\Ec_TopoI_data\Cov_depth_nodup\DSu_66_S53_ns_fm_ps_nodup.wig', 6984927]
                }
#Dict of wig files to be read with the number of reads indicated.
Output_path="F:\TopoI_ChIP-Seq\Ec_TopoI_data\\"
#Dict of wif files to be read - fold enrichments.
Wig_FE_input_dict={'Rep1_FE' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_Rep1_IP_div_Rep1_mock_no_psc.wig',
                   'Rep2_FE' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_Rep2_IP_div_Rep2_mock_no_psc.wig',
                   'Rif_Rep1_FE' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_Rif_Rep1_IP_div_Rif_Rep1_mock_no_psc.wig',
                   'Rif_Rep2_FE' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_Rif_Rep2_IP_div_Rif_Rep2_mock_no_psc.wig'
                   }

#Dict of wif files to be read - fold enrichments.
Wig_av_FE_input_dict={'Rep12_FE' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_average_Rep1_FE_Rep2_FE.wig',
                   'Rif_Rep12_FE' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_average_Rif_Rep1_FE_Rif_Rep2_FE.wig',
                   'Ded_Rep12_Rif12_FE' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_average_FE_Ded.wig'
                   }

#Path to the input gff annotation.
path_to_gff_annot="C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\Escherichia_coli_K_12_w3110_Mu_genes_annotation.gff"
#Path to the file with regions to be omitted (e.g. deletions).
Deletions="C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\Deletions_and_dps_region_w3110_G_Mu_SGS.broadPeak"

#Averaged FE over different groups of genes.
FE_input_dict_5000_1000={'Gene' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_Gene_width_5000bp_gb_1000bp.wig',
                         'rRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_rRNA_width_5000bp_gb_1000bp.wig',
                         'tRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_tRNA_width_5000bp_gb_1000bp.wig',
                         'ncRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_ncRNA_width_5000bp_gb_1000bp.wig'
                         }

FE_input_dict_5000_5000={'Gene' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_Gene_width_5000bp_gb_5000bp.wig',
                         'rRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_rRNA_width_5000bp_gb_5000bp.wig',
                         'tRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_tRNA_width_5000bp_gb_5000bp.wig',
                         'ncRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_ncRNA_width_5000bp_gb_5000bp.wig'
                         }

FE_input_dict_15000_5000={'Gene' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_Gene_width_15000bp_gb_5000bp.wig',
                         'rRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_rRNA_width_15000bp_gb_5000bp.wig',
                         'tRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_tRNA_width_15000bp_gb_5000bp.wig',
                         'ncRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Ded_Rep12_Rif12_FE_over_ncRNA_width_15000bp_gb_5000bp.wig'
                         }

Rif_noRif_FE_input_dict_5000_1000={'Gene' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_Gene_width_5000bp_gb_1000bp.wig',
                                   'Gene Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_Gene_width_5000bp_gb_1000bp.wig',
                                   'rRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_rRNA_width_5000bp_gb_1000bp.wig',
                                   'rRNA Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_rRNA_width_5000bp_gb_1000bp.wig',
                                   'tRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_tRNA_width_5000bp_gb_1000bp.wig',
                                   'tRNA Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_tRNA_width_5000bp_gb_1000bp.wig',
                                   'ncRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_ncRNA_width_5000bp_gb_1000bp.wig',
                                   'ncRNA Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_ncRNA_width_5000bp_gb_1000bp.wig'
                                   }  
Rif_noRif_FE_input_dict_5000_5000={'Gene' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_Gene_width_5000bp_gb_5000bp.wig',
                                   'Gene Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_Gene_width_5000bp_gb_5000bp.wig',
                                   'rRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_rRNA_width_5000bp_gb_5000bp.wig',
                                   'rRNA Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_rRNA_width_5000bp_gb_5000bp.wig',
                                   'tRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_tRNA_width_5000bp_gb_5000bp.wig',
                                   'tRNA Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_tRNA_width_5000bp_gb_5000bp.wig',
                                   'ncRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_ncRNA_width_5000bp_gb_5000bp.wig',
                                   'ncRNA Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_ncRNA_width_5000bp_gb_5000bp.wig'
                                   } 
Rif_noRif_FE_input_dict_15000_5000={'Gene' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_Gene_width_15000bp_gb_5000bp.wig',
                                   'Gene Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_Gene_width_15000bp_gb_5000bp.wig',
                                   'rRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_rRNA_width_15000bp_gb_5000bp.wig',
                                   'rRNA Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_rRNA_width_15000bp_gb_5000bp.wig',
                                   'tRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_tRNA_width_15000bp_gb_5000bp.wig',
                                   'tRNA Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_tRNA_width_15000bp_gb_5000bp.wig',
                                   'ncRNA' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rep12_FE_over_ncRNA_width_15000bp_gb_5000bp.wig',
                                   'ncRNA Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_av_Rif_Rep12_FE_over_ncRNA_width_15000bp_gb_5000bp.wig'
                                   } 
Annot_dict_inpath={'All_genes' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_genes_expression.txt',
                   'HEG_370' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_high_expression_genes_370.txt',
                   'LEG_370' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_low_expression_genes_370.txt',
                   'HEG_no_tRNA_rRNA_270' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_high_expression_genes_no_tRNA_no_rRNA_270.txt',
                   'LEG_270' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_low_expression_genes_270.txt',
                   'All_operons' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_operons_expression.txt',
                   'rRNA_operons' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_16S_rRNA_operons.txt',
                   'Long_and_active_operons_27' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_Top_long_and_active_expression_operons.txt',
                   'Short_and_active_operons_27' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_short_active_operons_27.txt',
                   'HEO_186' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_high_expression_operons_186.txt',
                   'LEO_186' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_low_expression_operons_186.txt',
                   'HEO_no_tRNA_rRNA_144' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_high_expression_operons_no_tRNA_rRNA_144.txt',
                   'LEO_144': 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_low_expression_operons_144.txt'
                   }

#######
#Opens and reads BED file with deletions coordinates.
#Example:
#GenomeID\tStart\tEnd
#NC_007779.1_w3110_Mu\t274500\t372148
#######

def deletions_info(del_path):
    del_ar=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        del_ar.append([int(line[1]), int(line[2])])
    filein.close()
    return del_ar

#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
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
    return NE_values

#######
#Positin-by-position division of two ars of equal length.
#######

def arrays_div(ar1, ar2):
    FE=[]
    if len(ar1)!=len(ar2):
        print("Arrays are not of equal length! Aborted!")
        return
    for i in range (len(ar1)):
        if ar1[i]!=0 and ar2[i]!=0:
            FE.append(ar1[i]/ar2[i])
        else:
            FE.append(0)
    return FE

#######
#Positin-by-position deduction of two ars of equal length.
#######

def array_ded(ar1, ar2):
    FE_ded=[]
    if len(ar1)!=len(ar2):
        print("Arrays are not of equal length! Aborted!")
        return    
    for i in range(len(ar1)):
        FE_ded.append(ar1[i]-ar2[i])
    return FE_ded

#######
#Write .wig file.
#######

def write_wig(ar, fileout_path, name):
    fileout=open(fileout_path, 'w')
    fileout.write(f'track type=wiggle_0 name="{name}" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom=NC_007779.1_w3110_Mu start=1 step=1\n')
    for point in ar:
        fileout.write(f'{point}\n')
    fileout.close()
    return

#########
##Compute correlation matrix and draw heatmaps.
#########

#Plot diagonal correlation matrix.
def correlation_matrix(df, cor_method, title, outpath):
    fig=plt.figure(figsize=(8,8), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    cax=ax1.imshow(df.corr(method=cor_method), interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1)
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    labels=list(df)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    #Add colorbar, make sure to specify tick locations to match desired ticklabels.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    fig.colorbar(cax, ticks=[-1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00], shrink=0.7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=400, figsize=(8, 8))
    plt.show()
    plt.close()
    return


#######
#Identify clusters in a corralation matrix (hierarchy clustering).
#Code stolen from https://github.com/TheLoneNut/CorrelationMatrixClustering/blob/master/CorrelationMatrixClustering.ipynb
#######

def Clustering(df):
    X = df.corr().values
    d = sch.distance.pdist(X)   # vector of pairwise distances
    L = sch.linkage(d, method='complete')
    ind = sch.fcluster(L, 0.5*d.max(), 'distance')
    columns = [df.columns.tolist()[i] for i in list((np.argsort(ind)))]
    df = df.reindex_axis(columns, axis=1)
    return df

#######
#Plot dendrogram on the basis of correlation matrix.
#######

def correlation_dendrogram(df, cor_method, title, outpath):
    corr_inv=1-df.corr(method=cor_method) #compute correlation and inverse to distance
    corr_inv_condensed=sch.distance.squareform(corr_inv) #convert to condensed
    z=sch.linkage(corr_inv_condensed, method='average')
    dendrogram=sch.dendrogram(z, labels=corr_inv.columns, leaf_rotation=90)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=400, figsize=(8, 8))
    plt.show()  
    plt.close()
    return


#######
#Performs normalization and returns +IP/mock fold enrichment.
#######

def norm_devide(Wig_input_dict, fileout_path):
    #WIG parsing
    reads_num_dict={}
    reads_num_ar=[]
    dict_of_wigs={}
    for name, file in Wig_input_dict.items():
        reads_num_dict[name]=file[1]
        reads_num_ar.append(file[1])
        dict_of_wigs[name]=wig_parsing(file[0])
        
    #Normalization on the total number of reads. 
    #Adds pseudocount (1) to avoid zero values
    Min_total_NE=min(reads_num_ar)
    print('Min_total_NE: ' + str(Min_total_NE))
    dict_of_wigs_norm={}
    for sample, wig_data in dict_of_wigs.items():
        dict_of_wigs_norm[sample]=[1.0 * x * Min_total_NE/reads_num_dict[sample] for x in wig_data]
        print('Passed!')
        
    #Pairwise division: +IP/mock
    FE_rep1=arrays_div(dict_of_wigs_norm['Rep1_IP'], dict_of_wigs_norm['Rep1_mock'])
    FE_rep2=arrays_div(dict_of_wigs_norm['Rep2_IP'], dict_of_wigs_norm['Rep2_mock'])
    FE_Rif_rep1=arrays_div(dict_of_wigs_norm['Rif_Rep1_IP'], dict_of_wigs_norm['Rif_Rep1_mock'])
    FE_Rif_rep2=arrays_div(dict_of_wigs_norm['Rif_Rep2_IP'], dict_of_wigs_norm['Rif_Rep2_mock'])
    FE_dataframe=pd.DataFrame({'FE_rep1': FE_rep1,
                               'FE_rep2': FE_rep2,
                               'FE_Rif_rep1': FE_Rif_rep1,
                               'FE_Rif_rep2': FE_Rif_rep2
                               })
    
    #Replica-wise deduction.
    Rep1_ded_Rep1_Rif=array_ded(FE_rep1, FE_Rif_rep1)
    Rep2_ded_Rep2_Rif=array_ded(FE_rep2, FE_Rif_rep2)
    Rep1_ded_Rep2_Rif=array_ded(FE_rep1, FE_Rif_rep2)
    Rep2_ded_Rep1_Rif=array_ded(FE_rep2, FE_Rif_rep1) 
    
    Ded_FE_dataframe=pd.DataFrame({'Rep1_ded_Rep1_Rif': Rep1_ded_Rep1_Rif,
                                   'Rep2_ded_Rep2_Rif': Rep2_ded_Rep2_Rif,
                                   'Rep1_ded_Rep2_Rif': Rep1_ded_Rep2_Rif,
                                   'Rep2_ded_Rep1_Rif': Rep2_ded_Rep1_Rif
                                   })    
    #Write fold enrichment data.
    #write_wig(FE_rep1, fileout_path+'Fold_enrichment\TopoA_Rep1_IP_div_Rep1_mock_no_psc.wig', 'Rep1_IP/Rep1_mock')
    #write_wig(FE_rep2, fileout_path+'Fold_enrichment\TopoA_Rep2_IP_div_Rep2_mock_no_psc.wig', 'Rep2_IP/Rep2_mock')
    #write_wig(FE_Rif_rep1, fileout_path+'Fold_enrichment\TopoA_Rif_Rep1_IP_div_Rif_Rep1_mock_no_psc.wig', 'Rif_Rep1_IP/Rif_Rep1_mock')
    #write_wig(FE_Rif_rep2, fileout_path+'Fold_enrichment\TopoA_Rif_Rep2_IP_div_Rif_Rep2_mock_no_psc.wig', 'Rif_Rep2_IP/Rif_Rep2_mock')
    #Write deduced fold enrichment data.
    #write_wig(Rep1_ded_Rep1_Rif, fileout_path+'Fold_enrichment\TopoA_Rep1_FE_ded_Rif_Rep1_FE.wig', 'Rep1_FE-Rif_Rep1_FE')
    #write_wig(Rep2_ded_Rep2_Rif, fileout_path+'Fold_enrichment\TopoA_Rep2_FE_ded_Rif_Rep2_FE.wig', 'Rep2_FE-Rif_Rep2_FE')
    #write_wig(Rep1_ded_Rep2_Rif, fileout_path+'Fold_enrichment\TopoA_Rep1_FE_ded_Rif_Rep2_FE.wig', 'Rep1_FE-Rif_Rep2_FE')
    #write_wig(Rep2_ded_Rep1_Rif, fileout_path+'Fold_enrichment\TopoA_Rep2_FE_ded_Rif_Rep1_FE.wig', 'Rep2_FE-Rif_Rep1_FE')    
    
    return FE_dataframe, Ded_FE_dataframe

#FE_and_ded_FE=norm_devide(Wig_input_dict, Output_path)
#correlation_matrix(FE_and_ded_FE[0], 'spearman', "Correlation between samples fold enrichment of\nTopoA ChIP-Seq experiment", Output_path+"Figures\\FE_cor_spearman.png")
#correlation_matrix(FE_and_ded_FE[0], 'pearson', "Correlation between samples fold enrichment of\nTopoA ChIP-Seq experiment", Output_path+"Figures\\FE_cor_pearson.png")
#correlation_matrix(FE_and_ded_FE[0], 'kendall', "Correlation between samples fold enrichment of\nTopoA ChIP-Seq experiment", Output_path+"Figures\\FE_cor_kendall.png")
#correlation_matrix(FE_and_ded_FE[1], 'spearman', "Correlation between samples fold enrichment pw deduced of\nTopoA ChIP-Seq experiment", Output_path+"Figures\\FE_ded_cor_spearman.png")
#correlation_matrix(FE_and_ded_FE[1], 'pearson', "Correlation between samples fold enrichment pw deduced of\nTopoA ChIP-Seq experiment", Output_path+"Figures\\FE_ded_cor_pearson.png")
#correlation_matrix(FE_and_ded_FE[1], 'kendall', "Correlation between samples fold enrichment pw deduced of\nTopoA ChIP-Seq experiment", Output_path+"Figures\\FE_ded_cor_kendall.png")
#correlation_dendrogram(FE_and_ded_FE[0], 'spearman', "Dendrogram of distances between samples forl enrichment of\nTopoA ChIP-Seq experiment", Output_path+"Figures\\FE_cor_spearman_dendrogram.png")
#correlation_dendrogram(FE_and_ded_FE[0], 'pearson', "Dendrogram of distances between samples forl enrichment of\nTopoA ChIP-Seq experiment", Output_path+"Figures\\FE_cor_pearson_dendrogram.png")
#correlation_dendrogram(FE_and_ded_FE[0], 'kendall', "Dendrogram of distances between samples forl enrichment of\nTopoA ChIP-Seq experiment", Output_path+"Figures\\FE_cor_kendall_dendrogram.png")
#correlation_dendrogram(FE_and_ded_FE[1], 'spearman', "Dendrogram of distances between samples forl enrichment\npw deduced of TopoA ChIP-Seq experiment", Output_path+"Figures\\FE_cor_ded_spearman_dendrogram.png")
#correlation_dendrogram(FE_and_ded_FE[1], 'pearson', "Dendrogram of distances between samples forl enrichment\npw deduced of TopoA ChIP-Seq experiment", Output_path+"Figures\\FE_cor_ded_pearson_dendrogram.png")
#correlation_dendrogram(FE_and_ded_FE[1], 'kendall', "Dendrogram of distances between samples forl enrichment\npw deduced of TopoA ChIP-Seq experiment", Output_path+"Figures\\FE_cor_ded_kendall_dendrogram.png")


#########
##Makes histograms: 1) Rep_FE, Rif_Rep; 2) Ded Rep_FE-Rif_Rep_FE.
#########

def plot_FE_dist_GW(ar1, name1, ar2, name2, ar12, name12, pathout):
    #Plot distribution of FE values.

    fig=plt.figure(figsize=(15, 5), dpi=100)
    plot0=plt.subplot2grid((1,2),(0,0), rowspan=1, colspan=1) 
    bins1=range(int(min(ar1+ar2)), int(max(ar1+ar2)), 1)
    plot0.hist(ar1, bins1, color='#ff878b', edgecolor='black', alpha=0.8, label=f'{name1}', zorder=10)
    plot0.hist(ar2, bins1, color='#7FCE79', edgecolor='black', alpha=0.5, label=f'{name2}', zorder=9)
    plot0.set_yscale('log')
    plot0.set_xlabel('Fold enrichment', size=17)
    plot0.set_ylabel('Number of positions', size=17)
    plot0.set_title('FE distributions', size=18)  
    plot0.legend(fontsize=22)
    
    bins2=range(int(min(ar12)), int(max(ar12)), 1)
    plot1=plt.subplot2grid((1,2),(0,1), rowspan=1, colspan=1)     
    plot1.hist(ar12, bins2, color='#ffce91', edgecolor='black', alpha=0.5, label=f'{name12}')
    plot1.set_yscale('log')
    plot1.set_xlabel('delta Fold enrichment', size=17)
    plot1.set_ylabel('Number of positions', size=17)
    plot1.set_title('delta FE distribution', size=18) 
    plot1.legend(fontsize=22)
    
    plt.tight_layout()
    plt.show()
    plt.savefig(pathout, dpi=300, figsize=(10, 4))
    plt.close() 
    return

#######
#Wrapper: read input data, compute FE, deduce FE, plot FE distributions.
#######

def FE_analysis(Wig_input_dict, pathout):
    #WIG parsing
    dict_of_wigs={}
    for name, file in Wig_input_dict.items():
        dict_of_wigs[name]=wig_parsing(file)
    
    #Replica-wise deduction.
    Rep1_ded_Rep1_Rif=array_ded(dict_of_wigs['Rep1_FE'], dict_of_wigs['Rif_Rep1_FE'])
    Rep2_ded_Rep2_Rif=array_ded(dict_of_wigs['Rep2_FE'], dict_of_wigs['Rif_Rep2_FE'])
    Rep1_ded_Rep2_Rif=array_ded(dict_of_wigs['Rep1_FE'], dict_of_wigs['Rif_Rep2_FE'])
    Rep2_ded_Rep1_Rif=array_ded(dict_of_wigs['Rep2_FE'], dict_of_wigs['Rif_Rep1_FE'])
    
    dict_of_wigs['Rep1_FE-Rif_Rep1_FE']=Rep1_ded_Rep1_Rif
    dict_of_wigs['Rep2_FE-Rif_Rep2_FE']=Rep2_ded_Rep2_Rif
    dict_of_wigs['Rep1_FE-Rif_Rep2_FE']=Rep1_ded_Rep2_Rif
    dict_of_wigs['Rep2_FE-Rif_Rep1_FE']=Rep2_ded_Rep1_Rif
    
    #Plot distributions of FE and FE deduction.
    #plot_FE_dist(dict_of_wigs['Rep1_FE'], 'Rep1_FE', dict_of_wigs['Rif_Rep1_FE'], 'Rif_Rep1_FE', Rep1_ded_Rep1_Rif, 'Rep1_FE-Rif_Rep1_FE', f'{pathout}Figures\Rep1_FE_vs_Rif_Rep1_FE.png')
    #plot_FE_dist(dict_of_wigs['Rep2_FE'], 'Rep2_FE', dict_of_wigs['Rif_Rep2_FE'], 'Rif_Rep2_FE', Rep2_ded_Rep2_Rif, 'Rep2_FE-Rif_Rep2_FE', f'{pathout}Figures\Rep2_FE_vs_Rif_Rep2_FE.png')
    #plot_FE_dist(dict_of_wigs['Rep1_FE'], 'Rep1_FE', dict_of_wigs['Rif_Rep2_FE'], 'Rif_Rep2_FE', Rep1_ded_Rep2_Rif, 'Rep1_FE-Rif_Rep2_FE', f'{pathout}Figures\Rep1_FE_vs_Rif_Rep2_FE.png')
    #plot_FE_dist(dict_of_wigs['Rep2_FE'], 'Rep2_FE', dict_of_wigs['Rif_Rep1_FE'], 'Rif_Rep1_FE', Rep2_ded_Rep1_Rif, 'Rep2_FE-Rif_Rep1_FE', f'{pathout}Figures\Rep2_FE_vs_Rif_Rep1_FE.png')
    return dict_of_wigs

#Replica-wise analysis
#FE_tracks_dict=FE_analysis(Wig_FE_input_dict, Output_path)


#######
#Parsing gff file and preparing gene annotation.
#######

def parse_gff_annotation(gff_inpath, deletions_inpath):
    #Parsing deletions
    deletions=deletions_info(deletions_inpath)
    
    filein=open(gff_inpath, 'r')
    genes_annotation={'Gene': {},
                      'rRNA': {},
                      'tRNA': {},
                      'ncRNA': {}
                      }
    data_source={}
    for line in filein:
        line=line.rstrip()
        if line[0]!='#':
            line=line.split('\t')
            #What occurs in the annotation:
            if line[1] not in data_source:
                data_source[line[1]]={line[2]: 1}
            else:
                if line[2] not in data_source[line[1]]:
                    data_source[line[1]][line[2]]=1
                else:
                    data_source[line[1]][line[2]]+=1
            #Classify genes:
            #Protein coding genes.
            if line[1]=='ena' and line[2]=='gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='protein_coding':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['Gene'][gene_name]=[gene_start, gene_end, gene_strand]
            #rRNA genes.
            elif line[1]=='ena' and line[2]=='rRNA_gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='rRNA':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['rRNA'][gene_name]=[gene_start, gene_end, gene_strand] 
            #tRNA genes.
            elif line[1]=='ena' and line[2]=='tRNA_gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='tRNA':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['tRNA'][gene_name]=[gene_start, gene_end, gene_strand]
            #Other non-coding RNAs.
            elif line[1]=='Rfam' and line[2]=='ncRNA_gene':
                gene_start=int(line[3])
                gene_end=int(line[4])
                gene_strand=line[6]
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='ncRNA':
                        gene_name=annot[1].split('=')[1]
                        genes_annotation['ncRNA'][gene_name]=[gene_start, gene_end, gene_strand]
    filein.close()            
    return genes_annotation, data_source

#genome_annotation=parse_gff_annotation(path_to_gff_annot, Deletions)

#print(genome_annotation[1])
#print('\n\n')
#print(len(genome_annotation[0]['Gene']))
#print(len(genome_annotation[0]['rRNA']))
#print(len(genome_annotation[0]['tRNA']))
#print(len(genome_annotation[0]['ncRNA']))

#######
#Makes binnnary annotation: wig-like track consists of 1 (gene on +), -1 (gene on -), 0 if intergenic region.
#######

def make_binary_annotation(annot_dict):
    wig_like_annot=[]
    for i in range(4647454):
        wig_like_annot.append(0)
    for gene_name, gene_info in annot_dict.items():
        gene_start=gene_info[0]
        gene_end=gene_info[1]
        if gene_info[2]=='+':
            for i in range(gene_end-gene_start):
                wig_like_annot[gene_start+i]=1
        elif gene_info[2]=='-':
            for i in range(gene_end-gene_start):
                wig_like_annot[gene_start+i]=-1            
    return wig_like_annot


#########
##Makes histogram for FE over TUs: US, GB, DS.
#########

def plot_FE_dist_UDB(ar0, name0, ar1, name1, ar2, name2, pathout):
    #Plot distribution of FE values.
    
    mean_FE0=round(np.mean(ar0),2)
    print(f'Mean FE in {name0}={mean_FE0}')
    fig=plt.figure(figsize=(15, 3), dpi=100)
    bins0=np.arange(float(-7), float(10), 0.25)
    plot0=plt.subplot2grid((1,3),(0,0), rowspan=1, colspan=1)
    plot0.hist(ar0, bins0, color='#ff878b', edgecolor='black', alpha=0.8, label=f'{name0}')
    plot0.annotate(f'Mean FE={mean_FE0}', xy=(0.60, 0.8), xycoords='axes fraction', size=15)
    plot0.set_yscale('log')
    plot0.set_xlabel('Fold enrichment', size=17)
    plot0.set_ylabel('Number of TUs', size=17)
    plot0.set_title('FE in upstreams', size=18)  
    #plot0.legend(fontsize=22)
       
    mean_FE1=round(np.mean(ar1),2)
    print(f'Mean FE in {name1}={mean_FE1}')
    bins1=np.arange(float(-7), float(10), 0.25)
    plot1=plt.subplot2grid((1,3),(0,1), rowspan=1, colspan=1)     
    plot1.hist(ar1, bins1, color='#ffce91', edgecolor='black', alpha=0.5, label=f'{name1}')
    plot1.annotate(f'Mean FE={mean_FE1}', xy=(0.60, 0.8), xycoords='axes fraction', size=15)
    plot1.set_yscale('log')
    plot1.set_xlabel('Fold enrichment', size=17)
    plot1.set_ylabel('Number of TUs', size=17)
    plot1.set_title('FE over TUs bodies', size=18) 
    #plot1.legend(fontsize=22)
    
    mean_FE2=round(np.mean(ar2),2)
    print(f'Mean FE in {name2}={mean_FE2}')
    bins2=np.arange(float(-7), float(10), 0.25)
    plot2=plt.subplot2grid((1,3),(0,2), rowspan=1, colspan=1) 
    plot2.hist(ar2, bins2, color='#7FCE79', edgecolor='black', alpha=0.5, label=f'{name2}')
    plot2.annotate(f'Mean FE={mean_FE2}', xy=(0.60, 0.8), xycoords='axes fraction', size=15)
    plot2.set_yscale('log')
    plot2.set_xlabel('Fold enrichment', size=17)
    plot2.set_ylabel('Number of TUs', size=17)
    plot2.set_title('FE in downstreams', size=18)  
    #plot2.legend(fontsize=22)    
    
    plt.tight_layout()
    plt.show()
    plt.savefig(pathout, dpi=300, figsize=(15, 3))
    plt.close() 
    return

#######
#Scale regions (gene bodies) to equal length: make long shorter and short longer.
#######

def scale_gene_body(ar, length):
    scaled=[]
    if len(ar)>length: #array should be shrinked
        #Determines positions to be taken (other positions will be discarded).
        positions_to_take=[]
        while len(positions_to_take)!=length:
            position=np.random.random_integers(0,len(ar)-1)
            if position not in positions_to_take:
                positions_to_take.append(position)
            else:
                continue
        positions_to_take.sort()
        for pos in positions_to_take:
            scaled.append(ar[pos])
    elif len(ar)<length:
        #Determine positions to be duplicated (other positions will be discarded).
        scaled=ar
        for i in range(length-len(ar)):
            position=np.random.random_integers(0,len(scaled))
            if position==0:
                scaled=scaled[:position+1]+scaled[position:position+1]+scaled[position+1:]
            else:
                scaled=scaled[:position]+scaled[position-1:position]+scaled[position:]        
    elif len(ar)==length:
        scaled=ar
            
    return scaled

#######
#Convert dictionary to array, discard keys.
#######

def dict_to_ar(dictionary):
    ar=[]
    for k,v in dictionary.items():
        ar.append(v[0]) 
    return ar

#######
#Write .tab file with FE info for genes US, GB, and DS.
#######

def write_genes_FE(dict1, dict2, dict3, path_out):
    fileout=open(path_out, 'w')
    fileout.write('Gene_name\tStart\tEnd\tStrand\tTopoA_FE_US\tTopoA_FE_GB\tTopoA_FE_DS\n')
    for k, v in dict1.items():
        fileout.write(f'{k}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[0]}\t{dict2[k][0]}\t{dict3[k][0]}\n')
    fileout.close()
    return

#######
#Returns FE or Ded FE over the set of genes (US, GB, DS) - for each gene separately.
#######

def genes_and_FE(gene_annotation, genes_set_name, FE_track, FE_track_name, out_path, deletions_inpath, win_width, length):
    #Parsing deletions
    deletions=deletions_info(deletions_inpath) 
    
    #Calculate FE over genes.
    gene_US=np.array([0.0]*win_width)
    gene_DS=np.array([0.0]*win_width)
    gene_B=[]
    gene_US_mean_dict={}
    gene_DS_mean_dict={}
    gene_B_mean_dict={}
    for gene_name, gene_info in gene_annotation.items():
        delited=0
        for deletion in deletions:
            if deletion[1]>=gene_info[0]>=deletion[0] or deletion[1]>=gene_info[1]>=deletion[0]:
                delited=1
        if delited==0:
            start=gene_info[0]
            end=gene_info[1]
            glen=len(FE_track)
            if gene_info[2]=='+':
                if start<win_width:
                    gene_US+=np.array(FE_track[glen-(win_width-start):] + FE_track[:start])
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[glen-(win_width-start):] + FE_track[:start]), start, end, gene_info[2]]
                else:
                    gene_US+=np.array(FE_track[start-win_width:start])
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[start-win_width:start]), start, end, gene_info[2]]
                if end+win_width>glen:
                    gene_DS+=np.array(FE_track[end:] + FE_track[:end+win_width-glen])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[end:] + FE_track[:end+win_width-glen]), start, end, gene_info[2]]
                else:
                    gene_DS+=np.array(FE_track[end:end+win_width])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[end:end+win_width]), start, end, gene_info[2]]
                gene_B.append(FE_track[start:end])
                gene_B_mean_dict[gene_name]=[np.mean(FE_track[start:end]), start, end, gene_info[2]]
            elif gene_info[2]=='-':
                if start<win_width:
                    gene_DS+=np.array(FE_track[:start][::-1] + FE_track[glen-(win_width-start):][::-1])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[:start][::-1] + FE_track[glen-(win_width-start):][::-1]), start, end, gene_info[2]]
                else:
                    gene_DS+=np.array(FE_track[start-win_width:start][::-1])
                    gene_DS_mean_dict[gene_name]=[np.mean(FE_track[start-win_width:start][::-1]), start, end, gene_info[2]]
                if end+win_width>glen:
                    gene_US+=np.array(FE_track[:end+win_width-glen][::-1] + FE_track[end:][::-1])  
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[:end+win_width-glen][::-1] + FE_track[end:][::-1]), start, end, gene_info[2]]
                else:
                    gene_US+=np.array(FE_track[end:end+win_width][::-1]) 
                    gene_US_mean_dict[gene_name]=[np.mean(FE_track[end:end+win_width][::-1]), start, end, gene_info[2]]
                gene_B.append(FE_track[start:end][::-1])
                gene_B_mean_dict[gene_name]=[np.mean(FE_track[start:end][::-1]), start, end, gene_info[2]]
    Num_genes=len(gene_annotation)
    gene_US=gene_US/Num_genes
    gene_DS=gene_DS/Num_genes
    print(f'FE over TUs computed...')
    
    #Scale GB length.
    print(f'GB scaling in progress, it takes some time...')
    gene_body=np.array([0.0]*length)
    for gene in gene_B:
        gene_scaled=np.array(scale_gene_body(gene, length))
        gene_body+=gene_scaled
    gene_B=gene_body/Num_genes
    
    #Write wig-like file with FE over US, GB, DS.
    print(f'Writing FE over TU, GB, DS...')
    write_wig(np.concatenate((gene_US, gene_B, gene_DS), axis=None), f'{Output_path}Fold_enrichment\\No_dps_TopoA_av_{FE_track_name}_over_{genes_set_name}_width_{win_width}bp_gb_{length}bp.wig', f'{win_width}_{length}')
    
    #Plot FE over US, GB, DS. 
    print(f'Plotting FE over TU, GB, DS...')
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    positions_DS=np.arange(1+length,win_width+1+length,1)
    positions_US=np.arange(-1-win_width,-1,1)  
    positions_Body=np.arange(0,length,1)
    plot1.plot(positions_US, gene_US, linestyle='-', color='#D17E7E', linewidth=1, label='Rep12') 
    plot1.plot(positions_DS, gene_DS, linestyle='-', color='#D17E7E', linewidth=1)
    plot1.plot(positions_Body, gene_B, linestyle='-', color='#D17E7E', linewidth=2)
    #plot1.set_ylim(-0.6, 0.5) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)       
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoA fold enrichment', size=20)
    plot1.set_title(f'{FE_track_name} over the {genes_set_name}s', size=20)     
    plt.savefig(f'{out_path}Figures\\No_dps_{FE_track_name}_over_{genes_set_name}s_{win_width}bp.png', dpi=400, figsize=(10, 6))   
    #plt.show()
    plt.close()      
    
    #Make ar from dict.
    gene_US_mean=dict_to_ar(gene_US_mean_dict)
    gene_DS_mean=dict_to_ar(gene_DS_mean_dict)
    gene_B_mean=dict_to_ar(gene_B_mean_dict)
    
    #Write table contains FE for US, GB, DS of TUs in a set.
    print(f'Writing FE for TUs\' TU, GB, DS...')
    write_genes_FE(gene_US_mean_dict, gene_B_mean_dict, gene_DS_mean_dict, f'{out_path}FE_of_genes\\No_dps_FE_{FE_track_name}_over_{genes_set_name}s_{win_width}bp.txt')
    
    #Plot distribution of mean TUs' FEs.
    print(f'Plotting FE distribution over TU, GB, DS...')
    plot_FE_dist_UDB(gene_US_mean, 'Upstream', gene_B_mean, 'TU body', gene_DS_mean, 'Downstream', f'{out_path}Figures\\No_dps_FE_distribution_{FE_track_name}_over_{genes_set_name}s_{win_width}bp.png')
    print(len(gene_US_mean), len(gene_DS_mean), len(gene_B_mean))
    return gene_US, gene_DS, gene_B, gene_US_mean, gene_DS_mean, gene_B_mean, gene_US_mean_dict, gene_DS_mean_dict, gene_B_mean_dict

#######
#Plots FE or Ded FE over US and DS.
#######

def plot_FE_over_genes(Rep1, Rep2, Rif_Rep1, Rif_Rep2, genes_set_name, out_path):
    win_width=5000
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    positions_DS=np.arange(1,win_width+1,1)
    positions_US=np.arange(-1-win_width,-1,1)   
    #Rep1
    plot1.plot(positions_US, Rep1[0], linestyle='-', color='#D17E7E', linewidth=1, label='Rep1-Rif_Rep1') 
    plot1.plot(positions_DS, Rep1[1], linestyle='-', color='#D17E7E', linewidth=1) 
    #Rep2
    plot1.plot(positions_US, Rep2[0], linestyle='-', color='#65A865', linewidth=1, label='Rep2-Rif_Rep2') 
    plot1.plot(positions_DS, Rep2[1], linestyle='-', color='#65A865', linewidth=1) 
    #Rif Rep1
    plot1.plot(positions_US, Rif_Rep1[0], linestyle='--', color='#D1A47E', linewidth=1, label='Rep1-Rif_Rep2') 
    plot1.plot(positions_DS, Rif_Rep1[1], linestyle='--', color='#D1A47E', linewidth=1) 
    #Rif Rep2
    plot1.plot(positions_US, Rif_Rep2[0], linestyle='--', color='#4B7E7E', linewidth=1, label='Rep2_Rif_Rep1') 
    plot1.plot(positions_DS, Rif_Rep2[1], linestyle='--', color='#4B7E7E', linewidth=1) 
    
    plot1.set_ylim(-0.6, 0.5) #(0.75, 1.6) for FE; (0.75, 1.6) for Ded FE
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoA fold enrichment', size=20)
    plot1.set_title(f'FE over the {genes_set_name}s', size=20)     
    plt.savefig(f'{out_path}Figures\\Ded_FE_over_{genes_set_name}s_{win_width}bp.png', dpi=400, figsize=(10, 6))   
    #plt.show()
    plt.close()
        
    return
  

#FE_Gene_Rep1=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', FE_tracks_dict['Rep1_FE'], 'Rep1_FE', Output_path)
#FE_Gene_Rep2=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', FE_tracks_dict['Rep2_FE'], 'Rep2_FE', Output_path)
#FE_Gene_Rif_Rep1=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', FE_tracks_dict['Rif_Rep1_FE'], 'Rif_Rep1_FE', Output_path)
#FE_Gene_Rif_Rep2=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', FE_tracks_dict['Rif_Rep2_FE'], 'Rif_Rep2_FE', Output_path)
#plot_FE_over_genes(FE_Gene_Rep1, FE_Gene_Rep2, FE_Gene_Rif_Rep1, FE_Gene_Rif_Rep2, 'Gene', Output_path)

#FE_Gene_Rep1_ded_Rif_Rep1=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', FE_tracks_dict['Rep1_FE-Rif_Rep1_FE'], 'Rep1_FE-Rif_Rep1_FE', Output_path)
#FE_Gene_Rep2_ded_Rif_Rep2=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', FE_tracks_dict['Rep2_FE-Rif_Rep2_FE'], 'Rep2_FE-Rif_Rep2_FE', Output_path)
#FE_Gene_Rep1_ded_Rif_Rep2=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', FE_tracks_dict['Rep1_FE-Rif_Rep2_FE'], 'Rep1_FE-Rif_Rep2_FE', Output_path)
#FE_Gene_Rep2_ded_Rif_Rep1=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', FE_tracks_dict['Rep2_FE-Rif_Rep1_FE'], 'Rep2_FE-Rif_Rep1_FE', Output_path)
#plot_FE_over_genes(FE_Gene_Rep1_ded_Rif_Rep1, FE_Gene_Rep2_ded_Rif_Rep2, FE_Gene_Rep1_ded_Rif_Rep2, FE_Gene_Rep2_ded_Rif_Rep1, 'Gene', Output_path)

#FE_rRNA_Rep1=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', FE_tracks_dict['Rep1_FE'], 'Rep1_FE', Output_path)
#FE_rRNA_Rep2=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', FE_tracks_dict['Rep2_FE'], 'Rep2_FE', Output_path)
#FE_rRNA_Rif_Rep1=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', FE_tracks_dict['Rif_Rep1_FE'], 'Rif_Rep1_FE', Output_path)
#FE_rRNA_Rif_Rep2=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', FE_tracks_dict['Rif_Rep2_FE'], 'Rif_Rep2_FE', Output_path)
#plot_FE_over_genes(FE_rRNA_Rep1, FE_rRNA_Rep2, FE_rRNA_Rif_Rep1, FE_rRNA_Rif_Rep2, 'rRNA', Output_path)

#FE_rRNA_Rep1_ded_Rif_Rep1=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', FE_tracks_dict['Rep1_FE-Rif_Rep1_FE'], 'Rep1_FE-Rif_Rep1_FE', Output_path)
#FE_rRNA_Rep2_ded_Rif_Rep2=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', FE_tracks_dict['Rep2_FE-Rif_Rep2_FE'], 'Rep2_FE-Rif_Rep2_FE', Output_path)
#FE_rRNA_Rep1_ded_Rif_Rep2=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', FE_tracks_dict['Rep1_FE-Rif_Rep2_FE'], 'Rep1_FE-Rif_Rep2_FE', Output_path)
#FE_rRNA_Rep2_ded_Rif_Rep1=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', FE_tracks_dict['Rep2_FE-Rif_Rep1_FE'], 'Rep2_FE-Rif_Rep1_FE', Output_path)
#plot_FE_over_genes(FE_rRNA_Rep1_ded_Rif_Rep1, FE_rRNA_Rep2_ded_Rif_Rep2, FE_rRNA_Rep1_ded_Rif_Rep2, FE_rRNA_Rep2_ded_Rif_Rep1, 'rRNA', Output_path)


#######
#Average two arrays position-by-position (for 2 biological replicas).
#######

def average_2_replicas(ar1, ar2):
    av_ar=[]
    for i in range(len(ar1)):
        av_ar.append((ar1[i]+ar2[i])/2)
    return av_ar

#######
#Average four arrays position-by-position (for 4 possible Ded FE combinations).
#######

def average_4_replicas(ar1, ar2, ar3, ar4):
    av_ar=[]
    for i in range(len(ar1)):
        av_ar.append((ar1[i]+ar2[i]+ar3[i]+ar4[i])/4)
    return av_ar

#Rep12=average_2_replicas(FE_tracks_dict['Rep1_FE'], FE_tracks_dict['Rep2_FE'])
#Rif_Rep12=average_2_replicas(FE_tracks_dict['Rif_Rep1_FE'], FE_tracks_dict['Rif_Rep2_FE'])
#Rep12_ded_Rif_Rep12=average_4_replicas(FE_tracks_dict['Rep1_FE-Rif_Rep1_FE'], FE_tracks_dict['Rep2_FE-Rif_Rep2_FE'], FE_tracks_dict['Rep1_FE-Rif_Rep2_FE'], FE_tracks_dict['Rep2_FE-Rif_Rep1_FE'])
#write_wig(Rep12, Output_path+'Fold_enrichment\TopoA_average_Rep1_FE_Rep2_FE.wig', 'av_Rep1_FE_Rep2_FE') 
#write_wig(Rif_Rep12, Output_path+'Fold_enrichment\TopoA_average_Rif_Rep1_FE_Rif_Rep2_FE.wig', 'av_Rif_Rep1_FE_Rif_Rep2_FE')
#write_wig(Rep12_ded_Rif_Rep12, Output_path+'Fold_enrichment\TopoA_average_FE_Ded.wig', 'av_no_Rif-Rif')

#######
#Reads dict of wig-files with averaged tracks.
#######

#FE averaged WIG parsing
#dict_of_wigs={}
#for name, file in Wig_av_FE_input_dict.items():
#    dict_of_wigs[name]=wig_parsing(file)  
    
#######
#Plot FE and Ded FE over US, GB, and DS.
####### 
    
def plot_averaged_FE_Dev_FE(Rep12, Rif_Rep12, Ded12, genes_set_name, output_path, win_width, length):
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    positions_DS=np.arange(1+length,win_width+1+length,1)
    positions_US=np.arange(-1-win_width,-1,1)  
    positions_Body=np.arange(0,length,1)
    #Rep12
    plot1.plot(positions_US, Rep12[0], linestyle='-', color='#D17E7E', linewidth=1, label='Rep12') 
    plot1.plot(positions_DS, Rep12[1], linestyle='-', color='#D17E7E', linewidth=1)
    plot1.plot(positions_Body, Rep12[2], linestyle='-', color='#D17E7E', linewidth=2)
    #Rif_Rep12
    plot1.plot(positions_US, Rif_Rep12[0], linestyle='-', color='#65A865', linewidth=1, label='Rif_Rep12') 
    plot1.plot(positions_DS, Rif_Rep12[1], linestyle='-', color='#65A865', linewidth=1) 
    plot1.plot(positions_Body, Rif_Rep12[2], linestyle='-', color='#65A865', linewidth=2) 
    
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoA fold enrichment', size=20)
    plot1.set_title(f'FE over the {genes_set_name}s', size=20)     
    plt.savefig(f'{output_path}Figures\\Averaged_FE_over_{genes_set_name}s_{win_width}bp_nd_with_body_{length}.svg', dpi=400, format='svg', figsize=(10, 6))   
    #plt.show()
    plt.close()    
    
    #Ded12
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)      
    plot1.plot(positions_US, Ded12[0], linestyle='--', color='#D1A47E', linewidth=1, label='Rep12-Rif_Rep12') 
    plot1.plot(positions_DS, Ded12[1], linestyle='--', color='#D1A47E', linewidth=1)
    plot1.plot(positions_Body, Ded12[2], linestyle='--', color='#D1A47E', linewidth=2)
    #plot1.set_ylim(-0.6, 0.5) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoA deduced fold enrichment', size=20)
    plot1.set_title(f'Deduced FE over the {genes_set_name}s', size=20)     
    plt.savefig(f'{output_path}Figures\\Averaged_Ded_FE_over_{genes_set_name}s_{win_width}bp_nd_with_body_{length}.svg', dpi=400, format='svg', figsize=(10, 6))   
    #plt.show()
    plt.close()      
    return


#Window_width=200
#Gene_body_length=1000

#Plot FE Rif and noRif.
#FE_Gene_Rep12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_Gene_Rif_Rep12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_Gene_Ded12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_rRNA_Rep12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_rRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_rRNA_Ded12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_tRNA_Rep12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_tRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_tRNA_Ded12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_ncRNA_Rep12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_ncRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_ncRNA_Ded12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#plot_averaged_FE_Dev_FE(FE_Gene_Rep12, FE_Gene_Rif_Rep12, FE_Gene_Ded12, 'Gene', Output_path, Window_width, Gene_body_length)
#plot_averaged_FE_Dev_FE(FE_rRNA_Rep12, FE_rRNA_Rif_Rep12, FE_rRNA_Ded12, 'rRNA', Output_path, Window_width, Gene_body_length)
#plot_averaged_FE_Dev_FE(FE_tRNA_Rep12, FE_tRNA_Rif_Rep12, FE_tRNA_Ded12, 'tRNA', Output_path, Window_width, Gene_body_length)
#plot_averaged_FE_Dev_FE(FE_ncRNA_Rep12, FE_ncRNA_Rif_Rep12, FE_ncRNA_Ded12, 'ncRNA', Output_path, Window_width, Gene_body_length)


#Window_width=5000
#Gene_body_length=5000

#Plot FE Rif and noRif.
#FE_Gene_Rep12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_Gene_Rif_Rep12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_Gene_Ded12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_rRNA_Rep12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_rRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_rRNA_Ded12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_tRNA_Rep12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_tRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_tRNA_Ded12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_ncRNA_Rep12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_ncRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_ncRNA_Ded12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#plot_averaged_FE_Dev_FE(FE_Gene_Rep12, FE_Gene_Rif_Rep12, FE_Gene_Ded12, 'Gene', Output_path, Window_width, Gene_body_length)
#plot_averaged_FE_Dev_FE(FE_rRNA_Rep12, FE_rRNA_Rif_Rep12, FE_rRNA_Ded12, 'rRNA', Output_path, Window_width, Gene_body_length)
#plot_averaged_FE_Dev_FE(FE_tRNA_Rep12, FE_tRNA_Rif_Rep12, FE_tRNA_Ded12, 'tRNA', Output_path, Window_width, Gene_body_length)
#plot_averaged_FE_Dev_FE(FE_ncRNA_Rep12, FE_ncRNA_Rif_Rep12, FE_ncRNA_Ded12, 'ncRNA', Output_path, Window_width, Gene_body_length)


#Window_width=15000
#Gene_body_length=5000

#Plot FE Rif and noRif.
#FE_Gene_Rep12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_Gene_Rif_Rep12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_Gene_Ded12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_rRNA_Rep12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_rRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_rRNA_Ded12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_tRNA_Rep12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_tRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_tRNA_Ded12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_ncRNA_Rep12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', dict_of_wigs['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_ncRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_ncRNA_Ded12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#plot_averaged_FE_Dev_FE(FE_Gene_Rep12, FE_Gene_Rif_Rep12, FE_Gene_Ded12, 'Gene', Output_path, Window_width, Gene_body_length)
#plot_averaged_FE_Dev_FE(FE_rRNA_Rep12, FE_rRNA_Rif_Rep12, FE_rRNA_Ded12, 'rRNA', Output_path, Window_width, Gene_body_length)
#plot_averaged_FE_Dev_FE(FE_tRNA_Rep12, FE_tRNA_Rif_Rep12, FE_tRNA_Ded12, 'tRNA', Output_path, Window_width, Gene_body_length)
#plot_averaged_FE_Dev_FE(FE_ncRNA_Rep12, FE_ncRNA_Rif_Rep12, FE_ncRNA_Ded12, 'ncRNA', Output_path, Window_width, Gene_body_length)



#######
#Parses WIG file with FE over TUs.
#######

def wig_FE_over_genes_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] in ['track']:
            ww_l=line[2].split('=')[1].rstrip('"').lstrip('"').split('_')
            win_width=int(ww_l[0])
            length=int(ww_l[1])
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    print(win_width, length)
    return NE_values, win_width, length

#######
#Returns smoothed tracks.
#######

def Smoothing(ends, window):
    smoothed=[]
    #Calculating the value for the first position
    sm=0.0
    window_float=float(window)
    sm+=np.mean(ends[:2*window])
    smoothed.append(sm)
    #Calculating values for the part of the array remains
    for i in range(len(ends)-2*window):
        sm+=(ends[i+2*window]-ends[i])/(2*window_float)
        smoothed.append(sm)
    return smoothed

#######
#Plot Ded FE for all groups of genes together.
#######

def plot_Ded_FE_all_gg(wig_in_dict, sm_window, output_path):
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=0
    length=0
    for name, file in wig_in_dict.items():
        print(name, file)
        data=wig_FE_over_genes_parsing(file)
        dict_of_wigs[name]=data[0]
        win_width=data[1]
        length=data[2]
    positions=np.arange(-win_width, win_width+length, 1)    
    print(win_width, length)
    
    #Plot FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    #Genes
    plot1.plot(positions, dict_of_wigs['Gene'], linestyle='-', color='#1411BB', linewidth=2, alpha=0.8, label='Gene (4043)', zorder=10)
    #rRNA
    plot1.plot(positions, dict_of_wigs['rRNA'], linestyle='-', color='#5FE400', linewidth=1, alpha=0.8, label='rRNA (22)', zorder=9)
    #tRNA
    plot1.plot(positions, dict_of_wigs['tRNA'], linestyle='-', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA (86)', zorder=8)
    #ncRNA
    plot1.plot(positions, dict_of_wigs['ncRNA'], linestyle='-', color='#FFC000', linewidth=1, alpha=0.8, label='ncRNA (18)', zorder=7)
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,1000).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(1000,win_width+1,1000).tolist()
    plot1.set_xticklabels(ticks_lables1)
    
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)
    
    plot1.set_yticks([0], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)
    
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoA fold enrichment', size=20)
    plot1.set_title(f'-Rif_FE - +Rif_FE over the different groups of genes', size=20)     
    plt.savefig(f'{output_path}Figures\\Averaged_Ded_FE_over_different_groups_of_genes_{win_width}bp_nd_with_body_{length}.png', dpi=400, figsize=(10, 6))   
    #plt.show()
    plt.close()     
    
    #Smoothing.
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+length-(sm_window-1), 1)  
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    #Genes
    plot1.plot(positions_sm, dict_of_wigs_sm['Gene'], linestyle='-', color='#1411BB', linewidth=2, alpha=0.8, label='Gene (4043)', zorder=10)
    #rRNA
    plot1.plot(positions_sm, dict_of_wigs_sm['rRNA'], linestyle='-', color='#5FE400', linewidth=1, alpha=0.8, label='rRNA (22)', zorder=9)
    #tRNA
    plot1.plot(positions_sm, dict_of_wigs_sm['tRNA'], linestyle='-', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA (86)', zorder=8)
    #ncRNA
    plot1.plot(positions_sm, dict_of_wigs_sm['ncRNA'], linestyle='-', color='#FFC000', linewidth=1, alpha=0.8, label='ncRNA (18)', zorder=7)    
    
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,1000).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(1000,win_width+1,1000).tolist()
    plot1.set_xticklabels(ticks_lables1)
    
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)   
    
    plot1.set_yticks([0], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)    
    
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoA fold enrichment', size=20)
    plot1.set_title(f'-Rif_FE - +Rif_FE over the different groups of genes', size=20)     
    plt.savefig(f'{output_path}Figures\\Averaged_Ded_FE_over_different_groups_of_genes_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(10, 6))   
    #plt.show()
    plt.close()    
    
    return

#plot_Ded_FE_all_gg(FE_input_dict_15000_5000, 50, Output_path)
#plot_Ded_FE_all_gg(FE_input_dict_15000_5000, 100, Output_path)
#plot_Ded_FE_all_gg(FE_input_dict_15000_5000, 200, Output_path)

#plot_Ded_FE_all_gg(FE_input_dict_5000_5000, 50, Output_path)
#plot_Ded_FE_all_gg(FE_input_dict_5000_5000, 100, Output_path)
#plot_Ded_FE_all_gg(FE_input_dict_5000_5000, 200, Output_path)

#plot_Ded_FE_all_gg(FE_input_dict_5000_1000, 50, Output_path)
#plot_Ded_FE_all_gg(FE_input_dict_5000_1000, 100, Output_path)
#plot_Ded_FE_all_gg(FE_input_dict_5000_1000, 200, Output_path)


#######
#Plot FE for all groups of genes together.
#######

def plot_FE_all_gg_Rif_no_Rif(wig_in_dict, sm_window, output_path):
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=0
    length=0
    for name, file in wig_in_dict.items():
        print(name, file)
        data=wig_FE_over_genes_parsing(file)
        dict_of_wigs[name]=data[0]
        win_width=data[1]
        length=data[2]
    positions=np.arange(-win_width, win_width+length, 1)    
    print(win_width, length)
    
    #Plot FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    #Genes
    plot1.plot(positions, dict_of_wigs['Gene'], linestyle='-', color='#1411BB', linewidth=2, alpha=0.8, label='Gene (4043)', zorder=10)
    plot1.plot(positions, dict_of_wigs['Gene Rif'], linestyle='--', color='#1411BB', linewidth=1, alpha=0.8, label='Gene Rif (4043)', zorder=9)
    #rRNA
    plot1.plot(positions, dict_of_wigs['rRNA'], linestyle='-', color='#EA003F', linewidth=1.5, alpha=0.8, label='rRNA (22)', zorder=8)
    plot1.plot(positions, dict_of_wigs['rRNA Rif'], linestyle='--', color='#EA003F', linewidth=0.8, alpha=0.8, label='rRNA Rif (22)', zorder=7)
    #tRNA
    plot1.plot(positions, dict_of_wigs['tRNA'], linestyle='-', color='#E692A9', linewidth=1, alpha=2, label='tRNA (86)', zorder=6)
    plot1.plot(positions, dict_of_wigs['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)
    #ncRNA
    plot1.plot(positions, dict_of_wigs['ncRNA'], linestyle='-', color='#FFC000', linewidth=1, alpha=1.5, label='ncRNA (18)', zorder=4)
    plot1.plot(positions, dict_of_wigs['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,5000).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(5000,win_width+1,5000).tolist()
    plot1.set_xticklabels(ticks_lables1)
    
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)
    
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)
    
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoA fold enrichment', size=20)
    plot1.set_title(f'-Rif_FE and +Rif_FE over tRNA and ncRNA', size=20)     
    plt.savefig(f'{output_path}Figures\\Averaged_noRif_and_Rif_FE_over_tRNA_and_ncRNA_{win_width}bp_nd_with_body_{length}.png', dpi=400, figsize=(10, 6))   
    #plt.show()
    plt.close()     
    
    #Smoothing.
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+length-(sm_window-1), 1)  
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)  
    #Genes
    plot1.plot(positions_sm, dict_of_wigs_sm['Gene'], linestyle='-', color='#1411BB', linewidth=2, alpha=0.8, label='Gene (4043)', zorder=10)
    plot1.plot(positions_sm, dict_of_wigs_sm['Gene Rif'], linestyle='--', color='#1411BB', linewidth=1, alpha=0.8, label='Gene Rif (4043)', zorder=9)
    #rRNA
    plot1.plot(positions_sm, dict_of_wigs_sm['rRNA'], linestyle='-', color='#EA003F', linewidth=1.5, alpha=0.8, label='rRNA (22)', zorder=8)
    plot1.plot(positions_sm, dict_of_wigs_sm['rRNA Rif'], linestyle='--', color='#EA003F', linewidth=0.8, alpha=0.8, label='rRNA Rif (22)', zorder=7)
    #tRNA
    plot1.plot(positions_sm, dict_of_wigs_sm['tRNA'], linestyle='-', color='#E692A9', linewidth=2, alpha=0.8, label='tRNA (86)', zorder=6)
    plot1.plot(positions_sm, dict_of_wigs_sm['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)
    #ncRNA
    plot1.plot(positions_sm, dict_of_wigs_sm['ncRNA'], linestyle='-', color='#FFC000', linewidth=1.5, alpha=0.8, label='ncRNA (18)', zorder=4)
    plot1.plot(positions_sm, dict_of_wigs_sm['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3) 
    
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,5000).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(5000,win_width+1,5000).tolist()
    plot1.set_xticklabels(ticks_lables1)
    
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)   
    
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)    
    
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoA fold enrichment', size=20)
    plot1.set_title(f'-Rif_FE and +Rif_FE over tRNA and ncRNA', size=20)     
    plt.savefig(f'{output_path}Figures\\Averaged_noRif_and_Rif_FE_over_tRNA_and_ncRNA_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(10, 6))   
    #plt.show()
    plt.close()    
    
    return


#plot_FE_all_gg_Rif_no_Rif(Rif_noRif_FE_input_dict_5000_1000, 50, Output_path)
#plot_FE_all_gg_Rif_no_Rif(Rif_noRif_FE_input_dict_5000_1000, 100, Output_path)
#plot_FE_all_gg_Rif_no_Rif(Rif_noRif_FE_input_dict_5000_1000, 200, Output_path)

#plot_FE_all_gg_Rif_no_Rif(Rif_noRif_FE_input_dict_5000_5000, 50, Output_path)
#plot_FE_all_gg_Rif_no_Rif(Rif_noRif_FE_input_dict_5000_5000, 100, Output_path)
#plot_FE_all_gg_Rif_no_Rif(Rif_noRif_FE_input_dict_5000_5000, 200, Output_path)

#plot_FE_all_gg_Rif_no_Rif(Rif_noRif_FE_input_dict_15000_5000, 50, Output_path)
#plot_FE_all_gg_Rif_no_Rif(Rif_noRif_FE_input_dict_15000_5000, 100, Output_path)
#plot_FE_all_gg_Rif_no_Rif(Rif_noRif_FE_input_dict_15000_5000, 200, Output_path)


#######
#Write .tab-file with FE (US, GB, DS), Rif FE (US, GB, DS), and Ded FE (US, GB, DS) information for each gene in a set.
#######

def write_genes_FE_Rif_FE_Ded_FE(dict1, dict2, dict3, dict4, dict5, dict6, dict7, dict8, dict9, path_out):
    fileout=open(path_out, 'w')
    fileout.write('Gene_name\tStart\tEnd\tStrand\tTopoA_noRif_FE_US\tTopoA_noRif_FE_GB\tTopoA_noRif_FE_DS\tTopoA_Rif_FE_US\tTopoA_Rif_FE_GB\tTopoA_Rif_FE_DS\tTopoA_Ded_FE_US\tTopoA_Ded_FE_GB\tTopoA_Ded_FE_DS\n')
    for k, v in dict1.items():
        fileout.write(f'{k}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[0]}\t{dict2[k][0]}\t{dict3[k][0]}\t{dict4[k][0]}\t{dict5[k][0]}\t{dict6[k][0]}\t{dict7[k][0]}\t{dict8[k][0]}\t{dict9[k][0]}\n')
    fileout.close()
    return

#######
#Wrapper: reads dict of .wig-files, plots distribution of FE.
#######

def FE_analysis_ComRep(Wig_input_dict, pathout):
    #WIG parsing
    dict_of_wigs={}
    for name, file in Wig_input_dict.items():
        dict_of_wigs[name]=wig_parsing(file)
    
    #Plot distributions of FE and FE deduction.
    #plot_FE_dist_GW(dict_of_wigs['Rep12_FE'], 'Rep12_FE', dict_of_wigs['Rif_Rep12_FE'], 'Rif_Rep12_FE', dict_of_wigs['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', f'{pathout}Figures\Rep12_FE_vs_Rif_Rep12_FE.png')
    return dict_of_wigs

#Combined replicas and Ded.
#FE_tracks_dict_CompRep=FE_analysis_ComRep(Wig_av_FE_input_dict, Output_path)
#Window_width=15000
#Gene_body_length=5000

#FE_Gene_Rep12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', FE_tracks_dict_CompRep['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_Gene_Rif_Rep12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', FE_tracks_dict_CompRep['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_Gene_Rif_Ded12=genes_and_FE(genome_annotation[0]['Gene'], 'Gene', FE_tracks_dict_CompRep['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_rRNA_Rep12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', FE_tracks_dict_CompRep['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_rRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', FE_tracks_dict_CompRep['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_rRNA_Rif_Ded12=genes_and_FE(genome_annotation[0]['rRNA'], 'rRNA', FE_tracks_dict_CompRep['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_tRNA_Rep12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', FE_tracks_dict_CompRep['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_tRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', FE_tracks_dict_CompRep['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_tRNA_Rif_Ded12=genes_and_FE(genome_annotation[0]['tRNA'], 'tRNA', FE_tracks_dict_CompRep['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#FE_ncRNA_Rep12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', FE_tracks_dict_CompRep['Rep12_FE'], 'Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_ncRNA_Rif_Rep12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', FE_tracks_dict_CompRep['Rif_Rep12_FE'], 'Rif_Rep12_FE', Output_path, Deletions, Window_width, Gene_body_length)
#FE_ncRNA_Rif_Ded12=genes_and_FE(genome_annotation[0]['ncRNA'], 'ncRNA', FE_tracks_dict_CompRep['Ded_Rep12_Rif12_FE'], 'Ded_Rep12_Rif12_FE', Output_path, Deletions, Window_width, Gene_body_length)

#write_genes_FE_Rif_FE_Ded_FE(FE_Gene_Rep12[6], FE_Gene_Rep12[8], FE_Gene_Rep12[7], FE_Gene_Rif_Rep12[6], FE_Gene_Rif_Rep12[8], FE_Gene_Rif_Rep12[7], FE_Gene_Rif_Ded12[6], FE_Gene_Rif_Ded12[8], FE_Gene_Rif_Ded12[7], f'{Output_path}FE_of_genes\\noRif_Rif_Ded_FE_over_Genes_{Window_width}bp.txt')
#write_genes_FE_Rif_FE_Ded_FE(FE_rRNA_Rep12[6], FE_rRNA_Rep12[8], FE_rRNA_Rep12[7], FE_rRNA_Rif_Rep12[6], FE_rRNA_Rif_Rep12[8], FE_rRNA_Rif_Rep12[7], FE_rRNA_Rif_Ded12[6], FE_rRNA_Rif_Ded12[8], FE_rRNA_Rif_Ded12[7], f'{Output_path}FE_of_genes\\noRif_Rif_Ded_FE_over_rRNAs_{Window_width}bp.txt')
#write_genes_FE_Rif_FE_Ded_FE(FE_tRNA_Rep12[6], FE_tRNA_Rep12[8], FE_tRNA_Rep12[7], FE_tRNA_Rif_Rep12[6], FE_tRNA_Rif_Rep12[8], FE_tRNA_Rif_Rep12[7], FE_tRNA_Rif_Ded12[6], FE_tRNA_Rif_Ded12[8], FE_tRNA_Rif_Ded12[7], f'{Output_path}FE_of_genes\\noRif_Rif_Ded_FE_over_tRNAs_{Window_width}bp.txt')
#write_genes_FE_Rif_FE_Ded_FE(FE_ncRNA_Rep12[6], FE_ncRNA_Rep12[8], FE_ncRNA_Rep12[7], FE_ncRNA_Rif_Rep12[6], FE_ncRNA_Rif_Rep12[8], FE_ncRNA_Rif_Rep12[7], FE_ncRNA_Rif_Ded12[6], FE_ncRNA_Rif_Ded12[8], FE_ncRNA_Rif_Ded12[7], f'{Output_path}FE_of_genes\\noRif_Rif_Ded_FE_over_ncRNAs_{Window_width}bp.txt')

#######
#Makes abs of all elements in an array.
#######

def make_abs(ar):
    ar_abs=[]
    for el in ar:
        ar_abs.append(abs(el))
    return ar_abs

#######
#Return mean FE for IG regions.
#######

def IG_FE(FE_wig, genes_wig):
    IG_values=[]
    for i in range(len(FE_wig)):
        if genes_wig[i]==0:
            IG_values.append(FE_wig[i])
    return IG_values

#######
#Concatenate ar of ars.
#######

def conc_ars(ar):
    conc=[]
    for in_ar in ar:
        conc+=in_ar
    return conc

#######
#Plot distributions of FE and Ded FE over IG and GB.
#######

def plot_GB_vs_IG(FE_IG, Rif_FE_IG, Ded_FE_IG, FE_GB, Rif_FE_GB, Ded_FE_GB, outpath):
    #Plot distribution of FE values over IG.
    fig=plt.figure(figsize=(15, 5), dpi=100)
    bins0=np.arange(float(min(FE_IG+FE_GB)), float(max(FE_IG+FE_GB)), 1)
    plot0=plt.subplot2grid((2,3),(0,0), rowspan=1, colspan=1)
    plot0.hist(FE_IG, bins0, color='#ff878b', edgecolor='black', alpha=0.8)
    plot0.set_yscale('log')
    plot0.set_xlabel('Fold enrichment', size=13)
    plot0.set_ylabel('Number of positions', size=13)
    plot0.set_title('FE in IG', size=18)  
    #Plot distribution of FE values over GB.
    plot0a=plt.subplot2grid((2,3),(1,0), rowspan=1, colspan=1)
    plot0a.hist(FE_GB, bins0, color='#ff878b', edgecolor='black', alpha=0.8)
    plot0a.set_yscale('log')
    plot0a.set_xlabel('Fold enrichment', size=13)
    plot0a.set_ylabel('Number of positions', size=13)
    plot0a.set_title('FE in GB', size=18)   
    
    #Plot distribution of Rif FE values over IG.  
    bins1=np.arange(float(min(Rif_FE_IG+Rif_FE_GB)), float(max(Rif_FE_IG+Rif_FE_GB)), 1)
    plot1=plt.subplot2grid((2,3),(0,1), rowspan=1, colspan=1)     
    plot1.hist(Rif_FE_IG, bins1, color='#ffce91', edgecolor='black', alpha=0.5)
    plot1.set_yscale('log')
    plot1.set_xlabel('Fold enrichment', size=13)
    plot1.set_ylabel('Number of positions', size=13)
    plot1.set_title('Rif FE in IG', size=18) 
    #Plot distribution of Rif FE values over GB.
    plot1a=plt.subplot2grid((2,3),(1,1), rowspan=1, colspan=1)     
    plot1a.hist(Rif_FE_GB, bins1, color='#ffce91', edgecolor='black', alpha=0.5)
    plot1a.set_yscale('log')
    plot1a.set_xlabel('Fold enrichment', size=13)
    plot1a.set_ylabel('Number of positions', size=13)
    plot1a.set_title('Rif FE in GB', size=18) 
    
    #Plot distribution of Ded FE values over IG.
    bins2=np.arange(float(min(Ded_FE_IG+Ded_FE_GB)), float(max(Ded_FE_IG+Ded_FE_GB)), 1)
    plot2=plt.subplot2grid((2,3),(0,2), rowspan=1, colspan=1) 
    plot2.hist(Ded_FE_IG, bins2, color='#7FCE79', edgecolor='black', alpha=0.5)
    plot2.set_yscale('log')
    plot2.set_xlabel('Delta Fold enrichment', size=13)
    plot2.set_ylabel('Number of positions', size=13)
    plot2.set_title('Ded FE in IG', size=18)  
    #Plot distribution of Ded FE values over GB.
    plot2=plt.subplot2grid((2,3),(1,2), rowspan=1, colspan=1) 
    plot2.hist(Ded_FE_GB, bins2, color='#7FCE79', edgecolor='black', alpha=0.5)
    plot2.set_yscale('log')
    plot2.set_xlabel('Delta Fold enrichment', size=13)
    plot2.set_ylabel('Number of positions', size=13)
    plot2.set_title('Ded FE in GB', size=18)      
    
    plt.tight_layout()
    #plt.show()
    plt.savefig(outpath, dpi=300, figsize=(15, 10))
    plt.close()     
    return

#######
#Wrapper: takes FE .wig-files, FE over gene bodies, computes FE over intergenic regions (IR), plots FE over GB and IR.
#######

def GB_vs_IG(Rep_FE, Rif_Rep_FE, Ded_FE, FE_GB_ar, Rif_FE_GB_ar, Ded_FE_GB_ar, set_name, genome_annot, outpath):
    #Makes wigs for each annotation.
    Gene_wig=make_binary_annotation(genome_annot[0]['Gene'])
    rRNA_wig=make_binary_annotation(genome_annot[0]['rRNA'])
    tRNA_wig=make_binary_annotation(genome_annot[0]['tRNA'])
    ncRNA_wig=make_binary_annotation(genome_annot[0]['ncRNA'])   
    All_genes_wig_abs=(np.array(make_abs(Gene_wig)) + np.array(make_abs(rRNA_wig)) + np.array(make_abs(tRNA_wig)) + np.array(make_abs(ncRNA_wig))).tolist()
    
    #Prepare FE of IG (arrays).
    FE_IG=IG_FE(Rep_FE, All_genes_wig_abs)
    Rif_FE_IG=IG_FE(Rif_Rep_FE, All_genes_wig_abs)
    Ded_FE_IG=IG_FE(Ded_FE, All_genes_wig_abs)
    
    #Prepare GB arrays:
    FE_GB=conc_ars(FE_GB_ar)
    Rif_FE_GB=conc_ars(Rif_FE_GB_ar)
    Ded_FE_GB=conc_ars(Ded_FE_GB_ar)
    
    #Plot distributions.
    plot_GB_vs_IG(FE_IG, Rif_FE_IG, Ded_FE_IG, FE_GB, Rif_FE_GB, Ded_FE_GB, f'{outpath}Figures\\FE_distributions_over_GB_and_IG_{set_name}.png')
    return

#GB_vs_IG(FE_tracks_dict_CompRep['Rep12_FE'], FE_tracks_dict_CompRep['Rif_Rep12_FE'], FE_tracks_dict_CompRep['Ded_Rep12_Rif12_FE'], 
#         FE_Gene_Rep12[2], FE_Gene_Rif_Rep12[2], FE_Gene_Rif_Ded12[2], 'Gene', genome_annotation, Output_path)



#######
#Reads annotation of particular set of genes .tab BroadPeak-like (determined on a basis of expression level).
#######

def parse_expression_annotation(annot_inpath):
    genes_annotation={}
    filein=open(annot_inpath, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        if line[0] not in ['GeneID', 'OperonID']:
            TU_name=line[1].lstrip('"').rstrip(';"')
            TU_start=int(line[2])
            TU_end=int(line[3])
            TU_strand=line[4]
            TU_expression=float(line[5].replace(',','.'))
            genes_annotation[TU_name]=[TU_start, TU_end, TU_strand, TU_expression]
    filein.close()            
    return genes_annotation

#######
#Wrapper: reads annotation, computes FE over annotations' TUs, plots distribution of FE, tracks of FE over TUs.
#######

def Wrap_expr_dep_analysis(annot_dict_inpath, deletions_inpath, dict_of_FE_tracks, Output_path, Window_width, Gene_body_length):
    #Parsing annotations of TU sets, prepare dict of annotations.
    Dict_of_annotations={}
    for set_name, annot_inpath in annot_dict_inpath.items():
        genes_annotation=parse_expression_annotation(annot_inpath)
        Dict_of_annotations[set_name]=genes_annotation
    
    #Computing FE over all TU sets supplied.
    FE_over_TUs={}
    for TU_set_name, TU_set in Dict_of_annotations.items():
        for FE_track_name, FE_track in dict_of_FE_tracks.items():
            print(f'Dealing with {FE_track_name} and {TU_set_name}')
            FE_of_TU_set=genes_and_FE(TU_set, TU_set_name, FE_track, FE_track_name, Output_path, deletions_inpath, Window_width, Gene_body_length)  
            FE_over_TUs[f'{TU_set_name}_{FE_track_name}']=FE_of_TU_set
    return

#Wrap_expr_dep_analysis(Annot_dict_inpath, Deletions, FE_tracks_dict_CompRep, Output_path, Window_width, Gene_body_length)

#######
#Scale a part of vector by scaling factor.
#######

def scale_part_of_seq(ar, start, end, scaling_factor):
    scaled_part=np.array(ar[start:end])/scaling_factor
    ar_sc=ar[:start] + scaled_part.tolist() + ar[end:]
    print(ar[start:start+10])
    print(ar_sc[start:start+10])
    return ar_sc

#######
#Plot FE for all groups of genes by expression together.
#######

def plot_FE_all_expression_gg_Rif_no_Rif(wig_in_dict, sm_window, output_path):
    #Number of genes within sets.
    TU_sets_v={'All_genes' : 4119, 'HEG_no_tRNA_rRNA_270' : 269, 'LEG_270' : 270, 'HEG_370' : 369, 'LEG_370' : 370,
               'All_genes_Rif' : 4119, 'HEG_no_tRNA_rRNA_270_Rif' : 269, 'LEG_270_Rif' : 270, 'HEG_370_Rif' : 369, 'LEG_370_Rif' : 370,
               'All_operons' : 2327, 'HEO_no_tRNA_rRNA_144' : 143, 'LEO_144' : 144, 'HEO_186' : 185, 'LEO_186' : 186,
               'All_operons_Rif' : 2327, 'HEO_no_tRNA_rRNA_144_Rif' : 143, 'LEO_144_Rif' : 144, 'HEO_186_Rif' : 185, 'LEO_186_Rif' : 186,
               'rRNA_operons' : 7}
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=0
    length=0
    for name, file in wig_in_dict.items():
        print(name, file)
        data=wig_FE_over_genes_parsing(file)
        win_width=data[1]
        length=data[2]
        restored_wig=scale_part_of_seq(data[0], win_width, win_width+length, TU_sets_v[name])
        dict_of_wigs[name]=restored_wig        
    positions=np.arange(-win_width, win_width+length, 1)    
    print(win_width, length)
       
    
    #Plot FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    ##Genes below
    #LEG_270
    #plot1.plot(positions, dict_of_wigs['LEG_270'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEG ({TU_sets_v["LEG_270"]})', zorder=6)
    #plot1.plot(positions, dict_of_wigs['LEG_270_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEG Rif ({TU_sets_v["LEG_270_Rif"]})', zorder=5)   
    #LEG_370
    #plot1.plot(positions, dict_of_wigs['LEG_370'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEG ({TU_sets_v["LEG_370"]})', zorder=6)
    #plot1.plot(positions, dict_of_wigs['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
    #All_genes
    #plot1.plot(positions, dict_of_wigs['All_genes'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All genes ({TU_sets_v["All_genes"]})', zorder=10)
    #plot1.plot(positions, dict_of_wigs['All_genes_Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All genes Rif ({TU_sets_v["All_genes_Rif"]})', zorder=9)
    #HEG_270
    #plot1.plot(positions, dict_of_wigs['HEG_no_tRNA_rRNA_270'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEG no rRNA, tRNA ({TU_sets_v["HEG_no_tRNA_rRNA_270"]})', zorder=8)
    #plot1.plot(positions, dict_of_wigs['HEG_no_tRNA_rRNA_270_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEG no rRNA, tRNA Rif ({TU_sets_v["HEG_no_tRNA_rRNA_270_Rif"]})', zorder=7)
    #HEG_370
    #plot1.plot(positions, dict_of_wigs['HEG_370'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEG ({TU_sets_v["HEG_370"]})', zorder=4)
    #plot1.plot(positions, dict_of_wigs['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)
    ##Operons below.
    #LEO_144
    plot1.plot(positions, dict_of_wigs['LEO_144'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEO ({TU_sets_v["LEO_144"]})', zorder=6)
    plot1.plot(positions, dict_of_wigs['LEO_144_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEO Rif ({TU_sets_v["LEO_144_Rif"]})', zorder=5)   
    #LEO_186
    #plot1.plot(positions, dict_of_wigs['LEO_186'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEO ({TU_sets_v["LEO_186"]})', zorder=6)
    #plot1.plot(positions, dict_of_wigs['LEO_186_Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
    #All_operons
    plot1.plot(positions, dict_of_wigs['All_operons'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All operons ({TU_sets_v["All_operons"]})', zorder=10)
    plot1.plot(positions, dict_of_wigs['All_operons_Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All operons Rif ({TU_sets_v["All_operons_Rif"]})', zorder=9)
    #HEO_144
    plot1.plot(positions, dict_of_wigs['HEO_no_tRNA_rRNA_144'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEO no rRNA, tRNA ({TU_sets_v["HEO_no_tRNA_rRNA_144"]})', zorder=8)
    plot1.plot(positions, dict_of_wigs['HEO_no_tRNA_rRNA_144_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEO no rRNA, tRNA Rif ({TU_sets_v["HEO_no_tRNA_rRNA_144_Rif"]})', zorder=7)
    #HEO_186 or rRNA operons
    #plot1.plot(positions, dict_of_wigs['rRNA_operons'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'rRNA operons ({TU_sets_v["rRNA_operons"]})', zorder=4)
    #plot1.plot(positions, dict_of_wigs['HEO_186'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEO ({TU_sets_v["HEO_186"]})', zorder=4)
    #plot1.plot(positions, dict_of_wigs['HEO_186_Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)    
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1) 
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoA fold enrichment', size=20)
    plot1.set_title(f'-Rif and +Rif over EDS of operons no dps', size=20)     
    plt.savefig(f'{output_path}Figures\\noRif_Rif_nodps_over_expression_og_{win_width}bp_nd_with_body_{length}.png', dpi=400, figsize=(10, 6))   
    #plt.show()
    plt.close()     
    
    #Smoothing.
    dict_of_wigs_sm={}
    for name, wig in dict_of_wigs.items():
        dict_of_wigs_sm[name]=Smoothing(wig, sm_window)
    positions_sm=np.arange(-win_width+sm_window, win_width+length-(sm_window-1), 1)  
        
    #Plot smoothed FE over genes.
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    ##Genes below
    #LEG_270
    #plot1.plot(positions_sm, dict_of_wigs_sm['LEG_270'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEG ({TU_sets_v["LEG_270"]})', zorder=6)
    #plot1.plot(positions_sm, dict_of_wigs_sm['LEG_270_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEG Rif ({TU_sets_v["LEG_270_Rif"]})', zorder=5)   
    #LEG_370
    #plot1.plot(positions_sm, dict_of_wigs_sm['LEG_370'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEG ({TU_sets_v["LEG_370"]})', zorder=6)
    #plot1.plot(positions_sm, dict_of_wigs_sm['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
    #All_genes
    #plot1.plot(positions_sm, dict_of_wigs_sm['All_genes'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All genes ({TU_sets_v["All_genes"]})', zorder=10)
    #plot1.plot(positions_sm, dict_of_wigs_sm['All_genes_Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All genes Rif ({TU_sets_v["All_genes_Rif"]})', zorder=9)
    #HEG_270
    #plot1.plot(positions_sm, dict_of_wigs_sm['HEG_no_tRNA_rRNA_270'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEG no rRNA, tRNA ({TU_sets_v["HEG_no_tRNA_rRNA_270"]})', zorder=8)
    #plot1.plot(positions_sm, dict_of_wigs_sm['HEG_no_tRNA_rRNA_270_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEG no rRNA, tRNA Rif ({TU_sets_v["HEG_no_tRNA_rRNA_270_Rif"]})', zorder=7)
    #HEG_370
    #plot1.plot(positions_sm, dict_of_wigs_sm['HEG_370'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEG ({TU_sets_v["HEG_370"]})', zorder=4)
    #plot1.plot(positions_sm, dict_of_wigs_sm['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)  
    ##Operons below.
    #LEO_144
    plot1.plot(positions_sm, dict_of_wigs_sm['LEO_144'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEO ({TU_sets_v["LEO_144"]})', zorder=6)
    plot1.plot(positions_sm, dict_of_wigs_sm['LEO_144_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEO Rif ({TU_sets_v["LEO_144_Rif"]})', zorder=5)   
    #LEO_186
    #plot1.plot(positions_sm, dict_of_wigs_sm['LEO_186'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEO ({TU_sets_v["LEO_186"]})', zorder=6)
    #plot1.plot(positions_sm, dict_of_wigs_sm['LEO_186_Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
    #All_operons
    plot1.plot(positions_sm, dict_of_wigs_sm['All_operons'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All operons ({TU_sets_v["All_operons"]})', zorder=10)
    plot1.plot(positions_sm, dict_of_wigs_sm['All_operons_Rif'], linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All operons Rif ({TU_sets_v["All_operons_Rif"]})', zorder=9)
    #HEO_144
    plot1.plot(positions_sm, dict_of_wigs_sm['HEO_no_tRNA_rRNA_144'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEO no rRNA, tRNA ({TU_sets_v["HEO_no_tRNA_rRNA_144"]})', zorder=8)
    plot1.plot(positions_sm, dict_of_wigs_sm['HEO_no_tRNA_rRNA_144_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEO no rRNA, tRNA Rif ({TU_sets_v["HEO_no_tRNA_rRNA_144_Rif"]})', zorder=7)
    #HEO_186 or rRNA operons
    #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA_operons'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'rRNA operons ({TU_sets_v["rRNA_operons"]})', zorder=4)
    #plot1.plot(positions_sm, dict_of_wigs_sm['HEO_186'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEO ({TU_sets_v["HEO_186"]})', zorder=4)
    #plot1.plot(positions_sm, dict_of_wigs_sm['HEO_186_Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)      
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=0.5)   
    plot1.set_yticks([1], minor='True')
    plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=0.5)    
    plot1.legend(fontsize=12)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel('TopoA fold enrichment', size=20)
    plot1.set_title(f'-Rif and +Rif over EDS of operons no dps', size=20)     
    plt.savefig(f'{output_path}Figures\\noRif_Rif_nodps_over_expression_og_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(10, 6))   
    #plt.show()
    plt.close()    
    return


############
####Genes:
############
noRif_FE_input_dict={'All_genes' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rep12_FE_over_All_genes_width_15000bp_gb_5000bp.wig',
                     'HEG_370' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rep12_FE_over_HEG_370_width_15000bp_gb_5000bp.wig',
                     'LEG_370' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rep12_FE_over_LEG_370_width_15000bp_gb_5000bp.wig',
                     'HEG_no_tRNA_rRNA_270' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rep12_FE_over_HEG_no_tRNA_rRNA_270_width_15000bp_gb_5000bp.wig',
                     'LEG_270' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rep12_FE_over_LEG_270_width_15000bp_gb_5000bp.wig',
                     }
#plot_FE_all_expression_gg_Rif_no_Rif(noRif_FE_input_dict, 200, Output_path)

Rif_FE_input_dict={'All_genes' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rif_Rep12_FE_over_All_genes_width_15000bp_gb_5000bp.wig',
                     'HEG_370' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rif_Rep12_FE_over_HEG_370_width_15000bp_gb_5000bp.wig',
                     'LEG_370' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rif_Rep12_FE_over_LEG_370_width_15000bp_gb_5000bp.wig',
                     'HEG_no_tRNA_rRNA_270' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rif_Rep12_FE_over_HEG_no_tRNA_rRNA_270_width_15000bp_gb_5000bp.wig',
                     'LEG_270' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rif_Rep12_FE_over_LEG_270_width_15000bp_gb_5000bp.wig',
                     }
#plot_FE_all_expression_gg_Rif_no_Rif(Rif_FE_input_dict, 200, Output_path)

nodps_noRif_FE_input_dict={'All_genes' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_All_genes_width_15000bp_gb_5000bp.wig',
                     'HEG_370' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_HEG_370_width_15000bp_gb_5000bp.wig',
                     'LEG_370' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_LEG_370_width_15000bp_gb_5000bp.wig',
                     'HEG_no_tRNA_rRNA_270' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_HEG_no_tRNA_rRNA_270_width_15000bp_gb_5000bp.wig',
                     'LEG_270' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_LEG_270_width_15000bp_gb_5000bp.wig',
                     }
#plot_FE_all_expression_gg_Rif_no_Rif(nodps_noRif_FE_input_dict, 200, Output_path)

nodps_Rif_FE_input_dict={'All_genes' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_All_genes_width_15000bp_gb_5000bp.wig',
                     'HEG_370' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_HEG_370_width_15000bp_gb_5000bp.wig',
                     'LEG_370' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_LEG_370_width_15000bp_gb_5000bp.wig',
                     'HEG_no_tRNA_rRNA_270' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_HEG_no_tRNA_rRNA_270_width_15000bp_gb_5000bp.wig',
                     'LEG_270' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_LEG_270_width_15000bp_gb_5000bp.wig',
                     }
#plot_FE_all_expression_gg_Rif_no_Rif(nodps_Rif_FE_input_dict, 200, Output_path)

nodps_noRif_and_Rif_FE_input_dict={'All_genes' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_All_genes_width_15000bp_gb_5000bp.wig',
                                   'HEG_no_tRNA_rRNA_270' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_HEG_no_tRNA_rRNA_270_width_15000bp_gb_5000bp.wig',
                                   'LEG_270' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_LEG_270_width_15000bp_gb_5000bp.wig',
                                   'All_genes_Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_All_genes_width_15000bp_gb_5000bp.wig',
                                   'HEG_no_tRNA_rRNA_270_Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_HEG_no_tRNA_rRNA_270_width_15000bp_gb_5000bp.wig',
                                   'LEG_270_Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_LEG_270_width_15000bp_gb_5000bp.wig',                                   
                                   }
#plot_FE_all_expression_gg_Rif_no_Rif(nodps_noRif_and_Rif_FE_input_dict, 200, Output_path)

############
####Operons:
############
noRif_FE_input_dict_op={'All_operons' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rep12_FE_over_All_operons_width_15000bp_gb_5000bp.wig',
                     'HEO_186' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rep12_FE_over_HEO_186_width_15000bp_gb_5000bp.wig',
                     'LEO_186' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rep12_FE_over_LEO_186_width_15000bp_gb_5000bp.wig',
                     'HEO_no_tRNA_rRNA_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rep12_FE_over_HEO_no_tRNA_rRNA_144_width_15000bp_gb_5000bp.wig',
                     'LEO_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rep12_FE_over_LEO_144_width_15000bp_gb_5000bp.wig',
                     }
#plot_FE_all_expression_gg_Rif_no_Rif(noRif_FE_input_dict_op, 200, Output_path)

Rif_FE_input_dict_op={'All_operons' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rif_Rep12_FE_over_All_operons_width_15000bp_gb_5000bp.wig',
                     'HEO_186' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rif_Rep12_FE_over_HEO_186_width_15000bp_gb_5000bp.wig',
                     'LEO_186' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rif_Rep12_FE_over_LEO_186_width_15000bp_gb_5000bp.wig',
                     'HEO_no_tRNA_rRNA_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rif_Rep12_FE_over_HEO_no_tRNA_rRNA_144_width_15000bp_gb_5000bp.wig',
                     'LEO_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\TopoA_av_Rif_Rep12_FE_over_LEO_144_width_15000bp_gb_5000bp.wig',
                     }
#plot_FE_all_expression_gg_Rif_no_Rif(Rif_FE_input_dict_op, 200, Output_path)

nodps_noRif_FE_input_dict_op={'All_operons' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_All_operons_width_15000bp_gb_5000bp.wig',
                     'HEO_186' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_HEO_186_width_15000bp_gb_5000bp.wig',
                     'LEO_186' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_LEO_186_width_15000bp_gb_5000bp.wig',
                     'HEO_no_tRNA_rRNA_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_HEO_no_tRNA_rRNA_144_width_15000bp_gb_5000bp.wig',
                     'LEO_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_LEO_144_width_15000bp_gb_5000bp.wig',
                     }
#plot_FE_all_expression_gg_Rif_no_Rif(nodps_noRif_FE_input_dict_op, 200, Output_path)

nodps_Rif_FE_input_dict_op={'All_operons' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_All_operons_width_15000bp_gb_5000bp.wig',
                            'HEO_186' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_HEO_186_width_15000bp_gb_5000bp.wig',
                            'LEO_186' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_LEO_186_width_15000bp_gb_5000bp.wig',
                            'HEO_no_tRNA_rRNA_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_HEO_no_tRNA_rRNA_144_width_15000bp_gb_5000bp.wig',
                            'LEO_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_LEO_144_width_15000bp_gb_5000bp.wig',
                            }
#plot_FE_all_expression_gg_Rif_no_Rif(nodps_Rif_FE_input_dict_op, 200, Output_path)

nodps_noRif_and_Rif_FE_input_dict_op={'All_operons' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_All_operons_width_15000bp_gb_5000bp.wig',
                                      'HEO_no_tRNA_rRNA_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_HEO_no_tRNA_rRNA_144_width_15000bp_gb_5000bp.wig',
                                      'LEO_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_LEO_144_width_15000bp_gb_5000bp.wig',
                                      'All_operons_Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_All_operons_width_15000bp_gb_5000bp.wig',
                                      'HEO_no_tRNA_rRNA_144_Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_HEO_no_tRNA_rRNA_144_width_15000bp_gb_5000bp.wig',
                                      'LEO_144_Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rif_Rep12_FE_over_LEO_144_width_15000bp_gb_5000bp.wig',                                   
                                      }
#plot_FE_all_expression_gg_Rif_no_Rif(nodps_noRif_and_Rif_FE_input_dict_op, 200, Output_path)

noRif_FE_input_dict_extended_op={'All_operons' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_All_operons_width_15000bp_gb_5000bp.wig',
                                 'rRNA_operons' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_rRNA_operons_width_15000bp_gb_5000bp.wig',
                                 'HEO_no_tRNA_rRNA_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_HEO_no_tRNA_rRNA_144_width_15000bp_gb_5000bp.wig',
                                 'LEO_144' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TU_sets_by_expression\\No_dps_TopoA_av_Rep12_FE_over_LEO_144_width_15000bp_gb_5000bp.wig',
                                 }
#plot_FE_all_expression_gg_Rif_no_Rif(noRif_FE_input_dict_extended_op, 200, Output_path)


##################################
##Analysis of TopoA FE over genes.
##################################

#######
#Read .tab-file with TopoA FE over TUs.
#######

def read_FE_tables_combine_together(Path_in, noRif_table_name, Rif_table_name, Ded_table_name, TU_set_name, Expression_data, path_out):
    #Read input dataframes.
    noRif_df=pd.read_csv(Path_in+noRif_table_name, sep='\t')
    Rif_df=pd.read_csv(Path_in+Rif_table_name, sep='\t')
    Ded_df=pd.read_csv(Path_in+Ded_table_name, sep='\t')
    #Combune into new dataframe.
    Data={'Gene_name' : noRif_df['Gene_name'], 'Strand' : noRif_df['Strand'], 'Start' : noRif_df['Start'], 'End' : noRif_df['End'],
          'TopoA_noRif_FE_US' : noRif_df['TopoA_FE_US'], 'TopoA_noRif_FE_GB' : noRif_df['TopoA_FE_GB'], 'TopoA_noRif_FE_DS' : noRif_df['TopoA_FE_DS'], 
          'TopoA_Rif_FE_US' : Rif_df['TopoA_FE_US'], 'TopoA_Rif_FE_GB' : Rif_df['TopoA_FE_GB'], 'TopoA_Rif_FE_DS' : Rif_df['TopoA_FE_DS'],
          'TopoA_Ded_FE_US' : Ded_df['TopoA_FE_US'], 'TopoA_Ded_FE_GB' : Ded_df['TopoA_FE_GB'], 'TopoA_Ded_FE_DS' : Ded_df['TopoA_FE_DS']
          }
    noRif_Rif_Ded_df=pd.DataFrame(Data, columns=['Gene_name', 'Strand', 'Start', 'End', 'TopoA_noRif_FE_US', 'TopoA_noRif_FE_GB', 'TopoA_noRif_FE_DS',
                                                 'TopoA_Rif_FE_US', 'TopoA_Rif_FE_GB', 'TopoA_Rif_FE_DS', 'TopoA_Ded_FE_US', 'TopoA_Ded_FE_GB', 'TopoA_Ded_FE_DS'])
    
    #Read expression data.
    Expression=pd.read_csv(Expression_data, sep='\t')
    #print(Expression)
    #Make some changes over expression data.
    for index, row in Expression.iterrows():
        #print(index, row)
        Genes_ID_changed=row['Gene_name'].rstrip(';')
        Expression_dot=row['Expression'].replace(',', '.')        
        Expression.at[index, 'Gene_name']=Genes_ID_changed
        Expression.at[index, 'Expression']=Expression_dot
    #Combine TopoA FE and expression data.
    Expression=Expression.astype({'Expression': np.float64})
    noRif_Rif_Ded_df_exp=pd.merge(noRif_Rif_Ded_df, Expression)
    #Find and mark tRNA TUs.
    tRNA_rRNA_dict_all={'Gene_name': [], 'tRNA': [], 'rRNA': []}
    tRNA_dict_local={'Gene_name': [], 'tRNA': [], 'TopoA_noRif_FE_GB': [], 'TopoA_Rif_FE_GB': [], 'TopoA_Ded_FE_GB': [], 'Expression': []}
    rRNA_dict_local={'Gene_name': [], 'rRNA': [], 'TopoA_noRif_FE_GB': [], 'TopoA_Rif_FE_GB': [], 'TopoA_Ded_FE_GB': [], 'Expression': []}
    rRNA_operons_genes_list=['16S ribosomal RNA', '23S ribosomal RNA', '5S ribosomal RNA']
    for index, row in noRif_Rif_Ded_df_exp.iterrows():
        tRNA_rRNA_dict_all['Gene_name'].append(row['Gene_name'])
        if (row['Gene_description'].find(' tRNA;')>=0 and row['Gene_description'].find('ribosomal')==-1) or row['Gene_description'].endswith(' tRNA')==True:
            tRNA_rRNA_dict_all['tRNA'].append(1)
            tRNA_dict_local['Gene_name'].append(row['Gene_name'])
            tRNA_dict_local['tRNA'].append(1)
            tRNA_dict_local['TopoA_noRif_FE_GB'].append(row['TopoA_noRif_FE_GB'])
            tRNA_dict_local['TopoA_Rif_FE_GB'].append(row['TopoA_Rif_FE_GB'])
            tRNA_dict_local['TopoA_Ded_FE_GB'].append(row['TopoA_Ded_FE_GB'])
            tRNA_dict_local['Expression'].append(row['Expression'])
        else:
            tRNA_rRNA_dict_all['tRNA'].append(0)
        if (row['Gene_description'].find(' tRNA;')>=0 and row['Gene_description'].find('16S')>=0) or row['Gene_description'] in rRNA_operons_genes_list:
            tRNA_rRNA_dict_all['rRNA'].append(1)
            rRNA_dict_local['Gene_name'].append(row['Gene_name'])
            rRNA_dict_local['rRNA'].append(1)
            rRNA_dict_local['TopoA_noRif_FE_GB'].append(row['TopoA_noRif_FE_GB'])
            rRNA_dict_local['TopoA_Rif_FE_GB'].append(row['TopoA_Rif_FE_GB'])
            rRNA_dict_local['TopoA_Ded_FE_GB'].append(row['TopoA_Ded_FE_GB'])
            rRNA_dict_local['Expression'].append(row['Expression']) 
        else:
            tRNA_rRNA_dict_all['rRNA'].append(0)
    #Combine with the whole dataset
    tRNA_rRNA_all_df=pd.DataFrame(tRNA_rRNA_dict_all, columns=('Gene_name', 'tRNA', 'rRNA'))
    noRif_Rif_Ded_df_exp_trRNA=pd.merge(noRif_Rif_Ded_df_exp, tRNA_rRNA_all_df)
    #Make local datasets for particular group of TUs.
    tRNA_dict_local_df=pd.DataFrame(tRNA_dict_local, columns=('Gene_name', 'tRNA', 'TopoA_noRif_FE_GB', 'TopoA_Rif_FE_GB', 'TopoA_Ded_FE_GB', 'Expression'))
    rRNA_dict_local_df=pd.DataFrame(rRNA_dict_local, columns=('Gene_name', 'rRNA', 'TopoA_noRif_FE_GB', 'TopoA_Rif_FE_GB', 'TopoA_Ded_FE_GB', 'Expression'))
    #Write new dataframe.
    noRif_Rif_Ded_df_exp_trRNA.to_csv(f'{path_out}FE_of_genes\\noRif_Rif_Ded_FE_Exp_tRNA_rRNA_over_{TU_set_name}_15000bp.txt', sep='\t', index=False)    
    return noRif_Rif_Ded_df_exp_trRNA, tRNA_dict_local_df, rRNA_dict_local_df

#######
#Unsorted/Sorted FE over TUs curve plotting.
#######

def plot_FE_ratio(NR_R_D_df, tRNA_data, rRNA_data, TU_set_name, pathout):
    #Scatter plot -Rif_FE vs +Rif_FE and correlation.
    noRif_Rif_corr=scipy.stats.pearsonr(NR_R_D_df['TopoA_noRif_FE_GB'], NR_R_D_df['TopoA_Rif_FE_GB'])
    fit=np.polyfit(NR_R_D_df['TopoA_noRif_FE_GB'], NR_R_D_df['TopoA_Rif_FE_GB'], 1)
    fit_fn=np.poly1d(fit)  
    
    fig=plt.figure(figsize=(8,8), dpi=100) 
    plot1=plt.subplot(111)    
    plot1.scatter(NR_R_D_df['TopoA_noRif_FE_GB'], NR_R_D_df['TopoA_Rif_FE_GB'], c='#74ffa0', alpha=1, linewidths=0.3, s=50, edgecolors='black', label='TopoA FE', zorder=7)
    plot1.scatter(tRNA_data['TopoA_noRif_FE_GB'], tRNA_data['TopoA_Rif_FE_GB'], c='#30006F', alpha=1, linewidths=0.3, s=50, edgecolors='black', label='TopoA FE tRNA', zorder=8)
    plot1.scatter(rRNA_data['TopoA_noRif_FE_GB'], rRNA_data['TopoA_Rif_FE_GB'], c='#E4FC00', alpha=1, linewidths=0.3, s=50, edgecolors='black', label='TopoA FE rRNA', zorder=9)
    plot1.annotate(f'Pearson correlation={round(noRif_Rif_corr[0],3)}\nSignificance={round(float(str(noRif_Rif_corr[1]).split("e")[0]),3)}e{str(noRif_Rif_corr[1]).split("e")[1]}', xy=(0.02, 0.55), xycoords='axes fraction', color='black', size=20) 
    plot1.plot(NR_R_D_df['TopoA_noRif_FE_GB'], fit_fn(NR_R_D_df['TopoA_noRif_FE_GB']), label=f'y={round(fit[0], 3)}x+{round(fit[1], 3)}', color='black', linestyle='--', linewidth=1, alpha=0.8, zorder=10)     
    plot1.set_xlabel('-Rif', size=30, labelpad=8, rotation=0)
    plot1.set_ylabel('+Rif', size=40, labelpad=8, rotation=90)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plot1.legend(fontsize=23) 
    #plot1.set_xticklabels(xticknames1)
    #plt.setp(plot1.set_xticklabels(xticknames1), rotation=0, fontsize=25)
    #plot1.set_yticklabels(yticknames1)
    #plt.setp(plot1.set_yticklabels(yticknames1), rotation=0, fontsize=25)
    plt.tight_layout()
    plt.show()
    fig.savefig(f'{pathout}Figures\Rif_noRif_over_GB_scatter_{TU_set_name}.png', dpi=400, figsize=(8, 8)) 
    plt.close()    


    #-Rif sorted.
    X=NR_R_D_df.index.tolist()
    print(scipy.stats.pearsonr(NR_R_D_df['TopoA_noRif_FE_GB'], NR_R_D_df['TopoA_Rif_FE_GB']))
    NR_R_D_df_sorted=NR_R_D_df.sort_values(by=['TopoA_noRif_FE_GB'])
    Y_noRif=NR_R_D_df_sorted['TopoA_noRif_FE_GB'].tolist()
    Y_Rif=NR_R_D_df_sorted['TopoA_Rif_FE_GB'].tolist()
    Y_Ded=NR_R_D_df_sorted['TopoA_Ded_FE_GB'].tolist()
    NR_R_D_df_sorted['New_index']=X
    #Prepare X_coords for tRNA and rRNA.
    tRNA_data_with_new_index=pd.merge(tRNA_data, NR_R_D_df_sorted)
    X_tRNA=tRNA_data_with_new_index.New_index.tolist()
    rRNA_data_with_new_index=pd.merge(rRNA_data, NR_R_D_df_sorted)
    X_rRNA=rRNA_data_with_new_index.New_index.tolist()
    #Plot data
    fig=plt.figure(figsize=(16,6), dpi=100) 
    plot1=plt.subplot(111)
    ticks1=[1]    
    plot1.set_yticks(ticks1, minor=True)
    plot1.grid(True, axis='y', which='minor', linewidth=2, linestyle='--', alpha=0.9)
    plot1.grid(True, axis='x', which='minor', linewidth=2, linestyle='--', alpha=0.9)
    plot1.scatter(X, Y_noRif, c='#74ffa0', alpha=1, linewidths=0.1, s=50, edgecolors='black', label='-Rif', zorder=8)
    plot1.scatter(X, Y_Ded, c='#8b74ff', alpha=1, linewidths=0.2, s=50, edgecolors='black', label='Ded', zorder=6)
    plot1.scatter(X_tRNA, tRNA_data['TopoA_noRif_FE_GB'], c='#30006F', alpha=1, linewidths=0.3, s=30, edgecolors='black', label='tRNA', zorder=9)
    plot1.scatter(X_rRNA, rRNA_data['TopoA_noRif_FE_GB'], c='#E4FC00', alpha=1, linewidths=0.3, s=30, edgecolors='black', label='rRNA', zorder=10)    
    #plot1.scatter(X, Y_Rif, c='#ff8074', alpha=1, linewidths=0.2, s=50, edgecolors='black', label='+Rif', zorder=9)
    xticknames1=np.arange(0, len(Y_noRif)+1, 400)
    #yticknames1=np.arange(0, max(Y_noRif)+1, 2).tolist()+[1.0]
    yticknames1=np.arange(int(min(Y_noRif+Y_Ded))-1, max(Y_noRif+Y_Ded)+1, 2).tolist()+[1.0]
    plot1.set_yticks(yticknames1, minor=False)
    plot1.set_xticks(xticknames1, minor=False)
    plot1.set_xlabel('Genes', size=33, labelpad=8, rotation=0)
    plot1.set_xticklabels(xticknames1)
    plt.setp(plot1.set_xticklabels(xticknames1), rotation=0, fontsize=25)
    plot1.set_yticklabels(yticknames1)
    plt.setp(plot1.set_yticklabels(yticknames1), rotation=0, fontsize=25)
    plot1.set_ylabel('-Rif FE', size=33, labelpad=8, rotation=90)
    lgnd=plot1.legend(fontsize=20) 
    lgnd.legendHandles[0]._sizes = [250]
    lgnd.legendHandles[1]._sizes = [250]   
    plt.tight_layout()
    plt.show()
    fig.savefig(f'{pathout}Figures\\noRif_Ded_over_GB_scatter_{TU_set_name}_sorted.png', dpi=400, figsize=(10, 5)) 
    plt.close()
    return


#######
#Expression and TopoA FE -Rif, +Rif.
#######

def plot_Exp_FE(NR_R_D_df, tRNA_data, rRNA_data, TU_set_name, pathout):
    #Scatter plot -Rif FE vs Expression and correlation.
    print(NR_R_D_df['TopoA_noRif_FE_GB'], '\n', NR_R_D_df['Expression'])
    FE_Exp_corr=scipy.stats.pearsonr(NR_R_D_df['TopoA_noRif_FE_GB'], NR_R_D_df['Expression'])
    Sig_sci="{:.2e}".format(FE_Exp_corr[1])
    print(Sig_sci)    
    print(FE_Exp_corr)
    
    fig=plt.figure(figsize=(8,8), dpi=100) 
    plot1=plt.subplot(111)    
    plot1.plot(NR_R_D_df['Expression'], NR_R_D_df['TopoA_noRif_FE_GB'], 'o', c='#74ffa0', alpha=1, markersize=10, markeredgecolor='black', markeredgewidth=0.1, label='TopoA +Rif', zorder=7)
    plot1.plot(tRNA_data['Expression'], tRNA_data['TopoA_noRif_FE_GB'], 'o', c='#30006F', alpha=1, markersize=10, markeredgecolor='black', markeredgewidth=0.1, label='TopoA +Rif tRNA', zorder=8)
    plot1.plot(rRNA_data['Expression'], rRNA_data['TopoA_noRif_FE_GB'], 'o', c='#E4FC00', alpha=1, markersize=10, markeredgecolor='black', markeredgewidth=0.1, label='TopoA +Rif rRNA', zorder=9)
    plot1.annotate(f'Pearson correlation={round(FE_Exp_corr[0],3)}\nSignificance={round(float(str(Sig_sci).split("e")[0]),3)}e{str(Sig_sci).split("e")[1]}', xy=(0.02, 0.65), xycoords='axes fraction', color='black', size=20, zorder=12) 
    plot1.set_xscale('log')
    plot1.set_xlabel('Expression', size=30, labelpad=8, rotation=0)
    plot1.set_ylabel('-Rif', size=40, labelpad=8, rotation=90)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=30)
    plot1.legend(fontsize=23, loc='upper left').set_zorder(11)
    plt.tight_layout()
    plt.show()
    fig.savefig(f'{pathout}Figures\\noRif_Exp_over_GB_scatter_{TU_set_name}.png', dpi=400, figsize=(8, 8)) 
    plt.close()   
    
    #Scatter plot +Rif FE vs Expression and correlation.
    FE_Exp_corr=scipy.stats.pearsonr(NR_R_D_df['TopoA_Rif_FE_GB'], NR_R_D_df['Expression'])  
    Sig_sci="{:.2e}".format(FE_Exp_corr[1])
    print(Sig_sci)
    print(FE_Exp_corr)
    
    fig=plt.figure(figsize=(8,8), dpi=100) 
    plot1=plt.subplot(111)    
    plot1.plot(NR_R_D_df['Expression'], NR_R_D_df['TopoA_Rif_FE_GB'], 'o', c='#74ffa0', alpha=1, markersize=10, markeredgecolor='black', markeredgewidth=0.1, label='TopoA +Rif', zorder=7)
    plot1.plot(tRNA_data['Expression'], tRNA_data['TopoA_Rif_FE_GB'], 'o', c='#30006F', alpha=1, markersize=10, markeredgecolor='black', markeredgewidth=0.1, label='TopoA +Rif tRNA', zorder=8)
    plot1.plot(rRNA_data['Expression'], rRNA_data['TopoA_Rif_FE_GB'], 'o', c='#E4FC00', alpha=1, markersize=10, markeredgecolor='black', markeredgewidth=0.1, label='TopoA +Rif rRNA', zorder=9)
    plot1.annotate(f'Pearson correlation={round(FE_Exp_corr[0],3)}\nSignificance={round(float(str(Sig_sci).split("e")[0]),3)}e{str(Sig_sci).split("e")[1]}', xy=(0.02, 0.65), xycoords='axes fraction', color='black', size=20, zorder=12) 
    plot1.set_xlabel('Expression', size=30, labelpad=8, rotation=0)
    plot1.set_ylabel('+Rif', size=40, labelpad=8, rotation=90)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=30)
    plot1.legend(fontsize=23, loc='upper left').set_zorder(11) 
    plot1.set_xscale('log')
    plt.tight_layout()
    plt.show()
    fig.savefig(f'{pathout}Figures\\Rif_Exp_over_GB_scatter_{TU_set_name}.png', dpi=400, figsize=(8, 8)) 
    plt.close()      
    return


#######
#Expression and TopoA Ded FE.
#######

def plot_Exp_Ded_FE(NR_R_D_df, tRNA_data, rRNA_data, TU_set_name, pathout):
    #Scatter plot Ded FE vs Expression and correlation.
    print(NR_R_D_df['TopoA_Ded_FE_GB'], '\n', NR_R_D_df['Expression'])
    FE_Exp_corr=scipy.stats.pearsonr(NR_R_D_df['TopoA_Ded_FE_GB'].abs(), NR_R_D_df['Expression'])
    Sig_sci="{:.2e}".format(FE_Exp_corr[1])
    print(Sig_sci)    
    print(FE_Exp_corr)
    
    fig=plt.figure(figsize=(8,8), dpi=100) 
    plot1=plt.subplot(111)    
    plot1.plot(NR_R_D_df['Expression'], NR_R_D_df['TopoA_Ded_FE_GB'].abs(), 'o', c='#74ffa0', alpha=1, markersize=10, markeredgecolor='black', markeredgewidth=0.1, label='TopoA |-Rif-+Rif|', zorder=7)
    plot1.plot(tRNA_data['Expression'], tRNA_data['TopoA_Ded_FE_GB'].abs(), 'o', c='#30006F', alpha=1, markersize=10, markeredgecolor='black', markeredgewidth=0.1, label='TopoA |-Rif-+Rif| tRNA', zorder=8)
    plot1.plot(rRNA_data['Expression'], rRNA_data['TopoA_Ded_FE_GB'].abs(), 'o', c='#E4FC00', alpha=1, markersize=10, markeredgecolor='black', markeredgewidth=0.1, label='TopoA |-Rif-+Rif| rRNA', zorder=9)
    plot1.annotate(f'Pearson correlation={round(FE_Exp_corr[0],3)}\nSignificance={round(float(str(Sig_sci).split("e")[0]),3)}e{str(Sig_sci).split("e")[1]}', xy=(0.02, 0.65), xycoords='axes fraction', color='black', size=20, zorder=12) 
    plot1.set_xscale('log')
    plot1.set_xlabel('Expression', size=30, labelpad=8, rotation=0)
    plot1.set_ylabel('-Rif-+Rif', size=40, labelpad=8, rotation=90)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=30)
    plot1.legend(fontsize=23, loc='upper left').set_zorder(1)
    plt.tight_layout()
    plt.show()
    fig.savefig(f'{pathout}Figures\\absDed_FE_Exp_over_GB_scatter_{TU_set_name}.png', dpi=400, figsize=(8, 8)) 
    plt.close()   
    
    #-Rif sorted.
    X=NR_R_D_df.index.tolist()
    NR_R_D_df_sorted=NR_R_D_df.sort_values(by=['TopoA_Ded_FE_GB'])
    Y_noRif=NR_R_D_df_sorted['TopoA_noRif_FE_GB'].tolist()
    Y_Rif=NR_R_D_df_sorted['TopoA_Rif_FE_GB'].tolist()
    Y_Ded=NR_R_D_df_sorted['TopoA_Ded_FE_GB'].tolist()
    NR_R_D_df_sorted['New_index']=X
    #Prepare X_coords for tRNA and rRNA.
    tRNA_data_with_new_index=pd.merge(tRNA_data, NR_R_D_df_sorted)
    X_tRNA=tRNA_data_with_new_index.New_index.tolist()
    rRNA_data_with_new_index=pd.merge(rRNA_data, NR_R_D_df_sorted)
    X_rRNA=rRNA_data_with_new_index.New_index.tolist()
    #Plot data
    fig=plt.figure(figsize=(16,6), dpi=100) 
    plot1=plt.subplot(111)
    ticks1=[0]    
    plot1.set_yticks(ticks1, minor=True)
    plot1.grid(True, axis='y', which='minor', linewidth=2, linestyle='--', alpha=0.9)
    plot1.grid(True, axis='x', which='minor', linewidth=2, linestyle='--', alpha=0.9)
    #plot1.scatter(X, Y_noRif, c='#74ffa0', alpha=1, linewidths=0.1, s=50, edgecolors='black', label='-Rif', zorder=8)
    plot1.scatter(X, Y_Ded, c='#74ffa0', alpha=1, linewidths=0.05, s=50, edgecolors='black', label='-Rif-+Rif', zorder=6)
    plot1.scatter(X_tRNA, tRNA_data['TopoA_Ded_FE_GB'], c='#30006F', alpha=1, linewidths=0.1, s=30, edgecolors='black', label='tRNA', zorder=9)
    plot1.scatter(X_rRNA, rRNA_data['TopoA_Ded_FE_GB'], c='#E4FC00', alpha=1, linewidths=0.1, s=30, edgecolors='black', label='rRNA', zorder=10)    
    #plot1.scatter(X, Y_Rif, c='#ff8074', alpha=1, linewidths=0.2, s=50, edgecolors='black', label='+Rif', zorder=9)
    xticknames1=np.arange(0, len(Y_Ded)+1, 400)
    yticknames1=np.arange(int(min(Y_Ded))-1, max(Y_Ded)+1, 2).tolist()+[0.0]
    plot1.set_yticks(yticknames1, minor=False)
    plot1.set_xticks(xticknames1, minor=False)
    plot1.set_xlabel('TUs', size=33, labelpad=8, rotation=0)
    plot1.set_xticklabels(xticknames1)
    plt.setp(plot1.set_xticklabels(xticknames1), rotation=0, fontsize=25)
    plot1.set_yticklabels(yticknames1)
    plt.setp(plot1.set_yticklabels(yticknames1), rotation=0, fontsize=25)
    plot1.set_ylabel('-Rif-+Rif', size=33, labelpad=8, rotation=90)
    lgnd=plot1.legend(fontsize=30, loc='lower right')
    lgnd.legendHandles[0]._sizes = [250] 
    lgnd.legendHandles[1]._sizes = [100] 
    lgnd.legendHandles[2]._sizes = [100]
    plt.tight_layout()
    plt.show()
    #fig.savefig(f'{pathout}Figures\\Ded_over_GB_scatter_{TU_set_name}_sorted.png', dpi=400, figsize=(10, 5)) 
    plt.close()
    return



#######
#Wrapper: reads and merges dataFrames with -Rif, +Rif, Ded data;
#Plot -Rif vs +Rif, computes correlation, sorts -Rif, plots again.
#######


def wrapper_compare_noRif_Rif_merge_plot(inpath, noRif_file, Rif_file, Ded_file, TU_set_name, Expression_data, outpath):
    noRif_Rif_Ded_data=read_FE_tables_combine_together(inpath, noRif_file, Rif_file, Ded_file, TU_set_name, Expression_data, outpath)
    #plot_FE_ratio(noRif_Rif_Ded_data[0], noRif_Rif_Ded_data[1], noRif_Rif_Ded_data[2], TU_set_name, outpath)
    #plot_Exp_FE(noRif_Rif_Ded_data[0], noRif_Rif_Ded_data[1], noRif_Rif_Ded_data[2], TU_set_name, outpath)
    plot_Exp_Ded_FE(noRif_Rif_Ded_data[0], noRif_Rif_Ded_data[1], noRif_Rif_Ded_data[2], TU_set_name, outpath)
    return

#wrapper_compare_noRif_Rif_merge_plot('F:\TopoI_ChIP-Seq\Ec_TopoI_data\FE_of_genes\\', 'FE_Rep12_FE_over_All_operonss_15000bp.txt', 
#                                    'FE_Rif_Rep12_FE_over_All_operonss_15000bp.txt', 'FE_Ded_Rep12_Rif12_FE_over_All_operonss_15000bp.txt', 
#                                    'All_operons', 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_operons_expression.txt', 
#                                    Output_path)

#wrapper_compare_noRif_Rif_merge_plot('F:\TopoI_ChIP-Seq\Ec_TopoI_data\FE_of_genes\\', 'FE_Rep12_FE_over_All_geness_15000bp.txt', 
#                                    'FE_Rif_Rep12_FE_over_All_geness_15000bp.txt', 'FE_Ded_Rep12_Rif12_FE_over_All_geness_15000bp.txt', 
#                                    'All_genes', 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_genes_expression.txt',
#                                    Output_path)




#######
#Average replicas of ChIP-Seq experiments.
#######

Replicas_input_dict={'Replic_1' : 'F:\E_coli_ChIPs\Cov_depth_nodup\SRR1946726_name_sorted_fm_ps_nd.wig',
                     'Replic_2' : 'F:\E_coli_ChIPs\Cov_depth_nodup\SRR1946727_name_sorted_fm_ps_nd.wig',
                     'Replic_3' : 'F:\E_coli_ChIPs\Cov_depth_nodup\SRR1946728_name_sorted_fm_ps_nd.wig'
                     }

def read_replicas_average_write(replicas_input_dict, ID, out_path):
    #Read data
    data_dict={}
    for name, inpath in replicas_input_dict.items():
        wig_data=wig_parsing(inpath)
        data_dict[name]=wig_data
    #Compute and plot correlation matrix.    
    wig_dataframe=pd.DataFrame(data_dict)    
    correlation_matrix(wig_dataframe, 'pearson', f'Correlation matrix between replicas of {ID} ChIP-Seq', f'{out_path}\Figures\Correlation_between_replicas_of_{ID}.png')
    #Average replicas.
    wig_av=[]
    for i in range(len(data_dict[name])):
        position_list=[]
        for replic_name, replic_data in data_dict.items():
            position_list.append(replic_data[i])
        position_av=np.mean(position_list)
        wig_av.append(position_av)
    #Write final wig track.      
    write_wig(wig_av, f'{out_path}\Cov_depth_nodup\{ID}_cov.wig', ID)
    return

#read_replicas_average_write(Replicas_input_dict, 'MukB_Nolivos', 'F:\E_coli_ChIPs')


#######
#Add new data for correlation analysis: RpoB, H-NS, Fis, GC%, etc.
#######

Special_data={'PolSofi' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\Pol_Sofi_LB_w3110_for_Mu.wig',
              'RpoB' : 'F:\E_coli_ChIPs\Cov_depth\RpoB_ME_Kahramanoglou.wig',
              'HNS' : 'F:\E_coli_ChIPs\Cov_depth\HNS_ME_Kahramanoglou.wig',
              'GC' : 'C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\E_coli_w3110_Mu_GC_133bp.wig',
              'Fis' : 'F:\E_coli_ChIPs\Cov_depth\Fis_ME_Kahramanoglou.wig',
              'MukB' : 'F:\E_coli_ChIPs\Cov_depth_nodup\MukB_Nolivos_cov.wig',
              'MatP' : 'F:\E_coli_ChIPs\Cov_depth_nodup\MatP_Nolivos_cov.wig',
              'TopoA -Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_average_Rep1_FE_Rep2_FE.wig',
              'TopoA +Rif' : 'F:\TopoI_ChIP-Seq\Ec_TopoI_data\Fold_enrichment\TopoA_average_Rif_Rep1_FE_Rif_Rep2_FE.wig',
              'Gyrase Cfx' : 'F:\Gyrase_Topo-Seq_data\Cfx_10mkM_average.wig',
              'TopoIV Cfx' : 'F:\TopoIV_Topo-Seq\Fold_enrichment\Cfx_average.wig',
              'RpoS' : 'F:\E_coli_ChIPs\Cov_depth\RpoS_Peano_av.wig',
              }

def read_correlate_data(replicas_input_dict, out_path):
    #Read data
    data_dict={}
    for name, inpath in replicas_input_dict.items():
        wig_data=wig_parsing(inpath)
        data_dict[name]=wig_data
    #Compute and plot correlation matrix.    
    wig_dataframe=pd.DataFrame(data_dict)    
    correlation_matrix(wig_dataframe, 'pearson', f'Correlation matrix between ChIP-Seq datasets', f'{out_path}\Figures\Correlation_between_E_coli_ChIP-Seq_datasets.png')
    #Perform hierarchy clustering and plot again.
    wig_dataframe_clusterized=Clustering(wig_dataframe)
    correlation_matrix(wig_dataframe_clusterized, 'pearson', f'Clusterized correlation matrix between ChIP-Seq datasets', f'{out_path}\Figures\Clusterized_correlation_between_E_coli_ChIP-Seq_datasets.png')
    return

read_correlate_data(Special_data, 'F:\TopoI_ChIP-Seq\Ec_TopoI_data')