###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes WIG files with information about the cumulative distribution of some 
#protein over transcription units (TUs). Plots this information.
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


#Path to the working directory.
PWD='C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs'
#Name of the signal to plotted (protein or smth.).
Signal_name='DRIP-Seq_combined_pictures'
#Half-window width will be used to smooth signal.
Sm_window=100
#Dictionary of pathes to input WIG data.
Wig_data_in_dict_operons={'All operons F'  : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_operons_no_rRNA_lacI_cat\\Signal_DRIP_CTD_minus_Rif_minus_fw_over_All_operons_no_rRNA_lacI_cat_width_15000bp_gb_5000bp.wig',
                          'All operons R'  : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_operons_no_rRNA_lacI_cat\\Signal_DRIP_CTD_minus_Rif_minus_rv_over_All_operons_no_rRNA_lacI_cat_width_15000bp_gb_5000bp.wig',
                          'All operons LF' : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_operons_no_rRNA_lacI_cat_L\\Signal_DRIP_CTD_minus_Rif_minus_fw_over_All_operons_no_rRNA_lacI_cat_L_width_15000bp_gb_5000bp.wig',
                          'All operons LR' : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_operons_no_rRNA_lacI_cat_L\\Signal_DRIP_CTD_minus_Rif_minus_rv_over_All_operons_no_rRNA_lacI_cat_L_width_15000bp_gb_5000bp.wig',
                          'All operons RF' : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_operons_no_rRNA_lacI_cat_R\\Signal_DRIP_CTD_minus_Rif_minus_fw_over_All_operons_no_rRNA_lacI_cat_R_width_15000bp_gb_5000bp.wig',
                          'All operons RR' : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_operons_no_rRNA_lacI_cat_R\\Signal_DRIP_CTD_minus_Rif_minus_rv_over_All_operons_no_rRNA_lacI_cat_R_width_15000bp_gb_5000bp.wig',
                          }
Wig_data_in_dict_genes={'All genes F'  : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_genes_no_rRNA_lacI_cat\\Signal_DRIP_CTD_minus_Rif_minus_fw_over_All_genes_no_rRNA_lacI_cat_width_15000bp_gb_5000bp.wig',
                        'All genes R'  : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_genes_no_rRNA_lacI_cat\\Signal_DRIP_CTD_minus_Rif_minus_rv_over_All_genes_no_rRNA_lacI_cat_width_15000bp_gb_5000bp.wig',
                        'All genes LF' : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_genes_no_rRNA_lacI_cat_L\\Signal_DRIP_CTD_minus_Rif_minus_fw_over_All_genes_no_rRNA_lacI_cat_L_width_15000bp_gb_5000bp.wig',
                        'All genes LR' : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_genes_no_rRNA_lacI_cat_L\\Signal_DRIP_CTD_minus_Rif_minus_rv_over_All_genes_no_rRNA_lacI_cat_L_width_15000bp_gb_5000bp.wig',
                        'All genes RF' : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_genes_no_rRNA_lacI_cat_R\\Signal_DRIP_CTD_minus_Rif_minus_fw_over_All_genes_no_rRNA_lacI_cat_R_width_15000bp_gb_5000bp.wig',
                        'All genes RR' : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Signal_over_TUs\Signal_wig\All_genes_no_rRNA_lacI_cat_R\\Signal_DRIP_CTD_minus_Rif_minus_rv_over_All_genes_no_rRNA_lacI_cat_R_width_15000bp_gb_5000bp.wig',
                        }
Wig_data_in_dict_transcripts={'All TUs CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\All_TUs_2173\Signal_TopA_CTD-Rif-_av_1_2_3_over_All_TUs_2173_width_15000bp_gb_5000bp.wig', 
                              'All TUs CTD-Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\All_TUs_2173\Signal_TopA_CTD-Rif+_av_1_2_3_over_All_TUs_2173_width_15000bp_gb_5000bp.wig', 
                              'All TUs CTD+Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\All_TUs_2173\Signal_TopA_CTD+Rif-_av_1_2_3_over_All_TUs_2173_width_15000bp_gb_5000bp.wig', 
                              'All TUs CTD+Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\All_TUs_2173\Signal_TopA_CTD+Rif+_av_2_3_over_All_TUs_2173_width_15000bp_gb_5000bp.wig',                               
                              'All TUs no tRNA, rRNA CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\All_TUs_no_tRNA_rRNA\Signal_TopA_CTD-Rif-_av_1_2_3_over_All_TUs_no_tRNA_rRNA_width_15000bp_gb_5000bp.wig', 
                              'All TUs no tRNA, rRNA CTD-Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\All_TUs_no_tRNA_rRNA\Signal_TopA_CTD-Rif+_av_1_2_3_over_All_TUs_no_tRNA_rRNA_width_15000bp_gb_5000bp.wig', 
                              'All TUs no tRNA, rRNA CTD+Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\All_TUs_no_tRNA_rRNA\Signal_TopA_CTD+Rif-_av_1_2_3_over_All_TUs_no_tRNA_rRNA_width_15000bp_gb_5000bp.wig', 
                              'All TUs no tRNA, rRNA CTD+Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\All_TUs_no_tRNA_rRNA\Signal_TopA_CTD+Rif+_av_2_3_over_All_TUs_no_tRNA_rRNA_width_15000bp_gb_5000bp.wig',                                                              
                              'HETU 323 CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HETU_323\Signal_TopA_CTD-Rif-_av_1_2_3_over_HETU_323_width_15000bp_gb_5000bp.wig',
                              'HETU 323 CTD-Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HETU_323\Signal_TopA_CTD-Rif+_av_1_2_3_over_HETU_323_width_15000bp_gb_5000bp.wig',
                              'HETU 323 CTD+Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HETU_323\Signal_TopA_CTD+Rif-_av_1_2_3_over_HETU_323_width_15000bp_gb_5000bp.wig',
                              'HETU 323 CTD+Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HETU_323\Signal_TopA_CTD+Rif+_av_2_3_over_HETU_323_width_15000bp_gb_5000bp.wig',                              
                              'HETU 321 no ompX CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_TopA_CTD-Rif-_av_1_2_3_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',
                              'HETU 321 no ompX CTD-Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_TopA_CTD-Rif+_av_1_2_3_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',    
                              'HETU 321 no ompX CTD+Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_TopA_CTD+Rif-_av_1_2_3_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',
                              'HETU 321 no ompX CTD+Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_TopA_CTD+Rif+_av_2_3_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',                                  
                              'LETU 245 CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LETU_245\Signal_TopA_CTD-Rif-_av_1_2_3_over_LETU_245_width_15000bp_gb_5000bp.wig',
                              'LETU 245 CTD-Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LETU_245\Signal_TopA_CTD-Rif+_av_1_2_3_over_LETU_245_width_15000bp_gb_5000bp.wig',
                              'LETU 245 CTD+Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LETU_245\Signal_TopA_CTD+Rif-_av_1_2_3_over_LETU_245_width_15000bp_gb_5000bp.wig',
                              'LETU 245 CTD+Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LETU_245\Signal_TopA_CTD+Rif+_av_2_3_over_LETU_245_width_15000bp_gb_5000bp.wig',                              
                              'LETU 244 no ybiI CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_TopA_CTD-Rif-_av_1_2_3_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig',
                              'LETU 244 no ybiI CTD-Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_TopA_CTD-Rif+_av_1_2_3_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig',
                              'LETU 244 no ybiI CTD+Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_TopA_CTD+Rif-_av_1_2_3_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig',
                              'LETU 244 no ybiI CTD+Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_TopA_CTD+Rif+_av_2_3_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig',                              
                              'LETU 248 no ybiI no appY CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LE_TU_248_no_ybiI_no_appY\Signal_TopA_CTD-Rif-_av_1_2_3_over_LE_TU_248_no_ybiI_no_appY_width_15000bp_gb_5000bp.wig',
                              'LETU 248 no ybiI no appY CTD-Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LE_TU_248_no_ybiI_no_appY\Signal_TopA_CTD-Rif+_av_1_2_3_over_LE_TU_248_no_ybiI_no_appY_width_15000bp_gb_5000bp.wig',
                              'LETU 248 no ybiI no appY CTD+Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LE_TU_248_no_ybiI_no_appY\Signal_TopA_CTD+Rif-_av_1_2_3_over_LE_TU_248_no_ybiI_no_appY_width_15000bp_gb_5000bp.wig',
                              'LETU 248 no ybiI no appY CTD+Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\LE_TU_248_no_ybiI_no_appY\Signal_TopA_CTD+Rif+_av_2_3_over_LE_TU_248_no_ybiI_no_appY_width_15000bp_gb_5000bp.wig',                                                        
                              'rRNA 7 CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\\rRNA_7\Signal_TopA_CTD-Rif-_av_1_2_3_over_rRNA_7_width_15000bp_gb_5000bp.wig',
                              'rRNA 7 CTD-Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\\rRNA_7\Signal_TopA_CTD-Rif+_av_1_2_3_over_rRNA_7_width_15000bp_gb_5000bp.wig', 
                              'rRNA 7 CTD+Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\\rRNA_7\Signal_TopA_CTD+Rif-_av_1_2_3_over_rRNA_7_width_15000bp_gb_5000bp.wig',
                              'rRNA 7 CTD+Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\\rRNA_7\Signal_TopA_CTD+Rif+_av_2_3_over_rRNA_7_width_15000bp_gb_5000bp.wig',                               
                              'tRNA 49 CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\\tRNA_49\Signal_TopA_CTD-Rif-_av_1_2_3_over_tRNA_49_width_15000bp_gb_5000bp.wig',
                              'tRNA 49 CTD-Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\\tRNA_49\Signal_TopA_CTD-Rif+_av_1_2_3_over_tRNA_49_width_15000bp_gb_5000bp.wig',
                              'tRNA 49 CTD+Rif-' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\\tRNA_49\Signal_TopA_CTD+Rif-_av_1_2_3_over_tRNA_49_width_15000bp_gb_5000bp.wig',
                              'tRNA 49 CTD+Rif+' : 'C:\Sutor\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\\tRNA_49\Signal_TopA_CTD+Rif+_av_2_3_over_tRNA_49_width_15000bp_gb_5000bp.wig'                              
                              }

'''
'HETU no tRNA, rRNA 317 CTD-Rif-' : 'D:\Sutormin_data\D_Sutormin\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HE_TU_317_no_tRNA_rRNA_ompX\Signal_TopoA -Rif_over_HE_TU_317_no_tRNA_rRNA_ompX_width_15000bp_gb_5000bp.wig', 
'HETU no tRNA, rRNA 317 CTD-Rif+' : 'D:\Sutormin_data\D_Sutormin\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HE_TU_317_no_tRNA_rRNA_ompX\Signal_TopoA +Rif_over_HE_TU_317_no_tRNA_rRNA_ompX_width_15000bp_gb_5000bp.wig', 
'HETU no tRNA, rRNA 317 CTD+Rif-' : 'D:\Sutormin_data\D_Sutormin\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HE_TU_317_no_tRNA_rRNA_ompX\Signal_TopoA -Rif_over_HE_TU_317_no_tRNA_rRNA_ompX_width_15000bp_gb_5000bp.wig', 
'HETU no tRNA, rRNA 317 CTD+Rif+' : 'D:\Sutormin_data\D_Sutormin\Science\Signal_over_TUs\All_combinations_for_TopA\Signal_of_TUs_wig\HE_TU_317_no_tRNA_rRNA_ompX\Signal_TopoA +Rif_over_HE_TU_317_no_tRNA_rRNA_ompX_width_15000bp_gb_5000bp.wig',
'''                              

Wig_data_in_dict_transcripts_RNApol={'All TUs CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_2173\Signal_PolSofi_over_All_TU_2173_width_15000bp_gb_5000bp.wig',         
                                     'HETU 321 no ompX CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_PolSofi_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',               
                                     'LETU 248 no ybiI no appY CTD-Rif-' : 'C:\Sutor\Science\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_PolSofi_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig'
                                     }

Wig_data_in_dict_transcripts_gyrase={'All TUs' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_2173\Signal_Gyrase Cfx_over_All_TU_2173_width_15000bp_gb_5000bp.wig', 
                                     'All TUs Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_2173\Signal_Gyrase Cfx +Rif_over_All_TU_2173_width_15000bp_gb_5000bp.wig', 
                                     'All TUs no tRNA, rRNA' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_249_no_tRNA_rRNA\Signal_Gyrase Cfx_over_All_TU_249_no_tRNA_rRNA_width_15000bp_gb_5000bp.wig', 
                                     'All TUs no tRNA, rRNA Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\All_TU_249_no_tRNA_rRNA\Signal_Gyrase Cfx +Rif_over_All_TU_249_no_tRNA_rRNA_width_15000bp_gb_5000bp.wig',  
                                     'HETU 321 no ompX' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_Gyrase Cfx_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',
                                     'HETU 321 no ompX Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\HE_TU_321_no_ompX\Signal_Gyrase Cfx +Rif_over_HE_TU_321_no_ompX_width_15000bp_gb_5000bp.wig',                              
                                     'LETU 244 no ybiI' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_Gyrase Cfx_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig',
                                     'LETU 244 no ybiI Rif' : 'F:\Signal_over_TUs\Transcript-based\Signal_of_TUs_wig\LE_TU_244_no_ybiI\Signal_Gyrase Cfx +Rif_over_LE_TU_244_no_ybiI_width_15000bp_gb_5000bp.wig',                           
                                     }
#Set type to choose plotting parameters: genes, operons or transcripts.
Set_type="operons"

#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
Out_path=f'{PWD}\Plots_together\{Signal_name}'
Dir_check_create(Out_path)


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
#Plot the signal for all groups of genes together.
#######


def plot_FE_all_expression_gg_Rif_no_Rif(wig_in_dict, sm_window, output_path, set_name, set_type):
    #Number of genes within sets.
    TU_sets_v={'All genes' : 4119, 'HEG 270' : 269, 'LEG 270' : 270, 'HEG 370' : 369, 'LEG 370' : 370,
               'All operons' : 2327, 'HEO 144' : 143, 'LEO 144' : 144, 'HEO 186' : 185, 'LEO 186' : 186, 'LAO 27' : 27, 'SAO 27' : 27, 'rRNA_operons' : 7, 
               'All TUs CTD-Rif-' : 1672, 'All TUs CTD-Rif+' : 1672, 'All TUs CTD+Rif-' : 1672, 'All TUs CTD+Rif+' : 1672,
               'All TUs no tRNA, rRNA CTD-Rif-' : 1634, 'All TUs no tRNA, rRNA CTD-Rif+' : 1634, 'All TUs no tRNA, rRNA CTD+Rif-' : 1634, 'All TUs no tRNA, rRNA CTD+Rif+' : 1634, 
               'HETU no tRNA, rRNA 317' : 202, 'LETU 249' : 205, 
               'HETU 323 CTD-Rif-' : 201, 'HETU 323 CTD-Rif+' : 201, 'HETU 323 CTD+Rif-' : 201, 'HETU 323 CTD+Rif+' : 201, 
               'HETU 321 no ompX CTD-Rif-' : 200, 'HETU 321 no ompX CTD-Rif+' : 200, 'HETU 321 no ompX CTD+Rif-' : 200, 'HETU 321 no ompX CTD+Rif+' : 200, 
               'LETU 245 CTD-Rif-' : 202, 'LETU 245 CTD-Rif+' : 202, 'LETU 245 CTD+Rif-' : 202, 'LETU 245 CTD+Rif+' : 202, 
               'LETU 244 no ybiI CTD-Rif-' : 201, 'LETU 244 no ybiI CTD-Rif+' : 201, 'LETU 244 no ybiI CTD+Rif-' : 201, 'LETU 244 no ybiI CTD+Rif+' : 201,
               'LETU 248 no ybiI no appY CTD-Rif-' : 200, 'LETU 248 no ybiI no appY CTD-Rif+' : 200, 'LETU 248 no ybiI no appY CTD+Rif-' : 200, 'LETU 248 no ybiI no appY CTD+Rif+' : 200, 
               'rRNA 7 CTD-Rif-' : 7, 'rRNA 7 CTD-Rif+' : 7, 'rRNA 7 CTD+Rif-' : 7, 'rRNA 7 CTD+Rif+' : 7, 
               'tRNA 49 CTD-Rif-' : 39, 'tRNA 49 CTD-Rif+' : 39, 'tRNA 49 CTD+Rif-' : 39,'tRNA 49 CTD+Rif+' : 39,
               'All operons F' : 2321, 'All operons R' : 2321, 'All operons LF' : 1091, 'All operons LR' : 1091, 'All operons RF' : 1230, 'All operons RR' : 1230,
               'All genes F' : 4088, 'All genes R' : 4088, 'All genes LF' : 1928, 'All genes LR' : 1928, 'All genes RF' : 2159, 'All genes RR' : 2159}
    #FE averaged WIG parsing
    dict_of_wigs={}
    win_width=15000
    length=5000
    for name, file in wig_in_dict.items():
        print(name, file)
        data=wig_FE_over_genes_parsing(file)
        win_width=data[1]
        length=data[2]
        dict_of_wigs[name]=data[0]        
    positions=np.arange(-win_width, win_width+length, 1)    
    print(win_width, length)
       
    
    #Plot FE over genes.
    plt.figure(figsize=(10, 6), dpi=100) #Def size - 10, 6; Closer look - 3, 6
    plot1=plt.subplot(111)
    ##Genes below
    if set_type=='genes':
        #LEG_270
        plot1.plot(positions, dict_of_wigs['All genes LF'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'All genes LF ({TU_sets_v["All genes LF"]})', zorder=6)
        plot1.plot(positions, dict_of_wigs['All genes LR'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'All genes LR ({TU_sets_v["All genes LR"]})', zorder=5)   
        #LEG_370
        #plot1.plot(positions, dict_of_wigs['LEG_370'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEG ({TU_sets_v["LEG_370"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
        #All_genes
        plot1.plot(positions, dict_of_wigs['All genes F'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All genes F ({TU_sets_v["All genes F"]})', zorder=10)
        plot1.plot(positions, dict_of_wigs['All genes R'], linestyle='--', color='#333738', linewidth=1.5, alpha=0.8, label=f'All genes R ({TU_sets_v["All genes R"]})', zorder=9)
        #HEG_270
        plot1.plot(positions, dict_of_wigs['All genes RF'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'All genes RF ({TU_sets_v["All genes RF"]})', zorder=8)
        plot1.plot(positions, dict_of_wigs['All genes RR'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'All genes RR ({TU_sets_v["All genes RR"]})', zorder=7)
        #HEG_370
        #plot1.plot(positions, dict_of_wigs['HEG_370'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEG ({TU_sets_v["HEG_370"]})', zorder=4)
        #plot1.plot(positions, dict_of_wigs['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)
    ##Operons below.
    elif set_type=="operons":
        #LEO_144
        #plot1.plot(positions, dict_of_wigs['LEO 144'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEO ({TU_sets_v["LEO 144"]})', zorder=6)
        #plot1.plot(positions, dict_of_wigs['LEO_144_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEO Rif ({TU_sets_v["LEO_144_Rif"]})', zorder=5)   
        #LEO_186 or SAO_27
        #plot1.plot(positions, dict_of_wigs['SAO 27'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'SAO ({TU_sets_v["SAO 27"]})', zorder=6)
        plot1.plot(positions, dict_of_wigs['All operons LF'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1.5, label=f'All operons LF ({TU_sets_v["All operons LF"]})', zorder=6)
        plot1.plot(positions, dict_of_wigs['All operons LR'], linestyle='--', color='#757d8b', linewidth=1, alpha=0.8, label=f'All operons LR ({TU_sets_v["All operons LR"]})', zorder=5)        
        #All_operons
        plot1.plot(positions, dict_of_wigs['All operons F'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All operons F ({TU_sets_v["All operons F"]})', zorder=10)
        plot1.plot(positions, dict_of_wigs['All operons R'], linestyle='--', color='#333738', linewidth=1.5, alpha=0.8, label=f'All operons R ({TU_sets_v["All operons R"]})', zorder=9)
        #HEO_144
        #plot1.plot(positions, dict_of_wigs['HEO 144'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEO no rRNA, tRNA ({TU_sets_v["HEO 144"]})', zorder=8)
        #plot1.plot(positions, dict_of_wigs['HEO_no_tRNA_rRNA_144_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEO no rRNA, tRNA Rif ({TU_sets_v["HEO_no_tRNA_rRNA_144_Rif"]})', zorder=7)
        #HEO_186 or rRNA operons or LAO_27
        #plot1.plot(positions, dict_of_wigs['LAO 27'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'LAO ({TU_sets_v["LAO 27"]})', zorder=4)
        #plot1.plot(positions, dict_of_wigs['rRNA_operons'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'rRNA operons ({TU_sets_v["rRNA_operons"]})', zorder=4)
        plot1.plot(positions, dict_of_wigs['All operons RF'], linestyle='-', color='#b08642', linewidth=1.5, alpha=1.5, label=f'All genes RF ({TU_sets_v["All genes RF"]})', zorder=4)
        plot1.plot(positions, dict_of_wigs['All operons RR'], linestyle='--', color='#b08642', linewidth=1, alpha=0.8, label=f'All genes RR ({TU_sets_v["All genes RR"]})', zorder=3)  
        ##Transcription units below.
    elif set_type=="transcripts":
        #Standard order of colors: #757d8b, #333738, #b08642, #e4d1b4.
        #LETU
        #plot1.plot(positions, np.array(dict_of_wigs['LETU 245 CTD-Rif-'])-0.07, linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LETU CTD-Rif- ({TU_sets_v["LETU 245 CTD-Rif-"]})', zorder=6) #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions, np.array(dict_of_wigs['LETU 245 CTD-Rif+'])-0.45, linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'LETU CTD-Rif+ ({TU_sets_v["LETU 245 CTD-Rif-"]})', zorder=5) #R123 -0.1; R12 -0; R3 -0.45     
        #plot1.plot(positions, np.array(dict_of_wigs['LETU 245 CTD+Rif-'])-0.15, linestyle='-', color='#b08642', linewidth=1.5, alpha=1, label=f'LETU CTD+Rif- ({TU_sets_v["LETU 245 CTD-Rif-"]})', zorder=6) #R123 -0.1; R123 -0.15; R123 -0.15
        #plot1.plot(positions, np.array(dict_of_wigs['LETU 245 CTD+Rif+'])-0.25, linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'LETU CTD+Rif+ ({TU_sets_v["LETU 245 CTD-Rif-"]})', zorder=5) #R23 -0.2; R23 -0.25; R23 -0.25             
        #LETU, no ybiI
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU 244 no ybiI CTD-Rif-'])-0.05, linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LETU CTD-/Rif- ({TU_sets_v["LETU 244 no ybiI CTD-Rif-"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU 244 no ybiI CTD-Rif+'])-0.1, linestyle='--', color='#757d8b', linewidth=1, alpha=1.5, label=f'LETU CTD-/Rif+ ({TU_sets_v["LETU 244 no ybiI CTD-Rif-"]})', zorder=5) #Def linewidth=1; #R123 -0.1; R12 +0; R3 -0.45    
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU 244 no ybiI CTD+Rif-'])-0.25, linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LETU CTD+Rif- ({TU_sets_v["LETU 244 no ybiI CTD-Rif-"]})', zorder=6) #Def linewidth=1.5; #R123 -0.1; R123 +0.15; R123 -0.15
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU 244 no ybiI CTD+Rif+'])-0.25, linestyle='-.', color='#757d8b', linewidth=1, alpha=1, label=f'LETU CTD+Rif+ ({TU_sets_v["LETU 244 no ybiI CTD-Rif-"]})', zorder=5) #Def linewidth=1; #R23 -0.2; R23 +0.25; R23 -0.25           
        #LETU, no ybiI, no appY
        plot1.plot(positions,  np.array(dict_of_wigs['LETU 248 no ybiI no appY CTD-Rif-'])-0.05, linestyle='-', color='#B6B8BD', linewidth=3, alpha=1, label=f'LETU CTD-/Rif- ({TU_sets_v["LETU 248 no ybiI no appY CTD-Rif-"]})', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU 248 no ybiI no appY CTD-Rif+'])-0.12, linestyle='--', color='#B6B8BD', linewidth=5, alpha=1.5, label=f'LETU CTD-/Rif+ ({TU_sets_v["LETU 248 no ybiI no appY CTD-Rif-"]})', zorder=5) #Def linewidth=0.8; #R123 -0.1; R12 +0; R3 -0.45    
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU 248 no ybiI no appY CTD+Rif-'])-0.25, linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LETU CTD+Rif- ({TU_sets_v["LETU 244 no ybiI no appY CTD-Rif-"]})', zorder=6) #Def linewidth=1.5; #R123 -0.1; R123 +0.15; R123 -0.15
        #plot1.plot(positions,  np.array(dict_of_wigs['LETU 248 no ybiI no appY CTD+Rif+'])-0.25, linestyle='-.', color='#757d8b', linewidth=1, alpha=1, label=f'LETU CTD+Rif+ ({TU_sets_v["LETU 244 no ybiI no appY CTD-Rif-"]})', zorder=5) #Def linewidth=1; #R23 -0.2; R23 +0.25; R23 -0.25         
        #All_TUs
        plot1.plot(positions, np.array(dict_of_wigs['All TUs CTD-Rif-'])-0.05, linestyle='-', color='#333738', linewidth=3, alpha=2, label=f'All TUs CTD-/Rif- ({TU_sets_v["All TUs CTD-Rif-"]})', zorder=10) #Def linewidth=2; #R123 -0.07; R12 -0.01; R3 -0.17
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs CTD-Rif+'])-0.12, linestyle='--', color='#333738', linewidth=5, alpha=0.8, label=f'All TUs CTD-/Rif+ ({TU_sets_v["All TUs CTD-Rif-"]})', zorder=9) #Def linewidth=0.8; #R123 -0.15; R12 +0; R3 -0.42
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs CTD+Rif-'])-0.27, linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'All TUs CTD+Rif- ({TU_sets_v["All TUs CTD-Rif-"]})', zorder=10) #Def linewidth=2; #R123 -0.25; R123 -0.25; R123 -0.27
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs CTD+Rif+'])-0.24, linestyle='--', color='#e4d1b4', linewidth=1, alpha=0.8, label=f'All TUs CTD+Rif+ ({TU_sets_v["All TUs CTD-Rif-"]})', zorder=9) #Def linewidth=1; #R23 -0.25; R23 -0.22; R23 -0.24
        #All_TUs no tRNA, rRNA.
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no tRNA, rRNA CTD-Rif-'])-0.05, linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All TUs no tRNA, rRNA CTD-/Rif- ({TU_sets_v["All TUs no tRNA, rRNA CTD-Rif-"]})', zorder=10)  #R123 -0.07; R12 -0.02; R3 -0.17
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no tRNA, rRNA CTD-Rif+'])-0.1, linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All TUs no tRNA, rRNA CTD-/Rif+ ({TU_sets_v["All TUs no tRNA, rRNA CTD-Rif-"]})', zorder=9)  #R123 -0.12; R12 -0.02; R3 -0.42
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no tRNA, rRNA CTD+Rif-'])-0.25, linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All TUs no tRNA, rRNA CTD+Rif- ({TU_sets_v["All TUs no tRNA, rRNA CTD-Rif-"]})', zorder=10)  #R123 -0.25; R123 -0.27; R123 -0.27
        #plot1.plot(positions, np.array(dict_of_wigs['All TUs no tRNA, rRNA CTD+Rif+'])-0.25, linestyle='-.', color='#333738', linewidth=1, alpha=0.8, label=f'All TUs no tRNA, rRNA CTD+Rif+ ({TU_sets_v["All TUs no tRNA, rRNA CTD-Rif-"]})', zorder=9)   #R23 -0.22; R23 -0.22; R23 -0.24   
        #HETU
        #plot1.plot(positions, np.array(dict_of_wigs['HETU 323 CTD-Rif-'])-0.25, linestyle='-', color='#757d8b', linewidth=1.5, alpha=0.8, label=f'HETU CTD-Rif- ({TU_sets_v["HETU 323 CTD-Rif-"]})', zorder=8) #R123 -0.05; R12 -0.05; R3 -0.25
        #plot1.plot(positions, np.array(dict_of_wigs['HETU 323 CTD-Rif+'])-0.32, linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'HETU CTD-Rif+ ({TU_sets_v["HETU 323 CTD-Rif-"]})', zorder=7) #R123 -0.05; R12 -0.03; R3 -0.32
        #plot1.plot(positions, np.array(dict_of_wigs['HETU 323 CTD+Rif-'])-0.28, linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HETU CTD+Rif- ({TU_sets_v["HETU 323 CTD-Rif-"]})', zorder=8) #R123 -0.20; R123 -0.25; R123 -0.28
        #plot1.plot(positions, np.array(dict_of_wigs['HETU 323 CTD+Rif+'])-0.24, linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'HETU CTD+Rif+ ({TU_sets_v["HETU 323 CTD-Rif-"]})', zorder=7) #R23 -0.15; R23 -0.20; R23 -0.24
        #HETU, no ompX
        plot1.plot(positions, np.array(dict_of_wigs['HETU 321 no ompX CTD-Rif-'])-0.05, linestyle='-', color='#b08642', linewidth=3, alpha=0.8, label=f'HETU CTD-/Rif- ({TU_sets_v["HETU 321 no ompX CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.05; R12 -0.05; R3 -0.25
        #plot1.plot(positions, np.array(dict_of_wigs['HETU 321 no ompX CTD-Rif+'])-0.12, linestyle='--', color='#b08642', linewidth=5, alpha=1, label=f'HETU CTD-/Rif+ ({TU_sets_v["HETU 321 no ompX CTD-Rif-"]})', zorder=7) #Def linewidth=0.8 #R123 -0.05; R12 -0.03; R3 -0.32
        #plot1.plot(positions, np.array(dict_of_wigs['HETU 321 no ompX CTD+Rif-'])-0.25, linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HETU CTD+Rif- ({TU_sets_v["HETU 321 no ompX CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.20; R123 -0.25; R123 -0.28
        #plot1.plot(positions, np.array(dict_of_wigs['HETU 321 no ompX CTD+Rif+'])-0.25, linestyle='-.', color='#b08642', linewidth=1, alpha=1, label=f'HETU CTD+Rif+ ({TU_sets_v["HETU 321 no ompX CTD-Rif-"]})', zorder=7) #Def linewidth=1 #R23 -0.15; R23 -0.20; R23 -0.24 
        #rRNA
        #plot1.plot(positions, dict_of_wigs['rRNA 7 CTD-Rif-'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=0.8, label=f'rRNA CTD-Rif- ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0; R12 -0; R3 -0
        #plot1.plot(positions, dict_of_wigs['rRNA 7 CTD-Rif+'], linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'rRNA CTD-Rif+ ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=7) #Def linewidth=1 #R123 -0; R12 -0; R3 -0
        #plot1.plot(positions, dict_of_wigs['rRNA 7 CTD+Rif-'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'rRNA CTD+Rif- ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0; R123 -0; R123 -0
        #plot1.plot(positions, dict_of_wigs['rRNA 7 CTD+Rif+'], linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'rRNA CTD+Rif+ ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=7) #Def linewidth=1 #R23 -0; R23 -0; R23 -0 
        #tRNA
        #plot1.plot(positions, np.array(dict_of_wigs['tRNA 49 CTD-Rif-'])-0.1, linestyle='-', color='#757d8b', linewidth=1.5, alpha=0.8, label=f'tRNA CTD-Rif- ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.1; R12 -0.1; R3 -0.1
        #plot1.plot(positions, np.array(dict_of_wigs['tRNA 49 CTD-Rif+'])-0.05, linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'tRNA CTD-Rif+ ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=7) #Def linewidth=1 #R123 -0.1; R12 -0.05; R3 -0.05
        #plot1.plot(positions, np.array(dict_of_wigs['tRNA 49 CTD+Rif-'])-0.25, linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'tRNA CTD+Rif- ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5 #R123 -0.2; R123 -0.25; R123 -0.25
        #plot1.plot(positions, np.array(dict_of_wigs['tRNA 49 CTD+Rif+'])-0.15, linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'tRNA CTD+Rif+ ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=7) #Def linewidth=1 #R23 -0.1; R23 -0.15; R23 -0.15   
                     
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks+list(range(-300, 200, 20))) #Start
    #plot1.set_xticks(ticks+list(range(4800, 5300, 20))) #End
    plot1.set_xticks(ticks, minor='False')
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1) 
    plot1.set_xticks([0, length], minor='True')
    #plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    #plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    #plot1.axvline(length, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.set_yticks([1], minor='True')
    #plot1.set_xlim(-300, 200) #Start
    #plot1.set_xlim(4800, 5300) #End
    #plot1.grid(axis='y', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)
    #plot1.axhline(1, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'{set_name} fold enrichment', size=20)
    plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_Rif_plus_CTD_minus_Rif_minus_CTD_minus_rtRNA_Start_{win_width}bp_nd_with_body_{length}_bp.png', dpi=400, figsize=(10, 6))  #Def size - 10, 6; Closer look - 3, 6
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
    if set_type=='genes':
        #LEG_270
        plot1.plot(positions_sm, dict_of_wigs_sm['All genes LF'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'All genes LF ({TU_sets_v["All genes LF"]})', zorder=6)
        plot1.plot(positions_sm, dict_of_wigs_sm['All genes LR'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'All genes LR ({TU_sets_v["All genes LR"]})', zorder=5)   
        #LEG_370
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEG_370'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'LEG ({TU_sets_v["LEG_370"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['tRNA Rif'], linestyle='--', color='#E692A9', linewidth=1, alpha=0.8, label='tRNA Rif (86)', zorder=5)        
        #All_genes
        plot1.plot(positions_sm, dict_of_wigs_sm['All genes F'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All genes F ({TU_sets_v["All genes F"]})', zorder=10)
        plot1.plot(positions_sm, dict_of_wigs_sm['All genes R'], linestyle='--', color='#333738', linewidth=1.5, alpha=0.8, label=f'All genes R ({TU_sets_v["All genes LR"]})', zorder=9)
        #HEG_270
        plot1.plot(positions_sm, dict_of_wigs_sm['All genes RF'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'All genes RF ({TU_sets_v["All genes RF"]})', zorder=8)
        plot1.plot(positions_sm, dict_of_wigs_sm['All genes RR'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'All genes RR ({TU_sets_v["All genes RR"]})', zorder=7)
        #HEG_370
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEG_370'], linestyle='-', color='#e4d1b4', linewidth=1, alpha=1.5, label=f'HEG ({TU_sets_v["HEG_370"]})', zorder=4)
        #plot1.plot(positions_sm, dict_of_wigs_sm['ncRNA Rif'], linestyle='--', color='#FFC000', linewidth=0.8, alpha=0.8, label='ncRNA Rif (18)', zorder=3)  
    ##Operons below.
    elif set_type=="operons":
        #LEO_144
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEO 144'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LEO ({TU_sets_v["LEO 144"]})', zorder=6)
        #plot1.plot(positions_sm, dict_of_wigs_sm['LEO_144_Rif'], linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LEO Rif ({TU_sets_v["LEO_144_Rif"]})', zorder=5)   
        #LEO_186 or SAO 27
        #plot1.plot(positions_sm, dict_of_wigs_sm['SAO 27'], linestyle='-', color='#bec1cb', linewidth=1, alpha=1.5, label=f'SAO ({TU_sets_v["SAO 27"]})', zorder=6)
        plot1.plot(positions_sm, dict_of_wigs_sm['All operons LF'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=1.5, label=f'All operons LF ({TU_sets_v["All operons LF"]})', zorder=6)
        plot1.plot(positions_sm, dict_of_wigs_sm['All operons LR'], linestyle='--', color='#757d8b', linewidth=1, alpha=0.8, label=f'All operons LR ({TU_sets_v["All operons LR"]})', zorder=5)        
        #All_operons
        plot1.plot(positions_sm, dict_of_wigs_sm['All operons F'], linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All operons F ({TU_sets_v["All operons F"]})', zorder=10)
        plot1.plot(positions_sm, dict_of_wigs_sm['All operons R'], linestyle='--', color='#333738', linewidth=1.5, alpha=0.8, label=f'All operons R ({TU_sets_v["All operons R"]})', zorder=9)
        #HEO_144
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEO 144'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HEO no rRNA, tRNA ({TU_sets_v["HEO 144"]})', zorder=8)
        #plot1.plot(positions_sm, dict_of_wigs_sm['HEO_no_tRNA_rRNA_144_Rif'], linestyle='--', color='#b08642', linewidth=1, alpha=1, label=f'HEO no rRNA, tRNA Rif ({TU_sets_v["HEO_no_tRNA_rRNA_144_Rif"]})', zorder=7)
        #HEO_186 or rRNA operons or LAO_27
        #plot1.plot(positions_sm, dict_of_wigs_sm['LAO 27'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'LAO ({TU_sets_v["LAO 27"]})', zorder=4)
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA_operons'], linestyle='-', color='#e4d1b4', linewidth=1.5, alpha=1, label=f'rRNA operons ({TU_sets_v["rRNA_operons"]})', zorder=4)
        plot1.plot(positions_sm, dict_of_wigs_sm['All operons RF'], linestyle='-', color='#b08642', linewidth=1.5, alpha=1.5, label=f'All operons RF ({TU_sets_v["All operons RF"]})', zorder=4)
        plot1.plot(positions_sm, dict_of_wigs_sm['All operons RR'], linestyle='--', color='#b08642', linewidth=1, alpha=0.8, label=f'All operons RR ({TU_sets_v["All operons RR"]})', zorder=3)  
    ##Transcription units below.
    elif set_type=="transcripts":
        #LETU
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['LETU 245 CTD-Rif-'])-0.07, linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LETU CTD-Rif- ({TU_sets_v["LETU 245 CTD-Rif-"]})', zorder=6)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['LETU 245 CTD-Rif+'])-0.45, linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'LETU CTD-Rif+ ({TU_sets_v["LETU 245 CTD-Rif-"]})', zorder=5)  
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['LETU 245 CTD+Rif-'])-0.15, linestyle='-', color='#b08642', linewidth=1.5, alpha=1, label=f'LETU CTD+Rif- ({TU_sets_v["LETU 245 CTD-Rif-"]})', zorder=6)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['LETU 245 CTD+Rif+'])-0.25, linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'LETU CTD+Rif+ ({TU_sets_v["LETU 245 CTD-Rif-"]})', zorder=5)             
        #LETU, no ybiI
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU 244 no ybiI CTD-Rif-'])-0.05, linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LETU CTD-/Rif- ({TU_sets_v["LETU 244 no ybiI CTD-Rif-"]})', zorder=6)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU 244 no ybiI CTD-Rif+'])-0.1, linestyle='--', color='#757d8b', linewidth=1, alpha=1, label=f'LETU CTD-/Rif+ ({TU_sets_v["LETU 244 no ybiI CTD-Rif-"]})', zorder=5)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU 244 no ybiI CTD+Rif-'])-0.25, linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LETU CTD+Rif- ({TU_sets_v["LETU 244 no ybiI CTD-Rif-"]})', zorder=6)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU 244 no ybiI CTD+Rif+'])-0.25, linestyle='-.', color='#757d8b', linewidth=1, alpha=1, label=f'LETU CTD+Rif+ ({TU_sets_v["LETU 244 no ybiI CTD-Rif-"]})', zorder=5)  
        #LETU, no ybiI, no appY
        plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU 248 no ybiI no appY CTD-Rif-'])-0.03, linestyle='-', color='#B6B8BD', linewidth=1.5, alpha=1, label=f'LETU CTD-/Rif- ({TU_sets_v["LETU 248 no ybiI no appY CTD-Rif-"]})', zorder=6)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU 248 no ybiI no appY CTD-Rif+'])-0.1, linestyle='--', color='#B6B8BD', linewidth=0.8, alpha=1, label=f'LETU CTD-/Rif+ ({TU_sets_v["LETU 248 no ybiI no appY CTD-Rif-"]})', zorder=5)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU 248 no ybiI no appY CTD+Rif-'])-0.25, linestyle='-', color='#757d8b', linewidth=1.5, alpha=1, label=f'LETU CTD+Rif- ({TU_sets_v["LETU 248 no ybiI no appy CTD-Rif-"]})', zorder=6)
        #plot1.plot(positions_sm,  np.array(dict_of_wigs_sm['LETU 248 no ybiI no appY CTD+Rif+'])-0.25, linestyle='-.', color='#757d8b', linewidth=1, alpha=1, label=f'LETU CTD+Rif+ ({TU_sets_v["LETU 248 no ybiI no appY CTD-Rif-"]})', zorder=5)          
        #All_TUs
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs CTD-Rif-'])-0.03, linestyle='-', color='#333738', linewidth=2, alpha=2, label=f'All TUs CTD-Rif- ({TU_sets_v["All TUs CTD-Rif-"]})', zorder=10)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs CTD-Rif+'])-0.1, linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All TUs CTD-Rif+ ({TU_sets_v["All TUs CTD-Rif-"]})', zorder=9)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs CTD+Rif-'])-0.27, linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'All TUs CTD+Rif- ({TU_sets_v["All TUs CTD-Rif-"]})', zorder=10)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs CTD+Rif+'])-0.24, linestyle='--', color='#e4d1b4', linewidth=1, alpha=0.8, label=f'All TUs CTD+Rif+ ({TU_sets_v["All TUs CTD-Rif-"]})', zorder=9)        
        #All_TUs no tRNA, rRNA.
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no tRNA, rRNA CTD-Rif-'])-0.05, linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All TUs no tRNA, rRNA CTD-/Rif- ({TU_sets_v["All TUs no tRNA, rRNA CTD-Rif-"]})', zorder=10)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no tRNA, rRNA CTD-Rif+'])-0.1, linestyle='--', color='#333738', linewidth=1, alpha=0.8, label=f'All TUs no tRNA, rRNA CTD-/Rif+ ({TU_sets_v["All TUs no tRNA, rRNA CTD-Rif-"]})', zorder=9)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no tRNA, rRNA CTD+Rif-'])-0.25, linestyle='-', color='#333738', linewidth=2, alpha=0.8, label=f'All TUs no tRNA, rRNA CTD+Rif- ({TU_sets_v["All TUs no tRNA, rRNA CTD-Rif-"]})', zorder=10)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['All TUs no tRNA, rRNA CTD+Rif+'])-0.25, linestyle='-.', color='#333738', linewidth=1, alpha=0.8, label=f'All TUs no tRNA, rRNA CTD+Rif+ ({TU_sets_v["All TUs no tRNA, rRNA CTD-Rif-"]})', zorder=9)        
        #HETU
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU 323 CTD-Rif-'])-0.25, linestyle='-', color='#757d8b', linewidth=1.5, alpha=0.8, label=f'HETU CTD-Rif- ({TU_sets_v["HETU 323 CTD-Rif-"]})', zorder=8)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU 323 CTD-Rif+'])-0.32, linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'HETU CTD-Rif+ ({TU_sets_v["HETU 323 CTD-Rif-"]})', zorder=7)  
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU 323 CTD+Rif-'])-0.28, linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HETU CTD+Rif- ({TU_sets_v["HETU 323 CTD-Rif-"]})', zorder=8)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU 323 CTD+Rif+'])-0.24, linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'HETU CTD+Rif+ ({TU_sets_v["HETU 323 CTD-Rif-"]})', zorder=7)          
        #HETU, no ompX
        plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU 321 no ompX CTD-Rif-'])-0.03, linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HETU CTD-/Rif- ({TU_sets_v["HETU 321 no ompX CTD-Rif-"]})', zorder=8)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU 321 no ompX CTD-Rif+'])-0.1, linestyle='--', color='#b08642', linewidth=0.8, alpha=1, label=f'HETU CTD-/Rif+ ({TU_sets_v["HETU 321 no ompX CTD-Rif-"]})', zorder=7)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU 321 no ompX CTD+Rif-'])-0.25, linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'HETU CTD+Rif- ({TU_sets_v["HETU 321 no ompX CTD-Rif-"]})', zorder=8)
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['HETU 321 no ompX CTD+Rif+'])-0.25, linestyle='-.', color='#b08642', linewidth=1, alpha=1, label=f'HETU CTD+Rif+ ({TU_sets_v["HETU 321 no ompX CTD-Rif-"]})', zorder=7)        
        #rRNA
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA 7 CTD-Rif-'], linestyle='-', color='#757d8b', linewidth=1.5, alpha=0.8, label=f'rRNA CTD-Rif- ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA 7 CTD-Rif+'], linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'rRNA CTD-Rif+ ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=7) #Def linewidth=1
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA 7 CTD+Rif-'], linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'rRNA CTD+Rif- ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions_sm, dict_of_wigs_sm['rRNA 7 CTD+Rif+'], linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'rRNA CTD+Rif+ ({TU_sets_v["rRNA 7 CTD-Rif-"]})', zorder=7) #Def linewidth=1        
        #tRNA
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['tRNA 49 CTD-Rif-'])-0.1, linestyle='-', color='#757d8b', linewidth=1.5, alpha=0.8, label=f'tRNA CTD-Rif- ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['tRNA 49 CTD-Rif+'])-0.05, linestyle='--', color='#333738', linewidth=1, alpha=1, label=f'tRNA CTD-Rif+ ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=7) #Def linewidth=1  
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['tRNA 49 CTD+Rif-'])-0.25, linestyle='-', color='#b08642', linewidth=1.5, alpha=0.8, label=f'tRNA CTD+Rif- ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=8) #Def linewidth=1.5
        #plot1.plot(positions_sm, np.array(dict_of_wigs_sm['tRNA 49 CTD+Rif+'])-0.15, linestyle='--', color='#e4d1b4', linewidth=1, alpha=1, label=f'tRNA CTD+Rif+ ({TU_sets_v["tRNA 49 CTD-Rif-"]})', zorder=7) #Def linewidth=1           
        
    #plot1.set_ylim(0.75, 1.6) #(0.75, 1.6) for FE; (-0.6, 0.5) for ded FE
    ticks=np.arange(-win_width,win_width+length+1,length).tolist()
    plot1.set_xticks(ticks)
    ticks_lables=ticks
    ticks_lables[ticks.index(0)]='TS'
    ticks_lables[ticks.index(length)]='TE'
    ticks_lables1=ticks_lables[:ticks_lables.index('TE')+1]+np.arange(length,win_width+1,length).tolist()
    plot1.set_xticklabels(ticks_lables1)
    plot1.set_xticks([0, length], minor='True')
    #plot1.grid(axis='x', which='minor', color='black', linestyle='--', alpha=0.7, linewidth=1)  
    #plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    #plot1.axvline(length, color='black', linestyle='--', alpha=0.7, linewidth=1)    
    plot1.set_yticks([1], minor='True')
    #plot1.grid(axis='y', which='minor', color='grey', linestyle='--', alpha=0.7, linewidth=1)   
    #plot1.axhline(1, color='black', linestyle='--', alpha=0.7, linewidth=1)
    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'{set_name} fold enrichment', size=20)
    plot1.set_title(f'{set_name} signal over {set_type}', size=20)     
    plt.savefig(f'{output_path}\\{set_name}_FE_over_{set_type}_Rif_plus_CTD_minus_Rif_minus_CTD_minus_smoothed_rtRNA_{win_width}bp_nd_with_body_{length}bp_smoothed_{2*sm_window}bp.png', dpi=400, figsize=(10, 6), transparent=True)   
    #plt.show()
    plt.close()    
    return


#DRIP-Seq.
plot_FE_all_expression_gg_Rif_no_Rif(Wig_data_in_dict_operons, Sm_window, Out_path, Signal_name, 'operons')
plot_FE_all_expression_gg_Rif_no_Rif(Wig_data_in_dict_genes, Sm_window, Out_path, Signal_name, 'genes')


#EcTopoI
#plot_FE_all_expression_gg_Rif_no_Rif(Wig_data_in_dict_transcripts, Sm_window, Out_path, Signal_name, Set_type)
#RNA-pol
#plot_FE_all_expression_gg_Rif_no_Rif(Wig_data_in_dict_transcripts_RNApol, Sm_window, Out_path, Signal_name, Set_type)
