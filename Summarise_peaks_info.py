###############################################
##Dmitry Sutormin, 2020##
##TopoI ChIP-Seq analysis##

#The returns signal (enrichment, score, GC%, other continuously distributed character) for
#sets of genomic intervals (Peaks, TUs, BIMEs-1, BIMEs-2, IHF sites, Fis sites, H-NS sites, MatP sites, etc.).
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import binom

#######
#Variables to be defined.
#######


#Path to the working directory with files containing info about intervals.
PWD_peaks="C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peaks_info\\"
#Input: Intervals (e.g. EcTopoI peaks) (tab-delimited with headers).
path_to_intervals_sets={'CTD-/Rif-'    : [PWD_peaks + "EcTopoI_noCTD_noRif_rep123_nm_0.001_peaks_width_GC_FE.csv", 'EcTopoI_noCTD_noRif_nm_0.001_peaks'],
                        'CTD-/Rif+'    : [PWD_peaks + "EcTopoI_noCTD_Rif_rep123_nm_0.001_peaks_width_GC_FE.csv", 'EcTopoI_noCTD_Rif_nm_0.001_peaks'],
                        'CTD+/Rif-'    : [PWD_peaks + "EcTopoI_CTD_noRif_rep123_nm_0.001_peaks_width_GC_FE.csv", 'EcTopoI_CTD_noRif_nm_0.001_peaks'],
                        'CTD+/Rif+'    : [PWD_peaks + "EcTopoI_CTD_Rif_rep12_nm_0.001_peaks_width_GC_FE.csv", 'EcTopoI_CTD_Rif_nm_0.001_peaks'],
                        'Shared peaks' : [PWD_peaks + "Shared_peaks_TopoA_nm_0.001_width_GC_FE_for_diffr_experiments.csv", 'Shared_peaks_TopoA_nm_0.001_peaks'],
                        }

#Path to the working directory with files containing continuous signal.
PWD_wig="C:\Sutor\Science\E_coli_ChIP-Seqs\All_tracks\\"
#Input: Continuous data (e.g. EcTopoI fold enrichment) (WIG).
path_to_cont_data={'HNS Kahramanoglou'   : PWD_wig + "Kahramanoglou_HNS_IP_ME.wig",
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

#Output: path to the dir to store output
Outputpath="C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peaks_info\\"
if not os.path.exists(Outputpath):
    os.makedirs(Outputpath)
    

#######
#Read intervals data.
#######

def peaks_data_pars(intervals_sets_path_dict):
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]] #Deletions in E. coli w3110 strain that corresponds to DY330 strain.
    #Read data.
    datasets_dict={}
    for set_name, filepath in intervals_sets_path_dict.items():
        filepath_name=filepath[1]
        dataset_df=pd.read_csv(filepath[0], sep='\t', header=0, index_col=False, dtype={'Start' : np.int64, 'End' : np.int64})
        dataset_df.Start=dataset_df.Start.round()
        dataset_df.End=dataset_df.End.round()
        #Remove items falling into delelted or masked regions.
        maska=(dataset_df['Start']>=0)
        for j in range(len(deletions)):
            maska=maska & (~((dataset_df['Start']<deletions[j][1]) & (dataset_df['Start']>deletions[j][0])) | ~((dataset_df['End']<deletions[j][1]) & (dataset_df['End']>deletions[j][0])))
        dataset_df_masked=dataset_df[maska]  
        #Keep data.            
        datasets_dict[set_name]=dataset_df_masked
        print("Number of " + str(set_name) + " regions: " + str(len(dataset_df_masked.index)))
    return datasets_dict


#######
#Parsing WIG file.
#######

def score_data_parser(inpath_dict):
    cont_char_dict={}
    for param_name, inpath in inpath_dict.items():
        param_file=open(inpath, 'r')
        ar=[]
        for line in param_file:
            line=line.rstrip().split(' ')
            if line[0] not in ['track', 'fixedStep']:
                ar.append(float(line[0]))
        param_file.close()
        print('Whole genome average ' + str(param_name) + ' : ' + str(sum(ar)/len(ar)))
        cont_char_dict[param_name]=ar 
    return cont_char_dict


#######
#Mask deletions.
#######

def mask_array(ar, regions_to_mask):
    #Mask deletions or smth in FE array.
    maska=[0]*len(ar)
    for deletion in regions_to_mask:
        del_len=deletion[1]-deletion[0]
        for i in range(del_len):
            maska[deletion[0]+i]=1
    ar_masked=np.ma.masked_array(ar, mask=maska)  
    ar_non_del_only=list(ar_masked[~ar_masked.mask])
    print(len(ar), len(ar_non_del_only))
    return ar_non_del_only


#######
#Return continous characteristics of intervals, mask deletions and intervals regions.
#######

def cont_char_of_intervals(intervals_sets_dict, cont_char_dict, name_of_intervals, name_of_cont_char):
    intervals_data=intervals_sets_dict[name_of_intervals]
    cont_char=cont_char_dict[name_of_cont_char]
    #Return enrichment of intervals.
    intervals_param_ar_tg=[]
    intervals_param_ar_sp=[]
    intervals=[]
    for index, row_data in intervals_data.iterrows():
        interval_start=int(row_data['Start'])
        interval_end=int(row_data['End'])
        intervals.append([interval_start, interval_end])
        intervals_param_ar_tg+=cont_char[interval_start:interval_end]   
        intervals_param_ar_sp.append(np.mean(cont_char[interval_start:interval_end]))
    
    #Prepare continuos character without deleted regions and regions defined by intervals. 
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]] #Deletions in E. coli w3110 strain that corresponds to DY330 strain.
    regions_to_mask=deletions+intervals
    cont_char_masked=mask_array(cont_char, regions_to_mask)
    return intervals_param_ar_tg, cont_char_masked, intervals_param_ar_sp


#######
#Functions wrapper.
#######

def func_wrapper(intervals_sets_path_dict, cont_char_path_dict, outpath):
    #Read intervals.
    intervals_sets_dict=peaks_data_pars(intervals_sets_path_dict)
    #Read continuous data.
    cont_char_dict=score_data_parser(cont_char_path_dict)
    #Process intervals, return signals of intervals.
    intervals_sets_dict_signals=intervals_sets_dict
    for intervals_set_name, intervals_set_data in intervals_sets_dict.items():
        for cont_char_name, cont_char_data in cont_char_dict.items():
            print('Now working with ' + intervals_set_name + ' and ' + cont_char_name)
            intervals_signal, cont_char_masked, intervals_signal_means=cont_char_of_intervals(intervals_sets_dict, cont_char_dict, intervals_set_name, cont_char_name)
            intervals_sets_dict_signals[intervals_set_name][cont_char_name]=intervals_signal_means
            
    #Write data.
    for intervals_set_name, intervals_set_data in intervals_sets_dict_signals.items():
        intervals_set_data.to_csv(outpath+intervals_sets_path_dict[intervals_set_name][1]+'_more_info.csv', sep='\t', index=False)
        
    return

func_wrapper(path_to_intervals_sets, path_to_cont_data, Outputpath)

print('Script ended its work succesfully!') 