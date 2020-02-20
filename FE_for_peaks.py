###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes narrowPeaks files with list of regions (e.g. ChIP-Seq peaks identified in different biological replicas),
#returns signal (fold enrichment) of them.
###############################################


#######
#Packages to be imported.
#######

import numpy as np
import Bio
from Bio import SeqIO


#Path to the working directory with NarrowPeak files.
PWD_peaks="C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\Reproducible_peaks\\"
#Dictionary of pathes to NarrowPeak file with peaks coordinates (MACS2 output).
Peaks_data={'TopoA_CTD_noRif_rep123_thr_3_nm_0.001_peaks' : PWD_peaks + "TopoA_CTD_noRif_rep123_thr_3_nm_0.001_peaks.narrowPeak",
            'TopoA_CTD_Rif_rep12_thr_2_nm_0.001_peaks' : PWD_peaks + "TopoA_CTD_Rif_rep12_thr_2_nm_0.001_peaks.narrowPeak",
            'TopoA_noCTD_noRif_rep123_thr_3_nm_0.001_peaks' : PWD_peaks + "TopoA_noCTD_noRif_rep123_thr_3_nm_0.001_peaks.narrowPeak",
            'TopoA_noCTD_Rif_rep123_thr_3_nm_0.001_peaks' : PWD_peaks + "TopoA_noCTD_Rif_rep123_thr_3_nm_0.001_peaks.narrowPeak",
            }

Peaks_data_shared=PWD_peaks + "TopoA_all_exp_shared_nm_0.001_peaks.narrowPeak"

#Path to the working directory with WIG files.
PWD_FE="C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Fold_enrichment\\"
#Dictionary of pathes to WIG file with fold enrichment tracks.
FE_data={'TopoA_CTD_noRif_rep123_thr_3_nm_0.001_peaks' : PWD_FE + "TopA_ChIP_CTD_plus_Rif_minus_average_FE.wig",
         'TopoA_CTD_Rif_rep12_thr_2_nm_0.001_peaks' : PWD_FE + "TopA_ChIP_CTD_plus_Rif_plus_average_FE.wig",
         'TopoA_noCTD_noRif_rep123_thr_3_nm_0.001_peaks' : PWD_FE + "TopA_ChIP_CTD_minus_Rif_minus_average_FE_1_2_3.wig",
         'TopoA_noCTD_Rif_rep123_thr_3_nm_0.001_peaks' : PWD_FE + "TopA_ChIP_CTD_minus_Rif_plus_average_FE_1_2_3.wig",
         }

#Path to the reference genome (e.g. E_coli_w3110_G_Mu.fasta).
Genome="C:\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\E_coli_w3110_G_Mu.fasta"


#######
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
        genome_id=record.name
    return len(genomefa), genomefa, genome_id

#######
#Opens and reads BED or narrowPeak files.
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
#Parses WIG file.
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
#Return FE for peaks.
#######

def return_peaks_FE(peaks_array, wig_FoldEnrich):
    #Iterate peaks, return FE of each region.
    peaks_array_FE=[]
    for i in range(len(peaks_array)):
        peak_start=peaks_array[i][0]
        peak_end=peaks_array[i][1]
        peak_signal=np.mean(wig_FoldEnrich[peak_start:peak_end])
        peaks_array_FE.append([peak_start, peak_end, peak_signal])
        
    #Check FE obtained.
    FE_list=[]
    for i in range(len(peaks_array_FE)):
        FE_list.append(peaks_array_FE[i][2])
    print(np.mean(FE_list), np.std(FE_list))
    return peaks_array_FE

#######
#Write peaks in broadPeak format.
#######   

def write_bed(peaks_ar, chrom_name, outpath):
    fileout=open(outpath, 'w')
    for i in range(len(peaks_ar)):
        fileout.write(chrom_name+'\t'+str(peaks_ar[i][0])+'\t'+str(peaks_ar[i][1])+'\tPeak_'+str(i)+'\t10\t.\t' + str(peaks_ar[i][2]) + '\t-1\t-1\n')
    fileout.close()
    return

#######
#Wrapper function.
#######

def wrapper(genome_path, Peaks_data, Peaks_data_shared, FE_data, outpath):
    #Read reference genome data.
    genome_length, genome_seq, chrom_name=read_genome(genome_path)
    
    ##Peaks lists correspond to WIG file with FE.
    #Read NarrowPeak.
    dict_of_peaks={}
    for peaks_name, peaks_path in Peaks_data.items():    
        dict_of_peaks[peaks_name]=deletions_info(peaks_path)
    #Read WIG.
    #Contains data of all replicas in separate arrays.
    dict_of_wigs={}
    for wig_name, wig_path in FE_data.items():
        dict_of_wigs[wig_name]=wig_parsing(wig_path)     
    #Return FE of peaks.
    for peaks_name, peaks_data in dict_of_peaks.items():
        print('Now working with ' + peaks_name)
        peaks_ar=peaks_data
        wig_FE=dict_of_wigs[peaks_name]
        peaks_ar_FE=return_peaks_FE(peaks_ar, wig_FE)
        dict_of_peaks[peaks_name]=peaks_ar_FE
    #Write data.
    for peaks_name, peaks_data in dict_of_peaks.items():
        write_bed(peaks_data, chrom_name, outpath+peaks_name+'_FE.narrowPeak')
    
    ##Peaks shared between all experimental conditions and replicas.
    Peaks_shared=deletions_info(Peaks_data_shared)
    #Return FE of shared peaks.
    Peaks_shared_dict={}
    for wig_name, wig_data in dict_of_wigs.items():
        print('Now working with shared peaks of ' + wig_name)
        Peaks_shared_dict[wig_name]=return_peaks_FE(Peaks_shared, wig_data) 
    #Write shared peaks data.
    for peaks_name, peaks_data in Peaks_shared_dict.items():
        write_bed(peaks_data, chrom_name, outpath+'TopoA_all_exp_shared_peaks_'+peaks_name+'_FE.narrowPeak')
        
    return


wrapper(Genome, Peaks_data, Peaks_data_shared, FE_data, PWD_peaks)