###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes NarrowPeaks with ChIP-Seq peaks,
#Returns sequences under them,
#Computes GC%, writes multifasta file.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import GC as GC_count
import pandas as pd


#Path to the working directory with NarrowPeak files.
PWD_peaks="C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\Peak_calling\Reproducible_peaks\\"
#Path to NarrowPeak file with peaks coordinates (MACS2 output).
Peaks_dict={'CTD-/Rif-' : [PWD_peaks + "TopoA_noCTD_noRif_rep123_thr_3_nm_0.001_peaks_FE.narrowPeak", 'EcTopoI_noCTD_noRif_rep123_nm_0.001_peaks'],
            'CTD-/Rif+' : [PWD_peaks + "TopoA_noCTD_Rif_rep123_thr_3_nm_0.001_peaks_FE.narrowPeak", 'EcTopoI_noCTD_Rif_rep123_nm_0.001_peaks'],
            'CTD+/Rif-' : [PWD_peaks + "TopoA_CTD_noRif_rep123_thr_3_nm_0.001_peaks_FE.narrowPeak", 'EcTopoI_CTD_noRif_rep123_nm_0.001_peaks'],
            'CTD+/Rif+' : [PWD_peaks + "TopoA_CTD_Rif_rep12_thr_2_nm_0.001_peaks_FE.narrowPeak", 'EcTopoI_CTD_Rif_rep12_nm_0.001_peaks'],
            'Shared peaks' : [PWD_peaks + "TopoA_all_exp_shared_peaks_TopoA_noCTD_noRif_rep123_thr_3_nm_0.001_peaks_FE.narrowPeak", 'TopoA_all_exp_shared_peaks_nm_0.001_peaks'],
            }

Shared_peaks_dict={'CTD-/Rif- S' : [PWD_peaks + "TopoA_all_exp_shared_peaks_TopoA_noCTD_noRif_rep123_thr_3_nm_0.001_peaks_FE.narrowPeak", 'Shared_peaks_TopoA_noCTD_noRif_nm_0.001'],
                   'CTD-/Rif+ S' : [PWD_peaks + "TopoA_all_exp_shared_peaks_TopoA_noCTD_Rif_rep123_thr_3_nm_0.001_peaks_FE.narrowPeak", 'Shared_peaks_TopoA_noCTD_Rif_nm_0.001'],
                   'CTD+/Rif- S' : [PWD_peaks + "TopoA_all_exp_shared_peaks_TopoA_CTD_noRif_rep123_thr_3_nm_0.001_peaks_FE.narrowPeak", 'Shared_peaks_TopoA_CTD_noRif_nm_0.001'],
                   'CTD+/Rif+ S' : [PWD_peaks + "TopoA_all_exp_shared_peaks_TopoA_CTD_Rif_rep12_thr_2_nm_0.001_peaks_FE.narrowPeak", 'Shared_peaks_TopoA_CTD_Rif_nm_0.001'],
                   }


#Path to the reference genome (e.g. E_coli_w3110_G_Mu.fasta).
Genome="C:\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\E_coli_w3110_G_Mu.fasta"
#Path to the working directory.
pwd="C:\Sutor\Science\TopoI-ChIP-Seq\Data_analysis\\"
#Outpath.
Path_out=pwd + "Seq_under_peaks\\"
if not os.path.exists(Path_out):
    os.makedirs(Path_out)
    

#######
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
    return genomefa

#######
#Opens and reads BED or NarrowPeak files.
#######

def deletions_info(del_path):
    del_ar=[]
    del_signal=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        del_ar.append([int(line[1]), int(line[2])])
        del_signal.append(float(line[6]))
    filein.close()
    return del_ar, del_signal

#######
#Return sub-sequences from reference genome by coordinates provided.
#######

def return_seqs(peaks_ar, genome):
    seq_dict={}
    for peak in peaks_ar:
        seq_dict[peak[0]]=genome[peak[0]:peak[1]]
    return seq_dict

#######
#Writes sequences as mfa file.
#######

def write_seqs(seqs_dict, pathout):
    fileout=open(pathout, 'w')
    for k in seqs_dict:
        fileout.write(f'>{k}\n{seqs_dict[k]}\n')
    fileout.close()
    return

#######
#Genome binning with average peak width, compute distribution of GC% of bins.
#######

def genome_bin_GC_dist(seq_dict, genome):
    #Calculate mean peaks length (width of a bin).
    #Keep GC% of peaks.
    peaks_gc_values=[]
    peaks_width=[]
    for k in seq_dict:
        peaks_gc_values.append(GC_count(seq_dict[k]))
        peaks_width.append(len(seq_dict[k]))
    median_len=int(np.median(peaks_width))-1
    #Bin genome and keep GC% of bins.
    number_of_genome_bin=int(len(genome)/median_len)-1
    genome_gc_values=[]
    for i in range(number_of_genome_bin):
        genome_gc_values.append(GC_count(genome[median_len*i:median_len*(i+1)]))
    
    GC_genome=GC_count(genome)
    print(f'GC% of the reference genome: {GC_genome}%') 
    print(f'GC% of the reference genome binned: {np.mean(genome_gc_values)}%')  
    print(f'GC% of the peaks regions: {np.mean(peaks_gc_values)}%')
    return GC_genome, median_len, genome_gc_values, peaks_gc_values, peaks_width

#######
#Plots distribution of peaks GC% in comparision to genome GC%.
#Plots distribution of peaks widths.
#Plots distribution of peaks fold enrichment (FE).
#######

def plot_distrib(set_name, set_data_ar, Min_GC, Max_GC, Min_width, Max_width, Min_FE, Max_FE, pathout):    
    #Plot GC% distributions, peaks width distributions.
    fig=plt.figure(figsize=(10, 3), dpi=100)
    ##GC%
    plot0=plt.subplot2grid((1,3),(0,0), rowspan=1, colspan=1)  
    genome_gc_values=set_data_ar[2]
    peaks_gc_values=set_data_ar[3]
    median_len=set_data_ar[1]    
    bins_GC=np.linspace(Min_GC, Max_GC, int(Max_GC-Min_GC))
    plot0.hist(genome_gc_values, bins_GC, color='#ff993a', edgecolor='black', alpha=0.5, density=True, label='Genome GC%') #Genome GC%
    plot0.hist(peaks_gc_values, bins_GC, color='#a826a6', edgecolor='black', alpha=0.3, density=True, label='Peaks GC%') #Peaks GC%
    plot0.annotate(f'Mean genome GC%={round(np.mean(genome_gc_values),2)}%', xy=(0.03, 0.9), xycoords='axes fraction', size=8)
    plot0.annotate(f'Mean peaks GC%={round(np.mean(peaks_gc_values),2)}%', xy=(0.03, 0.83), xycoords='axes fraction', size=8)
    plot0.annotate(f'Peaks median width={median_len}bp', xy=(0.03, 0.76), xycoords='axes fraction', size=8)
    plot0.set_xlabel('GC%', size=16)
    plot0.set_ylabel('Fraction of peaks(bins)', size=16)
    plot0.set_title('Peaks GC%\n'+set_name, size=14)
    plot0.legend(loc='center left', frameon=False, fontsize=8, handlelength=1)
    #Peaks width
    plot1=plt.subplot2grid((1,3),(0,1), rowspan=1, colspan=1)  
    peaks_width=set_data_ar[4]
    bins_width=np.linspace(Min_width, Max_width, 50)
    plot1.hist(peaks_width, bins_width, color='#3a41ff', edgecolor='black', alpha=0.8)
    plot1.annotate(f'Median peaks width={median_len}bp', xy=(0.25, 0.8), xycoords='axes fraction', size=8)
    plot1.annotate(f'Total number of peaks={len(peaks_width)}', xy=(0.25, 0.7), xycoords='axes fraction', size=8)
    plot1.set_xlabel('Peaks width', size=16)
    plot1.set_ylabel('Number of peaks', size=16)
    plot1.set_yscale('log', nonposy='clip')
    plot1.set_title('Peaks width\n'+set_name, size=14)  
    #Peaks FE
    peaks_FE=set_data_ar[5]
    plot2=plt.subplot2grid((1,3),(0,2), rowspan=1, colspan=1)  
    bins_FE=np.linspace(Min_FE, Max_FE, 50)
    plot2.hist(peaks_FE, bins_FE, color='#ffff3a', edgecolor='black', alpha=0.8)
    plot2.annotate(f'Median peaks FE={round(np.median(peaks_FE),3)} bp', xy=(0.25, 0.8), xycoords='axes fraction', size=8)
    plot2.set_xlabel('Peaks fold enrichment', size=16)
    plot2.set_ylabel('Number of peaks', size=16)
    plot2.set_yscale('log', nonposy='clip')
    plot2.set_title('Peaks FE\n'+set_name, size=14)        
    plt.tight_layout()
    plt.show()
    plt.savefig(pathout, dpi=300, figsize=(10, 5))
    plt.close()    
    return

#######
#Plots distribution of peaks fold enrichment (FE) only (for shared peaks).
#######

def plot_distrib_FE(dataset_name, peaks_FE, Min_FE, Max_FE, pathout):
    #Plot FE distribution.
    fig=plt.figure(figsize=(3, 3), dpi=100)
    plot0=plt.subplot2grid((1,1),(0,0), rowspan=1, colspan=1)     
    bins_FE=np.linspace(Min_FE, Max_FE, 50)
    plot0.hist(peaks_FE, bins_FE, color='#ffff3a', edgecolor='black', alpha=0.8)
    plot0.annotate(f'Median\npeaks FE={round(np.median(peaks_FE),2)} bp', xy=(0.25, 0.8), xycoords='axes fraction', size=11)
    plot0.set_xlabel('Peaks fold enrichment', size=16)
    plot0.set_ylabel('Number of peaks', size=16)
    plot0.set_yscale('log', nonposy='clip')
    plot0.set_title('Peaks FE\n'+dataset_name, size=14)        
    plt.tight_layout()
    plt.show()
    plt.savefig(pathout, dpi=300, figsize=(10, 5))
    plt.close()       
    return

#######
#Wrapper: reads NarrowPeak input, return sequences by coordinates provided from reference genome,
#writes sequences under peaks to mfa file, plots distributions of GC% of sequences under the peaks and peaks widths.
#######

def wrap_peaks(fasta_inpath, peaks_inpath_dict, shared_peaks_inpath_dict, outpath):
    ##Reads genome fasta.
    genome=read_genome(fasta_inpath)
    ##Reads peaks data.
    Peaks_data_dict={}
    for dataset_name, data_path in peaks_inpath_dict.items():    
        peaks_coords, peaks_FE=deletions_info(data_path[0])
        dataset_file_name=data_path[1]
        #Retrive sequences under peaks.
        seqs=return_seqs(peaks_coords, genome)
        #Keep data.
        Peaks_data_dict[dataset_name]=[peaks_coords, dataset_file_name, seqs, peaks_FE]
        #Write sequences.
        write_seqs(seqs, outpath+dataset_file_name+".fasta")
    
    ##Analyse sets of peaks.
    #Genome binning and GC% calculation.
    #Identify limits of plots axis.
    all_peaks_width=[]
    all_gc_values=[]
    all_FE_values=[]
    GC_width_FE_dict={}
    for peaks_dataset_name, peaks_data in Peaks_data_dict.items():   
        seq_dict=peaks_data[2]
        GC_genome, median_len, genome_gc_values, peaks_gc_values, peaks_width=genome_bin_GC_dist(seq_dict, genome)
        GC_width_FE_dict[peaks_dataset_name]=[GC_genome, median_len, genome_gc_values, peaks_gc_values, peaks_width, peaks_data[3]]
        all_peaks_width+=peaks_width
        all_gc_values=all_gc_values + genome_gc_values + peaks_gc_values
        all_FE_values+=peaks_data[3]
        #Create DataFrame and write file with peaks coordinates, FE, width, and GC%.
        peaks_coords=peaks_data[0]
        starts, ends = map(list, zip(*peaks_coords))
        Data_dictionary={'Start' : starts, 'End' : ends, 'Width' : peaks_width, 'GC%' : peaks_gc_values, 'FE' : peaks_data[3]}
        Data_DF=pd.DataFrame(Data_dictionary)
        Data_DF.to_csv(outpath+peaks_data[1]+'_width_GC_FE.csv', sep='\t', index=False)

    Min_GC=min(all_gc_values)
    Max_GC=max(all_gc_values)
    Min_width=min(all_peaks_width)
    Max_width=max(all_peaks_width)
    Min_FE=min(all_FE_values)
    Max_FE=max(all_FE_values)
    
    #Plot peaks GC%, width, and FE distributions.
    for peaks_dataset_name, set_data_ar in GC_width_FE_dict.items():
        name=Peaks_data_dict[peaks_dataset_name][1]
        plot_distrib(peaks_dataset_name, set_data_ar, Min_GC, Max_GC, Min_width, Max_width, Min_FE, Max_FE, outpath+name+"_GC_width.png")
   
    ##Analyse shared peaks.
    ##Reads peaks data.
    Shared_peaks_data_dict={}
    shared_all_FE_values=[]
    del Data_dictionary['FE']
    for dataset_name, data_path in shared_peaks_inpath_dict.items():    
        peaks_coords, peaks_FE=deletions_info(data_path[0])
        dataset_file_name=data_path[1]
        #Keep data.
        Shared_peaks_data_dict[dataset_name]=[dataset_file_name, peaks_FE] 
        shared_all_FE_values+=peaks_FE
        #Create DataFrame and write file with peaks coordinates, FE, width, and GC%.
        Data_dictionary[dataset_name+'_FE']=peaks_FE
    
    Data_DF=pd.DataFrame(Data_dictionary)
    Data_DF.to_csv(outpath+'Shared_peaks_TopoA_nm_0.001_width_GC_FE_for_diffr_experiments.csv', sep='\t', index=False)
    
    Shared_Min_FE=min(shared_all_FE_values)
    Shared_Max_FE=max(shared_all_FE_values)   
    
    #Plot shared peaks FE distributions.
    for dataset_name, set_data_ar in Shared_peaks_data_dict.items():
        plot_distrib_FE(dataset_name, set_data_ar[1], Shared_Min_FE, Shared_Max_FE, outpath+set_data_ar[0]+"_FE.png")

    return


wrap_peaks(Genome, Peaks_dict, Shared_peaks_dict, Path_out)