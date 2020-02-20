###############################################
##Dmitry Sutormin, 2020##
##TopoA ChIP-Seq analysis##

#Takes list of regions (BroadPeak, bed), return sequences of them,
#Identify k-meres and compare k-mere frequencis for regions against the whole genome sequence.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import Bio
from Bio import SeqIO
from matplotlib import pyplot as plt

#Path to file with regions.
Regions_path="C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Peaks\DRIP_Seq_CTD_minus_Rif_minus_forward_peaks_threshold_300_sequences_longer_20bp_no_rRNA_LacI.BroadPeak"
#Path to file with genome sequence.
Genome_path="C:\Sutor\Science\E_coli_DRIP-Seq\Scripts\Additional_genome_features\E_coli_w3110_G_Mu.fasta"
#Output path.
Output_path="C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Peaks_k_meres\\"
if not os.path.exists(Output_path):
    os.makedirs(Output_path)
    
    
    
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
#Return sequences for a set of regions.
#######

def return_regions_count_k_mers(regions_ar, genome_seq, K_mer_dict, k):
    regions_seq_ar=[]
    for region in regions_ar:
        region_seq=genome_seq[region[0]:region[1]]
        regions_seq_ar.append(region_seq)
        K_mer_dict=naive_k_mer_identification(region_seq, K_mer_dict, k)
    return K_mer_dict, regions_seq_ar

#######
#Naive methods for K-mers identification and counting.
#######

def naive_k_mer_identification(sequence, K_mer_dict, k):
    for i in range(len(sequence)-k):
        kmer=sequence[i:i+k]
        if kmer not in K_mer_dict:
            K_mer_dict[kmer]=1
        else:
            K_mer_dict[kmer]+=1
    return K_mer_dict

#######
#Calculate fraction of k-mers.
#######

def calc_k_mers_freq(k_mers_dict):
    #Calculate cumulative frequency of k-mers.
    cumul_freq=0
    for k_mer, freq in k_mers_dict.items():
        cumul_freq+=freq
    print('Total number of k-mers: ' + str(cumul_freq))
    #Calculate fraction of k-mers.
    for k_mer, freq in k_mers_dict.items():
        frequency=k_mers_dict[k_mer]
        k_mers_dict[k_mer]=frequency/float(cumul_freq)
    return k_mers_dict

#######
#Visualize k-mers frequency.
#######

def k_mers_frequency_distribution(k_mer_dict_1, k_mer_dict_2, name, output_path):
    #Sort k-mers by frequency.
    k_mer_dict_1_sorted={k: v for k, v in sorted(k_mer_dict_1.items(), key=lambda item: item[1], reverse=True)}
    #Keep k-mers changed their frequency dramatically.
    fileout_incr=open(output_path+name+'_k_mers_increased_freq.txt', 'w')
    fileout_decr=open(output_path+name+'_k_mers_decreased_freq.txt', 'w')
    fileout_incr.write('K-mer\tGenome frequency\tRegions frequency\tRatio of frequences\tRank for Genome\n')
    fileout_decr.write('K-mer\tGenome frequency\tRegions frequency\tRatio of frequences\tRank for Genome\n')
    #Convert to arrays.
    k_mers_ar_1=[]
    freq_ar_1=[]
    x_ar=[]
    freq_ar_2=[]
    freq_ratio=[]
    i=0
    for k_mer, freq in k_mer_dict_1_sorted.items():
        #For genome.
        k_mers_ar_1.append(k_mer)
        freq_ar_1.append(freq)  
        x_ar.append(i)
        #For regions.
        if k_mer in k_mer_dict_2:
            freq_ar_2.append(k_mer_dict_2[k_mer])
            kmers_freq_ratio=float(k_mer_dict_2[k_mer])/freq
            freq_ratio.append(kmers_freq_ratio)
            if kmers_freq_ratio>3:
                print(str(k_mer) + ' of rank ' + str(i) + ' out of ' + str(len(k_mer_dict_1_sorted)) + '; k-mers with frequency of regions ' + 
                      str(k_mer_dict_2[k_mer]) + ' and k-mers frequency in a genome ' + str(freq))
                fileout_incr.write(str(k_mer) + '\t' + str(freq) + '\t' + str(k_mer_dict_2[k_mer]) + '\t' + str(kmers_freq_ratio) + '\t' + str(i)+'/'+str(len(k_mer_dict_1_sorted)) + '\n')
            if kmers_freq_ratio<float(1/3):
                fileout_decr.write(str(k_mer) + '\t' + str(freq) + '\t' + str(k_mer_dict_2[k_mer]) + '\t' + str(kmers_freq_ratio) + '\t' + str(i)+'/'+str(len(k_mer_dict_1_sorted)) + '\n')
        else:
            freq_ar_2.append(-100) 
            freq_ratio.append(0)
        i+=1
    fileout_incr.close()
    fileout_decr.close()
    #Plot data.
    fig, plot_av=plt.subplots(1,1,figsize=(12,5), dpi=100)
    plot_av.plot(x_ar, freq_ar_1, 'ob', markersize=10, label='Genome k-mers frequences')
    plot_av.plot(x_ar, freq_ar_2, 'or', markersize=1, label='Region k-mers frequences')
    plot_av.plot(x_ar, freq_ratio, 'og', markersize=1, label='Ratio of k-mers frequences')
    plot_av.hlines(3, 0, len(x_ar), linestyles='dashed', linewidth=0.4)
    plot_av.hlines(1, 0, len(x_ar), linewidth=1)
    plot_av.hlines(float(1)/3, 0, len(x_ar), linestyles='dashed', linewidth=0.4)
    plot_av.set_xticks(x_ar, minor=False)
    plot_av.set_xticklabels(k_mers_ar_1, rotation=90, size=4)
    #plot_av.set_ylim([0,max(freq_ar_1+freq_ar_2+freq_ratio)])
    plot_av.set_yscale('log')
    plot_av.legend(loc='center right')  
    plt.show()
    plt.savefig(output_path+name+'_frequency_analysis.png', dpi=300, size=(15,5))
    return

#######
#Wrapper function.
#######

def wrapper_func(regions_path, genome_path, k_len, output_path):
    #Read genome.
    genome_sequence=read_genome(genome_path)
    #Read list of regions (peaks).
    regions_coords_ar, regions_signal_ar=deletions_info(regions_path)
    #Identify and count k-mers for genome.
    Genome_k_mers={}
    Genome_k_mers=naive_k_mer_identification(genome_sequence, Genome_k_mers, k_len)
    print('Number of different k-mers of length ' + str(k_len) + ' in a genome ' + str(len(Genome_k_mers)))
    #Identify and count k-mers for genomic regions.
    Regions_k_mers={}
    Regions_k_mers, Regions_seq_ar=return_regions_count_k_mers(regions_coords_ar, genome_sequence, Regions_k_mers, k_len)
    print('Number of different k-mers of length ' + str(k_len) + ' in regions ' + str(len(Regions_k_mers)))
    #Calculate fraction of k-mers.
    Genome_k_mers_fraction=calc_k_mers_freq(Genome_k_mers)
    Regions_k_mers_fraction=calc_k_mers_freq(Regions_k_mers) 
    #Plot the data.
    k_mers_frequency_distribution(Genome_k_mers_fraction, Regions_k_mers_fraction, 'DRIP_CTD_minus_Rif_minus_300_fw_sequences_longer_20bp_no_rRNA_LacI_k_mers_'+str(k_len), output_path)
    return

wrapper_func(Regions_path, Genome_path, 6, Output_path)