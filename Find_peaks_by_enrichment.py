###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes wig file with fold enrichment, finds regions with enrichment > threshold,
#returns BroadPeaks file with peaks.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import Bio
from Bio import SeqIO
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles

#Dictionary of pathes to wig file with fold enrichment.
Peaks_data={'DRIP_Seq_CTD_minus_Rif_minus_forward' : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Strand_average_scaled\DRIP_Seq_CTD_minus_Rif_minus_forward_av.wig',
            'DRIP_Seq_CTD_minus_Rif_minus_reverse' : 'C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Strand_average_scaled\DRIP_Seq_CTD_minus_Rif_minus_reverse_av.wig'}
#Threshold for reproducible peaks calling (must not exceed number of replicas).
Threshold=int(300)
#Path to the reference genome (e.g. E_coli_w3110_G_Mu.fasta).
Genome="C:\Sutor\Science\TopoI-ChIP-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\E_coli_w3110_G_Mu.fasta"
#Outpath.
Path_out='C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Peaks\\'
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
#Parsing WIG file.
#######

def score_data_parser(inpath, param_name):
    param_file=open(inpath, 'r')
    ar=[]
    for line in param_file:
        line=line.rstrip().split(' ')
        if line[0]=='fixedStep':
            chrom_name=line[1].split('=')[1]      
        if line[0] not in ['track', 'fixedStep']:
            ar.append(float(line[0]))
    param_file.close()
    print('Whole genome average ' + str(param_name) + ' : ' + str(sum(ar)/len(ar)))
    return ar, chrom_name


#######
#Find peaks in fold enrichment array, write sequences under peaks..
#######  

def Find_rep_peaks(genome_ar, thr, genome, outpath):
    #Output file with sequences.
    Output_fasta=open(outpath, 'w')
    #Keep data about peaks.
    peak=0
    peak_cumul_len=0
    peak_num=0
    rep_peaks_ar=[]
    for i in range(len(genome_ar)):
        if genome_ar[i]<thr and peak==0: #We are not in peak.
            continue
        elif genome_ar[i]>=thr and peak==0: #We are at left peak border.
            peak=1
            current_peak=[i]
            continue
        elif genome_ar[i]>=thr and peak==1: #We are within a peak.
            continue
        elif genome_ar[i]<thr and peak==1: #We are at the right peak border.
            peak=0
            current_peak.append(i)
            #Report only long enough peaks.
            if (current_peak[1]-current_peak[0])>20:
                rep_peaks_ar.append(current_peak)
                peak_cumul_len+=(current_peak[1]-current_peak[0])
                peak_sequence=genome[current_peak[0]:current_peak[1]]
                Output_fasta.write('>' + str(peak_num) + '\n' + str(peak_sequence) + '\n')
                peak_num+=1
            continue
    Output_fasta.close()
    print(f'Number of peaks found: {len(rep_peaks_ar)}')
    print(f'Cumulative length of peaks: {peak_cumul_len}')
    return rep_peaks_ar
        
#######
#Write reproducible peaks in broadPeak format.
#######   

def write_bed(rep_peaks_ar, chrom_name, outpath):
    fileout=open(outpath, 'w')
    for i in range(len(rep_peaks_ar)):
        fileout.write(chrom_name+'\t'+str(rep_peaks_ar[i][0])+'\t'+str(rep_peaks_ar[i][1])+'\tPeak_'+str(i)+'\t10\t.\t-1\t-1\t-1\n')
    fileout.close()
    return
    
#######
#Wrapper: takes wig file with fold enrichment, 
#finds regions with enrichment > threshold now are called peaks,
#returns BroadPeaks file with peaks.
#######    

def Wrapper(genome_path, reps_dict, thr, outpath):
    #Read genome.
    genome_seq=read_genome(genome_path)
    #Read intervals data.
    rep_data={}
    for name, rep_path in reps_dict.items():
        #Reads FE data.
        FE_data, chrom_name=score_data_parser(rep_path, name)
        #Finds peaks.
        Rep_peaks_array=Find_rep_peaks(FE_data, thr, genome_seq, f'{outpath}{name}_peaks_threshold_{thr}_sequences_longer_20bp.fasta')
        rep_data[name]=Rep_peaks_array
        #Write reproducible peaks.
        write_bed(Rep_peaks_array, chrom_name, f'{outpath}{name}_peaks_threshold_{thr}_sequences_longer_20bp.BroadPeak')
    return
            
Wrapper(Genome, Peaks_data, Threshold, Path_out)       
