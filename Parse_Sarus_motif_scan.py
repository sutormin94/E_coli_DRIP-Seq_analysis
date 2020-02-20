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

#######
#Variables to be defined.
#######

#Input: SARUS output, containing all scores for all positions in forward and reverse orientation.
Input_SARUS_output="F:\TopoI_ChIP-Seq\Ec_TopoI_data\Motif_scanning\ChIP-Munk\\Rif_Rep12_thr_0.001\\results.txt"

#Output path.
Output_path="F:\TopoI_ChIP-Seq\Ec_TopoI_data\Motif_scanning\ChIP-Munk\Rif_Rep12_thr_0.001\\"

#Name of a dataset.
Dataset_name="EcTopoI_Rif_motif_w3110_scanned"


#######
#Reads SARUS output and classify data on "+" (plus) strand, "-" (minus) strand and max of both.
#######

def read_sarus_parse(sarus_output):
    filein=open(sarus_output, 'r')
    ar_plus=[]
    ar_minus=[]
    ar_max_score=[]
    dict_plus_minus={}
    for line in filein:
        if line[0]=='>':
            line=line.lstrip('>').rstrip().split(' ')
            seq_id=line[0]
        else:
            line=line.rstrip().split('\t')
            score=line[0]
            position=line[1]
            strand=line[2]
            if strand=='+':
                ar_plus.append(line)
            elif strand=='-':
                ar_minus.append(line)
            if line[1] not in dict_plus_minus:
                dict_plus_minus[line[1]]=[line]
            else:
                dict_plus_minus[line[1]].append(line)
                max_score=max(dict_plus_minus[line[1]][0][0], dict_plus_minus[line[1]][1][0])
                ar_max_score.append([max_score, line[1], '.'])
    filein.close()
    print(f'File was parsed succesfully')
    return seq_id, ar_plus, ar_minus, ar_max_score


#######
#Write wig file.
#######

def write_wig(output_path, data, dataset_name, seq_id):
    wig_out=open(output_path, 'w')
    wig_out.write(f'track type=wiggle_0 name="{dataset_name}" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom={seq_id} start=1 step=1\n')
    for i in range(len(data)):
        wig_out.write(f'{data[i][0]}\n')    
    wig_out.close()
    return


#######
#Wrapper function.
#######

def wrapper(sarus_output, output_path, dataset_name):
    seq_id, ar_plus, ar_minus, ar_max_score=read_sarus_parse(sarus_output)
    write_wig(f'{output_path}{dataset_name}_plus.wig', ar_plus, f'{dataset_name}_plus', seq_id)
    print(f'Plus is done!')
    write_wig(f'{output_path}{dataset_name}_minus.wig', ar_minus, f'{dataset_name}_minus', seq_id)
    print(f'Minus is done!')
    write_wig(f'{output_path}{dataset_name}_both.wig', ar_max_score, f'{dataset_name}_both', seq_id)
    print(f'Both is done!')
    return

wrapper(Input_SARUS_output, Output_path, Dataset_name)