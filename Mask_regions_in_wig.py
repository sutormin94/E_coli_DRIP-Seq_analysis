###############################################
##Dmitry Sutormin, 2020##
##TopoA ChIP-Seq analysis##

#Takes wig file, masks regions (substitute values by 0).
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import Bio
from Bio import SeqIO
from matplotlib import pyplot as plt

#Path to the working directory.
PWD="C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\Strand_wise_FE_scaled\\"
#Dictionary of pathes to wig file with fold enrichment.
WIG_data={'DRIP_Seq_CTD_minus_Rif_minus_forward_av_FE' : PWD + 'DRIP_Seq_CTD_minus_Rif_minus_forward_av_FE.wig',
          'DRIP_Seq_CTD_minus_Rif_minus_reverse_av_FE' : PWD + 'DRIP_Seq_CTD_minus_Rif_minus_reverse_av_FE.wig',
          'DRIP_Seq_CTD_minus_Rif_plus_forward_av_FE' : PWD + 'DRIP_Seq_CTD_minus_Rif_plus_forward_av_FE.wig',
          'DRIP_Seq_CTD_minus_Rif_plus_reverse_av_FE' : PWD + 'DRIP_Seq_CTD_minus_Rif_plus_reverse_av_FE.wig',
          'DRIP_Seq_CTD_plus_Rif_minus_forward_av_FE' : PWD + 'DRIP_Seq_CTD_plus_Rif_minus_forward_av_FE.wig',
          'DRIP_Seq_CTD_plus_Rif_minus_reverse_av_FE' : PWD + 'DRIP_Seq_CTD_plus_Rif_minus_reverse_av_FE.wig',
          'DRIP_Seq_CTD_plus_Rif_plus_forward_av_FE' : PWD + 'DRIP_Seq_CTD_plus_Rif_plus_forward_av_FE.wig',
          'DRIP_Seq_CTD_plus_Rif_plus_reverse_av_FE' : PWD + 'DRIP_Seq_CTD_plus_Rif_plus_reverse_av_FE.wig'}
#List of regions to mask.
Masked=[[274500, 372148], [793800, 807500], [1199000, 1214000], [656600, 657860]]
#Name of a chromosome.
Chromosome_name="NC_007779.1_w3110_Mu"


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
##Mask regions.
#######

def Mask_regions(wig_ar, masked_list):
    for region in masked_list:
        zero_ar=[0]*(region[1]-region[0])
        print(len(wig_ar[region[0]:region[1]]), len(zero_ar))
        wig_ar[region[0]:region[1]]=zero_ar
    return wig_ar

#######
##Write wig.
#######

def write_wig(wig_data, Chromosome_name, name, outpath):
    #Write file.
    file_out=open(outpath, 'w')
    file_out.write('track type=wiggle_0 name="'+name+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
    for i in range(len(wig_data)):
        file_out.write(str(wig_data[i])+'\n')
    file_out.close()
    return

#######
##Wrapper.
#######

def Wrapper_function(wig_dict, masked_list, Chromosome_name, outpath):
    #Iterate wig files, mask regions.
    i=0
    for name, wig_path in wig_dict.items():
        wig_data=wig_parsing(wig_path)
        wig_data_masked=Mask_regions(wig_data, masked_list)
        write_wig(wig_data_masked, Chromosome_name, name+'_deletions_masked', outpath+name+'_deletions_masked.wig')
        print('Progress: ' + str(i) + '/' + str(len(wig_dict)))
        i+=1
    
    return

Wrapper_function(WIG_data, Masked, Chromosome_name, PWD)