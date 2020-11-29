###############################################
##Dmitry Sutormin, 2020##
##DRIP-Seq analysis##

####
#To subtract "+" and "-" strands of DRIP-Seq experiment (or any other strand-specific experiments).
####

###############################################

import numpy as np

#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\\"
#Path to the file with IP data
IP_path_dict={'1' :  PWD +        "WIG_av_str_subtr_scaled\DRIP_Seq_CTD_minus_Rif_minus_av_strands_subtr_scaled.wig",
              '2' :  PWD +        "WIG_av_str_subtr_scaled\DRIP_Seq_CTD_minus_Rif_plus_av_strands_subtr_scaled.wig",  
              '3' :  PWD +        "WIG_av_str_subtr_scaled\DRIP_Seq_CTD_plus_Rif_minus_av_strands_subtr_scaled.wig",
              '4' :  PWD +        "WIG_av_str_subtr_scaled\DRIP_Seq_CTD_plus_Rif_plus_av_strands_subtr_scaled.wig", 
              }

#Path to the file Mock control data
Mock_path_dict={'1' :  PWD +           "WIG_av_str_subtr_scaled\DRIP_Seq_CTD_minus_Rif_minus_RNAse_av_strands_subtr_scaled.wig",
                '2' :  PWD +           "WIG_av_str_subtr_scaled\DRIP_Seq_CTD_minus_Rif_plus_RNAse_av_strands_subtr_scaled.wig",
                '3' :  PWD +           "WIG_av_str_subtr_scaled\DRIP_Seq_CTD_plus_Rif_minus_RNAse_av_strands_subtr_scaled.wig",
                '4' :  PWD +           "WIG_av_str_subtr_scaled\DRIP_Seq_CTD_plus_Rif_plus_RNAse_av_strands_subtr_scaled.wig",                 
                }


#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  "DRIP_Seq_CTD_minus_Rif_minus_av_strands_subtr_scaled_control_subtr",
           '2' :  "DRIP_Seq_CTD_minus_Rif_plus_av_strands_subtr_scaled_control_subtr", 
           '3' :  "DRIP_Seq_CTD_plus_Rif_minus_av_strands_subtr_scaled_control_subtr",
           '4' :  "DRIP_Seq_CTD_plus_Rif_plus_av_strands_subtr_scaled_control_subtr", 
           }

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name_manual=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)

#Path to the output directory.
PWD_out=PWD + "WIG_av_str_subtr_scaled_control_subtr\\"
#Output path to the final file (fold enrichment).
FE_file_path_dict={'1' :  PWD_out +           "DRIP_Seq_CTD_minus_Rif_minus_av_strands_subtr_scaled_control_subtr.wig",
                   '2' :  PWD_out +           "DRIP_Seq_CTD_minus_Rif_plus_av_strands_subtr_scaled_control_subtr.wig",
                   '3' :  PWD_out +           "DRIP_Seq_CTD_plus_Rif_minus_av_strands_subtr_scaled_control_subtr.wig",
                   '4' :  PWD_out +           "DRIP_Seq_CTD_plus_Rif_plus_av_strands_subtr_scaled_control_subtr.wig", 
                   }


#######
#Parses WIG file.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    Dict_of_chromosomes_data={}
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0]=='fixedStep':
            chrom_name=line[1].split('=')[1]
            Dict_of_chromosomes_data[chrom_name]=[]
        if line[0] not in ['track', 'fixedStep']:
            Dict_of_chromosomes_data[chrom_name].append(float(line[0]))
    wigin.close()
    
    for Chromosome_name, data in Dict_of_chromosomes_data.items():
        data_array=np.array(data)
        data_mean=np.mean(data_array)
        print(f'Mean coverage of {Chromosome_name}: {data_mean}')
        #data_array_scaled=data_array/data_mean
        #Dict_of_chromosomes_data[Chromosome_name]=data_array_scaled
    return Dict_of_chromosomes_data


def read_files(input_dict):
    Data_dict={}
    for name, path in input_dict.items():
        Data_dict[name]=wig_parsing(path)
        print(f'Progress: {name}/{len(input_dict)}')
    return Data_dict

IP_dict=read_files(IP_path_dict)
Mock_dict=read_files(Mock_path_dict)


def subtract_write(IP_dict, Mock_dict, name_dict, Auto_or_manual, Chromosome_name_manual, FE_file_path_dict):
    for sample_name, sample_data in IP_dict.items():
        print(f'Now is processing: {sample_name}')
        print(f'Progress: {sample_name}/{len(IP_dict)}')
        FE_out=open(FE_file_path_dict[sample_name], 'w')
        #Write file with fold enrichment data.
        for Chromosome_name, data in sample_data.items():
            print(f'Average covarage of (+) strand: {np.mean(data)}')
            print(f'Average covarage of (-) strand: {np.mean(Mock_dict[sample_name][Chromosome_name])}')
            if Auto_or_manual==0:
                FE_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
            elif Auto_or_manual==1:
                FE_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name_manual+' start=1 step=1\n')
            for i in range(len(data)):
                FE_data_position=(IP_dict[sample_name][Chromosome_name][i])-(Mock_dict[sample_name][Chromosome_name][i])
                FE_out.write(str(FE_data_position)+'\n')
            
        FE_out.close()        
    return

subtract_write(IP_dict, Mock_dict, name_dict, Auto_or_manual, Chromosome_name_manual, FE_file_path_dict)
