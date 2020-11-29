###############################################
##Dmitry Sutormin, 2020##
##DRIP-Seq analysis##

####
#The only purpose - to subtract "+" and "-" strands of DRIP-Seq experiment (or any other strand-specific experiments).
####

###############################################

import numpy as np

#Path to the working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\\"
#Path to the files with data to be scaled.
IN_path_dict={'1' :  PWD +        "WIG_av_str_subtr\DRIP_Seq_CTD_minus_Rif_minus_av_strands_subtr.wig",
              '2' :  PWD +        "WIG_av_str_subtr\DRIP_Seq_CTD_minus_Rif_minus_RNAse_av_strands_subtr.wig",  
              '3' :  PWD +        "WIG_av_str_subtr\DRIP_Seq_CTD_minus_Rif_plus_av_strands_subtr.wig",
              '4' :  PWD +        "WIG_av_str_subtr\DRIP_Seq_CTD_minus_Rif_plus_RNAse_av_strands_subtr.wig",             
              }

#Dict with scaling coefficients.
#Scaling_dict={'1' : 0.782708832, '2' :  0.884843949,  '3' :  0.659649916, '4' :  0.5726933, '5' : 0.545898119, '6' :  0.524778634,
#              '7' : 0.739066282, '8' :  0.787327564,  '9' :  0.644613386, '10' : 0.596754035, '11' : 0.535745327,
#              '12' : 0.707823289, '13' : 0.814711565, '14' : 1.0, '15' : 0.616383167, '16' : 0.527514643, '17' : 0.538420646,
#              '18' : 0.615526028, '19' : 0.672435257, '20' : 0.543809175, '21' : 0.740944617, '22' : 0.573905487, '23' : 0.547143087}

#Scaling_dict={'1' : 0.770722056, '2' :  0.728681254,  '3' :  0.547790018, '4' :  0.566249681}

Scaling_dict={'1' : 0.840844951, '2' :  0.610590153,  '3' :  0.560772819, '4' :  0.620664397}

#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  "DRIP_Seq_CTD_minus_Rif_minus_av_strands_subtr_scaled",
           '2' :  "DRIP_Seq_CTD_minus_Rif_minus_RNAse_av_strands_subtr_scaled",
           '3' :  "DRIP_Seq_CTD_minus_Rif_plus_av_strands_subtr_scaled",
           '4' :  "DRIP_Seq_CTD_minus_Rif_plus_RNAse_av_strands_subtr_scaled",             
           }

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name_manual=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)

#Path to the output directory.
PWD_out=PWD + "WIG_av_str_subtr_scaled\\"
#Output path to the final file (fold enrichment).
FE_file_path_dict={'1' :  PWD_out + "DRIP_Seq_CTD_minus_Rif_minus_av_strands_subtr_scaled",
                   '2' :  PWD_out + "DRIP_Seq_CTD_minus_Rif_minus_RNAse_av_strands_subtr_scaled",
                   '3' :  PWD_out + "DRIP_Seq_CTD_minus_Rif_plus_av_strands_subtr_scaled",
                   '4' :  PWD_out + "DRIP_Seq_CTD_minus_Rif_plus_RNAse_av_strands_subtr_scaled",                    
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
        Dict_of_chromosomes_data[Chromosome_name]=data_array
    return Dict_of_chromosomes_data


def read_files(input_dict):
    Data_dict={}
    for name, path in input_dict.items():
        Data_dict[name]=wig_parsing(path)
        print(f'Progress: {name}/{len(input_dict)}')
    return Data_dict

IN_dict=read_files(IN_path_dict)

def scale_write(IN_dict, Scale_dict, name_dict, Auto_or_manual, Chromosome_name_manual, FE_file_path_dict):
    for sample_name, sample_data in IN_dict.items():
        print(f'Now is processing: {sample_name}')
        print(f'Progress: {sample_name}/{len(IN_dict)}')
        FE_for_out=open(FE_file_path_dict[sample_name]+'.wig', 'w')
        #Write file with scaled data.
        for Chromosome_name, data in sample_data.items():
            print(f'Average covarage before scaling: {np.mean(data)}')
            For_data_scaled=np.array(data)*Scale_dict[sample_name]
            print(f'Average covarage of after scaling: {np.mean(For_data_scaled)}')            
            if Auto_or_manual==0:
                FE_for_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
            elif Auto_or_manual==1:
                FE_for_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name_manual+' start=1 step=1\n')
            for i in range(len(data)):
                FE_for_out.write(str(For_data_scaled[i])+'\n')
        FE_for_out.close()      
    return

scale_write(IN_dict, Scaling_dict, name_dict, Auto_or_manual, Chromosome_name_manual, FE_file_path_dict)
