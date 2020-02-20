###############################################
##Dmitry Sutormin, 2020##
##DRIP-Seq analysis##

####
#The only purpose - to subtract "+" and "-" strands of DRIP-Seq experiment (or any other strand-specific experiments).
####

###############################################

import numpy as np

#Path to the working directory.
PWD="C:\Sutor\Science\E_coli_DRIP-Seq\WIG\\"
#Path to the file with IP data
IP_path_dict={'1' :  PWD +        "DRIP_CTD_1_S13_edt_forward_depth.wig",
              '2' :  PWD +        "DRIP_CTD_2_S15_edt_forward_depth.wig",  
              '3' :  PWD +        "DRIP_CTD_3_S17_edt_forward_depth.wig",
              '4' :  PWD +        "DRIP_CTD_r1_S19_edt_forward_depth.wig", 
              '5' :  PWD +        "DRIP_CTD_r2_S20_edt_forward_depth.wig",
              '6' :  PWD +        "DRIP_CTD_r3_S22_edt_forward_depth.wig",
              '7' :  PWD +        "DRIP_CTD_RNAse_1_S14_edt_forward_depth.wig",
              '8' :  PWD +        "DRIP_CTD_RNAse_2_S16_edt_forward_depth.wig",
              '9' :   PWD +       "DRIP_CTD_RNAse_3_S18_edt_forward_depth.wig",
              '10' :  PWD +       "DRIP_CTD_RNAse_r2_S21_edt_forward_depth.wig",  
              '11' :  PWD +       "DRIP_CTD_RNAse_r3_S23_edt_forward_depth.wig",
              '12' :  PWD +       "TopA-SPA_DRIP-Seq_1_S1_edt_forward_depth.wig", 
              '13' :  PWD +       "TopA-SPA_DRIP-Seq_2_S3_edt_forward_depth.wig",
              '14' :  PWD +       "TopA-SPA_DRIP-Seq_3_S5_edt_forward_depth.wig",
              '15' :  PWD +       "TopA-SPA_DRIP-Seq_1rif_S7_edt_forward_depth.wig",
              '16' :  PWD +       "TopA-SPA_DRIP-Seq_2rif_S9_edt_forward_depth.wig", 
              '17' :  PWD +       "TopA-SPA_DRIP-Seq_3rif_S11_edt_forward_depth.wig", 
              '18' :  PWD +       "TopA-SPA_DRIP-Seq_1H_S2_edt_forward_depth.wig",
              '19' :  PWD +       "TopA-SPA_DRIP-Seq_2H_S4_edt_forward_depth.wig", 
              '20' :  PWD +       "TopA-SPA_DRIP-Seq_3H_S6_edt_forward_depth.wig",
              '21' :  PWD +       "TopA-SPA_DRIP-Seq_1rH_S8_edt_forward_depth.wig",
              '22' :  PWD +       "TopA-SPA_DRIP-Seq_2rH_S10_edt_forward_depth.wig",
              '23' :  PWD +       "TopA-SPA_DRIP-Seq_3rH_S12_edt_forward_depth.wig",               
              }

#Path to the file Mock control data
Mock_path_dict={'1' :  PWD +      "DRIP_CTD_1_S13_edt_reverse_depth.wig",
                '2' :  PWD +      "DRIP_CTD_2_S15_edt_reverse_depth.wig",
                '3' :  PWD +      "DRIP_CTD_3_S17_edt_reverse_depth.wig",
                '4' :  PWD +      "DRIP_CTD_r1_S19_edt_reverse_depth.wig",   
                '5' :  PWD +      "DRIP_CTD_r2_S20_edt_reverse_depth.wig",
                '6' :  PWD +      "DRIP_CTD_r3_S22_edt_reverse_depth.wig",
                '7' :  PWD +      "DRIP_CTD_RNAse_1_S14_edt_reverse_depth.wig",
                '8' :  PWD +      "DRIP_CTD_RNAse_2_S16_edt_reverse_depth.wig",  
                '9' :   PWD +     "DRIP_CTD_RNAse_3_S18_edt_reverse_depth.wig",
                '10' :  PWD +     "DRIP_CTD_RNAse_r2_S21_edt_reverse_depth.wig",  
                '11' :  PWD +     "DRIP_CTD_RNAse_r3_S23_edt_reverse_depth.wig",
                '12' :  PWD +     "TopA-SPA_DRIP-Seq_1_S1_edt_reverse_depth.wig", 
                '13' :  PWD +     "TopA-SPA_DRIP-Seq_2_S3_edt_reverse_depth.wig",
                '14' :  PWD +     "TopA-SPA_DRIP-Seq_3_S5_edt_reverse_depth.wig",
                '15' :  PWD +     "TopA-SPA_DRIP-Seq_1rif_S7_edt_reverse_depth.wig",
                '16' :  PWD +     "TopA-SPA_DRIP-Seq_2rif_S9_edt_reverse_depth.wig",  
                '17' :  PWD +     "TopA-SPA_DRIP-Seq_3rif_S11_edt_reverse_depth.wig",
                '18' :  PWD +     "TopA-SPA_DRIP-Seq_1H_S2_edt_reverse_depth.wig",
                '19' :  PWD +     "TopA-SPA_DRIP-Seq_2H_S4_edt_reverse_depth.wig", 
                '20' :  PWD +     "TopA-SPA_DRIP-Seq_3H_S6_edt_reverse_depth.wig",
                '21' :  PWD +     "TopA-SPA_DRIP-Seq_1rH_S8_edt_reverse_depth.wig",
                '22' :  PWD +     "TopA-SPA_DRIP-Seq_2rH_S10_edt_reverse_depth.wig",
                '23' :  PWD +     "TopA-SPA_DRIP-Seq_3rH_S12_edt_reverse_depth.wig",                
                }

#Dict with scaling coefficients.
Scaling_dict={'1' : 0.782708832, '2' :  0.884843949,  '3' :  0.659649916, '4' :  0.5726933, '5' : 0.545898119, '6' :  0.524778634,
              '7' : 0.739066282, '8' :  0.787327564,  '9' :  0.644613386, '10' : 0.596754035, '11' : 0.535745327,
              '12' : 0.707823289, '13' : 0.814711565, '14' : 1.0, '15' : 0.616383167, '16' : 0.527514643, '17' : 0.538420646,
              '18' : 0.615526028, '19' : 0.672435257, '20' : 0.543809175, '21' : 0.740944617, '22' : 0.573905487, '23' : 0.547143087}

#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  "DRIP_Seq_CTD+/Rif-_Rep1_S13",
           '2' :  "DRIP_Seq_CTD+/Rif-_Rep2_S15",
           '3' :  "DRIP_Seq_CTD+/Rif-_Rep3_S17",
           '4' :  "DRIP_Seq_CTD+/Rif+_Rep1_S19",   
           '5' :  "DRIP_Seq_CTD+/Rif+_Rep2_S20",
           '6' :  "DRIP_Seq_CTD+/Rif+_Rep3_S22",
           '7' :  "DRIP_Seq_CTD+/Rif-_RNAseA_Rep1_S14",
           '8' :  "DRIP_Seq_CTD+/Rif-_RNAseA_Rep2_S16",  
           '9' :  "DRIP_Seq_CTD+/Rif-_RNAseA_Rep3_S18",
           '10' : "DRIP_Seq_CTD+/Rif+_RNAseA_Rep2_S21",  
           '11' : "DRIP_Seq_CTD+/Rif+_RNAseA_Rep3_S23",
           '12' : "DRIP_Seq_CTD-/Rif-_Rep1_S1", 
           '13' : "DRIP_Seq_CTD-/Rif-_Rep2_S3",
           '14' : "DRIP_Seq_CTD-/Rif-_Rep3_S5",
           '15' : "DRIP_Seq_CTD-/Rif+_Rep1_S7",
           '16' : "DRIP_Seq_CTD-/Rif+_Rep2_S9",  
           '17' : "DRIP_Seq_CTD-/Rif+_Rep3_S11",
           '18' : "DRIP_Seq_CTD-/Rif-_RNAseH_Rep1_S2",
           '19' : "DRIP_Seq_CTD-/Rif-_RNAseH_Rep2_S4", 
           '20' : "DRIP_Seq_CTD-/Rif-_RNAseH_Rep3_S6",
           '21' : "DRIP_Seq_CTD-/Rif+_RNAseH_Rep1_S8",
           '22' : "DRIP_Seq_CTD-/Rif+_RNAseH_Rep2_S10",
           '23' : "DRIP_Seq_CTD-/Rif+_RNAseH_Rep3_S12",            
           }

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name_manual=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)

#Path to the output directory.
PWD_out="C:\Sutor\Science\E_coli_DRIP-Seq\Data_analysis\WIG_scaled\\"
#Output path to the final file (fold enrichment).
FE_file_path_dict={'1' :  PWD_out + "DRIP_Seq_CTD_plus_Rif_minus_Rep1_S13",
                   '2' :  PWD_out + "DRIP_Seq_CTD_plus_Rif_minus_Rep2_S15",
                   '3' :  PWD_out + "DRIP_Seq_CTD_plus_Rif_minus_Rep3_S17",
                   '4' :  PWD_out + "DRIP_Seq_CTD_plus_Rif_plus_Rep1_S19",   
                   '5' :  PWD_out + "DRIP_Seq_CTD_plus_Rif_plus_Rep2_S20",
                   '6' :  PWD_out + "DRIP_Seq_CTD_plus_Rif_plus_Rep3_S22",
                   '7' :  PWD_out + "DRIP_Seq_CTD_plus_Rif_minus_RNAseA_Rep1_S14",
                   '8' :  PWD_out + "DRIP_Seq_CTD_plus_Rif_minus_RNAseA_Rep2_S16",  
                   '9' :  PWD_out + "DRIP_Seq_CTD_plus_Rif_minus_RNAseA_Rep3_S18",
                   '10' : PWD_out + "DRIP_Seq_CTD_plus_Rif_plus_RNAseA_Rep2_S21",  
                   '11' : PWD_out + "DRIP_Seq_CTD_plus_Rif_plus_RNAseA_Rep3_S23",
                   '12' : PWD_out + "DRIP_Seq_CTD_minus_Rif_minus_Rep1_S1", 
                   '13' : PWD_out + "DRIP_Seq_CTD_minus_Rif_minus_Rep2_S3",
                   '14' : PWD_out + "DRIP_Seq_CTD_minus_Rif_minus_Rep3_S5",
                   '15' : PWD_out + "DRIP_Seq_CTD_minus_Rif_plus_Rep1_S7",
                   '16' : PWD_out + "DRIP_Seq_CTD_minus_Rif_plus_Rep2_S9",  
                   '17' : PWD_out + "DRIP_Seq_CTD_minus_Rif_plus_Rep3_S11",
                   '18' : PWD_out + "DRIP_Seq_CTD_minus_Rif_minus_RNAseH_Rep1_S2",
                   '19' : PWD_out + "DRIP_Seq_CTD_minus_Rif_minus_RNAseH_Rep2_S4", 
                   '20' : PWD_out + "DRIP_Seq_CTD_minus_Rif_minus_RNAseH_Rep3_S6",
                   '21' : PWD_out + "DRIP_Seq_CTD_minus_Rif_plus_RNAseH_Rep1_S8",
                   '22' : PWD_out + "DRIP_Seq_CTD_minus_Rif_plus_RNAseH_Rep2_S10",
                   '23' : PWD_out + "DRIP_Seq_CTD_minus_Rif_plus_RNAseH_Rep3_S12",                 
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

IP_dict=read_files(IP_path_dict)
Mock_dict=read_files(Mock_path_dict)


def scale_write(IP_dict, Mock_dict, Scale_dict, name_dict, Auto_or_manual, Chromosome_name_manual, FE_file_path_dict):
    for sample_name, sample_data in IP_dict.items():
        print(f'Now is processing: {sample_name}')
        print(f'Progress: {sample_name}/{len(IP_dict)}')
        FE_for_out=open(FE_file_path_dict[sample_name]+'_forward.wig', 'w')
        FE_rev_out=open(FE_file_path_dict[sample_name]+'_reverse.wig', 'w')
        #Write file with scaled data.
        for Chromosome_name, data in sample_data.items():
            print(f'Average covarage of (+) strand before scaling: {np.mean(data)}')
            print(f'Average covarage of (-) strand before scaling: {np.mean(Mock_dict[sample_name][Chromosome_name])}')
            For_data_scaled=data*Scale_dict[sample_name]
            Rev_data_scaled=Mock_dict[sample_name][Chromosome_name]*Scale_dict[sample_name]
            print(f'Average covarage of (+) strand after scaling: {np.mean(For_data_scaled)}')
            print(f'Average covarage of (-) strand after scaling: {np.mean(Rev_data_scaled)}')            
            if Auto_or_manual==0:
                FE_for_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'_forward" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
                FE_rev_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'_reverse" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
            elif Auto_or_manual==1:
                FE_for_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'_forward" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name_manual+' start=1 step=1\n')
                FE_rev_out.write('track type=wiggle_0 name="'+name_dict[sample_name]+'_reverse" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name_manual+' start=1 step=1\n')
            for i in range(len(data)):
                FE_for_out.write(str(For_data_scaled[i])+'\n')
                FE_rev_out.write(str(Rev_data_scaled[i])+'\n')   
        FE_for_out.close()    
        FE_rev_out.close()    
    return

scale_write(IP_dict, Mock_dict, Scaling_dict, name_dict, Auto_or_manual, Chromosome_name_manual, FE_file_path_dict)
