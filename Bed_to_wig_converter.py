###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - convert bed-like file to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the working directory.
PWD="C:\Sutor\Science\E_coli_ChIP-Seqs\Singh_CRP\\"
#Path to the input file
filein_path_dict={'1' :  PWD + "Cov_depth\SRR5099114.bed",
                  '2' :  PWD + "Cov_depth\SRR5099115.bed",   
                  '3' :  PWD + "Cov_depth\SRR5099116.bed",
                  '4' :  PWD + "Cov_depth\SRR5099117.bed",    
                  '5' :  PWD + "Cov_depth\SRR5099118.bed",
                  '6' :  PWD + "Cov_depth\SRR5099119.bed",  
                  '7' :  PWD + "Cov_depth\SRR5099120.bed",
                  '8' :  PWD + "Cov_depth\SRR5099121.bed",  
                  '9' :   PWD + "Cov_depth_nodup\SRR5099114_name_sorted_fm_ps_nd.bed",
                  '10' :  PWD + "Cov_depth_nodup\SRR5099115_name_sorted_fm_ps_nd.bed",   
                  '11' :  PWD + "Cov_depth_nodup\SRR5099116_name_sorted_fm_ps_nd.bed",
                  '12' :  PWD + "Cov_depth_nodup\SRR5099117_name_sorted_fm_ps_nd.bed",    
                  '13' :  PWD + "Cov_depth_nodup\SRR5099118_name_sorted_fm_ps_nd.bed",
                  '14' :  PWD + "Cov_depth_nodup\SRR5099119_name_sorted_fm_ps_nd.bed",  
                  '15' :  PWD + "Cov_depth_nodup\SRR5099120_name_sorted_fm_ps_nd.bed",
                  '16' :  PWD + "Cov_depth_nodup\SRR5099121_name_sorted_fm_ps_nd.bed",                  
                  }

#Path to the output file.
fileout_path_dict={'1' :  PWD + "WIG\Singh_CRP_EE_IP_Rep1_SRR5099114.wig",
                   '2' :  PWD + "WIG\Singh_CRP_EE_IP_Rep2_SRR5099115.wig",   
                   '3' :  PWD + "WIG\Singh_CRP_EE_Mock_Rep1_SRR5099116.wig",
                   '4' :  PWD + "WIG\Singh_CRP_EE_Mock_Rep2_SRR5099117.wig",   
                   '5' :  PWD + "WIG\Singh_CRP_ME_IP_Rep1_SRR5099118.wig",
                   '6' :  PWD + "WIG\Singh_CRP_ME_IP_Rep2_SRR5099119.wig",   
                   '7' :  PWD + "WIG\Singh_CRP_ME_Mock_Rep1_SRR5099120.wig",
                   '8' :  PWD + "WIG\Singh_CRP_ME_Mock_Rep2_SRR5099121.wig",
                   '9' :   PWD + "WIG_nodup\Singh_CRP_EE_IP_Rep1_SRR5099114_nodup.wig",
                   '10' :  PWD + "WIG_nodup\Singh_CRP_EE_IP_Rep2_SRR5099115_nodup.wig",   
                   '11' :  PWD + "WIG_nodup\Singh_CRP_EE_Mock_Rep1_SRR5099116_nodup.wig",
                   '12' :  PWD + "WIG_nodup\Singh_CRP_EE_Mock_Rep2_SRR5099117_nodup.wig",   
                   '13' :  PWD + "WIG_nodup\Singh_CRP_ME_IP_Rep1_SRR5099118_nodup.wig",
                   '14' :  PWD + "WIG_nodup\Singh_CRP_ME_IP_Rep2_SRR5099119_nodup.wig",   
                   '15' :  PWD + "WIG_nodup\Singh_CRP_ME_Mock_Rep1_SRR5099120_nodup.wig",
                   '16' :  PWD + "WIG_nodup\Singh_CRP_ME_Mock_Rep2_SRR5099121_nodup.wig",                     
                    }

#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  "Singh_CRP_EE_IP_Rep1",
           '2' :  "Singh_CRP_EE_IP_Rep2",
           '3' :  "Singh_CRP_EE_Mock_Rep1",
           '4' :  "Singh_CRP_EE_Mock_Rep2", 
           '5' :  "Singh_CRP_ME_IP_Rep1",
           '6' :  "Singh_CRP_ME_IP_Rep2", 
           '7' :  "Singh_CRP_ME_Mock_Rep1",
           '8' :  "Singh_CRP_ME_Mock_Rep2", 
           '9' :   "Singh_CRP_EE_IP_Rep1_nodup",
           '10' :  "Singh_CRP_EE_IP_Rep2_nodup",
           '11' :  "Singh_CRP_EE_Mock_Rep1_nodup",
           '12' :  "Singh_CRP_EE_Mock_Rep2_nodup", 
           '13' :  "Singh_CRP_ME_IP_Rep1_nodup",
           '14' :  "Singh_CRP_ME_IP_Rep2_nodup", 
           '15' :  "Singh_CRP_ME_Mock_Rep1_nodup",
           '16' :  "Singh_CRP_ME_Mock_Rep2_nodup",            
           }

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)


def read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual):
    for sample_name, sample_path in filein_path_dict.items():
        print(f'Now is processing: {sample_path}')
        print(f'Progress: {sample_name}/{len(filein_path_dict)}')
        
        filein=open(filein_path_dict[sample_name], 'r')
        fileout=open(fileout_path_dict[sample_name], 'w')
        
        Ar_of_Cromosome_names=[]
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in Ar_of_Cromosome_names:
                if Auto_or_manual==0:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+str(line[0])+' start=1 step=1\n')
                elif Auto_or_manual==1:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
                Ar_of_Cromosome_names.append(line[0])
            else:
                fileout.write(line[2]+'\n')
            
        filein.close()
        fileout.close()    
    return


read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual)