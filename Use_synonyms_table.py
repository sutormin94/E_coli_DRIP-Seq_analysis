###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Reads table of E. coli W3110 genes synonyms, identifies synonyms o request.
#Additionally, returns information about membrane localisation of query gene.
###############################################

#######
#Packages to be imported.
#######

import manage_synonyms
from manage_synonyms import Synonyms as syn
import pandas as pd

#######
#Variables to be defined.
#######

#Path to the synonyms table.
Path_to_syns_table="F:\E_coli_membrane_proteins\Tables_of_synonyms\E_coli_W3110_based_table_of_synonyms_membrane_info.txt"

#List (table) of genes to be ajusted with synonyms table.
Path_to_genes_list="F:\Signal_over_TUs\Signal_of_TUs_tab_all\All_genes\Signal_over_TUs_All_genes.txt"

#Path to RegulonDB table of TFs sites.
TF_path="F:\RegulonDB_E_coli\BindingSiteSet.txt"


Dataframe_of_syns=pd.read_csv(Path_to_syns_table, sep='\t', index_col=0, header=None).T #Genes become columns, not rows.
Some_genes_list=pd.read_csv(Path_to_genes_list, sep='\t', header=0)


#######
#Add membrane assignment to some table.
#######

def add_membrane_info(some_genes_dataframe, dataframe_of_syns, path_out):
    Addit_info_dict={'Gene_name': [], 'Membrane_localization': [], 'Synonyms' : []}
    for index, row in some_genes_dataframe.iterrows():
        gene_name=row['Gene_name']
        mem_and_syns=syn.return_info(gene_name, dataframe_of_syns)
        if mem_and_syns!=-1:
            Addit_info_dict['Gene_name'].append(gene_name)
            Addit_info_dict['Membrane_localization'].append(mem_and_syns[1])
            Addit_info_dict['Synonyms'].append(mem_and_syns[2])
        else:
            Addit_info_dict['Gene_name'].append('-')
            Addit_info_dict['Membrane_localization'].append('-')
            Addit_info_dict['Synonyms'].append('-')            
    Addit_info_df=pd.DataFrame(Addit_info_dict, columns=['Gene_name', 'Membrane_localization', 'Synonyms'])
    
    List_and_add_info=pd.merge(some_genes_dataframe, Addit_info_df, how='left', on='Gene_name')
    List_and_add_info.to_csv(path_out, sep='\t', index=False)   
    return


#######
#Concatenate list.
#######

def concat_list(list_in, sep):
    str_out=''
    for element in list_in:
        str_out+=str(element)+sep
    return str_out

#######
#Read RegulonDB table for TF sites, prepare dataframe to be merged.
#######

def read_regulonDB_TF(path_in, dataframe_of_syns, sep, some_genes_dataframe, path_out):
    filein=open(path_in, 'r')
    Genes_TF_dict={}
    for line in filein:
        if line[0] not in ['#']:
            line=line.rstrip().split('\t')
            TF=line[1]
            TUs=line[7].split('-')
            TU_proximal=TUs[0]
            if len(TU_proximal)>4: #Not a single gen, but an operon.
                gen_name_basis=TU_proximal[:4]
            else:
                gen_name_basis=TU_proximal
            if len(line)==14:
                Evidence=line[13]
            else:
                Evidence='-'
            if gen_name_basis not in Genes_TF_dict:
                Genes_TF_dict[gen_name_basis]={'Number_of_TF_sites' : 1, 'TF_list' : [TF],  'Evidence_list' : [Evidence]}
            else:
                Genes_TF_dict[gen_name_basis]['Number_of_TF_sites']+=1
                Genes_TF_dict[gen_name_basis]['TF_list'].append(TF)
                Genes_TF_dict[gen_name_basis]['Evidence_list'].append(Evidence)
                    
    #Reshape TF data.
    Missed=0
    Genes_TF_data={'Gene_name_RDB' : [], 'Number_of_TF_sites': [], 'TF_list' : [], 'Evidence_list' : [], 'Gene_name' : []}
    for gene_name, gene_data in Genes_TF_dict.items():
        w3110_gene_name=syn.return_info(gene_name, dataframe_of_syns)
        if w3110_gene_name!=-1:
            Genes_TF_data['Gene_name_RDB'].append(gene_name)
            Genes_TF_data['Number_of_TF_sites'].append(gene_data['Number_of_TF_sites'])
            Genes_TF_data['TF_list'].append(concat_list(gene_data['TF_list'], sep))
            Genes_TF_data['Evidence_list'].append(concat_list(gene_data['Evidence_list'], sep))
            Genes_TF_data['Gene_name'].append(w3110_gene_name[0])
        else:
            print(f'{gene_name} was not found among E. coli W3110 genes, hence is dropped.')
            Missed+=1
    
    Genes_TF_dataframe=pd.DataFrame(Genes_TF_data, columns=['Gene_name_RDB', 'Number_of_TF_sites', 'TF_list', 'Evidence_list', 'Gene_name'])   
    print(Genes_TF_dataframe)
    print(f'{Missed}/{len(Genes_TF_dict)} genes were missed (do not correspond to any gene from E. coli W3110 annotation)')
    List_and_added_info=pd.merge(some_genes_dataframe, Genes_TF_dataframe, how='left', on='Gene_name')
    List_and_added_info.to_csv(path_out, sep='\t', index=False)       
    return List_and_added_info

Datframe_with_TF=read_regulonDB_TF(TF_path, Dataframe_of_syns, ';', Some_genes_list, 'F:\Signal_over_TUs\Signal_of_TUs_tab_all\All_genes\Signal_over_TUs_and_TF.txt')
add_membrane_info(Datframe_with_TF, Dataframe_of_syns, 'F:\Signal_over_TUs\Signal_of_TUs_tab_all\All_genes\Signal_over_TUs_and_TF_syns_Mem.txt')