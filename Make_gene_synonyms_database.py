###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Prepares database of gene synonyms for E. coli, based on w3110 genome annotation,
#GO terms annotation (driven by Ecocyc), Ecocyc database and Uniprot database.
#Returns a table with all the synonyms availible.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn3


#######
#Variables to be defined.
#######

#Path to the working directory.
PWD="F:\E_coli_membrane_proteins"

#Path to GO list of membrane proteins genes.
GO_membrane_input="F:\E_coli_membrane_proteins\Databases\GO_Ecocyc_CC_membrane_terms.txt" #+

#Path to the ecocyc GA_all database.
Ecocyc_GO_input="F:\E_coli_membrane_proteins\Databases\GO_Ecocyc_all_sym.txt" #+

#Path to Ecocyc database, contains info for all organism genes.
Ecocyc_input="F:\E_coli_membrane_proteins\Databases\Ecocyc_genes.col" #+

#Path to Uniprot database, contains info for all organism genes.
Uniprot_input="F:\E_coli_membrane_proteins\Databases\\Uniprot_ecoli_genes.txt" #+

#Path to E. coli W3110 genome annotation.
W3110_annotation="C:\Sutor\science\TopoI_Topo-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\Escherichia_coli_K_12_w3110_Mu_genes_annotation.gff" #+


#######
#Checks if directory exists and if not creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return

#Output path.
Out_path=PWD
Dir_check_create(Out_path)
Dir_check_create(PWD+'\Tables_of_synonyms\\')



#######
#Reads file with GO-listed membrane proteins genes. Returns unique genes.
#######

def read_GO_input(file_input):
    #Read data, make dict of ars.
    filein=open(file_input)
    Dict_of_GO_groups={}
    Dict_of_GO_groups_names=[]
    for line in filein:
        if len(Dict_of_GO_groups)==0 and line[0]=="#":
            ar=[]
            group_name=line[3:-4]
            Dict_of_GO_groups_names.append(group_name)
            Dict_of_GO_groups[group_name]=[]
        elif len(Dict_of_GO_groups)!=0 and line[0]=="#":
            Dict_of_GO_groups[group_name]=ar
            group_name=line[3:-4]
            Dict_of_GO_groups_names.append(group_name)
            Dict_of_GO_groups[group_name]=[]
            ar=[]
        elif line[0]!="#":
            line=line.rstrip(',\n').split('\t')[1].split(',')
            for gene in line:
                if gene not in ar:
                    ar.append(gene)
    Dict_of_GO_groups[group_name]=ar
    filein.close()
    
    #Plot Venn diagram shows overlapping between genes sets.
    if len(Dict_of_GO_groups_names)==2:
        venn2([set(Dict_of_GO_groups[Dict_of_GO_groups_names[0]]), set(Dict_of_GO_groups[Dict_of_GO_groups_names[1]])], 
              (Dict_of_GO_groups_names[0], Dict_of_GO_groups_names[1]))
        plt.show()    
    
    if len(Dict_of_GO_groups_names)==3:
        venn3([set(Dict_of_GO_groups[Dict_of_GO_groups_names[0]]), set(Dict_of_GO_groups[Dict_of_GO_groups_names[1]]), set(Dict_of_GO_groups[Dict_of_GO_groups_names[2]])], 
          (Dict_of_GO_groups_names[0], Dict_of_GO_groups_names[1], Dict_of_GO_groups_names[2]))
        plt.show()
    
    #Reorder data - make dict of genes: {'Gene_name' : 'Gene_set1; Gene_set2;'}
    Dict_of_genes={}
    for set_name, set_list in Dict_of_GO_groups.items():
        for gene_name in set_list:
            if gene_name not in Dict_of_genes:
                Dict_of_genes[gene_name]=set_name
            else:
                old_set_name=Dict_of_genes[gene_name]
                new_set_name=old_set_name+'; '+set_name
                Dict_of_genes[gene_name]=new_set_name
    return Dict_of_GO_groups, Dict_of_genes

Membrane_GO_Genes_groups_dict, Membrane_GO_Genes_dict=read_GO_input(GO_membrane_input) #Dict of ars and dict of genes names.
print(len(Membrane_GO_Genes_dict), 'GO Ecocyc membrane proteins')
#print('\n')

All_GO_Genes_groups_dict, All_GO_Genes_dict=read_GO_input(Ecocyc_GO_input) #Dict of ars and dict of genes names.
print(len(All_GO_Genes_dict), 'GO Ecocyc all genes')
#print('\n')


#######
#Parsing gff file and preparing gene annotation.
#######

def parse_gff_annotation(gff_inpath):
    #Read gff file.
    filein=open(gff_inpath, 'r')
    #Classify some genes.
    genes_annotation={'Gene': {},
                      'rRNA': {},
                      'tRNA': {},
                      'ncRNA': {}
                      }
    #All genes together.
    genes_annotation_tg={}
    #List all groups of genes and count them.
    data_source={}
    
    for line in filein:
        line=line.rstrip()
        if line[0]!='#':
            line=line.split('\t')
            #What occurs in the annotation:
            if line[1] not in data_source:
                data_source[line[1]]={line[2]: 1}
            else:
                if line[2] not in data_source[line[1]]:
                    data_source[line[1]][line[2]]=1
                else:
                    data_source[line[1]][line[2]]+=1
            #Classify genes:
            #Protein coding genes.
            if line[1]=='ena' and line[2]=='gene':
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='protein_coding':
                        gene_name=annot[1].split('=')[1]
                        gene_id=annot[4].split('=')[1]
                        genes_annotation['Gene'][gene_name]=gene_id
                        genes_annotation_tg[gene_name]=[gene_id]
            #rRNA genes.
            elif line[1]=='ena' and line[2]=='rRNA_gene':
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='rRNA':
                        gene_name=annot[1].split('=')[1]
                        gene_id=annot[3].split('=')[1]
                        genes_annotation['rRNA'][gene_name]=gene_id
                        genes_annotation_tg[gene_name]=[gene_id]
            #tRNA genes.
            elif line[1]=='ena' and line[2]=='tRNA_gene':
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='tRNA':
                        gene_name=annot[1].split('=')[1]
                        gene_id=annot[3].split('=')[1]
                        genes_annotation['tRNA'][gene_name]=gene_id
                        genes_annotation_tg[gene_name]=[gene_id]
            #Other non-coding RNAs.
            elif line[1]=='Rfam' and line[2]=='ncRNA_gene':
                annot=line[8].split(';')
                if annot[1][0]=="N":
                    if annot[2].split('=')[1]=='ncRNA':
                        gene_name=annot[1].split('=')[1]
                        gene_id=annot[3].split('=')[1]
                        genes_annotation['ncRNA'][gene_name]=gene_id
                        genes_annotation_tg[gene_name]=[gene_id]
    filein.close()            
    return genes_annotation_tg, genes_annotation, data_source

W3110_GFF_genome_annotation=parse_gff_annotation(W3110_annotation) #Dict.
print(len(W3110_GFF_genome_annotation[0]), 'W3110 GFF annotation')
#print(W3110_GFF_genome_annotation[0])
#print('\n')


#######
#Parsing Uniprot database for E. coli genes.
#######

def parse_Uniprot(input_path):
    up_db=open(input_path, 'r')
    dict_of_genes={}
    for line in up_db:
        if line[0] not in ['#']:
            line=line.rstrip().split('\t')
            names_ar=[]
            #Ordered locus name.
            locus_on=line[0].split('; ')
            if len(locus_on)==2:
                locus_on_b=locus_on[0]
                locus_on_jw=locus_on[1]
                names_ar+=[locus_on_b, locus_on_jw]
            elif len(locus_on)==1:
                locus_on_b_jw=locus_on[0]
                names_ar+=[locus_on_b_jw]
                #print(locus_on_b_jw)
            #Swiss-Prot Entry name.
            SP_entry_name=line[1]
            if SP_entry_name!='-':
                names_ar+=[SP_entry_name]
            #AC accession number.
            Gene_AC=line[2]
            if Gene_AC!='-':
                names_ar+=[Gene_AC]            
            #EcoGene.
            EcoGene_ID=line[3]
            if EcoGene_ID!='-':
                names_ar+=[EcoGene_ID]               
            #Gene name and synonyms.
            Gene_names=line[5].split('; ')
            if len(Gene_names)>1:
                Primary_name=Gene_names[0]
                Synonyms=Gene_names[1:]
                names_ar+=Gene_names
                #print(Primary_name, Synonyms)
            elif len(Gene_names)==1:
                Primary_name=Gene_names[0]
                names_ar+=Gene_names
                #print(Primary_name)
            #Combine into dict.
            if Primary_name!='-':
                dict_of_genes[Primary_name]=names_ar
    return dict_of_genes


Uniprot_data=parse_Uniprot(Uniprot_input)
print(len(Uniprot_data), 'Uniprot data')
#print(Uniprot_data)
#print('\n')



#######
#Parsing Ecocyc database for E. coli genes.
#######

def parse_Ecocyc(input_path):
    ec_db=open(input_path, 'r')
    dict_of_genes={}
    for line in ec_db: 
        if line[0] not in ['#']:
            line=line.rstrip().split('\t')
            names_ar=[]
            #Unique ID.
            Unique_ID=line[0]
            names_ar.append(Unique_ID)
            #Ordered locus name (BLATTNER-ID).
            locus_on_b=line[1]
            if locus_on_b!='':
                names_ar.append(locus_on_b)
            #Primary name.
            Primary_name=line[2]
            names_ar.append(Primary_name)
            #Swiss-Prot Entry name.
            SP_entry_name=line[4] 
            if SP_entry_name!='':
                names_ar.append(SP_entry_name)            
            #Synonyms.
            Synonyms=line[8:12] 
            Synonyms_no_empty=list(filter(None, Synonyms)) #With the help of https://stackoverflow.com/questions/3845423/remove-empty-strings-from-a-list-of-strings
            names_ar+=Synonyms_no_empty
            dict_of_genes[Primary_name]=names_ar
    return dict_of_genes

Ecocyc_data=parse_Ecocyc(Ecocyc_input)
print(len(Ecocyc_data), 'Ecocyc data')
#print(Ecocyc_data)
#print('\n')


#######
#Merge two lists.
#######

def merge_lists(list1, list2):
    merged_list=list1
    for item in list2:
        if item not in list1:
            merged_list.append(item)
    return merged_list

#######
#Write table with synonyms.
#######

def Write_dict_of_lists(dict_of_lists, pathout):
    fileout=open(pathout, 'w')
    for primary_name, data in dict_of_lists.items():
        fileout.write(f'{primary_name}\t')
        for i in range(len(data)-1):
            fileout.write(f'{data[i]}\t')
        fileout.write(f'{data[len(data)-1]}\n')
    return

#######
#Combine datasets with the W3110 GFF is a base.
#######

def combine_datasets(w3110_gff, ecocyc, uniprot, GO_ecocyc, GO_membrane, pathout):
    GFF_genes_listed=len(w3110_gff)
    Ecocyc_genes_listed=len(ecocyc)  
    Uniprot_genes_listed=len(uniprot) 
    
    #Merge Ecocyc and Uniprot: add Uniprot data to queries from Ecocyc.
    Ecocyc_Uniprot_merged={}
    count_eco_uni=0
    count_eco_uni_syn=0
    count_eco_uni_locus_tag=0
    for ecocyc_name, ecocyc_data in ecocyc.items():
        if ecocyc_name in uniprot:      #Look for primary gene name.
            eco_uni_combined=merge_lists(ecocyc_data, uniprot[ecocyc_name])
            count_eco_uni+=1
            Ecocyc_Uniprot_merged[ecocyc_name]=eco_uni_combined
            continue
        else:
            found_locus_tag=0
            for uniprot_name, uniprot_data in uniprot.items(): #If nothing found with primary name, search for locus tag.
                if ecocyc_data[1] in uniprot_data:
                    count_eco_uni_locus_tag+=1  
                    found_locus_tag=1
                    eco_uni_combined=merge_lists(ecocyc_data, uniprot_data)
                    Ecocyc_Uniprot_merged[ecocyc_name]=eco_uni_combined
                    continue 
            if found_locus_tag==0:
                found_syn=0
                for uniprot_name, uniprot_data in uniprot.items(): #If nothing returned with locus tag, look among synonyms.
                    if ecocyc_name in uniprot_data:
                        count_eco_uni_syn+=1      
                        found_syn=1
                        eco_uni_combined=merge_lists(ecocyc_data, uniprot_data)
                        Ecocyc_Uniprot_merged[ecocyc_name]=eco_uni_combined
                        continue          
                if found_syn==0:        #If nothing is still be found - just take quiery as it is.
                    eco_uni_combined=ecocyc_data
                    Ecocyc_Uniprot_merged[ecocyc_name]=eco_uni_combined
                    continue                    
    print('Ecocyc vs Uniprot by Primary name, among synonyms, and by locus tag', Ecocyc_genes_listed, Uniprot_genes_listed, 
          count_eco_uni, count_eco_uni_syn, count_eco_uni_locus_tag, count_eco_uni+count_eco_uni_syn+count_eco_uni_locus_tag, len(Ecocyc_Uniprot_merged))    
    
    #Merge Uniprot and Ecocyc: add unique subjects from Uniprot, didn't identify on a previous step of merging.
    count_uni_eco=0
    count_uni_eco_syn=0
    count_uni_eco_locus_tag=0
    for uniprot_name, uniprot_data in uniprot.items():
        if uniprot_name in ecocyc: #Look for primary gene name.
            count_uni_eco+=1
            continue
        else:
            found_locus_tag=0
            for ecocyc_name, ecocyc_data in ecocyc.items(): #If nothing found with primary name, search for locus tag.
                if uniprot_data[0] in ecocyc_data:
                    count_uni_eco_locus_tag+=1
                    continue
            if found_locus_tag==0:
                found_syn=0
                for ecocyc_name, ecocyc_data in ecocyc.items(): #If nothing returned with locus tag, look among synonyms.
                    if uniprot_name in ecocyc_data:
                        count_uni_eco_syn+=1     
                        found_syn=1
                        continue
                if found_syn==0: #If nothing is still be found - just take quiery as it is.
                    Ecocyc_Uniprot_merged[uniprot_name]=uniprot_data
                    continue
    print('Uniprot vs Ecocyc by Primary name, among synonyms, and by locus tag', Uniprot_genes_listed, Ecocyc_genes_listed, 
          count_uni_eco, count_uni_eco_syn, count_uni_eco_locus_tag, count_uni_eco+count_uni_eco_syn+count_uni_eco_locus_tag, len(Ecocyc_Uniprot_merged)) 
    
    Ecocyc_Uniprot_merged_genes_listed=len(Ecocyc_Uniprot_merged) 
    
    #Overlap W3110 GFF and Ecocyc merged with Uniprot.
    count_gff_eco=0
    count_gff_eco_syn=0
    W3110_based_list_of_syns={}
    for gff_name, gff_data in w3110_gff.items():
        if gff_name in Ecocyc_Uniprot_merged: #Look for primary gene name.
            count_gff_eco+=1
            gff_eco_uni_combined=merge_lists(gff_data, Ecocyc_Uniprot_merged[gff_name])
            W3110_based_list_of_syns[gff_name]=gff_eco_uni_combined
            continue
        else:
            found=0
            for ecocyc_name, ecocyc_data in Ecocyc_Uniprot_merged.items():
                if gff_name in ecocyc_data:     #If nothing returned with locus tag, look among synonyms.
                    count_gff_eco_syn+=1
                    gff_eco_uni_combined=merge_lists(gff_data, ecocyc_data)
                    W3110_based_list_of_syns[gff_name]=gff_eco_uni_combined
                    found=1
                    continue
            if found==0:
                W3110_based_list_of_syns[gff_name]=gff_data     #If nothing is still be found - just take quiery as it is.
                continue
    print('W3110 GFF vs Uniprot/Ecocyc by Primary name and among synonyms', GFF_genes_listed, Ecocyc_Uniprot_merged_genes_listed, 
          count_gff_eco, count_gff_eco_syn, count_gff_eco+count_gff_eco_syn, len(W3110_based_list_of_syns))
    
    #Overlap W3110 GFF and GO Membrane database.
    GO_membrane_genes_listed=len(GO_membrane)
    count_GO_all_ecouni=0
    count_GO_all_ecouni_syn=0
    W3110_based_list_of_syns_and_membrane_info={}
    for GO_gene, gene_type in GO_membrane.items():
        if GO_gene in W3110_based_list_of_syns: #Look for primary gene name.
            count_GO_all_ecouni+=1
            W3110_based_list_of_syns_and_membrane_info[GO_gene]=[gene_type]+W3110_based_list_of_syns[GO_gene]
            continue
        else:
            found=0
            for w3110_name, w3110_syns in W3110_based_list_of_syns.items():
                if GO_gene in w3110_syns:      #If nothing returned with locus tag, look among synonyms.
                    count_GO_all_ecouni_syn+=1
                    found=1
                    W3110_based_list_of_syns_and_membrane_info[w3110_name]=[gene_type]+W3110_based_list_of_syns[w3110_name]
                    continue
            if found==0:
                print(GO_gene, gene_type)
                continue
    #Add not membrane w3110 genes.
    for w3110_name, gene_syns in W3110_based_list_of_syns.items():
        if w3110_name not in W3110_based_list_of_syns_and_membrane_info:
            W3110_based_list_of_syns_and_membrane_info[w3110_name]=['-']+W3110_based_list_of_syns[w3110_name]
    print('GO membrane vs W3110 GFF powered with Uniprot/Ecocyc synonyms by Primary name and among synonyms', GO_membrane_genes_listed, GFF_genes_listed, 
          count_GO_all_ecouni, count_GO_all_ecouni_syn, count_GO_all_ecouni+count_GO_all_ecouni_syn)  
    
    #Write combined data.
    Write_dict_of_lists(W3110_based_list_of_syns, f'{pathout}\Tables_of_synonyms\\E_coli_W3110_based_table_of_synonyms.txt')
    Write_dict_of_lists(Ecocyc_Uniprot_merged, f'{pathout}\Tables_of_synonyms\\E_coli_Uniprot_Ecocyc_merged_table_of_synonyms.txt')  
    Write_dict_of_lists(W3110_based_list_of_syns_and_membrane_info, f'{pathout}\Tables_of_synonyms\\E_coli_W3110_based_table_of_synonyms_membrane_info.txt')  

    return

combine_datasets(W3110_GFF_genome_annotation[0], Ecocyc_data, Uniprot_data, All_GO_Genes_dict, Membrane_GO_Genes_dict, Out_path)