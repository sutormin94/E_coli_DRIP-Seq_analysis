###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Reads table of E. coli W3110 genes synonyms, identifies synonyms o request.
#Additionally, returns information about membrane localisation of query gene.
###############################################

#######
#Packages to be imported.
#######

import pandas as pd

class Synonyms(object):
    
    #######
    #Returns information for query: primary gene name (based on E. coli W3110 annotation), list of synonyms and membrane assignment.
    #######    

    def return_info(query, syns_dataframe):
        if query in syns_dataframe:
            Primary_name=query
            Membrane_assignment=syns_dataframe[query].dropna()[1]
            Synonyms=syns_dataframe[query].dropna()[1:].tolist()
            #print(f'{query} is found: {Primary_name}, {Synonyms}, {Membrane_assignment}')
            return Primary_name, Membrane_assignment, Synonyms   
        else:
            found=0
            for gene in syns_dataframe:
                Primary_name=gene
                Membrane_assignment=syns_dataframe[gene].dropna()[1]
                Synonyms=syns_dataframe[gene].dropna()[1:].tolist()
                if query in Synonyms:
                    #print(f'{query} is found: {Primary_name}, {Synonyms}, {Membrane_assignment}')
                    found=1
                    return Primary_name, Membrane_assignment, Synonyms                
            if found==0:
                print(f'{query} was not found in synonyms database!')       
                return -1

 