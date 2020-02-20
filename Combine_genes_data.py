###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes tab files contain cumulative average signal of differet factors
#over TUs (produced by FE_over_US_GB_DS.py script).
#Combine information about TUs bodies signal with expression data in a dataframe and write it.
#Additionally computes and plots correlation matrix for signals over TUs bodies.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import scipy
import pandas as pd
from pandas import DataFrame
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from matplotlib import cm as cm

#Path to the working directory.
PWD='F:\Signal_over_TUs'
#Dictionary of input TAB data files (information about the signal over TUs).
TUs_data_dict={'TopoA -Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\TopoA -Rif_over_All genes_15000bp.txt',
               'TopoA +Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\TopoA +Rif_over_All genes_15000bp.txt',
               'TopoIV' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\TopoIV Cfx_over_All genes_15000bp.txt',
               'Gyrase -Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\Gyrase Cfx_over_All genes_15000bp.txt',
               'Gyrase +Rif' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\Gyrase Cfx +Rif_over_All_genes_15000bp.txt',
               'MukB' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\MukB_over_All genes_15000bp.txt',
               'RpoB' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\RpoB_over_All genes_15000bp.txt',
               'H-NS' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\HNS_over_All genes_15000bp.txt',
               'Fis' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\Fis_over_All genes_15000bp.txt',
               'GC' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\GC_over_All genes_15000bp.txt',
               'Expression' : 'C:\Sutor\science\TopoI_Topo-Seq\Scripts\TopoA_ChIP-Seq\Additional_genome_features\DOOR_Mu_del_cor_genes_expression.txt',
               'MatP' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\MatP_over_All genes_15000bp.txt',
               'PolSofi' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\PolSofi_over_All genes_15000bp.txt',
               'RpoS' : 'F:\Signal_over_TUs\Signal_of_TUs_tab\All_genes\RpoS_over_All genes_15000bp.txt',
               }
#Main table (to be added to).
Main_table='F:\Signal_over_TUs\Signal_of_TUs_tab_all\All_genes\Old\\Signal_over_TUs_and_regulonDB_info_eq_len_All_genes_membrane_syns_TF.xlsx'
#Set of genes.
Set_of_genes='All_genes'

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
Dir_check_create(PWD+'\Signal_of_TUs_tab_all\\'+Set_of_genes)
Dir_check_create(PWD+'\Figures\Datasets_correlation\\'+Set_of_genes)


#########
##Compute correlation matrix and draw heatmaps.
#########

#Plot diagonal correlation matrix.
def correlation_matrix(df, cor_method, title, outpath):
    fig=plt.figure(figsize=(8,8), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    cax=ax1.imshow(df.corr(method=cor_method), interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1)
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    labels=list(df)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    #Add colorbar, make sure to specify tick locations to match desired ticklabels.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    fig.colorbar(cax, ticks=[-1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00], shrink=0.7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=400, figsize=(8, 8))
    plt.show()
    plt.close()
    return


#######
#Identify clusters in a corralation matrix (hierarchy clustering).
#Code was stolen from https://github.com/TheLoneNut/CorrelationMatrixClustering/blob/master/CorrelationMatrixClustering.ipynb
#with some modifications.
#######

def Clustering(df, cor_method):
    X = df.corr(method=cor_method).values
    d = sch.distance.pdist(X)   # vector of pairwise distances
    L = sch.linkage(d, method='complete')
    ind = sch.fcluster(L, 0.5*d.max(), 'distance')
    columns = [df.columns.tolist()[i] for i in list((np.argsort(ind)))]
    df = df.reindex_axis(columns, axis=1)
    return df


#######
#Read .tab-file with TopoA FE over TUs.
#######

def read_FE_tables_combine_together(data_dict, main_table, cor_method, set_of_genes_name, path_out):
    #Read input dataframes.
    Dict_of_full_dataset={}
    Df_of_critical_data=pd.DataFrame()
    #Initiate with Expression data.
    Expression_dataframe=pd.read_csv(data_dict['Expression'], sep='\t')
    for index, row in Expression_dataframe.iterrows():
        Expression_dot=row['Expression'].replace(',', '.')    
        Expression_dataframe.at[index, 'Expression']=Expression_dot
    Expression_dataframe=Expression_dataframe.astype({'Expression': np.float64})
    Dict_of_full_dataset['Expression']=Expression_dataframe
    Df_of_critical_data=Expression_dataframe
    #Read other dataframes.
    for dataset_name, dataset in data_dict.items():
        data_dataframe=pd.read_csv(dataset, sep='\t')
        if dataset_name!='Expression':
            Dict_of_full_dataset[dataset_name]=data_dataframe
            current_set_df=pd.DataFrame()
            current_set_df['Gene_name']=data_dataframe['Gene_name'] #Extracts column with gene names.
            current_set_df[f'{dataset_name}_FE_GB']=data_dataframe.iloc[: , 5:6] #Extracts column with signal over GB: FE_GB
            Df_of_critical_data=pd.merge(Df_of_critical_data, current_set_df)
            
    #Read main table to merge data to.
    main_data=pd.read_excel(main_table)
    main_and_critical_df=pd.merge(Df_of_critical_data, main_data, how='right', on=['Gene_name', 'GeneID', 'Start', 'End', 'Strand', 'OperonID'])
    
    #Write new dataframe.
    Df_of_critical_data.to_csv(f'{path_out}\Signal_of_TUs_tab_all\\{set_of_genes_name}\\Signal_over_TUs_{set_of_genes_name}.txt', sep='\t', index=False)    
    main_and_critical_df.to_csv(f'{path_out}\Signal_of_TUs_tab_all\\{set_of_genes_name}\\Signal_over_TUs_and_regulonDB_info_eq_len_{set_of_genes_name}.txt', sep='\t', index=False)        
    
    #Data cross-correlations.
    #Prepare only signal-containing columns.
    Data_to_correlate=Df_of_critical_data
    Data_to_correlate.drop(['GeneID', 'Gene_name', 'Start', 'End', 'Strand', 'Gene_description', 'OperonID'], 1, inplace=True)
    #Plot the correlation matrix.
    correlation_matrix(Data_to_correlate, cor_method, f'Correlation of signals over TUs for {set_of_genes_name}', 
                       f'{path_out}\Figures\Datasets_correlation\\{set_of_genes_name}\\Signal_over_TUs_correlation_for_{set_of_genes_name}.png')
    #Perform the hierarchial clustering
    Data_to_correlate_Clusterized=Clustering(Data_to_correlate, cor_method)
    #Plot the correlation matrix after clustering.
    correlation_matrix(Data_to_correlate_Clusterized, cor_method, f'Clusterized correlation of signals over TUs for {set_of_genes_name}', 
                       f'{path_out}\Figures\Datasets_correlation\\{set_of_genes_name}\\Signal_over_TUs_clusterized_correlation_for_{set_of_genes_name}.png')    
    
    return

read_FE_tables_combine_together(TUs_data_dict, Main_table, 'pearson', Set_of_genes, Out_path)