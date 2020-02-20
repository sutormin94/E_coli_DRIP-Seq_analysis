###############################################
##Dmitry Sutormin, 2019##
##TopoI ChIP-Seq analysis##

#The script tests sets of genomic intervals (Peaks, TUs, BIMEs-1, BIMEs-2, IHF sites, Fis sites, H-NS sites, MatP sites, etc.)
#for the enrichment with some continously distributed character (RNApol fold enrichment, score, GC%, etc.) (t-test).
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import binom
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

#######
#Variables to be defined.
#######

#Path to the working directory.
PWD="C:\Sutor\Science\E_coli_ChIP-Seqs\Data_for_ORFs\\"

#Input: Intervals (e.g. EcTopoI peaks) with additional info (FE, GC%, width, etc.).
path_to_intervals_data=PWD + "Genes_EP_expression_and_all_info_more_info.csv"

#Output: path to the dir to store output
Outputpath=PWD + "Interactions_between_factors\\"
if not os.path.exists(Outputpath):
    os.makedirs(Outputpath)
    

#######
#Compare RNApol occupation of BIMEs and overall RNApol occupation.
#######

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels)+1))
    ax.set_xticklabels(labels, size=35)
    ax.set_xlim(0.25, len(labels)+0.75)
    return

def Intervals_occupation(intervals_param_ar, cont_char_masked, name_of_intervals, name_of_cont_char, params, outpath):
    #Plot it.
    pos=[1,2]
    dataset=[intervals_param_ar, cont_char_masked]
    #print(len(pos))
    #print(len(dataset))
    #Violin plots for whole-genome signal vs signal of intervals (e.g. peaks).
    fig=plt.figure(figsize=(7,12), dpi=100)
    plt1=fig.add_subplot(1,1,1) 
    violins=plt1.violinplot(dataset, positions=pos, widths=0.9, showmeans=True, showmedians=True, points=200)
    #print(violins)
    for vio in violins['bodies']:
        vio.set_facecolor('#ff7762')
        vio.set_edgecolor('black')
        vio.set_alpha(1)
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    vbars=violins['cbars']
    vbars.set_linewidth(2)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    labels=[f'{name_of_intervals}\nsites', f'{name_of_cont_char}\nsites']
    set_axis_style(plt1, labels)
    yticknames1=np.arange(params[0], params[1], params[2]) #EcTopoI 0, 35, 1, GC% 0, 100, 10, RNApol 0, 21, 2
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylim(params[3], params[4]) #EcTopoI -0.2, 4, GC% 10, 77 RNApol -1, 21
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=35)   
    plt1.annotate(f'Mean {name_of_intervals}\nsignal={round(np.mean(intervals_param_ar),1)}', xy=(params[5], params[6]), xycoords='data', size=35, rotation=90) #EcTopoI 0.4, 3.8, GC% 0.3, 35, RNApol 0.4, 11
    plt1.annotate(f'Mean {name_of_cont_char}\nsignal={round(np.mean(cont_char_masked),1)}', xy=(params[7], params[6]), xycoords='data', size=35, rotation=90) #EcTopoI 1.4, 3.8, GC% 1.3, 35, RNApol 1.4, 11
    plt1.set_yscale('log')
    
    Intervals_stat=stats.ttest_ind(cont_char_masked, intervals_param_ar)
    print(f'\nT-test for all genome sites vs {name_of_intervals} sites {name_of_cont_char} FE means\n' + 'p-value=' + str(Intervals_stat[1]) +'\n' + 't-statistic=' + str(Intervals_stat[0]) + '\n')    
    #plt.show()
    plt.tight_layout()
    plt.savefig(f'{outpath}{name_of_cont_char}_signal_of_{name_of_intervals}_sites_and_genome.png', dpi=400, figsize=(7, 12)) 
    plt.close()    
    return

#######
#2D plot.
#######

def plot_2D(intervals_char_1_orig, intervals_char_2_orig, char_name_1, char_name_2, intervals_name, method, limit_1, limit_2, outpath):
    ##Filtering 1.
    intervals_char_1_int=intervals_char_1_orig[intervals_char_1_orig<limit_1]
    intervals_char_2_int=intervals_char_2_orig[intervals_char_1_orig<limit_1]
    ##Filtering 2.
    intervals_char_1=intervals_char_1_int[intervals_char_2_int<limit_2]
    intervals_char_2=intervals_char_2_int[intervals_char_2_int<limit_2]    
    
    ##Fitting.
    if method=='lin':
        #Linear fitting of linear data.
        fit=np.polyfit(intervals_char_1, intervals_char_2, 1)
        print(fit)
        fit_fn=np.poly1d(fit)  
        ##Pearson correlation.
        pearson_cor=scipy.stats.pearsonr(intervals_char_1, intervals_char_2)
        print(f'Paerson correlation ({char_name_1}, {char_name_2}) for {intervals_name} {pearson_cor}')  
        ##Plot data.
        fig=plt.figure(figsize=(5,4), dpi=100)
        ax=fig.add_subplot(111)
        ax.plot(intervals_char_1, intervals_char_2, 'ro')
        ax.plot(intervals_char_1, fit_fn(intervals_char_1), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3))) 
        ax.annotate(f'Pearson correlation=\n{round(pearson_cor[0], 3)}', xy=(0.6, 0.2), xycoords='axes fraction', size=9)
        ax.set_xlabel(char_name_1, fontsize=12)
        ax.set_ylabel(char_name_2, fontsize=12)
        ax.set_title(f'{char_name_1} vs {char_name_2} for\n {intervals_name} peaks', fontsize=12)
        plt.tight_layout()        
        char_name_1=char_name_1.replace("/", '_')
        plt.savefig(f'{outpath}\{char_name_1}_lin_vs_{char_name_2}_lin_for_{intervals_name}_peaks.png', dpi=300, figsize=(5,4))   
        
    if method=='semilogx':
        #Additional data filtration. Getting rid of zero values.
        intervals_char_1_fin=intervals_char_1[intervals_char_1>0]
        intervals_char_2_fin=intervals_char_2[intervals_char_1>0]
        #Linear fitting of linear data.
        fit=np.polyfit(np.log10(intervals_char_1_fin), intervals_char_2_fin, 1)
        print(fit)
        fit_fn=np.poly1d(fit)  
        ##Pearson correlation.
        pearson_cor=scipy.stats.pearsonr(np.log10(intervals_char_1_fin), intervals_char_2_fin)
        print(f'Paerson correlation ({char_name_1} log, {char_name_2}) for {intervals_name} {pearson_cor}')  
        ##Plot data.
        fig=plt.figure(figsize=(5,4), dpi=100)
        ax=fig.add_subplot(111)
        ax.plot(np.log10(intervals_char_1_fin), intervals_char_2_fin, 'go')
        ax.plot(np.log10(intervals_char_1_fin), fit_fn(np.log10(intervals_char_1_fin)), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3))) 
        ax.annotate(f'Pearson correlation=\n{round(pearson_cor[0], 3)}', xy=(0.6, 0.2), xycoords='axes fraction', size=9)
        ax.set_xlabel(f'log({char_name_1})', fontsize=12)
        ax.set_ylabel(char_name_2, fontsize=12)
        ax.set_title(f'log({char_name_1}) vs {char_name_2} for\n {intervals_name} peaks', fontsize=12)
        plt.tight_layout()        
        char_name_1=char_name_1.replace("/", '_')
        plt.savefig(f'{outpath}\{char_name_1}_log_vs_{char_name_2}_lin_for_{intervals_name}_peaks.png', dpi=300, figsize=(5,4))  
        
    elif method=='log':
        #Additional data filtration. Getting rid of zero values.
        intervals_char_1_fin=intervals_char_1[(intervals_char_1>0) & (intervals_char_2>0)]
        intervals_char_2_fin=intervals_char_2[(intervals_char_1>0) & (intervals_char_2>0)]      
        #Linnear fitting of log data.
        fit=np.polyfit(np.log10(intervals_char_1_fin), np.log10(intervals_char_2_fin), 1)
        print(fit)
        fit_fn=np.poly1d(fit) 
        ##Pearson correlation.
        pearson_cor=scipy.stats.pearsonr(np.log10(intervals_char_1_fin), np.log10(intervals_char_2_fin))
        print(f'Paerson correlation ({char_name_1} log, {char_name_2}) log for {intervals_name} {pearson_cor}')       
        ##Plot data.
        fig=plt.figure(figsize=(5,4), dpi=100)
        ax=fig.add_subplot(111)       
        ax.plot(np.log10(intervals_char_1_fin), np.log10(intervals_char_2_fin), 'bo')
        ax.plot(np.log10(intervals_char_1_fin), fit_fn(np.log10(intervals_char_1_fin)), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))   
        ax.annotate(f'Pearson correlation=\n{round(pearson_cor[0], 3)}', xy=(0.6, 0.2), xycoords='axes fraction', size=9)
        ax.set_xlabel(f'log({char_name_1})', fontsize=12)
        ax.set_ylabel(f'log({char_name_2})', fontsize=12)
        ax.set_title(f'log({char_name_1}) vs log({char_name_2}) for\n {intervals_name} peaks', fontsize=12)
        plt.tight_layout()        
        char_name_1=char_name_1.replace("/", '_')
        plt.savefig(f'{outpath}\{char_name_1}_log_vs_{char_name_2}_log_for_{intervals_name}_peaks.png', dpi=300, figsize=(5,4))        
    plt.close()
    return

#######
#Linear regression.
#######

def run_regression(Intervals_data_orig, limit):
    ##Filtering.
    Intervals_data=Intervals_data_orig[Intervals_data_orig['EcTopoI CTD-/Rif-']<limit]
    
    #Based on https://realpython.com/linear-regression-in-python/
    model = LinearRegression(fit_intercept=True, normalize=True, n_jobs=-1)
    
    features_list=['CsiR Aquino', 'Nac Aquino', 'NtrC Aquino', 'OmpR Aquino', 'Fur Beauchene', 'BolA Dressaire', 'NsrR Mehta', 'FNR Myers', 'Lrp Kroner',                    
                   'EcTopoI score', 'CRP Singh', 'RpoS Peano', 'HNS Kahramanoglou', 
                   'Cra Kim', 'ArcA Park', 'GadE Seo', 'GadW Seo', 'OxyR Seq', 'RpoS Seo', 'SoxR Seo', 'SoxS Seo', 'GadX Seo', 'RpoD Myers', 'Dps Antipov',         
                   'RpoB Borukhov', 'RpoB Kahramanoglou', 'Fis Kahramanoglou', 'MatP Nolivos', 'MukB Nolivos', 'RNA-Seq Sutormin', 'TopoIV Sayyed',       
                   'TopoIV Sutormin', 'Gyrase Sutormin', 'Gyrase Rif Sutormin', 'GC']
    model.fit(Intervals_data[features_list], Intervals_data['EcTopoI CTD-/Rif-'])
    r_sq = model.score(Intervals_data[features_list], Intervals_data['EcTopoI CTD-/Rif-'])          
    
    
    #model.fit(Intervals_data[['GC', 'RpoB Borukhov', 'RpoD Myers', 'RNA-Seq Sutormin', 'Gyrase Sutormin', 'EcTopoI score']], np.log(Intervals_data['EcTopoI CTD-/Rif-']))
    #r_sq = model.score(Intervals_data[['GC', 'RpoB Borukhov', 'RpoD Myers', 'RNA-Seq Sutormin', 'Gyrase Sutormin', 'EcTopoI score']], np.log(Intervals_data['EcTopoI CTD-/Rif-']))      
    #model.fit(Intervals_data[['GC', 'RpoB Borukhov', 'RpoD Myers', 'Fis Kahramanoglou', 'HNS Kahramanoglou', 
    #                          'RNA-Seq Sutormin', 'Gyrase Sutormin', 'EcTopoI score']], np.log(Intervals_data['EcTopoI CTD-/Rif-']))
    #r_sq = model.score(Intervals_data[['GC', 'RpoB Borukhov', 'RpoD Myers', 'Fis Kahramanoglou', 'HNS Kahramanoglou', 
    #                                   'RNA-Seq Sutormin', 'Gyrase Sutormin', 'EcTopoI score']], np.log(Intervals_data['EcTopoI CTD-/Rif-']))    
    #model.fit(Intervals_data[['GC', 'RpoB Borukhov', 'RpoD Myers', 'Fis Kahramanoglou', 'HNS Kahramanoglou', 'MatP Nolivos',
    #                          'MukB Nolivos', 'RNA-Seq Sutormin', 'Gyrase Sutormin', 'EcTopoI score', 'TopoIV Sayyed']], np.log(Intervals_data['EcTopoI CTD-/Rif-']))
    #r_sq = model.score(Intervals_data[['GC', 'RpoB Borukhov', 'RpoD Myers', 'Fis Kahramanoglou', 'HNS Kahramanoglou', 'MatP Nolivos',
    #                                   'MukB Nolivos', 'RNA-Seq Sutormin', 'Gyrase Sutormin', 'EcTopoI score', 'TopoIV Sayyed']], np.log(Intervals_data['EcTopoI CTD-/Rif-']))
    print('Coefficient of determination:', r_sq)
    print('Intercept:', model.intercept_)
    print('Slope:', model.coef_)    
    return


#######
#Principal component analysis.
#######

def PCA_analysis(Intervals_info, outpath, dataset_name):
    #Based on this example: https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60
    features=['CsiR Aquino', 'Nac Aquino', 'NtrC Aquino', 'OmpR Aquino', 'Fur Beauchene', 'BolA Dressaire', 'NsrR Mehta', 'FNR Myers', 'Lrp Kroner',                    
              'CRP Singh', 'RpoS Peano', 'HNS Kahramanoglou', 
              'Cra Kim', 'ArcA Park', 'GadE Seo', 'GadW Seo', 'OxyR Seq', 'RpoS Seo', 'SoxR Seo', 'SoxS Seo', 'GadX Seo', 'RpoD Myers', 'Dps Antipov',         
              'RpoB Borukhov', 'RpoB Kahramanoglou', 'Fis Kahramanoglou', 'MatP Nolivos', 'MukB Nolivos', 'RNA-Seq Sutormin', 'TopoIV Sayyed',       
              'TopoIV Sutormin', 'Gyrase Sutormin', 'Gyrase Rif Sutormin', 'GC', 'DRIP-Seq fw Sutormin', 'DRIP-Seq rv Sutormin']
    #Features.
    x=Intervals_info.loc[:, features].values
    #Standardizing the features.
    x=StandardScaler().fit_transform(x)
    #PCA model.
    pca=PCA(n_components=3)
    #PCA transformation.
    principalComponents=pca.fit_transform(x)   
    #Keep PCA data as dataframe.
    PCA_dataframe=pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2', 'principal component 3']) 
    #Add regions coordinates.
    PCA_dataframe_coord=pd.concat([PCA_dataframe, Intervals_info[['Start']]], axis = 1)
    PCA_dataframe_coord.to_csv(outpath+dataset_name+'_PCA_analysis.csv', sep='\t', header=True, index=True)
    
    #Performance of PCA.
    pca_expl=np.round(pca.explained_variance_ratio_, 3)
    print(pca_expl)    
    
    #Plot PCA results.
    fig = plt.figure(figsize = (12,4))
    ax = fig.add_subplot(1,3,1) 
    ax.set_xlabel(f'Principal Component 1 ({pca_expl[0]})', fontsize = 15)
    ax.set_ylabel(f'Principal Component 2 ({pca_expl[1]})', fontsize = 15)
    ax.set_title('3 component PCA: 1 vs 2', fontsize = 15)
    ax.scatter(PCA_dataframe['principal component 1'], PCA_dataframe['principal component 2'], c='red' , s=50)
    for index_num in PCA_dataframe_coord.index:
        coordinate=PCA_dataframe_coord.loc[index_num, 'Start']
        x_coord=PCA_dataframe_coord.loc[index_num, 'principal component 1']
        y_coord=PCA_dataframe_coord.loc[index_num, 'principal component 2']
        ax.annotate(coordinate, xy=(x_coord, y_coord), xycoords='data', size=4)
    ax.grid()    
    
    ax = fig.add_subplot(1,3,2) 
    ax.set_xlabel(f'Principal Component 1 ({pca_expl[0]})', fontsize = 15)
    ax.set_ylabel(f'Principal Component 3 ({pca_expl[2]})', fontsize = 15)
    ax.set_title('3 component PCA: 1 vs 3', fontsize = 15)
    ax.scatter(PCA_dataframe['principal component 1'], PCA_dataframe['principal component 3'], c='blue' , s=50)
    for index_num in PCA_dataframe_coord.index:
        coordinate=PCA_dataframe_coord.loc[index_num, 'Start']
        x_coord=PCA_dataframe_coord.loc[index_num, 'principal component 1']
        y_coord=PCA_dataframe_coord.loc[index_num, 'principal component 3']
        ax.annotate(coordinate, xy=(x_coord, y_coord), xycoords='data', size=4)    
    ax.grid() 
    
    ax = fig.add_subplot(1,3,3) 
    ax.set_xlabel(f'Principal Component 2 ({pca_expl[1]})', fontsize = 15)
    ax.set_ylabel(f'Principal Component 3 ({pca_expl[2]})', fontsize = 15)
    ax.set_title('3 component PCA: 2 vs 3', fontsize = 15)
    ax.scatter(PCA_dataframe['principal component 2'], PCA_dataframe['principal component 3'], c='green' , s=50)
    for index_num in PCA_dataframe_coord.index:
        coordinate=PCA_dataframe_coord.loc[index_num, 'Start']
        x_coord=PCA_dataframe_coord.loc[index_num, 'principal component 2']
        y_coord=PCA_dataframe_coord.loc[index_num, 'principal component 3']
        ax.annotate(coordinate, xy=(x_coord, y_coord), xycoords='data', size=4)    
    ax.grid()    
    
    plt.tight_layout()
    plt.savefig(outpath+dataset_name+'_PCA_analysis_annotated.png', dpi=300)
    
    
    #Plot dependences of features with PC.
    #plt.figure(figsize = (20,10))
    plt.matshow(pca.components_,cmap='viridis')
    plt.colorbar()
    plt.yticks([0,1,2], ['1st Comp','2nd Comp','3rd Comp'], fontsize=12, rotation=0)
    plt.ylim([-0.5, 2.5]) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    plt.xticks(range(len(features)), features, rotation=65, ha='left')
    plt.tight_layout()
    plt.show()    
    plt.savefig(outpath+dataset_name+'_PCA_analysis_features_impact.png', bbox_inches = "tight", dpi=300, figsize=(20,10))    
    
    return


#######
#Functions wrapper.
#######

def func_wrapper(intervals_sets_path, outpath):
    #Read data
    Intervals_info=pd.read_csv(intervals_sets_path, sep='\t', header=0, index_col=False, dtype={'Start' : np.int64, 'End' : np.int64})
    Intervals_info_plus=Intervals_info[Intervals_info['Strand']=='+']
    Intervals_info_minus=Intervals_info[Intervals_info['Strand']=='-']

    ##2D plots.
    limit_1=3000
    limit_2=1000
    #plot_2D(Intervals_info['Expression_E'], Intervals_info['RpoB Borukhov'], 'Expression_RNA-Seq', 'RpoB_ChIP-Seq_Borukhov', 'ORFs_all', 'lin', limit_1, limit_2, outpath)
    #DRIP with all genes.
    #plot_2D(Intervals_info['Expression_E'], Intervals_info['DRIP-Seq fw Sutormin'], 'Expression_RNA-Seq', 'DRIP-Seq_fw', 'ORFs_all', 'lin', limit_1, limit_2, outpath)
    #plot_2D(Intervals_info['Expression_E'], Intervals_info['DRIP-Seq rv Sutormin'], 'Expression_RNA-Seq', 'DRIP-Seq_rv', 'ORFs_all', 'lin', limit_1, limit_2, outpath)
    #DRIP with + genes.
    #plot_2D(Intervals_info_plus['Expression_E'], Intervals_info_plus['DRIP-Seq fw Sutormin'], 'Expression_RNA-Seq', 'DRIP-Seq_fw', 'ORFs_plus', 'lin', limit_1, limit_2, outpath)
    #plot_2D(Intervals_info_plus['Expression_E'], Intervals_info_plus['DRIP-Seq rv Sutormin'], 'Expression_RNA-Seq', 'DRIP-Seq_rv', 'ORFs_plus', 'lin', limit_1, limit_2, outpath)
    #DRIP with - genes.
    #plot_2D(Intervals_info_minus['Expression_E'], Intervals_info_minus['DRIP-Seq fw Sutormin'], 'Expression_RNA-Seq', 'DRIP-Seq_fw', 'ORFs_minus', 'lin', limit_1, limit_2, outpath)
    #plot_2D(Intervals_info_minus['Expression_E'], Intervals_info_minus['DRIP-Seq rv Sutormin'], 'Expression_RNA-Seq', 'DRIP-Seq_rv', 'ORFs_minus', 'lin', limit_1, limit_2, outpath)
    
    #DRIP with all genes.
    #plot_2D(Intervals_info['RpoB Borukhov'], Intervals_info['DRIP-Seq fw Sutormin'], 'RpoB_ChIP-Seq_Borukhov', 'DRIP-Seq_fw', 'ORFs_all', 'lin', limit_1, limit_2, outpath)
    #plot_2D(Intervals_info['RpoB Borukhov'], Intervals_info['DRIP-Seq rv Sutormin'], 'RpoB_ChIP-Seq_Borukhov', 'DRIP-Seq_rv', 'ORFs_all', 'lin', limit_1, limit_2, outpath)    
    #DRIP with + genes.
    #plot_2D(Intervals_info_plus['RpoB Borukhov'], Intervals_info_plus['DRIP-Seq fw Sutormin'], 'RpoB_ChIP-Seq_Borukhov', 'DRIP-Seq_fw', 'ORFs_plus', 'lin', limit_1, limit_2, outpath)
    #plot_2D(Intervals_info_plus['RpoB Borukhov'], Intervals_info_plus['DRIP-Seq rv Sutormin'], 'RpoB_ChIP-Seq_Borukhov', 'DRIP-Seq_rv', 'ORFs_plus', 'lin', limit_1, limit_2, outpath)    
    #DRIP with - genes.
    #plot_2D(Intervals_info_minus['RpoB Borukhov'], Intervals_info_minus['DRIP-Seq fw Sutormin'], 'RpoB_ChIP-Seq_Borukhov', 'DRIP-Seq_fw', 'ORFs_minus', 'lin', limit_1, limit_2, outpath)
    #plot_2D(Intervals_info_minus['RpoB Borukhov'], Intervals_info_minus['DRIP-Seq rv Sutormin'], 'RpoB_ChIP-Seq_Borukhov', 'DRIP-Seq_rv', 'ORFs_minus', 'lin', limit_1, limit_2, outpath)      

    
    #Run linear regression.
    #run_regression(Intervals_info, limit_1)
    #Run PCA analysis.
    #PCA_analysis(Intervals_info, outpath, 'Operons_factors_PCA')
    
    
    
    ##Violin plots.
    #DRIP-Seq signal on different replichores and strands.
    Ter_position=1593600
    Ori_position=3711900
    Intervals_info_right_fw=Intervals_info.loc[((Intervals_info['Start']>0) & (Intervals_info['Start']<Ter_position)) | (Intervals_info['Start']>Ori_position), 'DRIP-Seq fw Sutormin'].values
    Intervals_info_left_fw=Intervals_info.loc[((Intervals_info['Start']>Ter_position) & (Intervals_info['Start']<Ori_position)), 'DRIP-Seq fw Sutormin'].values
    Intervals_info_right_rv=Intervals_info.loc[((Intervals_info['Start']>0) & (Intervals_info['Start']<Ter_position)) | (Intervals_info['Start']>Ori_position), 'DRIP-Seq rv Sutormin'].values
    Intervals_info_left_rv=Intervals_info.loc[((Intervals_info['Start']>Ter_position) & (Intervals_info['Start']<Ori_position)), 'DRIP-Seq rv Sutormin'].values
    
    Intervals_occupation(Intervals_info_left_fw, Intervals_info_left_rv, 'DRIP_LF', 'DRIP_LR', [1, 3500, 500, 1, 3500, 0.35, 150, 1.35], outpath)   
    Intervals_occupation(Intervals_info_right_fw, Intervals_info_right_rv, 'DRIP_RF', 'DRIP_RR', [1, 3500, 500, 1, 3500, 0.35, 150, 1.35], outpath)   
    
    #RNA-Seq signal on different replichores and strands.
    Intervals_info_right_fwg=Intervals_info.loc[((((Intervals_info['Start']>0) & (Intervals_info['Start']<Ter_position)) | (Intervals_info['Start']>Ori_position)) & (Intervals_info['Strand']=='+')), 'Expression_E'].values
    Intervals_info_left_fwg=Intervals_info.loc[(((Intervals_info['Start']>Ter_position) & (Intervals_info['Start']<Ori_position))  & (Intervals_info['Strand']=='+')), 'Expression_E'].values
    Intervals_info_right_rvg=Intervals_info.loc[(((Intervals_info['Start']>0) & (Intervals_info['Start']<Ter_position)) | (Intervals_info['Start']>Ori_position) & (Intervals_info['Strand']=='-')), 'Expression_E'].values
    Intervals_info_left_rvg=Intervals_info.loc[(((Intervals_info['Start']>Ter_position) & (Intervals_info['Start']<Ori_position)) & (Intervals_info['Strand']=='-')), 'Expression_E'].values
    
    Intervals_occupation(Intervals_info_right_fwg, Intervals_info_right_rvg, 'RNA_RF', 'RNA_RR', [1, 75000, 1000, 1, 75000, 0.35, 1500, 1.35], outpath)   
    Intervals_occupation(Intervals_info_left_fwg, Intervals_info_left_rvg, 'RNA_LF', 'RNA_LR', [1, 75000, 1000, 1, 75000, 0.35, 1500, 1.35], outpath) 
    
    #Cont nuber of genes laying on different replichores and with different orientation.
    print(Intervals_info.loc[(((Intervals_info['Start']>0) & (Intervals_info['Start']<Ter_position)) | (Intervals_info['Start']>Ori_position)) & (Intervals_info['Strand']=='+')].shape)
    print(Intervals_info.loc[(((Intervals_info['Start']>0) & (Intervals_info['Start']<Ter_position)) | (Intervals_info['Start']>Ori_position)) & (Intervals_info['Strand']=='-')].shape)
    print(Intervals_info.loc[((Intervals_info['Start']>Ter_position) & (Intervals_info['Start']<Ori_position)) & (Intervals_info['Strand']=='+')].shape)
    print(Intervals_info.loc[((Intervals_info['Start']>Ter_position) & (Intervals_info['Start']<Ori_position)) & (Intervals_info['Strand']=='-')].shape)  
    
    #EcTopoI intervals.
    #Intervals_occupation(intervals_param_ar_tg_2, cont_char_masked_2, name_of_intervals_1, name_of_cont_char_2, [0, 21, 2, -1, 21, 0.35, 11, 1.35], outpath)   
    #RNApol intervals.
    #Intervals_occupation(intervals_param_ar_tg_2_3, cont_char_masked_2_3, name_of_intervals_2, name_of_cont_char_3, [0, 5, 1, -0.2, 4.5, 0.35, 2.1, 1.35], outpath)
    return

func_wrapper(path_to_intervals_data, Outputpath)

print('Script ended its work succesfully!') 