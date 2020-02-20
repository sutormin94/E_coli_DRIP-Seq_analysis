###############################################
##Dmitry Sutormin, 2019##
##TopoA ChIP-Seq analysis##

#Takes narrowPeaks files with ChIP-Seq peaks identified in different biological replicas,
#return narrowPeak only with reproducible peaks.
###############################################

#######
#Packages to be imported.
#######

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import stats
from scipy.stats import binom, poisson

#Path to the all-containing table.
Data_path="F:\Signal_over_TUs\Signal_of_TUs_tab_all\All_genes\Signal_over_TUs_and_TF_syns_Mem.txt"
#Path to the output histogram.
Outpath="F:\TopoI_ChIP-Seq\Ec_TopoI_data\Figures\Membrane_Promoter_complexity_and_FE\\"



#######
#Plot distributions.
#######

def plot_signal_dist(S1, S2, S3, S4, S1a, S2a, S3a, S4a, outpath):
    all_values=pd.concat([S1, S2, S3, S4, S1a, S2a, S3a, S4a], ignore_index=True)
    min_FE=min(all_values)
    max_FE=max(all_values)
    #Plot distribution of S1.
    fig=plt.figure(figsize=(15, 5), dpi=100)
    bins0=np.arange(min_FE, max_FE, 0.2)
    plot0=plt.subplot2grid((2,4),(0,0), rowspan=1, colspan=1)
    plot0.hist(S1, bins0, color='#ff878b', edgecolor='black', alpha=0.8)
    plot0.annotate(f'Mean FE={round(np.mean(S1),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot0.set_yscale('log')
    plot0.set_xlabel('Fold enrichment', size=13)
    plot0.set_ylabel('Number of genes', size=13)
    plot0.set_title('All genes -Rif', size=15)  
    #Plot distribution of S1a.
    plot0a=plt.subplot2grid((2,4),(1,0), rowspan=1, colspan=1)
    plot0a.hist(S1a, bins0, color='#ff878b', edgecolor='black', alpha=0.8)
    plot0a.annotate(f'Mean FE={round(np.mean(S1a),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot0a.set_yscale('log')
    plot0a.set_xlabel('Fold enrichment', size=13)
    plot0a.set_ylabel('Number of genes', size=13)
    plot0a.set_title('All genes +Rif', size=15)   
    
    #Plot distribution of S2.  
    plot1=plt.subplot2grid((2,4),(0,1), rowspan=1, colspan=1)     
    plot1.hist(S2, bins0, color='#ffce91', edgecolor='black', alpha=0.5)
    plot1.annotate(f'Mean FE={round(np.mean(S2),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot1.set_yscale('log')
    plot1.set_xlabel('Fold enrichment', size=13)
    plot1.set_ylabel('Number of genes', size=13)
    plot1.set_title('Membrane proteins genes \n-Rif', size=15) 
    #Plot distribution of S2a.
    plot1a=plt.subplot2grid((2,4),(1,1), rowspan=1, colspan=1)     
    plot1a.hist(S2a, bins0, color='#ffce91', edgecolor='black', alpha=0.5)
    plot1a.annotate(f'Mean FE={round(np.mean(S2a),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot1a.set_yscale('log')
    plot1a.set_xlabel('Fold enrichment', size=13)
    plot1a.set_ylabel('Number of genes', size=13)
    plot1a.set_title('Membrane proteins genes \n+Rif', size=15) 
    
    #Plot distribution of Ded FE values over IG.
    plot2=plt.subplot2grid((2,4),(0,2), rowspan=1, colspan=1) 
    plot2.hist(S3, bins0, color='#7FCE79', edgecolor='black', alpha=0.5)
    plot2.annotate(f'Mean FE={round(np.mean(S3),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot2.set_yscale('log')
    plot2.set_xlabel('Fold enrichment', size=13)
    plot2.set_ylabel('Number of genes', size=13)
    plot2.set_title('Complex promoter genes \n-Rif', size=15)  
    #Plot distribution of Ded FE values over GB.
    plot2a=plt.subplot2grid((2,4),(1,2), rowspan=1, colspan=1) 
    plot2a.hist(S3a, bins0, color='#7FCE79', edgecolor='black', alpha=0.5)
    plot2a.annotate(f'Mean FE={round(np.mean(S3a),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot2a.set_yscale('log')
    plot2a.set_xlabel('Fold enrichment', size=13)
    plot2a.set_ylabel('Number of genes', size=13)
    plot2a.set_title('Complex promoter genes \n+Rif', size=15)      
    
    #Plot distribution of Ded FE values over IG.
    plot3=plt.subplot2grid((2,4),(0,3), rowspan=1, colspan=1) 
    plot3.hist(S4, bins0, color='#7FCE79', edgecolor='black', alpha=0.5)
    plot3.annotate(f'Mean FE={round(np.mean(S4),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot3.set_yscale('log')
    plot3.set_xlabel('Fold enrichment', size=13)
    plot3.set_ylabel('Number of genes', size=13)
    plot3.set_title('Complex promoter genes \n-Rif', size=15)  
    #Plot distribution of Ded FE values over GB.
    plot3a=plt.subplot2grid((2,4),(1,3), rowspan=1, colspan=1) 
    plot3a.hist(S4a, bins0, color='#7FCE79', edgecolor='black', alpha=0.5)
    plot3a.annotate(f'Mean FE={round(np.mean(S4a),2)}', xy=(0.40, 0.8), xycoords='axes fraction', size=12)
    plot3a.set_yscale('log')
    plot3a.set_xlabel('Fold enrichment', size=13)
    plot3a.set_ylabel('Number of genes', size=13)
    plot3a.set_title('Complex promoter genes \n+Rif', size=15)      
        
    
    plt.tight_layout()
    plt.show()
    #plt.savefig(outpath+'Membrane_loc_promoter_complexity_and_TopoA_FE.png', dpi=300, figsize=(15, 10))
    #plt.close()     
    return

#########
#Analysis of the enrichment of membrane protein encoding genes and genes with complex promoters among 
#genes enriched with TopoA in -Rif conditions.
#########

def Gene_prop_stat(Depth_of_sampling, Confidence, Data_frame_all, stat_method, Feature, outpath):
    #Compare real distribution of data values with stat model.
    fig, ax=plt.subplots(1,1)
    #Real data.
    bins=np.linspace(min(Data_frame_all[Feature]), max(Data_frame_all[Feature]), 100)
    ax.hist(Data_frame_all[Feature], bins, density=True, histtype='stepfilled', alpha=0.4, label='Real data distibution')
    Prop_Expected=round(np.mean(Data_frame_all[Feature]),2)
    #Model distribution.
    if stat_method=='poisson':
        x=np.arange(poisson.ppf(0.01, Prop_Expected), poisson.ppf(0.99, Prop_Expected))
        ax.plot(x, poisson.pmf(x, Prop_Expected), 'r-', lw=3, alpha=1, label=f'Model distibution\nmu={Prop_Expected}')
    elif stat_method=='gamma':
        a, b, c=stats.gamma.fit(Data_frame_all[Feature])
        gamma_pdf=stats.gamma.pdf(bins, a, b, c)
        ax.plot(bins, gamma_pdf, 'r-', lw=3, alpha=1, label=f'Model distibution gamma\na={round(a,2)}, b={round(b,2)}, c={round(c,2)}')
    elif stat_method=='norm':
        m, s=stats.norm.fit(Data_frame_all[Feature])
        norm_pdf=stats.norm.pdf(bins, m, s)
        ax.plot(bins, norm_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution normal\nmean={round(m,2)}, std={round(s,2)}')
    elif stat_method=='f':
        a, b, c, d=stats.f.fit(Data_frame_all[Feature])
        f_pdf=stats.f.pdf(bins, a, b, c, d)
        ax.plot(bins, f_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution f\na={round(a,2)}, b={round(b,2)}, c={round(c,2)}, d={round(d)}')
    elif stat_method=='expon':
        a, b=stats.expon.fit(Data_frame_all[Feature])
        expon_pdf=stats.expon.pdf(bins, a, b)
        ax.plot(bins, expon_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution expon\na={round(a,4)}, b={round(b,4)}')
    elif stat_method=='genpareto':
        a, b, c=stats.genpareto.fit(Data_frame_all[Feature])
        genpareto_pdf=stats.genpareto.pdf(bins, a, b, c)
        ax.plot(bins, genpareto_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution genpareto\na={round(a,2)}, b={round(b,2)}, c={round(c,2)}')    
    elif stat_method=='gilbrat':
        a, b=stats.gilbrat.fit(Data_frame_all[Feature])
        gilbrat_pdf=stats.gilbrat.pdf(bins, a, b)
        ax.plot(bins, gilbrat_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution gilbrat\na={round(a,2)}, b={round(b,2)}')   
    elif stat_method=='halfcauchy':
        a, b=stats.halfcauchy.fit(Data_frame_all[Feature])
        halfcauchy_pdf=stats.halfcauchy.pdf(bins, a, b)
        ax.plot(bins, halfcauchy_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution halfcauchy\na={round(a,2)}, b={round(b,2)}')          
    elif stat_method=='halfgennorm':
        a, b, c=stats.halfgennorm.fit(Data_frame_all[Feature])
        halfgennorm_pdf=stats.halfgennorm.pdf(bins, a, b, c)
        ax.plot(bins, halfgennorm_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution halfgennorm\na={round(a,2)}, b={round(b,2)}, c={round(c,2)}')      
    elif stat_method=='lomax':
        a, b, c=stats.lomax.fit(Data_frame_all[Feature])
        lomax_pdf=stats.lomax.pdf(bins, a, b, c)
        ax.plot(bins, lomax_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution lomax\na={round(a,2)}, b={round(b,2)}, c={round(c,2)}')    
    elif stat_method=='pareto':
        a, b, c=stats.pareto.fit(Data_frame_all[Feature])
        pareto_pdf=stats.pareto.pdf(bins, a, b, c)
        ax.plot(bins, pareto_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution pareto\na={round(a,2)}, b={round(b,2)}, c={round(c,2)}')    
    elif stat_method=='reciprocal':
        a, b, c, d=stats.reciprocal.fit(Data_frame_all[Feature])
        reciprocal_pdf=stats.reciprocal.pdf(bins, a, b, c, d)
        ax.plot(bins, reciprocal_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution reciprocal\na={round(a,2)}, b={round(b,2)}, c={round(c,2)}, d={round(d,2)}')     
    elif stat_method=='wald':
        a, b=stats.wald.fit(Data_frame_all[Feature])
        wald_pdf=stats.wald.pdf(bins, a, b)
        ax.plot(bins, wald_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution wald\na={round(a,2)}, b={round(b,2)}')       
    elif stat_method=='genlogistic':
        a, b, c=stats.genlogistic.fit(Data_frame_all[Feature])
        genlogistic_pdf=stats.genlogistic.pdf(bins, a, b, c)
        ax.plot(bins, genlogistic_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution genlogistic\na={round(a,2)}, b={round(b,2)}, c={round(c,2)}')        
    elif stat_method=='mielke':
        a, b, c, d=stats.mielke.fit(Data_frame_all[Feature])
        mielke_pdf=stats.mielke.pdf(bins, a, b, c, d)
        ax.plot(bins, mielke_pdf, 'r-', lw=3, alpha=0.5, label=f'Model distibution mielke\na={round(a,2)}, b={round(b,2)}, c={round(c,2)}, d={round(d,2)}')     
    ax.legend(loc='best')
    plt.show()    
    #plt.savefig(f'{outpath}\{Feature}_distribution_fitting_with_{stat_method}.png', dpi=300)
    
    
    
    #Store stat data.
    Window, Expected, Low, High, Top, Bot, Stat=[], [], [], [], [], [], []
    #Inint statistics.
    if stat_method=='poisson':
        Prop_interval=poisson.interval(Confidence, Prop_Expected)
    elif stat_method=='norm':
        Prop_interval=stats.norm.interval(Confidence, m, s)    
    elif stat_method=='gamma':
        Prop_interval=stats.gamma.interval(Confidence, a, b, c)
    elif stat_method=='f':
        Prop_interval=stats.f.interval(Confidence, a, b, c, d)
    elif stat_method=='expon':
        Prop_interval=stats.expon.interval(Confidence, a, b)
    elif stat_method=='genpareto':
        Prop_interval=stats.genpareto.interval(Confidence, a, b, c)
    elif stat_method=='gilbrat':
        Prop_interval=stats.gilbrat.interval(Confidence, a, b)   
    elif stat_method=='halfcauchy':
        Prop_interval=stats.halfcauchy.interval(Confidence, a, b)  
    elif stat_method=='halfgennorm':
        Prop_interval=stats.halfgennorm.interval(Confidence, a, b, c)    
    elif stat_method=='lomax':
        Prop_interval=stats.lomax.interval(Confidence, a, b, c)    
    elif stat_method=='pareto':
        Prop_interval=stats.pareto.interval(Confidence, a, b, c) 
    elif stat_method=='reciprocal':
        Prop_interval=stats.reciprocal.interval(Confidence, a, b, c, d)      
    elif stat_method=='wald':
        Prop_interval=stats.wald.interval(Confidence, a, b)      
    elif stat_method=='genlogistic':
        Prop_interval=stats.genlogistic.interval(Confidence, a, b, c) 
    elif stat_method=='mielke':
        Prop_interval=stats.mielke.interval(Confidence, a, b, c, d)      
    #Calculate enrichment.
    for Sampling in range(1,Depth_of_sampling):
        Window.append(Sampling)
        Expected.append(Prop_Expected)
        Low.append(Prop_interval[0]), High.append(Prop_interval[1])
        
        Top_FE=Data_frame_all.head(Sampling)
        Prop_Top=np.mean(Top_FE[Feature])
        Top.append(Prop_Top)
        
        Bot_FE=Data_frame_all.tail(Sampling)
        Prop_Bot=np.mean(Bot_FE[Feature])
        Bot.append(Prop_Bot)
        
        if stat_method=='norm':
            print('Perform t-test to compare means.')
            T_test=stats.ttest_ind(Top_FE[Feature], Bot_FE[Feature])
            print(T_test)
            Stat.append(float(T_test[1]))
        else: 
            MWU=stats.mannwhitneyu(Top_FE[Feature], Bot_FE[Feature])
            print('Perform Mann-Whithney test to compare samples.')
            print(MWU)
            Stat.append(float(MWU[1]))
        
        print(f'Mean {Feature} level: {Prop_Expected}')
        print(f'{Feature} level of Top genes: {Prop_Top}\n{Feature} level of Bot genes: {Prop_Bot}')           
        

    Stat=np.array(Stat)
    return Window, Expected, Low, High, Top, Bot, Stat


#######
#Plot MG and CP stat analysis.
#######

def Mem_CP_among_Top_Bot(Window, Expected, Low, High, Top, Bot, confidence, char, outpath):
    fig=plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)
    plot1.plot(Window, Low, linestyle=":", color="grey", linewidth=0.5, zorder=7) #Lower confidential border.
    plot1.plot(Window, High, linestyle=":", color="grey", linewidth=0.5, zorder=8) #Upper confidential border.   
    plot1.plot(Window, Expected, linestyle="-", color="red", linewidth=2, label='Expected number of genes', zorder=9) #Expected (mean).   
    plot1.fill_between(Window, Low, High, facecolor='blue', alpha=0.3, zorder=10) #Fill confident interval.
    plot1.plot(Window, Top, linestyle="-", color="black", linewidth=1, label=f'Number of {char} genes among high-TopoA-signal', zorder=12) #Number of genes among Top.  
    plot1.plot(Window, Bot, linestyle="--", color="black", linewidth=1, label=f'Number of {char} genes among low-TopoA-signal',  zorder=11) #Number of genes among Bot. 
    plot1.annotate(f'Confidence={confidence}', xy=(0.02, 0.75), xycoords='axes fraction', size=15)
    plot1.set_xlabel('Number of genes', size=20)
    plot1.set_ylabel(f'Number of {char} genes', size=20)
    plot1.set_title(f'{char} genes and TopoA signal', size=30)        
    plot1.legend(fontsize=12)    
  
    plt.show()
    #plt.savefig(f'{outpath}Enrichment_of_{char}_genes_among_top_and_bot_{len(Window)}_of_+Rif_TopoA_signal_confidence_{confidence}.png', dpi=300, figsize=(10, 6))
    #plt.close()
    return


#######
#Plot properties (GC%, Expression, etc) of Top and Bot genes sorted by TopoA FE.
#######

def Prop_among_Top_Bot(Window, Expected, Low, High, Top, Bot, Stat, confidence, stat_method, char, outpath):
    fig=plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)   
    plot1.plot(Window, Expected, linestyle="-", color="red", linewidth=2, label=f'Mean {char}', zorder=9) #Expected (mean). 
    #plot1.plot(Window, Low, linestyle=":", color="grey", linewidth=0.5, zorder=7) #Lower confidential border.
    #plot1.plot(Window, High, linestyle=":", color="grey", linewidth=0.5, zorder=8) #Upper confidential border. 
    #plot1.fill_between(Window, Low, High, facecolor='blue', alpha=0.3, zorder=10) #Fill confident interval.
    Top_masked=np.ma.masked_where(Stat>((1-confidence)/2), Top)
    Bot_masked=np.ma.masked_where(Stat>((1-confidence)/2), Bot)
    plot1.fill_between(Window, Bot_masked, Top_masked, facecolor='orange', alpha=0.3, zorder=11) #Fill confident differences between Top and Bot.
    plot1.plot(Window, Top, linestyle="-", color="black", linewidth=1, label=f'Mean {char} of high-TopoA-signal genes', zorder=12) #Signal for Top (mean).  
    plot1.plot(Window, Bot, linestyle="--", color="black", linewidth=1, label=f'Mean {char} of low-TopoA-signal genes',  zorder=11) #Signal for Bot (mean).
    plot1.annotate(f'Confidence ({stat_method}) ={confidence}', xy=(0.43, 0.25), xycoords='axes fraction', size=15)
    plot1.set_xlabel('Number of genes', size=20)
    plot1.set_ylabel(f'{char} of genes', size=20)
    plot1.set_title(f'{char} of genes and TopoA signal', size=30)        
    plot1.legend(loc='best', fontsize=12)    
  
    plt.show()
    #plt.savefig(f'{outpath}{char}_of_genes_top_and_bot_{len(Window)}_of_+Rif_TopoA_conf_{confidence}_with_{stat_method}.png', dpi=300, figsize=(10, 6))
    #plt.close()
    return

#######
#Wrap data.
#######

def wrapper_gene_groups(pathin, file_type, pathout):
    #Read data.
    if file_type=='xlsx':
        Data_frame=pd.read_excel(pathin, header=0)
    elif file_type=='txt':
        Data_frame=pd.read_csv(pathin, sep='\t', header=0)
    #All genes.
    All_noRif=Data_frame[Data_frame['TopoA -Rif_FE_GB'].notnull()]['TopoA -Rif_FE_GB']
    All_Rif=Data_frame[Data_frame['TopoA +Rif_FE_GB'].notnull()]['TopoA +Rif_FE_GB']
    #Return data based on membrane localization.
    In_membrane_dataframe=Data_frame.loc[Data_frame['Membrane_localization']!='-']
    Membrane_TopoA_noRif=In_membrane_dataframe[In_membrane_dataframe['TopoA -Rif_FE_GB'].notnull()]['TopoA -Rif_FE_GB']
    Membrane_TopoA_Rif=In_membrane_dataframe[In_membrane_dataframe['TopoA +Rif_FE_GB'].notnull()]['TopoA +Rif_FE_GB']
    #Return data based on promoter complexity.
    TF_known_dataframe=Data_frame[Data_frame['Number_of_TF_sites'].notnull()]
    Complex_promoter_TopoA_noRif=TF_known_dataframe[TF_known_dataframe['TopoA -Rif_FE_GB'].notnull()]['TopoA -Rif_FE_GB']
    Complex_promoter_TopoA_Rif=TF_known_dataframe[TF_known_dataframe['TopoA +Rif_FE_GB'].notnull()]['TopoA +Rif_FE_GB']
    #Select data based on membrane localization & promoter complexity.
    In_mem_complex_prom_dataframe=In_membrane_dataframe[In_membrane_dataframe['Number_of_TF_sites'].notnull()]
    IMCP_TopoA_noRif=In_mem_complex_prom_dataframe[In_mem_complex_prom_dataframe['TopoA -Rif_FE_GB'].notnull()]['TopoA -Rif_FE_GB']
    IMCP_TopoA_Rif=In_mem_complex_prom_dataframe[In_mem_complex_prom_dataframe['TopoA +Rif_FE_GB'].notnull()]['TopoA +Rif_FE_GB']    
    
    Num_of_genes=len(All_noRif)
    Num_of_genes_M=len(Membrane_TopoA_noRif)
    Num_of_genes_CP=len(Complex_promoter_TopoA_noRif)
    Num_of_genes_MCP=len(IMCP_TopoA_noRif)
    
    print(f'Total number of genes: {Num_of_genes}\nNumber of genes encoding membrane proteins: {Num_of_genes_M}')
    print(f'Number of genes with complex promoters: {Num_of_genes_CP}')
    print(f'Number of genes encoding mem proteins and having complex promoters: {Num_of_genes_MCP}')
    plot_signal_dist(All_noRif, Membrane_TopoA_noRif, Complex_promoter_TopoA_noRif, IMCP_TopoA_noRif, 
                     All_Rif, Membrane_TopoA_Rif, Complex_promoter_TopoA_Rif, IMCP_TopoA_Rif, 
                     pathout)
    
    #Analysis of the enrichment of membrane protein encoding genes and genes with complex promoters among 
    #genes enriched with TopoA in -Rif conditions.
    Data_frame_NN=Data_frame[Data_frame['TopoA +Rif_FE_GB'].notnull()]
    Data_frame_NN_SNR=Data_frame_NN.sort_values('TopoA +Rif_FE_GB', axis=0, ascending=False, na_position='last')
    Depth_of_sampling=500
    Confidence=0.999    
    #Store stat data.
    Window=[]
    Expected_CP, Low_CP, High_CP, Top_CP, Bot_CP=[], [], [], [], []
    Expected_M, Low_M, High_M, Top_M, Bot_M=[], [], [], [], []
    #Calculate enrichment.
    for Sampling in range(1,Depth_of_sampling):
        Window.append(Sampling)
        Num_of_CP_Sampl_Expected=round(Sampling*(Num_of_genes_CP/Num_of_genes),0)
        Num_of_MG_Sampl_Expected=round(Sampling*(Num_of_genes_M/Num_of_genes),0)
        Expected_CP.append(Num_of_CP_Sampl_Expected), Expected_M.append(Num_of_MG_Sampl_Expected)
        
        Top_noRif_FE=Data_frame_NN_SNR.head(Sampling)
        Num_of_CP_Top=len(Top_noRif_FE[Top_noRif_FE['Number_of_TF_sites'].notnull()])
        Num_of_MG_Top=len(Top_noRif_FE.loc[Top_noRif_FE['Membrane_localization']!='-'])
        Top_CP.append(Num_of_CP_Top), Top_M.append(Num_of_MG_Top)
        
        Bot_noRif_FE=Data_frame_NN_SNR.tail(Sampling)
        Num_of_CP_Bot=len(Bot_noRif_FE[Bot_noRif_FE['Number_of_TF_sites'].notnull()])
        Num_of_MG_Bot=len(Bot_noRif_FE.loc[Bot_noRif_FE['Membrane_localization']!='-'])
        Bot_CP.append(Num_of_CP_Bot), Bot_M.append(Num_of_MG_Bot)
    
        interval_CP=binom.interval(Confidence, Sampling, Num_of_genes_CP/Num_of_genes)
        interval_M=binom.interval(Confidence, Sampling, Num_of_genes_M/Num_of_genes)
        Low_CP.append(interval_CP[0]), Low_M.append(interval_M[0]), High_CP.append(interval_CP[1]), High_M.append(interval_M[1])
    
        print(f'Lower conf interval genes CP: {interval_CP[0]}\nExpected num of CP genes: {Num_of_CP_Sampl_Expected}\nHigher conf interval genes CP: {interval_CP[1]}')
        print(f'Number of CP genes among Top: {Num_of_CP_Top}\nNumber of CP genes among Bot: {Num_of_CP_Bot}')
        print(f'Lower conf interval genes M: {interval_M[0]}\nExpected num of M genes: {Num_of_MG_Sampl_Expected}\nHigher conf interval genes M: {interval_M[1]}')
        print(f'Number of M genes among Top: {Num_of_MG_Top}\nNumber of M genes among Bot: {Num_of_MG_Bot}')     
        
    Mem_CP_among_Top_Bot(Window, Expected_CP, Low_CP, High_CP, Top_CP, Bot_CP, Confidence, 'CP', pathout)
    Mem_CP_among_Top_Bot(Window, Expected_M, Low_M, High_M, Top_M, Bot_M, Confidence, 'M', pathout)
    
    #Analysis of different parameters (GC%, Expression level, etc) among 
    #genes enriched with TopoA in -Rif conditions.   
    Confidence=0.999
    Depth_of_sampling=500
    Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp=Gene_prop_stat(Depth_of_sampling, Confidence, Data_frame_NN_SNR, 'wald', 'Expression', pathout)   
    Prop_among_Top_Bot(Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp, Confidence, 'Mann-Whitney U', 'Expression', pathout)
    Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp=Gene_prop_stat(Depth_of_sampling, Confidence, Data_frame_NN_SNR, 'genlogistic', 'GC_FE_GB', pathout)   
    Prop_among_Top_Bot(Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp, Confidence, 'T-test', 'GC_FE_GB', pathout) 
    Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp=Gene_prop_stat(Depth_of_sampling, Confidence, Data_frame_NN_SNR, 'wald', 'RpoB_FE_GB', pathout)   
    Prop_among_Top_Bot(Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp, Confidence, 'Mann-Whitney U', 'RpoB_FE_GB', pathout) 
    Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp=Gene_prop_stat(Depth_of_sampling, Confidence, Data_frame_NN_SNR, 'wald', 'PolSofi_FE_GB', pathout)   
    Prop_among_Top_Bot(Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp, Confidence, 'Mann-Whitney U', 'PolSofi_FE_GB', pathout)       
    #Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp=Gene_prop_stat(Depth_of_sampling, Confidence, Data_frame_NN_SNR, 'gamma', 'Gyrase -Rif_FE_GB', pathout)   
    #Prop_among_Top_Bot(Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp, Confidence, 'Mann-Whitney U', 'Gyrase -Rif_FE_GB', pathout)  
    Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp=Gene_prop_stat(Depth_of_sampling, Confidence, Data_frame_NN_SNR, 'gamma', 'Gyrase +Rif_FE_GB', pathout)   
    Prop_among_Top_Bot(Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp, Confidence, 'Mann-Whitney U', 'Gyrase +Rif_FE_GB', pathout)    
    Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp=Gene_prop_stat(Depth_of_sampling, Confidence, Data_frame_NN_SNR, 'gamma', 'TopoIV_FE_GB', pathout)   
    Prop_among_Top_Bot(Window, Expected_E, Low_E, High_E, Top_Exp, Bot_Exp, Stat_Exp, Confidence, 'Mann-Whitney U', 'TopoIV_FE_GB', pathout)     
    return

wrapper_gene_groups(Data_path, 'txt', Outpath)