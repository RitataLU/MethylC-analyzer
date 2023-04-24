# 2022.01.15

import sys
import subprocess, sys
import pandas as pd
import numpy as np
#import seaborn as sns
import os
import random
#import matplotlib
#if os.environ.get('DISPLAY','') == '':
#    print('no display found. Using non-interactive Agg backend')
#    matplotlib.use('Agg')

#import matplotlib.pyplot as plt
#from matplotlib import rcParams
#import matplotlib as mpl
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
import math
from math import log
from scipy import stats
import scipy.stats as ss
import time
import argparse
import glob
import pyBigWig
from scipy.stats import rankdata
from scipy import stats

parser = argparse.ArgumentParser()

parser.add_argument("-a",help="Name of experimental group",dest='group1',default='')
parser.add_argument("-b",help="Name of control group",dest='group2',default='')
parser.add_argument("-d",help="Minimum depth of sites. Default=4",dest='depth',default=4)
parser.add_argument("-r",help="Size of region. Default=500",dest='region',default=500)
parser.add_argument("-q",help="Minimum sites within a region. Default=4",dest='qualified',default=4)
parser.add_argument("-context",help="Context used. Default=CG",dest='context',default='CG')
parser.add_argument("-hc",help="Methylation cutoff of PCA & Heatmap. Default = 0.2",dest='heatmap_cutoff',default=0.2)
#parser.add_argument("-hcgc",help="PCA & Heatmap_CG_cutoff",dest='heatmap_cg_cutoff',default=0.2)
#parser.add_argument("-hchgc",help="PCA & Heatmap_CHG_cutoff",dest='heatmap_chg_cutoff',default=0.2)
#parser.add_argument("-hchhc",help="PCA & Heatmap_CHH_cutoff",dest='heatmap_chh_cutoff',default=0.2)
parser.add_argument("-dmrc",help="Methylation cutoff of DMR. Default = 0.1",dest='dmr_cutoff',default=0.1)
#parser.add_argument("-dmrcg",help="DMR_CG_cutoff",dest='dmr_cg_cutoff',default=0.1)
#parser.add_argument("-dmrchg",help="DMR_CHG_cutoff",dest='dmr_chg_cutoff',default=0.1)
#parser.add_argument("-dmrchh",help="DMR_CHH_cutoff",dest='dmr_chh_cutoff',default=0.1)
parser.add_argument("-test",help="DMR testing method. 0:TTest, 1:KS, 2:MWU. Default=0",dest='testMethod',default=0)
parser.add_argument("-pvalue",help="p-value cutoff for identifying DMR. Default = 0.05",dest='pvalue',default=0.05)
#parser.add_argument("-fdr",help="fdr for identifying DMR",dest='fdr',default=0.05)
parser.add_argument("-bs",help="Bin size of chrView and Metaplot. Default = 1000000",dest='bin_size',default=1000000)
parser.add_argument("-p",help="promoter_size",dest='promoter_size',default=2000)
parser.add_argument("command",help="commands of MethylC-Analyser")
parser.add_argument("samples_list",help="samples CGmap description")
parser.add_argument("input_gtf_file",help="path of gene annotation")
parser.add_argument("path_to_files",help="path of all target files")

args = parser.parse_args()

depth=int(args.depth)
region=int(args.region)
qualifiedSite=int(args.qualified)
pca_heat_cut=float(args.heatmap_cutoff)
#pca_heat_cg_cut=float(args.heatmap_cg_cutoff)
#pca_heat_chh_cut=float(args.heatmap_chh_cutoff)
#pca_heat_chg_cut=float(args.heatmap_chg_cutoff)
dmr_cut=float(args.dmr_cutoff)
#dmr_CG_cut=float(args.dmr_cg_cutoff)
#dmr_CHH_cut=float(args.dmr_chh_cutoff)
#dmr_CHG_cut=float(args.dmr_chg_cutoff)
pvalue = float(args.pvalue)
testmethod = int(args.testMethod)
context = args.context.upper() if args.context.upper() in ['CG','CHG','CHH'] else 'CG'
#fdr = float(args.fdr)
binSize=int(args.bin_size)
promoter_size=int(args.promoter_size)
#input_gene=pd.read_csv(str(args.input_gtf_file),sep='\t',header=None)
path_to_files=str(args.path_to_files)
samples_list=path_to_files+str(args.samples_list)
input_gtf_file=path_to_files+str(args.input_gtf_file)
input_gene_name=input_gtf_file
command=str(args.command)
path_to_files=str(args.path_to_files)



#Choose Tools
# Heatmap_PCA=input('Heatmap & PCA Analysis?  ([y]/n): ') or 'y'
# DMR=input('Identify DMR?  ([y]/n): ') or 'y'
# DMG=input('Identify DMG?  ([y]/n): ') or 'y'
# Fold_Enrichment=input('Use Fold Enrichment Analysis?  ([y]/n): ') or 'y'
# ChrView=input('Chromosome View Analysis?  ([y]/n): ') or 'y'
# Metaplot=input('Metaplot Analysis?  ([y]/n): ') or 'y'

DMR_exp=args.group1
DMR_ctrl=args.group2

#all functions
#transfer txt to bed
# def barplot:
#     df = pd.read_csv("CommonRegion_"+context+".txt",sep="\t")

def bed_form(inputxt):
    subprocess.call('''awk '{if (NR!=1) print $1"\t"$2"\t"$3}' %s > %s'''%(inputxt, inputxt+".bed"), shell=True)
    return inputxt+".bed"

# def Find_DMR(context, cutoff):
#     file1=pd.read_csv("CommonRegion_"+context+".txt",sep="\t",dtype =
#             {0:str,1:int,2:int},index_col=[0,1,2])
#     # set types for each coulmn 
#     #file1.iloc[:,0] = str(file1.iloc[:,0])
#     #file1.iloc[:,1:3] = file1.iloc[:,1:3].astype(int)
#     #add small value (0.001) to prevent meth_level is 0 for statistics (ttest)
#     file1.iloc[:,0] = file1.iloc[:,0:].astype(float)+0.0001 
#     file1.iloc[:,file1.shape[1]-1] = file1.iloc[:,file1.shape[1]-1].astype(float)+0.0001 

#     PVALUE =[]
#     for i in file1.index:
#         expGp = file1.loc[i,expgroup].dropna()
#         ctrlGp = file1.loc[i,ctrlgroup].dropna()

#         Delta = expGp.mean()-ctrlGp.mean()
#         pval =  stats.ttest_ind(expGp , ctrlGp)[1]
#         pKS = stats.kstest(expGp , ctrlGp)[1]
#         pMWU = stats.mannwhitneyu(expGp , ctrlGp)[1]
#         p = pd.concat([pd.DataFrame(expGp),pd.DataFrame(ctrlGp)]).T
#         p['Delta'] = Delta
#         p['pval'] = pval
#         p['pKS'] = pKS
#         p['pMWU'] = pMWU
#         PVALUE.append(p)

#     merge = pd.concat(PVALUE)
#     sig = merge[merge.pval <= pvalue]
#     sig_all = sig[(sig.Delta >= cutoff) | (sig.Delta <= -1*cutoff)]
#     sig_all.iloc[:,0] = sig_all.iloc[:,0:].astype(float)-0.0001 
#     sig_all.iloc[:,file1.shape[1]-1] = sig_all.iloc[:,file1.shape[1]-1].astype(float)-0.0001 
#     sig_all.to_csv('DMR_'+context+'_all_'+str(cutoff)+'.txt', sep='\t')
#     sig_all[sig_all.Delta >0].to_csv('DMR_'+context+'_hyper_'+str(cutoff)+'.txt',sep='\t')
#     sig_all[sig_all.Delta <0].to_csv('DMR_'+context+'_hypo_'+ str(cutoff)+'.txt',sep='\t')

def Find_DMR2(context, cutoff, test_method):
#    union=pd.read_csv('Unionsite.txt',sep='\t',na_values='-')
#    expgroup = samples[samples[2] == 'WT'][0].to_list()
#    ctrlgroup = samples[samples[2] == 'met1'][0].to_list()
    chrs = union['chr'].unique()
    data_holder = []
    positions = []
    meanMaths = []
    pvalues = []
    for chromosome in chrs:
        subset = union[(union['context'] == context) & (union['chr'] == chromosome)]
        if subset.shape[0] == 0: continue 
        maxPos = subset['pos'].max()
        bins = range(0,maxPos,region)
        groups = subset.groupby(pd.cut(subset['pos'], bins))
        for sRange,sValues in groups:
            minDepth = sValues.iloc[:,3:].count().min()
            if minDepth >= qualifiedSite:
                expValue = sValues.loc[:,expgroup]
                expValue2 = [x for sublist in expValue.values for x in sublist if math.isnan(x) == False]
                ctrlValue = sValues.loc[:,ctrlgroup]
                ctrlValue2 = [x for sublist in ctrlValue.values for x in sublist if math.isnan(x) == False]
                pKS = stats.kstest(expValue2,ctrlValue2)[1]
                pMWU = stats.mannwhitneyu(expValue2,ctrlValue2)[1]
                meanMeth = sValues.iloc[:,3:].astype(float).mean()
                meanMeth2 = [("%.3f" %x) for x in meanMeth]
                meanMeth3 = meanMeth + [random.random() * 0.00001 for x in meanMeth]
                pTTest = 1.0 if meanMeth.sum() == 0 else stats.ttest_ind(meanMeth3.loc[expgroup],meanMeth3.loc[ctrlgroup])[1]
                methMeth = sValues.iloc[:,3:].astype(float).mean().tolist()
                start = sRange.left
                end = sRange.right
                deltaMean = meanMeth.loc[expgroup].mean() - meanMeth.loc[ctrlgroup].mean()
                out = [chromosome,start,end] + meanMeth2  + [deltaMean,pTTest,pKS,pMWU]
                data_holder.append(out)

    merge = pd.DataFrame(data_holder,columns = ['Chr','Start','End'] + samples[0].tolist() + ['DeltaMean','pTTest','pKS','pMWU'])
    tests_methods = {0:'pTTest',1:'pKS',2:'pMWU'}
    pvals = merge.loc[:,tests_methods[test_method]]
    sig = merge[pvals <= pvalue]
    sig_all = sig[(sig.DeltaMean >= cutoff) | (sig.DeltaMean <= -1*cutoff)]
    sig_all.to_csv(path_to_files +'DMR_'+context+'_all_'+str(cutoff)+'.txt', sep='\t',index = False)
    sig_all[sig_all.DeltaMean >0].to_csv(path_to_files +'DMR_'+context+'_hyper_'+str(cutoff)+'.txt',sep='\t',index = False)
    sig_all[sig_all.DeltaMean <0].to_csv(path_to_files +'DMR_'+context+'_hypo_'+ str(cutoff)+'.txt',sep='\t',index = False)

def check_dmgempty(file):

    if os.path.getsize(file) > 0:
        df = pd.read_csv(file,header=None,sep='\t')
        value = df[3].nunique()
    else:
        value = 0

    return value

def DMR_DMGPlot(Context,cutoff):

    DMR_hyper = pd.read_csv(path_to_files+"DMR_"+Context+"_hyper_"+str(cutoff)+".txt", sep='\t')
    DMR_hypo = pd.read_csv(path_to_files+"DMR_"+Context+"_hypo_"+str(cutoff)+".txt", sep='\t')

    R_hyper = len(DMR_hyper)
    R_hypo = len(DMR_hypo)


    Gg_hyper = check_dmgempty(path_to_files +"DMG_"+Context+"_hyper_"+str(cutoff)+"_Genebody_list.txt")
    Gg_hypo  = check_dmgempty(path_to_files +"DMG_"+Context+"_hypo_"+str(cutoff)+"_Genebody_list.txt")

    Gp_hyper = check_dmgempty(path_to_files +"DMG_"+Context+"_hyper_"+str(cutoff)+"_Promoter_list.txt")
    Gp_hypo = check_dmgempty(path_to_files +"DMG_"+Context+"_hypo_"+str(cutoff)+"_Promoter_list.txt")

    # DMG_hyper_genebody = pd.read_csv("DMG_"+Context+"_hyper_"+str(cutoff)+"_Genebody_list.txt", header=None,sep='\t')
    # DMG_hypo_genebody = pd.read_csv("DMG_"+Context+"_hypo_"+str(cutoff)+"_Genebody_list.txt", header=None,sep='\t')

    # DMG_hyper_promoter = pd.read_csv("DMG_"+Context+"_hyper_"+str(cutoff)+"_Promoter_list.txt",header=None, sep='\t')
    # DMG_hypo_promoter = pd.read_csv("DMG_"+Context+"_hypo_"+str(cutoff)+"_Promoter_list.txt", header=None,sep='\t')


    # Gg_hyper = DMG_hyper_genebody[3].nunique()
    # Gg_hypo = DMG_hypo_genebody[3].nunique()

    # Gp_hyper = DMG_hyper_promoter[3].nunique()
    # Gp_hypo = DMG_hypo_promoter[3].nunique()

    data = { "Category": ['DMR', 'DMR', 'DMG(genebody)','DMG(genebody)','DMG(promoter)','DMG(promoter)'],
  "Number": [R_hyper, R_hypo, Gg_hyper,Gg_hypo,Gp_hyper,Gp_hypo],
  "Direction":['hyper','hypo','hyper','hypo','hyper','hypo']}
    df = pd.DataFrame(data)

    # set plot style: grey grid in the background:
    sns.set_style("whitegrid")
    # sns.set_context("talk")
    sns.set_context("notebook", font_scale=2, rc={"lines.linewidth": 2.5})
    colors = ["#007A78", "#FFC745"]
    # Set the figure size
    plt.figure(figsize=(10,6)) 
    # grouped barplot
    sns.barplot(x="Category", y="Number", hue="Direction", data=df, errorbar=None,palette=colors)
    
    plt.legend(fontsize = 15, 
               bbox_to_anchor= (1.01, 1), 
               title="", 
               shadow = False, 
               facecolor = 'white')
    plt.xlabel('')
    plt.ylabel('Numbers')

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams["axes.labelsize"] = 40
    plt.savefig(path_to_files +'Summary_DMR_DMG_numbers_'+ Context+'_'+str(cutoff)+'.pdf',dpi=300,bbox_inches="tight")


#enrichment
#1. overlap with dmr
def overlap(bed1,bed2):
        overlap=subprocess.check_output("bedtools intersect -a %s -b %s |awk '{size+=$3 - $2}END{print size}' "%(bed1,bed2), shell=True)
        if len(overlap)<2:
                return 1
        else:
                return int(overlap)
#enrichment plot
def Enrichment(tag, Cut):
    savefile=pd.DataFrame(columns=["dmr","feature","overlap_size","dmr_size","feature_size","genome_size","enrichment"])
    dmr=path_to_files +"DMR_"+tag+"_all_"+str(Cut)+".txt"
    for i in annotation_name:
        feature=str(i)
        overlap_size=overlap(input_gene_name +'_'+i+'_merge.bed',bed_form(dmr))
        dmrlen=pd.read_csv(dmr,sep='\t')
        if len(dmrlen)<=1:
            dmr_size=1000000000000
        else:
            dmr_size=int(subprocess.check_output('''awk '{size+=$3-$2}END{print size}' %s'''%(bed_form(dmr)), shell=True))

        feature_size=overlap(bed_form(path_to_files +"CommonRegion_"+tag+".txt"),bed_form(input_gene_name +'_'+i+'_merge.bed'))
        genome_size=int(subprocess.check_output('''awk '{size+=$3-$2}END{print size}' %s'''%(bed_form(path_to_files +"CommonRegion_"+tag+".txt")), shell=True))*1.0

        enrichment= math.log(((float(overlap_size)/float(dmr_size))/(float(feature_size)/float(genome_size))),2)
        out = [dmr,feature,overlap_size,dmr_size,feature_size,genome_size,enrichment]
        savefile.loc[len(savefile)]=out
        plt.style.use('ggplot')
        annotationname = ['Promoter','Genebody','Exon','Intron','5UTR','CDS','3UTR','IGR']
        annotationname_index = range(len(annotationname))

    Eplt=[savefile['enrichment']]
    for value in Eplt:
        fig = plt.figure(figsize=(8,6)) 

        ax1 = fig.add_subplot(1,1,1)
        ax1.bar(annotationname_index, value,align='center',color='darkblue')
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        plt.xticks(annotationname_index, annotationname, rotation=90,fontsize=15)
        plt.ylabel("Fold Enrichment (log2)")
        plt.title(tag+'_Fold_Enrichment')
        plt.grid()
        plt.savefig(path_to_files+tag+'_Fold_Enrichment'+'.pdf',dpi=300,bbox_inches='tight')
        #plt.close(fig)

###generating metaplot_delta files
def Delta_Meta(context):

    #sample_list=pd.read_csv("chrView_list.txt",header=None,sep='\t')
    
    expData = samples[samples[2] == metaplot_exp][0]
    ctrlData = samples[samples[2] == metaplot_ctrl][0]
    # read chrview files in two group 

    exp_df = []
    for expfile in expData:
        exp_df.append(pd.read_csv(path_to_files+expfile+"_"+ context +".matrix.gz" ,sep='\t', skiprows=1,compression ='gzip',na_values='-',header =None))

    ctrl_df = []
    for ctrlfile in ctrlData:
        ctrl_df.append(pd.read_csv(path_to_files+ctrlfile+"_"+ context +".matrix.gz" ,sep='\t', skiprows=1,compression ='gzip',na_values='-',header =None))
    
    #sum all valuew in each row in 3 context
    mexp,mctrl=0,0
    for j in range(0,len(expData)):
        mexp+=exp_df[j].iloc[:,6:86]
    
    for k in range(0,len(ctrlData)):
        mctrl+=ctrl_df[k].iloc[:,6:86]
        

   # get average M level in each row (gene)
    mexp/=len(expData)
    mctrl/=len(ctrlData)


    ## delta m level in 3 context
    deltaexp =0

    delta=mexp-mctrl
    
    #generate delta report
    exp_df[0].iloc[:,6:86]=delta

    report = exp_df[0]
    report.to_csv(path_to_files+'metaplot_delta_'+context+'.txt',sep='\t',index=False)

    pd.DataFrame([[metaplot_exp+' - '+metaplot_ctrl,path_to_files+'metaplot_delta.txt']]).to_csv(path_to_files+'metaplot_delta_list.txt',sep='\t',index=False,header=None)

     
###processing start, generating common regions

contexts=["CG","CHG","CHH"]
#unionsite

combined = pd.DataFrame()
samples = pd.read_csv(samples_list,header=None,sep="\t")
#define groups
expgroup = samples[samples[2] == DMR_exp][0].to_list()
ctrlgroup = samples[samples[2] == DMR_ctrl][0].to_list()


for sample in samples.itertuples():
    print("Now processing " + sample[2])
    CGmap = pd.read_csv(path_to_files+sample[2], header=None,sep="\t",dtype =
            {0:str,2:str,3:str,5:str,6:str,7:int},usecols=[0,2,3,5,6,7],index_col=[0,1,2],compression='gzip')

    CGmap = CGmap.loc[CGmap[7] >= depth,5:5]

    sample_col = [sample[1]]
    #sample_col = [sample[1] + '_meth',sample[1] + '_mC',sample[1] + '_C']
    CGmap.columns = sample_col
    if len(combined) == 0:
        combined = CGmap
    else:
        combined = pd.merge(combined,CGmap,left_index=True,right_index=True,how = 'outer')


combined.reset_index(inplace=True)
combined.rename(columns = {0:'chr', 2:'pos', 3: 'context'},inplace=True)
combined[['pos']] = combined[['pos']].astype(int)
#process chromosome name
#deleted

#combined[['chr']] = combined[['chr']].astype(int)

combined = combined.sort_values(['chr','pos'],ascending=[True,True])
combined.to_csv(path_to_files+'Unionsite.txt',sep = '\t',na_rep='-',index=False)


#common region

union=pd.read_csv(path_to_files+'Unionsite.txt',sep='\t',na_values='-')
#union=combined
#1.unionsite --> combined
chrs = union['chr'].unique()


#context one by one
for cxt in contexts:
        outfile = path_to_files+ 'CommonRegion_' + cxt + '.txt'
        with open(outfile,'w') as of:
                of.write('\t'.join(['chr','start','end'] + union.columns.values.tolist()[3:]) + "\n")
                for chromosome in chrs:
                        subset = union[(union['context'] == cxt) & (union['chr'] == chromosome)]
                        if subset.shape[0] == 0: continue 
                        maxPos = subset['pos'].max()
                        bins = range(0,maxPos,region)
                        groups = subset.groupby(pd.cut(subset['pos'], bins))
                        for sRange,sValues in groups:
                            #rows, columns = sValues.shape
                            #with a region find out minimun sites, less than depthcutoff
                                minDepth = sValues.iloc[:,3:].count().min()
                                if minDepth >= qualifiedSite:
                                        methMeth = sValues.iloc[:,3:].astype(float).mean().tolist()
                                        methMeth2 = [("%.3f" %x) for x in methMeth]
                                        start = sRange.left
                                        end = sRange.right
                                        out = map(str,[chromosome,start,end] + methMeth2)
                                        outStr = '\t'.join(out)
                                        of.write(outStr + "\n")

                                else:
                                        pass

##The Average Methylaion level
cg = pd.read_csv(path_to_files+"CommonRegion_CG.txt",sep="\t")
chg = pd.read_csv(path_to_files+"CommonRegion_CHG.txt",sep="\t")
chh = pd.read_csv(path_to_files+"CommonRegion_CHH.txt",sep="\t")

end = cg.shape[1]
cgdf = pd.DataFrame(cg.iloc[:,3:end].mean())
cgdf['context'] = 'CG'
chgdf = pd.DataFrame(chg.iloc[:,3:end].mean())
chgdf['context'] = 'CHG'
chhdf = pd.DataFrame(chh.iloc[:,3:end].mean())
chhdf['context'] = 'CHH'
merge1 = pd.concat([cgdf,chgdf])
merge2 = pd.concat([merge1,chhdf])
merge2['sample'] = merge2.index
merge2.columns =['methylation level', 'context', 'sample']
merge2['methylation level'] = merge2['methylation level']*100
merge2['group'] = merge2.index.map(samples.set_index(0)[2])

#ploting barplot 
import seaborn as sns
import matplotlib
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib as mpl
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

sns.set_style("whitegrid")
# sns.set_context("talk")
sns.set_context("notebook", font_scale=2, rc={"lines.linewidth": 2.5})

plt.figure(figsize=(8,6)) 

# fig.set_size_inches(8,6)
g = sns.catplot(
    data=merge2, kind="bar",
    x="context", y= 'methylation level', hue='group', palette="dark", 
    alpha=.6, height=6, legend = False)

plt.legend(fontsize = 15, 
               bbox_to_anchor= (1.3, 1), 
               title="", 
               shadow = False, 
               facecolor = 'white')
ax = g.facet_axis(0,0)
for p in ax.patches:
    ax.text(p.get_x() - 0.001, 
            p.get_height() * 1.05, 
           '{0:.1f}'.format(p.get_height()),   #Used to format it K representation
            color='black', 
            rotation='horizontal', 
            size='x-small')

g.set_axis_labels("", "Methylation level (%)")
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["axes.labelsize"] = 40

plt.savefig(path_to_files +'Average_methylation_levels.pdf',dpi=300,bbox_inches="tight")

#heatmap_PCA
if(command=='Heatmap_PCA' or command=='all'):
    print ("*------------------------*")
    print ("|generating Heatmap $ PCA|")
    print ("*------------------------*")
    if context == 'CG':
        subprocess.call('''Rscript --slave /MethylC-analyzer/scripts/heatmap_PCA_all.R %s %s %s'''%(path_to_files + "CommonRegion_CG.txt",pca_heat_cut,path_to_files), shell=True)
    elif context == 'CHG':
        subprocess.call('''Rscript --slave /MethylC-analyzer/scripts/heatmap_PCA_all.R %s %s %s'''%(path_to_files +"CommonRegion_CHG.txt",pca_heat_cut,path_to_files), shell=True)
    elif context == 'CHH':
        subprocess.call('''Rscript --slave /MethylC-analyzer/scripts/heatmap_PCA_all.R %s %s %s'''%(path_to_files +"CommonRegion_CHH.txt",pca_heat_cut,path_to_files ), shell=True)

else:
    pass

# Identify DMR
if(command=='DMR' or command=='all'):

    print ("*---------------*")
    print ("|Identifying DMR|")
    print ("*---------------*")

    expgroup = samples[samples[2] == DMR_exp][0]
    ctrlgroup = samples[samples[2] == DMR_ctrl][0]

    Find_DMR2(context,dmr_cut,testmethod) 

else:
    pass


# preprocessing for DMG, fold enrichment
# generated BED for each genomic region


gene = "Genebody"
exon = "exons"
intron = "introns"
utr3 = "3utr"
utr5 = "5utr"
cds = "cds"
promoter = "Promoter"
igr = "IGR"
bed12=[exon,intron,utr5,cds,utr3]

annotation_name=[promoter,gene,exon,intron,utr5,cds,utr3,igr]

subprocess.call('''python /MethylC-analyzer/scripts/extract_transcript_regions.py -i %s -o %s --gtf'''%(input_gtf_file,input_gtf_file), shell=True)

#Convert this blockbed (bed12) to bed6|
for i in bed12:
    subprocess.call('''cat %s | bed12ToBed6 -i stdin -n > %s'''%(input_gtf_file+'_'+i+'.bed',input_gtf_file+'_'+i+'_bed6.bed'),shell=True)
    # print (i)

subprocess.call('''rm %s*3utr.bed %s*5utr.bed %s*_cds.bed %s*_exons.bed %s*_introns.bed '''%(path_to_files,path_to_files,path_to_files,path_to_files,path_to_files),shell=True)
subprocess.call('''rm %s*_coding*.bed %s*noncoding*.bed %s*5utr_start.bed'''%(path_to_files,path_to_files,path_to_files),shell=True)

#find gene_body.bed
genes = pd.read_csv(str(input_gtf_file), header=None, sep="\t",dtype = {0 :str})
#input_gene=pd.read_csv(str(input_gtf_file),sep='\t',header=None)

#genes=genes[genes[0].str.len()<=5]
genes.columns=['chr','unknow', 'exon', 'g_str', 'g_end', 'g_score', 'g_dir','.', 'gene_name']
genes=genes[genes.exon=='exon']

genes.gene_name=' '+genes.gene_name

gene_col=genes['gene_name'].str.split(';', expand=True) 
#header=row1
gene_col.columns=gene_col.iloc[1,:]
#filter gene_col contains 'gene_id'
gene_id = gene_col.filter(regex='gene_id')
gene_id = gene_id.iloc[:,0].str.split(' ', expand=True)
#remove head and end's character
gene_id[2] = gene_id[2].map(lambda x: x.lstrip('"').rstrip('"'))
gene_id.columns=['num','g_name','gene_id']
gene_bed = genes.loc[:,['chr', 'g_str', 'g_end', 'g_score','g_dir']].join(gene_id['gene_id'])
#change order
gene_bed = gene_bed.loc[:,['chr', 'g_str', 'g_end', 'gene_id','g_score','g_dir']]
#keep only one g_str, g_end site, ex:
gene_bed=gene_bed.drop_duplicates(subset=['g_str','g_end'],keep='first')
gene_bed=gene_bed.sort_values(['chr','g_str'],ascending=[True,True])
gene_bed=gene_bed.drop_duplicates(subset=['g_str'],keep='last')
# combine gene exons, and keep mininum g_str and maximum g_end
gene_group=gene_bed.groupby(['chr','gene_id','g_score','g_dir']).agg({'g_str':'min', 'g_end':'max'}).reset_index()
gene_group = gene_group.drop_duplicates(subset=['g_str','g_end'],keep='first')
gene_body=gene_group.sort_values(['chr','g_str'],ascending=[True,True])
gene_body=gene_body.drop_duplicates(subset=['g_str'],keep='last')
#redo the order of columns
gene_body = gene_body.loc[:,['chr', 'g_str', 'g_end', 'gene_id','g_score','g_dir']]
gene_body=gene_body[~gene_body.gene_id.str.contains('MI')]
gene_body.to_csv(input_gtf_file+'_Genebody_bed6.bed', sep='\t',index=False, header=None)

#find promoter.bed
gene_body['pro_str'] = np.where(gene_body.g_dir == '+', gene_body.g_str - promoter_size, gene_body.g_end - 0)
gene_body['pro_end'] = np.where(gene_body.g_dir == '+', gene_body.g_str + 0, gene_body.g_end + promoter_size)
num = gene_body._get_numeric_data()
num[num<0]=0
gene_promoter = gene_body.loc[:, ['chr','pro_str','pro_end','gene_id','g_score','g_dir']]
gene_promoter.columns=['chr','g_str','g_end','gene_id','g_score','g_dir']

gene_promoter.to_csv(input_gtf_file+"_Promoter_bed6.bed", sep='\t',index=False, header=None)


#find igr.bed
gene_body['igr_str'] = gene_body['g_end'].shift(1).fillna(0).astype(int)+1
gene_body['igr_end'] = gene_body['g_str']-1
gene_body['igr_chr'] = gene_body['chr'].shift(1).fillna('chr1')
igrcol = gene_body.loc[:,['igr_chr','g_str','igr_str','igr_end']]
igrcol.columns = ['chr', 'g_str', 'igr_str', 'igr_end']
genecol=gene_body.loc[:,['chr','g_str','g_end','gene_id','g_score','g_dir']]
geneigr = pd.merge(genecol, igrcol, how='left', on=['chr', 'g_str'])
geneigr.igr_str=geneigr.igr_str.fillna(0).astype(int)
geneigr.igr_end=geneigr.igr_end.fillna(geneigr.g_str-1).astype(int)
geneigr = geneigr.loc[:,['chr', 'igr_str','igr_end', 'gene_id','g_score', 'g_dir']]
geneigr = geneigr[geneigr['igr_str']<geneigr['igr_end']]
geneigr=geneigr.drop_duplicates(subset=['igr_str','igr_end'],keep='first')
geneigr.to_csv(input_gene_name+"_IGR_bed6.bed", sep="\t", index=False, header=None)


for i in annotation_name:
        subprocess.call('''bedtools sort -i %s|bedtools merge -c 4,5,6 -o collapse,collapse,collapse >%s '''%(input_gene_name+'_'+i+'_bed6.bed',input_gene_name +'_'+i+'_merge.bed'),shell=True)

subprocess.call('''rm %s*3utr_bed6.bed %s*cds_bed6.bed %s*5utr_bed6.bed %s*exons_bed6.bed %s*gene_igr_bed6.bed %s*introns_bed6.bed '''%(path_to_files,path_to_files,path_to_files,path_to_files,path_to_files,path_to_files),shell=True)



#DMR enrichmet cal & plot
if(command=='Fold_Enrichment' or command=='all'):
    print ("*--------------------------------*")
    print ("|Applying DMR enrichment analysis|")
    print ("*--------------------------------*")

    Enrichment(context, dmr_cut)
#    Enrichment('CG', dmr_CG_cut)
#    Enrichment('CHG', dmr_CHG_cut)
#    Enrichment('CHH', dmr_CHH_cut)

else:
    pass

#DMG
def dmg(tag,dmrfile,direction,cutoff):
    subprocess.call('''bedtools intersect -a %s -b %s -wo >%s '''%(dmrfile,input_gene_name +'_Genebody_bed6.bed',path_to_files+"DMG_"+tag+"_"+direction+"_"+str(cutoff)+"_Genebody_list.txt"),shell=True)
    subprocess.call('''bedtools intersect -a %s -b %s -wo >%s '''%(dmrfile,input_gene_name +'_Promoter_bed6.bed',path_to_files+"DMG_"+tag+"_"+direction+"_"+str(cutoff)+"_Promoter_list.txt"),shell=True)
    
 


if(command=='DMG' or command=='all' ):

    print ("*---------------*")
    print ("|Identifying DMG|")
    print ("*---------------*")
    if context == 'CG':
        bed_form(path_to_files+"DMR_CG_all_"+str(dmr_cut)+".txt")
        bed_form(path_to_files+"DMR_CG_hyper_"+str(dmr_cut)+".txt")
        bed_form(path_to_files+"DMR_CG_hypo_"+str(dmr_cut)+".txt")
        # dmg('CG',"DMR_CG_all_"+str(dmr_CG_cut)+".txt.bed", dmr_CG_cut)
        dmg('CG',path_to_files+"DMR_CG_hyper_"+str(dmr_cut)+".txt.bed",'hyper', dmr_cut)
        dmg('CG',path_to_files+"DMR_CG_hypo_"+str(dmr_cut)+".txt.bed", 'hypo',dmr_cut)
        DMR_DMGPlot('CG', dmr_cut )
    elif context == 'CHG':
        bed_form(path_to_files+"DMR_CHG_all_"+str(dmr_cut)+".txt")
        bed_form(path_to_files+"DMR_CHG_hyper_"+str(dmr_cut)+".txt")
        bed_form(path_to_files+"DMR_CHG_hypo_"+str(dmr_cut)+".txt")
        dmg('CHG',path_to_files+"DMR_CHG_hyper_"+str(dmr_cut)+".txt.bed",'hyper', dmr_cut)
        dmg('CHG',path_to_files+"DMR_CHG_hypo_"+str(dmr_cut)+".txt.bed", 'hypo',dmr_cut)
        DMR_DMGPlot('CHG',dmr_cut)
    elif context == 'CHH':
        bed_form(path_to_files+"DMR_CHH_all_"+str(dmr_cut)+".txt")
        bed_form(path_to_files+"DMR_CHH_hyper_"+str(dmr_cut)+".txt")
        bed_form(path_to_files+"DMR_CHH_hypo_"+str(dmr_cut)+".txt")
        dmg('CHH',path_to_files+"DMR_CHH_hyper_"+str(dmr_CHH_cut)+".txt.bed",'hyper', dmr_cut)
        dmg('CHH',path_to_files+"DMR_CHH_hypo_"+str(dmr_CHH_cut)+".txt.bed", 'hypo',dmr_cut)
        DMR_DMGPlot('CHH',dmr_cut)

else:
    pass

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#####chrview
if(command=='ChrView' or command=='all'):

    chrview_exp= DMR_exp
    chrview_ctrl= DMR_ctrl


    print ("*--------------------------*")
    print ("|Generating ChrView figures|")
    print ("*--------------------------*")

    #combined = pd.DataFrame()
    #samples = pd.read_csv("samples_list.txt",header=None,sep="\t")
    #depth = 4
    #binSize = 100000000
    for sample in samples.itertuples():
        #print("Now reading " + sample[2])
        CGmap = pd.read_csv(path_to_files+sample[2], compression='gzip',header=None,sep="\t",dtype = {0 :str})
        chrs = CGmap[0].unique()
        Position = 0
        #print ("Position\tchromosome\tsRange\tmeanCG\tmeanCHG\tmeanCHH")
        with open(path_to_files+"chrView_list.txt",'a') as f:
            f.write(sample[1]+'\t'+sample[1]+"_"+str(binSize)+"_chrView.txt"+'\n')

        #savefile=pd.DataFrame(columns=["Position","chromosome","sRange","meanCG","meanCHG","meanCHH"])
        savefile=pd.DataFrame(columns=["Position","chromosome","start","end","meanCG","meanCHG","meanCHH"])
        #mylist=[]
        for chromosome in chrs:
            subset = CGmap[(CGmap[7] >= depth) & (CGmap[0] == chromosome) ]
            if subset.shape[0] == 0: continue 
            subset = subset._convert(numeric=True)
            maxPos = subset[2].max()
            bins = range(0,maxPos,binSize)
            groups = subset.groupby(pd.cut(subset[2], bins))



            for sRange,sValues in groups:
                Position += 1
                start = sRange.left
                end = sRange.right
                meanCG = sValues[sValues[3] == 'CG'][5].mean(skipna=True)
                meanCHG = sValues[sValues[3] == 'CHG'][5].mean(skipna=True)
                meanCHH = sValues[sValues[3] == 'CHH'][5].mean(skipna=True)
                out = map(str,[Position,chromosome,start,end,meanCG,meanCHG,meanCHH])

                out = pd.Series([Position,chromosome,start,end,meanCG,meanCHG,meanCHH],index=["Position","chromosome","start","end","meanCG","meanCHG","meanCHH"])
                savefile=savefile.append(out,ignore_index=True)
        savefile.to_csv(path_to_files+sample[1]+"_tmpchrView.txt",sep='\t',index=None)
    #sort chromosome
    chrlist=pd.read_csv(path_to_files+"chrView_list.txt",header=None,sep='\t')
    chrlist[2] = chrlist[0].map(samples.set_index(0)[2])
    chrlist.to_csv(path_to_files+"chrView_list.txt", sep ='\t',index=None, header=None)

    for sample in chrlist[0]:
        subprocess.call(''' (head -n 1 %s && tail -n +2 %s |sort -k2,2 -V -s) > %s '''%(path_to_files+sample+'_tmpchrView.txt',path_to_files+sample+'_tmpchrView.txt',path_to_files+sample+"_"+str(binSize)+'_chrView.txt'), shell = True)

    subprocess.call('''rm %s*_tmpchrView.txt'''%(path_to_files),shell=True)

    #plotting
    subprocess.call("Rscript --slave /MethylC-analyzer/scripts/chrView.R %s"%(path_to_files) , shell=True)
    
############chrView_delta################################################
    chrlist=pd.read_csv(path_to_files +"chrView_list.txt",header=None,sep='\t')
    
    expData = chrlist[chrlist.iloc[:,2] == chrview_exp][1]
    ctrlData = chrlist[chrlist.iloc[:,2] == chrview_ctrl][1]
    # read chrview files in two group 
    exp_df = []
    for expfile in expData:
        exp_df.append(pd.read_csv(path_to_files+expfile ,sep='\t'))

    ctrl_df = []
    for ctrlfile in ctrlData:
        ctrl_df.append(pd.read_csv(path_to_files+ctrlfile ,sep='\t'))
    
    #sum all valuew in each row in 3 context
    mcgexp,mchgexp,mchhexp,mcgctrl,mchgctrl,mchhctrl=0,0,0,0,0,0
    for j in range(0,len(expData)):
        mcgexp+=exp_df[j].iloc[:,4]
        mchgexp+=exp_df[j].iloc[:,5]
        mchhexp+=exp_df[j].iloc[:,6]
    
    for k in range(0,len(ctrlData)):
        mcgctrl+=ctrl_df[k].iloc[:,4]
        mchgctrl+=ctrl_df[k].iloc[:,5]
        mchhctrl+=ctrl_df[k].iloc[:,6]


   # get average M level in each row(region)
    mcgexp/=len(expData)
    mchgexp/=len(expData)
    mchhexp/=len(expData)
    mcgctrl/=len(ctrlData)
    mchgctrl/=len(ctrlData)
    mchhctrl/=len(ctrlData)

    ## delta m level in 3 context
    mcg,mchg,mchh=0,0,0

    mcg=mcgexp-mcgctrl
    mchg=mchgexp-mchgctrl
    mchh=mchhexp-mchhctrl

    #generate delta report
    exp_df[0].iloc[:,4]=mcg
    exp_df[0].iloc[:,5]=mchg
    exp_df[0].iloc[:,6]=mchh
    report = exp_df[0]
    report.columns = ['Position', 'chromosome', 'start', 'end', 'deltaCG', 'deltaCHG','deltaCHH']
    report.to_csv(path_to_files+'chrView_delta.txt',sep='\t',index=False)

    pd.DataFrame([[path_to_files+chrview_exp+' - '+chrview_ctrl,'chrView_delta.txt']]).to_csv(path_to_files+'chrView_delta_list.txt',sep='\t',index=False,header=None)


    #plotting chrview difference
    subprocess.call("Rscript --slave /MethylC-analyzer/scripts/chrView_delta.R %s"%(path_to_files), shell=True)
else:
    pass



#metaplot

if(command=='Metaplot' or command=='all'):
    metaplot_exp=DMR_exp
    metaplot_ctrl=DMR_ctrl
    
    metaplot_gene_feature=str('Genebody')

    print ("*--------------------*")
    print ("|Generating metaplots|")
    print ("*--------------------*")

    #union to bw
    #cgmap=pd.read_csv("union4.txt",sep='\t',dtype={0:str})
    cgmap=union
    genebodybed=input_gene_name+"_"+metaplot_gene_feature+"_merge.bed"
    bwheader=[]
    for chr in cgmap.iloc[:,0].unique():
        k=cgmap[cgmap.iloc[:,0]==chr]
        maxPos=k.iloc[:,1].max()
        bwheader.append((str(chr),int(maxPos)))

    #3 first sample
    for i in range(3,cgmap.shape[1]):
        #print i
        Samplename=cgmap.columns[i]
        cgoutfile=path_to_files+Samplename+"_CG.bw"
        chgoutfile=path_to_files+Samplename+"_CHG.bw"
        chhoutfile=path_to_files+Samplename+"_CHH.bw"

        bw_CG = pyBigWig.open(cgoutfile,"w")
        bw_CG.addHeader(bwheader)
        bw_CHG = pyBigWig.open(chgoutfile,"w")
        bw_CHG.addHeader(bwheader)
        bw_CHH = pyBigWig.open(chhoutfile,"w")
        bw_CHH.addHeader(bwheader)

        for chromosome in cgmap.iloc[:,0].unique():
            k=cgmap[cgmap.iloc[:,0]==chromosome]
            CG = k[(k.iloc[:,2] == 'CG') & (k.iloc[:,i] != np.nan)]
            CHG = k[(k.iloc[:,2] == 'CHG') & (k.iloc[:,i] != np.nan)]
            CHH = k[(k.iloc[:,2] == 'CHH') & (k.iloc[:,i] != np.nan)]

            bw_CG.addEntries(str(chromosome),[x-1 for x in CG.iloc[:,1].tolist()],values=CG.iloc[:,i].astype(float).tolist(), span=1)
            bw_CHG.addEntries(str(chromosome),[x-1 for x in CHG.iloc[:,1].tolist()],values=CHG.iloc[:,i].astype(float).tolist(), span=1)
            bw_CHH.addEntries(str(chromosome),[x-1 for x in CHH.iloc[:,1].tolist()],values=CHH.iloc[:,i].astype(float).tolist(), span=1)

        bw_CG.close()
        bw_CHG.close()
        bw_CHH.close()

        #catch warninig 
        devnull = open(os.devnull, 'w')
           
        subprocess.call('''computeMatrix scale-regions -S %s -R %s -a 2000 -b 2000 -m 4000 -bs 100 -out %s'''%(path_to_files+Samplename+"_CG.bw",genebodybed,path_to_files+Samplename+"_CG.matrix.gz"),shell=True,stdout=devnull, stderr=devnull)
        subprocess.call('''computeMatrix scale-regions -S %s -R %s -a 2000 -b 2000 -m 4000 -bs 100 -out %s'''%(path_to_files+Samplename+"_CHG.bw",genebodybed,path_to_files+Samplename+"_CHG.matrix.gz"),shell=True,stdout=devnull, stderr=devnull)
        subprocess.call('''computeMatrix scale-regions -S %s -R %s -a 2000 -b 2000 -m 4000 -bs 100 -out %s'''%(path_to_files+Samplename+"_CHH.bw",genebodybed,path_to_files+Samplename+"_CHH.matrix.gz"),shell=True,stdout=devnull, stderr=devnull)

    #metaplot

    subprocess.call("Rscript --slave /MethylC-analyzer/scripts/metaplot.R "+path_to_files+ ' '+metaplot_gene_feature,shell=True)

    ##generating delta files
    Delta_Meta('CG')
    Delta_Meta('CHG')
    Delta_Meta('CHH')
    #ploting delta meta
    subprocess.call("Rscript --slave /MethylC-analyzer/scripts/metaplot_delta.R " +metaplot_exp+' '+metaplot_ctrl + ' '+path_to_files,shell=True)

else:
    pass


#### output plotting log ###
fp = open("plot.log", "a")
print("#command used for plotting", file=fp)
fp.close()

if (command=='Heatmap_PCA' or command=='all'):
    fp = open(path_to_files +"plot.log", "a")
    print("#PCA & Heatmap:",file=fp)
    print ("Rscript --slave /MethylC-analyzer/scripts/heatmap_PCA_all.R "+ path_to_files +"CommonRegion_"+context+".txt "+ str(pca_heat_cut)+ " "+str(path_to_files), file=fp)
    fp.close()
    
if(command=='ChrView' or command=='all'):
    fp = open("plot.log", "a")
    print("#ChrView:",file=fp)
    print ("Rscript --slave /MethylC-analyzer/scripts/chrView.R", file=fp)
    print("#delta ChrView:",file=fp)
    print ("Rscript --slave /MethylC-analyzer/scripts/chrView_delta.R", file=fp)
    fp.close()

if(command=='Metaplot' or command=='all'):
    fp = open("plot.log", "a")
    print("#metaplot:",file=fp)
    print ("Rscript --slave /MethylC-analyzer/scripts/metaplot.R "+ str(metaplot_gene_feature), file=fp)
    print("#delta metaplot:",file=fp)
    print ("Rscript --slave /MethylC-analyzer/scripts/metaplot.R "+ str(metaplot_exp)+' '+str(metaplot_ctrl), file=fp)
    fp.close()
