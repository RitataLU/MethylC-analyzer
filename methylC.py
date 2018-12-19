import sys
import subprocess, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
from math import log
import scipy.stats as ss
import time
import argparse
import glob
import pyBigWig


parser = argparse.ArgumentParser()

parser.add_argument("-d",help="min site of #C+#T",dest='depth',default=4)
parser.add_argument("-r",help="size of region",dest='region',default=500)
parser.add_argument("-q",help="qualified site within a region",dest='qualified',default=4)
parser.add_argument("-hcgc",help="PCA & Heatmap_CG_cutoff",dest='heatmap_cg_cutoff',default=0.2)
parser.add_argument("-hchgc",help="PCA & Heatmap_CHG_cutoff",dest='heatmap_chg_cutoff',default=0.2)
parser.add_argument("-hchhc",help="PCA & Heatmap_CHH_cutoff",dest='heatmap_chh_cutoff',default=0.2)
parser.add_argument("-dcgc",help="DMR_CG_cutoff",dest='dmr_cg_cutoff',default=0.2)
parser.add_argument("-dchgc",help="DMR_CHG_cutoff",dest='dmr_chg_cutoff',default=0.2)
parser.add_argument("-dchhc",help="DMR_CHH_cutoff",dest='dmr_chh_cutoff',default=0.2)
parser.add_argument("-b",help="resolution of chrView and Metaplot",dest='bin_size',default=1000000)
parser.add_argument("-p",help="promoter_size",dest='promoter_size',default=2000)
parser.add_argument("input_gtf_file",help="path of gene annotation")
#parser.add_argument("sample_list",help="path of sample list")

args = parser.parse_args()

depth=int(args.depth)
region=int(args.region)
quaifiedSite=int(args.quaified)
pca_heat_cg_cut=float(args.heatmap_cg_cutoff)
pca_heat_chh_cut=float(args.heatmap_chh_cutoff)
pca_heat_chg_cut=float(args.heatmap_chg_cutoff)
dmr_cg_cut=float(args.dmr_cg_cutoff)
dmr_chh_cut=float(args.dmr_chh_cutoff)
dmr_chg_cut=float(args.dmr_chg_cutoff)
binSize=int(args.bin_size)
promoter=int(args.promoter_size)
input_gene=pd.read_csv(str(args.input_gtf_file),sep='\t',header=None)
input_gtf_file=str(args.input_gtf_file)

#Choose Tools
Heatmap_PCA=raw_input('Heatmap & PCA Analysis?  (y/n)')
Fold_Enrichment=raw_input('Use Fold Enrichment Analysis?  (y/n)')
DMG=raw_input('DMG Analysis?  (y/n)')
ChrView=raw_input('Chromosome View Analysis?  (y/n)')
Metaplot=raw_input('Metaplot Analysis?  (y/n)')

if(ChrView=='y'):
    chrview_lo=str(raw_input('enter experimental group of chrview'))
    chrview_ro=str(raw_input('enter control group of chrview'))
else:
    pass

if(Metaplot=='y'):
    metaplot_lo=str(raw_input('enter experimental group of metaplot'))
    metaplot_ro=str(raw_input('enter control group of metaplot'))
    metaplot_gene_feature=str(raw_input('enter one genomic feature of metaplot'))
else:
    pass


context=["CG","CHG","CHH"]
#unionsite

combined = pd.DataFrame()
samples = pd.read_csv("samples_list.txt",header=None,sep="\t")

for sample in samples.itertuples():
    print("Now processing " + sample[2])
    CGmap = pd.read_csv(sample[2], header=None,sep="\t",dtype =
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

combined[['chr']] = combined[['chr']].astype(int)

combined = combined.sort_values(['chr','pos'],ascending=[True,True])
combined.to_csv('Unionsite.txt',sep = '\t',na_rep='-',index=False)


#common region

combined=pd.read_csv('Unionsite.txt',sep='\t',na_values='-')
union=combined
#1.unionsite --> combined
chrs = union['chr'].unique()


#context one by one
for cxt in context:
	outfile = 'CommonRegion_' + cxt + '.txt'
	with open(outfile,'w') as of:
		of.write('\t'.join(['chr','start','end'] + union.columns.values.tolist()[3:]) + "\n")
		for chromosome in chrs:
			subset = union[(union['context'] == cxt) & (union['chr'] == chromosome)]
			maxPos = subset['pos'].max()
			bins = range(0,maxPos,500)
			groups = subset.groupby(pd.cut(subset['pos'], bins))
			for sRange,sValues in groups:
			    #rows, columns = sValues.shape
			    #with a region find out minimun sites, less than depthcutoff
				minDepth = sValues.iloc[:,3:].count().min()
				if minDepth >= quaifiedSite:
					methMeth = sValues.iloc[:,3:].astype(float).mean().tolist()
					methMeth2 = [("%.3f" %x) for x in methMeth]
					start = sRange.left
					end = sRange.right
					out = map(str,[chromosome,start,end] + methMeth2)
					outStr = '\t'.join(out)	
					of.write(outStr + "\n")				

				else:
					pass

					
#DMR
filename1=samples.iloc[0,2] #read group name via samples (read in upper step)
filename2_temp=samples
filename2=filename2_temp.iloc[len(filename2_temp.index)-1,2]

file1=pd.read_csv("CommonRegion_CG.txt",sep="\t",header=None)
file2=pd.read_csv("CommonRegion_CHG.txt",sep="\t",header=None)
file3=pd.read_csv("CommonRegion_CHH.txt",sep="\t",header=None)

CG_data=np.array(file1)
CHH_data=np.array(file2)
CHG_data=np.array(file3)

cutoff_CG=dmr_cg_cut
cutoff_CHH=dmr_chh_cut
cutoff_CHG=dmr_chg_cut

index1_CG=file1.iloc[0,:].str.contains(filename1,regex=False)
index2_CG=file1.iloc[0,:].str.contains(filename2,regex=False)
index1_CG=np.where(np.array(index1_CG)==True)[0]
index2_CG=np.where(np.array(index2_CG)==True)[0]

index1_CHH=file2.iloc[0,:].str.contains(filename1,regex=False)
index2_CHH=file2.iloc[0,:].str.contains(filename2,regex=False)
index1_CHH=np.where(np.array(index1_CHH)==True)[0]
index2_CHH=np.where(np.array(index2_CHH)==True)[0]

index1_CHG=file3.iloc[0,:].str.contains(filename1,regex=False)
index2_CHG=file3.iloc[0,:].str.contains(filename2,regex=False)
index1_CHG=np.where(np.array(index1_CHG)==True)[0]
index2_CHG=np.where(np.array(index2_CHG)==True)[0]


CG_delta=np.array(0.0)
CG_delta=np.append(dmr_cg_cut,CG_data[1:,index2_CG].astype(np.float).mean(axis=1)-CG_data[1:,index1_CG].astype(np.float).mean(axis=1)) #diffrence of average

CHH_delta=np.array(0.0)
CHH_delta=np.append(dmr_chh_cut,CHH_data[1:,index2_CHH].astype(np.float).mean(axis=1)-CHH_data[1:,index1_CHH].astype(np.float).mean(axis=1)) #diffrence of average

CHG_delta=np.array(0.0)
CHG_delta=np.append(dmr_chg_cut,CHG_data[1:,index2_CHG].astype(np.float).mean(axis=1)-CHG_data[1:,index1_CHG].astype(np.float).mean(axis=1)) #diffrence of average


CG_data[1:,index1_CG[0]]=CG_data[1:,index1_CG[0]].astype(np.float)+0.0001 #t-test avoid zero condition
CG_data[1:,index2_CG[0]]=CG_data[1:,index2_CG[0]].astype(np.float)+0.0001 #t-test avoid zero condition

CHH_data[1:,index1_CHH[0]]=CHH_data[1:,index1_CHH[0]].astype(np.float)+0.0001 #t-test avoid zero condition
CHH_data[1:,index2_CHH[0]]=CHH_data[1:,index2_CHH[0]].astype(np.float)+0.0001 #t-test avoid zero condition

CHG_data[1:,index1_CHG[0]]=CHG_data[1:,index1_CHG[0]].astype(np.float)+0.0001 #t-test avoid zero condition
CHG_data[1:,index2_CHG[0]]=CHG_data[1:,index2_CHG[0]].astype(np.float)+0.0001 #t-test avoid zero condition



df_a=pd.DataFrame(CG_data[1:,index1_CG]) #preproccess of merge dataframe (convert numpy_array to pandas_dataframe)
df_b=pd.DataFrame(CG_data[1:,index2_CG]) #preproccess of merge dataframe (convert numpy_array to pandas_dataframe)

df_c=pd.DataFrame(CHH_data[1:,index1_CHH]) #preproccess of merge dataframe (convert numpy_array to pandas_dataframe)
df_d=pd.DataFrame(CHH_data[1:,index2_CHH]) #preproccess of merge dataframe (convert numpy_array to pandas_dataframe)

df_e=pd.DataFrame(CHG_data[1:,index1_CHG]) #preproccess of merge dataframe (convert numpy_array to pandas_dataframe)
df_f=pd.DataFrame(CHG_data[1:,index2_CHG]) #preproccess of merge dataframe (convert numpy_array to pandas_dataframe)



idx1 = df_a.index.intersection(df_b.index) #establish index to use in ttest all element together
CG_ttest=np.array(1.1)
CG_ttest=np.append(CG_ttest,ss.stats.ttest_ind(df_a.loc[idx1].astype(np.float), df_b.loc[idx1].astype(np.float), axis=1)[1]) #t-test

idx2 = df_c.index.intersection(df_d.index) #establish index to use in ttest all element together
CHH_ttest=np.array(1.1)
CHH_ttest=np.append(CHH_ttest,ss.stats.ttest_ind(df_c.loc[idx2].astype(np.float), df_d.loc[idx2].astype(np.float), axis=1)[1]) #t-test

idx3 = df_e.index.intersection(df_f.index) #establish index to use in ttest all element together
CHG_ttest=np.array(1.1)
CHG_ttest=np.append(CHG_ttest,ss.stats.ttest_ind(df_e.loc[idx3].astype(np.float), df_f.loc[idx3].astype(np.float), axis=1)[1]) #t-test



CG_data[1:,index1_CG[0]]-=0.0001 #t-test avoid zero condition recovery
CG_data[1:,index2_CG[0]]-=0.0001 #t-test avoid zero condition recovery

CHH_data[1:,index1_CHH[0]]-=0.0001 #t-test avoid zero condition recovery
CHH_data[1:,index2_CHH[0]]-=0.0001 #t-test avoid zero condition recovery

CHG_data[1:,index1_CHG[0]]-=0.0001 #t-test avoid zero condition recovery
CHG_data[1:,index2_CHG[0]]-=0.0001 #t-test avoid zero condition recovery



CGpd=pd.DataFrame(pd.concat([pd.DataFrame(CG_data),pd.DataFrame(CG_delta),pd.DataFrame(CG_ttest)],axis=1,ignore_index=True)) #conbine dataframe
CGpd.columns=pd.concat([file1.iloc[0],pd.DataFrame(['delta','pvalue'])])[0] #fix header name
CGpd=CGpd.drop([0]) #drop 1st row (header name)

CHHpd=pd.DataFrame(pd.concat([pd.DataFrame(CHH_data),pd.DataFrame(CHH_delta),pd.DataFrame(CHH_ttest)],axis=1,ignore_index=True)) #conbine dataframe
CHHpd.columns=pd.concat([file2.iloc[0],pd.DataFrame(['delta','pvalue'])])[0] #fix header name
CHHpd=CHHpd.drop([0]) #drop 1st row (header name)

CHGpd=pd.DataFrame(pd.concat([pd.DataFrame(CHG_data),pd.DataFrame(CHG_delta),pd.DataFrame(CHG_ttest)],axis=1,ignore_index=True)) #conbine dataframe
CHGpd.columns=pd.concat([file3.iloc[0],pd.DataFrame(['delta','pvalue'])])[0] #fix header name
CHGpd=CHGpd.drop([0]) #drop 1st row (header name)


CG_DMRpd=CGpd.loc[np.abs(CGpd['delta'])>=cutoff_CG].loc[CGpd['pvalue']<0.05]

CHH_DMRpd=CHHpd.loc[np.abs(CHHpd['delta'])>=cutoff_CHH].loc[CHHpd['pvalue']<0.05]

CHG_DMRpd=CHGpd.loc[np.abs(CHGpd['delta'])>=cutoff_CHG].loc[CHGpd['pvalue']<0.05]


CG_DMRpd[CG_DMRpd['delta']>0].to_csv('DMR_CG_hyper.txt',sep='\t',index=None)
CG_DMRpd[CG_DMRpd['delta']<0].to_csv('DMR_CG_hypo.txt',sep='\t',index=None)
CG_DMRpd.to_csv('DMR_CG_all.txt',sep='\t',index=None)

CHH_DMRpd[CHH_DMRpd['delta']>0].to_csv('DMR_CHH_hyper.txt',sep='\t',index=None)
CHH_DMRpd[CHH_DMRpd['delta']<0].to_csv('DMR_CHH_hypo.txt',sep='\t',index=None)
CHH_DMRpd.to_csv('DMR_CHH_all.txt',sep='\t',index=None)

CHG_DMRpd[CHG_DMRpd['delta']>0].to_csv('DMR_CHG_hyper.txt',sep='\t',index=None)
CHG_DMRpd[CHG_DMRpd['delta']<0].to_csv('DMR_CHG_hypo.txt',sep='\t',index=None)
CHG_DMRpd.to_csv('DMR_CHG_all.txt',sep='\t',index=None)


#heatmap_PCA
if(Heatmap_PCA=='y'):
    subprocess.call('''Rscript --slave heatmap_PCA_all.R %s %s'''%("CommonRegion_CG.txt",pca_heat_cg_cut), shell=True)
    subprocess.call('''Rscript --slave heatmap_PCA_all.R %s %s'''%("CommonRegion_CHG.txt",pca_heat_chh_cut), shell=True)
    subprocess.call('''Rscript --slave heatmap_PCA_all.R %s %s'''%("CommonRegion_CHH.txt",pca_heat_chg_cut), shell=True)
else:
    pass



#transfer txt to bed
def bed_form(inputxt):
	subprocess.call('''awk '{if (NR!=1) print $1"\t"$2"\t"$3}' %s > %s'''%(inputxt, inputxt+".bed"), shell=True)
	return inputxt+".bed"

#input_gene=input_gene
#input_gene='arabidopsis.gtf'
#annotation name
gene = "gene_body"
exon = "exons"
intron = "introns"
utr3 = "3utr"
utr5 = "5utr"
cds = "cds"
promoter = "gene_promoter"
igr = "gene_igr"
annotation_name=[promoter,gene,exon,intron,utr5,cds,utr3,igr]

#subprocess.call('''gffread %s -T -o %s'''%(input_gene,input_gene+'.gtf'),shell=True)

#enrichment
print ("*--------------------------------*")
print ("|Extract UTR, exon, cds from gene|")
print ("*--------------------------------*")

subprocess.call('''python ./extract_transcript_regions.py -i %s -o %s --gtf'''%(input_gtf_file,input_gtf_file), shell=True)

print ("*-------------------------------------*")
print ("|Convert this blockbed (bed12) to bed6|")
print ("*-------------------------------------*")

for i in annotation_name:
    subprocess.call('''cat %s | bed12ToBed6 -i stdin -n > %s'''%(input_gtf_file+'_'+i+'.bed',input_gtf_file+'_'+i+'_bed6.bed'),shell=True)
    print (i)

#find gene_body.bed
genes = pd.read_csv(input_gene, header=None, sep="\t",dtype = {0 :str})
#genes=genes[genes[0].str.len()<=5]
genes.columns=['chr','unknow', 'exon', 'g_str', 'g_end', 'g_score', 'g_dir','.', 'gene_name']
genes=genes[genes.exon=='exon']

genes.gene_name=' '+genes.gene_name

gene_col=genes['gene_name'].str.split(';', expand=True)
gene_col.columns=gene_col.ix[1,:]
gene_id = gene_col.filter(regex='gene_id')
gene_id = gene_id.ix[:,0].str.split(' ', expand=True)
gene_id[2] = gene_id[2].map(lambda x: x.lstrip('"').rstrip('"'))
gene_id.columns=['num','g_name','gene_id']
gene_bed = genes.ix[:,['chr', 'g_str', 'g_end', 'g_score','g_dir']].join(gene_id['gene_id'])
gene_bed = gene_bed.ix[:,['chr', 'g_str', 'g_end', 'gene_id','g_score','g_dir']]
gene_bed=gene_bed.drop_duplicates(subset=['g_str','g_end'],keep='first')
gene_bed=gene_bed.sort_values(['chr','g_str'],ascending=[True,True])
gene_bed=gene_bed.drop_duplicates(subset=['g_str'],keep='last')
gene_group=gene_bed.groupby(['chr','gene_id','g_score','g_dir']).agg({'g_str':'min', 'g_end':'max'}).reset_index()
gene_group = gene_group.drop_duplicates(subset=['g_str','g_end'],keep='first')
gene_body=gene_group.sort_values(['chr','g_str'],ascending=[True,True])
gene_body=gene_body.drop_duplicates(subset=['g_str'],keep='last')
gene_body = gene_body.ix[:,['chr', 'g_str', 'g_end', 'gene_id','g_score','g_dir']]
#gene_body=gene_body[~gene_body.gene_id.str.contains('MI')]
gene_body.to_csv(input_gene_name+'_gene_body_bed6.bed', sep='\t',index=False, header=None)

#find promoter.bed
gene_body['pro_str'] = np.where(gene_body.g_dir == '+', gene_body.g_str - 2000, gene_body.g_end - 0)
gene_body['pro_end'] = np.where(gene_body.g_dir == '+', gene_body.g_str + 0, gene_body.g_end + 2000)
num = gene_body._get_numeric_data()
num[num<0]=0
gene_promoter = gene_body.ix[:, ['chr','pro_str','pro_end','gene_id','g_score','g_dir']]
gene_promoter.columns=['chr','g_str','g_end','gene_id','g_score','g_dir']

gene_promoter.to_csv(input_gene_name+"_gene_promoter_bed6.bed", sep='\t',index=False, header=None)


#find igr.bed
gene_body['igr_str'] = gene_body['g_end'].shift(1).fillna(0).astype(int)+1
gene_body['igr_end'] = gene_body['g_str']-1
gene_body['igr_chr'] = gene_body['chr'].shift(1).fillna('chr1')
igrcol = gene_body.ix[:,('igr_chr','g_str','igr_str','igr_end')]
igrcol.columns = ['chr', 'g_str', 'igr_str', 'igr_end']
genecol=gene_body.ix[:,['chr','g_str','g_end','gene_id','g_score','g_dir']]
geneigr = pd.merge(genecol, igrcol, how='left', on=['chr', 'g_str'])
geneigr.igr_str=geneigr.igr_str.fillna(0).astype(int)
geneigr.igr_end=geneigr.igr_end.fillna(geneigr.g_str-1).astype(int)
geneigr = geneigr.ix[:,('chr', 'igr_str','igr_end', 'gene_id','g_score', 'g_dir')]
geneigr = geneigr[geneigr['igr_str']<geneigr['igr_end']]
geneigr=geneigr.drop_duplicates(subset=['igr_str','igr_end'],keep='first')
geneigr.to_csv(input_gene_name+"_gene_igr_bed6.bed", sep="\t", index=False, header=None)


for i in annotation_name:
	subprocess.call('''bedtools sort -i %s|bedtools merge -c 4,5,6 -o collapse,collapse,collapse >%s '''%(input_gene_name+'_'+i+'_bed6.bed',input_gene_name +'_'+i+'_merge.bed'),shell=True)

#input_gene_exon_merge.bed



#enrichment
#1. overlap with dmr
def overlap(bed1,bed2):
	overlap=subprocess.check_output("bedtools intersect -a %s -b %s |awk '{size+=$3 - $2}END{print size}' "%(bed1,bed2), shell=True)
	if len(overlap)<2:
		return 1
	else:
		return int(overlap)
	

#DMR enrichmet cal & plot
if(Fold_Enrichment=='y'):
    for tag in context:
        savefile=pd.DataFrame(columns=["dmr","feature","overlap_size","dmr_size","feature_size","genome_size","enrichment"])
        for i in annotation_name:
            #dmr="DMR_"+tag+"_"+str(dmr_cut)+"_all.txt"
            dmr="DMR_"+tag+"_all.txt"
            feature=str(i)
            overlap_size=overlap(input_gene_name +'_'+i+'_merge.bed',bed_form(dmr))
            dmrlen=pd.read_csv(dmr,sep='\t')
            if len(dmrlen)<=1:
                dmr_size=1000000000000
            else:
                dmr_size=int(subprocess.check_output('''awk '{size+=$3-$2}END{print size}' %s'''%(bed_form(dmr)), shell=True))

            feature_size=overlap(bed_form("CommonRegion_"+tag+".txt"),bed_form(input_gene_name +'_'+i+'_merge.bed'))
            genome_size=int(subprocess.check_output('''awk '{size+=$3-$2}END{print size}' %s'''%(bed_form("CommonRegion_"+tag+".txt")), shell=True))*1.0

            enrichment=math.log(((float(overlap_size)/float(dmr_size))/(float(feature_size)/float(genome_size))),2)
            out = [dmr,feature,overlap_size,dmr_size,feature_size,genome_size,enrichment]
            savefile.loc[len(savefile)]=out

            savefile.to_csv("enrichment_"+tag+".txt",sep='\t',index=None)
            #plot enrichment
            #print "\n""Making enrichment plot""\n"
            plt.style.use('ggplot')
            annotationname = ['Promoter','Genebody','Exon','Intron','5UTR','CDS','3UTR','IGR']
            annotationname_index = range(len(annotationname))
        Eplt=[savefile['enrichment']]

        for value in Eplt:
            fig = plt.figure()
            ax1 = fig.add_subplot(1,1,1)
            ax1.bar(annotationname_index, value,align='center',color='darkblue')
            ax1.xaxis.set_ticks_position('bottom')
            ax1.yaxis.set_ticks_position('left')
            plt.xticks(annotationname_index, annotationname, rotation=90,fontsize='small')
            plt.ylabel("Fold Enrichment (log2)")
            plt.title(tag+'_Fold_Enrichment')
            plt.savefig(tag+'_Fold_Enrichment'+'.png',dpi=400,bbox_inches='tight')
            plt.close(fig)
else:
    pass
        

#DMG
def dmg(tag):
    #subprocess.call('''bedtools intersect -a %s -b %s >%s '''%(input_gene_name +'_gene_body_merge.bed',"DMR_"+tag+"_"+str(dmr_cut)+"_all.txt.bed","DMG_"+tag+"_gene_list.txt"),shell=True)
    #subprocess.call('''bedtools intersect -a %s -b %s >%s '''%(input_gene_name +'_gene_promoter_merge.bed',"DMR_"+tag+"_"+str(dmr_cut)+"_all.txt.bed","DMG_"+tag+"_promoter_list.txt"),shell=True)
    subprocess.call('''bedtools intersect -a %s -b %s >%s '''%(input_gene_name +'_gene_body_bed6.bed',"DMR_"+tag+"_all.txt.bed","DMG_"+tag+"_gene_list.txt"),shell=True)
    subprocess.call('''bedtools intersect -a %s -b %s >%s '''%(input_gene_name +'_gene_promoter_bed6.bed',"DMR_"+tag+"_all.txt.bed","DMG_"+tag+"_promoter_list.txt"),shell=True)

if(DMG=='y'):
    for k in context:
        dmg(k)
else:
    pass
        

#----------------------------------------------------------------------------------------------------------------------------------------------------------
#chrView
if(ChrView=='y'):
    #combined = pd.DataFrame()
    #samples = pd.read_csv("samples_list.txt",header=None,sep="\t")
    #minCoverage = 4
    #binSize = 1000000

    for sample in samples.itertuples():
        #print("Now reading " + sample[2])
        CGmap = pd.read_csv(sample[2], compression='gzip',header=None,sep="\t",dtype = {0 :str})
        chrs = CGmap[0].unique()
        Position = 0
        #print ("Position\tchromosome\tsRange\tmeanCG\tmeanCHG\tmeanCHH")
        with open("chrView_list.txt",'a') as f:
                f.write(sample[1]+'\t'+sample[1]+"_chrView.txt"+'\n')


        #savefile=pd.DataFrame(columns=["Position","chromosome","sRange","meanCG","meanCHG","meanCHH"])
        savefile=pd.DataFrame(columns=["Position","chromosome","start","end","meanCG","meanCHG","meanCHH"])
        #mylist=[]		
        for chromosome in chrs:
            subset = CGmap[(CGmap[7] >= depth) & (CGmap[0] == chromosome) ]
            maxPos = subset[2].max()
            bins = range(0,maxPos,binSize)
            groups = subset.groupby(pd.cut(subset[2], bins))

            #make list for R draw	
            #with open("chrView_list.txt",'a') as f:
                #f.write(sample[1]+'\t'+sample[1]+"_"+str(binSize)+"_chrView.txt"+'\n':)

            for sRange,sValues in groups:
                Position += 1
                start = sRange.left
                end = sRange.right
                meanCG = sValues[sValues[3] == 'CG'][5].mean()
                meanCHG = sValues[sValues[3] == 'CHG'][5].mean()
                meanCHH = sValues[sValues[3] == 'CHH'][5].mean()
                out = map(str,[Position,chromosome,start,end,meanCG,meanCHG,meanCHH])

                #out = map(str,[Position,chromosome,sRange,meanCG,meanCHG,meanCHH])
                #out = pd.Series([Position,chromosome,sRange,meanCG,meanCHG,meanCHH],index=["Position","chromosome","start","end,","meanCG","meanCHG","meanCHH"])
                out = pd.Series([Position,chromosome,start,end,meanCG,meanCHG,meanCHH],index=["Position","chromosome","start","end","meanCG","meanCHG","meanCHH"])

                savefile=savefile.append(out,ignore_index=True)
                #print(sample[1]+"_chrView_2.txt")	

        savefile.to_csv(sample[1]+"_chrView.txt",sep='\t',index=None)
    subprocess.call("Rscript --slave chrView.R" , shell=True)

    #chrView_mean
    chrlist=pd.read_csv("chrView_list.txt",header=None,sep='\t')

    filename1=samples.iloc[0,2] #read group name via samples (read in upper step)
    filename2_temp=samples
    filename2=filename2_temp.iloc[len(filename2_temp.index)-1,2]

    #check split index
    si=0
    for i in chrlist.iloc[:,1]:
        if filename1 in i:
            si+=1
        else:
            break

    #load data
    data=[]
    for i in chrlist.iloc[:,1]:
        data.append(pd.read_csv(i,sep='\t'))

    mcg1,mchg1,mchh1,mcg2,mchg2,mchh2=0,0,0,0,0,0
    for i in range(0,len(chrlist)):
        if(i<si):
            mcg1+=data[i].iloc[:,4]
            mchh1+=data[i].iloc[:,5]
            mchg1+=data[i].iloc[:,6]
        elif(i>=si):
            mcg2+=data[i].iloc[:,4]
            mchh2+=data[i].iloc[:,5]
            mchg2+=data[i].iloc[:,6]

    mcg1/=si
    mchg1/=si
    mchh1/=si
    mcg2/=len(chrlist)-si
    mchg2/=len(chrlist)-si
    mchh2/=len(chrlist)-si

    mcg,mchg,mchh=0,0,0
    if(chrview_lo==filename1 and chrview_ro==filename2):
        mcg=mcg1-mcg2
        mchg=mchg1-mchg2
        mchh=mchh1-mchh2
    elif(chrview_lo==filename2 and chrview_ro==filename1):
        mcg=mcg2-mcg1
        mchg=mchg2-mchg1
        mchh=mchh2-mchh1

    data[0].iloc[:,4]=mcg
    data[0].iloc[:,5]=mchg
    data[0].iloc[:,6]=mchh

    report=data[0]
    report.to_csv('chrView_mean.txt',sep='\t',index=False)

    pd.DataFrame([[chrview_lo+' - '+chrview_ro,'chrView_mean.txt']]).to_csv('chrView_mean_list.txt',sep='\t',index=False,header=None)

    subprocess.call("Rscript --slave chrView-mean.R" , shell=True)
else:
    pass

   
#metaplot
if(Metaplot=='y'):
    #union to bw
    #cgmap=pd.read_csv("union4.txt",sep='\t',dtype={0:str})
    cgmap=union
    genebodybed=input_gene_name+"_"+metaplot_gene_feature+"_merge.bed"
    bwheader=[]
    for chr in cgmap.iloc[:,0].unique():
        k=cgmap[cgmap.iloc[:,0]==chr]
        maxPos=k.iloc[:,1].max()
        bwheader.append((str(chr),maxPos))


    for i in range(3,cgmap.shape[1]):
        #print i
        Samplename=cgmap.columns[i]
        cgoutfile=Samplename+"_CG.bw"
        chgoutfile=Samplename+"_CHG.bw"
        chhoutfile=Samplename+"_CHH.bw"

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

            bw_CG.addEntries(str(chromosome),[x-1 for x in CG.ix[:,1].tolist()],values=CG.ix[:,i].astype(float).tolist(), span=1)
            bw_CHG.addEntries(str(chromosome),[x-1 for x in CHG.ix[:,1].tolist()],values=CHG.ix[:,i].astype(float).tolist(), span=1)
            bw_CHH.addEntries(str(chromosome),[x-1 for x in CHH.ix[:,1].tolist()],values=CHH.ix[:,i].astype(float).tolist(), span=1)

        bw_CG.close()
        bw_CHG.close()
        bw_CHH.close()

        subprocess.call('''computeMatrix scale-regions -S %s -R %s -a 2000 -b 2000 -m 4000 -bs 100 -out %s'''%(Samplename+"_CG.bw",genebodybed,Samplename+"_CG.matrix.gz"),shell=True)
        subprocess.call('''computeMatrix scale-regions -S %s -R %s -a 2000 -b 2000 -m 4000 -bs 100 -out %s'''%(Samplename+"_CHG.bw",genebodybed,Samplename+"_CHG.matrix.gz"),shell=True)
        subprocess.call('''computeMatrix scale-regions -S %s -R %s -a 2000 -b 2000 -m 4000 -bs 100 -out %s'''%(Samplename+"_CHH.bw",genebodybed,Samplename+"_CHH.matrix.gz"),shell=True)

    #metaplot
    subprocess.call("Rscript --slave metaplot.R "+metaplot_gene_feature,shell=True)
    subprocess.call("Rscript --slave metaplot-mean.R "+metaplot_gene_feature+' '+metaplot_lo+' '+metaplot_ro,shell=True)
    
else:
    pass
