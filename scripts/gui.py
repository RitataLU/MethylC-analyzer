import sys
from PyQt5 import QtWidgets as qw
from PyQt5 import QtGui as qg
from PyQt5 import QtCore as qc
import subprocess
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
#import pyBigWig
import threading


class MainWindow(qw.QMainWindow):
    def __init__(self):
        super(MainWindow,self).__init__()
        
        self.startbtn=qw.QPushButton("Start Analysis")
        self.vbbtn=qw.QPushButton("Browser")
        self.statuslabel=qw.QLabel("-")
        self.progressBar=qw.QProgressBar()
        self.statusBar().addPermanentWidget(self.statuslabel)
        self.statusBar().addPermanentWidget(self.progressBar)
        self.statusBar().addPermanentWidget(self.vbbtn)
        self.statusBar().addPermanentWidget(self.startbtn)
        
        self.analysis=RunThread(self.progressBar.value())
        self.analysis.count.connect(self.progressBar.setValue)
        self.analysis.state.connect(self.statuslabel.setText)
        self.startbtn.clicked.connect(self.start)
        self.vbbtn.clicked.connect(self.launchb)
        
        
    def start(self):
        self.analysis.start()
        
    def launchb(self):
        if(w.c==0):
            b.setCentralWidget(DashBoardDefault())
            b.show()
        else:
            b.setCentralWidget(DashBoard())
            b.showMaximized()

        
class Browser(qw.QMainWindow):
    def __init__(self):
        super(Browser,self).__init__()
        
        

class Panel(qw.QWidget):
    def __init__(self):
        super(Panel,self).__init__()
        #para
        self.spfpath=''
        self.cgmfpath=''
        self.igfpath=''
        self.wdpath=''
        
        self.c=0
        self.rs=''
        self.lvv=''
        self.cgv=''
        self.chhv=''
        self.chgv=''
        
        
        #1
        groupBox1 = qw.QGroupBox("Select FIles")

        self.label1 = qw.QLabel("-")
        btn1=qw.QPushButton("Select Sample List File")
        btn1.clicked.connect(self.samplelistfd)
        self.label3 = qw.QLabel("-")
        btn3=qw.QPushButton("Select Input Gene File")
        btn3.clicked.connect(self.inputgenefd)

        layout=qw.QGridLayout()
        layout.addWidget(self.label1,0,0)
        layout.addWidget(btn1,1,0)
        layout.addWidget(self.label3,0,1)
        layout.addWidget(btn3,1,1)
        
        groupBox1.setLayout(layout)

        #2
        groupBox2 = qw.QGroupBox("Setteing Parameters")
        
        label4 = qw.QLabel("Min Depth Value")
        self.tbn1=qw.QLineEdit("4")
        self.tbn1.setValidator(qg.QIntValidator())
        label5 = qw.QLabel("Region Size")
        label5.setToolTip("smaller region size may let analysis more accurate, but it also takes more time and performence, default is 500")
        self.tbn2=qw.QLineEdit("500")
        self.tbn2.setValidator(qg.QIntValidator())
        label6 = qw.QLabel("Region Quaified Value")
        self.tbn3=qw.QLineEdit("4")
        self.tbn3.setValidator(qg.QIntValidator())
        
        #dmr_layout
        dmr_layout=qw.QHBoxLayout()
        label7 = qw.QLabel("DMR_Cut_Off_Value")
        
        label_dmrcg=qw.QLabel("CG")
        label_dmrchh=qw.QLabel("CHH")
        label_dmrchg=qw.QLabel("CHG")
        self.tbn_dmrcg=qw.QLineEdit("0.2")
        self.tbn_dmrcg.setValidator(qg.QDoubleValidator())
        self.tbn_dmrchh=qw.QLineEdit("0.2")
        self.tbn_dmrchh.setValidator(qg.QDoubleValidator())
        self.tbn_dmrchg=qw.QLineEdit("0.2")
        self.tbn_dmrchg.setValidator(qg.QDoubleValidator())
        
        dmrcg_layout=qw.QHBoxLayout()
        dmrcg_layout.addWidget(label_dmrcg)
        dmrcg_layout.addWidget(self.tbn_dmrcg)
        dmrchh_layout=qw.QHBoxLayout()
        dmrchh_layout.addWidget(label_dmrchh)
        dmrchh_layout.addWidget(self.tbn_dmrchh)
        dmrchg_layout=qw.QHBoxLayout()
        dmrchg_layout.addWidget(label_dmrchg)
        dmrchg_layout.addWidget(self.tbn_dmrchg)
        
        dmr_layout.addLayout(dmrcg_layout)
        dmr_layout.addLayout(dmrchh_layout)
        dmr_layout.addLayout(dmrchg_layout)
        
        
        #heat_layout
        heat_layout=qw.QHBoxLayout()
        labelpcacut=qw.QLabel("Heatmap_Cut Off Value")
        
        label_heatcg=qw.QLabel("CG")
        label_heatchh=qw.QLabel("CHH")
        label_heatchg=qw.QLabel("CHG")
        self.tbn_heatcg=qw.QLineEdit("0.2")
        self.tbn_heatcg.setValidator(qg.QDoubleValidator())
        self.tbn_heatchh=qw.QLineEdit("0.2")
        self.tbn_heatchh.setValidator(qg.QDoubleValidator())
        self.tbn_heatchg=qw.QLineEdit("0.2")
        self.tbn_heatchg.setValidator(qg.QDoubleValidator())
        
        heatcg_layout=qw.QHBoxLayout()
        heatcg_layout.addWidget(label_heatcg)
        heatcg_layout.addWidget(self.tbn_heatcg)
        heatchh_layout=qw.QHBoxLayout()
        heatchh_layout.addWidget(label_heatchh)
        heatchh_layout.addWidget(self.tbn_heatchh)
        heatchg_layout=qw.QHBoxLayout()
        heatchg_layout.addWidget(label_heatchg)
        heatchg_layout.addWidget(self.tbn_heatchg)
        
        heat_layout.addLayout(heatcg_layout)
        heat_layout.addLayout(heatchh_layout)
        heat_layout.addLayout(heatchg_layout)
        
        
        labelbinsize=qw.QLabel("Bin Size")
        self.tbnbinsize=qw.QLineEdit("1000000")
        self.tbnbinsize.setValidator(qg.QIntValidator())
        labelpromotor=qw.QLabel("Promotor Size")
        self.tbnpromotor=qw.QLineEdit("2000")
        self.tbnpromotor.setValidator(qg.QIntValidator())
        
        layout=qw.QGridLayout()
        layout.addWidget(label4,0,0)
        layout.addWidget(self.tbn1,0,1)
        layout.addWidget(label5,1,0)
        layout.addWidget(self.tbn2,1,1)
        layout.addWidget(label6,2,0)
        layout.addWidget(self.tbn3,2,1)
        layout.addWidget(label7,3,0)
        layout.addLayout(dmr_layout,3,1)
        layout.addWidget(labelpcacut,6,0)
        layout.addLayout(heat_layout,6,1)
        layout.addWidget(labelbinsize,7,0)
        layout.addWidget(self.tbnbinsize,7,1)
        layout.addWidget(labelpromotor,8,0)
        layout.addWidget(self.tbnpromotor,8,1)
        
        groupBox2.setLayout(layout)

        #3
        groupBox3 = qw.QGroupBox("Choosing Tools")
        
        self.cb1 = qw.QCheckBox("HeatMap & PCA")
        self.cb2 = qw.QCheckBox("Enrichment")
        self.cb3 = qw.QCheckBox("DMG")
        self.cb4 = qw.QCheckBox("ChrView")
        self.chrbtn=qw.QPushButton("..")
        self.chrbtn.setEnabled(False)
        self.chrbtn.clicked.connect(self.launchchrset)
        self.cb5 = qw.QCheckBox("MetaPlot")
        self.metabtn=qw.QPushButton("..")
        self.metabtn.setEnabled(False)
        self.metabtn.clicked.connect(self.launchmetaset)
        
        layout=qw.QGridLayout()
        layout.addWidget(self.cb1,0,0)
        layout.addWidget(self.cb2,1,0)
        layout.addWidget(self.cb3,2,0)
        layout.addWidget(self.cb4,3,0)
        layout.addWidget(self.chrbtn,3,1)
        layout.addWidget(self.cb5,4,0)
        layout.addWidget(self.metabtn,4,1)
        
        groupBox3.setLayout(layout)
        
        #4
        groupBox4 = qw.QGroupBox("Output Setting")
        self.label10 = qw.QLabel("./")
        self.btn4=qw.QPushButton("Choose Work Directory")
        self.btn4.clicked.connect(self.workdirfd)
        self.btn4.setEnabled(False)
        self.cb6 = qw.QCheckBox("Auto shutdown PC")

        
        layout=qw.QGridLayout()
        layout.addWidget(self.label10,0,0)
        layout.addWidget(self.btn4,1,0)
        layout.addWidget(self.cb6,2,0)
        
        groupBox4.setLayout(layout)
        
        groupBox5 = qw.QGroupBox()
        layout2=qw.QGridLayout()
        layout2.addWidget(groupBox3,0,0)
        layout2.addWidget(groupBox4,0,1)
        groupBox5.setLayout(layout2)
        
        
        
        #layout
        grid1=qw.QGridLayout()
        grid1.addWidget(groupBox1,0,0)
        grid1.addWidget(groupBox2,1,0)
        grid1.addWidget(groupBox5,2,0)
        self.setLayout(grid1)
    
    def samplelistfd(self):
        self.spfpath, filetype = qw.QFileDialog.getOpenFileName(self,"Select File","./","All Files (*);;Text Files (*.txt)")
        self.label1.setText(self.spfpath)
        if(self.spfpath!=''):
            self.chrbtn.setEnabled(True)
            self.metabtn.setEnabled(True)
        else:
            self.chrbtn.setEnabled(False)
            self.metabtn.setEnabled(False)
        
    def inputgenefd(self):
        self.igfpath, filetype = qw.QFileDialog.getOpenFileName(self,"Select File","./","All Files (*);;gtf Files (*.gtf)")
        self.label3.setText(self.igfpath)
        
    def workdirfd(self):
        self.wdpath = qw.QFileDialog.getExistingDirectory(self,"Choose Work Directory")
        self.label10.setText(self.wdpath)
        
    def launchchrset(self):
        samples = pd.read_csv(self.label1.text(),header=None,sep="\t")
        obj=samples.iloc[:,2].unique()
        
        cc.combox1.clear()
        cc.combox2.clear()
        
        cc.combox1.addItems(obj)
        cc.combox2.addItems(obj)
        
        cc.show()
        
    def launchmetaset(self):
        samples = pd.read_csv(self.label1.text(),header=None,sep="\t")
        obj=samples.iloc[:,2].unique()

        mm.combox1.clear()
        mm.combox2.clear()
        
        mm.combox1.addItems(obj)
        mm.combox2.addItems(obj)
        
        mm.show()
        
        
        
class DashBoard(qw.QWidget):
    def __init__(self):
        super(qw.QWidget, self).__init__()
        self.layout = qw.QVBoxLayout(self)
 
        # Initialize tab screen
        self.tabs = qw.QTabWidget()
        
        #tab1
        self.tab1 = qw.QWidget()
        console1 = cgtab()
        console1.setMinimumSize(1500, 600)
        
        scroll = qw.QScrollArea()
        scroll.setWidget(console1)
        scroll.setAutoFillBackground(True)
        scroll.setWidgetResizable(False)
        vbox = qw.QVBoxLayout()
        vbox.addWidget(scroll)  
        self.tab1.setLayout(vbox)
        
        #tab2
        self.tab2 = qw.QWidget()
        console2 = chhtab()
        console2.setMinimumSize(1500, 600)
        
        scroll = qw.QScrollArea()
        scroll.setWidget(console2)
        scroll.setAutoFillBackground(True)
        scroll.setWidgetResizable(False)
        vbox = qw.QVBoxLayout()
        vbox.addWidget(scroll)  
        self.tab2.setLayout(vbox)
        
        #tab3
        self.tab3 = qw.QWidget()
        console3 = chgtab()
        console3.setMinimumSize(1500, 600)
        
        scroll = qw.QScrollArea()
        scroll.setWidget(console3)
        scroll.setAutoFillBackground(True)
        scroll.setWidgetResizable(False)
        vbox = qw.QVBoxLayout()
        vbox.addWidget(scroll)  
        self.tab3.setLayout(vbox)
        
        
        self.tabs.addTab(self.tab1,"CG")
        self.tabs.addTab(self.tab2,"CHH")
        self.tabs.addTab(self.tab3,"CHG")
        # Add tabs to widget
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
        
    
class RunThread(qc.QThread):
    count=qc.pyqtSignal(int)
    state=qc.pyqtSignal(str)
    
    def __init__(self, parent=None):
        super(RunThread, self).__init__()
 
    def __del__(self):
        self.wait()
        
    def run(self):
        mw.startbtn.setEnabled(False)
        if(w.spfpath=='' or w.igfpath==''):
            w.c=1
            print('Required Files are not enough')
            mw.startbtn.setEnabled(True)
        elif(w.rs=='t' or w.lvv=='t' or w.cgv=='t' or w.chhv=='t' or w.chgv=='t'):
            print('Parameters cannot be empty or iligal value')
            mw.startbtn.setEnabled(True)
        else:
            #TotalPipeline
            
            self.state.emit('Preparing..')
            
            path=w.wdpath #working path
            heatmap_check=w.cb1.checkState()
            enrichment_check=w.cb2.checkState()
            dmg_check=w.cb3.checkState()
            chrview_check=w.cb4.checkState()
            metaplot_check=w.cb5.checkState()
			
            autoshutdown_check=w.cb6.checkState()
            
            
            self.state.emit('Setting Parameters')
            
            depth=int(w.tbn1.text()) #number of C+T
            region=int(w.tbn2.text())
            qualifiedSite =int(w.tbn3.text()) #number of least valid region count
            pca_heat_cg_cut=float(w.tbn_heatcg.text()) #heatmap least valid max-min to display
            pca_heat_chh_cut=float(w.tbn_heatchh.text()) #heatmap least valid max-min to display
            pca_heat_chg_cut=float(w.tbn_heatchg.text()) #heatmap least valid max-min to display
            cutoff_CG=float(w.tbn_dmrcg.text()) #DMR_CG cutoff value
            cutoff_CHH=float(w.tbn_dmrchh.text()) #DMR_CHH cutoff value
            cutoff_CHG=float(w.tbn_dmrchg.text()) #DMR_CHG cutoff value
            promoter_value=2000 #range of promotor (In front of gene start)
            binSize=1000000 #accurucy of chrview and metaplot
            input_gene_name=w.label3.text() #input gene file name (obtain path)
            
            
            
            self.state.emit('Loading input gene file')
            self.count.emit(1)
            input_gene=pd.read_csv(w.label3.text(),sep='\t',header=None)

            #parser.add_argument('input_genome')

            #input_genome=input_genome
            #input_gene=input_gene


            #samtools faidx ref.fasta %out=fai
            #def genome_size(input_genome):
                #subprocess.call('''samtools faidx %s'''%(input_genome), shell=True)
                #genomefai=pd.read_csv(input_genome+".fai",sep='\t',header=None)
                #genome_size=genomefai[1].sum()*1.0
                #return genome_size

            context=["CG","CHG","CHH"]
            #unionsite
            
            combined = pd.DataFrame()
            
            self.state.emit('Reading sample list..')
            self.count.emit(2)
            samples = pd.read_csv(w.label1.text(),header=None,sep="\t")
            
            self.state.emit('Generating Unionsite Table')
            self.count.emit(4)
            
            for sample in samples.itertuples():
				#elapsed_time = time.time() - start_time
				print("Now processing " + sample[2])
				CGmap = pd.read_csv(sample[2], header=None,sep="\t",dtype =
						{0:str,2:str,3:str,5:str,6:str,7:int},usecols=[0,2,3,5,6,7],index_col=[0,1,2],compression='gzip')
				#elapsed_time = time.time() - start_time
				#print("File read " + time.strftime("%H:%M:%S",time.gmtime(elapsed_time)))

				CGmap = CGmap.loc[CGmap[7] >= depth,5:5]
				#elapsed_time = time.time() - start_time
				#print("Filtered " + time.strftime("%H:%M:%S",time.gmtime(elapsed_time)))
				sample_col = [sample[1]]
				#sample_col = [sample[1] + '_meth',sample[1] + '_mC',sample[1] + '_C']
				CGmap.columns = sample_col
				if len(combined) == 0:
					combined = CGmap
				else:
					combined = pd.merge(combined,CGmap,left_index=True,right_index=True,how = 'outer')
				#elapsed_time = time.time() - start_time
				#print("Combined " + time.strftime("%H:%M:%S",time.gmtime(elapsed_time)))

				combined.reset_index(inplace=True)
				combined.rename(columns = {0:'chr', 2:'pos', 3: 'context'},inplace=True)
				combined[['pos']] = combined[['pos']].astype(int)
			
            #process chromosome name
			#deleted
			
            combined = combined.sort_values(['chr','pos'],ascending=[True,True])
            #elapsed_time = time.time() - start_time
            #print("Sorted " + time.strftime("%H:%M:%S",time.gmtime(elapsed_time)))
            combined.to_csv('Unionsite.txt',sep = '\t',na_rep='-',index=False)
            
			
            #common region
            self.state.emit('Reading Unionsite File')
            self.count.emit(24)
            combined=pd.read_csv('Unionsite.txt',sep='\t')
            union=combined
            #1.unionsite --> combined
            chrs = union['chr'].unique()
            
            self.state.emit('Calculating CommonRegion..')
            self.count.emit(26)
            
            #context one by one
            for cxt in context:
                outfile = 'CommonRegion_' + cxt + '.txt'
                with open(outfile,'w') as of:
                    of.write('\t'.join(['chr','start','end'] + union.columns.values.tolist()[3:]) + "\n")
                    for chromosome in chrs:
                        subset = union[(union['context'] == cxt) & (union['chr'] == chromosome)]
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
            

            #DMR
            filename1=samples.iloc[0,2] #read group name via samples (read in upper step)
            filename2_temp=samples
            filename2=filename2_temp.iloc[len(filename2_temp.index)-1,2]
            
            self.state.emit('Loading CommonRegion Files')
            self.count.emit(46)
            file1=pd.read_csv("CommonRegion_CG.txt",sep="\t",header=None)
            file2=pd.read_csv("CommonRegion_CHH.txt",sep="\t",header=None)
            file3=pd.read_csv("CommonRegion_CHG.txt",sep="\t",header=None)

            CG_data=np.array(file1)
            CHH_data=np.array(file2)
            CHG_data=np.array(file3)


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
            CG_delta=np.append(cutoff_CG,CG_data[1:,index2_CG].astype(np.float).mean(axis=1)-CG_data[1:,index1_CG].astype(np.float).mean(axis=1)) #diffrence of average

            CHH_delta=np.array(0.0)
            CHH_delta=np.append(cutoff_CHH,CHH_data[1:,index2_CHH].astype(np.float).mean(axis=1)-CHH_data[1:,index1_CHH].astype(np.float).mean(axis=1)) #diffrence of average

            CHG_delta=np.array(0.0)
            CHG_delta=np.append(cutoff_CHG,CHG_data[1:,index2_CHG].astype(np.float).mean(axis=1)-CHG_data[1:,index1_CHG].astype(np.float).mean(axis=1)) #diffrence of average


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


            self.state.emit('Generating DMR Table')
            self.count.emit(48)
            CG_DMRpd[CG_DMRpd['delta']>0].to_csv('DMR_CG_hyper.txt',sep='\t',index=None)
            CG_DMRpd[CG_DMRpd['delta']<0].to_csv('DMR_CG_hypo.txt',sep='\t',index=None)
            CG_DMRpd.to_csv('DMR_CG_all.txt',sep='\t',index=None)

            CHH_DMRpd[CHH_DMRpd['delta']>0].to_csv('DMR_CHH_hyper.txt',sep='\t',index=None)
            CHH_DMRpd[CHH_DMRpd['delta']<0].to_csv('DMR_CHH_hypo.txt',sep='\t',index=None)
            CHH_DMRpd.to_csv('DMR_CHH_all.txt',sep='\t',index=None)

            CHG_DMRpd[CHG_DMRpd['delta']>0].to_csv('DMR_CHG_hyper.txt',sep='\t',index=None)
            CHG_DMRpd[CHG_DMRpd['delta']<0].to_csv('DMR_CHG_hypo.txt',sep='\t',index=None)
            CHG_DMRpd.to_csv('DMR_CHG_all.txt',sep='\t',index=None)

            
            if(heatmap_check==2):
                #heatmap_PCA
                self.state.emit('Graphing heatmap and PCA plot')
                self.count.emit(50)
                subprocess.call('''Rscript --slave heatmap_PCA_all.R %s %s'''%("CommonRegion_CG.txt",pca_heat_cg_cut), shell=True)
                subprocess.call('''Rscript --slave heatmap_PCA_all.R %s %s'''%("CommonRegion_CHG.txt",pca_heat_chh_cut), shell=True)
                subprocess.call('''Rscript --slave heatmap_PCA_all.R %s %s'''%("CommonRegion_CHH.txt",pca_heat_chg_cut), shell=True)
            else:
                self.count.emit(52)
            self.count.emit(52)


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
            
            self.state.emit('Generating enrichment Table')
            self.count.emit(52)
            subprocess.call('''./extract_transcript_regions.py -i %s -o %s --gtf'''%(input_gene_name,input_gene_name), shell=True)
            
            self.count.emit(54)
            

            print ("*-------------------------------------*")
            print ("|Convert this blockbed (bed12) to bed6|")
            print ("*-------------------------------------*")

            for i in annotation_name:
                subprocess.call('''cat %s | bed12ToBed6 -i stdin -n > %s'''%(input_gene_name+'_'+i+'.bed',input_gene_name+'_'+i+'_bed6.bed'),shell=True)

            self.count.emit(56)
                
            #find gene_body.bed
            genes = pd.read_csv(input_gene_name, header=None, sep="\t",dtype = {0 :str})
            #genes=genes[genes[0].str.len()<=5]
            genes.columns=['chr','unknow', 'exon', 'g_str', 'g_end', 'g_score', 'g_dir','.', 'gene_name']
            genes=genes[genes.exon=='exon']
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
            
            self.count.emit(57)
            
            
            #find promoter.bed
            gene_body['pro_str'] = np.where(gene_body.g_dir == '+', gene_body.g_str - promoter_value, gene_body.g_end - 0)
            gene_body['pro_end'] = np.where(gene_body.g_dir == '+', gene_body.g_str + 0, gene_body.g_end + promoter_value)
            num = gene_body._get_numeric_data()
            num[num<0]=0
            gene_promoter = gene_body.ix[:, ['chr','pro_str','pro_end','gene_id','g_score','g_dir']]
            gene_promoter.columns=['chr','g_str','g_end','gene_id','g_score','g_dir']
            gene_promoter.to_csv(input_gene_name+"_gene_promoter_bed6.bed", sep='\t',index=False, header=None)
            
            self.count.emit(58)


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
            self.count.emit(59)


            for i in annotation_name:
                subprocess.call('''bedtools sort -i %s|bedtools merge -c 4,5,6 -o collapse,collapse,collapse >%s '''%(input_gene_name+'_'+i+'_bed6.bed',input_gene_name +'_'+i+'_merge.bed'),shell=True)
            
            self.count.emit(60)
            #input_gene_exon_merge.bed


            #enrichment
            #1. overlap with dmr
            def overlap(bed1,bed2):
                overlap=subprocess.check_output("bedtools intersect -a %s -b %s |awk '{size+=$3 - $2}END{print size}' "%(bed1,bed2), shell=True)
                if len(overlap)<2:
                    return 1
                else:
                    return int(overlap)

            
            if(enrichment_check==2):
                #DMR enrichmet cal & plot
                self.state.emit('Graphing enrichment plot')
                self.count.emit(60)
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
                self.count.emit(63)
            self.count.emit(63)
            
            
            #DMG
            def dmg(tag):
                #subprocess.call('''bedtools intersect -a %s -b %s >%s '''%(input_gene_name +'_gene_body_merge.bed',"DMR_"+tag+"_"+str(dmr_cut)+"_all.txt.bed","DMG_"+tag+"_gene_list.txt"),shell=True)
                #subprocess.call('''bedtools intersect -a %s -b %s >%s '''%(input_gene_name +'_gene_promoter_merge.bed',"DMR_"+tag+"_"+str(dmr_cut)+"_all.txt.bed","DMG_"+tag+"_promoter_list.txt"),shell=True)
                subprocess.call('''bedtools intersect -a %s -b %s >%s '''%(input_gene_name +'_gene_body_merge.bed',"DMR_"+tag+"_all.txt.bed","DMG_"+tag+"_gene_list.txt"),shell=True)
                subprocess.call('''bedtools intersect -a %s -b %s >%s '''%(input_gene_name +'_gene_promoter_merge.bed',"DMR_"+tag+"_all.txt.bed","DMG_"+tag+"_promoter_list.txt"),shell=True)
            
            
            if(dmg_check==2):
                self.state.emit('Generating DMG Table')
                self.count.emit(63)
                for k in context:
                    dmg(k)
            else:
                self.count.emit(64)
            self.count.emit(64)

            #----------------------------------------------------------------------------------------------------------------------------------------------------------
            #chrView
            savefile=pd.DataFrame()
            
            if(chrview_check==2):
                self.state.emit('Generating chrView Table')
                self.count.emit(64)
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
                
                if(cc.cc1==filename1 and cc.cc2==filename2):
                    mcg=mcg1-mcg2
                    mchg=mchg1-mchg2
                    mchh=mchh1-mchh2
                elif(cc.cc1==filename2 and cc.cc2==filename1):
                    mcg=mcg2-mcg1
                    mchg=mchg2-mchg1
                    mchh=mchh2-mchh1

                data[0].iloc[:,4]=mcg
                data[0].iloc[:,5]=mchg
                data[0].iloc[:,6]=mchh

                report=data[0]
                report.to_csv('chrView_mean.txt',sep='\t',index=False)

                pd.DataFrame([[cc.cc1+' - '+cc.cc2,'chrView_mean.txt']]).to_csv('chrView_mean_list.txt',sep='\t',index=False,header=None)
                
                subprocess.call("Rscript --slave chrView-mean.R" , shell=True)
                
            else:
                self.count.emit(72)
            self.count.emit(72)

            
            
            #metaplot
            if(metaplot_check==2):
                #union to bw
                #cgmap=pd.read_csv("union4.txt",sep='\t',dtype={0:str})
                self.state.emit('Generating metaplot Files')
                self.count.emit(72)
                cgmap=union
                genebodybed=input_gene_name+"_"+mm.mm3+"_merge.bed"
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
                        CG = k[(k.iloc[:,2] == 'CG') & (k.iloc[:,i] != "-")]
                        CHG = k[(k.iloc[:,2] == 'CHG') & (k.iloc[:,i] != "-")]
                        CHH = k[(k.iloc[:,2] == 'CHH') & (k.iloc[:,i] != "-")]

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
                self.state.emit('Graphing metaplot')
                self.count.emit(96)
                subprocess.call("Rscript --slave metaplot.R "+mm.mm3,shell=True)
                subprocess.call("Rscript --slave metaplot-mean.R "+mm.mm3+' '+mm.mm1+' '+mm.mm2,shell=True)
                
            else:
                self.count.emit(82)
            self.count.emit(82)
            
            
            self.state.emit('Done')
            self.count.emit(100)
            
            w.c=1
			
            if(autoshutdown_check==2):
                subprocess.call("shutdown -h +3") #for linux
                subprocess.call("shutdown -s -t 180") #for windows
			
            mw.startbtn.setEnabled(True)
            #self.wait()

class cgtab(qw.QWidget):
    def __init__(self):
        super(cgtab,self).__init__()
        self.layout=qw.QGridLayout()
        
        self.label1 = qw.QLabel(self)
        self.pixmap1 = qg.QPixmap('heatmap_CG_0.2.png')
        self.pixmap1=self.pixmap1.scaled(qc.QSize(1280,800),qc.Qt.KeepAspectRatio)
        self.label1.setPixmap(self.pixmap1)
        self.label1.resize(self.pixmap1.width(),self.pixmap1.height())
        
        self.label2 = qw.QLabel(self)
        self.pixmap2 = qg.QPixmap('PCA_CG_0.2.jpeg')
        self.pixmap2=self.pixmap2.scaled(qc.QSize(1280,1280),qc.Qt.KeepAspectRatio)
        self.label2.setPixmap(self.pixmap2)
        self.label2.resize(self.pixmap2.width(),self.pixmap2.height())
        
        self.label3 = qw.QLabel(self)
        self.pixmap3 = qg.QPixmap('CG_Fold_Enrichment.png')
        self.pixmap3=self.pixmap3.scaled(qc.QSize(1280,1280),qc.Qt.KeepAspectRatio)
        self.label3.setPixmap(self.pixmap3)
        self.label3.resize(self.pixmap3.width(),self.pixmap3.height())
        
        self.label4 = qw.QLabel(self)
        self.pixmap4 = qg.QPixmap('chrViewq_CG.png')
        #self.pixmap4=self.pixmap4.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label4.setPixmap(self.pixmap4)
        self.label4.resize(self.pixmap4.width(),self.pixmap4.height())
        
        self.label6 = qw.QLabel(self)
        self.pixmap6 = qg.QPixmap('chrViewq_mean_CG.png')
        #self.pixmap6=self.pixmap6.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label6.setPixmap(self.pixmap6)
        self.label6.resize(self.pixmap6.width(),self.pixmap6.height())
        
        self.label5 = qw.QLabel(self)
        self.pixmap5 = qg.QPixmap('metaplot_CG.png')
        #self.pixmap5=self.pixmap5.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label5.setPixmap(self.pixmap5)
        self.label5.resize(self.pixmap5.width(),self.pixmap5.height())
        
        self.label7 = qw.QLabel(self)
        self.pixmap7 = qg.QPixmap('metaplot_mean_CG.png')
        #self.pixmap7=self.pixmap7.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label7.setPixmap(self.pixmap7)
        self.label7.resize(self.pixmap7.width(),self.pixmap7.height())
        
        self.layout = qw.QVBoxLayout(self)
        self.layout.addWidget(self.label1)
        self.layout.addWidget(self.label2)
        self.layout.addWidget(self.label3)
        self.layout.addWidget(self.label4)
        self.layout.addWidget(self.label6)
        self.layout.addWidget(self.label5)
        self.layout.addWidget(self.label7)
        self.setLayout(self.layout)
        
        
class chhtab(qw.QWidget):
    def __init__(self):
        super(chhtab,self).__init__()
        self.layout=qw.QGridLayout()
        
        self.label1 = qw.QLabel(self)
        self.pixmap1 = qg.QPixmap('heatmap_CHH_0.2.png')
        self.pixmap1=self.pixmap1.scaled(qc.QSize(1280,800),qc.Qt.KeepAspectRatio)
        self.label1.setPixmap(self.pixmap1)
        self.label1.resize(self.pixmap1.width(),self.pixmap1.height())
        
        self.label2 = qw.QLabel(self)
        self.pixmap2 = qg.QPixmap('PCA_CHH_0.2.jpeg')
        self.pixmap2=self.pixmap2.scaled(qc.QSize(1280,1280),qc.Qt.KeepAspectRatio)
        self.label2.setPixmap(self.pixmap2)
        self.label2.resize(self.pixmap2.width(),self.pixmap2.height())
        
        self.label3 = qw.QLabel(self)
        self.pixmap3 = qg.QPixmap('CHH_Fold_Enrichment.png')
        self.pixmap3=self.pixmap3.scaled(qc.QSize(1280,1280),qc.Qt.KeepAspectRatio)
        self.label3.setPixmap(self.pixmap3)
        self.label3.resize(self.pixmap3.width(),self.pixmap3.height())
        
        self.label4 = qw.QLabel(self)
        self.pixmap4 = qg.QPixmap('chrViewq_CHH.png')
        #self.pixmap4=self.pixmap4.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label4.setPixmap(self.pixmap4)
        self.label4.resize(self.pixmap4.width(),self.pixmap4.height())
        
        self.label6 = qw.QLabel(self)
        self.pixmap6 = qg.QPixmap('chrViewq_mean_CHH.png')
        #self.pixmap6=self.pixmap6.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label6.setPixmap(self.pixmap6)
        self.label6.resize(self.pixmap6.width(),self.pixmap6.height())
        
        self.label5 = qw.QLabel(self)
        self.pixmap5 = qg.QPixmap('metaplot_CHH.png')
        #self.pixmap5=self.pixmap5.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label5.setPixmap(self.pixmap5)
        self.label5.resize(self.pixmap5.width(),self.pixmap5.height())
        
        self.label7 = qw.QLabel(self)
        self.pixmap7 = qg.QPixmap('metaplot_mean_CHH.png')
        #self.pixmap7=self.pixmap7.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label7.setPixmap(self.pixmap7)
        self.label7.resize(self.pixmap7.width(),self.pixmap7.height())
        
        self.layout = qw.QVBoxLayout(self)
        self.layout.addWidget(self.label1)
        self.layout.addWidget(self.label2)
        self.layout.addWidget(self.label3)
        self.layout.addWidget(self.label4)
        self.layout.addWidget(self.label6)
        self.layout.addWidget(self.label5)
        self.layout.addWidget(self.label7)
        self.setLayout(self.layout)

class chgtab(qw.QWidget):
    def __init__(self):
        super(chgtab,self).__init__()
        self.layout=qw.QGridLayout()
        
        self.label1 = qw.QLabel(self)
        self.pixmap1 = qg.QPixmap('heatmap_CHG_0.2.png')
        self.pixmap1=self.pixmap1.scaled(qc.QSize(1280,800),qc.Qt.KeepAspectRatio)
        self.label1.setPixmap(self.pixmap1)
        self.label1.resize(self.pixmap1.width(),self.pixmap1.height())
        
        self.label2 = qw.QLabel(self)
        self.pixmap2 = qg.QPixmap('PCA_CHG_0.2.jpeg')
        self.pixmap2=self.pixmap2.scaled(qc.QSize(1280,1280),qc.Qt.KeepAspectRatio)
        self.label2.setPixmap(self.pixmap2)
        self.label2.resize(self.pixmap2.width(),self.pixmap2.height())
        
        self.label3 = qw.QLabel(self)
        self.pixmap3 = qg.QPixmap('CHG_Fold_Enrichment.png')
        self.pixmap3=self.pixmap3.scaled(qc.QSize(1280,1280),qc.Qt.KeepAspectRatio)
        self.label3.setPixmap(self.pixmap3)
        self.label3.resize(self.pixmap3.width(),self.pixmap3.height())
        
        self.label4 = qw.QLabel(self)
        self.pixmap4 = qg.QPixmap('chrViewq_CHG.png')
        #self.pixmap4=self.pixmap4.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label4.setPixmap(self.pixmap4)
        self.label4.resize(self.pixmap4.width(),self.pixmap4.height())
        
        self.label6 = qw.QLabel(self)
        self.pixmap6 = qg.QPixmap('chrViewq_mean_CHG.png')
        #self.pixmap6=self.pixmap6.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label6.setPixmap(self.pixmap6)
        self.label6.resize(self.pixmap6.width(),self.pixmap6.height())
        
        self.label5 = qw.QLabel(self)
        self.pixmap5 = qg.QPixmap('metaplot_CHG.png')
        #self.pixmap5=self.pixmap5.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label5.setPixmap(self.pixmap5)
        self.label5.resize(self.pixmap5.width(),self.pixmap5.height())
        
        self.label7 = qw.QLabel(self)
        self.pixmap7 = qg.QPixmap('metaplot_mean_CHG.png')
        #self.pixmap7=self.pixmap7.scaled(qc.QSize(1600,1600),qc.Qt.KeepAspectRatio)
        self.label7.setPixmap(self.pixmap7)
        self.label7.resize(self.pixmap7.width(),self.pixmap7.height())
        
        self.layout = qw.QVBoxLayout(self)
        self.layout.addWidget(self.label1)
        self.layout.addWidget(self.label2)
        self.layout.addWidget(self.label3)
        self.layout.addWidget(self.label4)
        self.layout.addWidget(self.label6)
        self.layout.addWidget(self.label5)
        self.layout.addWidget(self.label7)
        self.setLayout(self.layout)

            
class DashBoardDefault(qw.QWidget):
    def __init__(self):
        super(DashBoardDefault, self).__init__()
        
        self.layout = qw.QGridLayout(self)
        dtext=qw.QLabel('This will show anlysis report when analysis in done')
        
        self.layout.addWidget(dtext,0,0)
        self.setLayout(self.layout)

        
class chrsetting(qw.QWidget):
    def __init__(self):
        super(chrsetting, self).__init__()
        self.cc1=''
        self.cc2=''
        
        self.layout = qw.QGridLayout(self)
        
        self.combox1=qw.QComboBox()
        self.label=qw.QLabel(' - ')
        self.combox2=qw.QComboBox()
        self.subbtn=qw.QPushButton('Apply')
        self.subbtn.clicked.connect(self.savepara)
        
        self.layout.addWidget(self.combox1,0,0)
        self.layout.addWidget(self.label,0,1)
        self.layout.addWidget(self.combox2,0,2)
        self.layout.addWidget(self.subbtn,1,3)
        self.setLayout(self.layout)
        
    def savepara(self):
        self.cc1=self.combox1.currentText()
        self.cc2=self.combox2.currentText()
        self.close()
        
        
class metasetting(qw.QWidget):
    def __init__(self):
        super(metasetting, self).__init__()
        self.mm1=''
        self.mm2=''
        self.mm3=''
        
        self.layout = qw.QGridLayout(self)
        
        self.combox1=qw.QComboBox()
        self.label=qw.QLabel(' - ')
        self.combox2=qw.QComboBox()
        self.combox3=qw.QComboBox()
        self.subbtn=qw.QPushButton('Apply')
        self.subbtn.clicked.connect(self.savepara)
        
        gene = "gene_body"
        exon = "exons"
        intron = "introns"
        utr3 = "3utr"
        utr5 = "5utr"
        cds = "cds"
        promoter = "gene_promoter"
        igr = "gene_igr"
        annotation_name=[promoter,gene,exon,intron,utr5,cds,utr3,igr]
        
        self.combox3.addItems(annotation_name)
        
        self.layout.addWidget(self.combox1,0,0)
        self.layout.addWidget(self.label,0,1)
        self.layout.addWidget(self.combox2,0,2)
        self.layout.addWidget(self.combox3,0,3)
        self.layout.addWidget(self.subbtn,1,4)
        self.setLayout(self.layout)
    
    def savepara(self):
        self.mm1=self.combox1.currentText()
        self.mm2=self.combox2.currentText()
        self.mm3=self.combox3.currentText()
        self.close()
        
            
if __name__=="__main__":
    app=qw.QApplication(sys.argv)

    w=Panel()
    mw=MainWindow()
    mw.setCentralWidget(w)
    mw.show()
    
    b=Browser()
    
    cc=chrsetting()
    mm=metasetting()

    sys.exit(app.exec_())
    



