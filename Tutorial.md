# Tutorial


# <a name="SystemRequirements"></a>System Requirements
* Linux or Mac OS Environment
* CPU：No special restrictions, but CPU has 16 cores is more efficient

* MEM：12GB or higer (for plant sample) / 256GB or higher (for human sample)

* Python 3.9 and R>3.6
* [deepTools](https://deeptools.readthedocs.org/)
* [BEDtools](http://bedtools.readthedocs.org/)
* Python module
```
 numpy
 pandas
 math
 scipy
 matplolib
 argparse
 glob
 pyBigWig
 PyQt5
 seaborn
```
* R package
 ```
  ComplexHeatmap
  gplots
  ggplot2
  viridis
 ```
# Installation
   
1. Recommand to create a [conda](https://docs.conda.io/en/latest/miniconda.html) environment somewhere on your disk, and then activate it.
  
  ```
  $ conda create -n methylC_analyzer_env python=3.9
  $ conda activate Yout conda path methylC_analyzer_env

 ```
3. Download the source code and install the requirements.

  ```
  $ git clone https://github.com/RitataLU/MethylC-analyzer.git
  
 ```
4. Install Package - Run MethylC-analyzer/requirements/base.txt

  e.g.
  ```  
  sh MethylC-analyzer/requirements/base.txt
  ```


# Run demo 
1. go to a folder named 'Demo' where demo files are and decompress 'demo.gz' 

e.g.
```
$ cd MethylC-analyzer/Demo
$ tar -xvf demo.gz 
```
2. Make sure all scripts and demo files are in the same foder
  e.g. move scripts in the 'MethylC-analyzer/Demo/demo'
```
$ cd demo
$ cp ../../scripts/* ./

```
3. Download demo CG map files 
> 


1. Make a sample description file and name it as "samples_list.txt" in the location where methylc.py script. The file should be tab-delimited without a header.

> Sample Description File (tab-delimited, no header in the first line)
Sample list ( sample_name  CGmap_location  group )


**samples list format:**
>> sample_name  CGmap_location  group (tab-delimited, no header in the first line)
```    
MT1     MT1s.CGmap.gz   MT
MT2     MT2s.CGmap.gz   MT
WT1     WT1s.CGmap.gz   WT
WT2     WT2s.CGmap.gz   WT
```

**Input:**
1. gene annotation (GTF)

2. CGmap 

**Usage:**
```
$ python MethylC.py samples_list.txt TAR10_2_demo.gtf -a WT -b MT -d 4 -r 200 -q 2 -bs 200000 -pvalue 0.1

usage: MethylC_new.py [-h] [-a GROUP1] [-b GROUP2] [-d DEPTH] [-r REGION]
                      [-q QUALIFIED] [-context CONTEXT] [-hc HEATMAP_CUTOFF]
                      [-dmrc DMR_CUTOFF] [-test TESTMETHOD] [-pvalue PVALUE]
                      [-bs BIN_SIZE] [-p PROMOTER_SIZE]
                      samples_list input_gtf_file

positional arguments:
  samples_list        samples CGmap description
  input_gtf_file      path of gene annotation

optional arguments:
  -h, --help          show this help message and exit
  -a GROUP1           Name of group1
  -b GROUP2           Name of group2
  -d DEPTH            Minimum depth of sites. Default=4
  -r REGION           Size of region. Default=500
  -q QUALIFIED        Minimum sites within a region. Default=4
  -context CONTEXT    Context used. Default=CG
  -hc HEATMAP_CUTOFF  Methylation cutoff of PCA & Heatmap. Default = 0.2
  -dmrc DMR_CUTOFF    Methylation cutoff of DMR. Default = 0.1
  -test TESTMETHOD    DMR testing method. 0:TTest, 1:KS, 2:MWU. Default=0
  -pvalue PVALUE      p-value cutoff for identifying DMR. Default = 0.05
  -bs BIN_SIZE        Bin size of chrView and Metaplot. Default = 1000000
  -p PROMOTER_SIZE    promoter_size
    ```

** activate interface (choose analysis want to process)
```
Heatmap & PCA Analysis?  (y/n): y
Identify DMR?  (y/n): y
Identify DMG?  (y/n): y
Use Fold Enrichment Analysis?  (y/n): y
Chromosome View Analysis?  (y/n): y
Metaplot Analysis?  (y/n): y

enter experimental group name analysis: MT
enter control group name analysis: WT
```
**Output Figures**


* Average methylation level
* Heatmap & PCA for variable regions 
* PCA for variable regions 
* Identifying Differentially Methylated Regions (DMRs)
* Genomic regions fold enrichment analysis for DMRs 
* Identifying Differentially Methylated Genes (DMGs)
* The distribution fo DNA methylation on each chromosome
* Metaplot for each profile & comparison between groups 



