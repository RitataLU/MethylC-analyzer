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
1. Download [Demo](https://paoyang.ipmb.sinica.edu.tw/MethylC-analyzer/Demo.tar.gz) files and decompress 'Demo.tar.gz' file

```
$ wget --no-check-certificate https://paoyang.ipmb.sinica.edu.tw/MethylC-analyzer/Demo.tar.gz
$ tar -xvf Demo.tar.gz
```
> Demo contains 4 CGmap, one Samples_list.txt and one GTF file
```
WT1s.CGmap.gz
WT2s.CGmap.gz
MT1s.CGmap.gz
MT2s.CGmap.gz
samples_list.txt
TAR10_2_demo.gtf

```
2. Make sure all scripts and all demo files are in the same foder

> e.g. move scripts in the 'MethylC-analyzer/Demo'
```
$ cd Demo
$ cp ../../scripts/* ./

```

1. A sample description file "samples_list.txt" provides in Demo. The file is ***tab-delimited*** without a header.

> Sample Description File (sample_name  CGmap_location  group (tab-delimited, no header in the first line)

```    
MT1     MT1s.CGmap.gz   MT
MT2     MT2s.CGmap.gz   MT
WT1     WT1s.CGmap.gz   WT
WT2     WT2s.CGmap.gz   WT
```

**Input:**
1. gene annotation (GTF)

2. samples_list.txt (description file, tab-delimited)

**Usage:**
```
$ python MethylC.py samples_list.txt TAR10_2_demo.gtf -a MT -b WT -d 4 -r 200 -q 2 -bs 200000 -pvalue 0.1

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
  -a GROUP1           Name of experimental group
  -b GROUP2           Name of control group
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



