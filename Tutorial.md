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


1. Make a sample description file and name it as "samples_list.txt" in the location where methylc.py script. The file should be tab-delimited without a header.

> Sample Description File (tab-delimited, no header in the first line)
Sample list ( sample_name  CGmap_location  group )


**samples list format:**
>> sample_name  CGmap_location  group (tab-delimited, no header in the first line)
```    
wt1 wt1_demo.CGmap.gz   WT
wt2 wt2_demo.CGmap.gz   WT
met1_1  met1_1_demo.CGmap.gz    met1
met1_2  met1_2_demo.CGmap.gz    met1
```

**Input:**
1. gene annotation (GTF)

2. CGmap (post-alignment data by utilizing Bsseeker2)

**Usage:**
```
$ python MethylC.py samples_list.txt TAR10_2_demo.gtf -a WT -b MT -d 4 -r 200 -q 2 -bs 200000 -pvalue 0.1

usage: MethylC.py [-h] [-d DEPTH] [-r REGION] [-q QUALIFIED] [-hcgc HEATMAP_CG_CUTOFF] [-hchgc HEATMAP_CHG_CUTOFF]
                  [-hchhc HEATMAP_CHH_CUTOFF] [-dmrcg DMR_CG_CUTOFF] [-dmrchg DMR_CHG_CUTOFF] [-dmrchh DMR_CHH_CUTOFF] [-pvalue PVALUE]
                  [-b BIN_SIZE] [-p PROMOTER_SIZE]
                  samples_list input_gtf_file

positional arguments:
  samples_list          samples CGmap description
  input_gtf_file        path of gene annotation

optional arguments:

    -d, --DEPTH <INT> 
    minimum sites of methlated cytosine and unmethylated cytosine, default is 4
    
    -r, --REGION <INT>  
    size of region, default is 500 bp
    
    -q, --QUAIFIED <INT>
    qulified sites within a region, default is 4
    
    -hcgc, --HEATMAP_CG_CUTOFF <INT>
    PCA and Heatmap cutoff:
    CG Methylation difference between maximum and minimum regions , default is 0.2
    
    -hchgc, --HEATMAP_CHG_CUTOFF<INT>
    PCA and Heatmap cutoff:
    CHG Methylation difference between maximum and minimum regions , default is 0.2
    
    -hchhc, --HEATMAP_CHH_CUTOFF <INT>
    PCA and Heatmap cutoff:
    CHH Methylation difference between maximum and minimum regions , default is 0.2

    -dcgc, --DMR_CG_CUTOFF <INT>
    CG Methylation difference between 2 groups , default is 0.1
  
    -dchgc, --DMR_CHG_CUTOFF <INT>
    CHG Methylation difference between 2 groups , default is 0.1
    
    -dchhc, --DMR_CHH_CUTOFF <INT>
    CHH Methylation difference between 2 groups , default is 0.1
                        
    -b, --BIN_SIZE <INT>
    Cutoff of chrView and Metaplot:
    Seperate genome into several bins, and Size of bin, default is 1000000 bp
    
    -pvalue criteria for identfying DMR default is 0.05
    
    -p, --promoter <INT>
    Size of promoter, default is 2,000 bp before transcription start site
    ```

** activate interface (choose analysis want to process)
```
Heatmap & PCA Analysis?  (y/n): y
Identify DMR?  (y/n): y
Identify DMG?  (y/n): y
Use Fold Enrichment Analysis?  (y/n): y
Chromosome View Analysis?  (y/n): y
Metaplot Analysis?  (y/n): y
enter experimental group name analysis: met1
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



