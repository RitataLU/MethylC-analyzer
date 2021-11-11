
# MethylC-analyzer

MethylC-analyzer is a analyzer developing for analysing DNA methylation on WGBS and RRBS, it could utilize not only individual sample also do comparison between two groups.
 
MethylC-analyzer will produce 7 analysis and each analysis contains CG, CHG and CHH 3 context:

![MethylC-analyzer Flow](https://github.com/RitataLU/MethylC-analyzer/blob/master/MethylC-analyzer_main.png)

* Average methylation level
* Heatmap for variable regions 
* PCA for variable regions 
* Identifying Differentially Methylated Regions (DMRs)
* Genomic regions fold enrichment analysis for DMRs 
* Identifying Differentially Methylated Genes (DMGs)
* The distribution fo DNA methylation on each chromosome
* Metaplot for each profile & comparison between groups 

# <a name="SystemRequirements"></a>System Requirements
* Linux or Mac OS Environment
* CPU：No special restrictions, but CPU has 16 cores is more efficient

* MEM：12GB or higer (for plant sample) / 256GB or higher (for human sample)

* Python 3.9
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
  viridus
 ```


# Installation

1. Obtain Python 3.9
    
2. Recommand to create a [conda](https://docs.conda.io/en/latest/miniconda.html) environment somewhere on your disk, and then activate it.
  
  ```
  $ conda create -n methylC_analyzer_env python=3.9
  $ conda activate methylC_analyzer_env

 ```
3. Download the source code and install the requirements.

  ```
  $ git clone https://github.com/RitataLU/MethylC-analyzer.git
  
 ```
4. Install Package - Run MethylC-analyzer/requirements/base.txt

    ex: sh MethylC-analyzer/requirements/base.txt



# Tutorial 
Please follow the tutorial of example use case
[Tutorial](https://github.com/RitataLU/MethylC-analyzer/blob/master/Tutorial.md)

# Run MethylC-analyzer

1.   Make a sample list and name it as "samples_list.txt" in the location where methylc.py script

samples list format:
    sample_name  CGmap_location  group (seperate with a tab)
```    
wt1     wt1.CGmap.gz  WT
wt2     wt2.CGmap.gz  WT
wt3     wt3.CGmap.gz  WT
met1_1       met1_1.CGmap.gz    met1
met1_2       met1_2.CGmap.gz     met1
met1_3       met1_3.CGmap.gz    met1

```

**Input:**
1. gene annotation (GTF)

2. CGmap (post-alignment data by utilizing [Bs-seeker2] (https://github.com/BSSeeker/BSseeker2) and [Bs-seeker3] (https://github.com/khuang28jhu/bs3))

**Usage:**
```
$ python MethylC.py samples_list.txt TAR10.genes.gtf 


usage: MethylC.py [-h] [-d DEPTH] [-r REGION] [-q QUALIFIED] [-hcgc HEATMAP_CG_CUTOFF] [-hchgc HEATMAP_CHG_CUTOFF]
                  [-hchhc HEATMAP_CHH_CUTOFF] [-dmrcg DMR_CG_CUTOFF] [-dmrchg DMR_CHG_CUTOFF] [-dmrchh DMR_CHH_CUTOFF] [-pvalue PVALUE]
                  [-b BIN_SIZE] [-p PROMOTER_SIZE]
                  samples_list input_gtf_file

positional arguments:
  samples_list          samples CGmap description
  input_gtf_file        path of gene annotation
  
  ```
  
  ## Arguments
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

** activate interface (Users select analysis that want to process)
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

1. The average methylation in 3 context (CG, CHG, CHH)

 <img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/Average_methylation_levels.png" width="400">


2. PCA & Heatmap show variable region among samples

PCA: 

 <img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/PCA_CG_0.5.png" width="400">
 
Heatmap: 

 <img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/Heatmap_CG_0.5.png" width="400">


3. The distribution fo DNA methylation on each chromosome

<img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/chrView_CG.png" width="1000">

4. The distribution fo DNA methylation difference on each chromosome
<img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/chrView_delta_CG.png" width="1000">


5. Summary of dentifying Differentially Methylated Regions (DMRs) & Differentially Methylated Genes (DMGs)

<img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/Summary_DMR_DMG_numbers_CG_0.2.png" width="400">

   
6. Genomic regions fold enrichment analysis for DMRs 

<img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/CG_Fold_Enrichment.png" width="400">

7. The distribution of DNA methylation around gene body

<img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/metaplot_CG.png" width="400">

8. The distribution of DNA methylation difference around gene body
<img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/metaplot_delta_CG.png" width="400">



 
  
  
	

