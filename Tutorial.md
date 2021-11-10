# Tutorial

# Installation

1. Obtain Python 3
    
2. Recommand to create a [conda](https://docs.conda.io/en/latest/miniconda.html) environment somewhere on your disk, and then activate it.
  
  ```
  $ conda create -n methylC_analyzer_env 
  $ conda activate (Yout conda path) methylC_analyzer_env

 ```
3. Download the source code and install the requirements.

  ```
  $ git clone https://github.com/RitataLU/MethylC-analyzer.git
  
 ```
4. Install Package - Run MethylC-analyzer/requirements/base.txt

    ex: sudo sh MethylC-analyzer/requirements/base.txt



# Run demo 



```
$ cd ../demo
$ tar -xvf demo.gz 
```

```
$ mv ../scripts/* ./

```

1.   Make a sample list and name it as "samples_list.txt" in the location where MethylC.py scripts

samples list format:
    sample_name  CGmap_location  group (seperate with a tab)
```    
wt1 ./wt1_demo.CGmap.gz WT
wt2 ./wt2_demo.CGmap.gz WT
met1_1  ./met1_1_demo.CGmap.gz met1
met1_2  ./met1_2_demo.CGmap.gz  met1
```

**Input:**
1. gene annotation (GTF)

2. CGmap (post-alignment data by utilizing Bsseeker2)

**Usage:**
```
$ python MethylC.py samples_list.txt TAR10_demo.gtf -d 2 -r 10 -q 2 -hcgc 0.001 -hchhc 0.001 -hchgc 0.001 -b 2000 -pvalue 0.1


usage: MethylC.py [-h] [-d DEPTH] [-r REGION] [-q QUALIFIED] [-hcgc HEATMAP_CG_CUTOFF] [-hchgc HEATMAP_CHG_CUTOFF]
                  [-hchhc HEATMAP_CHH_CUTOFF] [-dmrcg DMR_CG_CUTOFF] [-dmrchg DMR_CHG_CUTOFF] [-dmrchh DMR_CHH_CUTOFF] [-pvalue PVALUE]
                  [-b BIN_SIZE] [-p PROMOTER_SIZE]
                  samples_list input_gtf_file

positional arguments:
  samples_list          samples CGmap description
  input_gtf_file        path of gene annotation

optional arguments:
  -h, --help            show this help message and exit
  -d DEPTH              min site of #C+#T
  -r REGION             size of region
  -q QUALIFIED          qualified site within a region
  -hcgc HEATMAP_CG_CUTOFF
                        PCA & Heatmap_CG_cutoff
  -hchgc HEATMAP_CHG_CUTOFF
                        PCA & Heatmap_CHG_cutoff
  -hchhc HEATMAP_CHH_CUTOFF
                        PCA & Heatmap_CHH_cutoff
  -dmrcg DMR_CG_CUTOFF  DMR_CG_cutoff
  -dmrchg DMR_CHG_CUTOFF
                        DMR_CHG_cutoff
  -dmrchh DMR_CHH_CUTOFF
                        DMR_CHH_cutoff
  -pvalue PVALUE        p-value for identifying DMR
  -b BIN_SIZE           resolution of chrView and Metaplot
  -p PROMOTER_SIZE      promoter_size
  
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
    CG Methylation difference between 2 groups , default is 0.2
  
    -dchgc, --DMR_CHG_CUTOFF <INT>
    CHG Methylation difference between 2 groups , default is 0.2
    
    -dchhc, --DMR_CHH_CUTOFF <INT>
    CHH Methylation difference between 2 groups , default is 0.2
                        
    -b, --BIN_SIZE <INT>
    Cutoff of chrView and Metaplot:
    Seperate genome into several bins, and Size of bin, default is 1000000 bp
    
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



