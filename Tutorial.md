# Tutorial

# Installation

1. Obtain Python 2.7 and virturalenv.

    MethylC-analyzer depends on [SAMtools](http://www.htslib.org/) and
    [BEDtools](http://bedtools.readthedocs.org/), so please make sure you
    already have them on your server.

    
2. Create a virtual environment somewhere on your disk, and then activate it.

  ```
  $ virtualenv --no-site-packages --python=python2.7 methylC_env
  $ source methylC_env/bin/activate

 ```
3. Download the source code and install the requirements.

  ```
  $ git clone https://github.com/RitataLU/MethylC-analyzer.git
  $ sudo sh MethylC-analyzer/requirements/base.txt
 ```

4. Add MethylC-anlyzer/script path to the PATH environment variable.
``` 
$ PATH=$PATH:(MethylC-analyzer/script file path)
$ source ~/.bash_profile
```

# Run demo
Download the demo input file in ATACgraph folder

```
$ cd MethylC-analyzer/demo
$ tar -xvf demo.tar.gz 
```

1.   Make a sample list and name it as "samples_list.txt" in the location where methylc.py script

samples list format:
sample_name  CGmap_location  group (seperate with a tab)
```    
wt1     ./wt1_demo.CGmap.gz  WT
wt2     ./wt2_demo.CGmap.gz  WT
met1_1       ./met1_1_demo.CGmap.gz    met1
met1_2       .met1_2_demo.CGmap.gz     met1
```

**Input:**
1. gene annotation
    
       gene annotation in GTF

2. CGmap

        CGmap files after mapping by utilizing Bsseeker2

**Usage:**
```
$ python MethylC.py samples_list.txt TAR10.genes.gtf


usage: methylC.py [-h] [-d DEPTH] [-r REGION] [-q QUAIFIED]
                  [-hcgc HEATMAP_CG_CUTOFF] [-hchhc HEATMAP_CHH_CUTOFF]
                  [-hchgc HEATMAP_CHG_CUTOFF] [-dcgc DMR_CG_CUTOFF]
                  [-dchhc DMR_CHH_CUTOFF] [-dchgc DMR_CHG_CUTOFF]
                  [-b BIN_SIZE] [-p PROMOTER_SIZE]
                  input_gtf_file sample_list

positional arguments:
  input_gtf_file        path of gene annotation
  sample_list           path of sample list

optional arguments:
  -h, --help            show this help message and exit
  -d DEPTH              min site of #C+#T
  -r REGION             size of region
  -q QUAIFIED           min number of region member
  -hcgc HEATMAP_CG_CUTOFF
                        Heatmap_CG_cutoff
  -hchhc HEATMAP_CHH_CUTOFF
                        Heatmap_CHH_cutoff
  -hchgc HEATMAP_CHG_CUTOFF
                        Heatmap_CHG_cutoff
  -dcgc DMR_CG_CUTOFF   DMR_CG_cutoff
  -dchhc DMR_CHH_CUTOFF
                        DMR_CHH_cutoff
  -dchgc DMR_CHG_CUTOFF
                        DMR_CHG_cutoff
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

**activate interface (choose analysis want to process)
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
**Output Figureq**


* Average methylation level

* Heatmap & PCA for variable regions 
* PCA for variable regions 
* Identifying Differentially Methylated Regions (DMRs)
* Genomic regions fold enrichment analysis for DMRs 
* Identifying Differentially Methylated Genes (DMGs)
* The distribution fo DNA methylation on each chromosome
* Metaplot for each profile & comparison between groups 



