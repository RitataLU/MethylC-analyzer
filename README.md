
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

* Python 2.7
* [SAMtools](http://www.htslib.org/)
* [deepTools](https://deeptools.readthedocs.org/)
* [BEDtools](http://bedtools.readthedocs.org/)
* R package
* ComplexHeatmap
* ggplot2

# Tutorial 
Please follow the tutorial of example use case
* [Tutorial](https://github.com/RitataLU/ATACgraph/blob/master/Tutorial.md)


# Installation
1. Obtain Python 2.7 and virturalenv

MethylC-analyzer depends on SAMtools and BEDtools, so please make sure you already have them on your server.

3. Create a virtual environment somewhere on your disk, and then activate it.
4. Download the source code and install the requirements.


```
$ git clone https://github.com/RitataLU/ATACgraph.git
$ cd ATACgraph
$ sudo sh ./base.txt

``` 
Python Modules 'Numpy', 'pandas', 'Metplotlib' and pyBigWig. To install the packages, use the following commands on an UNIX terminal:

```    
        $ pip install numpy
  	$ pip install pandas
  	$ pip install math
  	$ pip install scipy
  	$ pip install matplolib
  	$ pip install argparse
  	$ pip install glob
  	$ pip install pyBigWig
  	$ pip install collections
  	$ pip install gzip
  	$ pip install re
  	$ pip install PyQt5
```  


# Running MethylC-analyzer


> 1.   Make a sample list and name it as "samples_list.txt" in the location where methylc.py script

> > samples list format:
    sample name  CGmap location  group (seperate with a tab)
```    
      wt_1     ./wt_1.CGmap.gz	WT
      wt_2     ./wt_2.CGmap.gz 	WT
      met1_1       ./met1_1.CGmap.gz	met1
      met1_2       .met1_2.CGmap.gz	met
```
> 2. * Usage:

```
$ python MethylC.py samples_list.txt input_gtf_file


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
    
   
    
  ## Input

     1. gene annotation
     
        gene annotation in GTF

     2. CGmap

        CGmap files after mapping by utilizing Bsseeker2
  
  
  ## GUI interface
  
  The MethylC-analyzer also provide a user friendly GUI interface to let users who are not familiar programming.
  
  [GUI tutorial](https://github.com/RitataLU/MethylC-analyzer/blob/master/GUI_Tutorial.md)

  
  
  
	

