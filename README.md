
# MethylC-analyzer

MethylC-analyzer is a analyzer developing for analysing WGBS and RRBS, it can utilize not only individual sample also do comparison between groups.
 
MethylC-analyzer will produce 7 analysis and each analysis contains CG, CHG and CHH 3 context:

![MethylC-analyzer Flow](https://github.com/RitataLU/Methylpip/blob/master/Main%20flow.jpg)

* Heatmap 
* PCA
* Differentially Methylated Regions (DMRs)
* DMRs Fold Enrichment (DMRs location enerichment)
* Differentially Methylated Genes (DMGs)
* Whole genome chromosome View for each profile & comparison between groups
* Metaplot for each profile & comparison between groups 

# system requirement 
* CPU：No special restrictions, but CPU has 16 cores is more efficient

* MEM：12GB or higer (for plant sample) / 256GB or higher (for human sample)

* Python 2.7
* [SAMtools](http://www.htslib.org/)
* [deepTools](https://deeptools.readthedocs.org/)
* [BEDtools](http://bedtools.readthedocs.org/)

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


> 1.   Make a sample list and name it as "sample_list.txt" in the location where methylc.py script

> > sample list format:
    sample name  CGmap location  group (seperate with a tab)
```    
      Contrl_1     ./contrl_1.CGmap  Control
      Contrl_2     ./contrl_2.CGmap  Control
      Meta_1       ./Meta_1.CGmap    Meta
      Meta_2       ./Meta_2.CGmap    Meta
```
> 2. * Usage:

```
$ python MethylC.py input_gtf_file


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
     
        gene annotation in gtf or gff

     2. CGmap

        CGmap files after mapping by utilizing Bsseeker2
  
  
  ## GUI interface
  
  The MethylC-analyzer also provide a user friendly GUI interface to let users who are not familiar programming.
  
  [GUI tutorial](https://github.com/RitataLU/Methylpip/blob/master/GUI_Tutorial.docx/  "GUI tutorial")
  
  
  
	

