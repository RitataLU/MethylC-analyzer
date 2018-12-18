
# MethylC-analyzer

MethylC-analyzer is a analyzer developing for analysing WGBS and RRBS, it can utilize not only individual sample also do comparison between groups.
 
MethylC-analyzer will produce 7 analysis and each analysis contains CG, CHG and CHH 3 context:
* Heatmap 
* PCA
* Differentially Methylated Regions (DMRs)
* DMRs Fold Enrichment (DMRs location enerichment)
* Differentially Methylated Genes (DMGs)
* Whole genome chromosome View for each profie & comparison between groups
* Metaplot for each profie & comparison between groups 

# system requirement 
* CPU：No special restrictions, but CPU has 16 cores is more efficient

* MEM：12GB or higer (for plant sample) / 256GB or higher (for human sample)

* Python 2.7
* [SAMtools](http://www.htslib.org/)
* [deepTools](https://deeptools.readthedocs.org/)
* [BEDtools](http://bedtools.readthedocs.org/)

Python Modules 'Numpy', 'pandas', 'Metplotlib' and pyBigWig. To install the packages, use the following commands on an UNIX terminal:

	<pip install numpy>
	<pip install pandas>
	<pip install matplolib>
	<pip install pyBigWig>

# Running MethylC

* Usage:

> 1.   Make a sample list first
> > sample list format:
```
      sample name\tCGmap location
      Contrl_1\tcontrl_1.CGmap 
      Contrl_2\tcontrl_2.CGmap 
      Meta_1\tMeta_1.CGmap 
      Meta_2\tMeta_2.CGmap 
```

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
  -d DEPTH              min size of #C+#T
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
  
  
