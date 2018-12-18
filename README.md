
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
CPU：No special restrictions, but CPU has 16 cores is more efficient

MEM：12GB or higer (for plant sample) / 256GB or higher (for human sample)

Python 2.7 or above

 `[`SAMtools`]` `(`http://www.htslib.org/`)`
 
 [deepTools] (https://deeptools.readthedocs.org/"link")

 [BEDtools] (http://bedtools.readthedocs.org/"link")


Python Modules 'Numpy', 'pandas' and 'Metplotlib'. To install the packages, use the following commands on an UNIX terminal:
  
  `$ pip install numpy
  `$ pip install pandas`
  `$ pip install matplolib`
