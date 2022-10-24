
# MethylC-analyzer

MethylC-analyzer is a analyzer developing for analysing DNA methylation on WGBS and RRBS, it could utilize not only individual sample also do comparison between two groups.

MethylC-analyzer will produce 7 analysis and each analysis contains CG, CHG and CHH 3 context:

![MethylC-analyzer Flow](https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/Fig1_tmp.png)

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

###  From Docker 

:pushpin: Recommend using docker image to avoid enveioment conflict

```dockerfile
docker pull peiyulin/methylc
```

### From Github

1. Obtain Python 3.9
2. Recommend to create a [conda](https://docs.conda.io/en/latest/miniconda.html) environment somewhere on your disk, and then activate it.

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

# Input files


1. CGmap.gz file (need gzip compressed format) is the output of BS-Seeker2.(post-alignment data by utilizing [Bs-seeker2](https://github.com/BSSeeker/BSseeker2) and [Bs-seeker3](https://github.com/khuang28jhu/bs3))
> CGmap format
```
chr1    G       13538   CG      CG      0.6     6       10
chr1    G       13539   CHG     CC      0.0     0       9
chr1    G       13541   CHH     CA      0.0     0       9
chr1    G       13545   CHH     CA      0.0     0       8
```

> The methylation calling files from other aligners/callers, MethylC-analyzer provides a python script (methcalls2CGmap.py) to convert them to CGmap.gz, including CX report files generated by Bismark, the methylation calls generated by methratio.py in BSMAP (v2.73), and the TSV files exported from the methylation calling status with METHimpute.

```
usage: methcalls2CGmap.py [-h] [-n FILENAME]
                         [-f {bismark,bsmap,methimpute}]

optional arguments:
 -h, --help            show this help message and exit

Input format:
 -n FILENAME, --filename FILENAME
                       the file name that the users want to convert to CGMap
                       format
 -f {bismark,bsmap,methimpute}, --format {bismark,bsmap,methimpute}
                       the type of file to CGmap
```
>> Example for converting methylation calls to CGmap.gz:

```
# bismark to CGmap.gz
python methcalls2CGmap.py -n CX_report.txt.gz -f bismark

```

2.Gene annotation (GTF)

> gene annotation in GTF file: User can downloaded from [ensemble FTP](https://useast.ensembl.org/info/data/ftp/index.html)


# Tutorial 
Please follow the tutorial of example use case

[MethylC-analyzer docker tutorial](https://github.com/RitataLU/MethylC-analyzer/blob/master/Tutorial_Docker.md)  :mega:**Recommend**

[MethylC-analyzer command line tutorial](https://github.com/RitataLU/MethylC-analyzer/blob/master/Tutorial.md)



# Run MethylC-analyzer

Make a sample description file and name it as "samples_list.txt" in the location where methylc.py script. The file should be tab-delimited without a header.

> Sample Description File (***tab-delimited***, no header in the first line)
Sample list ( sample_name  CGmap_location  group )

```

wt1	wt1.CGmap.gz	WT
wt2	wt2.CGmap.gz	WT
wt3	wt3.CGmap.gz	WT
met1_1	met1_1.CGmap.gz	met1
met1_2	met1_2.CGmap.gz	met1
met1_3	met1_3.CGmap.gz	met1

```


**Usage:**
```
$ python MethylC.py samples_list.txt TAR10.genes.gtf 


usage: MethylC_new.py [-h] [-a GROUP1] [-b GROUP2] [-d DEPTH] [-r REGION]
                      [-q QUALIFIED] [-context CONTEXT] [-hc HEATMAP_CUTOFF]
                      [-dmrc DMR_CUTOFF] [-test TESTMETHOD] [-pvalue PVALUE]
                      [-bs BIN_SIZE] [-p PROMOTER_SIZE]
                      samples_list input_gtf_file

positional arguments:
  samples_list        samples CGmap description
  input_gtf_file      path of gene annotation

  

  
  ## Arguments
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
  

 ## activate interface (Users select analysis that want to process)

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



 


​	

