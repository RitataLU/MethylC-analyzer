# Docker for MethylC-analyzer



## Installation

pull the official docker image by running:

```dockerfile
docker pull peiyulin/methylc
```



## Launch MethylC-Analyze container 

To launch the MethylC-Analyzer Docker container, the folder has to be mounted from the host system into the container with the -v argument. The examples mount the current working directory under /app inside the container and then run the command:



```bash
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py -h

usage: MethylC.py [-h] [-a GROUP1] [-b GROUP2] [-d DEPTH] [-r REGION]
                  [-q QUALIFIED] [-context CONTEXT] [-hc HEATMAP_CUTOFF]
                  [-dmrc DMR_CUTOFF] [-test TESTMETHOD] [-pvalue PVALUE]
                  [-bs BIN_SIZE] [-p PROMOTER_SIZE]
                  command samples_list input_gtf_file path_to_files
```

for example:

```dockerfile
cd path_to_data/

docker run --rm -v $(pwd):/app peiyulin/methylc:0.1:latest python /MethylC-analyzer/scripts/MethylC.py command /app/input /app/output /app/ 
```



## Examples

Download the demo input files:

```bash
wget --no-check-certificate https://paoyang.ipmb.sinica.edu.tw/MethylC-analyzer/Demo.tar.gz
tar -xvf Demo.tar.gz
cd Demo
```



### Run it all in one

Run all the functions in MethylC-analyzer in one command by using `all`:

```dockerfile
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py all samples_list.txt TAR10_2_demo.gtf /app/ -a MT -b WT -d 4 -r 200 -q 2 -bs 200000 -pvalue 0.1
```

**Input**

* Samples CGmap files

* A GTF file

* A sample list txt file for sample description 

  * The file is ***tab-delimited*** without a header

    ```bash
    ### samples_list.txt in Demo
    MT1     MT1s.CGmap.gz   MT
    MT2     MT2s.CGmap.gz   MT
    WT1     WT1s.CGmap.gz   WT
    WT2     WT2s.CGmap.gz   WT
    ```

    > Format descriptions:
    >
    > (1) sample_name (2) CGmap_location (3) group



### Run it separately 

#### Heatmap & PCA Analysis

```dockerfile
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py Heatmap_PCA samples_list.txt TAR10_2_demo.gtf /app/ -a MT -b WT -d 4 -r 200 -q 2 -bs 200000 -pvalue 0.1
```

**Output**

* Text files

  * Common region of each context: CommonRegion_CG.txt, CommonRegion_CHG.txt, CommonRegion_CHH.txt

  * Methylation union site: Unionsite.txt

  * Log of plotting: plot.log


* Bed files

  * Merge bed files: gtf_introns_merge.bed, gtf_3utr_merge.bed,  gtf_5utr_merge.bed, gtf_Promoter_merge.bed, gtf_exons_merge.bed, gtf_cds_merge.bed, gtf_Genebody_merge.bed

  * bed6 files: gtf_Promoter_bed6.bed, gtf_Genebody_bed6.bed, .gtf_IGR_bed6.bed


* Figures

  * Heatmap (leaf): Heatmap_CG_0.2.pd

  * PCA (middle): PCA_CG_0.2.pdf

  * The average methylation for each context (right): Average_methylation_levels.pdf

<p align="center"><img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/fig1.png" width="280"><img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/fig2.png" width="280"><img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/fig3.png" width="280"></p>



#### Identify DMR

```dockerfile
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py DMR samples_list.txt TAR10_2_demo.gtf /app/ -a MT -b WT -d 4 -r 200 -q 2 -bs 200000 -pvalue 0.1
```

**Output**

* Text files
  * DMR for CG context: DMR_CG_all_0.1.txt, DMR_CG_hyper_0.1.txt, DMR_CG_hypo_0.1.txt             



#### Identify DMG

```dockerfile
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py DMG samples_list.txt TAR10_2_demo.gtf /app/ -a MT -b WT -d 4 -r 200 -q 2 -bs 200000 -pvalue 0.1
```

**Output**

* Text files

  * Hyper DMG list: DMG_CG_hyper_0.1_Genebody_list.txt, DMG_CG_hyper_0.1_Promoter_list.txt

  * Hypo DMG list: DMG_CG_hypo_0.1_Genebody_list.txt, DMG_CG_hypo_0.1_Promoter_list.txt

* Bed files
  * bed for DMR: DMR_CG_all_0.1.txt.bed, DMR_CG_hypo_0.1.txt.bed, DMR_CG_hyper_0.1.txt.bed

* Figure
  * Bar chart for DMG: Summary_DMR_DMG_numbers_CG_0.1.pdf

<p align="center"><img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/fig4.png" width="400"></p>



#### Fold Enrichment Analysis

```dockerfile
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py Fold_Enrichment samples_list.txt TAR10_2_demo.gtf /app/ -a MT -b WT -d 4 -r 200 -q 2 -bs 200000 -pvalue 0.1
```

* Bed files

  * for common regions: CommonRegion_CG.txt.bed

  * bed.bed files

* Figures
  * Fold_Enrichment plot: CG_Fold_Enrichment.pdf

<p align="center"><img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/fig5.png" width="400"></p>



#### Chromosome View Analysis

```dockerfile
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py ChrView samples_list.txt TAR10_2_demo.gtf /app/ -a MT -b WT -d 4 -r 200 -q 2 -bs 200000 -pvalue 0.1
```

**Output**

* Text files
  * Chromosome view of each sample: MT1_200000_chrView.txt, MT2_200000_chrView.txt, WT1_200000_chrView.txt, WT2_200000_chrView.txt
  * Chromosome view list: chrView_delta_list.txt, chrView_delta.txt, chrView_list.txt 

* Figures
  * Chromosome view of each context (left): chrView_CG.pdf, chrView_CHG.pdf, chrView_CHH.pdf   
  * Delta chromosome view of each context (right): chrView_delta_CG.pdf, chrView_delta_CHG.pdf, chrView_delta_CHH.pdf

<p align="center"><img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/fig6.png" width="400"><img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/fig7.png" width="400"></p>

#### Metaplot Analysis

```dockerfile
docker run --rm -v $(pwd):/app peiyulin/methylc:V1.0 python /MethylC-analyzer/scripts/MethylC.py Metaplot samples_list.txt TAR10_2_demo.gtf /app/ -a MT -b WT -d 4 -r 200 -q 2 -bs 200000 -pvalue 0.1
```

**Output**

* text files
  * metaplot_delta_CG.txt, metaplot_delta_CHG.txt, metaplot_delta_CHH.txt

* bigwig files
  * For each sample and each context: MT1_CG.bw, MT1_CHG.bw, MT1_CHH.bw, MT2_CG.bw, MT2_CHG.bw, MT2_CHH.bw, WT1_CG.bw, WT1_CHG.bw, WT1_CHH.bw, WT2_CG.bw, WT2_CHG.bw, WT2_CHH.bw
* Matrix
  * MT1_CHH.matrix.gz, MT2_CHH.matrix.gz, MT2_CHH.matrix.gz, MT2_CHH.matrix.gz, 

* Figures
  * Meta gene plots for each context: metaplot_CG.pdf, metaplot_CHG.pdf, metaplot_CHH.pdf
  * Meta gene plots for each context: metaplot_delta_CG.pdf, metaplot_delta_CHG.pdf, metaplot_delta_CHH.pdf

<p align="center"><img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/fig8.png" width="400"><img src="https://github.com/RitataLU/MethylC-analyzer/blob/master/Figures/fig9.png" width="400"></p>
