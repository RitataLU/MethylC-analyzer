
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


# Installation

1. Obtain Python 2.7 and virturalenv.

    MethylC-analyzer depends on `SAMtools <http://www.htslib.org/>`_ and
    `BEDtools <http://bedtools.readthedocs.org/>`_, so please make sure you
    already have them on your server.

    
2. Create a virtual environment somewhere on your disk, and then activate it.

  ```

  $ virtualenv --no-site-packages --python=python2.7 methylC_env
  $ source methylC_env/bin/activate

 ```
3. Download the source code and install the requirements.

  ```

  $ git clone https://github.com/RitataLU/MethylC-analyzer.git
  $ pip install -r MethylC-analyzer/requirements/base.txt
 ```

pip will install the following packages: 

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

# Tutorial 
Please follow the tutorial of example use case
* [Tutorial](https://github.com/RitataLU/ATACgraph/blob/master/Tutorial.md)

  
# GUI interface
  
The MethylC-analyzer also provides a user friendly GUI interface to let users who are not familiar programming.
please follow the tutorial   
[GUI tutorial](https://github.com/RitataLU/MethylC-analyzer/blob/master/GUI_Tutorial.md)

  
  
  
	

