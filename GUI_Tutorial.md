How to run GUI
==============

1. Run the GUI interface by command line

```
  python gui.py
```

![Alt text](https://github.com/RitataLU/MethylC-analyzer/blob/master/GUI_step1.png)

2. Choose the location of input Files (sample list & Gene annotation gtf file)
  
  This section let users to select input files, containing sample list (.txt) and gene annotation gtf file (.gtf)

![Alt text](https://github.com/RitataLU/MethylC-analyzer/blob/master/GUI_step2.png)

  
Sample list format should be a txt file and must contains 3 columns, sample name, CGmap location and sample’s group
e.g. 
```
  WT_1	C:/Users/usrname/project/CGmap/WT_1.CGmap.gz	WT
  WT_2	C:/Users/usrname/project/CGmap/WT_2.CGmap.gz	WT
  Mut_1	C:/Users/usrname/project/CGmap/Mut_1.CGmap.gz	Mutant
  Mut_2	C:/Users/usrname/project/CGmap/Mut_2.CGmap.gz	Mutant
```
Selecting gene annotation file’s path (*gene annotation must be a gtf format) 


3. Setting Parameters 

This section let users to set analysis parameters

![Alt text](https://github.com/RitataLU/MethylC-analyzer/blob/master/GUI_step3.png)

  (1)	Each parameter has its defalt setting and users can click the textbox to modify.
  (2)	GUI system provides parameter descriptions, the interface will display it when the mouse move beside the box. Users also can read the parameters description file to help you setting it.
  (3) The GUI has fool-proof design. If the parameter is integer only, the program can’t let you input float number or other illegal values. It can help users to set parameters more easily.
 
 
4. Analysis selection

![Alt text](https://github.com/RitataLU/MethylC-analyzer/blob/master/GUI_step4.png)
This section can users choose analysis
  (1) Check the tool’s box. If users choose ChrView or Metaplot, then users should provide advanced options (comparison setting). 
  (2) Advanced options of ChrView and Metaplot are selecting comparison groups. User should decide to how to comparison between 2 groups (The group name is load from sample list)

![Alt text](https://github.com/RitataLU/MethylC-analyzer/blob/master/GUI_group.png)

  (3)  The Metaplot advanced option needs to assign a genomic feature (5UTR, 3UTR, CDS, intron, exon, promoter, genebody)which will analysis plot in Metaplot.
  
![Alt text](https://github.com/RitataLU/MethylC-analyzer/blob/master/GUI_genomic.png)
