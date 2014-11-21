inovaGMWF
========================
Genomic and Molecular WorkFlows for the Inova Project.

Dependencies
------------------
Python 2.7
numpy 1.6.1

Files
-------------------
####genWF.py
General WorkFlows for the basic analysis done on genomic and molecular data.
Currently merges, checks and filters feature matrices & runs a pairwise analysis 
on metadata.  Includes the creation of detailed log for provenance.
Call file (ie $python2.7 genWF.py) for usage, or use -h (ie $python2.7 genWF.py -h)
for more information.

This information is limited, see the [google doc](https://docs.google.com/a/systemsbiology.org/document/d/1mLPYANWA1IHjjzHw22UNAsatM2zOZ7xyLE7g21mTfS4/edit#)
for more info of the workflow.
Its private, so you need to request sharing.
	
####genWF.cfg	
The configuration file for all parameters in genWF.py.  Very flexable, 
can start workflow at any point, see comments in the file for more info.

####statsUtil.py
Basic statistic utility methods used in genWF.py, mostly related to visualizing the data 
for statistically significant feature pairs.

####genUtil.py
General utility methods for the INOVA project:
converting/identifying sample ID formats, 
saving feature matrices (FMs),
converting individual based FMs to family based FMs.

####gnmcUtil.py
Utility methods for dealing with genomic data (sequencing data).
Reading the source files (like the CG manifest).
Soon to add:
Reading VCF, collection of QC measures, handling of transcript data, running pairwise on VCF.

To Dos
---------------------
Will improve the documentation in this file soon.

