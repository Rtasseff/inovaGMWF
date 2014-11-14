inovaGMWF
========================
Genomic and Molecular WorkFlows for the Inova Project.

Dependencies
------------------
Python 2.7
numpy 1.6.1

Files
-------------------
genWF.py
General WorkFlows for the basic analysis done on genomic and molecular data.
Currently merges, checks and filters feature matrices & runs a pairwise analysis 
on metadata.  Includes the creation of detailed log for provenance.
Call file (ie $python2.7 genWF.py) for usage, or use -h (ie $python2.7 genWF.py -h)
for more information.
	
genWF.cfg	
The configuration file for all parameters in genWF.py.  Very flexable, 
can start workflow at any point, see comments in the file for more info.

util.py
Basic utility methods used in genWF.py

To Dos
---------------------
Will improve the documentation in this file soon.

