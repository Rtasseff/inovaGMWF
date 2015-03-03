inovaGMWF
========================
Genomic and Molecular WorkFlows for the Inova Project.  
Currently, molecular workflows have been reduced in priority and metadata analysis has increased.

Dependencies
------------------
Python 2.7 \\
numpy 1.6.1 \\
ISB internal software: pairwise-2.0.0 (2.0.1 for full functionality)\\

Files
-------------------
####genWF.py
General WorkFlows for the basic analysis done on genomic and molecular data.
Currently formats (clinical), merges, checks and filters feature matrices; runs a pairwise analysis 
on metadata, runs pairwise on transcript data (and other sources generally) 
and creates summaries and plots all on a single workflow that is directed via a config file. 
Includes the creation of detailed log for provenance.
Call file (ie $python2.7 genWF.py) for usage, or use -h (ie $python2.7 genWF.py -h)
for more information.\\
typical call example:\\
python genWF.py genWF.cfg

This information is limited, see the [google doc](https://docs.google.com/a/systemsbiology.org/document/d/1mLPYANWA1IHjjzHw22UNAsatM2zOZ7xyLE7g21mTfS4/edit#)
for more info of the workflow.
It's private, so you need to request sharing.
	
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

####custWF.py
Truthfully, this is a bit of a workaround.  
We are using this file to house small workflows, or procedures, as methods
that we do not wish to incorporate with genWF.
This could be because they are prototypes, they are run 
in parallel, perhaps on different servers/systems (ie SGI), or too much time would 
be required to merge into genWF.

While it is a bit sloppy, and the whole structure of this code may need 
to be revisited, this allows us to move forward with important
analysis tasks without spending resources on redesigning the current code.

The main() method has several code snippets to run various procedures.
The comments should make it clear what they are and where they start/stop.
The easiest thing is to just comment out all but the needed procedure, 
and change the corresponding hard-coded paths/options.
We have tried to make it flow from top to bottom in the order of 
the typical workflow, but there are multiple 'in-between' steps 
(e.g. PBS job submissions or bash scripts) that prevented us from 
forming a single push-button code.

####bashShellScripts
This directory houses useful bash shell scripts that were used during 
the typical workflow.

pw\_region\_postproc.sh - offers some post processing of multiple small pairwise runs into one large one and sorts for FDR.  Used after calls like custWF::runPW\_INDvFAM() and before statsUtil::fdr\_bh\_filterSortFile()

runEigenstrat.pbs - bash code snippets to run eigenstrat analysis, last part (longest part) can be run via PBS job arrays.

####extras
Folder of extra code snippets, even less orered that custWF.py. 


