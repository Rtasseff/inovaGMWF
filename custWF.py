"""Methods to supplement the general workflow of genWF.py.
"""

import genWF
import numpy as np
import time
import sys
import os
import random

PWPATH_TITAN = '/titan/cancerregulome8/TCGA/scripts/pairwise-2.0.0-current'
PWPATH_SGI = '/isb/rkramer/bin/pairwise-2.0.1'
MAXQ = .1

def run2FMPWwList(FM1Path,FM2Path,pwOutPath,outDir,pwWhich=PWPATH_SGI,samples='',maxQ=MAXQ,testFListPath=''):
	"""Run pairwise between two existing FMs.
	FM1 and FM2 considered test and target features, respectively.
	Procedure creates tmp files in outDir then removes them. 
	PW output saved to pwOutPath. 
	Output is filtered by FDR q values with a max set in maxQ.
	Uses PW code at pwWhich.
	If samples == '', use FM1 header row to define sample list,
	else if samples is a path, load a tsv sample list from that path,
	otherwise samples must be a str list containing the sample ids.
	If specified a new line separated list of test features below 
	maxQ will be saved at testFListPath.
	"""

	FM1 = open(FM1Path)
	if samples=='':
		#samples must be taken from the FM
		samples = FM1.next().strip().split('\t')[1:]
	else: 
		# need to advacne the FM since header not used
		FM1.next()
		# something was passed, what is it
		if type(samples)==str:
			# must be a path
			samples = np.loadtxt(samples,dtype=str,delimiter='\t')
		# if not a str, must be some kind of list or np array (otherwise an error will be throw later
		

	# get feature names for test list
	testF = [line.strip().split('\t')[0] for line in FM1]
	FM1.close()
	# get feature names for target list
	FM2 = open(FM2Path)
	FM2.next()
	targF = [line.strip().split('\t')[0] for line in FM2]
	FM2.close()

	# create pair list file for pw
	pairListPath = outDir+'/tmpPairList_'+str(random.randrange(16**5))+'.tsv'
	genWF._mkPairFile(testF,targF,pairListPath)

	# join FMs for pairwise (could do it faster, but this has mulitple checks
	FMTmpPath = outDir+'/tmpFM_'+str(random.randrange(16**5))+'.tsv'
	genWF.catFM([FM1Path,FM2Path],samples,FMTmpPath,checkSampIDs=False)

	genWF.runPairwise(FMTmpPath,outDir,pwOutPath,pwWhich=pwWhich,pairFile=pairListPath,fdrFilter=maxQ)

	# get list of test variables
	if testFListPath!='':
		pwOut = open(pwOutPath)
		testFList = [line.strip().split('\t')[0] for line in pwOut]
		pwOut.close()
		# remove repeats:
		testFList = list(set(testFList))
		testFListFile = open(testFListPath,'w')
		for value in testFList:
			# ignore comment lines
			if value[0]!='#':testFListFile.write(value+'\n')
		testFListFile.close()
	
	os.remove(pairListPath)
	os.remove(FMTmpPath)

	


def main():
	outDir = '/isb/rtasseff/results/var_batch_20150121'
	testFMPath = '/isb/rtasseff/data/featureMatrices/DF5_MergedVCF_ForPairwise_ExcludingPO.fm'
	targFMPath = '/isb/rtasseff/data/featureMatrices/BATCH_GNMC_IND_20150109.fm'
	fullOutPath = '/isb/rtasseff/results/var_batch_20150121/fullOut.dat'
	listOutPath = '/isb/rtasseff/results/var_batch_20150121/var_batchFilter_featureNameList_20150121.dat'
	samplePath = '/isb/rtasseff/data/support/sampleIDList_ind_DF5_itmiFormat.dat'

	run2FMPWwList(testFMPath,targFMPath,fullOutPath,outDir,maxQ=MAXQ,testFListPath=listOutPath,samples=samplePath)		

if __name__ == '__main__':
	main()
