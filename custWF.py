"""Methods to supplement the general workflow of genWF.py.
"""

import genWF
import numpy as np
import time
import sys
import os
import random

PWPATH = '/titan/cancerregulome8/TCGA/scripts/pairwise-2.0.0-current'
MAXQ = .9

def run2FMPWwList(FM1Path,FM2Path,pwOutPath,outDir,pwWhich=PWPATH,samples=[],maxQ=MAXQ,testFListPath=''):
	"""Run pairwise between two existing FMs.
	FM1 and FM2 considered test and target features, respectively.
	Procedure creates tmp files in outDir then removes them. 
	PW output saved to pwOutPath. 
	Output is filtered by FDR q values with a max set in maxQ.
	Uses PW code at pwWhich.
	samples is a str list of sample ids to test on,
	if empty uses list in FM1.	
	If specified a new line separated list of test features below 
	maxQ will be saved at testFListPath.
	"""
	FM1 = open(FM1Path)
	if len(samples)==0:
		# get sample names
		samples = FM1.next().strip().split('\t')[1:]
	else: FM1.next()

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
	genWF.catFM([FM1Path,FM2Path],samples,FMTmpPath,studyID='101',allowedSuffix=['FAM'])

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
	outDir = '/isb/rtasseff'
	testFMPath = '/isb/rtasseff/tmp.fm'
	targFMPath = '/isb/rtasseff/tmp.fm'
	fullOutPath = '/isb/rtasseff/fullOut.dat'
	listOutPath = '/isb/rtasseff/testList.dat'

	run2FMPWwList(testFMPath,targFMPath,fullOutPath,outDir,pwWhich=PWPATH,maxQ=MAXQ,testFListPath=listOutPath)		

if __name__ == '__main__':
	main()
