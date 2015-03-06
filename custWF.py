"""Methods to supplement the general workflow of genWF.py.
"""

import genWF
import numpy as np
import time
import sys
import os
import random
import genUtil
import gzip
import statsUtil
import warnings


PWPATH_TITAN = '/titan/cancerregulome8/TCGA/scripts/pairwise-2.0.0-current'
PWPATH_SGI = '/isb/rkramer/bin/pairwise-2.0.1'

def run2FMPWwList(FM1Path,FM2Path,pwOutPath,outDir,pwWhich=PWPATH_SGI,samples='',maxQ=.1,testFListPath=''):
	"""Run pairwise between two existing FMs.
	FM1 and FM2 considered test and target features, respectively.
	Procedure creates tmp files in outDir then removes them. 
	PW output saved to pwOutPath. 
	If maxQ>0 output is filtered by FDR q values with a max set in maxQ.
	Uses PW code at pwWhich.
	If samples == '', use FM1 header row to define sample list,
	else if samples is a path, load a tsv sample list from that path,
	otherwise samples must be a str list containing the sample ids.
	If specified a new line separated list of test features below 
	maxQ will be saved at testFListPath.
	"""

	FM1 = genUtil.open2(FM1Path)

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
	FM2 = genUtil.open2(FM2Path)
	FM2.next()
	targF = [line.strip().split('\t')[0] for line in FM2]
	FM2.close()

	# create pair list file for pw
	pairListPath = outDir+'/tmpPairList_'+str(random.randrange(16**5))+'.tsv'
	genWF._mkPairFile(testF,targF,pairListPath)

	# join FMs for pairwise (could do it faster, but this has mulitple checks
	FMTmpPath = outDir+'/tmpFM_'+str(random.randrange(16**5))+'.tsv'
	genWF.catFM([FM1Path,FM2Path],samples,FMTmpPath,checkSampIDs=False)

	if maxQ < 0: maxQ = ''

	genWF.runPairwise(FMTmpPath,outDir,pwOutPath,pwWhich=pwWhich,pairFile=pairListPath,fdrFilter=maxQ)

	# get list of test variables
	if testFListPath!='':
		getSimpleList(pwOutPath,testFListPath)
		
	os.remove(pairListPath)
	os.remove(FMTmpPath)

def runPW_INDvFAM(FM1Path,FM2Path,pwOutPath,outDir,nbListPath,pwWhich=PWPATH_SGI,samples='',maxQ=''):
	"""Run pairwise between two existing FMs.
	FM1 and FM2 considered test and target features, respectively.
	Here we further assume that the test matrix is an IND-based FM
	and that FM2 is a FAM-based FM, we convert the former.
	Procedure creates tmp files in outDir then removes them. 
	PW output saved to pwOutPath. 
	If maxQ>0 output is filtered by FDR q values with a max set in maxQ.
	Uses PW code at pwWhich.
	If samples == '', use FM2 header row to define sample list,
	else if samples is a path, load a tsv sample list from that path,
	otherwise samples must be a str list containing the sample ids.
	These must be Family based IDs.

	This code has a very specific use, which is running a pairwise
	between an exiting FAM based target (phenotypes) against many small
	FMs that were gathered from a vcf.  Because of the large number of 
	test FMs it is not worth while to do a single IND 2 FAM conversion 
	for each and keep it.
	"""

	FM2 = genUtil.open2(FM2Path)

	if samples=='':
		#samples must be taken from the FM
		samples = FM2.next().strip().split('\t')[1:]
	else: 
		# need to advacne the FM since header not used
		FM2.next()
		# something was passed, what is it
		if type(samples)==str:
			# must be a path
			samples = np.loadtxt(samples,dtype=str,delimiter='\t')
		# if not a str, must be some kind of list or np array (otherwise an error will be throw later
		

	# get feature names for test list
	targF = [line.strip().split('\t')[0] for line in FM2]
	FM2.close()

	# convert the test matrix
	famFMPath = outDir+'/tmpFAM_'+str(random.randrange(16**5))+'.fm'
	nbList = genUtil.getRepNBList(nbListPath)
	genUtil.indFM2FamFM(FM1Path,famFMPath,nbList,samples)




	# get feature names for target list
	famFM = genUtil.open2(famFMPath)
	famFM.next()
	testF = [line.strip().split('\t')[0] for line in famFM]
	famFM.close()

	# create pair list file for pw
	pairListPath = outDir+'/tmpPairList_'+str(random.randrange(16**5))+'.tsv'
	genWF._mkPairFile(testF,targF,pairListPath)

	# join FMs for pairwise (could do it faster, but this has mulitple checks
	FMTmpPath = outDir+'/tmpFM_'+str(random.randrange(16**5))+'.fm'
	genWF.catFM([famFMPath,FM2Path],samples,FMTmpPath,checkSampIDs=False)

	if maxQ < 0: maxQ = ''

	genWF.runPairwise(FMTmpPath,outDir,pwOutPath,pwWhich=pwWhich,pairFile=pairListPath,fdrFilter=maxQ)
		
	os.remove(pairListPath)
	os.remove(FMTmpPath)
	os.remove(famFMPath)

	

def getSimpleList(pwOutPath,testFListPath):
	"""Make a simple non-redundant list
	of the test features (0th column of 
	the pairwise output at pwOutPath)
	and save it to testFListPath.
	"""
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


def _mkDir(outDir):
	if os.path.exists(outDir):
		warnings.warn("The dir "+outDir+"already exists, results may be overwritten.",UserWarning)	
	else:
		os.makedirs(outDir)

def _pw2VarTable(inPath,outPath):
	"""Take the standard pairwise output at inPath
	and put it into our 'results table' format 
	at outPath.
	"""
	pw = open(inPath)
	rt = open(outPath,'w')
	rt.write('chr\tpos\tgene\teffect\tp-value\n')
	for line in pw:
		tmp = line.strip().split('\t')
		fName = tmp[0].split(':')
		var = fName[-1].split('_')
		chro = var[0]
		pos = var[1]
		if fName[0]!='C':eff=tmp[3]
		else:eff='NA'
		gene='NA' # currently no easy way to grab this here
		p = 10**(-1*float(tmp[5]))
		pStr = '%05.4E' % (p)
		rt.write(chro+'\t'+pos+'\t'+gene+'\t'+eff+'\t'+pStr+'\n')
	pw.close()
	rt.close()



def parseResultTables(outDir,inputPaths,phenoCodes):
	"""Method will take paths to existing 
	result output and parse them into tables
	within a set directory structure under outDir.
	inputPaths is a directory that contains 
	all the paths to be parsed.
	These paths are actually the outputs of 
	several different tests.
	Possible keys:
	varPheno-FAM	Path to specific file
			variant pairwise all phenotypes for all family members
			results are initially in a single file to be sorted,
			assumes vars are 1st feature and pheno is second in 
			the standard pairwise output file.
	varBatch-FAM	Path to specific file
			variant pairwise for batch for all family members
			results are initially in a single file to be sorted,
			assumes vars are 1st feature and batch is second in 
			the standard pairwise output file.
	transPheno-FAM	Path to specific file
			transcript burden pairwise all phenotypes for all family members
			results are initially in a single file to be sorted
			assumes vars are 1st feature and pheno is second in 
			the standard pairwise output file
	transBatch-FAM	Path to specific file
			transcript burden pairwise for batch for all family members
			results are initially in a single file to be sorted
			assumes transcripts are 1st feature and batch is second in 
			the standard pairwise output file
	eigenstrat-M	Variable path, all results must be in same dir
			standard eigenstrat output for mother against all phenotypes
	eigenstrat-F	Variable path, all results must be in same dir
			standard eigenstrat output for father against all phenotypes
	eigenstrat-NB	Variable path, all results must be in same dir
			standard eigenstrat output for newborn against all phenotypes
	cifBat-FAM	Variable path, all results must be in same dir
			cifBath output for testing family trio on all phenotypes
	
	Variable path = used to simplify the path definitions, for these outputs
			phenotype targets are in separate files, but we can expect 
			them to be in the same dir with the same format for file
			name, but with different strings for each phenotype.
			Therefore use the string $PHENOCODE$ in place of the 
			phenotype string.  This code will search all phenotype
			strings in the list phenoCodes.
	If key is not found parsing for the corresponding result will be skipped.
	"""
	#NOTE: this is for readme documentation, probably better ways to do it
	lastEdit = '20150306'
	
	### setup
	
	# --set up the output dir
	 mkDir(outDir)
	
	# may not overwrite all tests, so just append info if DIR exits, let user sort it out
	readme = open(outDir+'/README.txt','a')

	readme.write('------------------------------------------------------------------------------')
	readme.write('Result Tables for Genomic Association Tests parsed on '+time.strftime("%c")+'.\n')
	readme.write('This version of parsing code last edited on '+lastEdit+'.\n')
	readme.write('First level dirs refer to individual tests, followed by family member.\n')
	readme.write('Tables are saved per target (phenotype) under the member directory.\n')

	
	### varPheno-FAM
	if inputPaths.has_key('varPheno-FAM'):
		path = inputPaths['varPheno-FAM']

		# --Create folder
		outDirTmp = outDir+'/varPheno'
		mkDir(outDirTmp)

		# over writing tests, so over write any readme.
		readmeTmp = open(tmpOutDir+'/README.txt','w')
		readmeTmp.write('--Parsing varPheno on '+time.strftime("%c")+'.\n')
		readmeTmp.write('Contains the pairwise test results for variant calls vs phenotypes.\n')
		readmeTmp.write('- source path:'+path+'\n')
		readmeTmp.write('- results dir:'+outDirTmp+'\n')

		# split the original file by members, 1 should set terms to members:
		membList = genWF.splitPWResults(path,outDirTmp,1)

		# gonna need those names later
		inFileName = os.path.splitext(os.path.basename(path))[0]

		# go through each family member 
		for memb in membList:
			path2 = outDirTmp+'/'+inFileName+'_subset_'+memb+'.dat'
			outDirTmp2 = outDirTmp+'/'+memb
			mkDir(outDirTmp2)
			# now split this even further into phenotypes
			phenoList = genWF.splitPWResults(path2,outDirTmp2,4,fName='target')
			# now parse each phenotype
			for pheno in phenoList:
				path3 = outDirTmp2+'/'+inFileName+'_subset_'+memb+'_subset_'+pheno+'.dat'
				_pw2VarTable(path3,outDirTmp2+'/resultTable_'+pheno+'.tsv')
				# remove this path, NOTE:comment if you want to keep split output
				os.remove(path3)
			os.remove(path2)

		# that should be each memb each pheno

			

			
		

		# put some notes in original readme
		readme.write('-varPheno DONE\n')
	else:
		readme.write('-varPheno SKIPPED\n')




	### varBatch:FAM

	### transPheno:FAM

	### transBatch:FAM

	### eigenstrat:M

	### eigenstrat:F

	### eigenstrat:NB

	### cifBat:FAM


def main():

#	###### run a single massive pairwise for getting the batch list -->
#	### NOTE:version 2.0.1 of pairwise will break over a certain size FM and did on DF5, 
#	### see below for alternate solutions RAT 20150202
#	outDir = '/isb/rtasseff/results/var_batch_20150128'
#	testFMPath = '/isb/rtasseff/data/featureMatrices/DF5_MergedVCF_ForPairwise_ExcludingPO.fm'
#	targFMPath = '/isb/rtasseff/data/featureMatrices/BATCH_GNMC_IND_20150109.fm'
#	fullOutPath = outDir+'/fullOut.dat'
#	listOutPath = outDir+'/var_batchFilter_featureNameList_20150128.dat'
#	samplePath = '/isb/rtasseff/data/support/sampleIDList_ind_DF5_itmiFormat.dat'
#
#	run2FMPWwList(testFMPath,targFMPath,fullOutPath,outDir,maxQ=.1,testFListPath=listOutPath,samples=samplePath)		
#	# <--
#
#	##### PRE run multiple small pairwise to get batch results, -->
#	# NOTE:designed to be called from PBS script for qsub -J
#	# in 20150203 on DF5 took < 5min per job
#	# after PBS run use the pw_region_postproc.sh shell script 
#	# and then the POST code segment below
#	tag = sys.argv[1]
#	outDir = '/isb/rtasseff/results/var_batch_20150204'
#	testFMPath = '/isb/rtasseff/data/featureMatrices/data_VCF_FM_regions_20150304/'+tag+'.Filtered.fm.gz'
#	targFMPath = '/isb/rtasseff/data/featureMatrices/BATCH_GNMC_IND_20150109.fm'
#	fullOutPath = outDir+'/fullPWOut_'+tag+'.dat'
#	samplePath = '/isb/rtasseff/data/support/sampleIDList_ind_DF5_itmiFormat.dat'
#	
#	run2FMPWwList(testFMPath,targFMPath,fullOutPath,outDir,maxQ=-1,samples=samplePath)		
#	# <--
	
#	### NOTE: must run pw_region_postproc.sh between PRE and POST!!!
#
#	##### POST run multiple small pairwise results
#	# assumes you did the PRE step above, and then the pw_region_postproc.sh,
#	# now you can run the last part to get FDR values and a unique list
#	# run from command line on ITMI SGI login took <30min
#	outDir = '/isb/rtasseff/results/var_batch_20150204'
#	inSortedPWPath =  outDir+'/fullPWOut_sorted.dat'
#	outFDRFilteredPWPath =  outDir+'/sigFDRPWOut.dat'
#	outBatchFList = outDir+'/fList_var_batch_20150204.dat'
#
#	statsUtil.fdr_bh_filterSortFile(inSortedPWPath,outFDRFilteredPWPath,alpha=.1,col=5,logTrans=True)
#	getSimpleList(outFDRFilteredPWPath,outBatchFList)


	##### PRE run multiple small pairwise to get phenotype results, -->
	# NOTE:designed to be called from PBS script for qsub -J
	# in 20150210 on DF5 took < 5min per job
	# after PBS run use the pw_region_postproc.sh shell script 
	# and then the POST code segment '##### POST run multiple small pairwise results'
	tag = sys.argv[1]
	outDir = '/isb/rtasseff/results/var_pheno_20150210'
	testFMPath = '/isb/rtasseff/data/featureMatrices/data_VCF_FM_regions_20150304/'+tag+'.Filtered.fm.gz'
	targFMPath = '/isb/rtasseff/data/featureMatrices/data_CLIN_Critical_Phenotype_20150203.fm'
	fullOutPath = outDir+'/fullPWOut_'+tag+'.dat'
	nbListPath = '/isb/rtasseff/data/support/repNBList.tsv'
	
	runPW_INDvFAM(testFMPath,targFMPath,fullOutPath,outDir,nbListPath)		
	# <--




if __name__ == '__main__':
	main()
