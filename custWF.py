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
import gnmcUtil
import statsUtil


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
	Specifically reading for variants, expecting 
	naming convention fName last field (description)
	to be <chr>_<position>
	"""
	pw = genUtil.open2(inPath)
	rt = open(outPath,'w')
	rt.write('chr\tpos\tgene\teffect\tp-value\n')
	for line in pw:
		if line[0]!='#':
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

def _pw2TransTable(inPath,outPath,regionDic={}):
	"""Take the standard pairwise output at inPath
	and put it into our 'results table' format 
	at outPath.
	Specifically reading for transcript burdens, expecting 
	naming convention fName last field (description)
	to be <burden measure>_<gene>_<transcript ID>
	Additional information in regionDic,
	which hold gnmcUtil.Region objects 
	with region names as the key,
	expecting region name to be <gene>_<transcript ID>,
	as used in split-transcript code.
	"""

	pw = genUtil.open2(inPath)
	rt = open(outPath,'w')
	rt.write('transID\tgene\tchr\tstartPos\tstopPos\tnVar\teffect\tp-value\n')

	for line in pw:
		if line[0]!='#':
			try:
				# NOTE: issue with some transcript names containing ':', 
				# skipping for now, but should really fix that.
				tmp = line.strip().split('\t')
				fName = tmp[0].split(':')
				trans = fName[-1].split('_')
				transID = '_'.join(trans[2:])
				gene = trans[1]
				# original region names were appended to a descriptor in the feature name (MiAC)
				# NOTE: this all requires that the original naming scheme was kept
				# in the future it would be better to standardize the meta info in a DB
				# beyond using our simple multi field (: sep) feature naming convention.
				regionName = '_'.join(trans[1:])
				if fName[0]!='C':eff=tmp[3]
				else:eff='NA'
				p = 10**(-1*float(tmp[5]))
				pStr = '%05.4E' % (p)
				if len(regionDic)==0:
					chro = 'NA'
					start = 'NA'
					stop = 'NA'
					nVar = 'NA'
				else:
					region = regionDic[regionName]
					chro = region.chrom
					start = str(region.startPos)
					stop = str(region.stopPos)
					nVar = str(region.nVar)

				rt.write(transID+'\t'+gene+'\t'+chro+'\t'+start+'\t'+stop+'\t'+nVar+'\t'+eff+'\t'+pStr+'\n')
			except:
				print 'issue with line: '+line+' in '+inPath
	pw.close()
	rt.close()

def _es2VarTable(inPath,outPath):
	"""Take the standard eigenstrat output at inPath
	and put it into our 'results table' format 
	at outPath.
	
	"""
	#NOTE: would be more flexible/consistent in long term to combine all
	# *2VarTable methods into one method, limiting test specific code to different file parsers;
	# however, in the interest of time and uncertainty of reuse, we separated them out for now.
	fin = genUtil.open2(inPath)
	rt = open(outPath,'w')
	rt.write('chr\tpos\tgene\teffect\tp-value\n')
	# first line is comment:
	fin.next()
	for line in fin:
		try:
			tmp = line.strip().split()
			chro = tmp[1]
			pos = tmp[2]
			eff=tmp[4]
			gene=tmp[-1] # currently no easy way to grab this here
			pStr = tmp[6]
			rt.write(chro+'\t'+pos+'\t'+gene+'\t'+eff+'\t'+pStr+'\n')
		except:
			print 'issue with line: '+line+' in '+inPath
	fin.close()
	rt.close()

def _cbat2VarTable(inPath,outPath):
	"""Take the standard cifBat output at inPath
	and put it into our 'results table' format 
	at outPath.
	
	"""
	#NOTE: would be more flexible/consistent in long term to combine all
	# *2VarTable methods into one method, limiting test specific code to different file parsers;
	# however, in the interest of time and uncertainty of reuse, we separated them out for now.
	fin = genUtil.open2(inPath)
	rt = open(outPath,'w')
	rt.write('chr\tpos\tgene\teffect\tp-value\n')
	# first line is comment:
	fin.next()
	for line in fin:
		try:
			tmp = line.strip().split()
			name = tmp[0].split(':')
			chro = name[0]
			pos = name[1]
			eff=tmp[16]
			gene=tmp[-1] # currently no easy way to grab this here
			pStr = tmp[17]
			rt.write(chro+'\t'+pos+'\t'+gene+'\t'+eff+'\t'+pStr+'\n')
		except:
			print 'issue with line: '+line+' in '+inPath
	fin.close()
	rt.close()

def filterVarBatch(inTablePath,outTablePath,varBatchList):
	"""Filter the table at inTablePath for the variant batch features
	in the varBatchList list and save it to outTablePath.
	"""
	inTable = open(inTablePath)
	outTable = open(outTablePath,'w')

	line = inTable.next()
	outTable.write(line)

	for line in inTable:
		tmp = line.split('\t')
		var = tmp[0]+'_'+tmp[1]
		if not var in varBatchList:
			outTable.write(line)
	inTable.close()
	outTable.close()

def filterTransBatch(inTablePath,outTablePath,transBatchList):
	"""Filter the table at inTablePath for the transcript batch features
	in the transBatchList list and save it to outTablePath.
	"""
	inTable = open(inTablePath)
	outTable = open(outTablePath,'w')

	line = inTable.next()
	outTable.write(line)

	for line in inTable:
		tmp = line.split('\t')
		trans = 'MiAC_'+tmp[1]+'_'+tmp[0]
		if not trans in transBatchList:
			outTable.write(line)
	inTable.close()
	outTable.close()

def pValueSummaries(inTablePath,outDir, pCol=4,qMax=.1,runFDR=True,runQQ=True):
	"""Given a table, at inTablePath, that has p-values in column pCol
	create summaries in outDir.  
	Including an abbreviated table of all rows with FDR<qMax
	and a q-q plot.
	"""
	if runFDR or runQQ:
		inTable = open(inTablePath)
		
		#skip header
		tmp = inTable.next()
		pValues = np.array([line.strip().split('\t')[pCol] for line in inTable],dtype=float)
		inTable.close()

		if runQQ:
			# create qq-plot
			statsUtil.qqPlot(pValues,outDir+'/qqPlot.png')

		if runFDR:

			# get FDR
			h,q,_ = statsUtil.fdr_bh(pValues,alpha=qMax)

			# Create new table with FDR
			outTable = open(outDir+'/topFDRTable.tsv','w')

			# add some summary lines
			outTable.write('# Hits < qMax ('+str(qMax)+'): '+str(np.sum(h))+'\n')
			outTable.write('# Min p-value: '+str(np.min(pValues))+', q='+str(np.min(q))+'\n')

			inTable = open(inTablePath)
			head = inTable.next()
			head = head.strip().split('\t')
			head.append('FDR_q')
			outTable.write('\t'.join(head)+'\n')

			i = 0
			for line in inTable:
				if q[i]<qMax:
					line = line.strip().split('\t')
					line.append(str(q[i]))
					outTable.write('\t'.join(line)+'\n')
				i = i+1

			inTable.close()
			outTable.close()

			
		




def parseResultTables(outDir,inputPaths,phenoCodes,transManifestPath='',runFDR=True,runQQ=True):
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
	transBatchList	The path to the list of feature names for transcripts that 
			were found to be batch related and should be removed.  If 
			no key, batch filtering will be skipped.
	varBatchList	The path to the list of feature names for variants that were 
			found to be batch related and should be removed.  If no key,
			batch filtering will be skipped.
	
	Variable path = used to simplify the path definitions, for these outputs
			phenotype targets are in separate files, but we can expect 
			them to be in the same dir with the same format for file
			name, but with different strings for each phenotype.
			Therefore use the string $PHENOCODE$ in place of the 
			phenotype string.  This code will search all phenotype
			strings in the list phenoCodes.
	If key is not found parsing for the corresponding result will be skipped.

	transManifestPath can be set to the path for the transcript region manifest to 
	allow for additional annotations to be added to results tables.
	"""
	#NOTE: this is for readme documentation, probably better ways to do it
	lastEdit = '20150321'
	
	### --setup-- ###

	
	# --set up the output dir
	_mkDir(outDir)
	
	# -- top level readme
	
	# may not overwrite all tests, so just append info if DIR exits, let user sort it out
	readme = open(outDir+'/README.txt','a')

	readme.write('------------------------------------------------------------------------------\n')
	readme.write('Result Tables for Genomic Association Tests parsed on '+time.strftime("%c")+'.\n')
	readme.write('This version of parsing code last edited on '+lastEdit+'.\n')
	readme.write('First level dirs refer to individual tests, followed by family member.\n')
	readme.write('Tables are saved per target (phenotype) under the member directory.\n')

	# --transcript annotations 
	# we want to get additional annotations out of the 
	# manifest file, if it is provided
	# doing it here rather than in table generation
	# because it can be used multiple times.
	# using a dictionary is the easiest, but memory intensive.
	regionDic = {}
	if transManifestPath!='':
		readme.write('-Transcript manifest used for annotations at: '+transManifestPath+'\n')
		regionReader = gnmcUtil.RegionManifestReader(transManifestPath)
		for region in regionReader:
			regionDic[region.name]=region
		regionReader.close()
			

	# --batch filter lists
	# if a path to the batch feature list is provided
	# it will be loaded into a dic (NOTE:very memory intensive)
	# and a batch filter feature list will be made in 
	# addition to the unfiltered table
	varBatchList = []
	if inputPaths.has_key('varBatchList'):
		readme.write('-Variant batch filter list at: '+inputPaths['varBatchList']+'\n')
		inList = open(inputPaths['varBatchList'])
		varBatchList = [line.strip().split(':')[-1] for line in inList]

	transBatchList = []
	if inputPaths.has_key('transBatchList'):
		readme.write('-Transcript batch filter list at: '+inputPaths['transBatchList']+'\n')
		inList = open(inputPaths['transBatchList'])
		transBatchList = [line.strip().split(':')[-1] for line in inList]

	
	### --varPheno-FAM-- ###
	if inputPaths.has_key('varPheno-FAM'):
		path = inputPaths['varPheno-FAM']

		# --Create folder
		outDirTmp = outDir+'/varPheno'
		_mkDir(outDirTmp)

		# over writing tests, so over write any readme.
		readmeTmp = open(outDirTmp+'/README.txt','w')
		readmeTmp.write('--Parsing varPheno on '+time.strftime("%c")+'.\n')
		readmeTmp.write('Contains the pairwise test results for variant calls vs phenotypes (defined on families).\n')
		readmeTmp.write('- source path:'+path+'\n')
		

		# split the original file by members, 1 should set terms to members:
		membList = genWF.splitPWResults(path,outDirTmp,1)

		# gonna need those names later
		inFileName = os.path.splitext(os.path.basename(path))[0]

		# go through each family member 
		for memb in membList:
			path2 = outDirTmp+'/'+inFileName+'_subset_'+memb+'.dat'
			outDirTmp2 = outDirTmp+'/'+memb
			_mkDir(outDirTmp2)
			# now split this even further into phenotypes
			phenoList = genWF.splitPWResults(path2,outDirTmp2,4,fName='target')
			# now parse each phenotype
			for pheno in phenoList:
				try:
					path3 = outDirTmp2+'/'+inFileName+'_subset_'+memb+'_subset_'+pheno+'.dat'
					outDirTmp3 = outDirTmp2+'/'+pheno
					_mkDir(outDirTmp3)
					_pw2VarTable(path3,outDirTmp3+'/resultTable.tsv')
					# remove this path, NOTE:comment if you want to keep split output
					os.remove(path3)

					if len(varBatchList)>0:
						filterVarBatch(outDirTmp3+'/resultTable.tsv',outDirTmp3+'/resultTable_batchFiltered.tsv',varBatchList)
						# make some additional summary stuff using the table:
						pValueSummaries(outDirTmp3+'/resultTable_batchFiltered.tsv',outDirTmp3, pCol=4,qMax=.1,runFDR=runFDR,runQQ=runQQ)

				except Exception as e:
					readmeTmp.write('-'+memb+'-'+pheno+' unexpected error:'+str(e)+'\n')
			os.remove(path2)

		# that should be each memb each pheno

			

			
		
		readmeTmp.close()
		# put some notes in original readme
		readme.write('-varPheno DONE\n')
		readme.write('--results dir:'+outDirTmp+'\n')
	else:
		readme.write('-varPheno SKIPPED\n')




	### --varBatch-FAM-- ###
	if inputPaths.has_key('varBatch-FAM'):
		path = inputPaths['varBatch-FAM']

		# --Create folder
		outDirTmp = outDir+'/varBatch'
		_mkDir(outDirTmp)

		# over writing tests, so over write any readme.
		readmeTmp = open(outDirTmp+'/README.txt','w')
		readmeTmp.write('--Parsing varBatch on '+time.strftime("%c")+'.\n')
		readmeTmp.write('Contains the pairwise test results for variant calls vs batch features (defined on individuals).\n')
		readmeTmp.write('- source path:'+path+'\n')
		

		#NOTE: no family members in batch for variant tests, these are done directly agains vcf and individually defined batch features!

		
		# split the original file by phenotypes, 3 should set terms to pheotypes for individual feature names:
		phenoList = genWF.splitPWResults(path,outDirTmp,3,fName='target')

		# gonna need those names later
		inFileName = os.path.splitext(os.path.basename(path))[0]

		  
		for pheno in phenoList:
			try:
				path2 = outDirTmp+'/'+inFileName+'_subset_'+pheno+'.dat'
				outDirTmp2 = outDirTmp+'/'+pheno
				_mkDir(outDirTmp2)
				_pw2VarTable(path2,outDirTmp2+'/resultTable.tsv')
				# remove this path, NOTE:comment next line if you want to keep the full pairwise output per phenotype
				os.remove(path2)

				# make some additional summary stuff using the table:
				pValueSummaries(outDirTmp2+'/resultTable.tsv',outDirTmp2, pCol=4,qMax=.1,runFDR=runFDR,runQQ=runQQ)

			except Exception as e:
				readmeTmp.write('-'+memb+'-'+pheno+' unexpected error:'+str(e)+'\n')

		readmeTmp.close()
		# put some notes in original readme
		readme.write('-varBatch DONE\n')
		readme.write('--results dir:'+outDirTmp+'\n')
	else:
		readme.write('-varBatch SKIPPED\n')



	### --transPheno-FAM-- ##
	if inputPaths.has_key('transPheno-FAM'):
		path = inputPaths['transPheno-FAM']

		# --Create folder
		outDirTmp = outDir+'/transPheno'
		_mkDir(outDirTmp)

		# over writing tests, so over write any readme.
		readmeTmp = open(outDirTmp+'/README.txt','w')
		readmeTmp.write('--Parsing transPheno on '+time.strftime("%c")+'.\n')
		readmeTmp.write('Contains the pairwise test results for transcript burden vs phenotypes (defined on families).\n')
		readmeTmp.write('Originally we expect the transcript burden to be the minor allele count (MiAC) of ppc variants;\n')
		readmeTmp.write('However, that level of detail is not automatically confirmed by this script.\n')
		readmeTmp.write('- source path:'+path+'\n')
		

		# split the original file by members, 1 should set terms to members:
		membList = genWF.splitPWResults(path,outDirTmp,1)

		# gonna need those names later
		inFileName = os.path.splitext(os.path.basename(path))[0]

		# go through each family member 
		for memb in membList:
			path2 = outDirTmp+'/'+inFileName+'_subset_'+memb+'.dat'
			outDirTmp2 = outDirTmp+'/'+memb
			_mkDir(outDirTmp2)
			# now split this even further into phenotypes
			phenoList = genWF.splitPWResults(path2,outDirTmp2,4,fName='target')
			# now parse each phenotype
			for pheno in phenoList:
				try:
					path3 = outDirTmp2+'/'+inFileName+'_subset_'+memb+'_subset_'+pheno+'.dat'
					outDirTmp3 = outDirTmp2+'/'+pheno
					_mkDir(outDirTmp3)
					_pw2TransTable(path3,outDirTmp3+'/resultTable.tsv',regionDic=regionDic)
					# remove this path, NOTE:comment if you want to keep split output
					os.remove(path3)
					if len(transBatchList)>0:
						filterTransBatch(outDirTmp3+'/resultTable.tsv',outDirTmp3+'/resultTable_batchFiltered.tsv',transBatchList)
						
						# make some additional summary stuff using the table:
						pValueSummaries(outDirTmp3+'/resultTable_batchFiltered.tsv',outDirTmp3, pCol=7,qMax=.1,runFDR=runFDR,runQQ=runQQ)

				except Exception as e:
					readmeTmp.write('-'+memb+'-'+pheno+' unexpected error:'+str(e)+'\n')
			os.remove(path2)

		# that should be each memb each pheno

			
		
		readmeTmp.close()
		# put some notes in original readme
		readme.write('-transPheno DONE\n')
		readme.write('--results dir:'+outDirTmp+'\n')
	else:
		readme.write('-transPheno SKIPPED\n')



	### --transBatch-FAM-- ###
	if inputPaths.has_key('transBatch-FAM'):
		path = inputPaths['transBatch-FAM']

		# --Create folder
		outDirTmp = outDir+'/transBatch'
		_mkDir(outDirTmp)

		# over writing tests, so over write any readme.
		readmeTmp = open(outDirTmp+'/README.txt','w')
		readmeTmp.write('--Parsing transBatch on '+time.strftime("%c")+'.\n')
		readmeTmp.write('Contains the pairwise test results for transcript burden vs batches (defined on families).\n')
		readmeTmp.write('Originally we expect the transcript burden to be the minor allele count (MiAC) of ppc variants;\n')
		readmeTmp.write('However, that level of detail is not automatically confirmed by this script.\n')
		readmeTmp.write('- source path:'+path+'\n')
		

		# split the original file by members, 1 should set terms to members:
		membList = genWF.splitPWResults(path,outDirTmp,1)

		# gonna need those names later
		inFileName = os.path.splitext(os.path.basename(path))[0]

		# go through each family member 
		for memb in membList:
			path2 = outDirTmp+'/'+inFileName+'_subset_'+memb+'.dat'
			outDirTmp2 = outDirTmp+'/'+memb
			_mkDir(outDirTmp2)
			# now split this even further into phenotypes
			phenoList = genWF.splitPWResults(path2,outDirTmp2,4,fName='target')
			# now parse each phenotype
			for pheno in phenoList:
				try:
					path3 = outDirTmp2+'/'+inFileName+'_subset_'+memb+'_subset_'+pheno+'.dat'
					outDirTmp3 = outDirTmp2+'/'+pheno
					_mkDir(outDirTmp3)
					_pw2TransTable(path3,outDirTmp3+'/resultTable.tsv',regionDic=regionDic)
					# remove this path, NOTE:comment if you want to keep split output
					os.remove(path3)

					# make some additional summary stuff using the table:
					pValueSummaries(outDirTmp3+'/resultTable.tsv',outDirTmp3, pCol=7,qMax=.1,runFDR=runFDR,runQQ=runQQ)

				except Exception as e:
					readmeTmp.write('-'+memb+'-'+pheno+' unexpected error:'+str(e)+'\n')
			os.remove(path2)

		# that should be each memb each pheno

			
		
		readmeTmp.close()
		# put some notes in original readme
		readme.write('-transBatch DONE\n')
		readme.write('--results dir:'+outDirTmp+'\n')
	else:
		readme.write('-transBatch SKIPPED\n')


	### --eigenstrat-M-- ###
	if inputPaths.has_key('eigenstrat-M'):
		path = inputPaths['eigenstrat-M']

		# --Create folder
		outDirTmp = outDir+'/eigenstrat'
		_mkDir(outDirTmp)
		outDirTmp = outDir+'/eigenstrat/M'
		_mkDir(outDirTmp)

		# over writing tests, so over write any readme.
		readmeTmp = open(outDirTmp+'/README.txt','w')
		readmeTmp.write('--Parsing eigenstrat on '+time.strftime("%c")+'.\n')
		readmeTmp.write('Contains the eigenstrat test results for variants in mother vs phenotypes (defined on families).\n')
		readmeTmp.write('Originally we expect the eigenstrat test to adjust for the first 10 PC;\n')
		readmeTmp.write('However, that level of detail is not automatically confirmed by this script.\n')
		readmeTmp.write('- source path:'+path+'\n')


		# different files for each phenotype:
		for pheno in phenoCodes:
			try:
				path2 = path.replace("$PHENOCODE$",pheno)
				outDirTmp2 = outDirTmp+'/'+pheno
				_mkDir(outDirTmp2)
				_es2VarTable(path2,outDirTmp2+'/resultTable.tsv')
				
				if len(varBatchList)>0:
					filterVarBatch(outDirTmp2+'/resultTable.tsv',outDirTmp2+'/resultTable_batchFiltered.tsv',varBatchList)
					# make some additional summary stuff using the table:
					pValueSummaries(outDirTmp2+'/resultTable_batchFiltered.tsv',outDirTmp2, pCol=4,qMax=.1,runFDR=runFDR,runQQ=runQQ)

			except Exception as e:
			    readmeTmp.write('-'+pheno+' unexpected error:'+str(e)+'\n')

		

		readmeTmp.close()
		# put some notes in original readme
		readme.write('-eigenstrat-M DONE\n')
		readme.write('--results dir:'+outDirTmp+'\n')
	else:
		readme.write('-eigenstrat-M SKIPPED\n')


	### --eigenstrat-F-- ###
	if inputPaths.has_key('eigenstrat-F'):
		path = inputPaths['eigenstrat-F']

		# --Create folder
		outDirTmp = outDir+'/eigenstrat'
		_mkDir(outDirTmp)
		outDirTmp = outDir+'/eigenstrat/F'
		_mkDir(outDirTmp)

		# over writing tests, so over write any readme.
		readmeTmp = open(outDirTmp+'/README.txt','w')
		readmeTmp.write('--Parsing eigenstrat on '+time.strftime("%c")+'.\n')
		readmeTmp.write('Contains the eigenstrat test results for variants in father vs phenotypes (defined on families).\n')
		readmeTmp.write('Originally we expect the eigenstrat test to adjust for the first 10 PC;\n')
		readmeTmp.write('However, that level of detail is not automatically confirmed by this script.\n')
		readmeTmp.write('- source path:'+path+'\n')


		# different files for each phenotype:
		for pheno in phenoCodes:
			try:
				path2 = path.replace("$PHENOCODE$",pheno)
				outDirTmp2 = outDirTmp+'/'+pheno
				_mkDir(outDirTmp2)
				_es2VarTable(path2,outDirTmp2+'/resultTable.tsv')

				if len(varBatchList)>0:
					filterVarBatch(outDirTmp2+'/resultTable.tsv',outDirTmp2+'/resultTable_batchFiltered.tsv',varBatchList)
					# make some additional summary stuff using the table:
					pValueSummaries(outDirTmp2+'/resultTable_batchFiltered.tsv',outDirTmp2, pCol=4,qMax=.1,runFDR=runFDR,runQQ=runQQ)

			except Exception as e:
			    readmeTmp.write('-'+pheno+' unexpected error:'+str(e)+'\n')


		readmeTmp.close()
		# put some notes in original readme
		readme.write('-eigenstrat-F DONE\n')
		readme.write('--results dir:'+outDirTmp+'\n')
	else:
		readme.write('-eigenstrat-F SKIPPED\n')


	### --eigenstrat-NB-- ###
	if inputPaths.has_key('eigenstrat-NB'):
		path = inputPaths['eigenstrat-NB']

		# --Create folder
		outDirTmp = outDir+'/eigenstrat'
		_mkDir(outDirTmp)
		outDirTmp = outDir+'/eigenstrat/F'
		_mkDir(outDirTmp)

		# over writing tests, so over write any readme.
		readmeTmp = open(outDirTmp+'/README.txt','w')
		readmeTmp.write('--Parsing eigenstrat on '+time.strftime("%c")+'.\n')
		readmeTmp.write('Contains the eigenstrat test results for variants in newborn vs phenotypes (defined on families).\n')
		readmeTmp.write('Originally we expect the eigenstrat test to adjust for the first 10 PC;\n')
		readmeTmp.write('However, that level of detail is not automatically confirmed by this script.\n')
		readmeTmp.write('- source path:'+path+'\n')


		# different files for each phenotype:
		for pheno in phenoCodes:
			try:
				path2 = path.replace("$PHENOCODE$",pheno)
				outDirTmp2 = outDirTmp+'/'+pheno
				_mkDir(outDirTmp2)
				_es2VarTable(path2,outDirTmp2+'/resultTable.tsv')

				if len(varBatchList)>0:
					filterVarBatch(outDirTmp2+'/resultTable.tsv',outDirTmp2+'/resultTable_batchFiltered.tsv',varBatchList)

					# make some additional summary stuff using the table:
					pValueSummaries(outDirTmp2+'/resultTable_batchFiltered.tsv',outDirTmp2, pCol=4,qMax=.1,runFDR=runFDR,runQQ=runQQ)

			except Exception as e:
			    readmeTmp.write('-'+pheno+' unexpected error:'+str(e)+'\n')

		readmeTmp.close()
		# put some notes in original readme
		readme.write('-eigenstrat-NB DONE\n')
		readme.write('--results dir:'+outDirTmp+'\n')
	else:
		readme.write('-eigenstrat-NB SKIPPED\n')


	### --cifBat-FAM-- ###
	if inputPaths.has_key('cifBat-FAM'):
		path = inputPaths['cifBat-FAM']

		# --Create folder
		outDirTmp = outDir+'/cifBat'
		_mkDir(outDirTmp)
		outDirTmp = outDir+'/cifBat/FAM'
		_mkDir(outDirTmp)

		# over writing tests, so over write any readme.
		readmeTmp = open(outDirTmp+'/README.txt','w')
		readmeTmp.write('--Parsing cifBat on '+time.strftime("%c")+'.\n')
		readmeTmp.write('Contains the cifBat test results for variants in trios vs phenotypes (defined on families).\n')
		readmeTmp.write('Taking results from the additive model.\n')
		readmeTmp.write('- source path:'+path+'\n')


		# different files for each phenotype:
		for pheno in phenoCodes:
			try:
				path2 = path.replace("$PHENOCODE$",pheno.replace('_',''))
				outDirTmp2 = outDirTmp+'/'+pheno
				_mkDir(outDirTmp2)
				_cbat2VarTable(path2,outDirTmp2+'/resultTable.tsv')

				if len(varBatchList)>0:
					filterVarBatch(outDirTmp2+'/resultTable.tsv',outDirTmp2+'/resultTable_batchFiltered.tsv',varBatchList)
					# make some additional summary stuff using the table:
					pValueSummaries(outDirTmp2+'/resultTable_batchFiltered.tsv',outDirTmp2, pCol=4,qMax=.1,runFDR=runFDR,runQQ=runQQ)
				
			except Exception as e:
			    readmeTmp.write('-'+pheno+' unexpected error:'+str(e)+'\n')

		readmeTmp.close()
		# put some notes in original readme
		readme.write('-cifBat-FAM DONE\n')
		readme.write('--results dir:'+outDirTmp+'\n')
	else:
		readme.write('-cifBat-FAM SKIPPED\n')
	
	readme.close()



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
#	outDir = '/isb/rtasseff/results/transcript_20150216'
#	inSortedPWPath =  outDir+'/genPW1_data_GNMC_Trans_PW_20150206_vs_BATCH_GNMC_20150109_out_sorted.dat'
#	outFDRFilteredPWPath =  outDir+'/sigFDRPWOut.dat'
#	outBatchFList = outDir+'/fList_trans_batch_20150320.dat'
#
#	statsUtil.fdr_bh_filterSortFile(inSortedPWPath,outFDRFilteredPWPath,alpha=.1,col=5,logTrans=True)
#	getSimpleList(outFDRFilteredPWPath,outBatchFList)

#
#	##### PRE run multiple small pairwise to get phenotype results, -->
#	# NOTE:designed to be called from PBS script for qsub -J
#	# in 20150210 on DF5 took < 5min per job
#	# after PBS run use the pw_region_postproc.sh shell script 
#	# and then the POST code segment '##### POST run multiple small pairwise results'
#	tag = sys.argv[1]
#	outDir = '/isb/rtasseff/results/var_pheno_20150210'
#	testFMPath = '/isb/rtasseff/data/featureMatrices/data_VCF_FM_regions_20150304/'+tag+'.Filtered.fm.gz'
#	targFMPath = '/isb/rtasseff/data/featureMatrices/data_CLIN_Critical_Phenotype_20150203.fm'
#	fullOutPath = outDir+'/fullPWOut_'+tag+'.dat'
#	nbListPath = '/isb/rtasseff/data/support/repNBList.tsv'
#	
#	runPW_INDvFAM(testFMPath,targFMPath,fullOutPath,outDir,nbListPath)		
#	# <--


	#### Parse results tables, -->
	# giong through tests outputs to create tables of 
	# all results within a set dir structure.
	# Uses format specifications as of 20150310
	
	# give paths in dic
	inputPaths = {}
	# added tags to break up jobs using the PBS job array
	tag = int(sys.argv[1])
	if tag==1:
		inputPaths['varPheno-FAM'] = '/isb/rtasseff/results/var_pheno_20150210/fullPWOut.dat'
	elif tag==2:
		inputPaths['varBatch-FAM'] = '/isb/rtasseff/results/var_batch_20150204/fullPWOut.dat'
	elif tag==3:
		inputPaths['transPheno-FAM'] = '/isb/rtasseff/results/transcript_20150216/genPW1_data_GNMC_Trans_PW_20150206_vs_data_CLIN_Critical_Phenotype_20150213_out.dat'
	elif tag==4:
		inputPaths['transBatch-FAM'] = '/isb/rtasseff/results/transcript_20150216/genPW1_data_GNMC_Trans_PW_20150206_vs_BATCH_GNMC_20150109_out.dat'
	elif tag==5:
		inputPaths['eigenstrat-M'] = '/isb/rtasseff/results/eig_m_20150211/OUTPUT/final$PHENOCODE$.out'
	elif tag==6:
		inputPaths['eigenstrat-F'] = '/isb/rtasseff/results/eig_f_20150211/output/final$PHENOCODE$.out'
	elif tag==7:
		inputPaths['eigenstrat-NB'] = '/isb/rtasseff/results/eig_nb_20150211/output/final$PHENOCODE$.out'
	elif tag==8:
		inputPaths['cifBat-FAM'] = '/bigdata0/users/vdhankani/CIFBAT/OUTPUT/$PHENOCODE$/$PHENOCODE$Autosomal.tdt.out.gz'

	# batch list paths:
	inputPaths['transBatchList'] = '/isb/rtasseff/results/transcript_20150216/fList_trans_batch_20150320.dat'
	inputPaths['varBatchList'] = '/isb/rtasseff/results/var_batch_20150204/fList_var_batch_20150204.dat'
	
	# define phenotype codes in list
	phenoCodes = ['1n2v4',
	'Hypertension_Related',
	'Incompetent_Cervix',
	'Preeclampsia',
	'Uterine_Related',
	'1v4',
	'IdiopathicNA',
	'Infection_Related',
	'Preterm',
	'History_PTB',
	'Immune_Related',
	'Placenta_Related',
	'Prom_Related',
	'TermCategory',
	'Gestational_Age_at_Delivery']

	outDir = '/isb/rtasseff/results/result_tables_20150326'

	transManifestPath = '/isb/rtasseff/data/transcripts_20141125/transcriptManifest_20141125.dat'

	parseResultTables(outDir,inputPaths,phenoCodes,transManifestPath=transManifestPath,runFDR=True,runQQ=True)

	#### <---, Parse results tables


if __name__ == '__main__':
	main()
