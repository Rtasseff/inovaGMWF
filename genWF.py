""" Methods for a standard feature matrix based workflow for use in ISB-INOVA projects.
Focus on creating a general, semi-automated, preprocessing and analysis workflow 
for study 101, and for reuse on similar studies.

- Methods to create FMs out of source files
- Methods to merge FMs, check consistency and filter for stats purposes
- Methods to run pairwise analysis and summarize results (currently for medium to small runs only, metadata)

The goal was to organize all of the analysis workflow that 
was validated on DF4 and reused on DF5 for study 101, while 
allowing for more general use when possible. However, for 
several reasons (including time constraints, prototyping, 
procedures run in parallel, source data on different servers)
we have limited this to only a subset of the analysis and preprocessing.
Some procedure must be run directly from custWF and gnmcUtil or even 
using by resources not included in this package.

View the README for more information, including links to extended 
documentation on the more complete workflow.

Designed around study 101,
Some limitations on this relate to primary phenotype being defined by family.
as in NB of FAM was PTB. For some sources, data on all family members is available  
but in others this is not true (PTB on NB only, molecular data on M only)
so we have organized this around the family unit (suffix -FAM).  
In cases with twins, only one is chosen as representative.   
We note that this limitation is introduced for systematic analysis 
of feature matrix like input data.
"""

import numpy as np
import subprocess
import time
import shutil
import statsUtil 
import logging
import sys
import ConfigParser
import argparse
import os
import genUtil
import gnmcUtil
import random

floatFmtStr = '%05.4E'
nanValues = ['NA','NaN','na','nan']

disc="python workflow for standard genomic/molecular data analysis, INOVA 101."
version='0.1.0'

def _checkPatientID(patientID,studyID='101',allowedSuffix=['FAM']):
	"""Check that ID is consistent with expectations"""
	tmp = patientID.split('-')
	if tmp[0] != studyID:
		raise ValueError ("ERR:0001 unexpected study id "+str(tmp[0]))
	elif (not allowedSuffix[0]=='any') and ('-'.join(tmp[2:]) not in allowedSuffix): 
		raise ValueError ("ERR:0002 unexpected member id "+'-'.join(tmp[2:])+'. Suffixes limited to '+str(allowedSuffix))


def getPatientOrder(fmPath,studyID='101',allowedSuffix=['FAM']):
	"""Get the sample order from a properly formated, pre-existing
	feature matrix.

	fmPath	str, path to standard feature matrix
	"""
	fin = open( fmPath )
	line = fin.next() # first line is header with info
	patients = line.strip().split('\t')[1:] # first position is '.' place holder
	for patient in patients:
		_checkPatientID(patient,studyID=studyID,allowedSuffix=allowedSuffix)	
	return patients



def _catOrderedFM(filenames,outfile,skipHeader=True):
	for fname in filenames:
		with open(fname) as infile:
			if skilpHeader:infile.next()
			for line in infile:
				outfile.write(line)

			
		


def filterFM(fmInPath,fmOutPath,maxMiss=.9,minObs=5):
	"""Filter a given feature matrix for statistically 
	poor features for use in pairwise code.
	fmInPath	str, path to the input feature matrix
	fmOutPath 	str, path for the output feature matrix
	maxMiss		flaot, maximum fraction of missing calls allowed
	minObs	int, min number of observations per catigorey
	"""
	fin = open(fmInPath)
	fout = open(fmOutPath,'w')
	
	fout.write('.')
	line = fin.next().strip().split('\t')[1:]
	n = len(line)
	for sampID in line:
		fout.write('\t'+sampID)
	fout.write('\n')


		
	for line in fin:
		tmp = line.strip().split('\t')
		fname = tmp[0].split(':')
		
		# skip if this does not meet stats preprocessing requierments
		skip = False	
		if len(tmp[1:])!=n: 
			skip = True
			logging.warning('removing '+':'.join(fname)+', in correct number of columns.')
		if np.sum(np.array(tmp[1:],dtype=str)=='NA')>n*maxMiss:
			skip = True
			logging.warning('removing '+':'.join(fname)+', too many missing values.')
		if fname[0]=='C':
			# check for size
			unique = list(set(tmp[1:]))
			if len(unique)>30:
				skip = True
				logging.warning('removing '+':'.join(fname)+',too many catigories.')
		if fname[0]=='B':
			# check for size
			unique = list(set(tmp[1:]))
			if ('NA' in unique and len(unique)>3) or ('NA' not in unique and len(unique)>2):
				skip = True
				logging.warning('removing '+':'.join(fname)+',too many catigories for binary.')

		if fname[0]=='B' or fname[0]=='C':
			unique = list(set(tmp[1:]))
			values = np.array(tmp[1:],dtype=str)
			for label in unique:
				if np.sum(label==values)<minObs:
					skip = True
					logging.warning('removing '+':'.join(fname)+', not enough samples for a group.')
					break
		if fname[0]=='N':
			values = statsUtil._getFloat(np.array(tmp[1:],dtype=str))
			uniqueValues = list(set(values[~np.isnan(values)]))
			if len(uniqueValues)==2:
				newValues = np.array(np.zeros(len(values))*np.nan,dtype=str)
				fname[0] = 'B'
				newValues[values==uniqueValues[0]]='0'
				newValues[values==uniqueValues[1]]='1'
				tmp[1:] = newValues[:]

			elif np.all(values[~np.isnan(values)] > -1E-52) and np.median(values[~np.isnan(values)]) < 1E-52:
				# if the matrix is mostly zeros change it to binary, zero or else
				fname[0] = 'B'
				newValues = np.array(np.zeros(len(values))*np.nan,dtype=str)
				newValues[values > 1E-52] = '1'
				newValues[values < 1E-52] = '0'
				tmp[1:] = newValues[:]
			

		if not skip:
			fout.write(':'.join(fname)+'\t'+'\t'.join(tmp[1:])+'\n')
	fout.close()
	fin.close()

	





def catFM(fMNames,sampleIDs,foutPath,studyID='101',allowedSuffix=['FAM']):
	"""Concatenate multiple feature matrices together.
	Assumes first line is header for sample labels,
	and use sampleID to determine columns and their order.
	nan will be used to fill in missing samples,
	sample in sampleID but not in a FM header.
	fMnames	list of str paths for feature matrices to be cat'ed
	sampleIDs	list of str sample ids
	"""
	fout = open(foutPath,'w')
	fout.write('.\t'+'\t'.join(sampleIDs)+'\n')
	for finName in fMNames:
		with open(finName) as fin:
			labels = np.array(fin.next().strip().split('\t')[1:],dtype='|S15')

			# doing some checking here to see labels are as expected
			for label in labels:
				 _checkPatientID(label,studyID,allowedSuffix)

# originally this was all hard coded to assume matrices were correct for 101, but not labled correctly
#				tmp = labels[i].split('-')
#				# some of the logic here has been hard coded
#				# to assume study 101
#				if tmp[0]!='101':
#					ValueError('ERR:0010, In feature matrix at '+finName+', study ID is not 101. Please change source code for use with other studies')
#				labels[i] = tmp[0]+'-'+tmp[1]+'-FAM' # by assumption, mapping all col to -FAM, hard coded for 101.

			for sampleID in sampleIDs:
				if not np.any(sampleID==labels):
					logging.warning('WARN0002: sample ID '+sampleID+' not found in feature matrix at '+ finName)
				if np.sum(sampleID==labels)>1:
					logging.warning('WARN0003: sample ID '+sampleID+' had multiple entries in feature matrix at '+ finName+'. Using only the first.')

			# start appending data
			for line in fin:
				data = line.strip().split('\t')
				fout.write(data[0])
				data = np.array(data[1:],dtype=str)
				if len(data)!=(len(labels)):
					raise ValueError('ERR:0011, In feature matrix at '+finName+', the rows do not have equal lengths')
				
				for sampleID in sampleIDs:
					if not np.any(sampleID==labels):fout.write('\tnan')
					else:fout.write('\t'+data[sampleID==labels][0])
				fout.write('\n')

	fout.close()
					
			
def _getFeatureNamesByGroup(fmPath,lookfor,feildInd):
	"""Going through feature names in fmPath,
	record names that contain strings in the list 
	lookfor in the corresponding field index, fieldInd.
	"""
	n = len(lookfor)
	if n!=len(feildInd):
		raise ValueError('ERR:00013, lookfor string list must be equal to feildInd')
	# perp list of lists
	nameLists = []
	for i in range(n):
		nameLists.append([])
	fin = open(fmPath)
	# skip header 
	line = fin.next()
	for line in fin:
	# go through each sample name:
		fname = line.strip().split('\t')[0].split(':')
		for i in range(n):
		# go through all lookfor strings
			if fname[feildInd[i]].find(lookfor[i])>=0:
				nameLists[i].append(':'.join(fname))
	return nameLists

def _mkPairFile(nameList1,nameList2,foutPath):
	"""Create a file at foutPath that contains
	all pairs from nameList1 and 2. This file 
	can be used to guide the pairwise analysis.
	"""
	with open(foutPath,'w') as fout:
		for name1 in nameList1:
			for name2 in nameList2:
				fout.write(name1+'\t'+name2+'\n')

def mkPSPairFile(fmPath,foutPath,pairType='batch'):
	"""Make a preset (PS) pair file for use in pairwise 
	analysis on the feature matrix at fmPath.
	Save it to foutPath.
	Current presets pairType values:
	batch	compare batch features to critical phenotypes to detect bias sampling
	QC 	compare QC features to critical phenotypes to detect bias in confounding factors
	QCBatch	compare batch features to QC features to detect process changes
	"""
	if pairType=='batch':
		feildInd = [3,3,3]
		lookfor = ['BATCH','Critical_Phenotype','summary_data']
		nameLists =  _getFeatureNamesByGroup(fmPath,lookfor,feildInd)
		markerList = nameLists[0]
		# change as of 20141222 to include summary_data features with Critical_Phenotype features
		predictorList = np.append(nameLists[1],nameLists[2])
		_mkPairFile(markerList,predictorList,foutPath)
	elif pairType=='QC':
		feildInd = [3,3,3]
		lookfor = ['QC','Critical_Phenotype','summary_data']
		nameLists =  _getFeatureNamesByGroup(fmPath,lookfor,feildInd)
		markerList = nameLists[0]
		# change as of 20141222 to include summary_data features with Critical_Phenotype features
		predictorList = np.append(nameLists[1],nameLists[2])
		_mkPairFile(markerList,predictorList,foutPath)
	elif pairType=='QCBatch':
		feildInd = [3,3]
		lookfor = ['QC','BATCH']
		nameLists =  _getFeatureNamesByGroup(fmPath,lookfor,feildInd)
		_mkPairFile(nameLists[0],nameLists[1],foutPath)

	else:
		raise ValueError('ERR:00014, pairType value is not in preset list.')




			

def runPairwise(fMPath,outDir,outPath,pwWhich='/titan/cancerregulome8/TCGA/scripts/pairwise-2.0.0-current',pairFile='',fdrFilter=''):
	"""Run the pairwise code found at pwWhich on feature matrix at
	fMName and save the output to outName.
	"""
	# construct the call to the pairwise script
	call = pwWhich
	# check to see if we have a file to indicate the pairs to calculate
	if pairFile!='':
		call = call+' --by-name '+pairFile

	# check to see if we have a setting for fdr filtered output
	if fdrFilter!='':
		call = call+' -q '+str(fdrFilter)


	call = call +' '+fMPath+' '+outPath
	# redirecting output to an info file
	with open(outDir+"/stdout.out",'w') as stdout:
		subprocess.check_call(call,shell=True,stdout=stdout)

def _replaceNANs(x):
	"""Given the string array, replace
	any usable na value with 'nan' and return.
	"""
	for i in range(len(x)):
		if x[i] in nanValues:
			x[i] = 'nan'
	return (x)


def writePWSumLongSepTar(pwPath, repPath, minLogQ=2.0,nTopPairs=2):
	"""Takes the pairwise output at pwPath 
	and parses it into a summary report at repPath
	by filtering out pairs with a log_10(FDR Q)<minLogQ
	and by dividing the ouput up by the feature name of the 
	target features.

	Assumes pairwise output format from pairwise-2.0.0
	& column 0 represent covariates of interest (enforced in workflow by suppling pairFile, see runPairWise)
	& column 1 represents target of interest (enforced in workflow by suppling pairFile, see runPairWise)
	& organizing by target feature name which is the full entry of column 1.
	& assumes pairwise output can fit in np array

	returns a list of top scoring pairs (nTopPairs) for each data source
	"""
	pw = np.loadtxt(pwPath,dtype=str,delimiter='\t')
	# get data sources
	n,m = pw.shape
	features = pw[:,1] 
	# get all data sources
	uniqueFeatures = list(set(features))
	# save top pairs
	topPairsList = []
	# prepare to write report 
	with open(repPath,'w') as rep:
		rep.write('Suammary output file for pairwise analysis run: '+ time.strftime("%c")+'.\n')
		rep.write('Organized by target feature names found in ** **.\n')
		rep.write('Informative log should be in same directory.\n')
		rep.write('Showing results for pairs with log_10(FDR Q)>'+str(minLogQ)+'.\n')
		# write each section 1 by 1
		for feature in uniqueFeatures:
			rep.write('\nResults for target feature: **'+feature+'**-\n')
			
			pwTmp = pw[features==feature]
			p = np.array(pwTmp[:,5],dtype=str)
			p = _replaceNANs(p)
			p = 10.0**(-1*np.array(pwTmp[:,5],dtype=float))
			_,q,_ = statsUtil.fdr_bh(p)
			q = -1*np.log10(q)
			
			if np.any(q>minLogQ):
				# filter if too low
				pwTmp = pwTmp[q>minLogQ]
				q = q[q>minLogQ]
				# sort
				ind = np.argsort(q)
				pwTmp = pwTmp[ind[-1::-1]]
				q = q[ind[-1::-1]]

				nSig = len(pwTmp)
				rep.write(str(nSig) + ' significant hits found:\n')

				# print individual results:
				rep.write('\tFName_1\tFName_2\tRho\tlog_10(p)\tlog_10(FDR Q)\n')
				for i in range(nSig):
					line = pwTmp[i]
					rep.write('\t'+line[0]+'\t'+line[1]+'\t'+line[3]+'\t'+line[5]+'\t'+str(q[i])+'\n')
					# record top pairs
					if i < nTopPairs:
						topPairsList.append([line[0],line[1]])

			else:
				rep.write('\t------None > min log(q)--------\n')
	return(topPairsList)


def writePWSumLong(pwPath, repPath, minLogQ=2.0,nTopPairs=2):
	"""Takes the pairwise output at pwPath 
	and parses it into a summary report at repPath
	by filtering out pairs with a log_10(FDR Q)<minLogQ
	and by dividing the ouput up by data source.

	Assumes pairwise output format from pairwise-2.0.0
	& column 0 represent covariates of interest (enforced in workflow by suppling pairFile, see runPairWise)
	& column 1 represents target of interest (enforced in workflow by suppling pairFile, see runPairWise)
	& organizing by data source which is index 2 (start at 0) field of the feature name for the covariates.
	& assumes pairwise output can fit in np array

	returns a list of top scoring pairs (nTopPairs) for each data source
	"""
	pw = np.loadtxt(pwPath,dtype=str,delimiter='\t')
	# get data sources
	n,m = pw.shape
	dataSources = np.array(np.zeros(n),dtype='|S15')
	for i in range(n):
		dataSources[i] = pw[i,0].split(':')[2]
	# get all data sources
	uniqueSources = list(set(dataSources))
	# save top pairs
	topPairsList = []
	# prepare to write report 
	with open(repPath,'w') as rep:
		rep.write('Suammary output file for pairwise analysis run: '+ time.strftime("%c")+'.\n')
		rep.write('Organized by data sources found in ** **.\n')
		rep.write('Informative log should be in same directory.\n')
		rep.write('Showing results for pairs with log_10(FDR Q)>'+str(minLogQ)+'.\n')
		# write each section 1 by 1
		for source in uniqueSources:
			rep.write('\nResults for data source: **'+source+'**-\n')
			rep.write('\tFName_1\tFName_2\tRho\tlog_10(p)\tlog_10(FDR Q)\n')
			pwTmp = pw[dataSources==source]
			p = np.array(pwTmp[:,5],dtype=str)
			p = _replaceNANs(p)
			p = 10.0**(-1*np.array(pwTmp[:,5],dtype=float))
			_,q,_ = statsUtil.fdr_bh(p)
			q = -1*np.log10(q)
			
			if np.any(q>minLogQ):
				# filter if too low
				pwTmp = pwTmp[q>minLogQ]
				q = q[q>minLogQ]
				# sort
				ind = np.argsort(q)
				pwTmp = pwTmp[ind[-1::-1]]
				q = q[ind[-1::-1]]
				for i in range(len(pwTmp)):
					line = pwTmp[i]
					rep.write('\t'+line[0]+'\t'+line[1]+'\t'+line[3]+'\t'+line[5]+'\t'+str(q[i])+'\n')
					# record top pairs
					if i < nTopPairs:
						topPairsList.append([line[0],line[1]])

			else:
				rep.write('\t------None > min log(q)--------\n')
	return(topPairsList)
	

def _run2FMPW(FM1Path,FM2Path,pwOutPath,pwRepPath,outDir,pwWhich,samples=[],minLogQ=2.0,nTopPairs=2):
	"""Run pairwise between two exisitng FMs, creates tmp files in outDir
	then removes them. PW output saved to pwOutPath. 
	Summary created via writePWSumLongSepTar and saved to pwRepPath
	Note: if samples not indicated uses samples in FM1, 
	will look for these labels in FM2 (ordernot need to 
	be exact and misisng labels will be set as nan).
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
	_mkPairFile(testF,targF,pairListPath)

	# join FMs for pairwise (could do it faster, but this has mulitple checks
	FMTmpPath = outDir+'/tmpFM_'+str(random.randrange(16**5))+'.tsv'
	catFM([FM1Path,FM2Path],samples,FMTmpPath,studyID='101',allowedSuffix=['FAM'])

	runPairwise(FMTmpPath,outDir,pwOutPath,pwWhich=pwWhich,pairFile=pairListPath)


	# create a summary and plots
	#####NOTE: potential mem issue as this summary code loads the whole pw output as np file object!
	topPairsList = writePWSumLongSepTar(pwOutPath, pwRepPath, minLogQ=minLogQ,nTopPairs=nTopPairs)
	# plot top pairs
	if len(topPairsList)>0:
		genUtil.plotter(topPairsList,FMTmpPath,outDir=outDir)
		logging.info("Plots of top scoring pairs saved in {}".format(outDir))

	os.remove(pairListPath)
	os.remove(FMTmpPath)

	


def _runPresetMetaPW(outDir,outName,pairType,pwWhich,fmPath,minLogQ=2.0,nTopPairs=2):
	pairFile = outDir+'/pairwise_'+outName+'_pairList.dat'
	mkPSPairFile(fmPath,pairFile,pairType=pairType)
	logging.info("Pairwise tests to run indicated in {}.".format(pairFile))
	# run pairwise
	logging.info("Running pairwise code, {}.".format(pwWhich))
	pwOutPath = outDir+'/pairwise_'+outName+'_fullOut.dat'
	runPairwise(fmPath,outDir,pwOutPath,pwWhich=pwWhich,pairFile=pairFile)
	logging.info("Full pairwise output saved at {}.".format(pwOutPath))
	pwSumPath = outDir+'/pairwise_'+outName+'_summary.dat'
	topPairsList = writePWSumLong(pwOutPath,pwSumPath,minLogQ=minLogQ,nTopPairs=nTopPairs)
	logging.info("Summary of pairwise output saved at {}.".format(pwSumPath))
	# plot top pairs
	if len(topPairsList)>0:
		prefix = 'pwPlot_'+outName
		genUtil.plotter(topPairsList,fmPath,outDir=outDir,prefix=prefix)
		logging.info("Plots of top scoring pairs saved in {} with prefix {}_.".format(outDir,prefix))

	


def _parse_CmdArgs(parser):
	parser.add_argument("config",help="Path to the configuration file that contains all settable parameters/options.")
	parser.add_argument("-ow","--over_write", help="Force to overwrite any file in output dir specified in config file, except log.out, which is appended.",
		action="store_true")
	return(parser.parse_args())

def main():

	logName = 'log.out'
	



	# --get the input arguments
	parser = argparse.ArgumentParser(description="Run the"+disc+' Version='+version+'.')
	args = _parse_CmdArgs(parser)

	
	#--setup the config parser
	config = ConfigParser.SafeConfigParser()
	config.read(args.config)


	# --set up the output dir
	# out directory 
	outDir = config.get('genWF','outDir')
	if os.path.exists(outDir):
		# only overwrite if option is set
		if args.over_write:
			logging.warning("Using an existing directory for output, previous files may be overwritten: {}".format(outDir))
		else:
			raise ValueError ("The specified output dir already exists and the -ow command was not used to force overwrite.".format(outDir))
	else:
		os.makedirs(outDir)

	
	# --setup logger and decrease level:
	logging.getLogger('').handlers = []
	logging.basicConfig(filename=outDir+'/'+logName, level=logging.INFO, format='%(asctime)s %(message)s')
	logging.captureWarnings(True)

	# --record some basic information:
	logging.info("--------------------------------------------------------------")
	logging.info("Running {}, {}, version={}...".format(sys.argv[0],disc,version))
	logging.info("--------------------------------------------------------------")

	# --general setup
	logging.info("--Getting some general information.")

	# put a copy of the config file in the output dir
	shutil.copyfile(args.config,outDir+'/wf.cfg')
	logging.info('A copy of the configuration file for this run was saved to '+outDir+'/wf.cfg')
	
	pwWhich = config.get('genWF','pwWhich')	

	
	# need the sample list 
	fmSampPath = config.get('genWF','fmSampPath')
	logging.info("Standardizing to sample list in FM at {}.".format(fmSampPath))
	samples = getPatientOrder(fmSampPath)

	# get some parameters
	# minimum Log Q counted on summary reports
	minLogQ = float(config.get('genWF','minLogQ'))
	# number of max hits to plot data for after summary reports
	nTopPairs = int(config.get('genWF','nTopPairs'))



	
	#########--Get GNMC Data--########
	# get the genomic data from source files
	if config.getboolean('genWF','runGetGNMCData'): # GNMC data ---->
		logging.info("--Getting GNMC Data--\n\tRunning scripts to get & format data from source.")
	
		# -general info-
		# path to tsv of the individual sample ids (ids in VCF header for example) 
		indSampListPath = config.get('genWF','indSampListPath')
		logging.info("List of individual samples to collect GNMC data on at {}.".format(indSampListPath))
		indSampList = np.loadtxt(indSampListPath,dtype=str,delimiter='\t')	

		# need the list of newborns for mapping to correct family
		repNBListPath = config.get('genWF','repNBListPath')
		logging.info("List of representative newborns at {}.".format(repNBListPath))
		nbList = genUtil.getRepNBList(repNBListPath)


		# -Batch-

		if config.getboolean('genWF','runGetGNMCBatch'):
			logging.info("-Getting GNMC BATCH FM-.")  

			# parse the manifest file for gnmc batch info
			manPath =  config.get('genWF','manPath')
			logging.info("Using CG manifest at {}.".format(manPath))
			gnmcIndBatchFMOutPath = config.get('genWF','gnmcIndBatchFMOutPath') 
			gnmcUtil.parseMkGenomeBatchFM(manPath,gnmcIndBatchFMOutPath,indSampList)
			logging.info("Individual GNMC BATCH FM saved at {}.".format(gnmcIndBatchFMOutPath))
			# move to family based FM
			gnmcFamBatchFMOutPath = config.get('genWF','gnmcFamBatchFMOutPath')
			genUtil.indFM2FamFM(gnmcIndBatchFMOutPath,gnmcFamBatchFMOutPath,nbList,samples)
			logging.info("Family based GNMC BATCH FM saved at {}.".format(gnmcFamBatchFMOutPath))

		# -QC-

		if config.getboolean('genWF','runGetGNMCQC'):
			logging.info("-Transforming GNMC QC FM-.")  
			# get the IND FM path:
			gnmcIndQCFMInPath =  config.get('genWF','gnmcIndQCFMInPath')
			# get the ouput FAM FM
			gnmcFamQCFMOutPath =  config.get('genWF','gnmcFamQCFMOutPath')
			genUtil.indFM2FamFM(gnmcIndQCFMInPath,gnmcFamQCFMOutPath,repNBListPath,samples)
			logging.info("The IND QC FM ("+gnmcIndQCFMInPath+") has been transformed to FAM format ("+gnmcFamQCFMOutPath+").")

		# -meta merge-

		if config.getboolean('genWF','runGetGNMCMeta'):
			logging.info("-Merge GNMC metadata-.")  
			# get the out FM path:
			gnmcMDFMOutPath =  config.get('genWF','gnmcMDFMOutPath')
			# merge them
			catFM([gnmcFamBatchFMOutPath,gnmcFamQCFMOutPath],samples,gnmcMDFMOutPath,studyID='101',allowedSuffix=['FAM'])



		# -Transform IND 2 FAM
		indFMPathList = config.get('genWF','indFMPathList')
		if not indFMPathList=='na':
			logging.info("-Transform IND 2 FAM.")
			# parse input paths 
			indFMPathList = indFMPathList.split(',')
			# parse output paths
			famFMPathList = config.get('genWF','famFMPathList').split(',')
			# parse filtered outputs
			filterFMPathList = config.get('genWF','filterFMPathList').split(',')
			if (len(indFMPathList)==len(famFMPathList)) and (len(indFMPathList)==len(filterFMPathList)):
				for i in range(len(indFMPathList)):

					logging.info("Transforming IND FM at {} to FAM FM at {}.".format(indFMPathList[i],famFMPathList[i]))
					genUtil.indFM2FamFM(indFMPathList[i],famFMPathList[i],repNBListPath,samples)

					if not filterFMPathList[i]=='na':
						logging.info("Filtering FAM FM for basic statistical tests.")
						filterFM(famFMPathList[i],filterFMPathList[i],minObs=0)
						logging.info("Filtered FM at {} and saved to {}.".format(famFMPathList[i],filterFMPathList[i]))


			else:
				logging.warning('GNMC Transform IND 2 FAM NOT run as the FM in/out/filter path lists are diffrent sizes.')



	# <----- GNMC data 


	###### --Run metadata analysis-- ######
	if config.getboolean('genWF','runMetadata'): # metadata analysis --->
		logging.info("--Metadata Analysis--\n\tRunning standard association tests on metadata.")

		# -collect data
		if config.getboolean('genWF','runCollectData'):
			# feature matrices to use:
			metaFMs = []
			# mi rna
			mirMDFMPath = config.get('genWF','mirMDFMPath')	
			if mirMDFMPath not in nanValues: metaFMs.append(mirMDFMPath)
			# rna seq
			rnasMDFMPath = config.get('genWF','rnasMDFMPath')
			if rnasMDFMPath not in nanValues: metaFMs.append(rnasMDFMPath)
			# methylation 
			methMDFMPath = config.get('genWF','methMDFMPath')
			if methMDFMPath not in nanValues: metaFMs.append(methMDFMPath)
			# genomic 
			gnmcMDFMPath = config.get('genWF','gnmcMDFMPath')
			if gnmcMDFMPath not in nanValues: metaFMs.append(gnmcMDFMPath)
			# extra data
			extraMDFMPaths = config.get('genWF','extraMDFMPaths')
			if extraMDFMPaths not in nanValues: 
				for path in extraMDFMPaths.strip().split(','):
					metaFMs.append(path)

			clinDataFMPath = config.get('genWF','clinDataFMPath')
			metaFMs.append(clinDataFMPath)

			logging.info("Considering metadata data stored at:")
			for fm in metaFMs:
				logging.info("\t{}.".format(fm))

			# cat matrices together
			combMDFMOutPath = config.get('genWF','combMDFMOutPath')
			catFM(metaFMs,samples,combMDFMOutPath,studyID='101')
			logging.info("Metadata stored in single FM at {}.".format(combMDFMOutPath))

		# filter for pairwise
		if config.getboolean('genWF','runFilterFM'):
			logging.info("Filtering metadata for basic statistical tests.")
			filteredFMInPath = config.get('genWF','filteredFMInPath')
			filteredFMOutPath = config.get('genWF','filteredFMOutPath')
			filterFM(filteredFMInPath,filteredFMOutPath,minObs=0)
			logging.info("Filtered metadata at {} and saved to {}.".format(filteredFMInPath,filteredFMOutPath))



		# -batch vs CP
		if config.getboolean('genWF','runBatchCP'):
			logging.info("-Batch vs. Critical Phenotypes, id bias in sampling.")
			outName = config.get('genWF','batchCPName')

			# get the list of tests
			fmPath = config.get('genWF','batchCPFMPath')
			pairType = 'batch'
			_runPresetMetaPW(outDir,outName,pairType,pwWhich,fmPath,minLogQ=minLogQ,nTopPairs=nTopPairs)

			

		# -QC vs CP
		if config.getboolean('genWF','runQCCP'):
			logging.info("-QC vs. Critical Phenotypes, id bias in unknown latent variables.")
			outName = config.get('genWF','qCCPName')

			# get the list of tests
			fmPath = config.get('genWF','qCCPFMPath')
			pairType = 'QC'
			_runPresetMetaPW(outDir,outName,pairType,pwWhich,fmPath,minLogQ=minLogQ,nTopPairs=nTopPairs)



		# -QC vs Batch

		if config.getboolean('genWF','runQCBatch'):
			logging.info("-QC vs. Batch, id unexpected changes in process.")
			outName = config.get('genWF','qCBatchName')

			# get the list of tests
			fmPath = config.get('genWF','qCBatchFMPath')
			pairType = 'QCBatch'
			_runPresetMetaPW(outDir,outName,pairType,pwWhich,fmPath,minLogQ=minLogQ,nTopPairs=nTopPairs)

	# <--- metadata analysis

	# 

	###### --Run general pairwise 1-- ######
	if config.getboolean('genWF','runGenPW1'): # general pairwise 1 --->
		logging.info("--General PW 1--\n\tRunnin basic pairwise on FM pairs.")
		# get the list of FMs
		genPW1FMList = config.get('genWF','genPW1FMList').split(',')
		# list of paired files must be even
		if (len(genPW1FMList) % 2 == 0):
			for i in range(0,len(genPW1FMList),2):
				# get file name 
				logging.info('-Running Pair '+str(i/2+1))
				logging.info('Test features at '+genPW1FMList[i])
				logging.info('Target features at '+genPW1FMList[i+1])
				pwPath = outDir+'/genPW1_'
				pwPath += genPW1FMList[i][genPW1FMList[i].rfind('/')+1:genPW1FMList[i].rfind('.')]
				pwPath += '_vs_'
				pwPath += genPW1FMList[i+1][genPW1FMList[i+1].rfind('/')+1:genPW1FMList[i+1].rfind('.')]
				pwOutPath = pwPath+'_out.dat'
				pwRepPath = pwPath+'_summary.dat'
				# run the paired FM workflow
				# NOTE: potential mem issue as the summary code in this wf loads the whole pw output as np file object!
				logging.info('Output saved to '+pwOutPath)
				logging.info("Summary of pairwise output saved at {}.".format(pwRepPath))
				_run2FMPW(genPW1FMList[i],genPW1FMList[i+1],pwOutPath,pwRepPath,outDir,pwWhich,samples=samples,minLogQ=minLogQ,nTopPairs=nTopPairs)
				


		else:
			logging.warning('General pairwise 1 NOT run as the FM list does not have complete pairs (odd number).')

	# <---- general pairwise 1

	logging.info("run completed.")
	logging.info("")

if __name__ == '__main__':
	main()
