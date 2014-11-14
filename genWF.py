""" Methods for a standard feature matrix based workflow on INOVA
Methods to create FMs out of source files
- Methods to merge FMs, check consistancy and filter for stats purposes
- Methods to run pairwise analysis and sumarize resutls (currently for medium to small runs only)

Desinged around study 101,
Some limitations on this realte to primary pheontype being defined by NB 
as in NB was PTB. For some sources data on all family members is avalibel
but in others this is not true (PTB on NB only, molecualr data on M only)
so we have organized this around the family unit (suffix -FAM).  
In cases with twins, only one is chosen as representtive.  
This avoids exact duplicated data (M molecualr data)
and near exact reprelicated data (similarity in genomes).  
Basically to allow for the indepedence assumption of most analysis. 
We note that this limitation is introduced for systimatic analysis 
of feature matrix like input data.
Thus, allowable mebers = [NB, NB-A, NB-B, M, F] where all with same family ID 
are mapped to FAM when appropriate based on a predefined list of representitive newborns.  
"""

import numpy as np
import subprocess
import time
import statsUtil as su
import logging
import sys

floatFmtStr = '%05.4E'
nanValues = ['NA','NaN','na','nan']

disc="Python workflow for standard genomic/molecular data analysis, INOVA 101."
version='0.1.0'

def _checkPatientID(patientID,studyID='101',allowedSuffix=['FAM']):
	"""Check that ID is consistent with expectations"""
	tmp = patientID.split('-')
	if tmp[0] != studyID:
		raise ValueError ("ERR:0001 unexpected study id "+str(tmp[0]))
	elif '-'.join(tmp[2:]) not in allowedSuffix: 
		raise ValueError ("ERR:0002 unexpected member id "+'-'.join(tmp[2:])+'. Suffixes limited to '+str(allowedSuffix))


def getPatientOrder(fmPath,studyID='101',allowedSuffix=['FAM']):
	"""Get the smaple order from a properly fromated, pre-existing
	feature matrix.

	fmPath	str, path to standard featuremetrix
	"""
	fin = open( fmPath )
	line = fin.next() # first line is header with info
	patients = line.strip().split('\t')[1:] # first position is '.' place holder
	for patient in patients:
		_checkPatientID(patient,studyID=studyID,allowedSuffix=allowedSuffix)	
	return patients

def getRepNBList(nbListPath,):
	"""Create a list (np str array) that contains
	all usable newborns for family centric analysis.
	Used in maping individuals to -FAM
	for family centric analsyis.
	Tab sep input list should have all and only
	newborns used to represent a family
	unit. An error will be returned if
	a family id xxx in 101-xxx-NB is 
	included twice.

	The list will be used for simplicity;
	effectivly -F -M and -NB are 
	just replaced with -FAM. However,
	for multi births only the one
	in the list will be used in the 
	-FAM family centric analysis.
	"""
	nbList = np.loadtxt(nbListPath,dtype=str)
	n = len(nbList)
	used = np.zeros(n,dtype=int)
	# check
	for i in range(n):
		nb = nbList[i]
		if not nb.split('-')[2]=='NB':
			raise ValueError('ERR0012:non newborn id '+nb+' found in nb list at '+nbListPath)
		famID = int(nb.split('-')[1])
		if np.any(famID==used):
			raise ValueError('ERR0013:duplicate family ids '+str(famID)+' found in nb list at '+nbListPath)
		used[i] = famID
	return(nbList)

def _itmiID2isbID(itmiID):
	"""ITMI puts member ids as prefixes, but ISB puts
	them as sufixes.  This will echange that order.
	Note that nan will be returned if member id is 
	not yet accounted for.  For example sample
	PO-101-657-22 has PO member id that is 
	not accounted for in any logic so far or
	analysis done. Note that the 3 sets 
	of number complicates maters with this ID.
	"""
	isbID = 'na'
	tmp = itmiID.split('-')
	# hard coded allowed prefixes 
	if tmp[0] in ['F','M','NB']:
		if len(tmp)==3:
			isbID = tmp[1]+'-'+tmp[2]+'-'+tmp[0]
		if len(tmp)==4:
			isbID = tmp[2]+'-'+tmp[3]+'-'+tmp[0]+'-'+tmp[1]
	return(isbID)

def _mapInd2FAM(indID,nbList):
	"""effectivly -F and -M are 
	just replaced with -FAM. However,
	for multi births only the one
	in the list will be used in the 
	-FAM family centric analysis.
	The nbList will be used to deterimine
	if the passed individual is nb and if so
	is the representitve nb from the list.
	Returns the new member id and new -FAM id.
	Non representitve nb's or unrecognized 
	member ids will return na.
	indID	str, isb fomrated id for individual smaple
	nbList	np str array of usable representitive new borns
	"""
	membID = 'na'
	famID = 'na'
	tmp = indID.split('-')
	if len(tmp)==3 and tmp[2] in ['F','M']:
		membID = tmp[2]
		famID = tmp[0]+'-'+tmp[1]+'-FAM'
	# only checking for multi births which have 4 '-'
	elif np.any(indID==nbList):
		famID = tmp[0]+'-'+tmp[1]+'-FAM'
		membID = 'NB'
	return(famID,membID)
		






	


def _parseCGSampID(itmiID,nbList,studyID='101'):
	"""Parse the CG file sample id (equivlant to itmi id)
	into a usable isb id.  Currently coded for family centric analysis 
	and returns the family id and seperatly the member (if id is usable, na otherwise).
	itmiID	str, itmi id from CG file
	nbList	np str array of usable representitive new borns
	studyID	str the expected study id, for checking
	"""
	# assumes <member>(optional,-<member part 2>)-<study id>-<trio id>
	# returns the format <study ID>-<family ID>-<representitive member [NB,NB-A]>
	# memb can be used to check if this is a valid member ie NB,NB-1,M,F 

	if studyID != '101': raise ValueError('ERR:0008, Study not indicated to be 101, logic for 101 maps all features to family centric feature matrix')
	isbID = _itmiID2isbID(itmiID)

	return(_mapInd2FAM(isbID,nbList))

def _parseVersion(versionStr):
	# parse the version string form the CG manifest file
	tmp = versionStr.strip().split('.')
	if len(tmp) != 3:
		raise ValueError("ERR:0003 unexpected version "+versionStr)
	return int(tmp[-1])

def _paresShipDate(shipDateStr):
	# parse the shipped date string form the CG manifest file
	tmp = shipDateStr.strip().split('/')
	if len(tmp) != 3:
		raise ValueError("ERR:0004 unexpected ship date "+shipDateStr)
	tmp = tmp[0]+tmp[1]+tmp[2]
	return int(tmp)


def parseMkGenomeBatchFM(manPath,fmPath,samples,nbList,studyID='101',genFeatName="GNMC:BATCH"):
	"""Parse the source files for the 
	genomic batch information, currently
	the manifest file only.
	Then make and save the batch feature matrix.
	Currently creating a family centric feautre matrix,
	relevent to study 101, so features not mapable to a 
	family unit (ie -FAM) will be ignored.
	manPath 	str, path to input manifest file
	fmPath 	str, path to output genomic batch feature matrix
	samples	list of sample ids as strings 
			related to primary feature matrix
	nbList	np str array of usable representitive new borns
	studyID	str the expected study id, for checking

	genFeatName 	str, the general name of the these batch features
	"""
	# get data from the manifest
	# mainifest is small, use whats easiest (in mem, np arrays):
	samples = np.array(samples,dtype=str)
	data = np.loadtxt(manPath,dtype=str,delimiter='\t')
	labels = data[0]
	data = data[1:]
	n,m = data.shape
	# create matrix, numerical features, order and version
	nFeat = 2 # version, ship data
	# get index for key col
	if np.sum(labels=='Cust. Subject ID')==1:
		sampInd = np.arange(m)[labels=='Cust. Subject ID'][0]
	else:
		raise ValueError('ERR:0005, Sample ID not found')

	if np.sum(labels=='Assembly Version')==1:
		versionInd = np.arange(m)[labels=='Assembly Version'][0]
	else:
		raise ValueError('ERR:0006, Assembly Version not found')

	if np.sum(labels=='Shipped Date')==1:
		dateInd = np.arange(m)[labels=='Shipped Date'][0]
	else:
		raise ValueError('ERR:0007, Shipped Date not found')

	# hard code the positon and name of the features...
	featHeader = np.array(['C:NB:'+genFeatName+':version','C:M:'+genFeatName+':version','C:F:'+genFeatName+':version','N:NB:'+genFeatName+':shipDate','N:M:'+genFeatName+':shipDate','N:F:'+genFeatName+':shipDate'],dtype=str)
	# hard code the possible number of mebers, corrisponding the the feautres above
	nMembers = 3

	nSamples = len(samples)

	# create blank FM
	fm = np.array((np.zeros((nFeat*nMembers,nSamples)) + np.nan),dtype='|S15') # making it an str allows us to control format as we go and use non numerical data types if needed later.

	# loop through all data
	for i in range(n):
		# get sample ID
		sampID, memb = _parseCGSampID(data[i,sampInd],nbList,studyID)
		# hard code the memb index to conform to order above in FAM centric FM
		membInd = -1
		if memb=='NB': membInd = 0
		elif memb=='M': membInd = 1
		elif memb=='F': membInd = 2
		# ignore members that do not map to predefined family unit
		if membInd >= 0:
			# version
			iFeat = 0
			row = iFeat * nMembers + membInd
			fm[row,sampID==samples] = '%i'%(_parseVersion(data[i,versionInd]))
			

			# ship date
			iFeat = 1
			row = iFeat * nMembers + membInd
			fm[row,sampID==samples] = '%i'%(_paresShipDate(data[i,dateInd]))

	# save the final matrix
	saveFM(fm,fmPath,featHeader,samples)

	
def saveFM(data,fmPath,featHeader,samples):
	"""Given np arrays, save the feature matrix in 
	standard format.
	"""
	fout = open (fmPath,'w')
	n,m = data.shape
	if n!=len(featHeader):
		raise ValueError('ERR:00009, number of features (rows) and headers not equal')
	if m!=len(samples):
		raise ValueError('ERR:00010, number of patients and col not equal')
	fout.write('.\t')
	fout.write('\t'.join(samples)+'\n')
	for i in range(n):
		fout.write(featHeader[i]+'\t'+'\t'.join(np.array(data[i],dtype=str))+'\n')

	fout.close()


#def _catUnorderedFM

def _catOrderedFM(filenames,outfile,skipHeader=True):
	for fname in filenames:
		with open(fname) as infile:
			if skilpHeader:infile.next()
			for line in infile:
				outfile.write(line)

			
		


def filterFM(fmInPath,fmOutPath,maxMiss=.9,minObs=5):
	"""Filter a given feature matrix for statistically 
	poor feauters for use in pairwise code.
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
			logging.info('removing '+':'.join(fname)+', in correct number of columns.')
		if np.sum(np.array(tmp[1:],dtype=str)=='NA')>n*maxMiss:
			skip = True
			logging.info('removing '+':'.join(fname)+', too many missing values.')
		if fname[0]=='C':
			# check for size
			unique = list(set(tmp[1:]))
			if len(unique)>30:
				skip = True
				logging.info('removing '+':'.join(fname)+',too many catigories.')
		if fname[0]=='B':
			# check for size
			unique = list(set(tmp[1:]))
			if ('NA' in unique and len(unique)>3) or ('NA' not in unique and len(unique)>2):
				skip = True
				logging.info('removing '+':'.join(fname)+',too many catigories for binary.')

		if fname[0]=='B' or fname[0]=='C':
			unique = list(set(tmp[1:]))
			values = np.array(tmp[1:],dtype=str)
			for label in unique:
				if np.sum(label==values)<minObs:
					skip = True
					logging.info('removing '+':'.join(fname)+', not enough samples for a group.')
					break

		if not skip:
			fout.write(':'.join(fname)+'\t'+'\t'.join(tmp[1:])+'\n')
	fout.close()
	fin.close()

	





def catFM(fMNames,sampleIDs,foutPath,studyID='101',allowedSuffix=['FAM']):
	"""Concatinate multiple feature matrices together.
	Assumes first line is header for sample lables,
	and use sampleID to determine columns and their order.
	nan will be used to fill in missing samples,
	sample in sampleID but not in a FM header.
	fMnames	list of str paths for feature matrices to be cated
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
					logging.info('WARN0002: sample ID '+sampleID+' not found in feature matrix at '+ finName)
				if np.sum(sampleID==labels)>1:
					logging.info('WARN0003: sample ID '+sampleID+' had multiple entries in feature matrix at '+ finName+'. Using only the first.')

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
	"""Going thourgh feature names in fmPath,
	record names that conatin strings in the list 
	lookfor in the corrisponding feild index, feildInd.
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
	QC 	compare QC features to critical phenotypes to detect bias in counfounding factors
	QCBatch	compare batch features to QC features to detect process changes
	"""
	if pairType=='batch':
		feildInd = [3,3]
		lookfor = ['BATCH','Critical_Phenotype']
		nameLists =  _getFeatureNamesByGroup(fmPath,lookfor,feildInd)
		_mkPairFile(nameLists[0],nameLists[1],foutPath)
	elif pairType=='QC':
		feildInd = [3,3]
		lookfor = ['QC','Critical_Phenotype']
		nameLists =  _getFeatureNamesByGroup(fmPath,lookfor,feildInd)
		_mkPairFile(nameLists[0],nameLists[1],foutPath)
	elif pairType=='QCBatch':
		feildInd = [3,3]
		lookfor = ['QC','BATCH']
		nameLists =  _getFeatureNamesByGroup(fmPath,lookfor,feildInd)
		_mkPairFile(nameLists[0],nameLists[1],foutPath)

	else:
		raise ValueError('ERR:00014, pairType value is not in preset list.')




			

def runPairwise(fMPath,outDir,outPath,pwWhich='/titan/cancerregulome8/TCGA/scripts/pairwise-2.0.0-current',pairFile=''):
	"""Run the pairwise code found at pwWhich on feature matrix at
	fMName and save the output to outName.
	"""
	# construct the call to the pairwise script
	call = pwWhich
	# check to see if we have a file to indicate the pairs to calculate
	if pairFile!='':
		call = call+' --by-name '+pairFile

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


def writePWSumLong(pwPath, repPath, minLogQ=2.0):
	"""Takes the pairwise output at pwPath 
	and parses it into a suammry report at repPath
	by filtering out paris with a log_10(FDR Q)<minLogQ
	and by dividing the ouput up by data source.

	Assumes pairwise output format from pairwise-2.0.0
	& column 0 represent covariates of intrest (enforeced in workflow by suppling pairFile, see runPairWise)
	& column 1 represents target of intrest (enforeced in workflow by suppling pairFile, see runPairWise)
	& orginizing by data source which is index 2 (start at 0) feild of the feautre name for the covariates.
	& assumes pairwise output can fit in np array
	"""
	pw = np.loadtxt(pwPath,dtype=str,delimiter='\t')
	# get data sources
	n,m = pw.shape
	dataSources = np.array(np.zeros(n),dtype='|S15')
	for i in range(n):
		dataSources[i] = pw[i,0].split(':')[2]
	# get all data sources
	uniqueSources = list(set(dataSources))
	# prepare to write report 
	with open(repPath,'w') as rep:
		rep.write('Suammary output file for pairwise analysis run: '+ time.strftime("%c")+'.\n')
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
			_,q,_ = su.fdr_bh(p)
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
			else:
				rep.write('\t------None > min log(p)--------\n')



def main():

	clinFM = '/titan/ITMI1/projects/gamcop/data/featureMatrices/data_CLIN_20141106.fm'
	gnmcBatchFM = '/titan/ITMI1/projects/gamcop/data/featureMatrices/metadata_GNMC_20141111.fm'
	rnaBatchFM = '/titan/ITMI1/projects/gamcop/data/featureMatrices/metadata_RNASeq_20141106.fm'
	miBatchFM = '/titan/ITMI1/projects/gamcop/data/featureMatrices/metadata_MIRSeq_20141111.fm'
	wrkDir = '/users/rtasseff/inova/inovaMWF'
	logName = 'log.out'
	nbListPath = wrkDir+'/repNBList.tsv'
	outFM = wrkDir+'/test.fm'
	manPath = wrkDir+'/P286_Data_Delivery_TOTAL_021214.txt'
	pwWhich='/titan/cancerregulome8/TCGA/scripts/pairwise-2.0.0-current'

	# --setup logger and decrease level:
	logging.getLogger('').handlers = []
	logging.basicConfig(filename=wrkDir+'/'+logName, level=logging.INFO, format='%(asctime)s %(message)s')

	# --record some basic information:
	logging.info("--------------------------------------------------------------")
	logging.info("Running {}, {}, version={}...".format(sys.argv[0],disc,version))
	logging.info("--------------------------------------------------------------")

	# --general setup
	logging.info("--Getting some general information.")
	# need the list of newborns for mapping to correct family
	logging.info("List of representative newborns at {}.".format(nbListPath))
	nbList = getRepNBList(nbListPath)

	# need the sample list 
	logging.info("Standardizing to sample list in FM at {}.".format(clinFM))
	samples = getPatientOrder(clinFM)

	# --parse the manifest file for gnmc batch info
	logging.info("--Getting GNMC BATCH FM.")

	logging.info("Using CG manifest at {}.".format(manPath))
	parseMkGenomeBatchFM(manPath,gnmcBatchFM,samples,nbList)
	logging.info("GNMC BATCH FM saved at {}.".format(gnmcBatchFM))

	# --Run metadata association tests
	logging.info("--Running standard association tests on metadata.")

	# feature matrices to use:
	metaFMs = [miBatchFM,gnmcBatchFM,rnaBatchFM,clinFM]
	logging.info("Considering metadata data stored at:")
	for fm in metaFMs:
		logging.info("\t{}.".format(fm))

	# cat matrices together
	catFM(metaFMs,samples,outFM,studyID='101')
	logging.info("Metadata stored in single FM at {}.".format(outFM))

	# filter for pairwise
	logging.info("Filtering metadata for basic statistical tests.")
	filteredFM = wrkDir+'/test_tmp.fm'
	filterFM(outFM,filteredFM)
	logging.info("Filtered metadata at {}.".format(filteredFM))



	# -batch vs CP
	logging.info("-Batch vs. Critical Phenotypes, id bias in sampling.")

	# get the list of tests
	pairFile = wrkDir+'/pairFileTestBatch.dat'
	mkPSPairFile(filteredFM,pairFile,pairType='batch')
	logging.info("Pairwise tests to run indicated in {}.".format(pairFile))

	# run pairwise
	pwOutName = 'pair_test_batch'
	logging.info("Running pairwise code, {}.".format(pwWhich))
	pwOutPath = wrkDir+'/'+pwOutName+'_out.dat'
	runPairwise(filteredFM,wrkDir,pwOutPath,pwWhich=pwWhich,pairFile=pairFile)
	logging.info("Full pairwise output saved at {}.".format(pwOutPath))
	pwSumPath = wrkDir+'/'+pwOutName+'_summary.dat'
	writePWSumLong(pwOutPath,pwSumPath)
	logging.info("Summary of pairwise output saved at {}.".format(pwSumPath))



	# -QC vs CP
	logging.info("-QC vs. Critical Phenotypes, id bias in unknown latent variables.")

	# get the list of tests
	pairFile = wrkDir+'/pairFileTestQC.dat'
	mkPSPairFile(filteredFM,pairFile,pairType='QC')
	logging.info("Pairwise tests to run indicated in {}.".format(pairFile))

	# run pairwise
	pwOutName = 'pair_test_qc'
	logging.info("Running pairwise code, {}.".format(pwWhich))
	pwOutPath = wrkDir+'/'+pwOutName+'_out.dat'
	runPairwise(filteredFM,wrkDir,pwOutPath,pwWhich=pwWhich,pairFile=pairFile)
	logging.info("Full pairwise output saved at {}.".format(pwOutPath))
	pwSumPath = wrkDir+'/'+pwOutName+'_summary.dat'
	writePWSumLong(pwOutPath,pwSumPath)
	logging.info("Summary of pairwise output saved at {}.".format(pwSumPath))


	# -QC vs Batch
	logging.info("-QC vs. Critical Phenotypes, id unexpected changes in process.")

	# get the list of tests
	pairFile = wrkDir+'/pairFileTestQCBatch.dat'
	mkPSPairFile(filteredFM,pairFile,pairType='QCBatch')
	logging.info("Pairwise tests to run indicated in {}.".format(pairFile))

	# run pairwise
	pwOutName = 'pair_test_qcbatch'
	logging.info("Running pairwise code, {}.".format(pwWhich))
	pwOutPath = wrkDir+'/'+pwOutName+'_out.dat'
	runPairwise(filteredFM,wrkDir,pwOutPath,pwWhich=pwWhich,pairFile=pairFile)
	logging.info("Full pairwise output saved at {}.".format(pwOutPath))
	pwSumPath = wrkDir+'/'+pwOutName+'_summary.dat'
	writePWSumLong(pwOutPath,pwSumPath)
	logging.info("Summary of pairwise output saved at {}.".format(pwSumPath))





	logging.info("run completed.")
	logging.info("")

if __name__ == '__main__':
	main()
