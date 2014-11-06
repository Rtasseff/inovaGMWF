""" Methods for a standard feature matrix based workflow on INOVA
Methods to create FMs out of source files
(IP) Methods to merge FMs
(IP) Methods to run pairwise analysis

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

Thus:
Allowable mebers = [NB, NB-A, M, F] where all with same family ID 
are mapped to FAM when appropriate.  
"""

import numpy as np



floatFmtStr = '%05.4E'


def _checkPatientID(patientID,studyID='101',allowedSuffix=['FAM']):
	"""Check that ID is consistent with expectations"""
	print patientID
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

def _parseCGSampID(sampID,studyID='101'):
	# assumes <member>(optional,-<member part 2>)-<study id>-<trio id>
	# returns the format <study ID>-<family ID>-<representitive member [NB,NB-A]>
	# memb can be used to check if this is a valid member ie NB,NB-1,M,F 

	if studyID != '101': raise ValueError('ERR:0008, Study not indicated to be 101, logic for 101 maps all features to family centric feature matrix')
	tmp = sampID.strip().split('-')
	membID = -1
	
	# the sample id is diffrent form ours + its more variable (diffren pre and suffixes)
	# map to family ID
	if len(tmp) == 4: 
		memb = tmp[0]+'-'+tmp[1]
		sampIDNew = tmp[2]+'-'+tmp[3]+'-FAM'
		# its possible this parsing does not work, but that is only 
		# in cases where the family member is not part of the 
		# 101 family centric study, and the sample will return 
		# an inappropriate member ID, to be used by caller to ignore sample ID
	else: 
		memb = tmp[0]
		sampIDNew = tmp[1]+'-'+tmp[2]+'-FAM' # if this is wrong it will be caught via memb

	return(sampIDNew, memb)

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


def parseMkGenomeBatchFM(manPath,fmPath,patients,genFeatName="GNMC:BATCH"):
	"""Parse the source files for the 
	genomic batch information, currently
	the manifest file only.
	Then make and save the batch feature matrix.
	manPath 	str, path to input manifest file
	fmPath 	str, path to output genomic batch feature matrix
	patients	list of patients ids as strings 
			related to primary feature matrix
	genFeatName 	str, the general name of the these batch features
	"""
	# get data from the manifest
	# mainifest is small, use whats easiest (in mem, np arrays):
	patients = np.array(patients,dtype=str)
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
	featHeader = np.array(['N:NB:'+genFeatName+':version','N:M:'+genFeatName+':version','N:F:'+genFeatName+':version','N:NB:'+genFeatName+':shipDate','N:M:'+genFeatName+':shipDate','N:F:'+genFeatName+':shipDate'],dtype=str)
	# hard code the possible number of mebers, corrisponding the the feautres above
	nMembers = 3

	nPatients = len(patients)

	# create blank FM
	fm = np.array((np.zeros((nFeat*nMembers,nPatients)) + np.nan),dtype='|S15') # making it an str allows us to control format as we go and use non numerical data types if needed later.

	# loop through all data
	for i in range(n):
		# get sample ID
		sampID, memb = _parseCGSampID(data[i,sampInd])
		# hard code the memb index to conform to order above in FAM centric FM
		membInd = -1
		if memb=='NB' or memb=='NB-A': membInd = 0
		elif memb=='M': membInd = 1
		elif memb=='F': membInd = 2
		# ignore members that do not map to predefined family unit
		if membInd >= 0:
			# version
			iFeat = 0
			row = iFeat * nMembers + membInd
			fm[row,sampID==patients] = '%i'%(_parseVersion(data[i,versionInd]))
			

			# ship date
			iFeat = 1
			row = iFeat * nMembers + membInd
			fm[row,sampID==patients] = '%i'%(_paresShipDate(data[i,dateInd]))

	# save the final matrix
	saveFM(fm,fmPath,featHeader,patients)

	
def saveFM(data,fmPath,featHeader,patients):
	"""Given np arrays, save the feature matrix in 
	standard format.
	"""
	fout = open (fmPath,'w')
	n,m = data.shape
	if n!=len(featHeader):
		raise ValueError('ERR:0009, number of features (rows) and headers not equal')
	if m!=len(patients):
		raise ValueError('ERR:00010, number of patients and col not equal')
	fout.write('.\t')
	fout.write('\t'.join(patients)+'\n')
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
					print 'WARN0002: sample ID '+sampleID+' not found in feature matrix at '+ finName
				if np.sum(sampleID==labels)>1:
					print 'WARN0003: sample ID '+sampleID+' had multiple entries in feature matrix at '+ finName+'. Using only the first.'

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
					
			

			





def main():
	fmInPath = '/Users/RTasseff/Projects/INOVA/DF5/basePipeLine/2014_10_17_genomic_clinical_hilevel_RT.fm'
	fmPath = '/Users/RTasseff/Projects/INOVA/DF5/basePipeLine/test.fm'
	manPath = '/Users/RTasseff/Projects/INOVA/DF5/basePipeLine/P286_Data_Delivery_TOTAL_021214.txt'
	trioList = getPatientOrder(fmInPath)
	parseMkGenomeBatchFM(manPath,fmPath,trioList)
	catFM(['/Users/RTasseff/Projects/INOVA/DF5/basePipeLine/test.fm','/Users/RTasseff/Projects/INOVA/DF5/basePipeLine/DF5_Methylation_Batch_Vars.txt'],trioList,'/Users/RTasseff/Projects/INOVA/DF5/basePipeLine/test_2.fm',studyID='101')


if __name__ == '__main__':
	main()
