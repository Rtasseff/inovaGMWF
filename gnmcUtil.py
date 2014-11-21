"""
Python2.7 module created for INOVA Project:
Utility methods for working for the genomic data,
including the VCF and other source materials from
CG and ITMI, as well as internal tools like results
from RK's split-transcripts.

RAT 20141121
"""

import numpy as np
import genUtil 


import warnings


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


def parseMkGenomeBatchFM(manPath,fmPath,samples,genFeatName="GNMC:BATCH"):
	"""Parse the source files for the 
	genomic batch information, currently
	the manifest file only.
	Then make and save the batch feature matrix for the individuals.
	Header and samples extracted will be determined by 'samples'.
	Converts all source sample ids to the same id format as 'samples',
	which is determined automatically.
	'samples' is expected to have this format already.
	manPath 	str, path to input manifest file
	fmPath 		str, path to output genomic batch feature matrix
	samples		if list, used as actual samples with ISB ids as strings 
			if a single str then we assume its a path to tab 
			sep file containing the sample ids 
	genFeatName 	str, the general name of the these batch features
			<Data Source>:<Feature Class>
	"""
	# get data from the manifest
	# mainifest is small, use whats easiest (in memory, np arrays):

	# get samples
	if type(samples)==str:samples = np.loadtxt(samples,dtype=str,delimiter='\t')
	else:samples = np.array(samples,dtype=str)


	# find out format sample ID type:
	sampIDType = genUtil.idListType(samples)
	if type(sampIDType)==list:sampIDType=sampIDType[0]
	if sampIDType=='na':
		raise ValueError ('ERR:0002 the passed in sample IDs do not match any known format')
	nSamples = len(samples)


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
	featHeader = np.array(['C:'+genFeatName+':version','N:'+genFeatName+':shipDate'],dtype=str)
	if len(featHeader)!=nFeat: raise RuntimeError ("RTE0001: Hard coded options for number of features is wrong")


	# create blank FM
	fm = np.array((np.zeros((nFeat,nSamples)) + np.nan),dtype='|S15') # making it an str allows us to control format as we go and use non numerical data types if needed later.

	# loop through all data
	for i in range(n):
		# get sample ID
		if sampIDType=='isb':sampID = genUtil.cgID2isb(data[i,sampInd])
		elif sampIDType=='itmi':sampID = genUtil.cgID2itmi(data[i,sampInd])
		elif sampIDType=='cg':sampID = data[i,sampInd]
		else: raise RuntimeError ('RTE:0003 for some reason you are getting back an unknown sample id type indicator')

		# version
		row = 0 
		fm[row,sampID==samples] = '%i'%(_parseVersion(data[i,versionInd]))
		
		# ship date
		row = 1
		fm[row,sampID==samples] = '%i'%(_paresShipDate(data[i,dateInd]))

	# save the final matrix
	genUtil.saveFM(fm,fmPath,featHeader,samples)


def main():
	# test functions
	wrkDir = '/Users/RTasseff/Projects/INOVA/DF5/basePipeLine'
	manPath = wrkDir+'/P286_Data_Delivery_TOTAL_021214.txt'
	fmPath = wrkDir+'/test_IND_20141121.fm'
	samples = wrkDir+'/sampIDList_ITMI_VCF_DF5.tsv'
	famFMPath = wrkDir+'/test_FAM_20141121.fm'
	famSample = wrkDir+'/famID_DF5.tsv'
	nbList = wrkDir+'/repNBList.tsv'

	parseMkGenomeBatchFM(manPath,fmPath,samples)

	genUtil.indFM2FamFM(fmPath,famFMPath,nbList,famSample)

if __name__ == '__main__':
	main()
