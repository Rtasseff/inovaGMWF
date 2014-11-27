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
import sys
import gzip

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


def getQCFM_VCF(vcfPath,fmOutPath,genFeatName="GNMC:QC",sampStartInd = 9):
	"""Get the predetermined QC features from vcf at vcfPath 
	(use gzip if path ends in .gz) and create and individual
	sample feature matrix at fmOutPath. genFeatName is a str
	of the general name of the these QC features
	<Data Source>:<Feature Class>
	current features:
	MIE_count - total number of annotated MIEs over genome
	VQLow_count - total number of annotated VQLow calls over genome
	Missing_count - total number of missing calls over genome
	homoz_count - total number of homozygous calls 
	heteroz_count - total number of heterozygous calls
	Note are parsing assumptions:
	- last comment line is sample line
	- sample line first entry #CHROM
	- sample line has sampStartInd entries before the IDs start
	- col of variant line correspond with sample line
	- once variant lines start they continue to eof
	- vcf is tab sep

	"""

	# to ensure consistency
	mieInd = 0
	vqlowInd = 1
	missInd = 2
	homoInd = 3
	hetroInd = 4
	nFeat = 5

	featName = np.array(np.zeros(nFeat),dtype='|S60')
	featName[miwInd] = 'N:'+genFeatName+':MIE_count'
	featName[vqlowInd] = 'N:'+genFeatName+':VQLow_count'
	featName[missInd] = 'N:'+genFeatName+':miss_count'
	featName[homoInd] = 'N:'+genFeatName+':homoz_count'
	featName[hetroInd] = 'N:'+genFeatName+':heteroz_count'

	
	
	if vcfPath[-2:]=='gz':
		vcf = gzip.open(vcfPath,"r")
	else:
		vcf = open(vcfPath,"r")
	# keep track of position, top header
	readLine = False
	for line in vcf: ### in vcf -->
		col = line.strip().split('\t')
		if not readLine:
			# check if this is needed comment 
			if col[0]=='#CHROM':
				sampleList = np.array(col[sampStartInd:],dtype=str)
				nSamp = len(sampleList)
				data = np.zeros((nFeat,nSamp))
				# verify these are known sample ids formats
				known=True
				for samp in sampleList:
					idType = genUtil.idType(samp)
					if idType=='na':
						raise ValueError("ERR:0008, found an unknown sample ID format in "+vcfPath+", vcf may be incorrect, format may need to be added to known formats, or sampStartInd may need to be changed.")

		else:  #### in readable section, variant line -->
			calls = col[sampStartInd:]
			# check len
			if len(calls)!=nSamp:
				raise ValueError("ERR:0009, number of calls and samples does not match in "+vcfPath
			for i in range(nSamp):
				err = False
				if calls[i].find(".") <> -1: 
					data[missInd,i] += 1
					err = True
				if calls[i].find("VQLOW") <> -1: 
					data[vqlowInd,i] += 1
					err = True
				if calls[i].find("MIE") <> -1: 
					data[mieInd,i] += 1
					err = True
				if not err:
					tmp = calls[i].split('/')
					if len(tmp)==2:
						if tmp[0]==tmp[1]:data[homoInd,i] += 1
						else: data[hetroInd,i] += 1

				


		#### <-- in readable section, variant line
	### <-- in vcf

	# done with vcf save feature matrix
	saveFM(data,fmOutPath,featName,sampleList)


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
