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
import struct
import warnings
import string

def _getNumericCode(genotype):
	"""Convert vcf call to number.
	Genotype Code: Ref homozygous is 0, heterozygous is 1, 
	non-ref homozygous is 2, missing genotype is NA.
	ported from vcfToFMx.py vdhankan@systemsbiology.org 20150304.
	"""
	if genotype == "0/0" or genotype == "0":
		return "0"
	elif genotype == "0/1" or genotype == "1/0" or genotype == "1":
 		return "1"
	elif genotype == "1/1":
		return "2"
	elif genotype == "./." or genotype == ".":
		return "NA"

def vcf2FM(vcfFilename,outFMFilename):
	"""Convert vcf file (path=vcfFilename) to a FM (path=outFMFilename).
	Genotype Code: Ref homozygous is 0, heterozygous is 1, 
	non-ref homozygous is 2, missing genotype is NA.
	ported from vcfToFMx.py vdhankan@systemsbiology.org 20150304.
	"""
	vcfFile = gzip.open(vcfFilename,"r")
        outFile = gzip.open(outFMFilename,"w")
        for line in vcfFile:
                columns = line.strip().split()
                outputColumns = []
                if line.startswith('chr'):
                        featureName = "C:GNMC:Data:"+columns[0]+"_"+columns[1]
                        outputColumns.append(featureName)
                        outputColumns.extend([getNumericCode(x) for x in columns[9:len(columns)-13]]) #last 13 samples are POs
                        print >>outFile, '\t'.join(outputColumns)
                elif line.startswith('#CHROM'):
                        featureName = "."
                        outputColumns.append(featureName)
                        outputColumns.extend(columns[9:len(columns)-13])
                        print >>outFile, '\t'.join(outputColumns)
        vcfFile.close()

def vcf2FM_regions():
	"""This script generates genomic feature matrix to use for pairwise analysis.
	INPUT: Reads in 634 vcf regions one at a time. The path to these region files is hard coded to be /bigdata0/users/dmauldin/DF5/CMS-EH-MV/
	OUTPUT: One genomic feature matrix per vcf region. The output FMs are written in the working directory.
	Genotype Code: Ref homozygous is 0, heterozygous is 1, non-ref homozygous is 2, missing genotype is NA.
	ported from vcfToFMx.py vdhankan@systemsbiology.org 20150304
	"""
	warnings.warn('While the original code (vcfToFMx.py) has been tested, this integrated method has not been; if run is successful please remove this warning')
	outDir = '/isb/rtasseff/data/featureMatrices/data_VCF_FM_regions_20150304'

	for i in range(1,635):
		vcfFilename = "/bigdata0/users/dmauldin/DF5/CMS-EH-MV/"+str(i)+".Filtered.vcf.gz"
		outFMFilename = outDir+'/'+str(i)+".Filtered.fm.gz"
		vcf2FM(vcfFilename,outFMFilename)

	

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
	featName[mieInd] = 'N:'+genFeatName+':MIE_count'
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
				readLine = True
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
				raise ValueError("ERR:0009, number of calls and samples does not match in "+vcfPath)

			#print calls[0]
			for i in range(nSamp):
				err = False
				
				if calls[i].find("VQLOW") <> -1: 
					data[vqlowInd,i] += 1
					err = True
				if calls[i].find("MIE") <> -1: 
					data[mieInd,i] += 1
					err = True

				callTmp = calls[i].split(':')[0]
				if callTmp.find(".") <> -1: 
					data[missInd,i] += 1
					err = True
				if not err:
					tmp = callTmp.split('/')
					if len(tmp)==2:
						if tmp[0]==tmp[1]:data[homoInd,i] += 1
						else: data[hetroInd,i] += 1

				


		#### <-- in readable section, variant line
	### <-- in vcf

	# done with vcf save feature matrix
	genUtil.saveFM(data,fmOutPath,featName,sampleList)

def summQCRegionFM(qcRegionDir,qcRegionName,numQCRegions,outQCFMPath):
	"""Used to sum up all individual qc regions into a single FM,
	very adhoc logic here:
	loops though all regions to sum up a single FM
	for QC variables, assumes that all files are 
	in the same dir with <qcRegionName>_i.fm 
	where i = 1,2...,numQCRegions
	"""

	for i in range(1,numQCRegions+1):
		try:
			fmName = qcRegionName+'_'+str(i)+'.fm'
			fm = np.loadtxt(qcRegionDir+'/'+fmName,dtype=str,delimiter='\t')
			samples_tmp = list(fm[0,1:])
			features_tmp = list(fm[1:,0])
			if i==1:
				samples = samples_tmp
				features = features_tmp
				data = np.array(fm[1:,1:],dtype=float)
			else:
				if (features_tmp != features) or (samples_tmp != samples):
					raise ValueError('FM labels not the same in '+qcRegionName+'_'+str(i)+'.fm when adding up QC FMs')
				data += np.array(fm[1:,1:],dtype=float)
		except:
			warnings.warn("Could not add QC region. File not found: "+fmName)
		
	genUtil.saveFM(data,outQCFMPath,features,samples)

def readBiMat(path,nRow):
	"""Reads a bianary format and spits out a matrix.
	assumes 8 bit unsigned char, will put it 
	an np array by dividing into nRows
	type = unsigned 8 bit integer
	Note that all '0' values refer to nan as
	is the convention in split-transcripts.
	"""
	# open and read the bianary file
	tmp = open(path,'rb').next()
	# one byte for each char 
	nData = len(tmp)
	# get col
	nCol = nData/nRow 

	# need to setup format for conversion
	fmt = str(nData)+'B'
	# save data as np array, uint8
	data = np.array(struct.unpack(fmt,tmp),dtype='uint8')
	# put into matrix 
	X = np.reshape(data,(nRow,nCol))
	
	return(X)

class RegionManifestReader:
	"""Reads the file at the given path as the mainfest of the 
	regions that were identifed by split-transcripts.  
	Can iterate through using next(), returns a region.
	Assumes the paths in the mainfest for the data are 
	relative to the same dir as the manifest,
	otherwise you can specify the relative dir path
	on initilization.
	Must close the reder like a file object when done.
	"""
	# setting for parsing the manifest file from split-transcripts
	delim = '\t'
	commentStr = '#'
	dataDirPrePath =''
	# index for the col of each value
	IND_NAME = 0
	IND_CHROM = 1
	IND_START = 2
	IND_STOP = 3
	IND_NUM = 4
	IND_DATA_PATH = 5
	def __init__(self,path,dataDirPrePath=''):
		self.fin=open(path)
		# the paths in the manifest are relative paths, lets deal with that here
		# check to see if a dataDirPrePath is specified or needs to be parsed
		if dataDirPrePath=='':
			tmp = path.split('/')
			if len(tmp)>1:self.dataDirPrePath = '/'.join(tmp[:-1])+'/'
		else:
			if dataDirPrePath[-1]!='/':dataDirPrePath=dataDirPrePath+'/'
			self.dataDirPrePath = dataDirPrePath
	def next(self):
		# search for the next readable region 
		regionLine = False
		while not regionLine:
			line = self.fin.next()
			if line[0] != self.commentStr:
				values = line.strip().split(self.delim)
				if len(values)==6:regionLine=True
		# create and pass back the region
		name = values[self.IND_NAME]
		chrom = values[self.IND_CHROM]
		startPos= values[self.IND_START]
		stopPos = values[self.IND_STOP]
		nVar = int(values[self.IND_NUM])
		dataPath = self.dataDirPrePath+values[self.IND_DATA_PATH]
		region = Region(name,chrom,startPos,stopPos,nVar,dataPath)
		return region


	def close(self):
		self.fin.close()

	def __iter__(self):
		return self



class Region:
	"""Represents a region with all relevent information
	including a funciton to load the corrisponding region data.
	"""
	name = ''
	chrom = ''
	startPos = ''
	stopPos = ''
	nVar = ''
	dataPath = ''
	def __init__(self,name,chrom,startPos,stopPos,nVar,dataPath):
		self.name = name
		self.chrom = chrom
		self.startPos = startPos
		self.stopPos = stopPos
		self.nVar = nVar
		self.dataPath = dataPath

	def getData(self):
		"""Get the corresponding data matrix for this transcript 
		from the binary matrix.
		"""
		X = readBiMat(self.dataPath,self.nVar)
		return X



def _removeMissing(X,maxMiss):
	"""Return a boolean 
	vector to indicate which samples
	to keep from X based on missing calls < maxMiss.
	Assumes a call of 0 is missing (convention 
	in split-transcripts)
	"""
	n,m = X.shape
	keep = np.array(np.ones(m/2),dtype=bool)
	missed = np.sum(X==0,0)
	for i in range(0,m,2):
		if ((missed[i]+missed[i+1])/(n*2.0)) > maxMiss: keep[i/2]=False
	return keep
		
def _isMinor(x,miAFrac=.49):
	"""Given a vector of calls, x, return an
	int vector in corresponding order that determines if 
	calls are minor (1) or not (0).
	Calls are considered to be minor alleles iff their 
	call fraction in this population is < miAFrac. 
	"""
	n = len(x)
	unique = list(set(x))
	minor = np.array(np.zeros(n),dtype=int)
	for value in unique:
		# ignore zeros, used as missing by convention
		if value > 0:
			# calculate call fraction
			frac = np.sum(x==value)/float(n)
			if frac<miAFrac:
				minor[x==value]=1
	return minor

def _one2twoCol(x):
	"""Since samples are represented by two columns 
	in our genotype matrix, it is sometimes 
	needed to create a corresponding 2x vector
	such that each adjacent column is a copy of
	the original.  This will provide the same value
 	for the same samples over both columns in the 
	genotype matrix.
	"""
	y = np.array(np.ones(len(x)*2),dtype=x.dtype)
	y[::2] = x
	y[1::2] = x
	return y


def getMiAC(X,miAFrac=.49,maxMiss=.8):
	"""given a matrix of variant calls, X, 
	find the minor allele count for each sample 
	and return a vector of the results.
	Assumes X is nxm where:
	n = number of variant positions
	m = number of samples x 2
	Two adjacent columns are given to each sample,
	which is consistent with the binary matrices from 
	split-transcript and from the VCF files as these 
	typically describe diploids (2 copies of each potion).
	Calls are considered to be minor alleles iff their 
	call fraction in this population is < miAFrac.
	Samples with more than maxMiss missing (or nan) calls
	in this region will be marked as nan.
	"""
	n,m = X.shape
	fMiAC = np.zeros(m/2)

	# first get rid of smaples with excess missing calls
	# NOTE: not sure if we should do this before or after 
	#	wec calculate fractions...
	keep = _removeMissing(X,maxMiss)
	keep2 = _one2twoCol(keep)
	for i in range(n):
		minor = _isMinor(X[i,keep2],miAFrac)
		fMiAC[keep] += minor[::2] + minor[1::2]
	fMiAC[~keep] = np.nan
	return fMiAC

def mkMissingRegionFM(regionManifestPath,sampleIDPath,outFMPath,fBaseName='N:GNMC:annotation_missing:',maxMiss=.8):
	"""Make region based annotation FM that records 
	True if a sample has too many missing calls (missing>maxMiss),
	on the individual sample level, given 
	the transcript manifest and corresponding 
	data as extracted from split-transcripts.
	It is critical that the sampleIDPath points to a sample
	list that matches the columns of the region matrices 
	in regionManifestPath, no testing can be done!!
	"""
	#must have minNum or more variants to count as region
	minNum = 2

	# get the manifest file 
	regionReader = RegionManifestReader(regionManifestPath)
	samples = np.loadtxt(sampleIDPath,dtype=str,delimiter='\t')
	features = []
	
	count = 0
	maxCount = 100000
	for region in regionReader:
		#try:
		if region.nVar>=minNum:
			X = region.getData()
			features.append(fBaseName+region.name)
			missing = ~_removeMissing(X,maxMiss)
			if count==0:data = missing
			else: data = np.vstack([data,missing])

		count = count+1 
		if count>maxCount: break


	genUtil.saveFM(data,outFMPath,features,samples)


def mkRegionFM(regionManifestPath,sampleIDPath,outFMPath,fBaseName='N:GNMC:data',miAFrac=.49,maxMiss=.8):
	"""Make region based data FM,
	on the individual sample level, given 
	the transcript manifest and corresponding 
	data as extracted from split-transcripts.
	Currently calculate minor allele count (MiAC) features.
	"""
	#must have minNum or more variants to count as region
	minNum = 2

	# get the manifest file 
	regionReader = RegionManifestReader(regionManifestPath)
	samples = np.loadtxt(sampleIDPath,dtype=str,delimiter='\t')
	features = []
	
	count = 0
	maxCount = 100000
	for region in regionReader:
		#try:
		if region.nVar>=minNum:
			X = region.getData()
			features.append(fBaseName+':MiAC_'+region.name)
			fMiAC = getMiAC(X,miAFrac,maxMiss)
			if count==0:data = fMiAC
			else: data = np.vstack([data,fMiAC])

		count = count+1 
		if count>maxCount: break


	genUtil.saveFM(data,outFMPath,features,samples)

	

def main():
#	# test functions
#	wrkDir = '/Users/RTasseff/Projects/INOVA/DF5/basePipeLine'
#	manPath = wrkDir+'/P286_Data_Delivery_TOTAL_021214.txt'
#	fmPath = wrkDir+'/test_IND_20141121.fm'
#	samples = wrkDir+'/sampIDList_ITMI_VCF_DF5.tsv'
#	famFMPath = wrkDir+'/test_FAM_20141121.fm'
#	famSample = wrkDir+'/famID_DF5.tsv'
#	nbList = wrkDir+'/repNBList.tsv'
#
#	parseMkGenomeBatchFM(manPath,fmPath,samples)
#
#	genUtil.indFM2FamFM(fmPath,famFMPath,nbList,famSample)

	# --QC feautres --> 

#	# -get the QC FM of regions from comand line ->
#	# **done seperatly
#	# in this case each region is done individually,
#	# need to call the command in a loop from terminal
#	vcfPath = sys.argv[1]
#	fmOutPath = sys.argv[2]
#	getQCFM_VCF(vcfPath,fmOutPath)
#	# <--

#	# -addup QC regions ->
#	# done seperatly 
#	# loops though all regions to sum up a single FM
#	# for QC variables, assumes that all files are 
#	# in the same dir with <qcRegionName>_i.fm 
#	# where i = 1,2...,numQCRegions
#	qcRegionDir = sys.argv[1]
#	qcRegionName = sys.argv[2]
#	numQCRegions = int(sys.argv[3])
#	outQCFMPath = sys.argv[4]
#	summQCRegionFM(qcRegionDir,qcRegionName,numQCRegions,outQCFMPath)
#	# <-


#	
	# --Feature Matrix for transcript level data -->
	# running on 1 core on SGI for 61K regions requiered 13.5 hours
	regionManifestPath='/isb/rtasseff/data/transcripts_20141125/transcriptManifest_20141125.dat'
	sampleIDPath='/isb/rtasseff/data/transcripts_20141125/sampIDList_ITMI_VCF_DF5.tsv'
	if sys.argv[1]=='MiAC':
		outFMPath='/isb/rtasseff/data/featureMatrices/data_GNMC_Trans_IND_20150206.fm'
		mkRegionFM(regionManifestPath,sampleIDPath,outFMPath,fBaseName='N:GNMC:data',miAFrac=.49,maxMiss=.8)
	# now the FM for the annotations on misisng to be used in other methods:
	if sys.argv[1]=='miss':
		outFMPath='/isb/rtasseff/data/featureMatrices/annotation_GNMC_missing_20150206.fm'
		mkMissingRegionFM(regionManifestPath,sampleIDPath,outFMPath,maxMiss=.8)
	# <--

if __name__ == '__main__':
	main()
