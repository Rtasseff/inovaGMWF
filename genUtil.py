"""
General python2.7 utility methods for use on the inova project.

RAT 20141121


"""

import warnings

import numpy as np
import statsUtil 

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


def cgID2isb(cgID):
	"""CG (as in manifest file) 
	puts member ids as prefixes, but ISB puts
	them as suffixes.  This will change that order.
	Note that 'na' will be returned if member id is 
	not yet accounted for.  For example, sample
	PO-101-657-22 has PO member id that is 
	not accounted for in any logic so far or
	analysis done. Note that the 3 sets 
	of numbers complicates matters with this ID.
	"""
	isbID = 'na'
	tmp = cgID.split('-')
	# hard coded allowed prefixes 
	if tmp[0] in ['F','M','NB']:
		if len(tmp)==3:
			isbID = tmp[1]+'-'+tmp[2]+'-'+tmp[0]
		if len(tmp)==4:
			isbID = tmp[2]+'-'+tmp[3]+'-'+tmp[0]+'-'+tmp[1]
	return(isbID)

def itmiID2isb(itmiID):
	"""ITMI (as in VCF) puts member ids as prefixes and no dash 
	for twin differences, but ISB puts
	them as suffixes with a dash.  This will change that order.
	Note that 'na' will be returned if member id is 
	not yet accounted for.  For example, sample
	PO-101-657-22 has PO member id that is 
	not accounted for in any logic so far or
	analysis done. Note that the 3 sets 
	of numbers complicates matters with this ID.
	"""
	isbID = 'na'
	tmp = itmiID.split('-')
	# hard coded allowed prefixes 
	if tmp[0] in ['F','M','NB']:
		if len(tmp)==3:
			isbID = tmp[1]+'-'+tmp[2]+'-'+tmp[0]

	if len(tmp[0])>2 and tmp[0][:2] =='NB':
		if len(tmp)==3:
			isbID = tmp[1]+'-'+tmp[2]+'-NB-'+tmp[0][2:]

	return(isbID)

def cgID2itmi(cgID):
	"""Transform the CG id to an itmi id.
	"""
	itmiID = cgID
	tmp = cgID.split('-')

	# only diff is to join the NB-<twin id> -> NB<twin id>
	if len(tmp)==4 and tmp[0]=='NB':
		itmiID = tmp[0]+tmp[1]+'-'+tmp[2]+'-'+tmp[3]
	return (itmiID)

def isbID2itmi(isbID):
	"""Convert from the ISB sample ID for an individual
	to the itmi ID (used in VCF eg).  The differences
	is to switch the member ID to a prefix and to 
	remove the separator in twin IDs.
	"""
	itmiID = 'na'
	tmp = isbID.split('-')
	if len(tmp)==3:
		itmiID = tmp[2]+'-'+tmp[0]+'-'+tmp[1]
	elif len(tmp)==4:
		itmiID = tmp[2]+tmp[3]+'-'+tmp[0]+'-'+tmp[1]

	return(itmiID)

def isbID2cg(isbID):
	"""Convert from the ISB sample ID for an individual
	to the CG ID (used in manifest eg).  The differences
	is to switch the member ID to a prefix.
	"""
	cgID = 'na'
	tmp = isbID.split('-')
	if len(tmp)==3:
		cgID = tmp[2]+'-'+tmp[0]+'-'+tmp[1]
	elif len(tmp)==4:
		cgID = tmp[2]+'-'+tmp[3]+'-'+tmp[0]+'-'+tmp[1]

	return(cgID)


def idType(sampID):
	""" Find out the type of sample id (cg, itmi, isb)
	sampID is, and return str indicator, 
	if there are multiple possibilities return all str 
	indicators in a list.
	'na' is returned if the format is not recognized 
	"""
	tmp = sampID.split('-')
	if tmp[0]=='PO' and len(tmp)==4:sampIDType = ['itmi','cg']
	elif tmp[0] in ['NB','M','F'] and len(tmp)==3:sampIDType = ['itmi','cg']
	elif tmp[0]=='NB' and len(tmp)==4 and tmp[1] in ['A','B','C','D']:sampIDType ='cg'
	elif len(tmp[0])>2 and tmp[0][:2] =='NB' and len(tmp)==3:sampIDType ='itmi'
	elif tmp[-1] in ['NB','M','F'] and len(tmp)==3:sampIDType = 'isb'
	elif tmp[-1] in ['A','B','C','D'] and len(tmp)==4:sampIDType = 'isb'
	else: sampIDType = 'na'

	return(sampIDType)

def idListType(smapIDList):
	""" Find out the types of sample id (cg, itmi, isb)
	sampIDList is, and return str indicator, 
	if there are multiple possibilities return all str 
	indicators in a list.
	Note: we expect all ids are consistent and stop looking 
	once one id is determinant.
	"""
	sampIDType = 'na'
	for sampID in smapIDList:
		sampIDType = idType(sampID)
		if type(sampIDType)==str:break

	return sampIDType

def convIDType(inID,outIDType,inIDType=''):
	"""Convert the sample id, inID, of type inIDType 
	to the outIDType, return new str id.
	If inIDType blank, detect automatically
	allowable types: 'isb', 'itmi', 'cg'.
	"""
	if inIDType=='':
		_inIDType = idType(inID)
		if type(_inIDType)==list:_inIDType=_inIDType[0]
	else: _inIDType=inIDType

	if _inIDType=='isb':
		if outIDType=='cg':
			outID = isbID2cg(inID)
		elif outIDType=='itmi':
			outID = isbID2itmi(inID)
		else:
			raise ValueError('ERR0021: output id type, '+ outIDType+', not recognized.')	
	elif _inIDType=='itmi':
		if outIDType=='isb':
			outID = itmiID2isb(inID)
		elif outIDType=='cg':
			raise ValueError('ERR0023: conversion from cg to itmi not implemented.')

		else:
			raise ValueError('ERR0022: output id type, '+ outIDType+', not recognized.')
	elif _inIDType=='cg':
		if outIDType=='itmi':
			outID = cgID2itmi(inID)
		elif outIDType=='isb':
			outID = cgID2isb(inID)
		else:
			raise ValueError('ERR0024: output id type, '+ outIDType+', not recognized.')
	else:
		raise ValueError('ERR0020: input id type, '+ _inIDType+', not recognized for sample: '+inID)
	return(outID)

def nb2indList(nbList,outIDType='isb'):
	"""Create a list of individual sample ids
	that corresponds with the new born list (str list nbList).
	The format, or id type for the sample ids is 
	specified by outIDType.
	Return str list of sample ids
	"""
	# get the type of ids
	inIDType = idListType(nbList)
	# if degenerative can choose any 
	if type(inIDType)==list:inIDType=inIDType[0]
	# if na we have issue 
	if inIDType == 'na':
		raise ValueError('ERR0025: input id type derived from list is not recognized.')

	# lets always start with ISB format:
	if inIDType=='isb':
		_nbList = nbList
	else:
		_nbList = []
		for value in nbList:
			_nbList.append(convIDType(value,'isb',inIDType=inIDType))
	# create individual list F M and NB
	indList = []
	for value in _nbList:
		tmp = value.split('-')
		# new born
		indList.append(value)
		# mother
		indList.append(tmp[0]+'-'+tmp[1]+'-M')
		# father
		indList.append(tmp[0]+'-'+tmp[1]+'-F')
	# convert if needed
	if outIDType!='isb':
		for i in range(len(indList)):
			indList[i] = convIDType(indList[i],outIDType,inIDType='isb')

	return(indList)

def mkNBDict(nbList):
	"""Make a dictionary of family ids linked to 
	representative newborns.
	"""
	nbDic={}
	for nb in nbList:
		tmp = nb.split('-')
		if len(tmp)==3:nbDic[tmp[1]]=tmp[2]
		if len(tmp)==4:nbDic[tmp[1]]=tmp[2]+'-'+tmp[3]
	return(nbDic)


def getRepNBList(nbListPath):
	"""Create a list (np str array) that contains
	all usable newborns for family centric analysis.
	Used in mapping individuals to -FAM
	for family centric analysis.
	Tab sep input list should have all and only
	newborns used to represent a family
	unit. An error will be returned if
	a family id xxx in 101-xxx-NB is 
	included twice.

	The list will be used for simplicity;
	effectively -F -M and -NB are 
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


def mapInd2FAM(indID,nbList):
	"""effectively -F and -M are 
	just replaced with -FAM. However,
	for multi births only the one
	in the list will be used in the 
	-FAM family centric analysis.
	The nbList will be used to determine
	if the passed individual is nb and if so
	is the representative nb from the list.
	Returns the new member id and new -FAM id.
	Non representative nb's or unrecognized 
	member ids will return na.
	indID	str, isb formated id for individual sample
	nbList	np str array of usable representative newborns
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

def mapFam2Ind(famID,nbDic,member):
	"""effectively -FAM is 
	just replaced with -member. However,
	for multi births (member='NB') only the one
	in the nbDic will be returned.
	The nbDic will be used to determine
	if the requested individual is nb and if so
	return the representative nb from the list.
	"""

	tmp = famID.split('-')
	if member in ['F','M']:
		indID = tmp[0]+'-'+tmp[1]+'-'+member
	# only checking for multi births which have 4 '-'
	elif member=='NB':
		indID = tmp[0]+'-'+tmp[1]+'-'+nbDic[tmp[1]]
	else: indID = 'na'
	return(indID)

	

def indFM2FamFM(indFMPath,famFMPath,nbList,famSample):
	"""Convert a given individual FM to a family based 
	FM.  New Family based FM will use the sample header 
	in 'famSample', and will chose representative newborns
	(ie from twins) using the complete list of family representing 
	newborns in nbList.  This will automatically identify
	the individual sample id formats in the header of indFMPath 
	(id formats: isb, itmi, cg) and will convert from these
	to ISB format (which should be consistent with famSample).
	indFMPath 	str, path to input individual feature matrix
	famFMPath 	str, path to output family based feature matrix
	nbList		if list, used as actual nbList
			if a single str then we assume its a path to tab 
			sep file containing the nbList
	famSample	if list, used as actual family samples with ISB ids as strings 
			if a single str then we assume its a path to tab 
			sep file containing the sample ids 
	"""
	if type(famSample)==str:famSample = np.loadtxt(famSample,dtype=str,delimiter='\t')
	else:famSample = np.array(famSample,dtype=str)

	if type(nbList)==str:nbList = getRepNBList(nbList)
	else:nbList = np.array(nbList,dtype=str)
	nbDic =  mkNBDict(nbList)

	# members we expect to have
	members=['NB','M','F']
	nMembers = len(members)
	famSufix = 'FAM'


	# open the individual FM
	indFM = open(indFMPath)
	# get the samples
	line = indFM.next()
	indSampleList = np.array(line.strip().split('\t')[1:],dtype=str)
	# make an index of the individual samples for latter reference:
	indSampleIndex = np.arange(len(indSampleList),dtype=int)

	# find out format of individual sample IDs:
	indSampIDType = idListType(indSampleList)
	if type(indSampIDType)==list:indSampIDType=indSampIDType[0]
	if indSampIDType=='na':
		raise ValueError ('ERR:0002 the individual sample IDs do not match any known format, ids in FM at '+indFMPath)



	# open family FM
	famFM = open(famFMPath,'w')
	famFM.write('.'+'\t'+'\t'.join(famSample)+'\n')
	

	# create temporary data matrix to hold set of family features
	famData = np.array((np.zeros(len(famSample)) + np.nan),dtype='|S15')
	# preallocate space for the ind data, be sure to include 1 extra spot at end for nan's to be referenced at
	indData = np.array((np.zeros(len(indSampleList)+1) + np.nan),dtype='|S15')

	# create a map, index of ind sample for each spot in FAM 
	# allocate
	fam2indSampleIndex = np.zeros((nMembers,len(famSample)),dtype=int)
	# search
	# loop through members
	for i in range(nMembers):
		member = members[i]
		# loop though all family ids
		for j in range(len(famSample)):
			# get the sample id
			indSampID = mapFam2Ind(famSample[j],nbDic,member)
			# convert this isb ID to correct format for comparison to individual ids
			if indSampIDType == 'itmi':indSampID = isbID2itmi(indSampID)
			elif indSampIDType == 'cg':indSampID = isbID2cg(indSampID)
			elif not indSampIDType == 'isb':raise RuntimeError ('RTE:0003 for some reason you are getting back an unknown sample id type indicator')

			# find the corresponding index
			if np.any(indSampID==indSampleList):
				ind = indSampleIndex[indSampID==indSampleList][0]
				fam2indSampleIndex[i,j] = ind
			else:
				# latter we will use this index as an 'nan'
				fam2indSampleIndex[i,j] = len(famSample)
				# warn if sample not found
				warnings.warn("The individual sample "+indSampID+", corresponding to desired family sample, "+famSample[j]+", was not found.",UserWarning)  





	# loop through individual matrix
	for line in indFM:
		tmp=line.strip().split('\t')
		fname = tmp[0].split(':')
		indData[:-1] = tmp[1:]
		for i in range(nMembers):
			member = members[i]
			famData = indData[fam2indSampleIndex[i]]
			# now we have the data time to write:
			# feature name
			famFM.write(fname[0]+':'+member+':'+':'.join(fname[1:]))
			# now data
			famFM.write('\t'+'\t'.join(famData)+'\n')

	famFM.close()
	indFM.close()


			
		

def plotter(pairList,fmPath,outDir='.',prefix='pwPlot'):
	""" Given the input feature matrix at fmPath,
	create a pairwise plot for each pair contained in
	pairList.  If pairList is a str, assume its a path
	The file should have a tab sep pair of
	feature names, one on each line. 
	otherwise assume its a list of 2 value lists 
	of the feautre names.  Plots saved as:
	<outDir>/<prefix>_<name 1>_vs_<name 2>.png.
	"""
	if type(pairList)==str:
		pairNameList = np.loadtxt(pairList,dtype=str,delimiter='\t')
	else:
		pairNameList = pairList

	n = len(pairNameList)
	fm = open(fmPath)
	# skip header 
	line = fm.next()
	data = line.strip().split('\t')
	# assuming pair list is short so we will allocate as we go
	# create mini matrix for more efficent pair finding
	fmMini = np.array(data,dtype=str)
	# create indexing matrix for fast recall
	pairInd = np.zeros((n,2),dtype=int) - 1
	fmMiniCount = 1 # skip header row
	for line in fm:
		data = line.strip().split('\t')
		fname = data[0]
		keep = False
		# check all pairs for same
		for i in range(n):
			# pairs can have same name so go through all of them
			if fname==pairNameList[i][0]:
				keep = True
				pairInd[i,0] = fmMiniCount

			elif fname==pairNameList[i][1]:
				keep = True
				pairInd[i,1] = fmMiniCount
		# added if needed, but dont add if already added
		if keep:
			fmMini = np.vstack([fmMini,data])
			fmMiniCount += 1
			
	fm.close()
	# now we have the mini fm lets do the plots
	for i in range(n):
		x = fmMini[pairInd[i,0]]
		y = fmMini[pairInd[i,1]]
		# feature names can make bad file names
		# replace : with -
		# and . or / with '?'
		outfile = outDir+'/'+prefix+'_'+x[0].replace('/','?').replace('.','?').replace(':','-')+'_vs_'+y[0].replace('/','?').replace('.','?').replace(':','-')+'.png'
		try:
			statsUtil.plotPairwise(x,y,outfile=outfile,varType=['',''],varName=['',''])
		except ValueError as ve:
			warnings.warn("Could not create the top scoring pairwise figure {}.\n\t Value error occurred:\n\t{}".format(outfile,ve))
			print outfile
			print ve
			
	#remove mini fm
	del fmMini
	del pairInd


def main():
	pairPath = 'batchPairs.tsv'
	fmPath = '/titan/ITMI1/projects/gamcop/data/featureMatrices/comb_CP_META_filtered_20141125.fm'
	plotter(pairPath,fmPath,outDir='.',prefix='pwPlot')

if __name__ == '__main__':
	main()
