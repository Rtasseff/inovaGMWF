import numpy as np

floatFmtStr = '%05.4E'


def getTrioOrder(fmPath,studyID='101'):
	"""Get the smaple order from the
	standard featurematrix.

	fmPath	str, path to standard featuremetrix
	"""
	fin = open( fmPath )
	line = fin.next() # first line is header with info
	samp = line.strip().split('\t')[1:] # first position is '.' place holder
	# assume family member id string is attached but not wanted
	for i in range( len(samp) ):
		tmp = samp[i].strip().split('-')
		if tmp[0] != studyID:
			raise ValueError ("ERR:0001 unexpected study id "+str(tmp[0]))
		elif tmp[2] != 'NB':
			raise ValueError ("ERR:0002 unexpected member id "+str(tmp[2]))
		samp[i] = tmp[1]
	return samp

def _parseSampID(sampID):
	# parse the typical smaple id form the CG manifest file
	# return the trioID and the member status
	# of NB = 0, M = 1, F = 2
	# assumes <member>(optional,-<member part 2>)-<study id>-<trio id>
	tmp = sampID.strip().split('-')
	membID = -1

	if tmp[0]=='NB':membID=0
	elif tmp[0]=='M':membID=1
	elif tmp[0]=='F':membID=2
	return(tmp[-1],membID)

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


def parseMkGenomeBatch(manPath,fmPath,trioList,genFeatName="BATCH:Genomic"):
	"""Parse the source files for the 
	genomic batch information, currently
	the manifest file only.
	Then make and save the batch feature matrix.
	manPath 	str, path to input manifest file
	fmPath 	str, path to output genomic batch feature matrix
	trioList	list of trio ids as strings
	genFeatName 	str, the general name of the these batch features
	"""
	# get data from the manifest
	# mainifest is small, use whats easiest:
	trioList = np.array(trioList,dtype=str)
	data = np.loadtxt(manPath,dtype=str,delimiter='\t')[1:]
	n,m = data.shape
	# create matrix, numerical features, order and version
	nFeat = 2 # version, ship data
	featHeader = np.array(['N:NB:'+genFeatName+':version','N:M:'+genFeatName+':version','N:F:'+genFeatName+':version','N:NB:'+genFeatName+':shipDate','N:M:'+genFeatName+':shipDate','N:F:'+genFeatName+':shipDate'],dtype=str)
	nMembers = 3
	nTrios = len(trioList)
	fm = np.array((np.zeros((nFeat*nMembers,nTrios)) + np.nan),dtype='|S15') # making it an str allows us to control format as we go and use non numerical data types if needed later.
	# loop through all data
	for i in range(n):
		# get sample 
		trioID,membID = _parseSampID(data[i,6])
		if membID >= 0:
			# version - 9
			iFeat = 0
			row = iFeat * nMembers + membID
			fm[row,trioID==trioList] = '%i'%(_parseVersion(data[i,9]))
			

			# ship date - 14
			iFeat = 1
			row = iFeat * nMembers + membID
			fm[row,trioID==trioList] = '%i'%(_paresShipDate(data[i,14]))

	# save the final matrix
	saveFM(fm,fmPath,featHeader,trioList,studyID='101')

	
def saveFM(data,fmPath,featHeader,trioHeader,studyID='101'):
	"""Given np arrays, save the feature matrix in 
	standard format.
	"""
	fout = open (fmPath,'w')
	n,m = data.shape
	if n!=len(featHeader):
		raise ValueError('ERR:0005, number of features and headers not equal')
	if m!=len(trioHeader):
		raise ValueError('ERR:0005, number of trios and headers not equal')
	fout.write('.\t')
	fout.write('\t'.join(trioHeader)+'\n')
	for i in range(n):
		fout.write(featHeader[i]+'\t'+'\t'+'\t'.join(np.array(data[i],dtype=str))+'\n')

	fout.close()



	
	





def main():
	fmInPath = '/Users/RTasseff/Projects/INOVA/DF5/basePipeLine/2014_09_29_new_clinical.fm'
	fmPath = '/Users/RTasseff/Projects/INOVA/DF5/basePipeLine/test.fm'
	manPath = '/Users/RTasseff/Projects/INOVA/DF5/basePipeLine/P286_Data_Delivery_TOTAL_021214.txt'
	trioList = getTrioOrder(fmInPath)
	parseMkGenomeBatch(manPath,fmPath,trioList)


if __name__ == '__main__':
	main()
