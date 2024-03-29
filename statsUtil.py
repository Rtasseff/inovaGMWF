#!/usr/bin/env python
#
# 
#     Copyright (C) 2003-2012 Institute for Systems Biology
#                             Seattle, Washington, USA.
# 
#     This library is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 2.1 of the License, or (at your option) any later version.
# 
#     This library is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public
#     License along with this library; if not, write to the Free Software
#     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
# 
# 20141114 RAT
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.stats import rankdata


nanValues = ['NA','NaN','na','nan']

def qqPlot(pValues,outpath):
	"""Make a qq plot of pValues
	saved to output path, path should have png extension.
	"""
	n = len(pValues)
	pExp = np.arange(1,0,-1.0/n)
	plt.plot(pExp,np.sort(pValues)[::-1])
	plt.plot([pExp[0],pExp[-1]],[pExp[0],pExp[-1]])
	plt.xlabel('expected -logP')
	plt.ylabel('observed -logP')
	plt.savefig(outpath,format='png')
	plt.clf()
	plt.close()



def fdr_bh(p_full,alpha=.05):
	"""Performs the Benjamini & Hochberg 1995
	multiple test correction for controlling
	the false discovery rate in familywise 
	analysis.  Tests must be independent or 
	positivly corrilated.
	p	original pvalues, np 1d array
	alpha 	threshold FDR, scalar float
	returns	h, regect or accept, 1d np bool array
	returus	p_adj, adjusted pvalues, 1d np array
	returns	pCrit, the critial p-value cut off
	"""
	n_full = len(p_full)
	isNAN = np.isnan(p_full)
	if np.any(isNAN):
		#print 'nan found, will be ignored'	
		p = p_full[~isNAN]
	else:
		p = p_full

	m = len(p)
	sortInd = np.argsort(p)
	pSort = p[sortInd]
	unsortInd = np.argsort(sortInd)
	pAdj = np.zeros(m)*np.nan
	gamma = (np.arange(m)+1)*(alpha/m)
	pTmp = m*pSort/(np.arange(m)+1)
	for i in range(m):
		pAdj[i] = np.min(pTmp[i:])

	pAdjUnsort = pAdj[unsortInd]
	rejSort = pSort<=gamma
	# find the largest value still under threshold
	# note they are sorted
	maxInd = np.sum(rejSort)-1
	if maxInd<0:
		pCrit = 0
	else:
		pCrit = pSort[maxInd]

	h = p<=pCrit
	
	h_full = np.zeros(n_full)
	pAdjUnsort_full = np.zeros(n_full)+np.nan
	h_full[~isNAN] = h
	pAdjUnsort_full[~isNAN] = pAdjUnsort

	return h_full,pAdjUnsort_full,pCrit


def fdr_bh_filterSortFile(inFilePath,outFilePath,alpha=.05,col=0,delimiter='\t',logTrans=False,nRows=''):
	"""Designed for large files that have been sorted from highest p-value 
	(least sig) to lowest (using eg linux sort command).
	Assuming 1 p-value per row in the same column of each row.
	Using the Benjamini & Hochberg 1995
	multiple test correction FDR
	We select the top n rows on inFilePath where FDR q < alpha and 
	print to a new file, outFilePath with FDR q-values appended to the last col.
	If rows have more entries than the p-values use col to indicate which delimiter separated column
	(starting at 0) holds the p-values.
	If pValues in inFilePath are -log transformed set logTrans=True (note numbers should be low to high in this case).
	If you know the number of rows you can speed the calculation by setting nRows.
	Cannot properly handle NA values, although if you set nRows to the number of non-NA rows
	it 'should' work.
	"""
	if nRows=='':
		fin = open(inFilePath)
		nRows = 0
		for line in fin:
			nRows+=1
	
	

	m = nRows

	fin = open(inFilePath)
	fout = open(outFilePath,'w')


	denom = float(m)
	q = 1
	for line in fin:
		p = line.strip().split(delimiter)[col]
		if p not in nanValues:
			p = float(p)
			if p<0: p=0
			if logTrans:p=10**(-1.0*p)
			qTmp = m*p/denom
			# every iteration denom drops
			denom += -1.0
			q = np.min([q,qTmp])

			if q<alpha:
				fout.write(delimiter.join(line.strip().split(delimiter))+delimiter+str(q)+'\n')

	fin.close()
	fout.close()


def getGroups(values,labels):
	"""assuming the labels correspond to the values 
	we return an np array of values for each label in
	a list of np arrays.
	Assumes values and labels are np arrays for easy indexing.
	"""
	unique = list(set(labels))
	groups = []
	cats = []
	for label in unique:
		if label not in nanValues:
			groups.append(values[labels==label])
			cats.append(label)

	return(groups,cats)

def makeConTable(x,y):
	"""Given two lists of categorical/binary observations
	(which can be strings, but nan will be treated as missing)
	we generate a 2-D contingency table.
	Here we assume that x and y are ordered the same such that 
	element x[i] and y[i] represent the same observation/sample.
	"""
	n = len(x)
	if n!=len(y):
		raise ValueError('the two data vectors must be the same length')

	xInt, xCat = cat2int(x)
	yInt, yCat = cat2int(y)
	xN = len(xCat)
	yN = len(yCat)

	conTable = np.zeros((xN,yN))

	for i in range(n):
		if xInt[i] >= 0 and yInt[i]>=0:
			conTable[xInt[i],yInt[i]]+=1

	return(conTable,xCat,yCat)

def cat2int(y):
	"""change a set of str category labels with n
	unique values into an int array with unique 
	values of numbers from 0 to n-1
	returns the new int array and list of cat labels.
		a list of categories corresponding to int value
	nan values (either 'nan' or the np object) will be preserved;
	however int arrays cannot have nan so preserved as -1.
	input array y must be an np array for indexing reasons.

	as of now integers cannot have nan values, setting this to -1
	"""
	unique = list(set(y))
	yNew = np.array(np.zeros(len(y)),dtype=int)
	count = 0
	cats = []
	for i in range(len(unique)):
		tmp = unique[i]
		missing=False
		if type(tmp)==np.string_:
			if tmp in nanValues:missing=True
			yNew[y==unique[i]] = -1
		else: 
			if np.isnan(tmp): missing=True
			yNew[np.isnan(y)] = -1

		if not missing:
			yNew[y==unique[i]] = count
			count += 1
			cats.append(unique[i])
			

	return(yNew,cats)

	
def plotPairwise(x,y,varType=['',''],varName=['',''],outfile=''):
	"""Diffrent combinations of variable types
	are visulized by diffrent plots, this funciton
	uses a basic visulization for the correct combinaiton.
	x	np str array, if type = N then this should be 
		convertable to float array
	y	same as x but for other variable
	varType	list with 2 values to indicate type of x and y 
		if left blank we assume x and y are split lines form a 
		feature matrix and the first entry is the label
	varName list with 2 values to indicate name of x and y
		if this *and* varType above is blank this will be taken as the 
		label (first entry) of x and y.

	"""
	if varType[0]=='':
		varType[0] = x[0].split(':')[0]
		if varName[0]=='':varName[0] = x[0]
		x = x[1:]
		varType[1] = y[0].split(':')[0]
		if varName[1]=='':varName[1] = y[0]
		y = y[1:]

	if not varType[0]=='N' and not varType[0]=='C' and not varType[0]=='B':
		raise ValueError( 'Variable type for x is unknwon: '+varType[0])
	if not varType[1]=='N' and not varType[1]=='C' and not varType[1]=='B':
		raise ValueError( 'Variable type for y is unknwon: '+varType[1])


	if varType[0] == 'B': varType[0] = 'C' # no diff here
	if varType[1] == 'B': varType[1] = 'C' # no diff here

	
	# check if both numerical:
	if varType[0]=='N' and varType[1]=='N':
		xFloat = _getFloat(x)
		yFloat = _getFloat(y)

		# must ignore all nan values
		keep = np.array(np.ones(len(xFloat)),dtype=bool)
		keep[np.isnan(xFloat)]=False
		keep[np.isnan(yFloat)]=False
		# non parametric tests used so lets consider ranks
		xRank = rankdata(xFloat[keep])
		yRank = rankdata(yFloat[keep])
		#xRank = np.argsort(np.argsort(xFloat[keep]))
		#yRank = np.argsort(np.argsort(yFloat[keep]))

		plt.plot(xRank,yRank,'o',label='Data')
		try:
			# scatter plot with line from least squares
			A = np.vstack([xRank, np.ones(len(xRank))]).T
			m, c = np.linalg.lstsq(A, yRank)[0]
			plt.plot(xRank, m*xRank + c, 'r', label='Fitted line')
			plt.legend()

		except ValueError as ve:
			warnings.warn("Could not fit line for "+varName[0]+" vs "+varName[1]+": "+str(ve),UserWarning)
		
		
		
		plt.xlabel(varName[0]+' ranks')
		plt.ylabel(varName[1]+' ranks')
		#plt.ylabel(varName[1].split(':')[-1])
	elif varType[0]=='C' and varType[1]=='C':
		# right now we are just going with a simple 
		# stacked bar plots, more advanced things
		# like mosics are avalible in R or pythons statsmodels

		conTable,xCat,yCat = makeConTable(x,y)
		varInd = 0
		# we want more bars less stacks so n shoudl be smallest
		if len(xCat)>len(yCat):
			tmpCat = yCat
			yCat = xCat
			xCat = tmpCat
			conTable = conTable.T
			varInd = 1
		n,m = conTable.shape
		# consider normalized contigency table
		conTable = conTable/np.sum(conTable,0)
		ind = np.arange(m)
		bottomTmp = np.zeros(m)
		colorList = ['b','g','r','c','y','m']
		nColor = len(colorList)
		colorInd = 0
		for i in range(n):
			plt.bar(ind,conTable[i],bottom=bottomTmp,color=colorList[colorInd],label=varName[varInd]+'::'+str(xCat[i]))
			bottomTmp += conTable[i]
			colorInd +=1
			if colorInd==nColor:colorInd=0

		plt.ylabel('Normalized Count')
		plt.xlabel(varName[abs(varInd-1)])
		plt.xticks(ind,np.array(yCat,dtype=str),rotation=10)
		plt.legend()
			


	else:
		#####
		#fig = plt.figure(figsize=(8,3))

		# number vs category 
		#set up variables
		if varType[0]=='N':
			varN = _getFloat(x)
			varC = y
			varNInd = 0
			varCInd = 1
		else:
			varN = _getFloat(y)
			varC = x
			varNInd = 1
			varCInd = 0

		groups,unique = getGroups(varN,varC)
		bp = plt.boxplot(groups)
		plt.xticks(range(1,(len(unique)+1)),unique,rotation=10)
		plt.xlabel(varName[varCInd])
		plt.ylabel(varName[varNInd])
		#plt.ylabel(varName[varNInd].split(':')[-1])
	
	if outfile=='':
		plt.show()
	else:
		plt.savefig(outfile)
		plt.clf()
		plt.close()



def _getFloat(x):
	"""Create a float of this string"""
	if x.dtype==np.dtype('float64'):xFloat = x
	else:
		# had issues with S2 dtypes cutting off the nan to na
		xFloat = np.array(x,dtype='|S15')
		for i in range(len(xFloat)):
			if xFloat[i] in nanValues:
				xFloat[i] = 'nan'
		xFloat = np.array(xFloat,dtype=float)
	return xFloat
		
