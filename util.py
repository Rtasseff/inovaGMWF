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


def getGroups(values,labels):
	"""assuming the labels correspond to the values 
	we return an np array of values for each label in
	a list of np arrays.
	Assumes values and labels are np arrays for easy indexing.
	"""
	unique = list(set(labels))
	groups = []
	for label in unique:
		groups.append(values[labels==label])

	return(groups,unique)

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
		if not np.isnan(xInt[i]) and not np.isnan(yInt[i]):
			conTable[xInt[i],yInt[i]]+=1

	return(conTable,xCat,yCat)


	
def plotPairwise(x,y,varType=['',''],varName=['',''],outfile=''):
	"""Diffrent combinations of variable types
	are visulized by diffrent plots, this funciton
	uses a basic visulization for the correct combinaiton.
	x	np str array, if type = N then this should be 
		convertable to float array
	y	same as x but for other variable
	varType	tuple with 2 values to indicate type of x and y 
		if left blank we assume x and y are lines form a 
		feature matrix and the first entry is the label
	"""
	if varType[0]=='':
		varType[0] = x[0].split(':')[0]
		varName[0] = x[0].split(':')[-1]
		x = x[1:]
		varType[1] = y[0].split(':')[0]
		varName[1] = y[0].split(':')[-1]
		y = y[1:]

	if not varType[0]=='N' and not varType[0]=='C' and not varType[0]=='B':
		raise ValueError( 'Variable type for x is unknwon: '+varType[0])
	if not varType[1]=='N' and not varType[1]=='C' and not varType[1]=='B':
		raise ValueError( 'Variable type for y is unknwon: '+varType[1])


	if varType[0] == 'B': varType[0] = 'C' # no diff here
	if varType[1] == 'B': varType[1] = 'C' # no diff here

	# check if both numerical:
	if varType[0]=='N' and varType[1]=='N':
		# scatter plot with line from least squares
		A = np.vstack([x, np.ones(len(x))]).T
		m, c = np.linalg.lstsq(A, y)[0]

		plt.plot(np.array(x,dtype=float),np.array(y,dtype=float),'o',label='Data')
		plt.plot(x, m*x + c, 'r', label='Fitted line')
		plt.legend()

		plt.xlabel(varName[0])
		plt.ylabel(varName[1])
	elif varType[0]=='C' and varType[1]=='C':
		# right now we are just going with a simple 
		# stacked bar plots, more advanced things
		# like mosics are avalible in R or pythons statsmodels

		conTable,xCat,yCat = makeConTable(x,y)
		n,m = conTable.shape
		ind = np.arange(m)
		bottomTmp = np.zeros(m)
		colorList = ['b','g','r','c','y','m']
		nColor = len(colorList)
		colorInd = 0
		for i in range(n):
			plt.bar(ind,conTable[i],bottom=bottomTmp,color=colorList[colorInd],label=varName[0]+'::'+str(xCat[i]))
			bottomTmp = conTable[i]
			colorInd +=1
			if colorInd==nColor:colorInd=0

		plt.ylabel('Count')
		plt.xlabel(varName[1])
		plt.xticks(ind,np.array(yCat,dtype=str),rotation=45)
		plt.legend()
			


	else:
		#set up variables
		if varType[0]=='N':
			varN = np.array(x,dtype=float)
			varC = y
			varNInd = 0
			varCInd = 1
		else:
			varN = np.array(y,dtype=float)
			varC = x
			varNInd = 1
			varCInd = 0

		groups,unique = getGroups(varN,varC)
		bp = plt.boxplot(groups)
		plt.xticks(range(1,(len(unique)+1)),unique,rotation=45)
		plt.xlabel(varName[varCInd])
		plt.ylabel(varName[varNInd])
	
	if outfile=='':
		plt.show()
	else:
		plt.savefig(outfile)
		plt.clf()
		plt.close()




