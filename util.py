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
		raise ValueError( 'Variable type for x is unknwon: '+varType[1])


	if varType[0] == 'B': varType[0] = 'C' # no diff here
	if varType[1] == 'B': varType[1] = 'C' # no diff here

	# check if both numerical:
	if varType[0]=='N' and varType[1]=='N':
		plt.plot(np.array(x,dtype=float),np.array(y,dtype=float))
		plt.xlabel(varName[0])
		plt.ylabel(varName[1])
	elif varType[0]=='N' or varType[1]=='N':
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
		plt.xticks(unique)
		plt.xlabel(varName[varCInd])
		plt.ylabel(varName[varNInd])
	
	if outfile=='':
		plt.show()
#	else:
#		need to put logic to save as file here





	
