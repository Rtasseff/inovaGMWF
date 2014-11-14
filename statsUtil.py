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
# 20120904 RAT
import numpy as np
import scipy.stats as stats
#import rpy2.robjects as robjects
import scipy.linalg
import operator as op
import itertools
import matplotlib.pyplot as plt

def iqr(x):
	"""returns the inter quantile distance for
	the variables in x.  
	if x is a matrix we assume the columns are
	different distributions and the iqr returned
	will be a vector corrisponding to the columns of x.
	"""
	n = len(x)
	if len(x.shape)==1:
		xSort = np.sort(x)
		q1 = xSort[int(.25*n)]
		q3 = xSort[int(.75*n)]
	elif len(x.shape)==2:
		n,m = x.shape
		q1 = np.zeros(m)
		q3 = np.zeros(m)
		for i in range(m):
			xSort = np.sort(x[:,i])
			q1[i] = xSort[int(.25*n)]
			q3[i] = xSort[int(.75*n)]

	return(q3-q1)

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

def fisherComb(p):
	"""Apply fisher method to combine 
	the pvaules in the np array p, to
	get a single pvalue and return"""
	n = len(p)
	c = -2*np.sum(np.log(p))
	pComb = 1-stats.chi.cdf(np.sqrt(c),2*n)
	return pComb


def mi(bins,calcP=False,v=False):
	"""calculates mutual information 
	\sum_i\sum_j p(i,j)*log_2[p(i,j)/(p(i)*p(j))]
	for 2 variable histogram defined by bins
	bins	int 2d np array
	returns MI
	"""
	tot = float(np.sum(bins))
	n,m = bins.shape
	p_x = np.sum(bins,1)/tot
	p_y = np.sum(bins,0)/tot
	null_chi = 0
	mi = 0
	for i in range(n):
		for j in range(m):
			p_xy = bins[i,j]/tot
			frac = p_xy/(p_x[i]*p_y[j])
			if frac>0:
				mi = mi + p_xy*np.log2(frac)
			elif v:
				print 'warning some bins are empty, skipping them'
			if calcP:
				E_xy = tot*(p_x[i]*p_y[j])
				if E_xy == 0 and v: print 'WARNING 0 expectation in pvalue calc'
				if E_xy<1 and v: print 'WARNING expectation < 1 in pvalue calc, not good'
				if E_xy<5 and v: print 'WARNING expectation < 5 in pvalue calc, Cochran criterion may be violated'
				null_chi = null_chi + (bins[i,j]-E_xy)**2/E_xy

	if calcP:
		v = (n)*(m)
		p = 1-stats.chi.cdf(np.sqrt(null_chi),v)
		return mi,p
	else: return mi


def sampleIndexWR(n):
	"""Given an array size, n, returns
	the indicies for random uniform sampling 
	of the array, with replacment.
	"""
	indFloat = sampleWR(range(n),n)
	ind = map(int,indFloat)
	return(ind)	

def sampleWR(pop,k=0):
	"""given a population, pop,
	return a new sample, size=k,
	by sampling pop with replacment
	"""
	if k<1:k=len(pop)
	n = len(pop)
	sel = np.zeros(k)
	for i in range(k):
		index = np.random.randint(0,n)
		sel[i] = pop[index]

	return sel

def cmds(D):
	"""Preform clasical multidimensional scaling.
	using D as the pair-wsie distances Matrix
	returns Y the column matrix defining the
	new corrdinates that retain maximum conserved 
	distance and E the eigan values that describe
	the contribution of each dimension (Y col).
	"""
	[n,m] = D.shape
	eps = 1E-21
	if n!=m:
		raise ValueError('Wrong size matrix')

	# Construct an n x n centering matrix
	# The form is P = I - (1/n) U where U is a matrix of all ones
	P = np.eye(n) - (1/float(n) * np.ones((n,n)))
	B = np.dot(np.dot(P,-.5*D**2),P)
	# Calculate the eigenvalues/vectors
	[E, V] = scipy.linalg.eig((B+B.T)/2) # may help with round off error??
	E = np.real(E) # these come out as complex but imaginary part should be ~eps
	V = np.real(V) # same argument here
	# sort that mo fo
	ind = np.argsort(E)[::-1]
	E = E[ind]
	V = V[:,ind]
	# lets now create our return matrix 
	if np.sum(E>eps)==0:
		Y = 0
	else:
		Y = V[:,E>eps]*np.sqrt(E[E>eps])
	return(Y,E)


def makeDesingMat(y):
	"""transfomr a class label vector into the design / indicator matrix Y
	y	class vector corresponding to observations
	return	matrix with observations on rows and
		classes on cols, with 1 indicating obs in class 
	"""
	if np.sum(np.unique(y) != np.array(range(len(np.unique(y)))))>0:
		raise ValueError('Class vector y must be a numeric vector, with values as 0,1,2,...')
	return(np.eye(len(np.unique(y)))[y,:])


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


def cat2int(y):
	"""change a set of str category labels with n
	unique values into an int array with unique 
	values of numbers from 0 to n-1
	returns the new int array and list of cat labels.
		a list of categories corresponding to int value
	nan values (either 'nan' or the np object) will be preserved.
	input array y must be an np array for indexing reasons.
	"""
	unique = list(set(y))
	yNew = np.array(np.zeros(len(y)),dtype=int)
	count = 0
	for i in range(len(unique)):
		tmp = unique[i]
		missing=False
		if type(tmp)==np.string_:
			if tmp.lower()=='nan':missing=True
		else: 
			if np.isnan(tmp): missing=True

		if missing==False:
			yNew[y==unique[i]] = count
			count += 1
		else:
			yNew[y==unique[i]] = np.nan

	return(yNew,unique)

def softThresh(X,thresh):
	"""Take X, and do soft 
	threshold on it and return that.
	"""
	tmp = np.abs(X) - thresh
	tmp[tmp<0] = 0
	return(tmp*np.sign(X))

def enrich(n_pos,n_draw,total_pos,total_all):
	"""Standard enrichment test usign hypergeometric 
	distribution.
	ie I got n_pos red balls out of n_draw draws from
	a bag with total_pos red balls and total all balls of 
	any color, calculate the probability of drawing n_pos
	red balls at random.
	"""
	
	p = stats.hypergeom.sf(n_pos-1,total_all,total_pos,n_draw)
	return(p)

def enrichList(posList,drawList,backList):
	"""Standard enrichment test usign hypergeometric 
	distribution.
	"""
	total_all = len(set(backList))
	total_pos = len(set(posList))
	n_draw = len(set(drawList))
	n_pos = len(set(posList).intersection(set(drawList)))
	p = enrich(n_pos,n_draw,total_pos,total_all)
	
	return(p)


def rankSum(x,y,forceExact=False):
	"""This is a two-sided Wilcoxon rank sum test,
	robust and nonparametric way to estimate if 
	two samples x and y have diffrent medians.
	Test assumes that two samples x,y are independent but 
	from distributions with the same median.  
	This is tested agains the alternate hyp that 
	the two samples are indipendent but come from 
	distributions with diffrent medians.
	Identical to Mann-Whitney U test.
	By defult the normal approximation is
	made for the test statistic IF the 
	sum of the samples has 15 or more obs;
	ortherwise, we enumerate the combinations
	to find the exact p-value.
	If forceExact==True then the exact p-value
	is calculated regardless of sample size.
	returns 
	p-value - sacalr
	z, approximate test statistic - scalar

	This function is similar to matlab ranksum.
	"""
	# find smallest 
	xLen = len(x)
	yLen = len(y)
	if xLen<=yLen:
		n = xLen
		sampSm = x
		sampLrg = y
	else:
		n = xLen
		sampSm = x
		sampLrg = y
	
	
	sampComb = np.append(sampSm,sampLrg)
	# check for no varriation
	if np.var(sampComb) == 0: return(np.nan,np.nan)
	# get ranks, possibly tied
	ranks,tieAdj = tiedRank(sampComb)
	ranks = ranks +1
	# get the stat
	xRank = ranks[:n]
	w = np.sum(xRank)
	# calculate the z statistic
	cor = 2*tieAdj/( (xLen+yLen)*(xLen+yLen-1) )
	wVar = xLen*yLen*((xLen+yLen+1) - cor)/12.
	wMean = n*(xLen+yLen+1)/2.
	wC = w-wMean
	z = (wC - .5 * np.sign(wC))/np.sqrt(wVar)

	# more complex methods exist, but I am using this 
	# simple way to deal with small samples
	if forceExact or (xLen+yLen)<15:
		allComb = chooseAllComb(ranks,n)
		allCombSum = np.sum(allComb,1)
		allCombLen = len(allCombSum)
		plo = np.sum(allCombSum<=w)/float(allCombLen)
		phi = np.sum(allCombSum>=w)/float(allCombLen)
		p = np.min([plo,phi]) 
		p = np.min([2*p,1])
	else:
		p = 2 * stats.norm.cdf(-np.abs(z))

	
		
	return(p,z)

def rankSum1S(x,y,forceExact=False):
	"""This is a one-sided Wilcoxon rank sum test,
	robust and nonparametric way to estimate if 
	the median of sample x > y.
	Test assumes that x,y are independent.  
	This is tested aginst the alternate hyp that 
	the two samples are indipendent but come from 
	distributions with the same medians.
	Identical to a one-sided Mann-Whitney U test.
	By defult the normal approximation is
	made for the test statistic IF the 
	sum of the samples has 15 or more obs;
	ortherwise, we enumerate the combinations
	to find the exact p-value.
	If forceExact==True then the exact p-value
	is calculated regardless of sample size.
	returns 
	p-value - sacalr
	z, approximate test statistic - scalar

	This function is similar to matlab ranksum.
	"""
	# find smallest 
	xLen = len(x)
	yLen = len(y)
	n = xLen
	
	sampComb = np.append(x,y)
	# check for no varriation
	if np.var(sampComb) == 0: return(np.nan,np.nan)
	# get ranks, possibly tied
	ranks,tieAdj = tiedRank(sampComb)
	ranks = ranks +1
	# get the stat
	xRank = ranks[:n]
	w = np.sum(xRank)
	# calculate the z statistic
	cor = 2*tieAdj/( (xLen+yLen)*(xLen+yLen-1) )
	wVar = xLen*yLen*((xLen+yLen+1) - cor)/12.
	wMean = n*(xLen+yLen+1)/2.
	wC = w-wMean
	z = (wC - .5 * np.sign(wC))/np.sqrt(wVar)

	# more complex methods exist, but I am using this 
	# simple way to deal with small samples
	if forceExact or (xLen+yLen)<15:
		allComb = chooseAllComb(ranks,n)
		allCombSum = np.sum(allComb,1)
		allCombLen = len(allCombSum)
		p = np.sum(allCombSum>=w)/float(allCombLen)
	else:
		p = stats.norm.cdf(-z)

	
		
	return(p,z)


def tiedRank(x):
	"""Calculates ranks for vector x while considering 
	ties.  If values are the same then the rank assigned
	is the average of the ranks they would span.
	Also returns the tie adjustment used in some tests.
	"""
	n = len(x)
	# sor the values 
	ind = np.argsort(x)
	xSort = x[ind]
	# get the nan values, at the end
	nNan = np.sum(np.isnan(x))
	xLen = n-nNan
	
		
	# find all ties
	tie = xSort[:xLen-1] == xSort[1:xLen]	
	tieInd = np.arange(xLen-1)[tie]
	tieInd = np.append(tieInd,xLen+2)
	nTies = len(tieInd)

	rankTmp = np.arange(float(n)) # the nan does not seem to work with ints 
	rankTmp[xLen:] = np.nan
	tieAdj = 0
	count = 0
	while count<nTies-1:
		tieStart = tieInd[count]
		nTied = 2
		# count number of ties, multiple ties will have consecuitive numbers
		while tieInd[count+1] == tieInd[count]+1:
			count = count + 1
			nTied = nTied + 1
			
		tieAdj = tieAdj + nTied*(nTied-1)*(nTied+1)/2.
		# compute average 
		rankTmp[tieStart:tieStart+nTied] = np.sum(rankTmp[tieStart:tieStart+nTied])/nTied
		count = count+1

	# put in order
	rank = np.zeros(n)
	rank[ind] = rankTmp
	return(rank,tieAdj)

	

def nck(n, k):
	"""Calculate n choose k
	N!/K!(N-K)!.  No warning for 
	large numbers that will take 
	forever!"""
	k = min(k, n-k)
	if k == 0: return 1
	numer = reduce(op.mul, xrange(n, n-k, -1))
	denom = reduce(op.mul, xrange(1, k+1))
	return numer//denom

def chooseAllComb(V,k):
	"""choose all possible combinations of 
	k items from the vector V
	blows up fast so be careful
	"""
	return(np.array(list(itertools.combinations(V,k))))

def runPairwise(x,y,xType, yType, obsMinWarn=5, obsMinError=1):
	"""Simple function to run a standard pairwise
	test on features of potentially different types
	to identify univariate statistical relationships between x and y.
	We have implemented robust methods when standard,
	well excepted, easy to implement options are available.
	x and y are np arrays of the same size (int float and str allowed,
	where str 'nan' and 'NAN' or other types np.nan is used for missing data.
	xType and yType holds a upper case string indicating type
	N,C or B for numerical (ordered) categorical or binary.
	Tests are determined based on the features being compared:
	N-N = Spearman rank, B-B = Fisher Exact, C-C = ChiSq,
	N-B = Ranks Sum, N-C = Kruskal Wallis.
	The test also determines what is reported in r (correlation, effect size, separability), 
	which, like p, is symmetric.
	N-N = spearman rho correlation coefficient (bounded by -1,1),
	B-B = phi coefficient which is the person coefficient analogue 
	C-C = cramer's V (generalization of phi), which is related to correlation for the chi sq test.
	N-C = as Kruskal Wallis has no simple single metric for effect size,
	we will use a measure of separability similar to multiclass LDA
	sqrt(variance of class means / variance of samples) (note this is nonrobust)
	which ranges from zero with no separability to 1
	N-B = effect size determined by z/sqrt(N), which is related to t-test correlation coeff,




	"""
	n = len(x)
	if n!=len(y): raise ValueError('x and y must be same size')

	obsMinWarnFlag=False
	obsMinErrorFlag=False

	if xType=='N' and yType=='N':
		# numerical - numerical, Spearman Rank
		# if its a string array, need to convert to float:
		if type(x[0])==np.string_:x = np.array(x,dtype=float)
		if type(y[0])==np.string_:y = np.array(y,dtype=float)
		r,p = stats.spearmanr(x,y)

	elif ((xType=='N' and (yType=='C' or yType=='B')) or ((xType=='C' or xType=='B') and yType=='N')):
		# numerical - categorical/binary, Kruskal Wallis
		# find groups:
		cat = False 
		if xType=='N':
			if type(x[0])==np.string_:x = np.array(x,dtype=float)
			values = x
			labels = y
			if yType=='C':cat=True
		else:
			if type(y[0])==np.string_:y = np.array(y,dtype=float)
			values = y
			labels = x
			if xType=='C':cat=True

		groups,_ = getGroups(values,labels)
		m = len(groups)
		for i in range(m):
			if len(groups[i])<obsMinWarn:obsMinWarnFlag==True
			if len(groups[i])<obsMinError:obsMinErrorFlag==True

		if obsMinErrorFlag==True:
			p = np.nan
			r = np.nan
		elif cat==True:
			# numerical - categorical, Kruskal Wallis
			h,p = stats.mstats.kruskalwallis(*groups)
			# calculate our ad hoc effect size
			classMeans = np.array([])
			for i in range(m):
				classMeans = np.append(classMeans,np.mean(groups[i]))
			r = np.sqrt(np.var(classMeans)/np.var(values))
		else:
			# numerical - binary, Rank Sum
			p,z = rankSum(groups[0],groups[1])
			r = z / np.sqrt(len(values))

		 

	elif (xType=='C' or xType=='B') and (yType=='C' or yType=='B'):
		# categorical variables
		conTable, xCat, yCat = makeConTable(x,y)
		# check table for limits on observation counts
		for i in range(len(xCat)):
			for j in range(len(yCat)):
				if conTable[i,j]<obsMinError:obsMinErrorFlag==True
				if conTable[i,j]<obsMinWarn:obsMinWarnFlag==True

		if xType=='B' and yType=='B':
			# special case of binary - binary, Fisher exact:
			odds,p = stats.fisher_exact(conTable)
			#r = np.log10(odds) no longer reporting odds ratio
			chi,_,_,_ = stats.chi2_contingency(conTable)
			r = np.sqrt(chi/n)
		else:
			# chi squared test
			chi,p,_,_ = stats.chi2_contingency(conTable)
			rows,cols = conTable.shape
			nMin = min([rows,cols])
			r = np.sqrt((chi/n)/(nMin-1))




	else:
		raise ValueError('prefix on labels is incorrect')

	return(r,p,obsMinErrorFlag)

	
	
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





	
