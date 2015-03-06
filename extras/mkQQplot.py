#import genWF

pwOutName='genPW1_data_GNMC_Trans_PW_20150206_vs_data_CLIN_Critical_Phenotype_20150213_out'
outDir = '/titan/ITMI1/projects/gamcop/analysis/gnmc/transcript_20150216'

#_=genWF.splitPWResults(outDir+'/'+pwOutName+'.dat',outDir,1)

import matplotlib.pyplot as plt
import numpy as np

famList = ['M','F','NB']




for fam in famList:
	finPath = outDir+'/'+pwOutName+'_subset_'+fam+'_sorted.dat'
	fin = open(finPath)
	pList = np.array([line.strip().split('\t')[5] for line in fin])
	n = len(pList)
	expPList = -1*np.log10(np.arange(1,0,-1.0/n))
	plt.plot(expPList,pList,'o',label=fam)
	fin.close()

plt.plot([expPList[0],expPList[-1]],[expPList[0],expPList[-1]])
plt.xlabel('expected -logP')
plt.ylabel('observed -logP')
plt.legend(loc=2)

plt.savefig(outDir+'/QQ_log_GNMC_Trans_PW_20150206.png',format='png')
plt.clf()
plt.close()



