import os
import sys
import numpy as np
import scipy.stats

def mkdir(path):
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False

def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s')

mkdir("./6-correlation_celltype_diff")

fout = open("./6-correlation_celltype_diff/SDOA",'w')
files = os.listdir("./5-result_conserved_TADs")			#get file list of dir
file_name = []  
for fl in files:
	file_name.append(fl)

	f = open("./5-result_conserved_TADs/"+str(fl))   		#open each file
	lines=f.readlines() 

	SDOA_diff = []
	exp_diff = []
	nrow = len(lines)					#get each line
	for i in range(len(lines)):
		L = lines[i].strip().split('\t')
		exp_diff.append(float(L[9]) - float(L[4]))
		SDOA_diff.append(float(L[8]) - float(L[3]))

	corr, pvalue = scipy.stats.spearmanr(exp_diff, SDOA_diff)
	#corr2, pvalue2 = scipy.stats.pearsonr(exp_diff, SDOA_diff)
	fout.write(fl+'\t'+str(corr)+'\t'+str(pvalue)+'\n')

	f.close()
fout.close()


fout = open("./6-correlation_celltype_diff/Dscore",'w')
files = os.listdir("./5-result_conserved_TADs_Dscore")			#get file list of dir
file_name = []  
for fl in files:
	file_name.append(fl)

	f = open("./5-result_conserved_TADs_Dscore/"+str(fl))   		#open each file
	lines=f.readlines() 

	SDOA_diff = []
	exp_diff = []
	nrow = len(lines)					#get each line
	for i in range(len(lines)):
		L = lines[i].strip().split('\t')
		exp_diff.append(float(L[9]) - float(L[4]))
		SDOA_diff.append(float(L[8]) - float(L[3]))

	corr, pvalue = scipy.stats.spearmanr(exp_diff, SDOA_diff)
	#corr2, pvalue2 = scipy.stats.pearsonr(exp_diff, SDOA_diff)
	fout.write(fl+'\t'+str(corr)+'\t'+str(pvalue)+'\n')

	f.close()
fout.close()


fout = open("./6-correlation_celltype_diff/AP",'w')
files = os.listdir("./5-result_conserved_TADs_AP")			#get file list of dir
file_name = []  
for fl in files:
	file_name.append(fl)

	f = open("./5-result_conserved_TADs_AP/"+str(fl))   		#open each file
	lines=f.readlines() 

	SDOA_diff = []
	exp_diff = []
	nrow = len(lines)					#get each line
	for i in range(len(lines)):
		L = lines[i].strip().split('\t')
		exp_diff.append(float(L[9]) - float(L[4]))
		SDOA_diff.append(float(L[8]) - float(L[3]))

	corr, pvalue = scipy.stats.spearmanr(exp_diff, SDOA_diff)
	#corr2, pvalue2 = scipy.stats.pearsonr(exp_diff, SDOA_diff)
	fout.write(fl+'\t'+str(corr)+'\t'+str(pvalue)+'\n')

	f.close()
fout.close()