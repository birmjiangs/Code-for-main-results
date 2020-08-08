import os
import sys
import numpy as np
import scipy.stats as stats
def mkdir(path):
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False

def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s')
################

mkdir("./overlapped")
mkdir("./1-result_TAD_Dscore_exp")

celltypes = ["GM12878","IMR90","K562","HUVEC"]
resolutions = ["5kb","10kb","25kb","50kb"]
ct_2_del_ncol = {"GM12878":0,"IMR90":4,"K562":6,"HUVEC":2}
out = []
for ct in celltypes:
	for res in resolutions:
		print("doing",ct,res)


		log2meanTPM = []
		Dscore = []
		os.system("python3 overlap.py ../data_Dscore/"+ct+"_"+res+".txt ../TSS_TPM/ 0 0 "+"./overlapped/"+ct+"_"+res)  
		
		fout = open("./1-result_TAD_Dscore_exp/"+ct+"_"+res,'w')
		f = open("./overlapped/"+ct+"_"+res+"/"+ct+"_"+res+".txt_gene_TPM_chart.txt")   	 		
		lines=f.readlines() 
		nrow = len(lines)					
		for i in range(len(lines)):
			L = lines[i].strip().split('\t')
			Dscore.append(float(L[3]))
			TAD = "\t".join(L[0:5])
			if len(L) == 5:
				log2meanTPM.append(0)
				fout.write(TAD+'\t0.0\n')
				continue
			mean_TPM = 0
			count_TPM = 0
			for j in range(10,len(L),13):
				mean_TPM += float(L[j+ct_2_del_ncol[ct]]) + float(L[j+ct_2_del_ncol[ct]+1])
				count_TPM += 2

			mean_TPM /= count_TPM
			log2meanTPM.append(np.log(mean_TPM+1)/np.log(2))
			fout.write(TAD+'\t'+str(log2meanTPM[-1])+'\n')
		fout.close()
		f.close()
		out.append(ct+'\t'+res+'\t'+str(stats.spearmanr(log2meanTPM,Dscore)[0])+'\t'+str(stats.spearmanr(log2meanTPM,Dscore)[1]))#+'\t'+str(stats.pearsonr(log2meanTPM,SDOA)[0]))

savetxt("./1-result_spearman_corr_Dscore_vs_exp",out)


		

				