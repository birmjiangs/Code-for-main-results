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


ct_2_TPM_col = {"GM12878":5,"IMR90":9,"K562":11,"HUVEC":7}
celltypes = ["GM12878","IMR90","K562","HUVEC"]
resolution = ["5kb","10kb","25kb","50kb"]
mkdir("./7-result_corr")


mkdir("./7-result_ct_diff_each_gene_SDOA/")
fout0 = open("./7-result_corr/SDOA",'w')
for a in range(len(celltypes)):
	for b in range(a):
		for res in resolution:
			ct1 = celltypes[a]
			ct2 = celltypes[b]
			gene_2_TPM_SDOA_ct1 = {}
			gene_2_TPM_SDOA_ct2 = {}

			f = open("./3.1-get_spearman_corr_SDOA_vs_exp/overlapped/"+ct1+"_"+res+"/"+ct1+"_gene_TPM_chart.txt")   		
			lines=f.readlines() 
			nrow = len(lines)					
			for i in range(len(lines)):
				L = lines[i].strip().split('\t')	
				SDOA = float(L[5])
				for j in range(6,len(L),13):
					TPM = (float(L[j+ct_2_TPM_col[ct1]])+float(L[j+ct_2_TPM_col[ct1]+1])) / 2
					gene = L[j+4]
					gene_2_TPM_SDOA_ct1[gene] = str(SDOA)+'\t'+str(TPM)
			f.close()
					

			f = open("./3.1-get_spearman_corr_SDOA_vs_exp/overlapped/"+ct2+"_"+res+"/"+ct2+"_gene_TPM_chart.txt")   		
			lines=f.readlines() 
			nrow = len(lines)					
			for i in range(len(lines)):
				L = lines[i].strip().split('\t')	
				SDOA = float(L[5])
				for j in range(6,len(L),13):
					TPM = (float(L[j+ct_2_TPM_col[ct2]])+float(L[j+ct_2_TPM_col[ct2]+1])) / 2
					gene = L[j+4]
					gene_2_TPM_SDOA_ct2[gene] = str(SDOA)+'\t'+str(TPM)
			f.close()


			SDOA_diffs = []
			log2_TPM_fcs = []
			fout = open("./7-result_ct_diff_each_gene_SDOA/"+ct1+"_"+ct2+"_"+res,'w')
			for gene in gene_2_TPM_SDOA_ct1:
				if gene not in gene_2_TPM_SDOA_ct2:
					continue
				SDOA_TPM_ct1 = list(map(float,gene_2_TPM_SDOA_ct1[gene].split("\t")))
				SDOA_TPM_ct2 = list(map(float,gene_2_TPM_SDOA_ct2[gene].split("\t")))
				log2_TPM_fc = np.log((SDOA_TPM_ct2[1]+1)/(SDOA_TPM_ct1[1]+1)) / np.log(2)
				SDOA_diff = SDOA_TPM_ct2[0] - SDOA_TPM_ct1[0]
				fout.write(gene+'\t'+str(SDOA_diff)+'\t'+str(log2_TPM_fc)+'\n')
				SDOA_diffs.append(SDOA_diff)
				log2_TPM_fcs.append(log2_TPM_fc)
			fout.close()
			corr,pv = stats.spearmanr(SDOA_diffs,log2_TPM_fcs)
			fout0.write(ct1+'\t'+ct2+'\t'+res+'\t'+str(corr)+'\t'+str(pv)+'\n')
fout0.close()





mkdir("./7-result_ct_diff_each_gene_Ds/")
fout0 = open("./7-result_corr/Ds",'w')
for a in range(len(celltypes)):
	for b in range(a):
		for res in resolution:
			ct1 = celltypes[a]
			ct2 = celltypes[b]
			gene_2_TPM_SDOA_ct1 = {}
			gene_2_TPM_SDOA_ct2 = {}

			f = open("./3.2-get_spearman_corr_Ds_vs_exp/overlapped/"+ct1+"_"+res+"/"+ct1+"_"+res+".txt_gene_TPM_chart.txt")   		
			lines=f.readlines() 
			nrow = len(lines)					
			for i in range(len(lines)):
				L = lines[i].strip().split('\t')	
				SDOA = float(L[3])
				for j in range(5,len(L),13):
					TPM = (float(L[j+ct_2_TPM_col[ct1]])+float(L[j+ct_2_TPM_col[ct1]+1])) / 2
					gene = L[j+4]
					gene_2_TPM_SDOA_ct1[gene] = str(SDOA)+'\t'+str(TPM)
			f.close()
					

			f = open("./3.2-get_spearman_corr_Ds_vs_exp/overlapped/"+ct2+"_"+res+"/"+ct2+"_"+res+".txt_gene_TPM_chart.txt")   		 		
			lines=f.readlines() 
			nrow = len(lines)					
			for i in range(len(lines)):
				L = lines[i].strip().split('\t')	
				SDOA = float(L[3])
				for j in range(5,len(L),13):
					TPM = (float(L[j+ct_2_TPM_col[ct2]])+float(L[j+ct_2_TPM_col[ct2]+1])) / 2
					gene = L[j+4]
					gene_2_TPM_SDOA_ct2[gene] = str(SDOA)+'\t'+str(TPM)
			f.close()


			SDOA_diffs = []
			log2_TPM_fcs = []
			fout = open("./7-result_ct_diff_each_gene_Ds/"+ct1+"_"+ct2+"_"+res,'w')
			for gene in gene_2_TPM_SDOA_ct1:
				if gene not in gene_2_TPM_SDOA_ct2:
					continue
				SDOA_TPM_ct1 = list(map(float,gene_2_TPM_SDOA_ct1[gene].split("\t")))
				SDOA_TPM_ct2 = list(map(float,gene_2_TPM_SDOA_ct2[gene].split("\t")))
				log2_TPM_fc = np.log((SDOA_TPM_ct2[1]+1)/(SDOA_TPM_ct1[1]+1)) / np.log(2)
				SDOA_diff = SDOA_TPM_ct2[0] - SDOA_TPM_ct1[0]
				fout.write(gene+'\t'+str(SDOA_diff)+'\t'+str(log2_TPM_fc)+'\n')
				SDOA_diffs.append(SDOA_diff)
				log2_TPM_fcs.append(log2_TPM_fc)
			fout.close()
			corr,pv = stats.spearmanr(SDOA_diffs,log2_TPM_fcs)
			fout0.write(ct1+'\t'+ct2+'\t'+res+'\t'+str(corr)+'\t'+str(pv)+'\n')

fout0.close()



mkdir("./7-result_ct_diff_each_gene_AP/")
fout0 = open("./7-result_corr/AP",'w')
for a in range(len(celltypes)):
	for b in range(a):
		for res in resolution:
			ct1 = celltypes[a]
			ct2 = celltypes[b]
			gene_2_TPM_SDOA_ct1 = {}
			gene_2_TPM_SDOA_ct2 = {}

			f = open("./3.3-get_spearman_corr_AP_vs_exp/overlapped/"+ct1+"_"+res+"/"+ct1+"_"+res+".txt_gene_TPM_chart.txt")		
			lines=f.readlines() 
			nrow = len(lines)					
			for i in range(len(lines)):
				L = lines[i].strip().split('\t')	
				SDOA = float(L[3])
				for j in range(5,len(L),13):
					TPM = (float(L[j+ct_2_TPM_col[ct1]])+float(L[j+ct_2_TPM_col[ct1]+1])) / 2
					gene = L[j+4]
					gene_2_TPM_SDOA_ct1[gene] = str(SDOA)+'\t'+str(TPM)
			f.close()
					

			f = open("./3.3-get_spearman_corr_AP_vs_exp/overlapped/"+ct2+"_"+res+"/"+ct2+"_"+res+".txt_gene_TPM_chart.txt")	
			lines=f.readlines() 
			nrow = len(lines)					
			for i in range(len(lines)):
				L = lines[i].strip().split('\t')	
				SDOA = float(L[3])
				for j in range(5,len(L),13):
					TPM = (float(L[j+ct_2_TPM_col[ct2]])+float(L[j+ct_2_TPM_col[ct2]+1])) / 2
					gene = L[j+4]
					gene_2_TPM_SDOA_ct2[gene] = str(SDOA)+'\t'+str(TPM)
			f.close()

			SDOA_diffs = []
			log2_TPM_fcs = []
			fout = open("./7-result_ct_diff_each_gene_AP/"+ct1+"_"+ct2+"_"+res,'w')
			for gene in gene_2_TPM_SDOA_ct1:
				if gene not in gene_2_TPM_SDOA_ct2:
					continue
				SDOA_TPM_ct1 = list(map(float,gene_2_TPM_SDOA_ct1[gene].split("\t")))
				SDOA_TPM_ct2 = list(map(float,gene_2_TPM_SDOA_ct2[gene].split("\t")))
				log2_TPM_fc = np.log((SDOA_TPM_ct2[1]+1)/(SDOA_TPM_ct1[1]+1)) / np.log(2)
				SDOA_diff = SDOA_TPM_ct2[0] - SDOA_TPM_ct1[0]
				#omit gene where AP = 0 in either cell type
				if SDOA_TPM_ct1[0] == 0 or SDOA_TPM_ct2[0] == 0:
					continue
				fout.write(gene+'\t'+str(SDOA_diff)+'\t'+str(log2_TPM_fc)+'\n')
				SDOA_diffs.append(SDOA_diff)
				log2_TPM_fcs.append(log2_TPM_fc)
			fout.close()
			corr,pv = stats.spearmanr(SDOA_diffs,log2_TPM_fcs)
			fout0.write(ct1+'\t'+ct2+'\t'+res+'\t'+str(corr)+'\t'+str(pv)+'\n')
fout0.close()

