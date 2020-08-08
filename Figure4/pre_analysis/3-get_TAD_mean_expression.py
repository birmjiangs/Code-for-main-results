import os
import sys
import numpy as np


def mkdir(path):
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False

def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s')

def sort_list(list_in): #use when using '\t' to seperate items
    list_out = sorted(list_in, key=lambda items: int(items.split('\t')[1]))
    list_out = sorted(list_out, key=lambda items: items.split('\t')[0])
    return list_out
  
####### 1. get gene position

mm9_gene = []
mm9_tss = []
gene_2_metadata = {}
f = open("./data/refGene.txt")   		
lines=f.readlines() 
nrow = len(lines)					
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    gene = L[12]
    if gene not in gene_2_metadata:
        gene_2_metadata[gene] = L[2]+'\t'+L[4]+'\t'+L[5]+'\t'+L[3]+'\t'+L[12]+'\t'+L[1]
        mm9_gene.append(gene_2_metadata[gene])
        if L[3] == "+":
            start = int(L[4]) - 1000
            end = int(L[4]) + 1000
            mm9_tss.append(L[2]+'\t'+str(start)+'\t'+str(end)+'\t'+L[3]+'\t'+L[12]+'\t'+L[1])
        if L[3] == "-":
            start = int(L[5]) - 1000
            end = int(L[5]) + 1000
            mm9_tss.append(L[2]+'\t'+str(start)+'\t'+str(end)+'\t'+L[3]+'\t'+L[12]+'\t'+L[1])

f.close()

mm9_gene = sort_list(mm9_gene)
mm9_tss = sort_list(mm9_tss)
mkdir("./3-mm9_refGene/")

savetxt("./3-mm9_refGene/refGene.tsv", mm9_gene)
savetxt("./3-mm9_refGene/refGene_tss.tsv", mm9_tss)
####### 2. get_TAD_overlap_gene
out_TADs = []
f = open("1-chart_8_stages.tsv")           
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    out_TADs.append("\t".join(L[0:3]))
f.close()
savetxt("./3-used_TADs", out_TADs)

os.system("python ./overlap.py 3-used_TADs ./3-mm9_refGene/ 0 0 ./3-overlapped")

TAD_2_genes = {}
f = open("./3-overlapped/3-used_TADs_refGene_tss.tsv.txt")
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    TAD = "\t".join(L[0:3])
    

    if len(L)<=3:
        TAD_2_genes[TAD] = ""
        continue

    for j in range(7,len(L),6):
        if TAD not in TAD_2_genes:
            TAD_2_genes[TAD] = L[j]
        else:
            TAD_2_genes[TAD] = TAD_2_genes[TAD] + ',' + L[j]

f.close()

####### 3. read expression data and get TAD mean expression of each stage
celltypes = ["HSC","MPP","CLP","ETP","DN2","DN3","DN4","DP"]

TAD_2_stages_mean_exp = {}
for ct in celltypes:
    gene_2_fpkm = {}

    f = open("./data/RNA-seq/"+ct+".txt")           
    lines=f.readlines() 
    nrow = len(lines)                   
    for i in range(len(lines)):
        L = lines[i].strip().split('\t')

        gene = L[4].split("|")[2]
        gene_2_fpkm[gene] = float(L[1])

    f.close()

    ## comvert FPKM to TPM
    #@@sum_FPKM = 0
    #@@for gene in gene_2_fpkm:
    #@@    sum_FPKM += gene_2_fpkm[gene]
    #@@for gene in gene_2_fpkm:
    #@@    gene_2_fpkm[gene] = gene_2_fpkm[gene] / sum_FPKM * 1000000.0
    #########

    for TAD in TAD_2_genes:
        if TAD_2_genes[TAD] == "":
            if TAD not in TAD_2_stages_mean_exp:
                TAD_2_stages_mean_exp[TAD] = TAD + "\tnan"
            else:
                TAD_2_stages_mean_exp[TAD] = TAD_2_stages_mean_exp[TAD] + "\tnan"
            continue

        #print(TAD_2_genes[TAD])
        exp_sum = 0
        count_gene = 0
        genes_in_TAD = TAD_2_genes[TAD].split(",")
        for gene in genes_in_TAD:
            if gene not in gene_2_fpkm:
                continue
            exp_sum += np.log(gene_2_fpkm[gene]+1)/np.log(2)  #use 1 as pseudocount
            count_gene += 1
            
        if count_gene == 0:  #all genes in TAD does not has expression data avaliable
            exp_mean = "nan"
        else:
            exp_mean = exp_sum / count_gene  #log mean exp (pseudocount = 1)
        if TAD not in TAD_2_stages_mean_exp:
            TAD_2_stages_mean_exp[TAD] = TAD + '\t' + str(exp_mean) 
        else:
            TAD_2_stages_mean_exp[TAD] = TAD_2_stages_mean_exp[TAD] + '\t' + str(exp_mean) 


out = []
for TAD in TAD_2_stages_mean_exp:
    out.append(TAD_2_stages_mean_exp[TAD])

out = sort_list(out)

savetxt("./3-result_TAD_mean_exp_8stages_FPKM", out)
#savetxt("./3-result_TAD_mean_exp_8stages", out)   #if do FPKM-TPM conversion