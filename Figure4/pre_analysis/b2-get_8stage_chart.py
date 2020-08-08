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

def sort_list(list_in): #use when using '\t' to seperate items
    list_out = sorted(list_in, key=lambda items: int(items.split('\t')[1]))
    list_out = sorted(list_out, key=lambda items: items.split('\t')[0])
    return list_out
  
################
cell_stages = ["HSC","MPP", "CLP","ETP","DN2","DN3","DN4","DP"]


##### get common TADs：
TAD_2_count_celltype = {}
for stage in cell_stages:
    f = open("./b1-inter_TAD_HiC_mean_count/"+stage)        
    lines=f.readlines() 
    nrow = len(lines)                   
    for i in range(len(lines)):
        L = lines[i].strip().split('\t')    
        TAD = L[0] + '\t' + L[1] + '\t' + L[2]
        if TAD not in TAD_2_count_celltype:
            TAD_2_count_celltype[TAD] = 1
        else:
            TAD_2_count_celltype[TAD] += 1

    f.close()

#############################################

#### read inter TAD count file
TAD_2_chart_SDOC = {}
TAD_2_chart_interTAD_contact = {}
for c in range(len(cell_stages)):
    stage = cell_stages[c]
    f = open("./b1-inter_TAD_HiC_mean_count/"+stage)   		
    lines=f.readlines() 
    nrow = len(lines)					
    for i in range(len(lines)):
        L = lines[i].strip().split('\t')    
        
        TAD = L[0] + '\t' + L[1] + '\t' + L[2]
        if TAD_2_count_celltype[TAD] != 8:  # not common TAD in all cell stages
            continue


        SDOC = L[3]
        mean_inter_TAD_contact = L[4]

        if c == 0:
            TAD_2_chart_SDOC[TAD] = SDOC
            TAD_2_chart_interTAD_contact[TAD] = mean_inter_TAD_contact
        else:
            print(i,len(TAD_2_chart_SDOC))
            TAD_2_chart_SDOC[TAD] = TAD_2_chart_SDOC[TAD] + '\t' + SDOC
            TAD_2_chart_interTAD_contact[TAD] = TAD_2_chart_interTAD_contact[TAD] + '\t' + mean_inter_TAD_contact
    
    f.close()
#################
### make chart

chart_SDOC = []
chart_interTAD_contact = []

for TAD in TAD_2_chart_SDOC:
    chart_SDOC.append(TAD+'\t'+TAD_2_chart_SDOC[TAD])
    chart_interTAD_contact.append(TAD+'\t'+TAD_2_chart_interTAD_contact[TAD])

chart_SDOC = sort_list(chart_SDOC)
chart_interTAD_contact = sort_list(chart_interTAD_contact)

chart_interTAD_contact_zscore = []
for i in range(len(chart_SDOC)):
    L = chart_interTAD_contact[i].split("\t")
    chart_interTAD_contact_zscore.append(list(map(float,L[3:])))

chart_interTAD_contact_zscore = np.array(chart_interTAD_contact_zscore)
for i in range(8):
    chart_interTAD_contact_zscore[:,i] = stats.zscore(chart_interTAD_contact_zscore[:,i])

chart_interTAD_contact_zscore_out = []
for i in range(len(chart_interTAD_contact_zscore)):
    L = list(map(str,chart_interTAD_contact_zscore[i,:].tolist()))


    line = chart_SDOC[i].split("\t")
    TAD = line[0] +'\t' +line[1] + '\t' + line[2]
    chart_interTAD_contact_zscore_out.append(TAD+'\t'+"\t".join(L))

mkdir("b2-charts")
savetxt("./b2-charts/chart_SDOC.txt", chart_SDOC)
savetxt('./b2-charts/chart_interTAD_contact.txt', chart_interTAD_contact)
savetxt('./b2-charts/chart_interTAD_contact_zscore.txt', chart_interTAD_contact_zscore_out)


######## get pearson correlation of each TAD
TAD_correlation = []
correlation = []
for i in range(len(chart_SDOC)):
    SDOC = list(map(float,chart_SDOC[i].split("\t")[3:]))
    TAD = chart_SDOC[i].split("\t")[0:3]
    interTAD_contact = list(map(float,chart_interTAD_contact_zscore_out[i].split("\t")[3:]))

    corr, pvalue = stats.pearsonr(SDOC, interTAD_contact)
    correlation.append(corr)
    TAD_correlation.append("\t".join(TAD)+'\t'+str(corr))

mkdir("b2-correlation_result")
savetxt("./b2-correlation_result/pearson_correlation.txt", TAD_correlation)


#### do pearson correlation only in TADs with interTAD contact alteration zscore > 2:
inter_TAD_contact_peak_valley_diff = []
for i in range(len(chart_interTAD_contact_zscore_out)):
    interTAD_contact = list(map(float,chart_interTAD_contact_zscore_out[i].split("\t")[3:]))
    diff = np.max(interTAD_contact) - np.min(interTAD_contact)
    inter_TAD_contact_peak_valley_diff.append(diff)


correlation_high_diff = []
TAD_correlation_high_diff = []
inter_TAD_contact_peak_valley_diff = stats.zscore(inter_TAD_contact_peak_valley_diff)
for i in range(len(inter_TAD_contact_peak_valley_diff)):
    SDOC = list(map(float,chart_SDOC[i].split("\t")[3:]))
    TAD = chart_SDOC[i].split("\t")[0:3]
    interTAD_contact = list(map(float,chart_interTAD_contact_zscore_out[i].split("\t")[3:]))

    if inter_TAD_contact_peak_valley_diff[i] > 2:
        corr, pv = stats.pearsonr(SDOC, interTAD_contact)
        correlation_high_diff.append(corr)
        TAD_correlation_high_diff.append("\t".join(TAD)+'\t'+str(corr))

savetxt("./b2-correlation_result/pearson_correlation_high_diff.txt", TAD_correlation_high_diff)


