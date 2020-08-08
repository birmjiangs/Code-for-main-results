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
  
################
celltypes = ["HSC","MPP","CLP","ETP","DN2","DN3","DN4","DP"]

chart = []
TAD_ct_2_SDOC = {}
TAD_ct_2_DHS = {}
TAD_ct_2_TAD_volume = {}
TAD_2_count = {}
for ct in celltypes:
    f = open("./4-QNed_SDOC/"+ct)   		
    lines=f.readlines() 
    nrow = len(lines)					
    for i in range(len(lines)):
        L = lines[i].strip().split('\t')
        TAD = L[0]+'\t'+L[1]+'\t'+L[2]

        TAD_ct_2_DHS[TAD+'\t'+ct] = L[3]
        TAD_ct_2_TAD_volume[TAD+'\t'+ct] = L[4]
        TAD_ct_2_SDOC[TAD+'\t'+ct] = L[5]
        if TAD not in TAD_2_count:
            TAD_2_count[TAD] = 1
        else:
            TAD_2_count[TAD] += 1
    f.close()

for TAD in TAD_2_count:
    if TAD_2_count[TAD] != len(celltypes):
        continue

    chart.append(TAD)

    for ct in celltypes:
        chart[-1] = chart[-1] + '\t' + TAD_ct_2_DHS[TAD+'\t'+ct]

    for ct in celltypes:
        chart[-1] = chart[-1] + '\t' + TAD_ct_2_TAD_volume[TAD+'\t'+ct]
        
    for ct in celltypes:
        chart[-1] = chart[-1] + '\t' + TAD_ct_2_SDOC[TAD+'\t'+ct]
        
chart = sort_list(chart)

#normalize DHS and TAD volume
#for i in range(3,3+3*len(celltypes))

savetxt("./1-chart_8_stages.tsv",chart)





######## do all QNed chart:
chart = []
TAD_ct_2_SDOC = {}
TAD_ct_2_DHS = {}
TAD_ct_2_TAD_volume = {}
TAD_2_count = {}
for ct in celltypes:
    f = open("./4-QNed_all/"+ct)          
    lines=f.readlines() 
    nrow = len(lines)                   
    for i in range(len(lines)):
        L = lines[i].strip().split('\t')
        TAD = L[0]+'\t'+L[1]+'\t'+L[2]

        TAD_ct_2_DHS[TAD+'\t'+ct] = L[3]
        TAD_ct_2_TAD_volume[TAD+'\t'+ct] = L[4]
        TAD_ct_2_SDOC[TAD+'\t'+ct] = L[5]
        if TAD not in TAD_2_count:
            TAD_2_count[TAD] = 1
        else:
            TAD_2_count[TAD] += 1
    f.close()

for TAD in TAD_2_count:
    if TAD_2_count[TAD] != len(celltypes):
        continue

    chart.append(TAD)

    for ct in celltypes:
        chart[-1] = chart[-1] + '\t' + TAD_ct_2_DHS[TAD+'\t'+ct]

    for ct in celltypes:
        chart[-1] = chart[-1] + '\t' + TAD_ct_2_TAD_volume[TAD+'\t'+ct]
        
    for ct in celltypes:
        chart[-1] = chart[-1] + '\t' + TAD_ct_2_SDOC[TAD+'\t'+ct]
        
chart = sort_list(chart)

#normalize DHS and TAD volume
#for i in range(3,3+3*len(celltypes))

savetxt("./1-chart_8_stages_all_QNed.tsv",chart)