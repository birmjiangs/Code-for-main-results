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
celltypes = ["GM12878", "IMR90", "K562", "HUVEC"]
for ct in celltypes:

    # do only once
    os.system("python3 ./overlap.py ./data/"+ct+"_10kb ./0-input_compartment/ 0 0 ./0-overlapped")
    
    out_TAD_compartment = []
    f = open("./0-overlapped/"+ct+"_10kb_0-PC1_res100000_"+ct+".txt.txt")   		
    lines=f.readlines() 
    nrow = len(lines)					
    for i in range(len(lines)):
        L = lines[i].strip().split('\t')
        if len(L) == 6:
        	out_TAD_compartment.append(L[0]+'\t'+L[1]+'\t'+L[2]+'\tno_data')
        	continue
        
        A_length = 0
        B_length = 0
        for j in range(9,len(L),6):
        	overlapped_length = min(float(L[2]),float(L[j-1])) - max(float(L[1]),float(L[j-2]))
        	if L[j] == "A":
        		A_length += overlapped_length
        	elif L[j] == "B":
        		B_length += overlapped_length
    
        if A_length + B_length == 0:
        	out_TAD_compartment.append(L[0]+'\t'+L[1]+'\t'+L[2]+'\tno_data')
        	continue
    
        if A_length == 0:
        	out_TAD_compartment.append(L[0]+'\t'+L[1]+'\t'+L[2]+'\tB')
        	continue
        if B_length == 0:
        	out_TAD_compartment.append(L[0]+'\t'+L[1]+'\t'+L[2]+'\tA')
        	continue
    
        if A_length > B_length:
        	out_TAD_compartment.append(L[0]+'\t'+L[1]+'\t'+L[2]+'\tmixed_A')
        if B_length >= A_length:
        	out_TAD_compartment.append(L[0]+'\t'+L[1]+'\t'+L[2]+'\tmixed_B')
    
    f.close()
    
    savetxt("./0-TAD_compartment_"+ct+".tsv", out_TAD_compartment)
    


os.system("python3 ./overlap.py ./data/GM12878_10kb ./data/super_enhancers/GM12878_SE 1 0 ./0-overlapped_SE")
os.system("python3 ./overlap.py ./data/IMR90_10kb ./data/super_enhancers/IMR90_SE 1 0 ./0-overlapped_SE")
os.system("python3 ./overlap.py ./data/K562_10kb ./data/super_enhancers/K562_SE 1 0 ./0-overlapped_SE")
os.system("python3 ./overlap.py ./data/HUVEC_10kb ./data/super_enhancers/HUVEC_SE 1 0 ./0-overlapped_SE")