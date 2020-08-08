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

celltypes = ["GM12878","IMR90","K562","HUVEC"] 
resolution = ["5kb","10kb","25kb","50kb"]

mkdir("./0-AP_combined_sorted")
################
for ct in celltypes:
	for res in resolution:
		out = []
		for n in range(23):
			chrN = "chr"+str(n)
			if n == 0:
				chrN = "chrX"
			try:
				f = open("../0-calculate_all_APs/"+res+"/3-result_TAD_AP/"+ct+"_"+chrN)  
			except:
				continue 		
			lines=f.readlines() 
			nrow = len(lines)					
			for i in range(len(lines)):
				L = lines[i].strip().split('\t')		
				out.append("\t".join(L[0:5]))

			f.close()
		out = sort_list(out)
		savetxt("./0-AP_combined_sorted/"+ct+"_"+res+".txt", out)
