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

    
celltypes = ["GM12878","IMR90","K562","HUVEC"]
resolution = ["5kb","10kb","25kb","50kb"]
intersect_ratio_thres = 0.7

mkdir("5-tmp")
################
mkdir("./5-tmp/5kb")
mkdir("./5-tmp/10kb")
mkdir("./5-tmp/25kb")
mkdir("./5-tmp/50kb")

########### part1 SDOA ############

for ct in celltypes:
	for res in resolution:
		os.system("cp ./3.1-get_spearman_corr_SDOA_vs_exp/1-result_TAD_SDOA_exp/"+ct+"_"+res+" ./5-tmp/"+res+"/")


for ct in celltypes:
	for res in resolution:
		os.system("python3 overlap.py ./5-tmp/"+res+"/"+ct+"_"+res+" ./5-tmp/"+res+"/ 0 0")


mkdir("./5-result_conserved_TADs/")


for a in range(len(celltypes)):
	for b in range(a):
		for res in resolution:
			ct1 = celltypes[a]
			ct2 = celltypes[b]

			fout = open("./5-result_conserved_TADs/"+res+"_"+ct1+"_"+ct2,'w')
			f = open("./overlapped/"+ct1+"_"+res+"_"+ct2+"_"+res+".txt")
			lines=f.readlines()
			nrow = len(lines)
			for i in range(len(lines)):
				L = lines[i].strip().split('\t')	
				if len(L) == 7:
					continue
				start1 = int(L[1])
				end1 = int(L[2])
				TAD1_len = end1 - start1
				for j in range(7,len(L)-1,7):
					start2 = int(L[j+1])
					end2 = int(L[j+2])
					intersect_len = min(end1,end2) - max(start1,start2)
					TAD2_len = end2 - start2
					if float(intersect_len) / TAD1_len >= intersect_ratio_thres and float(intersect_len) / TAD2_len >= intersect_ratio_thres:
						fout.write(L[0]+'\t'+L[1]+'\t'+L[2]+'\t'+L[5]+'\t'+L[6]+'\t'+L[j]+'\t'+L[j+1]+'\t'+L[j+2]+'\t'+L[j+5]+'\t'+L[j+6]+'\n')
						continue
			f.close()
			fout.close()


######## part2 Dscore #######
for ct in celltypes:
	for res in resolution:
		os.system("cp ./3.2-get_spearman_corr_Ds_vs_exp/1-result_TAD_Dscore_exp/"+ct+"_"+res+" ./5-tmp/"+res+"/")


for ct in celltypes:
	for res in resolution:
		os.system("python3 overlap.py ./5-tmp/"+res+"/"+ct+"_"+res+" ./5-tmp/"+res+"/ 0 0")

mkdir("./5-result_conserved_TADs_Dscore")
for a in range(len(celltypes)):
	for b in range(a):
		for res in resolution:
			ct1 = celltypes[a]
			ct2 = celltypes[b]

			fout = open("./5-result_conserved_TADs_Dscore/"+res+"_"+ct1+"_"+ct2,'w')
			f = open("./overlapped/"+ct1+"_"+res+"_"+ct2+"_"+res+".txt")
			lines=f.readlines()
			nrow = len(lines)
			for i in range(len(lines)):
				L = lines[i].strip().split('\t')	
				if len(L) == 6:
					continue
				start1 = int(L[1])
				end1 = int(L[2])
				TAD1_len = end1 - start1
				for j in range(6,len(L)-1,6):
					start2 = int(L[j+1])
					end2 = int(L[j+2])
					intersect_len = min(end1,end2) - max(start1,start2)
					TAD2_len = end2 - start2
					if float(intersect_len) / TAD1_len >= intersect_ratio_thres and float(intersect_len) / TAD2_len >= intersect_ratio_thres:
						fout.write(L[0]+'\t'+L[1]+'\t'+L[2]+'\t'+L[3]+'\t'+L[5]+'\t'+L[j]+'\t'+L[j+1]+'\t'+L[j+2]+'\t'+L[j+3]+'\t'+L[j+5]+'\n')
						continue
			f.close()
			fout.close()



########## part3 AP ########
for ct in celltypes:
	for res in resolution:
		os.system("cp ./3.3-get_spearman_corr_AP_vs_exp/1-result_TAD_AP_exp/"+ct+"_"+res+" ./5-tmp/"+res+"/")


for ct in celltypes:
	for res in resolution:
		os.system("python3 overlap.py ./5-tmp/"+res+"/"+ct+"_"+res+" ./5-tmp/"+res+"/ 0 0")


mkdir("./5-result_conserved_TADs_AP")
for a in range(len(celltypes)):
	for b in range(a):
		for res in resolution:
			ct1 = celltypes[a]
			ct2 = celltypes[b]

			fout = open("./5-result_conserved_TADs_AP/"+res+"_"+ct1+"_"+ct2,'w')
			f = open("./overlapped/"+ct1+"_"+res+"_"+ct2+"_"+res+".txt")
			lines=f.readlines()
			nrow = len(lines)
			for i in range(len(lines)):
				L = lines[i].strip().split('\t')	
				if len(L) == 6:
					continue
				start1 = int(L[1])
				end1 = int(L[2])
				TAD1_len = end1 - start1
				for j in range(6,len(L)-1,6):
					start2 = int(L[j+1])
					end2 = int(L[j+2])
					intersect_len = min(end1,end2) - max(start1,start2)
					TAD2_len = end2 - start2
					if float(intersect_len) / TAD1_len >= intersect_ratio_thres and float(intersect_len) / TAD2_len >= intersect_ratio_thres:
						fout.write(L[0]+'\t'+L[1]+'\t'+L[2]+'\t'+L[3]+'\t'+L[5]+'\t'+L[j]+'\t'+L[j+1]+'\t'+L[j+2]+'\t'+L[j+3]+'\t'+L[j+5]+'\n')
						continue
			f.close()
			fout.close()
