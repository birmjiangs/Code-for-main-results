import numpy as np

def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s')

TAD_out = []
extended_TAD_out = [] #for profile
f = open("./1-dataframe.tsv")           
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    if i == 0:
        continue
    L = lines[i].strip().split('\t')
    chrN = L[1][1:]
    start = int(L[2])
    end = int(L[3][:-1])
    TAD = chrN + '\t' + str(start) + '\t' + str(end)
    TAD_out.append(TAD)

    extended_start = 2*start - end
    extended_end =  2*end - start 
    extended_TAD = chrN + '\t' + str(extended_start) + '\t' + str(extended_end)
    extended_TAD_out.append(extended_TAD)
f.close()

savetxt("./pre1-TADs.tsv", TAD_out)
