import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s')

################
#get data
SDOA = []
DS = []
AP = []

f = open("./3.1-get_spearman_corr_SDOA_vs_exp/1-result_spearman_corr_SDOA_vs_exp")   		
lines=f.readlines() 
nrow = len(lines)					
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    SDOA.append(float(L[2]))
f.close()
f = open("./3.2-get_spearman_corr_Ds_vs_exp/1-result_spearman_corr_Dscore_vs_exp")   		
lines=f.readlines() 
nrow = len(lines)					
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    DS.append(float(L[2]))
f.close()

f = open("./3.3-get_spearman_corr_AP_vs_exp/1-result_spearman_corr_AP_vs_exp")   		
lines=f.readlines() 
nrow = len(lines)					
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    AP.append(float(L[2]))
f.close()

#plot
colors =[[0.4,0.65,0.85],
[0.4,0.8,0.65],
[0.9,0.55,0.4]]


plt.figure(figsize=(7, 2), dpi=300)

plt.rc('font',family='Arial')
plt.rc('font',size = 9)


plt.bar([i*3-0.7 for i in range(16)], SDOA, color = colors[0],  width = 0.7)
plt.bar([i*3 for i in range(16)], DS, color = colors[1], width = 0.7)
plt.bar([i*3+0.7 for i in range(16)], AP, color = colors[2], width = 0.7)

plt.title("Correlation of each parameter to gene expression level in each dataset")
plt.ylabel("Spearman correlation")
plt.savefig("./4-barplot.png")
plt.savefig("./4-barplot.eps")

print("SDOA:",SDOA)
print("DS:",DS)
print("AP:",AP)