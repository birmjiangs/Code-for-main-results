import os
import sys
import numpy as np
import matplotlib.pyplot as plt


colors =[[0.4,0.65,0.85],
[0.4,0.8,0.65],
[0.9,0.55,0.4]]

def mkdir(path):
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False
################
ticks = []
SDOA = []
Ds = []
AP = []
ct1_ct2_res_2_3results = {}
f = open("7-MANUAL_chart.txt")   		
lines=f.readlines() 
nrow = len(lines)					
for i in range(len(lines)):
	L = lines[i].strip().split('\t')
	#ticks.append(L[0]+'-'+L[1]+"("+L[2]+")")
	ticks.append(L[0]+'-'+L[1])
	SDOA.append(float(L[3]))
	Ds.append(float(L[4]))
	AP.append(float(L[5]))
	ct1_ct2_res_2_3results[L[0]+L[1]+L[2]] = [float(L[3]),float(L[4]),float(L[5])]

f.close()

plt.figure(figsize=(8, 3), dpi=600)
plt.rc('font',family='Arial')
plt.rc('font',size = 7)

plt.bar([i*3-0.7 for i in range(24)], SDOA, color = colors[0], width = 0.7)
plt.bar([i*3 for i in range(24)], Ds, color = colors[1], width = 0.7)
plt.bar([i*3+0.7 for i in range(24)], AP, color = colors[2], width = 0.7)

plt.xticks([i*3 for i in range(24)], ticks,rotation = 45)
plt.ylabel("Spearman correlation")
plt.title("Correlation of change in each parameter to gene expression alteration between cell types")
plt.tight_layout()

plt.savefig("./8-barplot.png")
plt.savefig("./8-barplot.eps")
plt.close()

print("SDOA:",SDOA)
print("DS:",Ds)
print("AP:",AP)


