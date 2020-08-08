
import os
import sys
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from random import *

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


df = pd.read_csv("./GM12878.tsv", sep = "\t")
print(df.head())
adjacent_diff = []
random_diff = []

for i in range(200):
	ind1 = int(len(df)*random())
	ind2 = int(len(df)*random())
	if ind1 == ind2:
		continue
	random_diff.append(abs(df.loc[ind1, "QNed_SDOCs"] - df.loc[ind2, "QNed_SDOCs"]))
for i in range(len(df) - 1):
	if df.loc[i,"TADs"].split("\t")[2] != df.loc[i+1,"TADs"].split("\t")[1]:
		continue
	adjacent_diff.append(abs(df.loc[i, "QNed_SDOCs"] - df.loc[i+1, "QNed_SDOCs"]))

print(np.mean(adjacent_diff), np.mean(random_diff))


color_bar = [[.6,.6,.6],[.6,.6,.6]]

plt.figure(figsize=(2, 3), dpi=600)
plt.rc('font',family='Arial')
plt.rc('font',size = 8)

boxprops = {'linewidth':1.5 , 'color': 'black'}
whiskerprops = {'linestyle': '--', 'linewidth':1.1, 'color': 'black'}
medianprops = {'linewidth':1 , 'color': 'black', 'drawstyle' : 'steps'}
bp4 = plt.boxplot([adjacent_diff, random_diff], positions = [i for i in range(2)] ,showfliers=False, meanline =False, patch_artist = True, whiskerprops = whiskerprops, medianprops = medianprops,widths = 0.68,whis = [5, 95])
for patch, color in zip(bp4['boxes'], color_bar):
    patch.set_facecolor(color)

plt.xticks([0,1], ["Adjacent TADs", "Random TADs"])
plt.ylabel("Difference in SDOC")
plt.tight_layout()
plt.savefig("./random_vs_adjacent.png")
plt.savefig("./random_vs_adjacent.pdf")
plt.close()

print(stats.ttest_ind(adjacent_diff, random_diff))