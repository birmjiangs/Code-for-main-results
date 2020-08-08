
import os
import sys
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression 
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as stats

matplotlib.rcParams['pdf.fonttype'] = 42
np.random.seed(5)

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

def draw_heatmap(data,xticks,ylabels):
 # cmap = cm.afmhot.reverse()
    figure=plt.figure(facecolor='w',figsize=(4,4),dpi=600)
    ax=figure.add_subplot(1,1,1,position=[0.1,0.15,0.8,0.8])
    ax.set_yticks([i for i in range(len(ylabels))])
    ax.set_yticklabels(ylabels)
    ax.set_xticks([i for i in range(len(xticks))])
    ax.set_xticklabels(xticks)
    vmax=data[0][0]
    vmin=data[0][0]
    for i in data:
        for j in i:
            if j>vmax:
                vmax=j
            if j<vmin:
                vmin=j
    map=ax.imshow(data,interpolation='nearest',cmap=cmap_SDOC,aspect='auto',vmin=-3,vmax=3)
    cb=plt.colorbar(mappable=map,cax=None,ax=None,aspect=9,shrink=0.3,ticks=[-3,3])


def draw_corrmap(data,xticks,yticks):
 # cmap = cm.afmhot.reverse()
    figure=plt.figure(facecolor='w',figsize=(3,2.7),dpi=600)
    ax=figure.add_subplot(1,1,1)
    ax.set_xticks([i for i in range(len(xticks))])
    ax.set_yticks([i for i in range(len(yticks))])
    ax.set_yticklabels([i for i in yticks], rotation = 45)
    ax.set_xticklabels([i for i in xticks], rotation = 45)
    vmax=data[0][0]
    vmin=data[0][0]
    for i in data:
        for j in i:
            if j>vmax:
                vmax=j
            if j<vmin:
                vmin=j
    map=ax.imshow(data,interpolation='nearest',cmap="RdBu_r",aspect='auto',vmin=-1,vmax=1)  #magma
    cb=plt.colorbar(mappable=map,cax=None,ax=None,aspect=10,shrink=0.8,ticks=[-1,1])
    plt.tight_layout()

    for y in range(data.shape[0]): #add text
        for x in range(data.shape[1]):
            plt.text(x, y, '%.2f' % data[y, x],
            horizontalalignment='center',
            verticalalignment='center',
            color = "white",
            )


def draw_boxplot_8stages(data_boxplot,file,ylim = [],color_bar = [[0.7,0.7,0.7] for i in range(8)]):
    plt.figure(figsize=(2.4, 3.2), dpi=600)
    plt.rc('font',family='Arial')
    plt.rc('font',size = 7)
    
    boxprops = {'linewidth':1.5 , 'color': 'black'}
    whiskerprops = {'linestyle': '--', 'linewidth':1.1, 'color': 'black'}
    medianprops = {'linewidth':1 , 'color': 'black', 'drawstyle' : 'steps'}
    bp = plt.boxplot(data_boxplot, showfliers=False, meanline =False, patch_artist = True, whiskerprops = whiskerprops, medianprops = medianprops,widths = 0.77,whis = [5, 95])
    for patch, color in zip(bp['boxes'], color_bar):
        patch.set_facecolor(color)
    
    plt.xticks([i+1 for i in range(8)],["HSC","MPP","CLP","ETP","DN2","DN3","DN4","DP"], rotation = 0)
    if ylim != []:
        plt.ylim(ylim)

    #do paired t test for mean expression of TADs in each stage
    pvalues = [0 for i in range(7)]
    for i in range(1,8):
        stat, pvalue = stats.ttest_rel(data_boxplot[:,i],data_boxplot[:,0])
        pvalues.append(pvalue)

        mark = "N.S."
        if pvalue <= 0.05:
            mark = "*"
        if pvalue < 0.001:
            mark = "**"
        if pvalue < 0.00001:
            mark = "***"
        plt.text(i+0.7,np.percentile(data_boxplot[:,i],96),mark)
        #print(data_boxplot[:,i])
        
    for i in range(8):
        plt.scatter(i+1,np.mean(data_boxplot[:,i]),s = 5, linewidth = 0, color = [0.4,1,1], zorder = 10)
    plt.savefig(file)
    plt.close()
    return pvalues

def do_boxplot_exp_fc(d1,d2,d3,d4,d5,d6,file):

    plt.figure(figsize=(2, 3.5), dpi=600)
    plt.rc('font',family='Arial')
    plt.rc('font',size = 8)
    
    boxprops = {'linewidth':1.5 , 'color': 'black'}
    whiskerprops = {'linestyle': '--', 'linewidth':1.1, 'color': 'black'}
    medianprops = {'linewidth':1 , 'color': 'black', 'drawstyle' : 'steps'}
    color_bar_this = [[0.9,0.9,0.9],
    [0.75,0.75,0.75],
    [0.6,0.6,0.6],
    [0.7,0.85,1],
    [0.5,0.7,0.9],
    [0.3,0.55,0.8]]
    bp = plt.boxplot([d1,d2,d3,d4,d5,d6], showfliers=False, meanline =True, patch_artist = True, whiskerprops = whiskerprops, medianprops = medianprops,widths = 0.77,whis = [5, 95])
    for patch, color in zip(bp['boxes'], color_bar_this):
        patch.set_facecolor(color)
    
    plt.xticks([i+1 for i in range(6)],["" for i in range(6)], rotation = 0)

    print("t test results for "+file+":")
    stat, pvalue = stats.mannwhitneyu(d1,d4)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d2,d5)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d3,d6)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d2,d3)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d5,d6)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d1,d3)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d4,d5)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d4,d6)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d5,d6)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d3,d5)
    print(pvalue)
    print("t test end")

    mkdir("./1-figures/exp_fc_boxplot")
    plt.savefig("./1-figures/exp_fc_boxplot/"+file+".png")
    plt.savefig("./1-figures/exp_fc_boxplot/"+file+".pdf")
    plt.close() 

def do_boxplot_exp_fc_SDOC_increase(d1,d2,d3,d4,d5,d6,file):

    plt.figure(figsize=(2, 3.5), dpi=600)
    plt.rc('font',family='Arial')
    plt.rc('font',size = 8)
    
    boxprops = {'linewidth':1.5 , 'color': 'black'}
    whiskerprops = {'linestyle': '--', 'linewidth':1.1, 'color': 'black'}
    medianprops = {'linewidth':1 , 'color': 'black', 'drawstyle' : 'steps'}
    color_bar_this = [[0.9,0.9,0.9],
    [0.75,0.75,0.75],
    [0.6,0.6,0.6],
    [1,0.9,0.6],
    [0.9,0.75,0.45],
    [0.8,0.65,0.35]]
    bp = plt.boxplot([d1,d2,d3,d4,d5,d6], showfliers=False, meanline =True, patch_artist = True, whiskerprops = whiskerprops, medianprops = medianprops,widths = 0.77,whis = [5, 95])
    for patch, color in zip(bp['boxes'], color_bar_this):
        patch.set_facecolor(color)
    
    plt.xticks([i+1 for i in range(6)],["" for i in range(6)], rotation = 0)
   
    print("t test results for "+file+":")
    stat, pvalue = stats.mannwhitneyu(d1,d4)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d2,d5)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d3,d6)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d2,d3)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d5,d6)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d1,d3)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d4,d5)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d4,d6)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d5,d6)
    print(pvalue)
    stat, pvalue = stats.mannwhitneyu(d3,d5)
    print(pvalue)
    print("t test end")

    mkdir("./1-figures/exp_fc_boxplot")
    plt.savefig("./1-figures/exp_fc_boxplot/"+file+".png")
    plt.savefig("./1-figures/exp_fc_boxplot/"+file+".pdf")
    plt.close() 



def draw_violinplot_8stages(data_boxplot,file,ylim = [],color_bar = [[0.7,0.7,0.7] for i in range(8)]):
    plt.figure(figsize=(4, 1.7), dpi=600)
    plt.rc('font',family='Arial')
    plt.rc('font',size = 8)
    
    violin_parts = plt.violinplot(data_boxplot, positions = [i+1 for i in range(8)] ,showmeans=False, showmedians = False, showextrema = False, widths =0.8, points = 1000)
    for vp in violin_parts['bodies']:
        vp.set_facecolor([0.6,0.6,0.6])
        vp.set_edgecolor([0,0,0])
        vp.set_linewidth(0)
        vp.set_alpha(1)
    
    plt.xticks([i+1 for i in range(8)],["HSC","MPP","CLP","ETP","DN2","DN3","DN4","DP"], rotation = 0)
    if ylim != []:
        plt.ylim(ylim)

    #do paired t test for mean expression of TADs in each stage
    pvalues = [0 for i in range(7)]
    #for i in range(1,8):
    #    stat, pvalue = stats.ttest_rel(data_boxplot[:,i],data_boxplot[:,0])
    #    pvalues.append(pvalue)
#
    #    mark = "N.S."
    #    if pvalue <= 0.05:
    #        mark = "*"
    #    if pvalue < 0.001:
    #        mark = "**"
    #    if pvalue < 0.00001:
    #        mark = "***"
    #    plt.text(i+0.7,np.percentile(data_boxplot[:,i],96),mark)
    #    #print(data_boxplot[:,i])
    #    
    for i in range(8):
        #plt.scatter(i+1,np.mean(data_boxplot[:,i]),s = 5, linewidth = 0, color = [0.4,1,1], zorder = 10)
        plt.plot([i+0.8,i+1.2],[np.median(data_boxplot[:,i]),np.median(data_boxplot[:,i])], linewidth = 1.5, color = "black")
    plt.savefig(file)
    plt.close()
    return pvalues


def draw_comparative_stacked_bar_8stages(data, file, ylim = []):

    plt.figure(figsize=(3, 3), dpi=600)
    plt.rc('font',family='Arial')
    plt.rc('font',size = 8)

    sum_len = len(data)
    print(data, np.shape(data))
    print(np.sum(data == 0, axis = 0),  np.shape(np.sum(data == 0, axis = 0)))
    

    data_count_compA = np.sum(data == 0, axis = 0)
    data_count_compB = np.sum(data == 1, axis = 0)
    data_count_A_maj = np.sum(np.abs(data - 0.25) < 0.25, axis = 0)
    data_count_B_maj = np.sum(np.abs(data - 0.75) <= 0.25, axis = 0) - data_count_compB

    stat, pv = stats.ttest_ind(data_count_compA[0:5], data_count_compA[5:8])
    print(pv)

    plt.bar([i for i in range(8)], [data_count_compB[i]/sum_len for i in range(8)], width = 0.75, linewidth = 0, color = [0.0,0.4,0.7],zorder = 0)
    plt.bar([i for i in range(8)], [(data_count_B_maj[i]+data_count_compB[i])/sum_len for i in range(8)], width = 0.75, linewidth = 0, color = [0,1,1],zorder = -1)
    plt.bar([i for i in range(8)], [(data_count_B_maj[i]+data_count_compB[i]+data_count_A_maj[i])/sum_len for i in range(8)], width = 0.75, linewidth = 0, color = [1,0.6,0.3],zorder = -2)
    plt.bar([i for i in range(8)], [1 for i in range(8)], width = 0.75, linewidth = 0, color = [0.7,0,0],zorder = -3)

    plt.xticks([i for i in range(8)],["HSC","MPP","CLP","ETP","DN2","DN3","DN4","DP"], rotation = 0)
    
    if ylim != []:
        plt.ylim(ylim)

    plt.tight_layout()
    plt.savefig(file)
    plt.close()

#--------------------------------------------------------------------------------------------------

mkdir("./1-out/")

celltypes = ["HSC","MPP","CLP","ETP","DN2","DN3","DN4","DP"]


data = []  #main dataset
TADs = []  #all TADs
## read DHS, TAD volume and SDOC (all QNed)
f = open("../pre_analysis/1-chart_8_stages_all_QNed.tsv")           
lines=f.readlines() 
nrow = len(lines)                    
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    L[3:] = list(map(float,L[3:]))
    data.append(["\t".join(L[0:3])]+L[3:])
    TADs.append("\t".join(L[0:3]))
f.close()
## read mean expression 
f = open("../pre_analysis/3-result_TAD_mean_exp_8stages_FPKM")           
lines=f.readlines() 
nrow = len(lines)                    
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    L[1:] = list(map(float,L[1:]))
    data[i] = data[i] + L[3:]
f.close()
## read inter TAD contact (z-scored)
f = open("../pre_analysis/b2-charts/chart_interTAD_contact_zscore.txt")           
lines=f.readlines() 
nrow = len(lines)                    
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    L[1:] = list(map(float,L[1:]))
    data[i] = data[i] + L[3:]
f.close()


#construct titles
title = ["TAD"]
for i in range(8):
    title.append("DHS_"+str(i))
for i in range(8):
    title.append("TAD_volume_"+str(i))
for i in range(8):
    title.append("SDOC_"+str(i))
for i in range(8):
    title.append("mean_expression_"+str(i))
for i in range(8):
    title.append("adjacent_TAD_contact_"+str(i))

#create dataframe of main dataset
df = pd.DataFrame(data,columns = title)
#inspect dataframe
print(df.head())
print(df.describe())
print("print#0", df["TAD"].head())
# do clustering(cluster lineages) using SDOC of all TADs
#pass, use R to get clustering tree

# select high diff SDOC TADs to do heatmap


# read TAD-TAD contact loess normalized chart
TAD_pair_2_contact = {}
f = open("../pre_analysis/b5-TADpair_contact_charts/data_for_clustering.txt")           
lines=f.readlines() 
nrow = len(lines)                    
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    TAD_pair = "\t".join(L[0:6])
    data = list(map(float,L[7:]))
    TAD_pair_2_contact[TAD_pair] = data
f.close()
print("done read TAD contact chart")

df["SDOC_std"] = np.std(np.array(df.loc[:,"SDOC_0":"SDOC_7"]),axis = 1)

mkdir("./1-figures")

#plot the high diff SDOC TADs 
SDOC_stddev_cutoff = 0.3
SDOC_stddev = np.array(sorted(df["SDOC_std"]))
colors = np.array([[0.7,0.7,0.7] for i in range(len(SDOC_stddev))])
colors[SDOC_stddev>SDOC_stddev_cutoff] = [1,0,0]

plt.scatter([i for i in range(len(colors))],sorted(df["SDOC_std"]), color = colors,s=4)
plt.plot([0,len(df)],[SDOC_stddev_cutoff,SDOC_stddev_cutoff],'--',color = 'black')
plt.savefig("./1-figures/SDOC_stddev_sorted.png")
plt.savefig("./1-figures/SDOC_stddev_sorted.pdf")
plt.close()

#extract high diff SDOC TADs from df:
df_high_SDOC_std = df.loc[df["SDOC_std"]>SDOC_stddev_cutoff,]
df_high_SDOC_std.to_csv("./1-dataframe_high_SDOC_std.tsv",sep = "\t")
#print(df_high_SDOC_std)
#print(type(df),type(df_high_SDOC_std))

# do clustering
data_for_clustering = df_high_SDOC_std.loc[:,"SDOC_0":"SDOC_7"]
corr_mat = np.corrcoef(data_for_clustering)

pred = KMeans(n_clusters=3, random_state=100).fit_predict(corr_mat)
clustering_result = pred

df_high_SDOC_std["cluster"] = clustering_result
df_high_SDOC_std = df_high_SDOC_std.sort_values(by = "cluster")

SDOC_to_plot = df_high_SDOC_std.loc[:,"SDOC_0":"SDOC_7"]


### draw heatmap
cdict_SDOC = {'blue': ((0.0, 0.9, 0.9),
                 (0.3, 0.8, 0.8),
                 (0.5, 0.0, 0.0),
                 (0.8, 0.0, 0.0),
                 (1.0, 0, 0)),
         'green': ((0.0, 0.4, 0.4),
                 (0.3, 0.35, 0.35),
                 (0.5, 0, 0),
                 (0.8, 0.8, 0.8),
                 (1.0, 1, 1)),
         'red': ((0.0, 0.2, 0.2),
                 (0.3, 0.18, 0.18),
                 (0.5, 0, 0),
                 (0.8, 0.8, 0.8),
                 (1.0, 1, 1))}

cdict_expression = {'blue': ((0.0, 0.7, 0.7),
                 (0.3, 0.8, 0.8),
                 (0.5, 1, 1),
                 (0.8, 0.0, 0.0),
                 (1.0, 0, 0)),
         'green': ((0.0, 0.4, 0.4),
                 (0.3, 0.5, 0.5),
                 (0.5, 1, 1),
                 (0.8, 0.3, 0.3),
                 (1.0, 0, 0)),
         'red': ((0.0, 0, 0),
                 (0.3, 0, 0),
                 (0.5, 1,1),
                 (0.8, 0.9, 0.9),
                 (1.0, 0.7, 0.7))}

cmap_SDOC = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict_SDOC,256)
cmap_exp = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict_expression,256)


#plot SDOC
plt.figure(figsize=(1, 1), dpi=600)
plt.rc('font',family='Arial')
plt.rc('font',size = 8)
heatmap_data = SDOC_to_plot.to_numpy()
for i in range(len(heatmap_data)):
    heatmap_data[i,:] = stats.zscore(heatmap_data[i,:])
draw_heatmap(heatmap_data,["HSC","MPP","CLP","ETP","DN2","DN3","DN4","DP"],[])
plt.savefig("./1-figures/SDOC_heatmap.png")
plt.savefig("./1-figures/SDOC_heatmap.pdf")
plt.close() 

mkdir("./1-results")
savetxt("./1-results/SDOC_heatmaps_data",SDOC_to_plot)



#1. get all SDOC decrease TADs from the clustering result
SDOC_decrease_TADs = df_high_SDOC_std[df_high_SDOC_std["cluster"] == 1]["TAD"].values.tolist()
SDOC_increase_TADs = df_high_SDOC_std[df_high_SDOC_std["cluster"] == 2]["TAD"].values.tolist()
SDOC_decrease_TADs = sort_list(SDOC_decrease_TADs)

print("total amount of SDOC decreasing TADs:", len(SDOC_decrease_TADs))

#2. boxplot pair-wise loess normalized contact in 8 stages
TAD_pair_contact_8stages = []
TAD_pairs_SDOC_decreasing = []
TAD_pairs_SDOC_increasing = [] ###?
for i in range(len(SDOC_decrease_TADs)):
    for j in range(i):
        chr1 = SDOC_decrease_TADs[i].split("\t")[0]
        chr2 = SDOC_decrease_TADs[j].split("\t")[0]
        if chr1 != chr2:
            continue

        TAD_pair = SDOC_decrease_TADs[i]+'\t'+SDOC_decrease_TADs[j]
        if TAD_pair not in TAD_pair_2_contact:
            continue

        TAD_pair_contact_8stages.append(TAD_pair_2_contact[TAD_pair])
        TAD_pairs_SDOC_decreasing.append(TAD_pair)


for i in range(len(SDOC_increase_TADs)):
    for j in range(i):
        chr1 = SDOC_increase_TADs[i].split("\t")[0]
        chr2 = SDOC_increase_TADs[j].split("\t")[0]
        if chr1 != chr2:
            continue

        TAD_pair = SDOC_increase_TADs[i]+'\t'+SDOC_increase_TADs[j]
        if TAD_pair not in TAD_pair_2_contact:
            continue

        TAD_pairs_SDOC_increasing.append(TAD_pair)



#3. randomly select TADs:
a = [i for i in range(len(TADs))]
b = [i for i in range(len(TADs))]
#remove SDOC decreasing TADs from the sampling pool:
SDOC_decreasing_TAD_index = df[df["TAD"].isin(SDOC_decrease_TADs)].index
SDOC_increasing_TAD_index = df[df["TAD"].isin(SDOC_increase_TADs)].index
a = np.setdiff1d(a,SDOC_decreasing_TAD_index)
b = np.setdiff1d(b,SDOC_increasing_TAD_index)
print('print#3.1',len(a),len(SDOC_decrease_TADs))
print('print#3.2',len(b),len(SDOC_increase_TADs))

random_index = np.random.choice(a, len(SDOC_decrease_TADs), replace=False)
random_index_SDOC_increasing_TAD = np.random.choice(b, len(SDOC_increase_TADs), replace=False)

random_TADs = np.array(TADs)[sorted(random_index)].tolist()  #sort index, or bugged(only half of TAD pair can be collected in TAD_pairs_random)
random_TADs_except_SDOC_increasing_TADS = np.array(TADs)[sorted(random_index_SDOC_increasing_TAD)].tolist()  #sort index, or bugged(only half of TAD pair can be collected in TAD_pairs_random)
sort_list(random_TADs)
sort_list(random_TADs_except_SDOC_increasing_TADS)

random_TAD_pair_contact_8stages = []
TAD_pairs_random = []
for i in range(len(random_TADs)):
    for j in range(i): 
        chr1 = random_TADs[i].split("\t")[0]
        chr2 = random_TADs[j].split("\t")[0]
        if chr1 != chr2:
            continue

        TAD_pair = random_TADs[i]+'\t'+random_TADs[j]
        if TAD_pair not in TAD_pair_2_contact:
            continue

        random_TAD_pair_contact_8stages.append(TAD_pair_2_contact[TAD_pair])
        TAD_pairs_random.append(TAD_pair)

TAD_pairs = TAD_pairs_SDOC_decreasing + TAD_pairs_random

random_TAD_pair_contact_8stages_SDOC_increasing_TADs = []
TAD_pairs_random_SDOC_increasing_TADs = []
for i in range(len(random_TADs_except_SDOC_increasing_TADS)):
    for j in range(i): 
        chr1 = random_TADs_except_SDOC_increasing_TADS[i].split("\t")[0]
        chr2 = random_TADs_except_SDOC_increasing_TADS[j].split("\t")[0]
        if chr1 != chr2:
            continue

        TAD_pair = random_TADs_except_SDOC_increasing_TADS[i]+'\t'+random_TADs_except_SDOC_increasing_TADS[j]
        if TAD_pair not in TAD_pair_2_contact:
            continue

        random_TAD_pair_contact_8stages_SDOC_increasing_TADs.append(TAD_pair_2_contact[TAD_pair])
        TAD_pairs_random_SDOC_increasing_TADs.append(TAD_pair)

TAD_pairs_2 = TAD_pairs_SDOC_increasing + TAD_pairs_random_SDOC_increasing_TADs


#### investigate compartment of low, mid and high TADs:
# add B compartment fraction of each TAD in 8 stages to the main dataframe
f = open("./pre2-out/result_TAD_8stage_B_proportion.tsv")  
B_comp_proportion_data_8stages = [[] for i in range(8)]         
lines=f.readlines() 
nrow = len(lines)                    
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    L[3:] = list(map(float,L[3:]))
    for j in range(8):
        B_comp_proportion_data_8stages[j].append(L[3+j])
f.close()

for i in range(8):
    stage_this = celltypes[i]
    col_name = "B_compartment_fraction_"+stage_this
    df[col_name] = B_comp_proportion_data_8stages[i]



#draw SDOC deceasing TAD expression and SDOC increasing TAD expression
exp_inc = df[df["TAD"].isin(SDOC_increase_TADs)].loc[:,"mean_expression_0":"mean_expression_7"]
exp_dec = df[df["TAD"].isin(SDOC_decrease_TADs)].loc[:,"mean_expression_0":"mean_expression_7"]

draw_boxplot_8stages(exp_inc.dropna(how = 'any').to_numpy(), "./1-figures/exp_SDOC_increase.png",[-0.4,8],[[1,0.9,0.6] for i in range(8)])
draw_boxplot_8stages(exp_dec.dropna(how = 'any').to_numpy(), "./1-figures/exp_SDOC_decrease.png",[-0.4,8],[[0.5,0.7,0.9] for i in range(8)])
draw_boxplot_8stages(exp_inc.dropna(how = 'any').to_numpy(), "./1-figures/exp_SDOC_increase.pdf",[-0.4,8],[[1,0.9,0.6] for i in range(8)])
draw_boxplot_8stages(exp_dec.dropna(how = 'any').to_numpy(), "./1-figures/exp_SDOC_decrease.pdf",[-0.4,8],[[0.5,0.7,0.9] for i in range(8)])





############## try using a simpler metric to replace TCI (better do it at last to avoid the tremendous time cost ) ########################
# 1. get data matrix
# 2. do clustering
# 3. draw heatmap and bar, do enrichment analysis


# 1.get data matrix for each chromosome:
# read TAD-TAD 8 stage contact file and get the slope 
print("123")
TAD_pair_2_contact_slope = {}

if os.path.isfile("./1-out/TAD_pair_2_slope"):
    f = open("./1-out/TAD_pair_2_slope")           
    lines=f.readlines() 
    nrow = len(lines)                   
    for i in range(len(lines)):
        L = lines[i].strip().split(':')
        TAD_pair = L[0]
        slope = float(L[1])
        TAD_pair_2_contact_slope[TAD_pair] = slope
    f.close()
else:
    out = []
    for TAD_pair in TAD_pair_2_contact:
         x = np.asarray([[i/7] for i in range(8)])
         y = np.asarray([[TAD_pair_2_contact[TAD_pair][i]] for i in range(8)])
         reg = LinearRegression().fit(x, y)
         slope = reg.coef_[0][0]
         TAD_pair_2_contact_slope[TAD_pair] = slope
         out.append(TAD_pair+':'+str(slope))
    savetxt("./1-out/TAD_pair_2_slope", out)

### 12.5 try change the defination of TCI to "mean slope of TAD-TAD contact"
SDOC_decrease_TADs_2_TCI = {}  #TCI: mean slope
SDOC_decrease_TADs_2_count_pairs = {}
for TAD_pair in TAD_pairs:
    L = TAD_pair.split("\t")
    TAD1 = "\t".join(L[0:3])
    TAD2 = "\t".join(L[3:6])
    if TAD1 not in SDOC_decrease_TADs or TAD2 not in SDOC_decrease_TADs:
        continue
    if TAD1 not in SDOC_decrease_TADs_2_TCI:
        SDOC_decrease_TADs_2_TCI[TAD1] = TAD_pair_2_contact_slope[TAD_pair]
        SDOC_decrease_TADs_2_count_pairs[TAD1] = 1
    else:
        SDOC_decrease_TADs_2_TCI[TAD1] += TAD_pair_2_contact_slope[TAD_pair]
        SDOC_decrease_TADs_2_count_pairs[TAD1] += 1
    if TAD2 not in SDOC_decrease_TADs_2_TCI:
        SDOC_decrease_TADs_2_TCI[TAD2] = TAD_pair_2_contact_slope[TAD_pair]
        SDOC_decrease_TADs_2_count_pairs[TAD2] = 1
    else:
        SDOC_decrease_TADs_2_TCI[TAD2] += TAD_pair_2_contact_slope[TAD_pair]
        SDOC_decrease_TADs_2_count_pairs[TAD2] += 1

count_SDOC_decrease_TAD = 0
for TAD in SDOC_decrease_TADs_2_TCI:
    SDOC_decrease_TADs_2_TCI[TAD] /= SDOC_decrease_TADs_2_count_pairs[TAD]  # got TCI for each SDOC decreasing TAD
    count_SDOC_decrease_TAD += 1
print("count SDOC decrease TADs:", count_SDOC_decrease_TAD)


random_TADs_2_TCI = {}  #TCI: mean slope
random_TADs_2_count_pairs = {}
for TAD_pair in TAD_pairs:
    L = TAD_pair.split("\t")
    TAD1 = "\t".join(L[0:3])
    TAD2 = "\t".join(L[3:6])
    if TAD1 not in random_TADs or TAD2 not in random_TADs:
        continue
    if TAD1 not in random_TADs_2_TCI:
        random_TADs_2_TCI[TAD1] = TAD_pair_2_contact_slope[TAD_pair]
        random_TADs_2_count_pairs[TAD1] = 1
    else:
        random_TADs_2_TCI[TAD1] += TAD_pair_2_contact_slope[TAD_pair]
        random_TADs_2_count_pairs[TAD1] += 1
    if TAD2 not in random_TADs_2_TCI:
        random_TADs_2_TCI[TAD2] = TAD_pair_2_contact_slope[TAD_pair]
        random_TADs_2_count_pairs[TAD2] = 1
    else:
        random_TADs_2_TCI[TAD2] += TAD_pair_2_contact_slope[TAD_pair]
        random_TADs_2_count_pairs[TAD2] += 1

for TAD in random_TADs_2_TCI:
    random_TADs_2_TCI[TAD] /= random_TADs_2_count_pairs[TAD]  # got TCI for each SDOC decreasing TAD


# plot TCI:
TCIs_SDOC_decreasing = sorted(list(SDOC_decrease_TADs_2_TCI.values()))
TCIs_random = sorted(list(random_TADs_2_TCI.values()))

mkdir("./1-figures/TCI(mean_slope)_plot")

plt.figure(figsize=(3, 2.5), dpi=600)
plt.rc('font',family='Arial')
plt.rc('font',size = 8)
plt.plot(TCIs_random, color = [0.8,0.8,0.8])
plt.plot(TCIs_SDOC_decreasing, color = [0.5,0.7,0.9])
plt.ylabel("TAD clustering index")
plt.tight_layout()
plt.savefig("./1-figures/TCI(mean_slope)_plot/TCI_plot.png")
plt.savefig("./1-figures/TCI(mean_slope)_plot/TCI_plot.pdf")
plt.close()
stat, pvalue = stats.ks_2samp(TCIs_SDOC_decreasing, TCIs_random)
print("K-S test results for TCI:", stat, pvalue)
stat, pvalue = stats.ttest_ind(TCIs_SDOC_decreasing, TCIs_random)
print("t test results for TCI:", stat, pvalue)
#also do hist plot:
plt.figure(figsize=(3, 2.5), dpi=600)
plt.rc('font',family='Arial')
plt.rc('font',size = 8)
plt.hist(TCIs_SDOC_decreasing, bins = int((max(TCIs_SDOC_decreasing)-min(TCIs_SDOC_decreasing))*10), density=False, color = [0.5,0.7,0.9])
plt.hist(TCIs_random, bins = int((max(TCIs_random)-min(TCIs_random))*10), density=False, color = [0.6,0.6,0.6], alpha = 0.6)
plt.plot([np.median(TCIs_SDOC_decreasing),np.median(TCIs_SDOC_decreasing)],[0,55], "--", linewidth = 1, color = [0,0.4,0.8])
plt.plot([np.median(TCIs_random),np.median(TCIs_random)],[0,55], "--", linewidth = 1, color = "black")
plt.text(1,49,"p = 1.5e-45")
plt.xlabel("TAD clustering index")
plt.ylabel("Frequency")
plt.tight_layout()
plt.ylim(0,55)
plt.xlim(-2,2)
plt.savefig("./1-figures/TCI(mean_slope)_plot/TCI_plot.png")
plt.savefig("./1-figures/TCI(mean_slope)_plot/TCI_plot.pdf")
plt.close()




#####  3.2.1 boxplot of SDOC decreasing TADs (SDOC, expression of 8 stages for 3 TCi group)
high_proportion_TADs = []  #>0.25
mid_proportion_TADs = [] #0-0.25
low_proportion_TADs = [] #<0

for TAD in SDOC_decrease_TADs_2_TCI:
    if SDOC_decrease_TADs_2_TCI[TAD] <= 0:
        low_proportion_TADs.append(TAD)
    elif SDOC_decrease_TADs_2_TCI[TAD] < 0.25:
        mid_proportion_TADs.append(TAD)
    elif SDOC_decrease_TADs_2_TCI[TAD] >= 0.25:
        high_proportion_TADs.append(TAD)


#@@expression_low_proportion_TADs = df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any')
#@@expression_mid_proportion_TADs = df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any')
#@@expression_high_proportion_TADs = df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any')
#@@#print results
#@@print("expression SDOC decreasing!!!!!!")
#@@print(np.mean(expression_low_proportion_TADs))
#@@print(np.mean(expression_mid_proportion_TADs))
#@@print(np.mean(expression_high_proportion_TADs))
#@@print("expression done")
#@@#plot results
#@@mkdir("./1-figures/3.2-boxplot_SDOC_decrease_TAD")
#@@
#@@draw_boxplot_8stages(expression_low_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_SDOC_decrease_TAD/expression_low_proportion_TADs.png", [-0.3,6.5])
#@@draw_boxplot_8stages(expression_mid_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_SDOC_decrease_TAD/expression_mid_proportion_TADs.png", [-0.3,6.5])
#@@draw_boxplot_8stages(expression_high_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_SDOC_decrease_TAD/expression_high_proportion_TADs.png", [-0.3,6.5])
#@@boxplot_1_dataset_3 = expression_low_proportion_TADs.to_numpy()
#@@boxplot_1_dataset_4 = expression_mid_proportion_TADs.to_numpy()
#@@boxplot_1_dataset_5 = expression_high_proportion_TADs.to_numpy()
#@@
#@@#get mean SDOC of 3 groups of TADs:
#@@SDOC_low_proportion_TADs = df[df["TAD"].isin(low_proportion_TADs)].loc[:,"SDOC_0":"SDOC_7"].dropna(how = 'any')
#@@SDOC_mid_proportion_TADs = df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"SDOC_0":"SDOC_7"].dropna(how = 'any')
#@@SDOC_high_proportion_TADs = df[df["TAD"].isin(high_proportion_TADs)].loc[:,"SDOC_0":"SDOC_7"].dropna(how = 'any')
#@@
#@@#get compartment B fraction of each TAD
#@@low_TCI_B_comp_fractions_SDOC_decreasing = df.loc[df["TAD"].isin(low_proportion_TADs), "B_compartment_fraction_HSC":"B_compartment_fraction_DP"].to_numpy()
#@@mid_TCI_B_comp_fractions_SDOC_decreasing = df.loc[df["TAD"].isin(mid_proportion_TADs), "B_compartment_fraction_HSC":"B_compartment_fraction_DP"].to_numpy()
#@@high_TCI_B_comp_fractions_SDOC_decreasing = df.loc[df["TAD"].isin(high_proportion_TADs), "B_compartment_fraction_HSC":"B_compartment_fraction_DP"].to_numpy()
B_comp_fractions_all_SDOC_decreasing = df.loc[df["TAD"].isin(SDOC_decrease_TADs), "B_compartment_fraction_HSC":"B_compartment_fraction_DP"].to_numpy()
#@@
#@@#print results
#@@
#@@#plot results
#@@draw_boxplot_8stages(SDOC_low_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_SDOC_decrease_TAD/SDOC_low_proportion_TADs.png",[-2.3,3])
#@@draw_boxplot_8stages(SDOC_mid_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_SDOC_decrease_TAD/SDOC_mid_proportion_TADs.png",[-2.3,3])
#@@draw_boxplot_8stages(SDOC_high_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_SDOC_decrease_TAD/SDOC_high_proportion_TADs.png",[-2.3,3])

# add mean exp and exp fold change boxplot for 6 conditions:

#@@exp_fc_low_proportion_TADs_SDOC_decreased = np.log(np.mean(df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_5":"mean_expression_7"].dropna(how = 'any').to_numpy()+1, axis = 1) / np.mean(df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_4"].dropna(how = 'any').to_numpy()+1, axis = 1)) / np.log(2)
#@@exp_fc_mid_proportion_TADs_SDOC_decreased = np.log(np.mean(df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_5":"mean_expression_7"].dropna(how = 'any').to_numpy()+1, axis = 1) / np.mean(df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_4"].dropna(how = 'any').to_numpy()+1, axis = 1)) / np.log(2)
#@@exp_fc_high_proportion_TADs_SDOC_decreased = np.log(np.mean(df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_5":"mean_expression_7"].dropna(how = 'any').to_numpy()+1, axis = 1) / np.mean(df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_4"].dropna(how = 'any').to_numpy()+1, axis = 1)) / np.log(2)

exp_fc_low_proportion_TADs_SDOC_decreased = np.log((df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1)) / np.log(2)
exp_fc_mid_proportion_TADs_SDOC_decreased = np.log((df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1))/ np.log(2)
exp_fc_high_proportion_TADs_SDOC_decreased = np.log((df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1)) / np.log(2)

exp_fc_low_proportion_TADs_SDOC_decreased = exp_fc_low_proportion_TADs_SDOC_decreased[exp_fc_low_proportion_TADs_SDOC_decreased != 0]
exp_fc_mid_proportion_TADs_SDOC_decreased = exp_fc_mid_proportion_TADs_SDOC_decreased[exp_fc_mid_proportion_TADs_SDOC_decreased != 0]
exp_fc_high_proportion_TADs_SDOC_decreased = exp_fc_high_proportion_TADs_SDOC_decreased[exp_fc_high_proportion_TADs_SDOC_decreased != 0]

#@@mean_exp_low_proportion_TADs_SDOC_decreased = np.mean(df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any'), axis = 1)
#@@mean_exp_mid_proportion_TADs_SDOC_decreased = np.mean(df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any'), axis = 1)
#@@mean_exp_high_proportion_TADs_SDOC_decreased = np.mean(df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any'), axis = 1)

########## 3.2.2 boxplot of random TADs (SDOC, expression of 8 stages for 3 TCi group)
high_proportion_TADs = []  #>0.25
mid_proportion_TADs = [] #0-0.25
low_proportion_TADs = [] #<0

for TAD in random_TADs_2_TCI:
    if random_TADs_2_TCI[TAD] <= 0:
        low_proportion_TADs.append(TAD)
    elif random_TADs_2_TCI[TAD] < 0.25:
        mid_proportion_TADs.append(TAD)
    elif random_TADs_2_TCI[TAD] >= 0.25:
        high_proportion_TADs.append(TAD)

#@@
#@@expression_low_proportion_TADs = df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any')
#@@expression_mid_proportion_TADs = df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any')
#@@expression_high_proportion_TADs = df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any')
#@@#print results
#@@print("expression SDOC decreasing!!!!!!")
#@@print(np.mean(expression_low_proportion_TADs))
#@@print(np.mean(expression_mid_proportion_TADs))
#@@print(np.mean(expression_high_proportion_TADs))
#@@print("expression done")
#@@#plot results
#@@
#@@mkdir("./1-figures/3.2-boxplot_random_TAD")
#@@draw_boxplot_8stages(expression_low_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_random_TAD/expression_low_proportion_TADs.png", [-0.3,6.5])
#@@draw_boxplot_8stages(expression_mid_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_random_TAD/expression_mid_proportion_TADs.png", [-0.3,6.5])
#@@draw_boxplot_8stages(expression_high_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_random_TAD/expression_high_proportion_TADs.png", [-0.3,6.5])
#@@boxplot_1_dataset_0 = expression_low_proportion_TADs.to_numpy()
#@@boxplot_1_dataset_1 = expression_mid_proportion_TADs.to_numpy()
#@@boxplot_1_dataset_2 = expression_high_proportion_TADs.to_numpy()
#@@
#@@
#@@#get mean SDOC of 3 groups of TADs:
#@@SDOC_low_proportion_TADs = df[df["TAD"].isin(low_proportion_TADs)].loc[:,"SDOC_0":"SDOC_7"].dropna(how = 'any')
#@@SDOC_mid_proportion_TADs = df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"SDOC_0":"SDOC_7"].dropna(how = 'any')
#@@SDOC_high_proportion_TADs = df[df["TAD"].isin(high_proportion_TADs)].loc[:,"SDOC_0":"SDOC_7"].dropna(how = 'any')
#@@
#@@#get compartment B fraction of each TAD
#@@low_TCI_B_comp_fractions_Random = df.loc[df["TAD"].isin(low_proportion_TADs), "B_compartment_fraction_HSC":"B_compartment_fraction_DP"].to_numpy()
#@@mid_TCI_B_comp_fractions_Random = df.loc[df["TAD"].isin(mid_proportion_TADs), "B_compartment_fraction_HSC":"B_compartment_fraction_DP"].to_numpy()
#@@high_TCI_B_comp_fractions_Random = df.loc[df["TAD"].isin(high_proportion_TADs), "B_compartment_fraction_HSC":"B_compartment_fraction_DP"].to_numpy()
B_comp_fractions_all_Random = df.loc[df["TAD"].isin(random_TADs), "B_compartment_fraction_HSC":"B_compartment_fraction_DP"].to_numpy()
#@@
#@@#print results
#@@
#@@#plot results
#@@draw_boxplot_8stages(SDOC_low_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_random_TAD/SDOC_low_proportion_TADs.png",[-2.3,3])
#@@draw_boxplot_8stages(SDOC_mid_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_random_TAD/SDOC_mid_proportion_TADs.png",[-2.3,3])
#@@draw_boxplot_8stages(SDOC_high_proportion_TADs.to_numpy(), "./1-figures/3.2-boxplot_random_TAD/SDOC_high_proportion_TADs.png",[-2.3,3])


#@@
#@@# do combined boxplot
#@@
#@@color_bar0 = [get_light_color([190,190,190]) for i in range(8)]
#@@color_bar1 = [get_light_color([140,140,140]) for i in range(8)]
#@@color_bar2 = [get_light_color([80,80,80]) for i in range(8)]
#@@color_bar3 = [get_light_color([0,238,238]) for i in range(8)]
#@@color_bar4 = [get_light_color([0,191,255]) for i in range(8)]
#@@color_bar5 = [get_light_color([65,105,225]) for i in range(8)]
#@@color0 = np.array([190,190,190])/255
#@@color1 = np.array([140,140,140])/255
#@@color2 = np.array([80,80,80])/255
#@@color3 = np.array([0,238,238])/255
#@@color4 = np.array([0,191,255])/255
#@@color5 = np.array([65,105,225])/255
#@@
#@@plt.figure(figsize=(5, 4), dpi=600)
#@@plt.rc('font',family='Arial')
#@@plt.rc('font',size = 8)
#@@
#@@
#@@linewidth_set = 1.4
#@@boxprops = {'linewidth':linewidth_set , 'color': color0}
#@@whiskerprops = {'linestyle': '-', 'linewidth':0.6, 'color': color0}
#@@medianprops = {'linewidth':0.5 , 'color': color0, 'drawstyle' : 'steps'}
#@@bp = plt.boxplot(boxplot_1_dataset_0, positions = [i*7-1 for i in range(8)], showfliers=False, meanline =False, patch_artist = True,  boxprops = boxprops, whiskerprops = whiskerprops, medianprops = medianprops, showcaps = False, widths = 0.75,whis = [10, 90])
#@@for patch, color in zip(bp['boxes'], color_bar0):
#@@    patch.set_facecolor(color)
#@@
#@@linewidth_set = 1.4
#@@boxprops = {'linewidth':linewidth_set , 'color': color1}
#@@whiskerprops = {'linestyle': '-', 'linewidth':0.6, 'color': color1}
#@@medianprops = {'linewidth':0.5 , 'color': color1, 'drawstyle' : 'steps'}
#@@bp = plt.boxplot(boxplot_1_dataset_1, positions = [i*7 for i in range(8)], showfliers=False, meanline =False, patch_artist = True,  boxprops = boxprops, whiskerprops = whiskerprops, medianprops = medianprops, showcaps = False, widths = 0.75,whis = [10, 90])
#@@for patch, color in zip(bp['boxes'], color_bar1):
#@@    patch.set_facecolor(color)
#@@
#@@boxprops = {'linewidth':linewidth_set , 'color': color2}
#@@whiskerprops = {'linestyle': '-', 'linewidth':0.6, 'color': color2}
#@@medianprops = {'linewidth':0.5 , 'color': color2, 'drawstyle' : 'steps'}
#@@bp1 = plt.boxplot(boxplot_1_dataset_2, positions = [i*7+1 for i in range(8)], showfliers=False, meanline =False, patch_artist = True,  boxprops = boxprops, whiskerprops = whiskerprops, medianprops = medianprops,showcaps = False, widths = 0.75,whis = [10, 90])
#@@for patch, color in zip(bp1['boxes'], color_bar2):
#@@    patch.set_facecolor(color)
#@@
#@@whiskerprops = {'linestyle': '-', 'linewidth':0.6, 'color': color3}
#@@boxprops = {'linewidth':linewidth_set , 'color': color3}
#@@medianprops = {'linewidth':0.5 , 'color': color3, 'drawstyle' : 'steps'}
#@@bp2 = plt.boxplot(boxplot_1_dataset_3, positions = [i*7+2 for i in range(8)], showfliers=False, meanline =False, patch_artist = True,  boxprops = boxprops, whiskerprops = whiskerprops, medianprops = medianprops,showcaps = False, widths = 0.75,whis = [10, 90])
#@@for patch, color in zip(bp2['boxes'], color_bar3):
#@@    patch.set_facecolor(color)
#@@
#@@boxprops = {'linewidth':linewidth_set , 'color': color4}
#@@whiskerprops = {'linestyle': '-', 'linewidth':0.6, 'color': color4}
#@@medianprops = {'linewidth':0.5 , 'color': color4, 'drawstyle' : 'steps'}
#@@bp3 = plt.boxplot(boxplot_1_dataset_4, positions = [i*7+3 for i in range(8)], showfliers=False, meanline =False, patch_artist = True, boxprops = boxprops, whiskerprops = whiskerprops, medianprops = medianprops,showcaps = False, widths = 0.75,whis = [10, 90])
#@@for patch, color in zip(bp3['boxes'], color_bar4):
#@@    patch.set_facecolor(color)
#@@
#@@boxprops = {'linewidth':linewidth_set , 'color': color5}
#@@whiskerprops = {'linestyle': '-', 'linewidth':0.6, 'color': color5}
#@@medianprops = {'linewidth':0.5 , 'color': color5, 'drawstyle' : 'steps'}
#@@bp4 = plt.boxplot(boxplot_1_dataset_5, positions = [i*7+4 for i in range(8)], showfliers=False, meanline =False, patch_artist = True, boxprops = boxprops, whiskerprops = whiskerprops, medianprops = medianprops,showcaps = False, widths = 0.75,whis = [10, 90])
#@@for patch, color in zip(bp4['boxes'], color_bar5):
#@@    patch.set_facecolor(color)
#@@
#@@#plt.xticks([i+1 for i in range(8)],["HSC","MPP","CLP","ETP","DN2","DN3","DN4","DP"], rotation = 0)
#@@ 
#@@for i in range(8):
#@@    plt.scatter(i*7-1,np.mean(boxplot_1_dataset_0[:,i]),s = 1, linewidth = 0, color = color0, zorder = 10)
#@@    plt.scatter(i*7,np.mean(boxplot_1_dataset_1[:,i]),s = 1, linewidth = 0, color = color1, zorder = 10)
#@@    plt.scatter(i*7+1,np.mean(boxplot_1_dataset_2[:,i]),s = 1, linewidth = 0, color = color2, zorder = 10)
#@@    plt.scatter(i*7+2,np.mean(boxplot_1_dataset_3[:,i]),s = 1, linewidth = 0, color = color3, zorder = 10)
#@@    plt.scatter(i*7+3,np.mean(boxplot_1_dataset_4[:,i]),s = 1, linewidth = 0, color = color4, zorder = 10)
#@@    plt.scatter(i*7+4,np.mean(boxplot_1_dataset_5[:,i]),s = 1, linewidth = 0, color = color5, zorder = 10)
#@@
#@@for i in range(7):
#@@    plt.plot([i*7-1,i*7+6],[np.mean(boxplot_1_dataset_0[:,i]),np.mean(boxplot_1_dataset_0[:,i+1])], linewidth = 1.5, color = color0, zorder = 10)   
#@@    plt.plot([i*7,i*7+7],[np.mean(boxplot_1_dataset_1[:,i]),np.mean(boxplot_1_dataset_1[:,i+1])], linewidth = 1.5, color = color1, zorder = 10)
#@@    plt.plot([i*7+1,i*7+8],[np.mean(boxplot_1_dataset_2[:,i]),np.mean(boxplot_1_dataset_2[:,i+1])], linewidth = 1.5, color = color2, zorder = 10)
#@@    plt.plot([i*7+2,i*7+9],[np.mean(boxplot_1_dataset_3[:,i]),np.mean(boxplot_1_dataset_3[:,i+1])], linewidth = 1.5, color = color3, zorder = 10)
#@@    plt.plot([i*7+3,i*7+10],[np.mean(boxplot_1_dataset_4[:,i]),np.mean(boxplot_1_dataset_4[:,i+1])], linewidth = 1.5, color = color4, zorder = 10)
#@@    plt.plot([i*7+4,i*7+11],[np.mean(boxplot_1_dataset_5[:,i]),np.mean(boxplot_1_dataset_5[:,i+1])], linewidth = 1.5, color = color5, zorder = 10)
#@@plt.ylim(0,7)
#@@plt.xlim(-2,54)
#@@mkdir("./1-figures/3.2-combined_boxplot")
#@@plt.savefig("./1-figures/3.2-combined_boxplot/expression_low_mid_high.png")
#@@plt.savefig("./1-figures/3.2-combined_boxplot/expression_low_mid_high.pdf")
#@@plt.close()

# add mean exp and exp fold change boxplot for 6 conditions:
#$$exp_fc_low_proportion_TADs = np.log(np.mean(df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_5":"mean_expression_7"].dropna(how = 'any').to_numpy()+1, axis = 1) / np.mean(df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_4"].dropna(how = 'any').to_numpy()+1, axis = 1)) / np.log(2)
#$$exp_fc_mid_proportion_TADs = np.log(np.mean(df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_5":"mean_expression_7"].dropna(how = 'any').to_numpy()+1, axis = 1) / np.mean(df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_4"].dropna(how = 'any').to_numpy()+1, axis = 1)) / np.log(2)
#$$exp_fc_high_proportion_TADs = np.log(np.mean(df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_5":"mean_expression_7"].dropna(how = 'any').to_numpy()+1, axis = 1) / np.mean(df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_4"].dropna(how = 'any').to_numpy()+1, axis = 1)) / np.log(2)

exp_fc_low_proportion_TADs = np.log((df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1)) / np.log(2)
exp_fc_mid_proportion_TADs = np.log((df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1))/ np.log(2)
exp_fc_high_proportion_TADs = np.log((df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1)) / np.log(2)

exp_fc_low_proportion_TADs = exp_fc_low_proportion_TADs[exp_fc_low_proportion_TADs != 0]
exp_fc_mid_proportion_TADs = exp_fc_mid_proportion_TADs[exp_fc_mid_proportion_TADs != 0]
exp_fc_high_proportion_TADs = exp_fc_high_proportion_TADs[exp_fc_high_proportion_TADs != 0]


#@@mean_exp_low_proportion_TADs = np.mean(df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any'), axis = 1)
#@@mean_exp_mid_proportion_TADs = np.mean(df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any'), axis = 1)
#@@mean_exp_high_proportion_TADs = np.mean(df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_0":"mean_expression_7"].dropna(how = 'any'), axis = 1)
#do combined boxplot (exp fold change)
do_boxplot_exp_fc(exp_fc_low_proportion_TADs,
exp_fc_mid_proportion_TADs,
exp_fc_high_proportion_TADs,
exp_fc_low_proportion_TADs_SDOC_decreased,
exp_fc_mid_proportion_TADs_SDOC_decreased,
exp_fc_high_proportion_TADs_SDOC_decreased, "boxplot_exp_fc")

#@@do_boxplot_exp_fc(mean_exp_low_proportion_TADs,
#@@mean_exp_mid_proportion_TADs,
#@@mean_exp_high_proportion_TADs,
#@@mean_exp_low_proportion_TADs_SDOC_decreased,
#@@mean_exp_mid_proportion_TADs_SDOC_decreased,
#@@mean_exp_high_proportion_TADs_SDOC_decreased, "boxplot_mean_exp")
#@@print("print #funny\n",np.mean(df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_5":"mean_expression_7"].dropna(how = 'any').to_numpy()+1, axis = 1),"\n",exp_fc_mid_proportion_TADs_SDOC_decreased)







########################################################################################################################
#----------------------------------------------------------------------------------------------------------------------#
#--------------------------------  BIB response - add SDOC increasing TAD TCI analysis  -------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
########################################################################################################################

print("TAD increase SDOC count:", len(list(set(SDOC_increase_TADs))))

SDOC_increase_TADs_2_TCI = {}  #TCI: mean slope
SDOC_increase_TADs_2_count_pairs = {}
for TAD_pair in TAD_pairs_2:
    L = TAD_pair.split("\t")
    TAD1 = "\t".join(L[0:3])
    TAD2 = "\t".join(L[3:6])
    if TAD1 not in SDOC_increase_TADs or TAD2 not in SDOC_increase_TADs:
        continue
    if TAD1 not in SDOC_increase_TADs_2_TCI:
        SDOC_increase_TADs_2_TCI[TAD1] = TAD_pair_2_contact_slope[TAD_pair]
        SDOC_increase_TADs_2_count_pairs[TAD1] = 1
    else:
        SDOC_increase_TADs_2_TCI[TAD1] += TAD_pair_2_contact_slope[TAD_pair]
        SDOC_increase_TADs_2_count_pairs[TAD1] += 1
    if TAD2 not in SDOC_increase_TADs_2_TCI:
        SDOC_increase_TADs_2_TCI[TAD2] = TAD_pair_2_contact_slope[TAD_pair]
        SDOC_increase_TADs_2_count_pairs[TAD2] = 1
    else:
        SDOC_increase_TADs_2_TCI[TAD2] += TAD_pair_2_contact_slope[TAD_pair]
        SDOC_increase_TADs_2_count_pairs[TAD2] += 1

count_TAD_increase_SDOC = 0
for TAD in SDOC_increase_TADs_2_TCI:
    count_TAD_increase_SDOC+=1
    SDOC_increase_TADs_2_TCI[TAD] /= SDOC_increase_TADs_2_count_pairs[TAD]  # got TCI for each SDOC increasing TAD
print("TAD increase SDOC count:", count_TAD_increase_SDOC)


random_TADs_2_TCI = {}  #TCI: mean slope
random_TADs_2_count_pairs = {}
for TAD_pair in TAD_pairs_2:
    L = TAD_pair.split("\t")
    TAD1 = "\t".join(L[0:3])
    TAD2 = "\t".join(L[3:6])
    if TAD1 not in random_TADs_except_SDOC_increasing_TADS or TAD2 not in random_TADs_except_SDOC_increasing_TADS:
        continue
    if TAD1 not in random_TADs_2_TCI:
        random_TADs_2_TCI[TAD1] = TAD_pair_2_contact_slope[TAD_pair]
        random_TADs_2_count_pairs[TAD1] = 1
    else:
        random_TADs_2_TCI[TAD1] += TAD_pair_2_contact_slope[TAD_pair]
        random_TADs_2_count_pairs[TAD1] += 1
    if TAD2 not in random_TADs_2_TCI:
        random_TADs_2_TCI[TAD2] = TAD_pair_2_contact_slope[TAD_pair]
        random_TADs_2_count_pairs[TAD2] = 1
    else:
        random_TADs_2_TCI[TAD2] += TAD_pair_2_contact_slope[TAD_pair]
        random_TADs_2_count_pairs[TAD2] += 1

for TAD in random_TADs_2_TCI:
    random_TADs_2_TCI[TAD] /= random_TADs_2_count_pairs[TAD]  # got TCI for each SDOC increasing TAD
# plot TCI:
TCIs_SDOC_increasing = sorted(list(SDOC_increase_TADs_2_TCI.values()))
TCIs_random = sorted(list(random_TADs_2_TCI.values()))

mkdir("./1-figures/TCI(mean_slope)_plot")

plt.figure(figsize=(3, 2.5), dpi=600)
plt.rc('font',family='Arial')
plt.rc('font',size = 8)
plt.plot(TCIs_random, color = [0.8,0.8,0.8])
plt.plot(TCIs_SDOC_increasing, color = [0.9,0.75,0.45])
plt.ylabel("TAD clustering index")
plt.tight_layout()
plt.savefig("./1-figures/TCI(mean_slope)_plot/TCI_plot(SDOC_increasing).png")
plt.savefig("./1-figures/TCI(mean_slope)_plot/TCI_plot(SDOC_increasing).pdf")
plt.close()
stat, pvalue = stats.ks_2samp(TCIs_SDOC_increasing, TCIs_random)
print("K-S test results for TCI:", stat, pvalue)
stat, pvalue = stats.ttest_ind(TCIs_SDOC_increasing, TCIs_random)
print("t test results for TCI:", stat, pvalue)

#also do hist plot:
plt.figure(figsize=(3, 2.5), dpi=600)
plt.rc('font',family='Arial')
plt.rc('font',size = 8)
plt.hist(TCIs_SDOC_increasing, bins = int((max(TCIs_SDOC_increasing)-min(TCIs_SDOC_increasing))*10), density=False, color = [0.9,0.75,0.45])
plt.hist(TCIs_random, bins = int((max(TCIs_random)-min(TCIs_random))*10), density=False, color = [0.6,0.6,0.6], alpha = 0.6)
plt.plot([np.median(TCIs_SDOC_increasing),np.median(TCIs_SDOC_increasing)],[0,55], "--", linewidth = 1, color = [0.7,0.55,0])
plt.plot([np.median(TCIs_random),np.median(TCIs_random)],[0,55], "--", linewidth = 1, color = "black")
plt.text(1,49,"p = 0.0025")
plt.xlabel("TAD clustering index")
plt.ylabel("Frequency")
plt.tight_layout()
plt.ylim(0,55)
plt.xlim(-2,2)
plt.savefig("./1-figures/TCI(mean_slope)_plot/TCI_hist(SDOC_increasing).png")
plt.savefig("./1-figures/TCI(mean_slope)_plot/TCI_hist(SDOC_increasing).pdf")
plt.close()


#####  3.2.1 boxplot of SDOC increasing TADs (SDOC, expression of 8 stages for 3 TCi group)
high_proportion_TADs = []  #>0.25
mid_proportion_TADs = [] #0-0.25
low_proportion_TADs = [] #<0

for TAD in SDOC_increase_TADs_2_TCI:
    if SDOC_increase_TADs_2_TCI[TAD] <= 0:
        low_proportion_TADs.append(TAD)
    elif SDOC_increase_TADs_2_TCI[TAD] < 0.25:
        mid_proportion_TADs.append(TAD)
    elif SDOC_increase_TADs_2_TCI[TAD] >= 0.25:
        high_proportion_TADs.append(TAD)

# add mean exp and exp fold change boxplot for 6 conditions:
exp_fc_low_proportion_TADs_SDOC_increased = np.log((df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1)) / np.log(2)
exp_fc_mid_proportion_TADs_SDOC_increased = np.log((df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1))/ np.log(2)
exp_fc_high_proportion_TADs_SDOC_increased = np.log((df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1)) / np.log(2)

exp_fc_low_proportion_TADs_SDOC_increased = exp_fc_low_proportion_TADs_SDOC_increased[exp_fc_low_proportion_TADs_SDOC_increased != 0]
exp_fc_mid_proportion_TADs_SDOC_increased = exp_fc_mid_proportion_TADs_SDOC_increased[exp_fc_mid_proportion_TADs_SDOC_increased != 0]
exp_fc_high_proportion_TADs_SDOC_increased = exp_fc_high_proportion_TADs_SDOC_increased[exp_fc_high_proportion_TADs_SDOC_increased != 0]

########## 3.2.2 boxplot of random TADs (SDOC, expression of 8 stages for 3 TCi group)
high_proportion_TADs = []  #>0.25
mid_proportion_TADs = [] #0-0.25
low_proportion_TADs = [] #<0

for TAD in random_TADs_2_TCI:
    if random_TADs_2_TCI[TAD] <= 0:
        low_proportion_TADs.append(TAD)
    elif random_TADs_2_TCI[TAD] < 0.25:
        mid_proportion_TADs.append(TAD)
    elif random_TADs_2_TCI[TAD] >= 0.25:
        high_proportion_TADs.append(TAD)


#print results
exp_fc_low_proportion_TADs = np.log((df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(low_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1)) / np.log(2)
exp_fc_mid_proportion_TADs = np.log((df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(mid_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1))/ np.log(2)
exp_fc_high_proportion_TADs = np.log((df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_7"].dropna(how = 'any').to_numpy()+0.1) / (df[df["TAD"].isin(high_proportion_TADs)].loc[:,"mean_expression_0"].dropna(how = 'any').to_numpy()+0.1)) / np.log(2)

exp_fc_low_proportion_TADs = exp_fc_low_proportion_TADs[exp_fc_low_proportion_TADs != 0]
exp_fc_mid_proportion_TADs = exp_fc_mid_proportion_TADs[exp_fc_mid_proportion_TADs != 0]
exp_fc_high_proportion_TADs = exp_fc_high_proportion_TADs[exp_fc_high_proportion_TADs != 0]

#do combined boxplot (exp fold change)
do_boxplot_exp_fc_SDOC_increase(exp_fc_low_proportion_TADs,
exp_fc_mid_proportion_TADs,
exp_fc_high_proportion_TADs,
exp_fc_low_proportion_TADs_SDOC_increased,
exp_fc_mid_proportion_TADs_SDOC_increased,
exp_fc_high_proportion_TADs_SDOC_increased, "boxplot_exp_fc(SDOC_increasing)")

mkdir("./1-figures/compartment_B_fraction/")
draw_comparative_stacked_bar_8stages(B_comp_fractions_all_SDOC_decreasing, "./1-figures/compartment_B_fraction/Fraction_of_compB_SDOC_decreasing_barplot.png")
draw_comparative_stacked_bar_8stages(B_comp_fractions_all_SDOC_decreasing, "./1-figures/compartment_B_fraction/Fraction_of_compB_SDOC_decreasing_barplot.pdf")