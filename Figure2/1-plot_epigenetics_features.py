#1. test n/V vs log2(n)/V

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats

from sklearn.linear_model import LinearRegression

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


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

def mod_zscore(x):
	median_val = np.median(x)
	std = 0
	for i in range(len(x)):
		std += (x[i]-median_val)**2
	std /= len(x)
	std = np.sqrt(std)
	return np.array([(i-median_val)/std for i in x])


##########

#######  1. build dataframe  


def build_dataframe(celltype):
	TADs = []
	DHS_peaks = []
	SDOCs = []
	log2n_SDOCs = []
	exp = []
	TAD_len = []
	DHS_peak_len = []

	
	
	f = open("./data/DNase_peak_len_10kb/"+celltype)   		
	lines=f.readlines() 
	nrow = len(lines)					
	for i in range(len(lines)):
	    L = lines[i].strip().split('\t')
	    DHS_peak_len.append(float(L[3]))
	f.close()
	
	f = open("./data/"+celltype+"_10kb")   		
	lines=f.readlines() 
	nrow = len(lines)					
	for i in range(len(lines)):
	    L = lines[i].strip().split('\t')
	    TADs.append("\t".join(L[0:3]))
	    SDOCs.append(float(L[3])/float(L[4]))
	    log2n_SDOCs.append(np.log(float(L[3])+1)/np.log(2)/float(L[4]))
	    exp.append(float(L[5]))
	    TAD_len.append((float(L[2]) - float(L[1]))/1000)
	f.close()
	
	mkdir("./1-tmp/")
	savetxt("./1-tmp/"+celltype+"_10kb_SDOCs", SDOCs)
	savetxt("./1-tmp/"+celltype+"_10kb_log2n_SDOCs", log2n_SDOCs)

	os.system("Rscript ./quantile_normalize.r "+ celltype) 

	### get QNed SDOC:
	
	QNed_SDOCs = []
	QNed_log2n_SDOCs = []
	
	f = open("./1-QNed_SDOC/"+celltype+"_10kb_SDOCs")  
	lines=f.readlines() 
	nrow = len(lines)					
	for i in range(len(lines)):
	    SDOC = float(lines[i].strip())
	    QNed_SDOCs.append(SDOC)
	f.close()
	
	f = open("./1-QNed_SDOC/"+celltype+"_10kb_log2n_SDOCs")  
	lines=f.readlines() 
	nrow = len(lines)					
	for i in range(len(lines)):
	    SDOC = float(lines[i].strip())
	    QNed_log2n_SDOCs.append(SDOC)
	f.close()
	
	
	df = pd.DataFrame()
	
	df["TADs"] = TADs
	df["SDOCs"] = SDOCs
	df["log2n_SDOCs"] = log2n_SDOCs
	df["QNed_SDOCs"] = QNed_SDOCs
	df["QNed_log2n_SDOCs"] = QNed_log2n_SDOCs  #not good
	df["log2_mean_expression"] = exp
	df["TAD_length_kb"] = TAD_len
	
	df["mod_Zscore_SDOC"] = mod_zscore(df["SDOCs"])
	df["Zscore_SDOC"] = stats.zscore(df["SDOCs"])

	df["group"] = [0 for i in range(len(df))]
	df.loc[df["QNed_SDOCs"] < -1.5 ,"group"] = 0
	df.loc[np.abs(df["QNed_SDOCs"] + 1) < 0.5 ,"group"] = 1
	df.loc[np.abs(df["QNed_SDOCs"]) < 0.5 ,"group"] = 2
	df.loc[np.abs(df["QNed_SDOCs"] - 1) < 0.5 ,"group"] = 3
	df.loc[df["QNed_SDOCs"] > 1.5,"group"] = 4
	


	print(df.head())

	mkdir("./1-dataframes")
	df.to_csv("./1-dataframes/"+celltype, sep = "\t")
	return df

###### done 


df_GM12878 = build_dataframe("GM12878")
df_IMR90 = build_dataframe("IMR90")
df_HUVEC = build_dataframe("HUVEC")
df_K562 = build_dataframe("K562")

#save dataframe to file:
mkdir("./1-dataframes/")
df_GM12878.to_csv("./1-dataframes/GM12878.tsv", sep = '\t', index = False)


def get_boxplot_data_by_Zscore(group_by_column, data_from_column):
	groups_GM12878 = [
	df_GM12878.loc[df_GM12878[group_by_column] < -1.5 ,data_from_column],
	df_GM12878.loc[np.abs(df_GM12878[group_by_column] + 1) < 0.5 ,data_from_column],
	df_GM12878.loc[np.abs(df_GM12878[group_by_column]) < 0.5 ,data_from_column],
	df_GM12878.loc[np.abs(df_GM12878[group_by_column] - 1) < 0.5 ,data_from_column],
	df_GM12878.loc[df_GM12878[group_by_column] > 1.5,data_from_column]]
	
	groups_IMR90 = [
	df_IMR90.loc[df_IMR90[group_by_column] < -1.5 ,data_from_column],
	df_IMR90.loc[np.abs(df_IMR90[group_by_column] + 1) < 0.5 ,data_from_column],
	df_IMR90.loc[np.abs(df_IMR90[group_by_column]) < 0.5 ,data_from_column],
	df_IMR90.loc[np.abs(df_IMR90[group_by_column] - 1) < 0.5 ,data_from_column],
	df_IMR90.loc[df_IMR90[group_by_column] > 1.5,data_from_column]]
	
	groups_HUVEC = [
	df_HUVEC.loc[df_HUVEC[group_by_column] < -1.5 ,data_from_column],
	df_HUVEC.loc[np.abs(df_HUVEC[group_by_column] + 1) < 0.5 ,data_from_column],
	df_HUVEC.loc[np.abs(df_HUVEC[group_by_column]) < 0.5 ,data_from_column],
	df_HUVEC.loc[np.abs(df_HUVEC[group_by_column] - 1) < 0.5 ,data_from_column],
	df_HUVEC.loc[df_HUVEC[group_by_column] > 1.5,data_from_column]]
	
	groups_K562 = [
	df_K562.loc[df_K562[group_by_column] < -1.5 ,data_from_column],
	df_K562.loc[np.abs(df_K562[group_by_column] + 1) < 0.5 ,data_from_column],
	df_K562.loc[np.abs(df_K562[group_by_column]) < 0.5 ,data_from_column],
	df_K562.loc[np.abs(df_K562[group_by_column] - 1) < 0.5 ,data_from_column],
	df_K562.loc[df_K562[group_by_column] > 1.5,data_from_column]]
	return groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC



def barplot_10_group_4_datasets(data1,data2,data3,data4,file,ylim = [],color_bar = [[0.7,0.7,0.7] for i in range(5)]):

    plt.figure(figsize=(5, 3), dpi=600)
    plt.rc('font',family='Arial')
    plt.rc('font',size = 8)

    sumlen_data1 = np.sum([len(data1[i]) for i in range(5)])
    sumlen_data2 = np.sum([len(data2[i]) for i in range(5)])
    sumlen_data3 = np.sum([len(data3[i]) for i in range(5)])
    sumlen_data4 = np.sum([len(data4[i]) for i in range(5)])

    plt.bar([i*5 for i in range(5)], [len(data1[i])/sumlen_data1 for i in range(5)], width = 0.95, linewidth = 0, color = [0.4, 0.4, 0.4])
    plt.bar([i*5+1 for i in range(5)], [len(data2[i])/sumlen_data2 for i in range(5)], width = 0.95, linewidth = 0, color = [0.55, 0.55, 0.55])
    plt.bar([i*5+2 for i in range(5)], [len(data3[i])/sumlen_data3 for i in range(5)], width = 0.95, linewidth = 0, color = [0.7, 0.7, 0.7])
    plt.bar([i*5+3 for i in range(5)], [len(data4[i])/sumlen_data4 for i in range(5)], width = 0.95, linewidth = 0, color = [0.85, 0.85, 0.85])
    
    plt.xticks([i*5+1.5 for i in range(5)],["<-1.5σ","-1.5σ~-0.5σ","-0.5σ~0.5σ","0.5σ~1.5σ",">1.5σ"], rotation = 45)
    
    if ylim != []:
        plt.ylim(ylim)

    plt.tight_layout()
    plt.savefig(file)
    plt.close()

def boxplot_10_group_4_datasets(data1,data2,data3,data4,file,ylim = [],color_bar = [[0.7,0.7,0.7] for i in range(5)]):

    plt.figure(figsize=(5, 3), dpi=600)
    plt.rc('font',family='Arial')
    plt.rc('font',size = 8)
    
    boxprops = {'linewidth':1.5 , 'color': 'black'}
    whiskerprops = {'linestyle': '--', 'linewidth':1.1, 'color': 'black'}
    medianprops = {'linewidth':1 , 'color': 'black', 'drawstyle' : 'steps'}


    color_bar1 = [[0.4, 0.4, 0.4] for i in range(5)]
    color_bar2 = [[0.55, 0.55, 0.55] for i in range(5)]
    color_bar3 = [[0.7, 0.7, 0.7] for i in range(5)]
    color_bar4 = [[0.85, 0.85, 0.85] for i in range(5)]

    bp1 = plt.boxplot(data1, positions = [i*5 for i in range(5)] ,showfliers=False, meanline =False, patch_artist = True, whiskerprops = whiskerprops, medianprops = medianprops,widths = 0.77,whis = [5, 95])
    for patch, color in zip(bp1['boxes'], color_bar1):
        patch.set_facecolor(color)

    bp2 = plt.boxplot(data2, positions = [i*5+1 for i in range(5)] ,showfliers=False, meanline =False, patch_artist = True, whiskerprops = whiskerprops, medianprops = medianprops,widths = 0.77,whis = [5, 95])
    for patch, color in zip(bp2['boxes'], color_bar2):
        patch.set_facecolor(color)

    bp3 = plt.boxplot(data3, positions = [i*5+2 for i in range(5)] ,showfliers=False, meanline =False, patch_artist = True, whiskerprops = whiskerprops, medianprops = medianprops,widths = 0.77,whis = [5, 95])
    for patch, color in zip(bp3['boxes'], color_bar3):
        patch.set_facecolor(color)

    bp4 = plt.boxplot(data4, positions = [i*5+3 for i in range(5)] ,showfliers=False, meanline =False, patch_artist = True, whiskerprops = whiskerprops, medianprops = medianprops,widths = 0.77,whis = [5, 95])
    for patch, color in zip(bp4['boxes'], color_bar4):
        patch.set_facecolor(color)
    

    plt.scatter([i*5 for i in range(5)], [np.mean(data1[i]) for i in range(5)], s = 5, color = "cyan", zorder = 5)
    plt.scatter([i*5+1 for i in range(5)], [np.mean(data2[i]) for i in range(5)], s = 5, color = "cyan", zorder = 5)
    plt.scatter([i*5+2 for i in range(5)], [np.mean(data3[i]) for i in range(5)], s = 5, color = "cyan", zorder = 5)
    plt.scatter([i*5+3 for i in range(5)], [np.mean(data4[i]) for i in range(5)], s = 5, color = "cyan", zorder = 5)


    plt.xticks([i*5+1.5 for i in range(5)],["<-1.5σ","-1.5σ~-0.5σ","-0.5σ~0.5σ","0.5σ~1.5σ",">1.5σ"], rotation = 45)
    if ylim != []:
        plt.ylim(ylim)

    plt.tight_layout()
    plt.savefig(file)
    plt.close()

mkdir("./1-figures/3-boxplots/")
mkdir("./1-figures/3-barplots")

groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC = get_boxplot_data_by_Zscore("mod_Zscore_SDOC", "SDOCs")
barplot_10_group_4_datasets(groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC, "./1-figures/3-barplots/rawSDOC_quantiles.png")
boxplot_10_group_4_datasets(groups_GM12878,
	groups_K562,
	groups_IMR90,
	groups_HUVEC,
	"./1-figures/3-boxplots/rawSDOC_quantiles.png",
	)


groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC = get_boxplot_data_by_Zscore("mod_Zscore_SDOC", "mod_Zscore_SDOC")
barplot_10_group_4_datasets(groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC, "./1-figures/3-barplots/modified_zscore_quantiles.png", [0,0.8])
boxplot_10_group_4_datasets(groups_GM12878,
	groups_K562,
	groups_IMR90,
	groups_HUVEC,
	"./1-figures/3-boxplots/modified_zscore_quantiles.png",
	)


groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC = get_boxplot_data_by_Zscore("mod_Zscore_SDOC", "log2_mean_expression")
barplot_10_group_4_datasets(groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC, "./1-figures/3-barplots/modified_zscore_quantiles_expression.png", [0,0.8])
boxplot_10_group_4_datasets(groups_GM12878,
	groups_K562,
	groups_IMR90,
	groups_HUVEC,
	"./1-figures/3-boxplots/modified_zscore_quantiles_expression.png",
	)

groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC = get_boxplot_data_by_Zscore("QNed_SDOCs", "QNed_SDOCs")
boxplot_10_group_4_datasets(groups_GM12878,
	groups_K562,
	groups_IMR90,
	groups_HUVEC,
	"./1-figures/3-boxplots/QNed_SDOC_quantiles.png",
	)


for i in range(5):
	print(len(groups_K562[i]))


groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC = get_boxplot_data_by_Zscore("Zscore_SDOC", "log2_mean_expression")
barplot_10_group_4_datasets(groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC, "./1-figures/3-barplots/Zscore_quantiles_expression.png", [0,0.8])
boxplot_10_group_4_datasets(groups_GM12878,
	groups_K562,
	groups_IMR90,
	groups_HUVEC,
	"./1-figures/3-boxplots/Zscore_quantiles_expression.png",
	)
print("\n")
for i in range(5):
	print(len(groups_K562[i]))



groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC = get_boxplot_data_by_Zscore("QNed_SDOCs", "log2_mean_expression")
barplot_10_group_4_datasets(groups_GM12878, groups_K562, groups_IMR90, groups_HUVEC, "./1-figures/3-barplots/QNed_SDOC_quantiles_expression.png", [0,0.8])
boxplot_10_group_4_datasets(groups_GM12878,
	groups_K562,
	groups_IMR90,
	groups_HUVEC,
	"./1-figures/3-boxplots/QNed_SDOC_quantiles_expression.png",
	)

for i in range(5):
	print(len(groups_K562[i]))


######## 4 analysis compartment along with other epigenetic features
	
def epi_boxplot_5groups(df, feature_to_plot, file):

    plt.figure(figsize=(1.5, 3), dpi=600)
    plt.rc('font',family='Arial')
    plt.rc('font',size = 8)
    boxprops = {'linewidth':1.5 , 'color': 'black'}
    whiskerprops = {'linestyle': '--', 'linewidth':1.1, 'color': 'black'}
    medianprops = {'linewidth':1 , 'color': 'black', 'drawstyle' : 'steps'}


    datas = [df.loc[df["group"] == i, feature_to_plot] for i in range(5)]
    bp1 = plt.boxplot(datas, positions = [i for i in range(5)] ,showfliers=False, meanline =False, patch_artist = True, whiskerprops = whiskerprops, medianprops = medianprops,widths = 0.77,whis = [5, 95])
    cmap = mpl.cm.get_cmap("RdYlBu_r")
    color_bar = [cmap((i+1)/5) for i in range(5)]
    for patch, color in zip(bp1['boxes'], color_bar):
        patch.set_facecolor(color)
    #plt.scatter([i for i in range(6)], [np.mean(data1[i]) for i in range(6)], s = 5, color = "blue", zorder = 5)

    plt.xticks([i for i in range(5)],["<-1.5σ","-1.5σ~-0.5σ","-0.5σ~0.5σ","0.5σ~1.5σ",">1.5σ"], rotation = 45)
    #if ylim != []:
    #    plt.ylim(ylim)

    plt.tight_layout()
    plt.savefig(file+'.png')
    plt.close()




#sort epigenetic chart:
f = open("./data/epigenetics/2-sorted_TAD_epi_chart.txt")   
out = []		
lines=f.readlines() 
nrow = len(lines)					
for i in range(len(lines)):
    out.append(lines[i].strip())
f.close()
out[1:] = sort_list(out[1:])
savetxt("./data/epigenetics/sorted_epi_chart.tsv", out)

tmp_pd_GM12878_proporties = pd.read_csv('./data/epigenetics/sorted_epi_chart.tsv', sep='\t')
print(tmp_pd_GM12878_proporties.head())

df_GM12878["TAD_length"] = tmp_pd_GM12878_proporties["TAD_length"]
df_GM12878["Alu"] = tmp_pd_GM12878_proporties["Alu"]
df_GM12878["GM12878_CTCF"] = tmp_pd_GM12878_proporties["GM12878_CTCF"]
df_GM12878["GM12878_RNAseq"] = tmp_pd_GM12878_proporties["GM12878_RNAseq"]
df_GM12878["GM12878_H3K27me3"] = tmp_pd_GM12878_proporties["GM12878_H3K27me3"]
df_GM12878["GM12878_DHS"] = tmp_pd_GM12878_proporties["GM12878_DHS"]
df_GM12878["GM12878_H3K4me3"] = tmp_pd_GM12878_proporties["GM12878_H3K4me3"]
df_GM12878["GM12878_MRE"] = tmp_pd_GM12878_proporties["GM12878_MRE"]

# get group number of GM12878 TADs by QNed SDOC
df_GM12878["group"] = [0 for i in range(len(df_GM12878))]
df_GM12878.loc[df_GM12878["QNed_SDOCs"] < -1.5 ,"group"] = 0
df_GM12878.loc[np.abs(df_GM12878["QNed_SDOCs"] + 1) < 0.5 ,"group"] = 1
df_GM12878.loc[np.abs(df_GM12878["QNed_SDOCs"]) < 0.5 ,"group"] = 2
df_GM12878.loc[np.abs(df_GM12878["QNed_SDOCs"] - 1) < 0.5 ,"group"] = 3
df_GM12878.loc[df_GM12878["QNed_SDOCs"] > 1.5,"group"] = 4

for i in range(5):
	print(len(df_GM12878.loc[df_GM12878["group"] == i,"group"]))



mkdir("./1-figures/4-boxplot_6group_features")
fratures_to_plot = ["TAD_length","Alu","GM12878_CTCF","GM12878_RNAseq","GM12878_H3K27me3","GM12878_DHS","GM12878_H3K4me3","GM12878_MRE"]
for feature in fratures_to_plot:
	epi_boxplot_5groups(df_GM12878, feature, "./1-figures/4-boxplot_6group_features/"+feature)

compartments = []
f = open("0-TAD_compartment_GM12878.tsv")   		
lines=f.readlines() 
nrow = len(lines)					
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    compartments.append(L[3])
f.close()

df_GM12878["compartment"] = compartments
		

### get AB compartment ratio:

def compartment_ratio(df_do, celltype):
	compartments = []
	f = open("0-TAD_compartment_"+celltype+".tsv")   		
	lines=f.readlines() 
	nrow = len(lines)					
	for i in range(len(lines)):
	    L = lines[i].strip().split('\t')
	    compartments.append(L[3])
	f.close()
	
	df_do["compartment"] = compartments

	A_ratio = np.array([np.sum(df_do.loc[df_do["group"].isin([i]), "compartment"] == "A") for i in range(5)]) / len(df_do)
	B_ratio = np.array([np.sum(df_do.loc[df_do["group"].isin([i]), "compartment"] == "B") for i in range(5)]) / len(df_do)
	mixed_A_ratio = np.array([np.sum(df_do.loc[df_do["group"].isin([i]), "compartment"] == "mixed_A") for i in range(5)]) / len(df_do)
	mixed_B_ratio = np.array([np.sum(df_do.loc[df_do["group"].isin([i]), "compartment"] == "mixed_B") for i in range(5)]) / len(df_do)
	sum_TADs = np.copy(A_ratio + B_ratio + mixed_A_ratio + mixed_B_ratio)
	
	A_ratio /= sum_TADs
	B_ratio /= sum_TADs
	mixed_A_ratio /= sum_TADs
	mixed_B_ratio /= sum_TADs
	
	#do barplot:
	plt.figure(figsize=(2, 3), dpi=400)
	plt.rc('font',family='Arial')
	plt.rc('font',size = 9)
	
	plt.bar([i for i in range(5)], [1 for i in range(5)], color = [0.7,0,0], width = 0.77, linewidth = 0, edgecolor = 'black')
	plt.bar([i for i in range(5)], B_ratio + mixed_B_ratio + mixed_A_ratio, color = [1,0.6,0.3], width = 0.77)
	plt.bar([i for i in range(5)], B_ratio + mixed_B_ratio, color = [0,1,1], width = 0.77)
	plt.bar([i for i in range(5)], B_ratio, color = [0.0,0.4,0.7], width = 0.77)
	plt.xticks([i for i in range(5)],["<-1.5σ","-1.5σ~-0.5σ","-0.5σ~0.5σ","0.5σ~1.5σ",">1.5σ"], rotation = 45)
	plt.tight_layout()
	plt.savefig("./1-figures/compartment"+celltype+".png")
	plt.close()


compartment_ratio(df_GM12878 ,"GM12878")
compartment_ratio(df_IMR90 ,"IMR90")
compartment_ratio(df_K562 ,"K562")
compartment_ratio(df_HUVEC ,"HUVEC")


## analysis SE count per TAD

#get mean super-enhancer per TAD
def plot_SD_count_per_TAD(ct, df_ct):
	count_super_enhancer = []
	TAD_length = []
	f = open("./0-overlapped_SE/"+ct+"_10kb_"+ct+".txt")           
	lines=f.readlines() 
	nrow = len(lines)                   
	for i in range(len(lines)):
	    L = lines[i].strip().split('\t')
	    count = (len(L) - 6) / 5
	    count_super_enhancer.append(count)
	    TAD_length.append(float(L[2]) - float(L[1]))
	f.close()
	
	df_ct["SE"] = count_super_enhancer
	df_ct["TAD_length"] = TAD_length

	print(df_ct["group"].isin([i]))
	SE_density = [np.sum(df_ct.loc[df_ct["group"].isin([i]), "SE"]) / np.sum(df_ct.loc[df_ct["group"].isin([i]),"TAD_length"]) * 10 ** 7 for i in range(5)]
	print(SE_density)
	
	cmap = mpl.cm.get_cmap("RdYlBu_r")
	color_bar = [cmap((i+1)/5) for i in range(5)]
	
	mkdir("./1-figures/4-super-enhancers")
	plt.figure(figsize=(2, 3), dpi=600)
	plt.rc('font',family='Arial')
	plt.rc('font',size = 9)
	plt.bar([i for i in range(5)], SE_density, color = color_bar, width = 0.8, edgecolor = 'black', linewidth = 1)
	plt.ylabel("Super-enhancer density (counts per BP)")
	plt.xticks([i for i in range(5)],["<-1.5σ","-1.5σ~-0.5σ","-0.5σ~0.5σ","0.5σ~1.5σ",">1.5σ"], rotation = 45)
	plt.tight_layout()
	plt.savefig("./1-figures/4-super-enhancers/"+ct+".png")
	plt.savefig("./1-figures/4-super-enhancers/"+ct+".pdf")
	plt.close()

plot_SD_count_per_TAD("GM12878", df_GM12878)
plot_SD_count_per_TAD("IMR90", df_IMR90)
plot_SD_count_per_TAD("K562", df_K562)
plot_SD_count_per_TAD("HUVEC", df_HUVEC)


