import os
import sys
import numpy as np
import scipy.stats as stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import AgglomerativeClustering

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
	list_out = sorted(list_in, key=lambda items: items[2])
	return list_out

def sort_list2(list_in): #use when using '\t' to seperate items
	list_out = sorted(list_in, key=lambda items: items[2].split("\t")[0])
	return list_out

def do_boxplot(data,colorBar,plot_name,add_significance = False):
	mkdir("./b5-boxplot/")

	plt.figure(figsize=(2.2, 4), dpi=600)
	plt.rc('font',family='Arial')
	plt.rc('font',size = 9)
	ax = plt.subplot(111)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	boxprops = {'linewidth':1.5 , 'color': 'black'}
	whiskerprops = {'linestyle': '--', 'linewidth':1.1, 'color': 'black'}
	medianprops = {'linewidth':1 , 'color': 'black', 'drawstyle' : 'steps'}
	bp = ax.boxplot(data, showfliers=False, meanline =False, patch_artist = True, whiskerprops = whiskerprops, medianprops = medianprops,widths = 0.77,whis = [5, 95])
	for patch, color in zip(bp['boxes'], colorBar):
		patch.set_facecolor(color)
	
	if add_significance != False:
		if add_significance == True:
			add_significance = [0,1]
		pos1 = add_significance[0]+1
		pos2 = add_significance[1]+1
		#do mann whitney u test
		print(np.mean(data[0]),np.mean(data[1]))
		try:
			st, pv = stats.mannwhitneyu(data[0], data[1], use_continuity=True, alternative=None)
			print(pv,"rep1,1-2",marker,species)
		except:
			pass
		# plot significance:
		highest = np.max([np.percentile(data[0],95), np.percentile(data[1],95)]) * 1.05
		ax.plot([pos1,pos2], [highest,highest], "-", color = "black", linewidth = 1)
		plt.ylim([None,highest*1.15])
		ax.text(-0.05+pos1,highest*1.03,'p = {:0.1e}'.format(pv))

	plt.xticks([i+1 for i in range(len(data))],["Dec SDOC", "random"],rotation = 0)
	plt.title(plot_name)
	plt.tight_layout()
	plt.savefig("./b5-boxplot/"+plot_name+".png")
	plt.savefig("./b5-boxplot/"+plot_name+".pdf")
	plt.savefig("./b5-boxplot/"+plot_name+".eps")

def do_boxplot_all(data,colorBar,plot_name,add_significance = False):
	mkdir("./b5-boxplot/")

	plt.figure(figsize=(6, 4), dpi=600)
	plt.rc('font',family='Arial')
	plt.rc('font',size = 9)
	ax = plt.subplot(111)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	boxprops = {'linewidth':1.5 , 'color': 'black'}
	whiskerprops = {'linestyle': '--', 'linewidth':1.1, 'color': 'black'}
	medianprops = {'linewidth':1 , 'color': 'black', 'drawstyle' : 'steps'}
	bp = ax.boxplot(data, showfliers=False, meanline =False, patch_artist = True, whiskerprops = whiskerprops, medianprops = medianprops,widths = 0.77,whis = [5, 95])
	for patch, color in zip(bp['boxes'], colorBar):
		patch.set_facecolor(color)
	
	if add_significance != False:
		if add_significance == True:
			add_significance = [0,1]
		pos1 = add_significance[0]+1
		pos2 = add_significance[1]+1
		#do mann whitney u test
		print(np.mean(data[0]),np.mean(data[1]))
		try:
			st, pv = stats.mannwhitneyu(data[0], data[1], use_continuity=True, alternative=None)
			print(pv,"rep1,1-2",marker,species)
		except:
			pass
		# plot significance:
		highest = np.max([np.percentile(data[0],95), np.percentile(data[1],95)]) * 1.05
		ax.plot([pos1,pos2], [highest,highest], "-", color = "black", linewidth = 1)
		plt.ylim([None,highest*1.15])
		ax.text(-0.05+pos1,highest*1.03,'p = {:0.1e}'.format(pv))

	#plt.xticks([i+1 for i in range(len(data))],["Dec SDOC", "random"],rotation = 0)
	plt.title(plot_name)
	plt.tight_layout()
	plt.savefig("./b5-boxplot/"+plot_name+".png")
	plt.savefig("./b5-boxplot/"+plot_name+".pdf")
	plt.savefig("./b5-boxplot/"+plot_name+".eps")
	plt.close()
################
lineages = ["HSC","MPP","CLP","ETP","DN2","DN3","DN4","DP"]
resolution = 20000
TAD_distance_low_cutoff = 2000000/resolution  # 3Mb distance at least
TAD_distance_high_cutoff = 50000000/resolution

TADpair_2_count_stages = {}
TADpair_2_8stages_contact = {}

for n in range(len(lineages)):
	ct = lineages[n]
	print("reading data form", ct)
	f = open("./b4-TAD_pair_contact_each_stage/" + ct)   		
	lines=f.readlines() 
	nrow = len(lines)					
	for i in range(len(lines)):
		L = lines[i].strip().split('\t')
		TAD_pair_coord = "\t".join(L[0:6])
		distance = int(L[8])  #int(L[6])
		contact = float(L[9])	#float(L[7])

		if TAD_pair_coord not in TADpair_2_8stages_contact:
			TADpair_2_8stages_contact[TAD_pair_coord] = [0 for i in range(8)]
			TADpair_2_8stages_contact[TAD_pair_coord][n] = contact
		else:
			TADpair_2_8stages_contact[TAD_pair_coord][n] = contact

		if TAD_pair_coord not in TADpair_2_count_stages:
			TADpair_2_count_stages[TAD_pair_coord] = 1
		else:
			TADpair_2_count_stages[TAD_pair_coord] += 1

	   #@@ if TAD_pair_coord not in TADpair_2_8stages_SDOC1:
	   #@@ 	TADpair_2_8stages_SDOC1[TAD_pair_coord] = [0 for i in range(8)]
	   #@@ 	TADpair_2_8stages_SDOC1[TAD_pair_coord][n] = SDOC[0]
	   #@@ else:
	   #@@ 	TADpair_2_8stages_SDOC1[TAD_pair_coord][n] = SDOC[0]

	   #@@ if TAD_pair_coord not in TADpair_2_8stages_SDOC2:
	   #@@ 	TADpair_2_8stages_SDOC2[TAD_pair_coord] = [0 for i in range(8)]
	   #@@ 	TADpair_2_8stages_SDOC2[TAD_pair_coord][n] = SDOC[1]
	   #@@ else:
	   #@@ 	TADpair_2_8stages_SDOC2[TAD_pair_coord][n] = SDOC[1]
	
	f.close()

chart_contact_8stages = []
data = []
TAD_pair_coords = []
TAD_distance = []
for TAD_pair_coord in TADpair_2_8stages_contact:
	if TADpair_2_count_stages[TAD_pair_coord] != 8:  #missing TAD data in at least one stage
		continue
	if sum(TADpair_2_8stages_contact[TAD_pair_coord]) == 0: #no contact in all stages
		#print("zero contact in all stages, omit this pair")
		continue
	L = TAD_pair_coord.split("\t")
	TAD_distance.append(str(abs((float(L[1])+float(L[2])) / 2 - (float(L[4])+float(L[5])) / 2)))
	chart_contact_8stages.append(TAD_pair_coord+'\t'+TAD_distance[-1]+'\t'+"\t".join(list(map(str,TADpair_2_8stages_contact[TAD_pair_coord]))))
	TAD_pair_coords.append(TAD_pair_coord)
	data.append(TADpair_2_8stages_contact[TAD_pair_coord])

data = np.array(data)

for i in range(8):
	data[:,i] = stats.zscore(data[:,i])


mkdir("./b5-TADpair_contact_charts")
#@@savetxt("./b5-TADpair_contact_charts/chart.txt", chart_contact_8stages)
#@@

## takes approx 25s

###### part 2: analysis time-resolved TAD pair contact between TADs in specific TAD clusters
#### 1. get TADs in specific clusters
#@@TAD_2_increasing_cluster = {}
#@@TAD_2_decreasing_cluster = {}
#@@
#@@### cluster 7: increasing SDOC, cluster 2,3: decreasing SDOC
#@@f = open("./b3-data_out/TAD_cluster")   		
#@@SDOC_increasing_clusters = ["7"]
#@@SDOC_decreasing_clusters = ["2","3"]
#@@
#@@lines=f.readlines() 
#@@nrow = len(lines)					
#@@for i in range(len(lines)):
#@@	L = lines[i].strip().split('\t')
#@@	TAD = "\t".join(L[0:3])
#@@	cluster = L[3]
#@@	if cluster in SDOC_increasing_clusters:
#@@		TAD_2_increasing_cluster[TAD] = int(cluster)  #7
#@@	if cluster in SDOC_decreasing_clusters:
#@@		TAD_2_decreasing_cluster[TAD] = int(cluster) #2,3
#@@
#@@
#@@# get chart of TAD pairs formed by TADs in specific cluster(SDOC increasing or decreasing)
#@@chart_TAD_pairs_increasing_SDOC = []
#@@chart_TAD_pairs_decreasing_SDOC = []
#@@chart_TAD_pairs_random = []
#@@
#@@out_chart_TAD_pairs_increasing_SDOC = []
#@@out_chart_TAD_pairs_decreasing_SDOC = []
#@@out_chart_TAD_pairs_random = []
#@@
#@@TAD_distance_low_cutoff = 2000000/resolution  # 3Mb distance at least
#@@TAD_distance_high_cutoff = 50000000/resolution
#@@
#@@for i in range(len(chart_zscored_contact)):
#@@	line = chart_zscored_contact[i]
#@@	L = line.split("\t")
#@@	TAD1 = "\t".join(L[0:3])
#@@	TAD2 = "\t".join(L[3:6])
#@@	distance = int(float(L[6])/resolution)
#@@	contacts = list(map(float,L[7:15]))
#@@
#@@
#@@	if distance <= TAD_distance_low_cutoff:
#@@		continue
#@@	if distance >= TAD_distance_high_cutoff:
#@@		continue
#@@
#@@	if TAD1 in TAD_2_increasing_cluster and TAD2 in TAD_2_increasing_cluster:
#@@		chart_TAD_pairs_increasing_SDOC.append(contacts)
#@@		out_chart_TAD_pairs_increasing_SDOC.append(line)
#@@	elif TAD1 in TAD_2_decreasing_cluster and TAD2 in TAD_2_decreasing_cluster:
#@@		chart_TAD_pairs_decreasing_SDOC.append(contacts)
#@@		out_chart_TAD_pairs_decreasing_SDOC.append(line)
#@@	else:#elif i % 500 == 1:
#@@		chart_TAD_pairs_random.append(contacts)
#@@		out_chart_TAD_pairs_random.append(line)
#@@
#@@
#@@chart_TAD_pairs_increasing_SDOC = np.array(chart_TAD_pairs_increasing_SDOC)
#@@chart_TAD_pairs_decreasing_SDOC = np.array(chart_TAD_pairs_decreasing_SDOC)
#@@chart_TAD_pairs_random = np.array(chart_TAD_pairs_random)
#@@
#@@print("increasing_SDOC:", np.mean(chart_TAD_pairs_increasing_SDOC, axis = 0))
#@@print("decreasing_SDOC:", np.mean(chart_TAD_pairs_decreasing_SDOC, axis = 0))
#@@print("random:", np.mean(chart_TAD_pairs_random, axis = 0))
#@@
#@@mkdir("./b5-data_out")
#@@savetxt("./b5-data_out/increasing_SDOC_contacts", out_chart_TAD_pairs_increasing_SDOC)  #useless, change to TAD_coords + data
#@@savetxt("./b5-data_out/decreasing_SDOC_contacts", out_chart_TAD_pairs_decreasing_SDOC)
#@@savetxt("./b5-data_out/random_TAD_contacts",out_chart_TAD_pairs_random)
#@@
#@@#do mann-whitney U test
#@@stat,pv = stats.mannwhitneyu(chart_TAD_pairs_decreasing_SDOC[:,-1], chart_TAD_pairs_random[:,-1], use_continuity=True, alternative=None)
#@@print(pv)
#@@stat,pv = stats.mannwhitneyu(chart_TAD_pairs_increasing_SDOC[:,-1], chart_TAD_pairs_random[:,-1], use_continuity=True, alternative=None)
#@@print(pv)
#@@
#@@
#@@#do violin plot:
#@@
#@@colorbar = [[0.65,0,0],[0.65,0.65,0.65]]
#@@data_plot = [chart_TAD_pairs_increasing_SDOC[:,-1],chart_TAD_pairs_random[:,-1]]
#@@plot_name = "Increasing SDOC vs random"
#@@do_boxplot(data_plot,colorbar,plot_name,1)
#@@
#@@colorbar = [[0,0.3,0.65],[0.65,0.65,0.65]]
#@@data_plot = [chart_TAD_pairs_decreasing_SDOC[:,-1],chart_TAD_pairs_random[:,-1]]
#@@plot_name = "Decreasing SDOC vs random"
#@@do_boxplot(data_plot,colorbar,plot_name,1)
#@@
#@@colorbar = [[0,0.3,0.65],[0.65,0,0],[0.65,0.65,0.65],
#@@[0,0.3,0.65],[0.65,0,0],[0.65,0.65,0.65],
#@@[0,0.3,0.65],[0.65,0,0],[0.65,0.65,0.65],
#@@[0,0.3,0.65],[0.65,0,0],[0.65,0.65,0.65],
#@@[0,0.3,0.65],[0.65,0,0],[0.65,0.65,0.65],
#@@[0,0.3,0.65],[0.65,0,0],[0.65,0.65,0.65],
#@@[0,0.3,0.65],[0.65,0,0],[0.65,0.65,0.65],
#@@[0,0.3,0.65],[0.65,0,0],[0.65,0.65,0.65]]
#@@
#@@data_plot = []
#@@for i in range(8):
#@@	data_plot = data_plot + [chart_TAD_pairs_decreasing_SDOC[:,i],chart_TAD_pairs_increasing_SDOC[:,i],chart_TAD_pairs_random[:,i]]
#@@plot_name = "Inter TAD contact"
#@@do_boxplot_all(data_plot,colorbar,plot_name)


##################
#part3 test if significant rised(decreased) TAD-TAD contact are more likely to occur on TADs with high SDOC change
mkdir("./b5-loess_plot")

### 1. put TAD pairs into distance groups
#index_TAD_pairs:
TAD_pair_2_index = {}
index_2_TAD_pairs = {}
for i in range(len(chart_contact_8stages)):
	L = chart_contact_8stages[i].split("\t")
	TAD_pair_coord = "\t".join(L[0:6])
	TAD_pair_2_index[TAD_pair_coord] = i
	index_2_TAD_pairs[i] = TAD_pair_coord
######

normalized_TAD_pair_contact_chart = []
test_not_normalized_TAD_pair_contact_chart = []
for n in range(len(lineages)):
	min_distance = 1000
	max_distance = 1000
	ct = lineages[n]

	#construct dataset for loess fit
	dataset_for_loess = []
	distance_2_data = {}

	for i in range(len(chart_contact_8stages)):
		L = chart_contact_8stages[i].split("\t")
		distance_bin = int(float(L[6])/resolution)

		mean_contact = float(L[7+n])	
		TAD_pair_index = TAD_pair_2_index["\t".join(L[0:6])]
		dataset_for_loess.append([distance_bin,mean_contact,TAD_pair_index])  #x,y,index

		if distance_bin not in distance_2_data:
			distance_2_data[distance_bin] = [mean_contact]
		else:
			distance_2_data[distance_bin].append(mean_contact)

	dataset_for_loess = sorted(dataset_for_loess)
	dataset_for_loess = np.array(dataset_for_loess)
	#print("loess datasets:",dataset_for_loess[10:20])


	### get stdvar for each distance_bin, only keep distance with more than 2 values
	var_data_for_loess = []
	for distance_bin in distance_2_data:
		if len(distance_2_data[distance_bin]) < 3:
			continue

		min_distance = min(min_distance,distance_bin)  #for estimation of loess fitted variance
		max_distance = max(max_distance,distance_bin)  #for estimation of loess fitted variance

		var_data_for_loess.append([distance_bin,np.std(distance_2_data[distance_bin])])

	var_data_for_loess = sorted(var_data_for_loess)
	var_data_for_loess = np.array(var_data_for_loess)

	#print("loess datasets(variance):",var_data_for_loess[10:20])

	########### do lowess on mean
	lowess_mean = sm.nonparametric.lowess(var_data_for_loess[:,1], var_data_for_loess[:,0], frac=0.01, it = 0)

	#save lowess result into distance_2_loess_var, adding missing values(small distances with less than 3 pairs ,middle, large distance with less than 3 pairs)
	distance_2_loess_var = {}
	for i in range(len(lowess_mean[:,0])-1):
		distance_this = lowess_mean[i,0]
		distance_next = lowess_mean[i+1,0]
		for j in range(int(distance_this+1), int(distance_next)):
			distance_2_loess_var[j] = (lowess_mean[i,1]+lowess_mean[i+1,1])/2
		distance_2_loess_var[distance_this] = lowess_mean[i,1]
		distance_2_loess_var[distance_next] = lowess_mean[i+1,1]

	for i in range(0,min_distance):
		distance_2_loess_var[i] = distance_2_loess_var[min_distance]

	print(max_distance,max(dataset_for_loess[:,0])+1)
	for i in range(int(max_distance),int(max(dataset_for_loess[:,0]))+1):   #variance of distance larger than the max calculated fitted variance are set to the calculated fitted variance of the largest distance
		distance_2_loess_var[i] = distance_2_loess_var[max_distance]

	############### do lowess on variance
	lowess_var = sm.nonparametric.lowess(dataset_for_loess[:,1], dataset_for_loess[:,0], frac=0.01, it = 0)
	distance_2_loess_mean = {}
	for i in range(len(lowess_var[:,0])):  #for each lowess fitted point
		distance_this = lowess_var[i,0]
		fitted_mean = lowess_var[i,1]
		distance_2_loess_mean[distance_this] = fitted_mean


	###### save loess result on each distances with any TAD-TAD contact
	out_lowess_result = []
	for distance in distance_2_loess_mean:
		out_lowess_result.append([distance,distance_2_loess_mean[distance],distance_2_loess_var[distance]])
	out_lowess_result = sorted(out_lowess_result)

	mkdir("./b5-tmp_data")
	savetxt("./b5-tmp_data/"+ct+"_lowess_mean_var.txt", out_lowess_result)

	####### get TAD-TAD pairs in each stage with a "high" TAD-TAD contact frequency, and plot with red color
	high_contact_pair_index = []
	for i in range(len(dataset_for_loess)):
		distance = dataset_for_loess[i,0]

		if dataset_for_loess[i,1] > distance_2_loess_mean[distance] + 2 * distance_2_loess_var[distance]:
			high_contact_pair_index.append(i)
	high_contact_pair_index = np.array(high_contact_pair_index)



	fig = plt.figure(figsize=(4, 4), dpi=600)
	plt.rc('font',family='Arial')
	plt.rc('font',size = 9)
	ax = fig.add_subplot(1, 1, 1)
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.scatter(var_data_for_loess[:,0], var_data_for_loess[:,1], s = 5, alpha = 0.5, color = [0.7,0.7,0.7], linewidth = 0)
	ax.plot(lowess_mean[:, 0], lowess_mean[:, 1], 'black')
	plt.ylim([1e-2,1e1])
	plt.savefig("./b5-loess_plot/"+str(ct)+"_var.png")
	plt.close()
	

	fig = plt.figure(figsize=(4, 4), dpi=600)
	plt.rc('font',family='Arial')
	plt.rc('font',size = 9)
	ax = fig.add_subplot(1, 1, 1)
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.scatter(dataset_for_loess[:,0], dataset_for_loess[:,1], s = 0.5, alpha = 0.02, color = [0.7,0.7,0.7], linewidth = 0)
	ax.scatter(dataset_for_loess[high_contact_pair_index,0], dataset_for_loess[high_contact_pair_index,1], s = 0.5, alpha = 1, color = [1,0.0,0.0], linewidth = 0)
	ax.plot(lowess_var[:, 0], lowess_var[:, 1], 'black')
	plt.ylim([1e-2,1e2])
	plt.savefig("./b5-loess_plot/"+str(ct)+"_mean.png")
	plt.close()


	##### calculate normalized TAD-TAD comtact (normalized by distance) using fitted mean and variance
	dataset_for_loess_sorted = np.array(sort_list(np.copy(dataset_for_loess)))
	#print(dataset_for_loess_sorted[1:10])
	if n == 0: #first lineage, out put TAD coordinates and distance
		for i in range(len(dataset_for_loess_sorted)):
			distance_bin = dataset_for_loess_sorted[i,0]
			TAD_coord = index_2_TAD_pairs[dataset_for_loess_sorted[i,2]]
			normalized_contact = (dataset_for_loess_sorted[i,1] - distance_2_loess_mean[distance_bin]) / distance_2_loess_var[distance_bin]
			normalized_TAD_pair_contact_chart.append(TAD_coord+'\t'+str(distance_bin)+'\t'+ str(normalized_contact))
			test_not_normalized_TAD_pair_contact_chart.append(TAD_coord+'\t'+str(distance_bin)+'\t'+ str(dataset_for_loess_sorted[i,1]))
	else:
		for i in range(len(dataset_for_loess_sorted)):
			distance_bin = dataset_for_loess_sorted[i,0]
			normalized_contact = (dataset_for_loess_sorted[i,1] - distance_2_loess_mean[distance_bin]) / distance_2_loess_var[distance_bin]
			normalized_TAD_pair_contact_chart[i] = normalized_TAD_pair_contact_chart[i] + '\t' + str(normalized_contact)
			test_not_normalized_TAD_pair_contact_chart[i] = test_not_normalized_TAD_pair_contact_chart[i] +'\t'+ str(dataset_for_loess_sorted[i,1])
			
savetxt("./b5-TADpair_contact_charts/loess_normalized.txt", normalized_TAD_pair_contact_chart)
savetxt("./b5-TADpair_contact_charts/loess_test_not_normalized.txt", test_not_normalized_TAD_pair_contact_chart)

#### takes ~159s

###### get significantly altered TAD-TAD pairs: 

TAD_pair_for_clustering = []
TAD_pair_for_clustering_raw = []
data_for_clustering = []
data_raw = []
for i in range(len(normalized_TAD_pair_contact_chart)):
	L1 = test_not_normalized_TAD_pair_contact_chart[i].split("\t")
	data_raw = list(map(float,L1[7:15]))

	L = normalized_TAD_pair_contact_chart[i].split("\t")
	distance = int(float(L[6]))
	if distance < TAD_distance_low_cutoff or distance > TAD_distance_high_cutoff:
		continue

	L[7:] = list(map(float,L[7:]))

	diff = max(L[7:]) - min(L[7:])
	if diff < 2 or max(L[7:]) < 2 or np.max(data_raw) < 0.1:  #filter out unqualified TAD-TAD pairs
		continue
	TAD_pair_for_clustering.append(normalized_TAD_pair_contact_chart[i])
	TAD_pair_for_clustering_raw.append(test_not_normalized_TAD_pair_contact_chart[i])
	data_for_clustering.append(L[7:15])

f.close()
print("total TAD pairs for clustering:", len(data_for_clustering))

mkdir("b5-TAD_pair_clustering_result")
savetxt("./b5-TADpair_contact_charts/data_for_clustering.txt", normalized_TAD_pair_contact_chart)#@@
