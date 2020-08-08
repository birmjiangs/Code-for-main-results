import os
import sys
import numpy as np
import scipy.stats as stats

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
# get GSM number to celltype

mkdir("./pre2-out")

GSM_2_celltype = {}
celltypes = ["HSC","MPP","CLP","ETP","DN2","DN3","DN4","DP"]
f = open("./data_PC1_from_GSE79422/accession_number_to_dataset.tsv")           
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    number = L[0]
    celltype = L[1].split("_")[1]
    GSM_2_celltype[number] = celltype
f.close()


# get common loci wiith PC1 data
loci_2_PC1_value = {}
loci_2_count = {}
celltypes_each_dataset = []

files = os.listdir("./data_PC1_from_GSE79422/extracted")     
for fl in files:
    f = open("./data_PC1_from_GSE79422/extracted/"+str(fl)) 
    GSM_number = fl.split("_")[0]
    celltype = GSM_2_celltype[GSM_number]   
    celltypes_each_dataset.append(celltype)

    lines=f.readlines() 
    nrow = len(lines)               
    for i in range(len(lines)):
        L = lines[i].strip().split()
        loci = L[0]+'\t'+L[1]+'\t'+L[2]
        PC1 = float(L[3])

        if loci not in loci_2_count:
            loci_2_count[loci] = 1
            loci_2_PC1_value[loci] = [PC1]
        else:
            loci_2_count[loci] += 1
            loci_2_PC1_value[loci].append(PC1)
    f.close()

count_datasets = len(files)
common_loci_2_PC1s = {}
for loci in loci_2_count:
    if loci_2_count[loci] != count_datasets:
        continue
    else:
        common_loci_2_PC1s[loci] = loci_2_PC1_value[loci]

#######################
# merge PC1 data within the same stage:
celltypes_each_dataset = np.array(celltypes_each_dataset)

common_loci_2_mean_PC1s_8stages = {}
common_loci_2_PC1s_8stages = {}  #specific PC1 data for AVONA test

for ct in celltypes:
    wanted_indeces = (celltypes_each_dataset == ct)
    print(wanted_indeces)
    for loci in common_loci_2_PC1s:

        PC1s = np.array(common_loci_2_PC1s[loci])[wanted_indeces].tolist()
        mean_PC1 = np.mean(PC1s)

        if loci not in common_loci_2_mean_PC1s_8stages:
            common_loci_2_mean_PC1s_8stages[loci] = [mean_PC1]
            common_loci_2_PC1s_8stages[loci] = [PC1s]
        else:
            common_loci_2_mean_PC1s_8stages[loci].append(mean_PC1)
            common_loci_2_PC1s_8stages[loci].append(PC1s)

out_loci_8stage_mean_PC1 = []
for loci in common_loci_2_mean_PC1s_8stages:
    out_loci_8stage_mean_PC1.append(loci+'\t'+"\t".join(list(map(str,common_loci_2_mean_PC1s_8stages[loci]))))
out_loci_8stage_mean_PC1 = sort_list(out_loci_8stage_mean_PC1)

savetxt("./pre2-out/8stages_mean_PC1.tsv", out_loci_8stage_mean_PC1)
mkdir("./pre2-tmp")
savetxt("./pre2-tmp/8stages_mean_PC1.tsv", out_loci_8stage_mean_PC1)


# do ANOVA test for each loci:
loci_2_significant_change = {}
for loci in common_loci_2_mean_PC1s_8stages:
    PC1s_all_stages = common_loci_2_PC1s_8stages[loci]
    stat, pvalue = stats.f_oneway(PC1s_all_stages[0],PC1s_all_stages[1],PC1s_all_stages[2],PC1s_all_stages[3],PC1s_all_stages[4],PC1s_all_stages[5],PC1s_all_stages[6],PC1s_all_stages[7])
    if pvalue < 0.05:
        loci_2_significant_change[loci] = 1
    else:
        loci_2_significant_change[loci] = 0


# check if all replicate has a A-B or B-A in at least one stage:
loci_2_all_reps_same_flip = {}
for loci in common_loci_2_PC1s_8stages:
    loci_2_all_reps_same_flip[loci] = 0 #init 

    PC1s_all_stages = common_loci_2_PC1s_8stages[loci]

    sign_last_stage = 0
    same_sign_last_stage = 0

    for i in range(8): #each stage
        PC1s_of_all_reps_in_this_stage = PC1s_all_stages[i]

        sign = 0
        same_sign = 1
        for j in range(len(PC1s_of_all_reps_in_this_stage)):  #go through each replicate in this stage to see if they share the same sign of PC1s
            if j == 0:
                if PC1s_of_all_reps_in_this_stage[j] >= 0:
                    sign = 1
                if PC1s_of_all_reps_in_this_stage[j] < 0:
                    sign = -1
                continue
            
            if PC1s_of_all_reps_in_this_stage[j] * sign < 0:
                same_sign = 0
                break

        if i != 0:
            #compare with last stage if last stage has same sign PC1 for all replicates
            if same_sign_last_stage == 0 or same_sign == 0:
                pass
            elif same_sign_last_stage == 1 and same_sign == 1:
                if sign_last_stage * sign == -1:  # this loci has an all-rep-A2B-or-B2A-change
                    loci_2_all_reps_same_flip[loci] += 1

        sign_last_stage = sign
        same_sign_last_stage = same_sign # 0: not same sign, 1: same sign



# identify 4 types of loci: stable A, stable B, A-B, B-A and transient flips
loci_2_type = []
loci_2_all_rep_PC1 = []
for i in range(len(out_loci_8stage_mean_PC1)):
    L = out_loci_8stage_mean_PC1[i].split("\t")
    loci = L[0]+'\t'+L[1]+'\t'+L[2]

    loci_2_all_rep_PC1.append([loci]+common_loci_2_PC1s_8stages[loci])

    all_rep_same_flip = loci_2_all_reps_same_flip[loci]

    loci_2_type.append("")

    stages_PC1 = list(map(float, L[3:]))
    count_A_2_B = 0
    count_B_2_A = 0
    for j in range(7):
        if stages_PC1[j] <= 0 and stages_PC1[j+1] > 0:
            count_B_2_A = 1
        elif stages_PC1[j] > 0 and stages_PC1[j+1] <= 0:
            count_A_2_B = 1
    if count_B_2_A + count_A_2_B == 0:
        if np.mean(stages_PC1) > 0:
            loci_2_type[-1] = loci+'\t'+"stable A" + '\t' + str(all_rep_same_flip)
        else:
            loci_2_type[-1] =  loci+'\t'+"stable B" + '\t' + str(all_rep_same_flip)
    elif count_A_2_B == 1 and count_B_2_A == 0:
        if loci_2_significant_change[loci] == 1:
            loci_2_type[-1] =  loci+'\t'+"A to B" + '\t' + str(all_rep_same_flip)
        else:
            loci_2_type[-1] =  loci+'\t'+"A to B N.S." + '\t' + str(all_rep_same_flip)
    elif count_B_2_A == 1 and count_A_2_B == 0:
        if loci_2_significant_change[loci] == 1:
            loci_2_type[-1] =  loci+'\t'+"B to A" + '\t' + str(all_rep_same_flip)
        else:
            loci_2_type[-1] =  loci+'\t'+"B to A N.S." + '\t' + str(all_rep_same_flip)
    else:
        if loci_2_significant_change[loci] == 1:
            loci_2_type[-1] =  loci+'\t'+"transient flips" + '\t' + str(all_rep_same_flip)
        else:
            loci_2_type[-1] =  loci+'\t'+"transient flips N.S." + '\t' + str(all_rep_same_flip)

savetxt("./pre2-out/loci_PC1_types.tsv", loci_2_type)
savetxt("./pre2-out/loci_all_reps_PC1s.tsv", loci_2_all_rep_PC1)

### get A/B compartment proportion for each TAD in each stage:
os.system("python ./overlap.py ./pre1-TADs.tsv ./pre2-tmp/ 0 0 ./pre2-overlapped")


TAD_2_count_A_bin_8stages = {}
TAD_2_count_B_bin_8stages = {}
TAD_2_B_proportion_8stages = {}
f = open("./pre2-overlapped/pre1-TADs.tsv_8stages_mean_PC1.tsv.txt")           
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    TAD = L[0] + '\t' + L[1] + '\t' + L[2]
    TAD_2_count_A_bin_8stages[TAD] = [0 for n in range(8)]
    TAD_2_count_B_bin_8stages[TAD] = [0 for n in range(8)]
    if len(L) == 3:
        TAD_2_B_proportion_8stages[TAD] = [float("nan") for n in range(8)]
        continue

    TAD_start = L[1]
    TAD_end = L[2]
    for j in range(6,len(L),11):
        if L[j-2] == TAD_end or L[j-1] == TAD_start:  #not really overlapped
            continue

        for n in range(8):
            PC1 = float(L[j+n])
            if PC1 >= 0:
                TAD_2_count_A_bin_8stages[TAD][n] += 1
            elif PC1 < 0:
                TAD_2_count_B_bin_8stages[TAD][n] += 1


    TAD_2_B_proportion_8stages[TAD] = [TAD_2_count_B_bin_8stages[TAD][n]/(TAD_2_count_B_bin_8stages[TAD][n]+TAD_2_count_A_bin_8stages[TAD][n]) for n in range(8)]

out_TAD_8stage_B_proportion = []
for TAD in TAD_2_B_proportion_8stages:
    out_TAD_8stage_B_proportion.append(TAD+'\t'+"\t".join(list(map(str,TAD_2_B_proportion_8stages[TAD]))))

savetxt("./pre2-out/result_TAD_8stage_B_proportion.tsv", out_TAD_8stage_B_proportion)