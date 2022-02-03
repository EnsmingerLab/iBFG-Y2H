from __future__ import division
import os
import argparse
import datetime
import math
import numpy as np
import pandas as pd
import logging
import logging.config
import itertools
import noz_score
import param
import evaluation
import plot
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import matthews_corrcoef
#from noz_score import *


def freq(matrix):
        
    # sum of matrix
    total = matrix.values.sum()
    freq_df = matrix / total

    return freq_df


def marginal_freq(matrix):
    # sum of matrix
    total = matrix.values.sum()
    col_freq = matrix.sum(axis=0)/total
    row_freq = matrix.sum(axis=1)/total
    return row_freq, col_freq


def load_YI1(yi):
    "Load gold standard for yi1"
    yi = pd.read_table(yi, sep="\t")
    yi.columns = ["AD", "DB"]
    yi["Interactions"] = yi[['AD', 'DB']].apply(lambda x: '_'.join(x), axis=1)
    return yi


def rename(df):
    """remove _BC-* from orf names"""
    return df.rename(columns = lambda x : str(x)[:-5], index = lambda x : str(x)[:-5])


def calculate_freq(GFP_pre, GFP_high, GFP_med):
    
    # for GFP_pre
    row_freq, col_freq = marginal_freq(GFP_pre)
    # for GFP_med
    med_freq = freq(GFP_med)
    # for GFP_high
    high_freq = freq(GFP_high)

    med_freq_renamed = rename(freq(GFP_med))

    AD_NAMES = list(med_freq_renamed.index)
    DB_NAMES = list(med_freq_renamed)

    #heat_freq(med_freq)
    
    return row_freq, col_freq, med_freq, high_freq, AD_NAMES, DB_NAMES


def get_mcc(dicts, gold_st, AD_freq, DB_freq, row_cut, col_cut):
    # dicts contains intereactions and scores for each index 
    # compare to gold standard and plot PRC
    MCC = pd.DataFrame({}, columns=["precision","recall","mcc","rank"])
    AD_GOLD = gold_st.AD.tolist()
    DB_GOLD = gold_st.DB.tolist()
    for j in dicts.keys():
        i = dicts[j]
        #i["Interaction"] = i.index
        i = i.reset_index()
        # compare gold with our set
        i["ishit_interaction"] = i.Interaction.isin(gold_st.Interactions).astype(int)
        
        ## build screen set
        # AD has to appear in AD_GOLD at least once
        # DB has to appear in DB_GOLD at least once
        i["is_AD"] = i.AD_x.isin(AD_GOLD)
        i["is_DB"] = i.DB_x.isin(DB_GOLD)
         
        # frequencies should be greater than floow
        i = pd.merge(i, AD_freq, on="AD_name")
        i = pd.merge(i, DB_freq, on="DB_name")
        
        i.rename(columns={"0_x":"AD_mfreq", "0_y":"DB_mfreq"}, inplace=True)
        i["AD_cut"] = i.AD_mfreq > row_cut
        i["DB_cut"] = i.DB_mfreq > col_cut
        # get union
        
        i["screen"] = (i.is_AD & i.is_DB & i.AD_cut & i.DB_cut).astype(int)
        i = i.sort_values(by="Score", ascending=False)  
        # hit
        y_pred = i.ishit_interaction.tolist()
        # screen 
        y_network = i.loc[i.screen == 1].ishit_interaction.tolist()
        # scores
        y_score = i.Score.tolist()
        MAXMCC = prcmcc(y_network, 1000)
        MAXMCC["rank"] = j
        
        MCC = MCC.append(MAXMCC)
    MCC = MCC.reset_index(drop=True)

    return MCC


def prcmcc(label, test_range):
    # pred: predicted labels
    # label: actual labels
    PRCMCC = []
    total_screen = len(label)
    if total_screen < test_range:
        test_range = total_screen-1

    for i in range(1,test_range+1):

        test_screen = label[:i]
        test_nonscreen = label[i:]
        # true positive: condition positive and predicted pos
        TP = sum(test_screen)
        # true negative: condition negative and predicted neg
        TN = test_nonscreen.count(0)
        # false positive: condition negative and predicted pos
        FP = test_screen.count(0)
        # false negative: condition positive and predicted neg
        FN = sum(test_nonscreen)
        # precision = TP / (TP+FP)
        precision = sum(test_screen)/i * 100
        # recall = TP / (TP+FN)
        recall = sum(test_screen)/sum(label) * 100
        try:
            MCC= (TP*TN-FP*FN)/math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))*100 
            PRCMCC.append([precision, recall, MCC])
        except Exception:
            PRCMCC.append([np.nan, np.nan, np.nan])

    df = pd.DataFrame(PRCMCC, columns=["precision", "recall", "mcc"])

    #MAXid = df.MCC.idxmax()
    #print MAXid
    #MAXMCC = df.loc[MAXid].tolist()
    #all_mcc.append(MAXMCC) 
    #MAXMCC.insert(0, total_screen)
    return df


def get_norm_score(weight, high_freq, med_freq, pre_freq):
    IS = ((weight * high_freq) + med_freq) / pre_freq

    # replace with nan
    IS = IS.replace(0, np.nan)

    # log2 median of DB
    log_med = np.log2(IS.median(axis=0))

    # log2 of matrix
    log_is = np.log2(IS)

    # subtract
    IS_normed = log_is.sub(log_med)

    return IS_normed


def test_rank(IS_normed, all_pairs, mix_index):

    dicts = {"is_{}".format(i): {} for i in mix_index}
    all_pairs = pd.DataFrame(all_pairs, columns=['AD', 'DB'])
    all_pairs['Interaction'] = all_pairs.AD.str.cat(all_pairs.DB, sep="_")

    transform = IS_normed.unstack().reset_index()
    transform.columns = ['DB_name', 'AD_name', 'Score']
    # drop nan score
    #transform = transform.dropna(how='any')
    # split cols
    transform[['DB', 'DB_BC']] = transform['DB_name'].str.split('_', expand=True)
    transform[['AD', 'AD_BC']] = transform['AD_name'].str.split('_', expand=True)
    # merge cols
    transform['Interaction'] = transform.AD.str.cat(transform.DB, sep="_")
    merged = pd.merge(all_pairs, transform, on="Interaction")
    g = merged.sort_values(["Score"], ascending=False).groupby(['Interaction'])
    

    for i in mix_index:
        d_name = "is_{}".format(i)
        dicts[d_name] = g.nth(i).dropna(how='any')
    
    return dicts


def score_main(GFP_pre, GFP_high, GFP_med, weights, floor_perc, gold_st):
        
    # to find optimum set of parameters
    # we test all possible combinations
    # generates a list of tuples containing all possible combinations
    l = [weights, floor_perc]
    comb = list(itertools.product(*l))
    
    # get frequencies

    row_freq, col_freq, med_freq, high_freq, AD_NAMES, DB_NAMES = calculate_freq(GFP_pre, GFP_high, GFP_med)
    
    shape = GFP_pre.shape
    total_rows = shape[0] #AD
    total_cols = shape[1] #DB

    row_sorted = sorted(row_freq.tolist())
    col_sorted = sorted(col_freq.tolist())
    
    # get gold standard
    AD_GOLD = gold_st.AD.tolist()
    DB_GOLD = gold_st.DB.tolist()

    # find intersection
    AD_intersect = list(set(AD_GOLD) & set(AD_NAMES))
    DB_intersect = list(set(DB_GOLD) & set(DB_NAMES))
    
    all_pairs = list(itertools.product(AD_intersect, DB_intersect))
    
    mix_index = [0,1,2,3]
    
    # optimization
 
    output_csv = pd.DataFrame({}, columns=["precision","recall","mcc","rank","weight", "floor"])
    for s in comb:
        weight = s[0]
        floor = s[1]
        
        row_cut = row_sorted[int(round(total_rows/floor))]
        col_cut = col_sorted[int(round(total_cols/floor))]

        AD_freq = row_freq.where(row_freq > row_cut, row_cut)
        DB_freq = col_freq.where(col_freq > col_cut, col_cut)

        # rebuild matrix from two vectors
        freq_mx = np.outer(AD_freq, DB_freq)
        
        pre_freq = pd.DataFrame(data = freq_mx, columns = DB_freq.index.tolist(), index = AD_freq.index.tolist())
    
        IS_normed = get_norm_score(weight, high_freq, med_freq, pre_freq)
        
        # test which set of bc we want to use
        # in this case, we will have max 4 min 1 score(s) for each orf pair
        
        dicts = test_rank(IS_normed, all_pairs, mix_index)
        
        AD_freq = AD_freq.to_frame()
        DB_freq = DB_freq.to_frame()
        
        AD_freq["AD_name"] = AD_freq.index
        DB_freq["DB_name"] = DB_freq.index
         
        mcc_list = get_mcc(dicts, gold_st, AD_freq, DB_freq, row_cut, col_cut)
        mcc_list["weight"] = weight
        mcc_list["floor"] = floor
        output_csv = output_csv.append(mcc_list)
        #break 
    return output_csv, row_freq, col_freq, med_freq, high_freq, AD_NAMES, DB_NAMES


def load_summary(mcc_sum):
    mcc_summary = pd.read_csv(mcc_sum)
    MAX = mcc_summary.loc[mcc_summary["mcc"].idxmax()]
    max_weight = MAX["weight"]
    max_rank = MAX["rank"]
    max_floor = MAX["floor"]
    max_mcc = mcc_summary[(mcc_summary["weight"] == max_weight) & (mcc_summary["rank"] == max_rank) & (mcc_summary["floor"] == max_floor)].reset_index(drop=True)
    # plot max mcc
    title = max_rank+";w="+str(round(max_weight, 1))+";f="+str(round(max_floor, 1))
    #plot.plot_prc(max_mcc.precision, max_mcc.recall, "./dk_prc_curve_opt.png", title)
    #plot.plot_prcmcc(max_mcc, "./dk_prcmcc_curve_opt.png", title)
    #print max_mcc
    #print "plots made"
    return max_weight, max_rank, max_floor


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='BFG-Y2H Scores')
    parser.add_argument("--sample", help="Path for read counts files")

    args = parser.parse_args()
    counts = args.sample

    sample_name = os.path.basename(counts)    
    
    os.chdir(counts)
    # find the optimized parameters
    yi = load_YI1(param.GOLD)
    for f in os.listdir(counts):
        if not f.endswith("_combined_counts.csv"):
            continue
        fname = "./"+f        
        if "pre" in f:
            GFP_pre = pd.read_table(fname, sep=",", index_col=0)
        elif "med" in f:
            GFP_med = pd.read_table(fname, sep =",", index_col=0)
        elif "high" in f:
            GFP_high = pd.read_table(fname, sep =",", index_col=0)
    
    output_csv, row_freq, col_freq, med_freq, high_freq, AD_NAMES, DB_NAMES = score_main(GFP_pre, GFP_high, GFP_med, param.weights, param.floor_perc, yi)
    output_csv.to_csv("DK_mcc_summary_yi1.csv", index=False)
    max_weight, max_rank, max_floor = load_summary("DK_mcc_summary_yi1.csv")
    maxmcc = pd.DataFrame({"max_weight": [max_weight], "max_rank": [max_rank], "max_floor": [max_floor]})
    print "max param found"
    # evaluation
    litbm = evaluation.load_litbm(param.litBM13)
    evaluation.dk_main(litbm, max_weight, max_rank, max_floor, high_freq, med_freq, row_freq, col_freq, AD_NAMES, DB_NAMES, "dk_mcc_summary_litbm13")
    
    noz_main, raw_scores = noz_score.main(GFP_pre, GFP_med, GFP_high, yi)
#    df = raw_scores.copy()
#    df = df.unstack().reset_index()
#    df.to_csv("noz_raw_score.csv", index=False)
    noz_main.to_csv("noz_mcc_summary_yi1.csv", index=False)
    max_rho, max_rank = noz_score.load_summary("noz_mcc_summary_yi1.csv")
    maxmcc["max_rho"] = [max_rho]
    maxmcc["noz_max_rank"] = [max_rank]
    maxmcc.to_csv("max_parameters.csv", index=False)
    # evaluation    
    evaluation.noz_main(litbm, raw_scores, max_rank, max_rho, "noz_mcc_summary_litbm13")
