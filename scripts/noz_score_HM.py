#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 20:17:22 2019

@author: harleymount
"""
counts='/Users/harleymount/Desktop/V6/'
scripts=counts+'scripts/'
import os
os.chdir(counts)
os.chdir(scripts)

from __future__ import division
import pandas as pd
import numpy as np
import os
import datetime
import math
import plot
import param
import score
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# calculating scores based on nozomu's paper
# load the file in the same way as score.py
# detail about how to calculate the score on google doc

def load_YI1(yi):
    "Load gold standard for yi1"
    yi = pd.read_table(yi, sep="\t")
    yi.columns = ["AD", "DB"]
    yi["Interactions"] = yi[['AD', 'DB']].apply(lambda x: '_'.join(x), axis=1)
    return yi

def pre_freq(rc_pre):
    """
    frequency for GFP pre. same as non-selective condition 
    in Nozomu's paper: use marginal frequencies
    """
    total_reads = rc_pre.values.sum()
    col_freq = rc_pre.sum(axis=0)/total_reads
    row_freq = rc_pre.sum(axis=1)/total_reads
    return row_freq, col_freq


def freq(rc_matrix):
    """
    freq = read counts for this pair / total reads in this matrix
    """
    total_reads = rc_matrix.values.sum()
    # normalize the score by +1
    rc_matrix += 1
    freq = rc_matrix / total_reads
    return freq


def get_score(plus_his_freq,minus_his_freq):
    
    """
    calculate raw score
    """
    # convert pre_freq to dataframe
    plus_his_freq = pd.DataFrame(data = plus_his_freq, columns=minus_his_freq.columns.tolist(), index=minus_his_freq.index.tolist())

    s = (minus_his_freq)/plus_his_freq
    return s

def norm_score(raw_s, q):
    """
    Normalize the raw score matrix s
    q: float, quantile 
    """
    s = raw_s.copy()
    # get median of all DB scores 
    med = s.median(axis=0)
    beta = s[s>med]-med
    
    beta_q = beta.quantile(q)
    # if s_ij - med(s_i) < beta_i
    m1 = (s-med) <= beta_q
    # if s_ij - med(s_i) >= beta_i
    m2 = (s-med) > beta_q

    s[m1] = 1
    s[m2] = (s-med)/beta_q
    s = s.replace(np.inf, np.nan) 
    return s


def get_rank(norm_s, rank=range(0,9)):  # note i changed here to 9 because we have 3x3 barcodes in any cases
    """
    In total we have 4 ranks for each protein pair
    aBC1-bBC1, aBC1-bBC2, aBC2-bBC1, aBC2-bBC2 
    norm_s: normalized score s'
    rank: range(0, 4)
    """
    output={}
    # unstack the matrix
    transform = norm_s.unstack().reset_index()
    # rename col
    transform.columns = ['DB', 'AD', 's_prime']
    # split cols
    transform[['DB', 'DB_BC']] = transform['DB'].str.split('_', expand=True)
    transform[['AD', 'AD_BC']] = transform['AD'].str.split('_', expand=True)
    
    # merge cols
    transform['Interaction'] = transform.AD.str.cat(transform.DB, sep="_")
    
    # sort by group
    sort = transform.sort_values(["s_prime"], ascending=False).groupby("Interaction")

    for i in rank:
        d_name="rank_{}".format(i)
        scores = sort.nth(i).dropna(how="any").sort_values(["s_prime"], ascending=False)
        #output[d_name] = scores[scores.s_prime > 1]

        output[d_name] = scores
    
    return output


def get_screen(dict_s, gold_st):

    """
    Get screening set
    """
    #dict_s=output_rank
    prcmcc = pd.DataFrame({}, columns=["precision","recall","mcc","rank"])
    AD_GOLD = gold_st.AD.tolist()
    DB_GOLD = gold_st.DB.tolist()
    #key='rank_0'
    for key in dict_s.keys():
        s_prime = dict_s[key]
        s_prime = s_prime.reset_index()
        s_prime["is_hit"] = s_prime.Interaction.isin(gold_st.Interactions).astype(int)
        
        s_prime["is_AD"] = s_prime.AD.isin(AD_GOLD)
        s_prime["is_DB"] = s_prime.DB.isin(DB_GOLD)

        s_prime["screen"] = (s_prime.is_AD & s_prime.is_DB).astype(int)
        s_prime = s_prime[s_prime["s_prime"]>1]
        network_orfs = s_prime.loc[s_prime.screen == 1].is_hit.tolist()
        #print sum(network_orfs)
        #print sum(s_prime.is_hit.tolist())

        MAXMCC = prcmcc_fun(network_orfs, 1000) # have to make sure score.prcmcc is same as used in DK script
        MAXMCC["rank"] = key
        prcmcc = prcmcc.append(MAXMCC)
    
    prcmcc = prcmcc.reset_index(drop=True)
    return prcmcc

def main(plus_his, minus_his, gold_st):

    #calculate GFP_pre freq
    row_freq, col_freq = pre_freq(plus_his)
    plus_his_freq = np.outer(row_freq, col_freq)
    minus_his_freq = freq(minus_his)


    # test plot
    # should be very little corr
    # freq_corr(high_freq, med_freq)
    
    # get raw scores
    raw_s = get_score(plus_his_freq, minus_his_freq)
    df = raw_s.copy()
    df = df.unstack().reset_index()
    df.to_csv("noz_raw_score.csv", index=False)
    #test = s.values.flatten()
    #test.sort()
    #plot_diff(test)
    
    percentile = np.arange(0.1, 0.8, 0.05)
    mcc_summary = pd.DataFrame({}, columns=["precision","recall","mcc","rank","rho"])
    # test different rho to optimize mcc
    #p=0.75
    
    for p in percentile:
        # get normalized scores
        norm_s = norm_score(raw_s, p)
        sample_name = "rho_"+ str(p)
        # test corr of score and normed score
        # norm_score_corr(sample_name, s, norm_s)
        output_rank = get_rank(norm_s, rank=range(0,4))
        mcc = get_screen(output_rank, gold_st)
        mcc["rho"] = p
        mcc_summary = mcc_summary.append(mcc)
        #print datetime.datetime.now()
        #print mcc_summary
        #break
    return mcc_summary, raw_s

def load_summary(mcc_sum):
    mcc_summary = pd.read_csv(mcc_sum)
    MAX = mcc_summary.loc[mcc_summary["mcc"].idxmax()]
    max_rho = MAX.rho
    max_rank = MAX["rank"]
    max_mcc = mcc_summary[(mcc_summary.rho == max_rho) & (mcc_summary["rank"] == max_rank)].reset_index(drop=True)
    # plot max mcc
    title = max_rank+";rho="+str(round(max_rho, 1))
#    plot.plot_prc(max_mcc.precision, max_mcc.recall, "./noz_prc_curve_opt.png", title)
#    plot.plot_prcmcc(max_mcc, "./noz_prcmcc_curve_opt.png", title)
#    print "plots made"
    return max_rho, max_rank



def prcmcc_fun(y_network, test_range):
    # pred: predicted labels
    # label: actual labels
    PRCMCC = []
    total_screen = len(y_network)
    if total_screen < test_range:
        test_range = total_screen-1
#i=100
    for i in range(1,test_range+1):

        test_screen = y_network[:i]
        test_nonscreen = y_network[i:]
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
        recall = sum(test_screen)/sum(y_network) * 100
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



#=========================Self written code

def main_final(plus_his, minus_his, gold_st):

    #calculate GFP_pre freq
    row_freq, col_freq = pre_freq(plus_his)
    plus_his_freq = np.outer(row_freq, col_freq)
    minus_his_freq = freq(minus_his)


    # test plot
    # should be very little corr
    # freq_corr(high_freq, med_freq)
    
    # get raw scores
    raw_s = get_score(plus_his_freq, minus_his_freq)
    df = raw_s.copy()
    df = df.unstack().reset_index()
    df.to_csv("noz_raw_score.csv", index=False)
    #test = s.values.flatten()
    #test.sort()
    #plot_diff(test)
    
    percentile = ho
    p=percentile
    mcc_summary = pd.DataFrame({}, columns=["precision","recall","mcc","rank","rho"])
    # test different rho to optimize mcc
    #p=0.75
    

    # get normalized scores
    norm_s = norm_score(raw_s, p)
    sample_name = "rho_"+ str(p)
    # test corr of score and normed score
    # norm_score_corr(sample_name, s, norm_s)
    output_rank = get_rank_final(norm_s, max_rank[-1])
    mcc = get_screen(output_rank, gold_st)
    mcc["rho"] = p
    mcc_summary = mcc_summary.append(mcc)
    #print datetime.datetime.now()
    #print mcc_summary
    #break
    return mcc_summary, raw_s, output_rank






def noz_main(litbm, score_matrix, max_rank, max_rho, fname):
    # get normalized score by using max_pho
    norm_s = norm_score(score_matrix, max_rho)
    # save the normalized scores to file
    df = norm_s.unstack().reset_index()
    df.to_csv("noz_norm_score.csv", index=False)
    # get optimized rank by using max_rank
    r = [int(max_rank.split("_")[1])]
    output_rank = get_rank(norm_s, r)
    mcc = get_screen(output_rank, litbm)
    # store mcc to file
    mcc.to_csv(fname+".csv", index=False)
    # plot prcmcc



#rank=max_rank[-1]

def get_rank_final(norm_s, rank):  # note i changed here to 9 because we have 3x3 barcodes in any cases
    """
    In total we have 4 ranks for each protein pair
    aBC1-bBC1, aBC1-bBC2, aBC2-bBC1, aBC2-bBC2 
    norm_s: normalized score s'
    rank: range(0, 4)
    """
    output={}
    # unstack the matrix
    transform = norm_s.unstack().reset_index()
    # rename col
    transform.columns = ['DB', 'AD', 's_prime']
    # split cols
    transform[['DB', 'DB_BC']] = transform['DB'].str.split('_', expand=True)
    transform[['AD', 'AD_BC']] = transform['AD'].str.split('_', expand=True)

    # merge cols
    transform['Interaction'] = transform.AD.str.cat(transform.DB, sep="_")
    
    # sort by group
    sort = transform.sort_values(["s_prime"], ascending=False).groupby("Interaction")

    d_name="rank_{}".format(rank)
    scores = sort.nth(int(rank)).dropna(how="any").sort_values(["s_prime"], ascending=False)
    #output[d_name] = scores[scores.s_prime > 1]

    output[d_name] = scores

    return output







def plot_prcmcc(df, name, title):

    plt.plot(df.precision.tolist(), ".", label="precision", markersize=0.8)
    plt.plot(df.recall.tolist(), ".", label="recall", markersize=0.8)
    plt.plot(df.mcc.tolist(), ".", label="mcc", markersize=0.8)
    plt.xlabel("ORF pairs ranked by IS")
    plt.legend(loc="upper right")
    plt.title(title)
    plt.savefig(name)
    plt.close()


plot_prcmcc(prcmcc, 'test', 'test')
#=====================================



if __name__ == "__main__":
    # test on yAD4 DB1

    
    os.chdir(test_dir)
    for f in os.listdir(test_dir):
        if not f.endswith("_combined_counts.csv"):
            continue
        fname = "./"+f
        if "plus" in f:
            plus_his = pd.read_table(fname, sep=",", index_col=0)
        elif "minus" in f:
            minus_his = pd.read_table(fname, sep =",", index_col=0)
        
    gold_st = load_YI1(param.GOLD)
    mcc_summary, s = main(plus_his, minus_his, gold_st)            
    mcc_summary.to_csv("noz_mcc_summary_opt.csv", index=False)
    ho, max_rank = load_summary("noz_mcc_summary_opt.csv")

    noz_main, raw_scores = main(plus_his, minus_his, gold_st)            
    noz_main.to_csv("noz_mcc_summary_yi1.csv", index=False)
    max_rho, max_rank = load_summary("noz_mcc_summary_yi1.csv")
    maxmcc["max_rho"] = [max_rho]
    maxmcc["noz_max_rank"] = [max_rank]
    maxmcc.to_csv("max_parameters.csv", index=False)
    # evaluation    
    noz_main(litbm, raw_scores, max_rank, max_rho, "noz_mcc_summary_litbm13")