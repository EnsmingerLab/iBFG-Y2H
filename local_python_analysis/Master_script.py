#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 16:02:38 2019

@author: harleymount
"""

#I will compile all of the codes from Roth lab into a single python code to run analysis on output files


#------------------------------------------------------------------------------------------------------------------------
#first need to declare the local directory that will be used for analysis
#inside this directory there should be a scripts directory with all of the roth lab python scripts (most importantly param.py), should also be a ref and summary directory
counts='/Users/harleymount/Desktop/V6/'
scripts=counts+'scripts/'
import os
os.chdir(counts)
os.chdir(scripts)



#------------------------------------------------------------------------------------------------------------------------
from __future__ import division
import argparse
import datetime
import math
import numpy as np
import pandas as pd
import logging
import logging.config
import itertools
#excluding load noz_score because I will put it in this script as its own function, no inter python file functions 
import param
#import plot #need to fix plot to have python3 style print statements
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import matthews_corrcoef
import matplotlib.pyplot as plt
os.chdir(counts)
#============================================================================================================================================
#============================================================================================================================================
#Functions from score_HM_working.py
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


def rename(df):
    """remove _BC-* from orf names"""
    return df.rename(columns = lambda x : str(x)[:-5], index = lambda x : str(x)[:-5])


def calculate_freq(plus_his, minus_his):
    
    # for plus_his
    row_freq, col_freq = marginal_freq(plus_his)
    # for GFP_med
    minus_his_freq = freq(minus_his)
    #data_freq_renamed = rename(freq(data)) # not sure if I need this line after changing my summary files to V3 and updating locus col info
    AD_NAMES = list(minus_his_freq.index)
    #AD_NAMES=[x.upper() for x in AD_NAMES] # making names in AD Names uppercase
    DB_NAMES = list(minus_his_freq)
    #DB_NAMES=[x.upper() for x in DB_NAMES] # making names in DB Names uppercase
    data_freq_renamed = rename(freq(plus_his))
    AD_NAMES = list(data_freq_renamed.index)
    DB_NAMES = list(data_freq_renamed)

    #heat_freq(med_freq)
    return row_freq, col_freq, minus_his_freq, AD_NAMES, DB_NAMES



def get_norm_score(weight, minus_his_freq, pre_freq):
    IS = ((weight * minus_his_freq)) / pre_freq

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
#transform[[]]
    transform = IS_normed.unstack().reset_index()
    transform.columns = ['DB_name', 'AD_name', 'Score']
    # drop nan score
    #transform = transform.dropna(how='any')
    # split cols
    transform[['DB', 'DB_BC']] = transform['DB_name'].str.split('_', expand=True)[[0,1]] #------------------THIS IS WHERE PROBLEM IS ,I addwd [[0,1]]
    transform[['AD', 'AD_BC']] = transform['AD_name'].str.split('_', expand=True)[[0,1]]
    # merge cols
    transform['Interaction'] = transform.AD.str.cat(transform.DB, sep="_")
    merged = pd.merge(all_pairs, transform, on="Interaction")
    g = merged.sort_values(["Score"], ascending=False).groupby(['Interaction'])#up to this point all positive controls appear to beshowing up
    
#
    for i in mix_index:
        d_name = "is_{}".format(i)
        dicts[d_name] = g.nth(i).dropna(how='any')
    
    return dicts


#------my own stuff lookign at high index PPI
#g = transform.sort_values(["Score"], ascending=False).groupby(['Interaction'])#up to this point all positive controls appear to beshowing up

#g.nth(0).dropna(how='any')

#===============

def get_mcc(dicts, gold_st, AD_freq, DB_freq, row_cut, col_cut):
    # dicts contains intereactions and scores for each index 
    # compare to gold standard and plot PRC
    MCC = pd.DataFrame({}, columns=["precision","recall","mcc","rank"])
    AD_GOLD = gold_st.AD.tolist()
    DB_GOLD = gold_st.DB.tolist()
#j='is_2' # for troubleshooting for loop below
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


def prcmcc(y_network, test_range):
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


def score_main(plus_his, minus_his, weights, floor_perc, gold_st): #I will make a new score_main function that just uses one sample instead of 3
        
    # to find optimum set of parameters
    # we test all possible combinations
    # generates a list of tuples containing all possible combinations
    l = [weights, floor_perc]
    comb = list(itertools.product(*l))
    
    # get frequencies

    row_freq, col_freq, minus_his_freq, AD_NAMES, DB_NAMES = calculate_freq(plus_his, minus_his)#for this line I will need to change it to only do teh calculations for my one combined counts sample
    
    
    
    shape = plus_his.shape
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
    
    #mix_index = [0,1,2,3] # have to change here to test matrices 1-8
    mix_index = [0,1,2,3,4] # have to change here to test matrices 1-8
    
    # optimization
 
    #s=(1,2.1)
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
    #--------------------------------------
    
    
    
    
        IS_normed = get_norm_score(weight, minus_his_freq, pre_freq)
        
        # test which set of bc we want to use
        # in this case, we will have max 4 min 1 score(s) for each orf pair
        
        dicts = test_rank(IS_normed, all_pairs, mix_index) ###CODE FAILS HERE PROB DUE TO NAMES OF ROW AND COL
        
        AD_freq = AD_freq.to_frame()
        DB_freq = DB_freq.to_frame()
        
        AD_freq["AD_name"] = AD_freq.index
        DB_freq["DB_name"] = DB_freq.index
         
        mcc_list = get_mcc(dicts, gold_st, AD_freq, DB_freq, row_cut, col_cut)
        mcc_list["weight"] = weight
        mcc_list["floor"] = floor
        output_csv = output_csv.append(mcc_list)
        #break 
    return output_csv, row_freq, col_freq, minus_his_freq, AD_NAMES, DB_NAMES
    


def load_summary(mcc_sum):
    mcc_summary = pd.read_csv(mcc_sum)
    MAX = mcc_summary.loc[mcc_summary["mcc"].idxmax()]
    max_weight = MAX["weight"]
    max_rank = MAX["rank"]
    max_floor = MAX["floor"]
    max_mcc = mcc_summary[(mcc_summary["weight"] == max_weight) & (mcc_summary["rank"] == max_rank) & (mcc_summary["floor"] == max_floor)].reset_index(drop=True)
    # plot max mcc
    title = max_rank+";w="+str(round(max_weight, 1))+";f="+str(round(max_floor, 1))
    plot_prc(max_mcc.precision, max_mcc.recall, "./dk_prc_curve_opt.png", title)
    plot_prcmcc(max_mcc, "./dk_prcmcc_curve_opt.png", title)
    #print max_mcc
    #print "plots made"
    return max_weight, max_rank, max_floor
#============================================================================================================================================
#============================================================================================================================================
#code is from evaluation.py


def load_litbm(lit):
    """
    load litbm to pandas dataframe
    input cols: AD DB
    output cols: AD DB Interactions
    """
    litbm = pd.read_table(lit)
    litbm["Interactions"] = litbm[['AD', 'DB']].apply(lambda x: '_'.join(x), axis=1)
    return litbm



def dk_main(litbm, max_weight, max_rank, max_floor, minus_his_freq, row_freq, col_freq, AD_NAMES, DB_NAMES, fname):
    
    shape = minus_his_freq.shape
    total_rows = shape[0] #AD
    total_cols = shape[1] #DB

    AD_intersect = list(set(litbm.AD) & set(AD_NAMES))
    DB_intersect = list(set(litbm.DB) & set(DB_NAMES))
    
    all_pairs = list(itertools.product(AD_intersect, DB_intersect))
    
    row_sorted = sorted(row_freq.tolist())
    col_sorted = sorted(col_freq.tolist())
    row_cut = row_sorted[int(round(total_rows/max_floor))]
    col_cut = col_sorted[int(round(total_cols/max_floor))]

    AD_freq = row_freq.where(row_freq > row_cut, row_cut)
    DB_freq = col_freq.where(col_freq > col_cut, col_cut)

    freq_mx = np.outer(AD_freq, DB_freq)
    AD_freq = AD_freq.to_frame()
    DB_freq = DB_freq.to_frame()
        
    AD_freq["AD_name"] = AD_freq.index
    DB_freq["DB_name"] = DB_freq.index

    pre_freq = pd.DataFrame(data = freq_mx, columns = DB_freq.index.tolist(), index = AD_freq.index.tolist())
    
    IS_normed = get_norm_score(max_weight, minus_his_freq, pre_freq)
    df = IS_normed.unstack().reset_index()
    df.to_csv("DK_norm_score.csv", index=False) #this line outputs the entire list of scores, but I also want to output to optimal list based on MCC stuff

#=================HM custom code
    #to output the optimal list I will
    #add column names to output 
    df.columns = ['DB_name', 'AD_name', 'Score']
    #adjust the columns to have a BC # column 
    df[['DB', 'DB_BC']] = df['DB_name'].str.split('_', expand=True)
    df[['AD', 'AD_BC']] = df['AD_name'].str.split('_', expand=True)#[[0,1]]
    #add interaction column to describe PPI
    df['Interaction'] = df.AD.str.cat(df.DB, sep="_")
    #now sort interactions by score, and then group into interaction PPI groups
    g_all=df.sort_values(["Score"], ascending=False).groupby(['Interaction'])
#determine numerical value for best rank and use for scoring
    best_rank_int=int(max_rank.split('_')[1])
    g_all.nth(best_rank_int).dropna(how='any').sort_values(["Score"], ascending=False).to_csv("DK_IS_scores_Best_Rank_Max_Floor.csv", index=False)#line to get the best ranked interaction score matrix based on MCC by changing .nth(x)
#=================HM custom code ends
    r = [int(max_rank.split("_")[1])] 
    dicts = test_rank(IS_normed, all_pairs, r)
    mcc = get_mcc(dicts, litbm, AD_freq, DB_freq, row_cut, col_cut)
    mcc.to_csv(fname+".csv", index=False)



#============================================================================================================================================
#============================================================================================================================================
#code from noz_score.py

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


def freq_NOZ(rc_matrix): #going to change the name of this function to prevent confusion, used to be called freq, now called freq_NOZ
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
    minus_his_freq = freq_NOZ(minus_his)


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



def load_summary_NOZ(mcc_sum):
    mcc_summary = pd.read_csv(mcc_sum)
    MAX = mcc_summary.loc[mcc_summary["mcc"].idxmax()]
    max_rho = MAX.rho
    max_rank = MAX["rank"]
    max_mcc = mcc_summary[(mcc_summary.rho == max_rho) & (mcc_summary["rank"] == max_rank)].reset_index(drop=True)
    # plot max mcc
    title = max_rank+";rho="+str(round(max_rho, 1))
    plot_prc(max_mcc.precision, max_mcc.recall, "./noz_prc_curve_opt.png", title)
    plot_prcmcc(max_mcc, "./noz_prcmcc_curve_opt.png", title)
#    print "plots made"
    return max_rho, max_rank






def noz_main(litbm, score_matrix, max_rank, max_rho, fname):
    # get normalized score by using max_pho
    norm_s = norm_score(score_matrix, max_rho)
    # save the normalized scores to file
    df = norm_s.unstack().reset_index()
    df.to_csv("noz_norm_score.csv", index=False)
    # get optimized rank by using max_rank
    r = [int(max_rank.split("_")[1])]
    output_rank = get_rank(norm_s, r)
    #=================HM custom code ============================================================
    df.columns = ['DB_name', 'AD_name', 'Score']
    #adjust the columns to have a BC # column 
    df[['DB', 'DB_BC']] = df['DB_name'].str.split('_', expand=True)
    df[['AD', 'AD_BC']] = df['AD_name'].str.split('_', expand=True)#[[0,1]]
    #add interaction column to describe PPI
    df['Interaction'] = df.AD.str.cat(df.DB, sep="_")
    #now sort interactions by score, and then group into interaction PPI groups
    g_all=df.sort_values(["Score"], ascending=False).groupby(['Interaction'])
#determine numerical value for best rank and use for scoring
    g_all.nth(r).dropna(how='any').sort_values(["Score"], ascending=False).to_csv("NOZ_IS_scores_Best_Rank_max_rho.csv", index=False)#line to get the best ranked interaction score matrix based on MCC by changing .nth(x)
    
    #====================
    mcc = get_screen(output_rank, litbm)
    # store mcc to file
    mcc.to_csv(fname+".csv", index=False)
    # plot prcmcc


#============================================================================================================================================
#============================================================================================================================================
    #Plotting scripts    

def plot_prcmcc(df, name, title):

    plt.plot(df.precision.tolist(), ".", label="precision", markersize=0.8)
    plt.plot(df.recall.tolist(), ".", label="recall", markersize=0.8)
    plt.plot(df.mcc.tolist(), ".", label="mcc", markersize=0.8)
    plt.xlabel("ORF pairs ranked by IS")
    plt.legend(loc="upper right")
    plt.title(title)
    plt.savefig(name)
    plt.close()





def plot_prc(precision, recall, output_file, title):

    plt.step(recall, precision, color='b', alpha=0.6,where='post')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 105])
    plt.xlim([0.0, 105])
    plt.title(title)
    plt.savefig(output_file)
    plt.close()




#============================================================================================================================================
#============================================================================================================================================
#code to be run to analyze files using the functions above    


#load the gold stadnard set: 
yi = load_YI1(param.GOLD)

#load count matrices:
for f in os.listdir(counts): #this runs fine, load the plus and minus his combined counts files
    if not f.endswith("_combined_counts.csv"):
        continue
    fname = "./"+f         
    if "plus" in f:
        plus_his = pd.read_table(fname, sep=",", index_col=0) #plus his counts file has plus in filename, manually added
    elif "minus" in f:
        minus_his = pd.read_table(fname, sep =",", index_col=0) #minus his is the same as plus his but with minus in name


#-----DK analkysis method
output_csv, row_freq, col_freq, minus_his_freq, AD_NAMES, DB_NAMES = score_main(plus_his, minus_his, param.weights, param.floor_perc, yi)
output_csv.to_csv("DK_mcc_summary_yi1.csv", index=False)
max_weight, max_rank, max_floor = load_summary("DK_mcc_summary_yi1.csv")
maxmcc = pd.DataFrame({"max_weight": [max_weight], "max_rank": [max_rank], "max_floor": [max_floor]})
litbm = load_litbm(param.litBM13)
dk_main(litbm, max_weight, max_rank, max_floor, minus_his_freq, row_freq, col_freq, AD_NAMES, DB_NAMES, "dk_mcc_summary_litbm_full_ranks.csv")




#----------Nozomu analysis method
gold_st = load_YI1(param.GOLD)
noz_main_values, raw_scores = main(plus_his, minus_his, gold_st)   
noz_main_values.to_csv("noz_mcc_summary_yi1.csv", index=False)         
max_rho, max_rank = load_summary_NOZ("noz_mcc_summary_yi1.csv")
maxmcc["max_rho"] = [max_rho]
maxmcc["noz_max_rank"] = [max_rank]
maxmcc.to_csv("max_parameters.csv", index=False)
noz_main(litbm, raw_scores, max_rank, max_rho, "noz_mcc_summary_litbm13")



