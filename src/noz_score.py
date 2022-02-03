from __future__ import division
import pandas as pd
import numpy as np
import os
import datetime
import math
import plot
import param
import score
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


def get_score(pre_freq, med_freq, high_freq):
    
    """
    calculate raw score
    """
    # convert pre_freq to dataframe
    pre_freq = pd.DataFrame(data = pre_freq, columns=med_freq.columns.tolist(), index=med_freq.index.tolist())

    s = (med_freq + high_freq)/pre_freq
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


def get_rank(norm_s, rank=range(0,4)):
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
    prcmcc = pd.DataFrame({}, columns=["precision","recall","mcc","rank"])
    AD_GOLD = gold_st.AD.tolist()
    DB_GOLD = gold_st.DB.tolist()
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

        MAXMCC = score.prcmcc(network_orfs, 1000)
        MAXMCC["rank"] = key
        prcmcc = prcmcc.append(MAXMCC)
    
    prcmcc = prcmcc.reset_index(drop=True)
    return prcmcc

def main(GFP_pre, GFP_med, GFP_high, gold_st):

    #calculate GFP_pre freq
    row_freq, col_freq = pre_freq(GFP_pre)
    GFP_pre_freq = np.outer(row_freq, col_freq)
    med_freq = freq(GFP_med)
    high_freq = freq(GFP_high)

    # test plot
    # should be very little corr
    # freq_corr(high_freq, med_freq)
    
    # get raw scores
    raw_s = get_score(GFP_pre_freq, med_freq, high_freq)
    df = raw_s.copy()
    df = df.unstack().reset_index()
    df.to_csv("noz_raw_score.csv", index=False)
    #test = s.values.flatten()
    #test.sort()
    #plot_diff(test)
    
    percentile = np.arange(0.1, 0.8, 0.05)
    mcc_summary = pd.DataFrame({}, columns=["precision","recall","mcc","rank","rho"])
    # test different rho to optimize mcc
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

if __name__ == "__main__":
    # test on yAD4 DB1
    test_dir = "/home/rothlab/rli/02_dev/08_bfg_y2h/181109_test/yAD1DB3/"
    
    os.chdir(test_dir)
    for f in os.listdir(test_dir):
        if not f.endswith("_combined_counts.csv"):
            continue
        fname = "./"+f
        if "pre" in f:
            GFP_pre = pd.read_table(fname, sep=",", index_col=0)
        elif "med" in f:
            GFP_med = pd.read_table(fname, sep =",", index_col=0)
        elif "high" in f:
            GFP_high = pd.read_table(fname, sep =",", index_col=0)
    
    gold_st = load_YI1(param.GOLD)
    mcc_summary, s = main(GFP_pre, GFP_med, GFP_high, gold_st)            
    mcc_summary.to_csv("noz_mcc_summary_opt.csv", index=False)
    ho, max_rank = load_summary("noz_mcc_summary_opt.csv")
