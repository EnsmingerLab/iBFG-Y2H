import score
import noz_score
import param
import itertools
import pandas as pd
import numpy as np

def load_litbm(lit):
    """
    load litbm to pandas dataframe
    input cols: AD DB
    output cols: AD DB Interactions
    """
    litbm = pd.read_table(lit)
    litbm["Interactions"] = litbm[['AD', 'DB']].apply(lambda x: '_'.join(x), axis=1)
    return litbm

def noz_main(litbm, score_matrix, max_rank, max_rho, fname):
    # get normalized score by using max_pho
    norm_s = noz_score.norm_score(score_matrix, max_rho)
    # save the normalized scores to file
    df = norm_s.unstack().reset_index()
    df.to_csv("noz_norm_score.csv", index=False)
    # get optimized rank by using max_rank
    r = [int(max_rank.split("_")[1])]
    output_rank = noz_score.get_rank(norm_s, r)
    mcc = noz_score.get_screen(output_rank, litbm)
    # store mcc to file
    mcc.to_csv(fname+".csv", index=False)
    # plot prcmcc
    

def dk_main(litbm, max_weight, max_rank, max_floor, high_freq, med_freq, row_freq, col_freq, AD_NAMES, DB_NAMES, fname):
    
    shape = med_freq.shape
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
    
    IS_normed = score.get_norm_score(max_weight, high_freq, med_freq, pre_freq)
    df = IS_normed.unstack().reset_index()
    df.to_csv("DK_norm_score.csv", index=False)
    
    r = [int(max_rank.split("_")[1])] 
    dicts = score.test_rank(IS_normed, all_pairs, r)
    mcc = score.get_mcc(dicts, litbm, AD_freq, DB_freq, row_cut, col_cut)
    mcc.to_csv(fname+".csv", index=False)

if __name__ == "__main__":
    
    test_dir = "/home/rothlab/rli/02_dev/08_bfg_y2h/rerun_analysis/yAD1DB1/"
    litbm=load_litbm(param.litBM13)
    print litbm
    #noz_main()

