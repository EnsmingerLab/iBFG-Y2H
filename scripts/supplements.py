import pandas as pd

def read_summary(AD_sum, DB_sum, AD_group="G0", DB_group="G0"):
    """
    Read from AD and DB summary files. 
    Grep gene names based on AD and DB group
    If group == G0, grep all
    """
    
    AD_summary = pd.read_table(AD_sum, sep="\t")
    DB_summary = pd.read_table(DB_sum, sep="\t")
    
    # grep group 
    if AD_group!="G0":
        if AD_group == 'GM':
            AD_summary = AD_summary[AD_summary.Plate.str.contains("Miha")]
        else:
            AD_summary = AD_summary[AD_summary.Group==AD_group]
    
    if DB_group!="G0":
        if DB_group == 'GM':
            DB_summary = DB_summary[DB_summary.Plate.str.contains("Miha")]
        else:
            DB_summary = DB_summary[DB_summary.Group==DB_group]
    
    # grep gene names
    AD_genes = AD_summary.Locus.tolist()
    DB_genes = DB_summary.Locus.tolist()

    return AD_genes, DB_genes


def parse_ds_ref(fasta):
    """
    separate dayag's fasta file
    """

    with open(fasta, "r") as ref, open("ds_AD_ref.fasta", "w") as ad, open("ds_DB_ref.fasta", "w") as db:
        c = ref.readlines()
        line =0
        while line in range(len(c)):
        #for line in c:
            print line
            if ">c" in c[line]: 
                line+=2 
                continue
            
            if ">AD" in c[line]:
                ad.write(c[line])
                ad.write(c[line+1])
            else:
                db.write(c[line])
                db.write(c[line+1])
            line+=2


def get_pair_counts(AD, DB, f):
    # get counts from combined counts file based
    # on AD and DB
    df = pd.read_csv(f, index_col=0)
    #df.set_index("0")
    #print df
    df = df.stack().reset_index()
    df.columns = ["AD", "DB","c"]
    count = df[(df.AD.str.contains(AD))& (df.DB.str.contains(DB))]
    #print df[df.AD.str.contains(DB)]
    print count
if __name__ == "__main__":
    fasta = "./ds_ref/barcodes.fasta"
    #parse_ds_ref(fasta)
    f = "/home/rothlab/rli/02_dev/08_bfg_y2h/181109_test/yAD3DB3/yAD3DB3_med_combined_counts.csv"
    get_pair_counts("YNL032W", "YNL099C", f)



