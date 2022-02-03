import numpy as np

###################################

# If you want to do alignments
# Set this to True
ALIGNMENT = True

# if you want to do read counts
# Set this to True
READ_COUNT = True

##################################

# for score optimization

# gold standard to use
GOLD="/home/rothlab/rli/02_dev/08_bfg_y2h/summary/YI_1.txt"

# weights to test
weights = np.arange(0, 2.4, 0.2)
floor_perc = np.arange(5,25,2.5)

# in this case we test all 4 indexes

###################################

# for calculating MCC

litBM13="/home/rothlab/rli/02_dev/08_bfg_y2h/summary/litbm_13.txt"

###################################
# summary files are used to grep gene names, group information 
# and to create fasta reference files

# in the summary files, the following columns must exist: 
# summary for AD (all the genes and group)
AD_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180627_byORFeome_AD.csv"
# summary for DB (all the genes and group)
DB_summary = "/home/rothlab/rli/02_dev/08_bfg_y2h/summary/20180627_byORFeome_DB_AA.csv"

REF_PATH = "/home/rothlab/rli/02_dev/08_bfg_y2h/ref/"

###################################

# program path

BOWTIE2 = "/home/rothlab/rli/lib/bowtie2-2.3.4.1-linux-x86_64/bowtie2 "
BOWTIE2_BUILD = "/home/rothlab/rli/lib/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build "
SAMTOOLS = "/home/rothlab/rli/lib/samtools-1.4.1/samtools "

###################################

# Padding sequences used 

# DB down tags
DB_Dn1 = "TCGATAGGTGCGTGTGAAGG"
DB_Dn2 = "CCTCAGTCGCTCAGTCAAG"
# DB up tags
DB_Up1 = "CCATACGAGCACATTACGGG"
DB_Up2 = "CTAACTCGCATACCTCTGATAAC"

# AD down tags
AD_Dn1 = "CTCCAGGGTTAGGCAGATG"
AD_Dn2 = "CAATCGCACTATCCCGCTG"
# AD up tags
AD_Up1 = "CCCTTAGAACCGAGAGTGTG"
AD_Up2 = "CACTCCGTTCGTCACTCAATAA"

###################################
