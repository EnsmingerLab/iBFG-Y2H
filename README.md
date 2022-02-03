# iBFG-Y2H
iBFG-Y2H analysis code



<!DOCTYPE html>
<html>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" href="https://www.w3schools.com/w3css/4/w3.css">
<link rel='stylesheet' href='https://fonts.googleapis.com/css?family=Roboto'>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">
<meta name="keywords" content="Harley, O'Connor, Mount, Bioinformatics, BFG, BFG-Y2H, reference, assembly, genomics, proteomics, sequencing, nanopore, analysis, Ensminger, research, legionella, pneumophila">

BFG-Y2H analysis using Roth Lab Oligo Pairs
        <p>Journal and walkthrough of how to analyze effector interactions using the Roth lab BFG-Y2H staggered oligo set</p>
  <h3 id='Request'>Directory Structure</h3>
  <p>On the guru cluster login and make a directory called 08_bfg_y2h. Within this directory several sub directories are made including src, summary, and ref.in the summary directory add your summary files as outlined in the next section. Ref should be left empty, but will hold the output fasta and bowtie index files made from your DB and AD barcodes. Src should contain all of the python scripts included in this github repository: <b>ADD LINK HERE</b>.
  <h3 id='Request'>Summary Files</h3>
  <p>To analyze data you will need several summary files in tab delimited text format. For some reason the barcode summary file should be saved with the .csv file ending. These files should appear like so: 
<b>LitBM13.txt example below: </b>
<pre>
AD	DB
PSMD2	PSMC1
EXOSC5	EXOSC1
CDKN2D	CDK4
MORF4L1	MRFAP1
PSMC2	PSMD5
CDK4	CCND2
BNIP3L	BNIP3
SH3GLB2	SH3GLB1
S100B	S100A1
BIK	BCL2L2
LSM3	LSM2
MRFAP1L1	MORF4L1
S100A1	S100P
ARL2BP	ARL2
IKZF1	IKZF5
CEP55	TSG101
AP2B1	LDLRAP1
CNOT7	BTG1
CCND3	CDK6
</pre>

<b>YI_1.txt example below</b>
<pre>
AD_ORF_ID	DB_ORF_ID
PSMD2	PSMC1
EXOSC5	EXOSC1
CDKN2D	CDK4
MORF4L1	MRFAP1
PSMC2	PSMD5
CDK4	CCND2
BNIP3L	BNIP3
SH3GLB2	SH3GLB1
S100B	S100A1
BIK	BCL2L2
LSM3	LSM2
MRFAP1L1	MORF4L1
S100A1	S100P
ARL2BP	ARL2
IKZF1	IKZF5
CEP55	TSG101
AP2B1	LDLRAP1
CNOT7	BTG1
CCND3	CDK6
....
....
....
</pre>

<b>AD_Barcodes_V4.csv (which is actually in tab delim format). DB Is identical. These gene names should not have _ seperators if possible, I think it will causes issues. </b>
<pre>
Index	Plate	Well384	Group	Kiloseq Plate	Well384	Well96	Name	Locus	ORF_BC_count	ORF_Dominant%	ORF_Second/Dominant%	ORF_BC_index	DnTag_Sequence	DnTag_Dominant%	DnTag_Second/Dominant%	UpTag_Sequence	UpTag_Dominant%	UpTag_Second/Dominant%
1242	AD_Platinum_4	P22	G1	AD_Platinum_4	P22	NA	ADEV_1	ADEV_BC-1	3	100	1	1	GTGCGTTCACTAGGGAAGTTCCTGT	100	1	CGTTGAGTAAGCAAGGGCCGTATCC	100	1
1243	AD_Platinum_4	P23	G1	AD_Platinum_4	P23	NA	ADEV_2	ADEV_BC-2	3	100	1	1	TCATAGGAGTTTTGAGAGCAGGTTA	100	1	ATAAGTTGTTTACCCAGTGAATATA	100	1
1244	AD_Platinum_4	P24	G1	AD_Platinum_4	P24	NA	ADEV_3	ADEV_BC-3	3	100	1	1	CTTGGTTTGACCCTCTTTTACTCCC	100	1	TCGATGATTGACTTGCCATCGTCTT	100	1
1189	AD_Platinum_4	B13	G1	AD_Platinum_4	B13	NA	AP2B1	AP2B1_BC-1	3	100	1	1	ATATTAATGACTTCAAAAATAAAGG	100	1	AGTGATTTCAGGTTTTTATGGCTTC	100	1
1190	AD_Platinum_4	B14	G1	AD_Platinum_4	B14	NA	AP2B1	AP2B1_BC-2	3	100	1	2	TCTATTCTTTTGGTCCTTCAGTGAC	100	1	TACGAGCCGGGAGGTCAAGGCGAAA	100	1
1191	AD_Platinum_4	B15	G1	AD_Platinum_4	B15	NA	AP2B1	AP2B1_BC-3	3	100	1	3	AAGTGTTTCTCAGTAGAGGAACATA	100	1	GGGCCAAAACTTTACTAACTCTCTG	100	1
1177	AD_Platinum_4	B01	G1	AD_Platinum_4	B01	NA	ARL2BP	ARL2BP_BC-1	3	100	1	1	TTATCCTAGTAAGACAACGGAGACC	100	1	ACTCCATATAATACACTTATGACAG	100	1
1178	AD_Platinum_4	B02	G1	AD_Platinum_4	B02	NA	ARL2BP	ARL2BP_BC-2	3	100	1	2	CAATGACCGAAATCCATGTGAGTTA	100	1	TAGGATTCTTGCTAATGGTAGGATC	100	1
1179	AD_Platinum_4	B03	G1	AD_Platinum_4	B03	NA	ARL2BP	ARL2BP_BC-3	3	100	1	3	CTAAAGAATAGAGGAATTCATACCC	100	1	TATGTTACCGAAGCATTGGCAGAGG	100	1
</pre>
  </p>
Adjusting Param.py
Now you will have to tell the param.py script what your directory structure looks like by editing it with vi param.py from within the roth lab cluster.


#import numpy as np


# If you want to do alignments
# Set this to True
ALIGNMENT = True

# if you want to do read counts
# Set this to True
READ_COUNT = True


# for score optimization

# gold standard to use
GOLD="/home/rothlab/hmount/02_dev/08_bfg_y2h/summary/YI_1.txt"

# weights to test
weights = np.arange(0, 2.4, 0.2)
floor_perc = np.arange(5,25,2.5)

# in this case we test all 4 indexes

###################################

# for calculating MCC

litBM13="/home/rothlab/hmount/02_dev/08_bfg_y2h/summary/litbm_13.txt"

###################################
# summary files are used to grep gene names, group information 
# and to create fasta reference files

# in the summary files, the following columns must exist: 
# summary for AD (all the genes and group)
AD_summary = "/home/rothlab/hmount/02_dev/08_bfg_y2h/summary/AD_barcodes_V4.csv"
# summary for DB (all the genes and group)
DB_summary = "/home/rothlab/hmount/02_dev/08_bfg_y2h/summary/DB_barcodes_V4.csv"

REF_PATH = "/home/rothlab/hmount/02_dev/08_bfg_y2h/ref/"

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
</pre>

<p>After modifying param.py, you can begin read alignment and counting using the code below. First you need to generate the reference files by running the code below while within 08_bfg_y2h. This failed for me until I installed argparse on my user account in guru. If you get an error pertaining to argparse fix it like so: 
pip install --user argparse. I had to modify main.py to not isolate group identities from file name, but instead I just put both AD and DB as G1 for group1 arbitrarily.</p>
<pre>
	/home/rothlab/rli/py/bin/python2.7 ./src/main.py --pfasta /home/rothlab/hmount/02/dev/08_bfg_y2h/ref
</pre>

<p>This code aborts after the reference files are made and will throw an error because fastq files havent been supplied yet. Make sure the ref file you generate roughly resembles that in my github with the same series of AD and DB files. Once these files are created the alignment and counting code can be run. Sge.sh uses sge_sub.sh to run the main.py script.
  </br>
  <pre>
  	./sge.sh -f /home/rothlab/hmount/01_ngsdata/181219_fastq/LP/ -o /home/rothlab/hmount/02_dev/08_bfg_y2H/test
  </pre>

  <p>Once this code finishes running you will be left with a directory that contains the following files for each pair of fastq files. For me it was two directories, one for plus his, one for minus.</p>
  <pre>
*_bc_corr.png

*_combined_counts.csv
*_dntag_rawcounts.csv
*_uptag_rawcounts.csv

*_AD_BC_bowtie.log
*_AD_BC_noh.csv
*_AD_BC_sorted.sam
*_DB_BC_bowtie.log
*_DB_BC_noh.csv
*_DB_BC_sorted.sam

main.log

</pre>

The *_combined_counts.csv files are used for scoring. For Python analysis will need to combine counts for the selective (-HIS) and non-selective condition (+HIS). I do the scoring of these combined counts locally in python. Using spyder I then manually run each line of the score_HM.py or noz_score_HM.py to analyze the counts files. First change to your analysis directory which should contain your combined_counts.csv files for both selective and non-selective conditions. Within this directory should be a subdirectory with the .src files like param.py. Change into this directory and import all of the necessary modules. This is important, if you don't do this the code won't run. Also the param.py file that you use locally will be different from that you use on the cluster. Make sure to adjust it accordingly. Inside your scoring directory make a directory that has the summary files. Also make a ref directory as well. Also for combined counts files for the plus his condition you need to add plus somewhere in the file name and minus for the minus his samples.

For Rstudio local analysis, which was actually used in the manuscript you can keep the up and down tag files seperate but they will need to be modified. This requires running the cluster_output_reformat_V2.R script before running analyze_cluster_output.R R script. The count files will also be named appropriately in a format like the following: 
<pre>
1_-His_1_UpUp_QC_counts.tsv
1_-His_1_DnDn_QC_counts.tsv
1_+His_1_UpUp_QC_counts.tsv
1_+His_1_DnDn_QC_counts.tsv
</pre>

Download all four counts.csv files (-HIS up, -HIS down, +HIS up and +HIS down) samples to a directory that will be used for scoring. Then just follow along with the R script and files in the directory included. You need to hhave the working directory be set to the directory with counts files, and you need to have a seperate directory with the barcodes_V3.csv, lib directory, and res directory. See examples provided. The analyze_cluster_output.R describes which files are required. 
	



 


<div style="margin-top: 100px; font-size: 8pt;text-align:center;">
    <p><a href="http://creativecommons.org/licenses/by-sa/4.0/">CC-BY-SA</a> The Ensminger Lab, University of Toronto, 2022</p>
  </div>




</body>
</html>

<!DOCTYPE html>





