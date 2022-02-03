# Ensminger_BFG
Files for BFG analysis in the Ensminger Lab


This set of files was adapted from the BFG-Y2H analysis pipeline written by Jochen Weile in the Roth Lab. The original protocol can be found here: http://dalai.mshri.on.ca/~jweile/projects/bfg_pipeline/doc/manual.html

I have adapted several of these scripts to optimize analysis for Ensminger Lab samples using the inducible barcode fusion genetic approach. 

I will highlight some key processing steps and differences here. 


#================================================================================================================
Often BFG runs will be performed with PhiX reads present as well, or in the presence of nextera libraries which are automatically demutliplexed. If PhiX was used it will need to be seperated from the BFG reads. To do so we use the following approach using flexbar: 


flexbar -r RCP_R1.fastq -p RCP_R2.fastq -b plate_adapter_fwd.fa -b2 plate_adapter_rev.fa -bk -u 9 -n 8

#------------------------------
>plate_adapter_fwd_BC1
NNNNNNNNNCCCTTAGAACCGAGAGTGTG
>plate_adapter_fwd_BC2
NNNNNNNNNCTCCAGGGTTAGGCAGATG


>plate_adapter_rev_BC1
NNNNNNNNNGTTATCAGAGGTATGCGAGTTAG
>plate_adapter_rev_BC2
NNNNNNNNNCTTGACTGAGCGACTGAGG
#------------------------------



This produces 4 files with reads, the BC1 reads in the fw and reverse direction, the BC2 reads in the fwd and reverse direction. The BC1 fwd and BC2 fwd reads are then catted together with in the same order for both 

cat flexbarOut_barcode_plate_adapter_fwd_BC1-plate_adapter_rev_BC1_1.fastq flexbarOut_barcode_plate_adapter_fwd_BC2-plate_adapter_rev_BC2_1.fastq > BFG_R1.fastq
cat flexbarOut_barcode_plate_adapter_fwd_BC1-plate_adapter_rev_BC1_2.fastq flexbarOut_barcode_plate_adapter_fwd_BC2-plate_adapter_rev_BC2_2.fastq > BFG_R2.fastq

#================================================================================================================

Files should then be uploaded to the scinet cluster like so


scp BFG_R1.fastq.gz username@niagara.computecanada.ca:~


next login to the Roth lab cluster and download the files to the desired directory by running this command

scp mounthar@niagara.computecanada.ca:~/BFG_R1.fastq.gz . 


#================================================================================================================

If specifying muxtags for specific projects the file format for the muxtags.csv file is like the following: 

bar2num plate index combo, sample identifier, selection, sample number, date, projectID
P01-P02,LP,+His,1,N1,8-15-2018,V2
P03-P04,LP,-His,1,N2,8-15-2018,V2
P05-P06,CBU,+His,1,N3,8-15-2018,EL
P07-P08,CBU,-His,1,N4,8-15-2018,EL





#================================================================================================================
to properly begin analysis cd into the directory with your BFG reads and copy this directory like so: 


If your barcodes were generated through RCP-PCR analysis you will likely have a file with barcode sequence, gene identity and plate position. This can be used with the python code barcode_csv_maker.py to generate a barcode.csv file.  
The input file should have the following columns (X implies unimportant):
X,X,LPG,X,BC1,X,BC2,X,X,PLATE(i.e.AD_PLAT_1),WELL (i.e. B09),



With a successfully created barcodes.csv file, the barcodes.fasta file can also be generated using the barcode_fasta_maker.py script which takes the barcodes.csv file as an input.

These files come pre-generated for the BFG V2 platinum set in this repository.




git clone 

