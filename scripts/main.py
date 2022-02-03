import argparse
import logging
import logging.config
import os
import glob
import re
import param
import create_fasta
import alignment
import supplements
import read_counts
import plot

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='BFG-Y2H')
    
    # make fasta from summary file (AD and DB)
    # -- create: creates fasta file from summary.csv
    # -- build: if the fasta files already exist, build index files for them
    parser.add_argument('--pfasta', help="Path to fasta file")
#    parser.add_argument("--adgroup", help="AD group number", required=True)
#    parser.add_argument("--dbgroup", help="DB group number", required=True)

    # parameters for cluster
    parser.add_argument("--fastq", help="Path to all fastq files you want to analyze")
    parser.add_argument("--output", help="Output path for sam files")
    
    # for analysis
    parser.add_argument("--r1", help=".sam file for read one")
    parser.add_argument("--r2", help=".sam file for read two")

    args = parser.parse_args()
    
    # required arguments
    AD_GROUP = 'G1'
    DB_GROUP = 'G1'

    #AD_REF = REF_PATH+"y_AD_"+AD_GROUP
    #DB_REF = REF_PATH+"y_DB_"+DB_GROUP

    # processing fasta file
    fasta_output = args.pfasta
         
    if fasta_output is not None:
        # example of create fasta for AD1DB4
        create_fasta.create_fasta(param.AD_summary, param.DB_summary, fasta_output, group_spec=True, AD=AD_GROUP, DB=DB_GROUP)#adjusted code here to add param.AD_Summary and param.DB_Summary

        # example of create fasta for all 
#        create_fasta(AD_summary, DB_summary, fasta_output)

        list_fasta = os.listdir(fasta_output)

        for fasta in list_fasta:
            create_fasta.build_index(os.path.join(fasta_output,fasta), fasta_output)
        
    ##########################################################
    ###################### Alignment #########################

    output = args.output

    # input fastq is always R1
    ad = args.fastq

    # grep dir name and basename
    dir_name = os.path.dirname(ad)
    basename = os.path.basename(ad)

    # list of files in that dir
    list_fastq = os.listdir(dir_name)
        
    # get base name (sample name)
    ad_base = basename.split("_R1")[0]
        
    # find corresponding R2
    db = [i for i in list_fastq if "R2" in i and i.split("_R2")[0]==ad_base][0]
    db = os.path.join(dir_name, db)

    db_base = db.split("_R2")[0]

    m = re.match(r"yAD([1-9]|M)DB([1-9]|M)", ad_base)
    
    #print AD_GROUP
    #print DB_GROUP

    AD_REF = param.REF_PATH+"y_AD_"+AD_GROUP
    DB_REF = param.REF_PATH+"y_DB_"+DB_GROUP


    if output is None: 
        exit(0)
    else:
        output_dir = os.path.join(output, ad_base.split("_")[0]+"/")
        
        if not os.path.isdir(output_dir):
            os.system("mkdir -p "+output_dir)
        
    output = output_dir
        
    os.chdir(output)
        
    logging.config.fileConfig("/home/rothlab/hmount/02_dev/08_bfg_y2h/src/logging.conf", disable_existing_loggers=False)  
    log = logging.getLogger("root")

    if param.ALIGNMENT:
        r1_sam = alignment.bowtie_align(ad, AD_REF, output)
        r2_sam = alignment.bowtie_align(db, DB_REF, output)
        
        # check if sam files exist
        if not os.path.isfile(r1_sam) or not os.path.isfile(r2_sam):
            # log error
            log.error("SAM FILE DOES NOT EXIST - CHECK ALIGNMENTS")
            exit(0)

        dir_name = os.path.dirname(r1_sam)
        r1_basename = os.path.basename(r1_sam)
        r2_basename = os.path.basename(r2_sam)

        # sort r1_sam
        log.info("Sorting sam files..")
        sorted_r1 = os.path.join(dir_name, r1_basename.replace(".sam", "_sorted.sam"))
        sort_r1 = param.SAMTOOLS+"sort -n -o "+sorted_r1+" "+r1_sam
        # sort r2_sam
        sorted_r2 = os.path.join(dir_name, r2_basename.replace(".sam", "_sorted.sam"))
        sort_r2 = param.SAMTOOLS+"sort -n -o "+sorted_r2+" "+r2_sam
        
        # remove headers
        r1 = os.path.join(dir_name, r1_basename.replace(".sam", "_noh.sam")) 
        r2 = os.path.join(dir_name, r2_basename.replace(".sam", "_noh.sam")) 
        
        os.system(sort_r1)
        os.system(sort_r2)

        # remove original sam file
        os.system("rm "+r1_sam)
        os.system("rm "+r2_sam)

        log.info("Sam files are sorted")

        log.info("Removing header..")
        os.system("grep -v \"^@\" "+sorted_r1+" > "+r1)
        os.system("grep -v \"^@\" "+sorted_r2+" > "+r2)
        
        r1_csv = os.path.join(dir_name, r1.replace(".sam", ".csv"))
        r2_csv = os.path.join(dir_name, r2.replace(".sam", ".csv"))

        log.info("Cut file....")
        os.system("cut -f 1-5 "+r1+" > "+ r1_csv)
        os.system("cut -f 1-5 "+r2+" > "+ r2_csv)

        # remove no header sam file
        os.system("rm "+r1)
        os.system("rm "+r2)
        log.info("File shrinked")
    
    else:
        for f in os.listdir(output_dir):
            if ad_base in f:
                if "R1" in f and ".csv" in f:
                    r1_csv = f
                if "R2" in f and ".csv" in f:
                    r2_csv = f

    if param.READ_COUNT:
        
        log.info("Counting reads for %s, %s", r1_csv, r2_csv)

        AD_genes, DB_genes = supplements.read_summary(param.AD_summary, param.DB_summary, AD_group=AD_GROUP, DB_group=DB_GROUP)
        up_matrix, dn_matrix = read_counts.RCmain(r1_csv, r2_csv, AD_genes, DB_genes)
        
        uptag_file = "./"+ad_base+"_uptag_rawcounts.csv"
        dntag_file = "./"+ad_base+"_dntag_rawcounts.csv"
        
        dn_matrix.to_csv(dntag_file)
        up_matrix.to_csv(uptag_file)

        combined = up_matrix + dn_matrix

        combined.to_csv("./"+ad_base+"_combined_counts.csv")
        
        basename = r1_csv.split("R1")[0]
        # plot up and dn corr
        # sample_bc_corr.png
        plot.bc_corr(basename, up_matrix, dn_matrix)

        log.info("Barcode counts corr plot")
