import param
import os
import logging
import logging.config

#logging.config.fileConfig("/home/rothlab/rli/02_dev/08_bfg_y2h/src/logging.conf")
align_log = logging.getLogger("alignments")

def bowtie_align(fastq, ref, output):
    """
    Align r1 and r2 to reference
    Log bowtie output 
    """
    align_log.info("Starts aligning... %s", fastq)
    align_log.info("Reference: %s", ref)
    align_log.info("Output: %s", output)

    basename = os.path.basename(fastq)

    if "R1" in fastq:
        params = "-q --norc --local --very-sensitive-local -t -p 23 --reorder "
        sam_file = basename.replace('.fastq.gz','_AD_BC.sam')
    elif "R2" in fastq:
        params = "-q --nofw --local --very-sensitive-local -t -p 23 --reorder "
        sam_file = basename.replace('.fastq.gz','_DB_BC.sam')

    input_f = "-x " + ref + " -U " + fastq + " -S " + os.path.join(output, sam_file)
    
    log_f = os.path.join(output, sam_file.replace(".sam", "_bowtie.log"))
    
    command = param.BOWTIE2 + params + input_f + " 2> " + log_f

    os.system(command)

    align_log.info("Alignment finished for %s", fastq)

    return os.path.join(output, sam_file)
