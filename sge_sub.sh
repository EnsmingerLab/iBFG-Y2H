#!/bin/bash
#$ -cwd
#$ -V

/home/rothlab/rli/py/bin/python2.7 ./src/main.py --fastq $1 --output $2 
