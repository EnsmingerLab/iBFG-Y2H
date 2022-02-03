#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 13:38:38 2018

@author: harleymount
""" 




#-----------------------------------AD PRS CODE----------------------------------------------------------------------------------------------------------------
sample=open('PLAT_AD_PRS_ONLY.csv', 'r')
line=sample.readline()
line=sample.readline()



prev_lpg=''
count=1
global_count=0

import csv 



with open('barcodes_AD_PRS.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        #writer.writerow(x)   
        while line != '':
            split=line.strip().split(',')
            current_lpg=split[2]
            if prev_lpg != current_lpg:      
                global_count+=1
                count=1
                writer.writerow(['AD', 'PRS', 'ADPRS'+str(global_count), split[2],'BC'+str(count), split[9].split('_')[0]+split[9].split('_')[1]+split[9].split('_')[2]+split[10], split[4], split[6]])
                prev_lpg=current_lpg
            elif prev_lpg == current_lpg:
                count+=1
                writer.writerow(['AD', 'PRS', 'ADPRS'+str(global_count), split[2],'BC'+str(count), split[9].split('_')[0]+split[9].split('_')[1]+split[9].split('_')[2]+split[10], split[4], split[6]])            
            line=sample.readline()

#--------------------------------------------------------------------------------------------------------------------------------------------
            
            
            
            
            
#-----------------------------------DB PRS CODE----------------------------------------------------------------------------------------------------------------
sample=open('PLAT_DB_PRS_ONLY.csv', 'r')
line=sample.readline()
line=sample.readline()



prev_lpg=''
count=1
global_count=0

import csv 



with open('barcodes_DB_PRS.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        #writer.writerow(x)   
        while line != '':
            split=line.strip().split(',')
            current_lpg=split[2]
            if prev_lpg != current_lpg:      
                global_count+=1
                count=1
                writer.writerow(['DB', 'PRS', 'DBPRS'+str(global_count), split[2],'BC'+str(count), split[9].split('_')[0]+split[9].split('_')[1]+split[9].split('_')[2]+split[10], split[4], split[6]])
                prev_lpg=current_lpg
            elif prev_lpg == current_lpg:
                count+=1
                writer.writerow(['DB', 'PRS', 'DBPRS'+str(global_count), split[2],'BC'+str(count), split[9].split('_')[0]+split[9].split('_')[1]+split[9].split('_')[2]+split[10], split[4], split[6]])            
            line=sample.readline()

#--------------------------------------------------------------------------------------------------------------------------------------------
            
            
            
#-----------------------------------DB LPG CODE----------------------------------------------------------------------------------------------------------------
sample=open('PLAT_DB_LPG_ONLY.csv', 'r')
line=sample.readline()
line=sample.readline()



prev_lpg=''
count=1
global_count=0

import csv 



with open('barcodes_DB_LPG.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        #writer.writerow(x)   
        while line != '':
            split=line.strip().split(',')
            current_lpg=split[2]
            if prev_lpg != current_lpg:      
                global_count+=1
                count=1
                writer.writerow(['DB', 'LPG', 'DBLPG'+str(global_count), split[2],'BC'+str(count), split[9].split('_')[0]+split[9].split('_')[1]+split[9].split('_')[2]+split[10], split[4], split[6]])
                prev_lpg=current_lpg
            elif prev_lpg == current_lpg:
                count+=1
                writer.writerow(['DB', 'LPG', 'DBLPG'+str(global_count), split[2],'BC'+str(count), split[9].split('_')[0]+split[9].split('_')[1]+split[9].split('_')[2]+split[10], split[4], split[6]])            
            line=sample.readline()

#--------------------------------------------------------------------------------------------------------------------------------------------
            
            
            
            
#-----------------------------------DB LPG CODE----------------------------------------------------------------------------------------------------------------
sample=open('PLAT_AD_LPG_ONLY.csv', 'r')
line=sample.readline()
line=sample.readline()



prev_lpg=''
count=1
global_count=0

import csv 



with open('barcodes_AD_LPG.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        #writer.writerow(x)   
        while line != '':
            split=line.strip().split(',')
            current_lpg=split[2]
            if prev_lpg != current_lpg:      
                global_count+=1
                count=1
                writer.writerow(['AD', 'LPG', 'ADLPG'+str(global_count), split[2],'BC'+str(count), split[9].split('_')[0]+split[9].split('_')[1]+split[9].split('_')[2]+split[10], split[4], split[6]])
                prev_lpg=current_lpg
            elif prev_lpg == current_lpg:
                count+=1
                writer.writerow(['AD', 'LPG', 'ADLPG'+str(global_count), split[2],'BC'+str(count), split[9].split('_')[0]+split[9].split('_')[1]+split[9].split('_')[2]+split[10], split[4], split[6]])            
            line=sample.readline()

#--------------------------------------------------------------------------------------------------------------------------------------------