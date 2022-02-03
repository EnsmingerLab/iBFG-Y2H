# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:01:24 2017

@author: Harley
"""


def rev_complementer(x):
    reverse_sequence=x[::-1]
    rev_com=''
    for i in reverse_sequence:
        if i =='G' or i=='g':
            rev_com+='C'
        elif i =='C' or i=='c':
            rev_com+='G'
        elif i =='A' or i=='a':
            rev_com+='T'
        elif i =='T' or i=='t':
            rev_com+='A'
    return rev_com


            

file=open('barcodes_944.csv', 'r')
line=file.readline().strip().split(',')
newlist=[]


while line != '':
    newlist.append(('>',line[5],'-UPTAG',line[6]))
    newlist.append(('>',line[5],'-DNTAG',line[7]))
    line=file.readline().strip().split(',')
    
    
output=open('barcodes.fa', 'w')

for i in newlist:
    output.write(i[0]+i[1]+i[2]+'\n'+i[3]+'\n')
    output.write(i[0]+'c'+i[1]+i[2]+'\n'+rev_complementer(i[3])+'\n')