# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 14:41:29 2022

@author: dalofa

A script for generating a csv table of protocluster predicted by antiSMASH.
The script takes a genbank-file annotated by antiSMASH and outputs a csv-table
from it.
"""

import sys

#Get command line inputs
argumentList = sys.argv[1:]

if len(argumentList)<2:
    print("Missing either input filename or output filename. Provide filename first and then output filename.")
    exit()
print("arguments provided", argumentList)
file_path=argumentList[0]
output_name=argumentList[1]


from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


#ollect all protoclysters from antiSMASH annotations
all_region=[]

contig_n =1
genbank = SeqIO.parse(file_path,format = "genbank")
for record in genbank:
    record.features = [f for f in record.features if f.type == "protocluster"]
    

    #rewrite to fix
    for feature in record.features:
        if feature.qualifiers["tool"]==['antismash']:
            c_num = "Protocluster number " + str(contig_n)+ "." +  str(feature.qualifiers["protocluster_number"])[2:-1]
            r_product = str(feature.qualifiers["product"])[2:-1].replace("'","")
            r_st=int(feature.location.start)+1 #because indexing from 0 in python
            r_end=int(feature.location.end)
            r_contig = record.id
            r_info=[c_num, r_product, r_st, r_end, r_contig]
            all_region.append(r_info)

    contig_n+=1

gbk_df = pd.DataFrame(all_region,columns=["protocluster","product","gbk_start","gbk_end","contig name"])       

#write gbk df to file
gbk_df.to_csv(output_name)