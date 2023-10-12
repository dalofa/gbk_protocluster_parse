# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 14:41:29 2022

@author: David Faurdal

Updated: 12-10-2023 to include region information in addition to protoclusters.

A script for generating a csv table of protocluster predicted by antiSMASH.
The script takes a genbank-file annotated by antiSMASH and outputs a csv-table
from it.
"""
import sys

# Get cmd line inputs
argumentList = sys.argv[1:]

if len(argumentList)<2:
    print("Missing either input filename or output filename. Provide filename first and then output filename.")
    exit()
print("arguments provided", argumentList)
file_path=argumentList[0]
output_name=argumentList[1]

from Bio import SeqIO
import pandas as pd

# Collect all protoclusters and all regions from antiSMASH annotations
all_protoclusters=[]
all_regions=[]
contig_n =1

genbank = SeqIO.parse(file_path,format = "genbank")

# get informations for regions and protoclusters annotated by antiSMASH
for record in genbank:
    
    # get information from regions annotated by antiSMASH
    regions = [f for f in record.features if f.type == "region" and f.qualifiers["tool"]==['antismash']]
    for region in regions:
        r_num = "Region number " + str(contig_n)+ "." +  str(region.qualifiers["region_number"])[2:-1]
        r_product = str(region.qualifiers["product"])[2:-1].replace("'","")
        r_st=int(region.location.start)+1 #because indexing from 0 in python
        r_end=int(region.location.end)
        r_contig = record.id
        r_info=[r_num, r_product, r_st, r_end, r_contig]
        all_regions.append(r_info)
        
        # get information from protoclusters annotated by antiSMASH
    protoclusters = [f for f in record.features if f.type == "protocluster" and f.qualifiers["tool"]==['antismash']]
    for cluster in protoclusters:
        c_num = "Protocluster number " + str(contig_n)+ "." +  str(cluster.qualifiers["protocluster_number"])[2:-1]
        r_product = str(cluster.qualifiers["product"])[2:-1].replace("'","")
        r_st=int(cluster.location.start)+1 #because indexing from 0 in python
        r_end=int(cluster.location.end)
        r_contig = record.id
        r_info=[c_num, r_product, r_st, r_end, r_contig]
        all_protoclusters.append(r_info)
    contig_n+=1

# match protoclusters to regions
all_total=[]

for region in all_regions:
    r_num = region[0]
    r_product = region[1]
    r_start = region[2]
    r_stop= region[3]
    r_contig = region[4]
    
    for proto in all_protoclusters:
        c_start = proto[2]
        c_stop = proto[2]
        
        if c_start >= r_start and c_stop <= r_stop:
            new_row = proto
            new_row.append(r_num)
            new_row.append(r_product)
            new_row.append(r_start)
            new_row.append(r_stop)
            new_row.append(r_contig)
            
            all_total.append(new_row)

# construct dataframe
total_df = pd.DataFrame(all_total,columns=["protocluster","product","gbk_start","gbk_end","contig name","region","product","gbk_start","gbk_end","contig name"])

# write collected dataframe to file
total_df.to_csv(output_name)