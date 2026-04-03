
# Description: Process bedpe file and select cis_reads to Higashi format
# Author: genger
# Data: 2022-07-04-Monday


from tqdm import *
import re
import pandas as pd

import os
import sys
import getopt

## change the path 
path = "./"
os.chdir(path)

## usage
def usage():
    print("")
    print("Process bedpe file and select cis_reads to Higashi format.\nDefault set is put cis reads by order(pos).\nDefault set is throw away chr XYM!")
    print("uage: python %s -option <argument>" %sys.argv[0])
    print(" -h/--help ")
    print(" --input_stat=<STRING> stat input file.")
    print(" --input_bed=<STRING> bedpe input file.")
    print(" --output=<STRING> higashi required file.")
    print(" --stat=<STRING> stat file.")

## deal with options
try: 
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help","input_stat=","input_bed=","output=","stat="])
except getopt.GetoptError:
    print("ERROR: Get option error.\nYou can contact the author through wechat 13958598285.")
    usage()
    sys.exit(2)

for opt, val in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit(1)
    else:
        if opt in ("--input_stat"):
            input_stat_url = val 
        if opt in ("--input_bed"):
            input_bed_url = val
        if opt in ("--output"):
            output_url = val
        if opt in ("--stat"):
            stat_url = val

try:
    print("Your input_stat is:",input_stat_url)
    print("Your input_bed is:",input_bed_url)
    print("Your output is:",output_url)
    print("Your stat is:",stat_url)
except:
    print("ERROR: Missing option.")
    usage()
    sys.exit(2)

raw_stat=pd.read_csv(input_stat_url)

bed=pd.read_table(input_bed_url,header=None)
bed.columns=["chrom1","start1","end1","chrom2","start2","end2","info","quality","strand1","strand2"]
bed["pixel"]=bed["info"].apply(lambda x:re.findall("_(\w+_\w+)",x)[0])
bed["mid1"]=(bed["start1"]+bed["end1"])/2
bed["mid2"]=(bed["start2"]+bed["end2"])/2

# all valid contact
stat=bed["pixel"].value_counts()
raw_stat["valid_contact"]=raw_stat["pixel"].apply(lambda x: stat[x] if x in stat else 0)

# trans contact
stat=bed.loc[bed["chrom1"]!=bed["chrom2"]]["pixel"].value_counts()
raw_stat["trans_contact"]=raw_stat["pixel"].apply(lambda x: stat[x] if x in stat else 0)

# cis contact
stat=bed.loc[bed["chrom1"]==bed["chrom2"]]["pixel"].value_counts()
raw_stat["cis_contact"]=raw_stat["pixel"].apply(lambda x: stat[x] if x in stat else 0)

# cis >= 10kb contact
stat=bed.loc[(bed["chrom1"]==bed["chrom2"]) & (abs(bed["mid1"]-bed["mid2"])>=10000)]["pixel"].value_counts()
raw_stat["cis_more_10kb"]=raw_stat["pixel"].apply(lambda x: stat[x] if x in stat else 0)

# cis < 10kb contact
stat=bed.loc[(bed["chrom1"]==bed["chrom2"]) & (abs(bed["mid1"]-bed["mid2"])<10000)]["pixel"].value_counts()
raw_stat["cis_less_10kb"]=raw_stat["pixel"].apply(lambda x: stat[x] if x in stat else 0)

# sex contact
stat=bed.loc[(bed["chrom1"]=="chrX") | (bed["chrom1"]=="chrY") | (bed["chrom1"]=="chrM") | (bed["chrom2"]=="chrX") | (bed["chrom2"]=="chrY") | (bed["chrom2"]=="chrM")]["pixel"].value_counts()
raw_stat["sex_contact"]=raw_stat["pixel"].apply(lambda x: stat[x] if x in stat else 0)

raw_stat.to_csv(stat_url,index=None)

out=bed[["pixel","chrom1","mid1","chrom2","mid2"]].copy()
out["mid1"]=round(out["mid1"]).astype("int")
out["mid2"]=round(out["mid2"]).astype("int")
out["count"]=1
out.columns=["barcode","chrom1","pos1","chrom2","pos2","count"]
out.to_csv(output_url,sep="\t",index=None)

