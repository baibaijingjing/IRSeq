#!/usr/bin/env python
#-*- coding=utf-8 -*-
import os
import sys
import time
import pandas as pd
import numpy as np
from collections import Counter
from argparse import ArgumentParser
sys.path.append("/mnt/X500/farmers/chensi/bak_sunhb/NL200/pipeline/irseq/v1.2/lib/lib")
import pypinyin,retrying,repoze.lru
def get_args():
    parser = ArgumentParser()
    parser.add_argument("-i",help="Input sample.list",required=True)
    parser.add_argument("-d",help="Filter sample dir, deault: 3.filter", required=True)
    parser.add_argument("-r",help="Random size,default = 1, unit is M",default=1,type=float)
    parser.add_argument("-f",help="subffix of random file, default: random", default="random")
    parser.add_argument("-p",help="Out path.list of random result file")
    parser.add_argument("-a",help="Repeat sampling if total reads less than threshold", action="store_true")
    return parser.parse_args()
def sampling(vector,size,replace=True):
    """This function is to normalization of vector
    -size   The upper limit of depths"""
    vector_range = []
    vector_return = []
    for ix,term in enumerate(vector):
        tmp = [ix] * int(term)
        vector_range.extend(tmp)
    np.random.shuffle(vector_range)
    choice = np.random.choice(vector_range, size, replace=replace)
    del vector_range
    count = Counter(choice)
    for ix,term in enumerate(vector):
        if count.get(ix,False):
            vector_return.append(count[ix])
        else:
            vector_return.append(0)
    return vector_return
def main():
    opt = get_args()
    opt.d = os.path.abspath(opt.d)
    inSampleList = os.path.abspath(opt.i)
    size = int(opt.r * 1000000)
    # outPathList = os.path.join(os.path.dirname(inSampleList),opt.p)
    suffix = opt.p if opt.p else time.strftime("%Y%m%d")
    outPathList = "out_path_random_" + suffix + ".list"
    outPathList = os.path.join(opt.d, outPathList)
    fw = open(outPathList,"w")
    print "random sampling size: %s" %size
    df = pd.read_table(inSampleList)

    for ix,line in df.iterrows():
        sampleID = line["sampleid"]
        sPath = "/".join([opt.d,sampleID,"%s.tsv"%sampleID])
        prefix,subfix = os.path.splitext(sPath)
        df = pd.read_table(sPath)
        outFile = prefix + "_" + opt.f + subfix
        if df["cloneCount"].sum() < size and not opt.a:
            os.system("ln -s %s %s"%(sPath, outFile))
            fw.write("\t".join([sampleID, outFile]) + "\n")
            continue
        random_sample = sampling(df["cloneCount"],size)
        df["cloneCount"] = random_sample
        df["cloneFraction"] = df["cloneCount"] / df["cloneCount"].sum()
        df = df[df["cloneCount"] > 0]
        # outfile = sPath
        df.to_csv(outFile,sep="\t",index=None)
        fw.write("\t".join([sampleID,outFile]) + "\n")
if __name__ == "__main__":
    main()
