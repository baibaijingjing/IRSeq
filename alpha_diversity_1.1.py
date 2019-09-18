#!/usr/bin/env python
#-*- coding=utf-8 -*-
import os
import sys
import pandas as pd
import numpy as np
import argparse
from datetime import datetime
sys.path.append("/mnt/X500/farmers/chensi/bak_sunhb/NL200/pipeline/irseq/v1.2/lib/lib")
import pypinyin,retrying,repoze.lru
#sys.path.append("/mnt/NL200/xiongdk/script/irseq/script/")
#sys.path.append("/mnt/NL200/sunhb/pipeline/irseq/v1.0/lib/")
sys.path.append('/mnt/X500/farmers/xiongdk/pipeline/irseq')
import irtools
def usage():
    year, month, day = datetime.now().year, datetime.now().month, datetime.now().day
    time = "".join([str(x)for x in year, month, day])
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required = True, help = "Input sample.list")
    parser.add_argument("-d", required = True, help = "Input dir of filter result *.csv, example, 3.filter")
    parser.add_argument("-s", action = "store_true", help = "stat tsv file not random_tsv file")
    parser.add_argument("-o", default = "4.personal_%s"%time, help = "Out put dir, default=4.personal_***")
    parser.add_argument("-t", default = time, help = "time of today")
    opt = parser.parse_args()
    return opt
#cloneCount      cloneFraction
def shannon(df):
    """https://en.wikipedia.org/wiki/Diversity_index#Shannon_index"""
    # total_reads = df.cloneCount.sum()
    # df.cloneFraction = df.cloneCount / total_reads
    return -(sum(df.cloneFraction * np.log(df.cloneFraction)))
def evenness(df):
    """S = H(obs)/H(max)
    Where H(obs) is the number derived from the Shannon diversity index and H(max) is the maximum possible value of H(obs) (if every species was equally likely)
    https://en.wikipedia.org/wiki/Species_evenness"""
    obs = shannon(df)
    count = df.shape[0]
    max_freq = 1.0 / count
    max_vector = np.repeat(max_freq,count)
    pre = -(sum(max_vector * np.log(max_vector)))
    return obs / pre
def clonality2(df):
    """ S = 1 − H(obs)/H(max)  ref: DOI: 10.1126/scitranslmed.3010760 """
    return 1 - evenness(df)

def main():
    opt = usage()
    inSample = opt.i
    inDir = os.path.abspath(opt.d)
    if not os.path.exists(opt.o):
        os.mkdir(opt.o)
    outFile = os.path.join(opt.o, "alpha_diversity_%s.txt"%opt.t)
    subfix = opt.s
    # outSamplePath = "out_smaple_path.list"
    # inpath = sys.argv[1].strip()
    # outpath = os.path.dirname(inpath)
    # outfile = sys.argv[2] if len(sys.argv) == 3 else os.path.join(outpath,"alpha_index.txt")
    head = ("sampleid,total_reads,max_freq,count,simpsonIndex,shannonIndex,d50Index,"
            "inverseSimpsonIndex,giniSimpsonIndex,evenness,richness,clonality").split(",")
    #outfile = os.path.join(outpath, "alpha_index.txt")
    fw = open(outFile,"w")
    fw.write("\t".join(head) + "\n")
    # fwp = open(outSamplePath, "w")
    dfs = pd.read_table(inSample)
    for ix, line in dfs.iterrows():
        sampleID = line["sampleid"]
        sampleID_1 = sampleID +".tsv"
        sampleID_2 = sampleID +"_random.tsv"
        samplePath_1 = os.path.join(os.path.join(inDir, sampleID), sampleID_1)
        samplePath_2 = os.path.join(os.path.join(inDir, sampleID), sampleID_2)
        if subfix:
            samplePath = samplePath_1
        elif os.path.exists(samplePath_2):
            samplePath = samplePath_2
        elif os.path.exists(samplePath_1):
            samplePath = samplePath_1
        else:
            print "Sample path is not exists: %s"%sampleID
            break
        # fwp.write("\t".join([sampleID, samplePath]) + "\n")
        df = pd.read_table(samplePath)
        total_reads = df.cloneCount.sum()
        df.cloneFraction = df.cloneCount / total_reads
        max_freq = df.cloneFraction.max()
        count = df.shape[0]
        simpson = irtools.diversity.simpson(df)
        shannon = irtools.diversity.shannon(df)
        d50Index = irtools.diversity.d50Index(df)
        inverseSimpsonIndex = irtools.diversity.inverseSimpson(df)
        giniSimpsonIndex = irtools.diversity.giniSimpson(df)
        evenness = irtools.diversity.evenness(df)
        richness = irtools.diversity.richness(df)
        clonality = irtools.diversity.clonality(df)
        # clonality_2 = clonality2(df)
        alpha_index = [sampleID,total_reads,max_freq,count,simpson,shannon,d50Index,inverseSimpsonIndex,giniSimpsonIndex,evenness,richness,clonality]
        alpha_index = map(str, alpha_index)
        fw.write("\t".join(alpha_index) + "\n")
    fw.close()
    # fwp.close()
if __name__ == "__main__":
    main()
