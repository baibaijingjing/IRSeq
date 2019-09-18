#!/usr/bin/env python
import os
import sys
import pandas as pd
from datetime import datetime
from argparse import ArgumentParser
sys.path.append("/mnt/X500/farmers/chensi/bak_sunhb/NL200/pipeline/irseq/v1.2/lib/lib")
import pypinyin,retrying,repoze.lru
sys.path.append("/mnt/X500/farmers/xiongdk/pipeline/irseq/")
import irtools
def get_args():
    parser = ArgumentParser()
    parser.add_argument("-i", help="Input sample.list", required=True)
    parser.add_argument("-d", help="Input workdir, example 2.mixcr", required=True)
    parser.add_argument("-o", help="outdir workdir, default=3.filter", default="3.filter")
    parser.add_argument("-f", help="Filter by frequence, conflict with -r, default=0.00001", type=float, default=0.00001)
    parser.add_argument("-r", help="Filter by reads, conflict with -f, recommand 2", type=int)
    parser.add_argument("-p", help="Out path.list of filter result file", default="out_path.list")
    parser.add_argument("-n", help="not plot 3D VJplot", action="store_true")
    parser.add_argument("-subfix",help="output file subfix", default=".tsv")
    parser.add_argument("-type",default='TRB',choices=['TRB', 'IGH'], help='default:TRB')
    return parser.parse_args()

def filter_freq(sub_df, min_freq):
    sub_df = sub_df[sub_df['cloneFraction'] >= min_freq]
    total_seq = sub_df['cloneCount'].sum()
    sub_df['cloneFraction'] = sub_df['cloneCount'] / total_seq
    return sub_df
def filter_reads(sub_df, min_reads):
    sub_df = sub_df[sub_df['cloneCount'] >= min_reads]
    total_seq = sub_df['cloneCount'].sum()
    sub_df['cloneFraction'] = sub_df['cloneCount'] / total_seq
    return sub_df
def main():
    opt = get_args()
    smapleList = opt.i
    barplot_3d = "/mnt/NL200/sunhb/pipeline/irseq/v1.1/lib/3dbarplot.py"
    df_sample = pd.read_table(smapleList,encoding="utf-8")
    opt.o = os.path.abspath(opt.o)
    year, month, day = datetime.now().year, datetime.now().month, datetime.now().day
    time = "".join([str(x)for x in year, month, day])
    if not os.path.exists(opt.o):
        os.makedirs(opt.o)
    outPathList = os.path.join(opt.o, opt.p)
    fw = open(outPathList, "w")
    res = []
    for ix,line in df_sample.iterrows():
        sampleID = line.sampleid
        sampleFile = os.path.abspath(opt.d) + "/" + sampleID + "/" + sampleID + opt.subfix
        outDir = os.path.join(opt.o, sampleID)
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        outFile = os.path.join(outDir, "%s.tsv"%sampleID)
        df = pd.read_table(sampleFile)
        total = df.cloneCount.sum()
        sub_df = df[~df['aaSeqCDR3'].str.contains('_|\*')]
        fltByAA = sub_df.cloneCount.sum()
        if opt.r:
            sub_df = filter_reads(sub_df, opt.r)
        else:
            sub_df = filter_freq(sub_df, opt.f)
        fw.write("\t".join([sampleID, outFile]) + "\n")
        fltByReadsOrFreq = sub_df.cloneCount.sum()
        sub_df.to_csv(outFile,index=False,sep="\t")
        res.append([sampleID, total, fltByAA, fltByReadsOrFreq])
        #code = "cd %s\n"%outDir
        #code += "python %s -indir %s -sc %s %s" %(barplot_3d, outDir, opt.type, sampleID)
        outFile = os.path.join(outDir, "%s.VJParing.pdf"%sampleID)
        if opt.n:
            continue
        irtools.plot.barplot3d(sub_df, outFile, sampleID)
        #os.system(code)
        # print "\t".join([sampleID,outFile])
    fw.close()
    res  =pd.DataFrame(res)
    res.columns = ["sampleid", "total", "fltByAA", "fltByReadsOrFreq"]
    out_stat = os.path.join(opt.o, "stat_filter_%s.txt"%time)
    res.to_csv(out_stat, sep = "\t", index=False)
if __name__ == "__main__":
    main()
