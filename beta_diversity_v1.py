#!/usr/bin/python
import os
import sys
import pandas as pd
import argparse
from datetime import datetime
#sys.path.append("/mnt/NL200/sunhb/pipeline/irseq/v1.0/lib/")
sys.path.append("/mnt/X500/farmers/xiongdk/pipeline/irseq/")
sys.path.append("/mnt/X500/farmers/chensi/bak_sunhb/NL200/pipeline/irseq/v1.2/lib/lib")
import pypinyin,retrying,repoze.lru
import irtools

def usage():
    year, month, day = datetime.now().year, datetime.now().month, datetime.now().day
    time = "".join([str(x)for x in year, month, day])
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required = True, help = "Input sample_pired.list")
    parser.add_argument("-d", required = True, help = "Input dir of filter result *.csv, example, 3.filter")
    parser.add_argument("-s", action = "store_true", help = "stat tsv file not random_tsv file")
    parser.add_argument("-o", default = "4.personal_%s"%time, help = "Out put dir, default=4.personal_***")
    parser.add_argument("-t", default = time, help = "time of today")
    opt = parser.parse_args()
    return opt

def main():
    opt = usage()
    infile = opt.i
    indir = opt.d
    odir = opt.o
    time = opt.t
    if not os.path.exists(odir):
        os.makedirs(odir)
    assert os.path.exists(indir), "Input dir is not exists"
    assert os.path.exists(infile), "Input sample_paired.list is not exists"
    sinfo = pd.read_table(infile)
    assert sinfo.shape[1] > 1, "sample_paired.list at least two columns"
    subfix = ".tsv" if opt.s else "_random.tsv"
    overlap = []
    for _, line in sinfo.iterrows():
        sid_a = line.ID1
        sid_b = line.ID2
        fpath_a = os.path.join(indir, sid_a, "%s%s"%(sid_a, subfix))
        fpath_b = os.path.join(indir, sid_b, "%s%s"%(sid_b, subfix))
        df_a = pd.read_table(fpath_a)
        df_b = pd.read_table(fpath_b)
        mh = irtools.overlap.MHoverlap(df_a, df_b)
        mm = irtools.overlap.MMoverlap(df_a, df_b)
        sm = irtools.overlap.SMoverlap(df_a, df_b)
        rcl = irtools.overlap.RCL(df_a, df_b)
        icc = irtools.overlap.ICC(df_a, df_b)
        overlap.append([sid_a, sid_b, mh, mm, sm, rcl, icc])
    overlap = pd.DataFrame(overlap)
    overlap.columns = ["ID1","ID2","MHoverlap","MMoverlap","SMoverlap","RCL","ICC"]
    overlap.to_csv("%s/overlap_%s.csv"%(odir, time), index=False)


if __name__ == "__main__":
    main()
