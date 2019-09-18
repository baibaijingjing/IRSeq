#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import sys
import glob
from datetime import datetime
import numpy as np
import pandas as pd
from pathlib import Path
from argparse import ArgumentParser
#sys.path.append("/mnt/X500/farmers/xiongdk/pipeline/irseq")
sys.path.append("/mnt/X500/farmers/chensi/software/untidy/geneplus/irseq")
sys.path.append("/mnt/X500/farmers/xiongdk/pipeline/irseqOnLine/lib")
sys.path.append("/mnt/X500/farmers/chensi/bak_sunhb/NL200/pipeline/irseq/v1.2/lib/lib")
import pypinyin,retrying,repoze.lru
import irtools
import indexes

def get_argv():
    parser = ArgumentParser()
    parser.add_argument("-i", required=True, help="Input sample list, at least one column of lib id, header is sampleid")
    parser.add_argument("-o", default="2.mixcr", help="Output dir, default=2.mixcr")
    opt = parser.parse_args()
    return opt

def main():
    opt = get_argv()
    infile = Path(opt.i).absolute()
    odir = Path(opt.o).absolute()
    assert infile.is_file(), "Your input samplelist is not a file, or is not exist!"
    if not odir.exists():
        odir.mkdir()
    sinfo = pd.read_table(infile)
    sid_fail = []
    sid_stat_df = []
    for sid in sinfo.sampleid:
        res = irtools.db.get_fpath(sid)
        if len(res) == 0:
            sid_fail.append(sid)
            continue
        res = res[-1]
        sdir = odir.joinpath(sid)
        if not sdir.exists():
            sdir.mkdir()
        os.system("ln -s %s %s"%(res, str(sdir)))
        sid_stat = irtools.db.get_runstat(sid)
        sid_stat_df.append(sid_stat)
    sid_stat_df = pd.concat(sid_stat_df)
    year, month, day = datetime.now().year, datetime.now().month, datetime.now().day
    time = "".join([str(x)for x in year, month, day])
    sid_stat_df.to_csv("%s/statistic_%s.csv"%(str(odir), time), index=False)
    sid_fail = "\n".join(sid_fail)
    with open("%s/sample_notin_db%s.txt"%(str(odir), time), "w") as fw:
        fw.write(sid_fail)
if __name__ == "__main__":
    main()
