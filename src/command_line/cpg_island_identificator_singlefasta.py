#!/usr/bin/python3
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq, reverse_complement
from Bio.Data import IUPACData
import pandas as pd
pd.options.mode.chained_assignment = None
# input parameters
ap = argparse.ArgumentParser()
ap.add_argument("-gc", "--gc_content",  required=False, default=50, help="min GC content(support for S and W nucleotides).Default= 50")
ap.add_argument("-ratio", "--gc_ratio", required=False, default=0.6, help="min ratio of the Obs/Exp value, type = float. Default= 0.6")
ap.add_argument("-step", "--step_size", required=True, help="step size for CpG identification, type = integer")
ap.add_argument("-win", "--window_size", required=True, help="window size for CpG identification, type = integer")
args = vars(ap.parse_args())
# calculate obs value
def obs(seq):
    return seq.count('CG')
# calculate Exp value
def exp(seq):
    return round(seq.count('C') * seq.count('G') / int(args['window_size']), 2)
# calculate gc content
def gc_content(seq):
    gc = sum(seq.count(x) for x in ["G", "C", "S"])
    return round(gc * 100 / sum(seq.count(x) for x in ["A", "T", "G", "C", "S", "W"]), 2)
# main
gcobs = []
gcexp = []
headers = [] # setup empty lists
# import each fasta file from the working directory
for filename in sorted(os.listdir(str(os.getcwd()))):
    if filename.endswith(".fa") or filename.endswith(".fasta"):
        record = SeqIO.read(filename, "fasta")
        rev = reverse_complement(record.seq)
        for i in range(0, len(record.seq) - int(args['window_size']) + 1, int(args['step_size'])):
            if gc_content(record.seq[i:i + int(args['window_size'])]) > float(args['gc_content']):
                gcobs.append(obs(record.seq[i:i + int(args['window_size'])]))
                gcexp.append(exp(record.seq[i:i + int(args['window_size'])]))
                headers.append(i)
            if gc_content(rev[i:i + int(args['window_size'])]) > float(args['gc_content']):
                gcobs.append(obs(rev[i:i + int(args['window_size'])]))
                gcexp.append(exp(rev[i:i + int(args['window_size'])]))
                headers.append(i -len(record.seq))
        # create data frame
        df = pd.DataFrame()
        df['start'] = headers
        df['obs'] = gcobs
        df['exp'] = gcexp
        df['obs/exp'] = round(df['obs']/df['exp'], 2)
        df = df[df['exp'] > 0]
        df = df[df['obs/exp'] > float(args['gc_ratio'])]
        df = df.sort_values(by=['obs/exp'], ascending=False)
        df['id'] = record.id
        start = df.iloc[:,0]
        end = start.astype(int) + int(args['window_size'])
        df1 = pd.DataFrame()
        df1 = df[["id"]]
        df1['start'] = start
        df1['end'] = end
        df1['obs'] = df[["obs"]]
        df1['exp'] = df[["exp"]]
        df1['obs/exp'] = df[["obs/exp"]]
        # export
        with open(''.join([filename.split(".")[0],"_CpG_islands_","w",str(args['window_size']),"_","s",str(args['step_size']),"_","r",str(args['gc_ratio']),"_","gc",str(args['gc_content']),".txt"]), 'a') as f:
            f.write(
                df1.to_csv(header = True, index = False, sep = '\t', doublequote= False, line_terminator= '\n')
            )
        # remove lists and dataframe to be used for the next file
        gcobs.clear(); gcexp.clear(); headers.clear(); del df; del start; del end; del df1


