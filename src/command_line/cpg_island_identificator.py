#!/usr/bin/python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq, reverse_complement
from Bio.Data import IUPACData
import pandas as pd
# input parameters
ap = argparse.ArgumentParser(description="identify CpG islands on one or many sequences based on the Gardiner-Garden and Frommer (1987) method")
ap.add_argument("-in", "--input", required=True, help="input single or multi-fasta file")
ap.add_argument("-gc", "--gc",  required=False, default=50, help="min GC content(support for S and W nucleotides). Default is 50")
ap.add_argument("-ratio", "--ratio", required=False, default=0.6, help="min ratio of the Obs/Exp value, type = float. Default is 0.6")
ap.add_argument("-step", "--step", required=True, help="step size for CpG identification, type = integer")
ap.add_argument("-win", "--window", required=True, help="window size for CpG identification, type = integer")
ap.add_argument("-out", "--output", required=True,  help="output txt file")
args = vars(ap.parse_args())
# calculate obs value
def obs(seq):
    return seq.count('CG')
# calculate Exp value
def exp(seq):
    return round(seq.count('C') * seq.count('G') / int(args['window']), 2)
# calculate gc content
def gc(seq):
    gc = sum(seq.count(x) for x in ["G", "C", "S"])
    return round(gc * 100 / sum(seq.count(x) for x in ["A", "T", "G", "C", "S", "W"]), 2)
# main
gcobs = []
gcexp = []
start = []
end = []
headers = [] # setup empty lists
# import multi-fasta
for record in SeqIO.parse(args['input'], "fasta"):
    rever = reverse_complement(record.seq)
# forward sequence
    for i in range(0, len(record.seq) - int(args['window']) + 1, int(args['step'])):
        if gc(record.seq[i:i + int(args['window'])]) > float(args['gc']):
            gcobs.append(obs(record.seq[i:i + int(args['window'])]))
            gcexp.append(exp(record.seq[i:i + int(args['window'])]))
            start.append(i + 1) # fix python index
            end.append(i + int(args['window'])) # retrieve the end position of each putative CpG island
            headers.append(record.id)
# reverse complement
        if gc(rever[i:i + int(args['window'])]) > float(args['gc']):
            gcobs.append(obs(rever[i:i + int(args['window'])]))
            gcexp.append(exp(rever[i:i + int(args['window'])]))
            start.append(i -len(record.seq))
            end.append(i -len(record.seq) + int(args['window'])) # retrieve the end position of each putative CpG island
            headers.append(record.id)
# create data frame
df = pd.DataFrame()
df['id'] = headers
df['start'] = start
df['end'] = end
df['obs'] = gcobs
df['exp'] = gcexp
df['obs/exp'] = round(df['obs']/df['exp'], 2) # calculate obs/exp ratio and filter by a user spesified threshold
df = df[df['exp'] > 0]
df = df[df['obs/exp'] > float(args['ratio'])]
df = df.sort_values(by=['obs/exp'], ascending=False) # sort by biggest obs/exp ratio
# export
with open(args['output'], 'a') as f:
    f.write(
        df.to_csv(header = True, index = False, sep = '\t', doublequote= False, line_terminator= '\n')
    )
