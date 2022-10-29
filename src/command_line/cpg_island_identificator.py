# python3
import argparse
from pyfaidx import Fasta
import pandas as pd
# input parameters
ap = argparse.ArgumentParser(description="identify CpG islands on one or many sequences based on the Gardiner-Garden and Frommer (1987) method")
ap.add_argument("-in", "--input", required=True, help="input single or multi-fasta file")
ap.add_argument("-gc", "--gc",  required=False, type=int, default=50, help="min GC content(support for S and W nucleotides), integer. Default is 50")
ap.add_argument("-ratio", "--ratio", required=False, type=float, default=0.6, help="min ratio of the Obs/Exp value, float. Default is 0.6")
ap.add_argument("-step", "--step", required=False, type=int, default=50, help="step size for CpG identification, integer. Default is 50")
ap.add_argument("-win", "--window", required=False, type=int,default=200, help="window size for CpG identification, integer. Default is 200")
ap.add_argument("-out", "--output", required=True, help="output txt file")
args = vars(ap.parse_args())
# calculate obs value
def obs(seq):
    return seq.count('CG')
# calculate Exp value
def ratio(seq):
    obs = seq.count('CG')
    exp = seq.count('C') * seq.count('G') / int(args['window'])
    return round(obs/exp, 2)
# calculate gc content
def gc(seq):
    gc = sum(seq.count(x) for x in ["G", "C", "S"])
    return round(gc * 100 / sum(seq.count(x) for x in ["A", "T", "G", "C", "S", "W"]), 2)
# main
gcobs = []
gccont = []
start = []
end = []
strand = []
headers = []
gcratio = [] # setup empty lists
# import multi-fasta
features = Fasta(args['input'])
for key in features.keys():
# forward sequence
    for i in range(0, len(features[key]) - args['window'] + 1, args['step']):
        if gc(features[key][i:i + args['window']].seq) >= args['gc'] and ratio(features[key][i:i + args['window']].seq) >= args['ratio']:
            gcobs.append(obs(features[key][i:i + args['window']].seq))
            gccont.append(gc(features[key][i:i + args['window']].seq))
            gcratio.append(ratio(features[key][i:i + args['window']].seq))
            start.append(i + 1) # fix python index
            end.append(i + args['window']) # retrieve the end position of each putative CpG island
            headers.append(key)
            strand.append('+')
# reverse complement
        if gc(features[key][i:i + args['window']].reverse.complement.seq) >= args['gc'] and ratio(features[key][i:i + args['window']].reverse.complement.seq) >= args['ratio']:
            gcobs.append(obs(features[key][i:i + args['window']].reverse.complement.seq))
            gccont.append(gc(features[key][i:i + args['window']].reverse.complement.seq))
            gcratio.append(ratio(features[key][i:i + args['window']].reverse.complement.seq))
            start.append(-1*(i + args['window'] - features[key][:].end))
            end.append(-1*(i - features[key][:].end)) # retrieve the end position of each putative CpG island
            headers.append(key)   
            strand.append('-')
# create data frame
df = pd.DataFrame()
df['id'] = headers
df['start'] = start
df['end'] = end
df['strand'] = strand
df[''.join(['%','GC'])] = gccont
df['obs'] = gcobs
df['obs/exp'] = gcratio
df = df.sort_values(by=['obs/exp'], ascending=False) # sort by biggest obs/exp ratio
# export
with open(args['output'], 'a') as f:
    f.write(
        df.to_csv(header = True, index = False, sep = '\t', doublequote= False, lineterminator= '\n')
    )
