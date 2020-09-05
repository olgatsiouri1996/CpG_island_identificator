# python3
import argparse
from Bio import SeqIO
import pandas as pd
# input parameters
ap = argparse.ArgumentParser()
ap.add_argument("-in", "--input_file", required=True, help="input fasta file")
ap.add_argument("-gc", "--gc_content", required=True, help="min GC content")
ap.add_argument("-ratio", "--gc_ratio", required=True, help="min ratio of the Obs/Exp value, type = float")
ap.add_argument("-win", "--window_size", required=True, help="window size for CpG identification")
ap.add_argument("-out", "--output_file", required=True, help="output txt file")
args = vars(ap.parse_args())
# calculate obs value
def obs(seq):
  return round(seq.count('CG'))
# calculate Exp value
def exp(seq):
  return round(seq.count('C') * seq.count('G') / int(args['window_size']))
# calculate gc content
def gc_content(seq):
  return round((seq.count('C') + seq.count('G')) / (seq.count('C') + seq.count('G') + seq.count('A') + seq.count('T')) * 100)
# main
gcobs = []
gcexp = []
headers = [] # setup empty lists
for record in SeqIO.parse(args['input_file'], "fasta"):
    if gc_content(record.seq) > int(args['gc_content']):
       gcobs.append(obs(record.seq))
       gcexp.append(exp(record.seq))
       headers.append(record.id)
# create data frame
df = pd.DataFrame()
df['id'] = headers
df['obs'] = gcobs
df['exp'] = gcexp
df['obs/exp'] = df['obs']/df['exp']
df = df[df['obs/exp'] > float(args['gc_ratio'])]
df = df.sort_values(by=['obs/exp'], ascending=False)
# export
with open(args['output_file'], 'a') as f:
    f.write(
        df.to_string(header = True, index = False)
    )

