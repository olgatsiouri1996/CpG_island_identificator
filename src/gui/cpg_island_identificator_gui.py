#!/usr/bin/python3
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq, reverse_complement
from Bio.Data import IUPACData
import pandas as pd
pd.options.mode.chained_assignment = None
# input parameters
@Gooey(required_cols=4, program_name='CpG island identificator', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="identify CpG islands on a sequence based on the Gardiner-Garden and Frommer (1987) method")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input single-fasta file")
    ap.add_argument("-gc", "--gc",  required=False, default=50, help="min GC content(support for S and W nucleotides).Default= 50")
    ap.add_argument("-ratio", "--ratio", required=False, default=0.6, help="min ratio of the Obs/Exp value, type = float. Default= 0.6")
    ap.add_argument("-step", "--step", required=True, help="step size for CpG identification, type = integer")
    ap.add_argument("-win", "--window", required=True, help="window size for CpG identification, type = integer")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output txt file")
    args = vars(ap.parse_args())
    # calculate obs value
    def obs(seq):
        return seq.count('CG')
    # calculate Exp value
    def exp(seq):
        return round(seq.count('C') * seq.count('G') / int(args['window']), 2)
    # calculate gc content
    def gc_content(seq):
        gc = sum(seq.count(x) for x in ["G", "C", "S"])
        return round(gc * 100 / sum(seq.count(x) for x in ["A", "T", "G", "C", "S", "W"]), 2)
    # main
    gcobs = []
    gcexp = []
    headers = [] # setup empty lists
    record = SeqIO.read(args['input'], "fasta")
    rev = reverse_complement(record.seq)
    for i in range(0, len(record.seq) - int(args['window']) + 1, int(args['step'])):
        if gc_content(record.seq[i:i + int(args['window'])]) > float(args['gc']):
            gcobs.append(obs(record.seq[i:i + int(args['window'])]))
            gcexp.append(exp(record.seq[i:i + int(args['window'])]))
            headers.append(i)
        if gc_content(rev[i:i + int(args['window'])]) > float(args['gc']):
            gcobs.append(obs(rev[i:i + int(args['window'])]))
            gcexp.append(exp(rev[i:i + int(args['window'])]))
            headers.append(i -len(record.seq))
    # create data frame
    df = pd.DataFrame()
    df['start'] = headers
    df['obs'] = gcobs
    df['exp'] = gcexp
    df['obs/exp'] = round(df['obs']/df['exp'], 2)
    df = df[df['exp'] > 0]
    df = df[df['obs/exp'] > float(args['ratio'])]
    df = df.sort_values(by=['obs/exp'], ascending=False)
    df['id'] = record.id
    start = df.iloc[:,0]
    end = start.astype(int) + int(args['window'])
    df1 = pd.DataFrame()
    df1 = df[["id"]]
    df1['start'] = start
    df1['end'] = end
    df1['obs'] = df[["obs"]]
    df1['exp'] = df[["exp"]]
    df1['obs/exp'] = df[["obs/exp"]]
    # export
    with open(args['output'], 'a') as f:
        f.write(
            df1.to_csv(header = True, index = False, sep = '\t', doublequote= False, line_terminator= '\n')
        )

if __name__ == '__main__':
    main()
