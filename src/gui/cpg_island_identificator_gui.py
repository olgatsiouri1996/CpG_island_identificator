# python3
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq, reverse_complement
from Bio.Data import IUPACData
import pandas as pd
# input parameters
@Gooey(required_cols=4, program_name='CpG island identificator', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="identify CpG islands on one or many sequences based on the Gardiner-Garden and Frommer (1987) method")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input single or multi-fasta file")
    ap.add_argument("-gc", "--gc",  required=False, type=int, default=50, help="min GC content(support for S and W nucleotides), integer")
    ap.add_argument("-ratio", "--ratio", required=False, type=float, default=0.6, help="min ratio of the Obs/Exp value, float")
    ap.add_argument("-step", "--step", required=True, type=int, help="step size for CpG identification, integer")
    ap.add_argument("-win", "--window", required=True, type=int, help="window size for CpG identification, integer")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output txt file")
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
    for record in SeqIO.parse(args['input'], "fasta"):
        rever = reverse_complement(record.seq)
    # forward sequence
        for i in range(0, len(record.seq) - args['window'] + 1, args['step']):
            if gc(record.seq[i:i + args['window']]) >= args['gc'] and ratio(record.seq[i:i + args['window']]) >= args['ratio']:
                gcobs.append(obs(record.seq[i:i + args['window']]))
                gccont.append(gc(record.seq[i:i + args['window']]))
                gcratio.append(ratio(record.seq[i:i + args['window']]))
                start.append(i + 1) # fix python index
                end.append(args['window']) # retrieve the end position of each putative CpG island
                headers.append(record.id)
                strand.append('+')
    # reverse complement
            if gc(rever[i:i + args['window']]) >= args['gc'] and ratio(rever[i:i + args['window']]) >= args['ratio']:
                gcobs.append(obs(rever[i:i + args['window']]))
                gccont.append(gc(rever[i:i + args['window']]))
                gcratio.append(ratio(rever[i:i + args['window']]))
                start.append(-1*(i + args['window']  - len(record.seq)))
                end.append(-1*(i - len(record.seq))) # retrieve the end position of each putative CpG island
                headers.append(record.id)
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
            df.to_csv(header = True, index = False, sep = '\t', doublequote= False, line_terminator= '\n')
        )

if __name__ == '__main__':
    main()