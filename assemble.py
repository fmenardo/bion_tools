from Bio import SeqIO
import argparse
import random
import re
from Bio.Seq import Seq


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])


def complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna])



parser = argparse.ArgumentParser()

parser.add_argument('INFILE',type=str,help='path to the fasta file')
parser.add_argument('-start','--start', metavar='', default='1', help='start position, default = 1',nargs='*')
parser.add_argument('-stop','--stop', metavar='', default='', help='stop position, default = last position',nargs='*')
parser.add_argument('-name','--name', metavar='name', default='', help='suffix for output files', type =str, required=True)
parser.add_argument('-rc','--reverse_complement', default=False,help='reverse and complement sequence output',action='store_true')
parser.add_argument('-r','--reverse', default=False,help='reverse (invert) sequence output',action='store_true')
parser.add_argument('-c','--complement', default=False,help='complement sequence output',action='store_true')


arguments = parser.parse_args()




fasta={}
out={}

for record in SeqIO.parse(arguments.INFILE, "fasta"):

        fasta.update({str(record.id):str(record.seq)})

        
for name, seq in fasta.items():
 
    length = len(seq)

    if (arguments.stop == ''):
	arguments.stop = (int(length))
    else: arguments.stop=arguments.stop[0]

       
    SEQ= seq[int(arguments.start[0])-1:int(arguments.stop)]
    #print(SEQ)
    SEQ = Seq(SEQ)
    if (arguments.reverse_complement):
        SEQ= SEQ.reverse_complement()
    if (arguments.complement):
        SEQ= SEQ.complement()
    if (arguments.reverse):
        SEQ= SEQ[::-1]


            
    out.update({str(name):str(SEQ)})
    file_out = str(arguments.name) + ".fasta"   
    with open(file_out, 'w') as f_out:
        for seq_name, sequence in out.items():
            f_out.write(">"+ str(seq_name) +"\n" + str(sequence) + "\n")

    break

