
from Bio import SeqIO
import argparse
import re








parser = argparse.ArgumentParser()

parser.add_argument('INFILE',type=str,help='path to the multi fasta alignment')






arguments = parser.parse_args()
tot={}
fasta={}
OUT={}

for record in SeqIO.parse(arguments.INFILE, "fasta"):
	length = len(str(record.seq))
	good_seq_length=len(re.findall('[ATCGatcg]', str(record.seq)))
	ratio= float(good_seq_length)/float(length)
	print str(record.id) + " " +str(ratio)









