
from Bio import SeqIO
import argparse
import re

############################################################			define arg type float 0 < X > 1		###############################################################

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x


parser = argparse.ArgumentParser()

parser.add_argument('INFILE',type=str,help='path to the multi fasta alignment')
parser.add_argument('-s','--same_length', default= False, help='check if all sequences have same length: if yes print the length, if not gives warning', action='store_true')

arguments = parser.parse_args()

fasta={}

flag=0
old_length = "lala"
for record in SeqIO.parse(arguments.INFILE, "fasta"):
	length = len(record.seq)
	if (arguments.same_length):
		if (old_length != length):
			if (old_length != "lala"):	
				flag=1

	else:
		print length
	old_length=length

if (flag==1):
	print "The sequences have different lengths"
else:
	print length	
	





