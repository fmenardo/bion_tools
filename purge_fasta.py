
from Bio import SeqIO
import argparse
import re




parser = argparse.ArgumentParser()

parser.add_argument('INFILE',type=str,help='path to the multi fasta alignment')
parser.add_argument('-l','--list_p', metavar='', help='list with seq to purge (or to keep, if -k is present', type =str,nargs=1)
parser.add_argument('-k','--keep', default=False,help='output fasta with taxa in the list',action='store_true')

arguments = parser.parse_args()

fasta={}

with open(arguments.list_p[0]) as f:
    lines = f.read().splitlines()



for record in SeqIO.parse(arguments.INFILE, "fasta"):
	flag=0
	length = len(record.seq)
	for line in lines:
		if (line == str(record.id)):
			flag=1

	if arguments.keep:
		if (flag==1):
			print ">" +str(record.id)
			print str(record.seq)

	else:
		if (flag ==0):
	
			print ">" +str(record.id)
			print str(record.seq)







