
from Bio import SeqIO
import argparse
import re


############################################################			define arg type float 0 < X > 1		###############################################################

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x

##########################################		find if a position is variable		#######################
def find_variable(n):
	out={}
	message=""
	#print "nnnnnnnnnnnnnn"

	pos=""
	for name, seq in fasta.items():
		#print str(name) + " " +str(arguments.outgroup)
		if (str(name) != str(arguments.outgroup)):
			pos += seq[n]
	seq_length=len(pos)
	good_seq_length=len(re.findall('[ATCGatcg]', pos))
	ratio= float(good_seq_length)/float(seq_length)
	bad_seq_length=seq_length - good_seq_length

	if ((ratio >= arguments.missing) or (arguments.missing_absolute>= bad_seq_length)) :
		A=len(re.findall('[Aa]', pos))		
		T=len(re.findall('[Tt]', pos))
		G=len(re.findall('[Gg]', pos))
		C=len(re.findall('[Cc]', pos))
		if (A > 0):
			A=1
		if (T > 0):
			T=1
		if (C > 0):
			C=1
		if (G > 0):
			G=1
		if ((A+C+G+T) > 1):
			message = "variable"
			
			#print "due"
			#print A + T + G +C
			for name_out, seq_out in fasta.items():
				OUT[name_out] += seq_out[n]
		else:
			if A == 1:
				message="A"
			if T == 1:
				message="T"
			if C == 1:
				message="C"
			if G == 1:
				message="G"

	else:
		message="missing"

	#print out
	return (OUT,message)






parser = argparse.ArgumentParser()

parser.add_argument('INFILE',type=str,help='path to the multi fasta alignment')
parser.add_argument('-m','--missing', metavar='0-1', default='1', help='minimum percentage of NOT missing data  -m and -ma can coexist', type =restricted_float,nargs='?')
parser.add_argument('-ma','--missing_absolute', metavar='0 -  n_seq', default='0', help='max number of sequences with missing data, -m and -ma can coexist', type =int,nargs='?')
parser.add_argument('-o','--outgroup', metavar='OUT', default="very_unlikely_name",help='name of outgroup sequence not to count in the missing data but to output ',type =str, nargs='?')
parser.add_argument('-c','--count', default=False,help='count the numbero of variable / not variable / with more missing than accepted',action='store_true')




arguments = parser.parse_args()






tot={}
fasta={}
OUT={}

for record in SeqIO.parse(arguments.INFILE, "fasta"):
	length = len(record.seq)
	fasta.update({str(record.id):str(record.seq)})
	OUT.update({str(record.id):""})
n=0
count_variable=0
count_invariant=0
count_missing=0
count_invA=0
count_invT=0
count_invC=0
count_invG=0
list_pos=[]


while n < length:
	
	(OUT,message)=find_variable(n)
	n=n+1

	if message == "variable":
		count_variable =int(count_variable)+1
		list_pos.append(n)

	if message == "missing":
		count_missing=int(count_missing) +1	
	if message == "A":	
		count_invariant = int(count_invariant)+1
		count_invA = int(count_invA) + 1	
	if message == "T":	
		count_invariant = int(count_invariant)+1
		count_invT = int(count_invT) + 1
	if message == "C":	
		count_invariant = int(count_invariant)+1
		count_invC = int(count_invC) + 1
	if message == "G":	
		count_invariant = int(count_invariant)+1
		count_invG = int(count_invG) + 1



result={}

if arguments.count:

	F=open(arguments.INFILE+"_count_invariant","w")
	

	F.write (str(count_variable) + "	" + "variable site\n")
	F.write (str(count_invariant) + "	" + "invariant sites of which\n")
	F.write (str(count_invA) + "	" + "A\n")
	F.write (str(count_invT) + "	" + "T\n")
	F.write (str(count_invC) + "	" + "C\n")
	F.write (str(count_invG) + "	" + "G\n")
	F.write (str(count_missing) + "	" + "sites failed the missing data filter\n")
	F.write (str(length) + "	" + "in total\n")
	F.close()


for name_out, seq_out in OUT.items():
		print ">" + name_out
		print seq_out
	


F=open(arguments.INFILE+"_pos","w")
F.write (n)
F.close()



