# bion_tools
Basic python and perl tools for bioinformatic 

Python tools depend on BioPhyton, and should compatible with both python 2 and 3
Perl tools depende on the Menafmod.pm which is incuded, if you place bion_tools into your home they should work straight away, otherwise you have to modify the path at the beginning of the script

**assemble.py:** cut, complement and/or reverse fasta, works on file with a single sequence, if there are multipe sequences will consider only the first

**check_length_fasta.py:** check the length of every sequence in a fasta file

**check_missing_fasta.py:** check the proportion of missing data (gaps, N, X, anything that is not actgACTG) for each sequence in a fasta file

**purge_fasta.py:** extract or delete sequences from a multi fasta file

**transpose:** transpose a file, needs two arguments, infile and separator

**convert_ali:** convert alignments between different formats (FASTA NEXUS PHYLIP)

**freq_count:** tally occurrences of lines in a file

