# bion_tools
### Python and perl tools for bioinformatic 

Python tools depend on BioPhyton, and should compatible with both python 2 and 3.

Perl tools depend on the module menafmod.pm which is included, if you place bion_tools into your home they should work straight away, otherwise you have to modify the path at the beginning of each script.

**base_env.yml** file to creat virtual environment where all tools should work, useful also to generate custom kernels for jupyter notebooks.

Short instructions:

if you never used anaconda create a .condarc file in your home and edit to determine the default location of conda environments and pkgs eg.

    envs_dirs:
      - /data/fmenar/conda/envs
    pkgs_dirs:
      - /data/fmenar/conda/pkgs

Than:

    conda env create -f base_env.yml
    source activate base_env

if you want to add your environment to the kernel list


    ipython kernel install --user --name jupyter_nb_base
    conda deactivate


**assemble.py:** cut, complement and/or reverse sequence from a fasta, works on files with a single sequence, if there are multiple sequences will consider only the first

**check_length_fasta.py:** check the length of every sequence in a fasta file

**check_missing_fasta.py:** check the proportion of missing data (gaps, N, X, anything that is not actgACTG) for each sequence in a fasta file

**purge_fasta.py:** extract or delete sequences from a multi fasta file by name

**transpose:** transpose a file, needs two arguments, infile and separator

**convert_ali:** convert alignments between different formats (FASTA NEXUS PHYLIP)

**freq_count:** tally occurrences of lines in a file

