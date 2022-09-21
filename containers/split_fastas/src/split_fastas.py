#!/usr/bin/python3

import sys
from Bio import SeqIO
import os



multifasta_inpath = sys.argv[1]
split_fasta_outfolder = sys.argv[2]


for record in SeqIO.parse(multifasta_inpath, "fasta"):
    outfile_name = record.id + '.fasta'
    with open(os.path.join(split_fasta_outfolder, outfile_name), 'w') as splitfile:
        SeqIO.write(record, splitfile, "fasta")
