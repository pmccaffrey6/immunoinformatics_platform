#!/usr/bin/python

import sys
from Bio Import SeqIO
import os



multifasta_inpath = sys.argv[1]
split_fasta_outfolder = sys.argv[2]


for record in SeqIO.parse(multifasta_inpath, "fasta"):
    with open(os.path.join(f"{split_fasta_outfolder}", f"{record.id}.fasta"), 'w') as splitfile:
        SeqIO.write(record, splitfile, "fasta")
