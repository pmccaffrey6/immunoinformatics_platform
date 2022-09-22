#!/usr/bin/python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq
import os



multifasta_inpath = sys.argv[1]
split_fasta_outfolder = sys.argv[2]


for record in SeqIO.parse(multifasta_inpath, "fasta"):
    outfile_name = record.id.replace('.','_') + '.fasta'
    record.seq = Seq(str(record.seq).replace('B','').replace('X','').replace('J','').replace('U',''))
    with open(os.path.join(split_fasta_outfolder, outfile_name), 'w') as splitfile:
        SeqIO.write(record, splitfile, "fasta")
