#!/usr/bin/env python3

from Bio import SeqIO

print(sys.version)

read_file = snakemake.input[0]
db_file = snakemake.output[0]

read_index = SeqIO.index_db(db_file,
                            read_file,
                            'fastq')
