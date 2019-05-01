#!/usr/bin/env python3

from Bio import SeqIO
import csv
import warnings

r1_db_file = snakemake.input['r1_idx']
r2_db_file = snakemake.input['r2_idx']
sam_file = snakemake.input['sam']
r1_out = snakemake.output['r1']
r2_out = snakemake.output['r2']

# get the ids
with open(sam_file, 'rt') as f:
    reader = csv.reader(f, delimiter='\t')
    sam_query_ids = sorted(set(
        x[0] for x in reader if not x[0].startswith('@')))

# open the dbs
r1_db = SeqIO.index_db(r1_db_file)
r2_db = SeqIO.index_db(r2_db_file)

# read the reads, handle exceptions
r1_reads = []
r2_reads = []
for read in sam_query_ids:
    try:
        r1_reads.append(r1_db[read])
        r2_reads.append(r2_db[read])
    except KeyError:
            warnings.warn(f'{read} not in database')
            pass

# write output
SeqIO.write(r1_reads,
            r1_out,
            'fastq')
SeqIO.write(r2_reads,
            r2_out,
            'fastq')
