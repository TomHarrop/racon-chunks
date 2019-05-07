#!/usr/bin/env python3

import Bio
from Bio import SeqIO
import csv
import logging
import os
import platform
import sqlite3
import sys

# set up log
logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.DEBUG)

# debug biopython issue
logging.debug('sys.version')
logging.debug(sys.version)
logging.debug('sqlite3.version')
logging.debug(sqlite3.version)
logging.debug('platform.python_implementation()')
logging.debug(platform.python_implementation())
logging.debug('platform.platform()')
logging.debug(platform.platform())
logging.debug('Bio.__version__')
logging.debug(Bio.__version__)
logging.debug('os.environ')
logging.debug(os.environ)


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
            logging.warning(f'{read} not in database')
            pass

# write output
SeqIO.write(r1_reads,
            r1_out,
            'fastq')
SeqIO.write(r2_reads,
            r2_out,
            'fastq')
