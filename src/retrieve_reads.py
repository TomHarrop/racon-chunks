#!/usr/bin/env python3

import logging
import os
from Bio import SeqIO

# set up log
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    filename=snakemake.log[0],
    level=logging.DEBUG)

# catch files from snakemake
read_id_list = snakemake.input['read_ids']
read_file = snakemake.input['fastq']
outdir = snakemake.params['outdir']
read_no = snakemake.params['read_no']

# dev
# read_id_list = ['output/040_read-chunks/chunk_87.txt',
#                 'output/040_read-chunks/chunk_999.txt']
# read_file = 'test/r1.fq'
# outdir = 'test'

# dict of chunk to outfile
chunk_to_outfile = {os.path.basename(x).rstrip('.txt'):
                    os.path.join(outdir,
                                 os.path.basename(x).rstrip('.txt') +
                                 f'_r{read_no}.fq')
                    for x in read_id_list}

# initialise a dict of read_id to chunk
read_to_chunk = dict()

# loop over read_id_list 
for read_id_file in read_id_list:
    my_chunk_id = os.path.basename(read_id_file).rstrip('.txt')
    logging.info(f'Processing {my_chunk_id}')
    # get the list of read ids
    with open(read_id_file, 'rt') as f:
        ids = [x.rstrip('\n') for x in f.readlines()]
    # check if id is already indexed
    for id in ids:
        if id in read_to_chunk:
            read_to_chunk[id] = read_to_chunk[id].append(my_chunk_id)
        else:
            read_to_chunk[id] = [my_chunk_id]

# initialise a dict of chunk to seqrecords
chunk_to_seqrecords = {key: [] for key in set(os.path.basename(x).rstrip('.txt')
                                              for x in read_id_list)}

# read the fastq into memory
logging.info(f'Reading {read_file}')
for seq_rec in SeqIO.parse(read_file, 'fastq'):
    try:
        write_chunks = read_to_chunk[seq_rec.id]
        # I don't think we can append a fastq record to a file, so we have to
        # hold this in memory.
        for chunk in write_chunks:
            chunk_to_seqrecords[chunk].append(seq_rec)
    except KeyError as e:
        logging.debug(f'Read {seq_rec.id} not in dict(read_to_chunk)')

# write the output
logging.info(f'Writing output to {outdir}')
for chunk in chunk_to_seqrecords:
    logging.info(f'Writing {chunk}')
    SeqIO.write(chunk_to_seqrecords[chunk],
                chunk_to_outfile[chunk],
                'fastq')
