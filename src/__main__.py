#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging
import snakemake
from pkg_resources import resource_filename


#############
# FUNCTIONS #
#############

def parse_arguments():
    # parse arguments
    parser = argparse.ArgumentParser(
        prog='racon_chunks')
    parser.add_argument(
        '--reads',
        required=True,
        help='Reads in [interleaved,] uncompressed fastq',
        type=str,
        dest='reads')
    parser.add_argument(
        '--assembly',
        required=True,
        help='Genome assembly in uncompressed fasta',
        type=str,
        dest='assembly')
    parser.add_argument(
        '--outdir',
        required=True,
        help='Output directory',
        type=str,
        dest='outdir')
    parser.add_argument(
        '--output_filename',
        required=False,
        help='Filename for the polished fasta file',
        type=str,
        dest='output_filename',
        default='racon.fasta')
    default_threads = 1
    parser.add_argument(
        '--threads',
        help=('Number of threads. Default: %i' % default_threads),
        metavar='int',
        type=int,
        dest='threads',
        default=default_threads)
    default_chunks = 2000
    parser.add_argument(
        '--chunks',
        help=('Number of chunks. Default: %i' % default_chunks),
        metavar='int',
        type=int,
        dest='chunks',
        default=default_chunks)
    default_wait = 60
    parser.add_argument(
        '--wait_min',
        help=(('Number of minutes to wait for each racon job. '
               'Default: %i') % default_wait),
        metavar='int',
        type=int,
        dest='wait',
        default=default_wait)
    default_fraction = 1.0
    parser.add_argument(
        '--fraction',
        help=('Fraction of fastq reads to use. Default: %f' % default_fraction),
        metavar='float',
        type=float,
        dest='fraction',
        default=default_fraction)

    args = vars(parser.parse_args())
    return(args)


########
# MAIN #
########

def main():
    # set up log
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.DEBUG)

    # get the snakefile
    snakefile = resource_filename(__name__, 'Snakefile') # SET THIS UP IN setup.py
    logging.debug(f'Using snakefile {snakefile}')

    # get args
    args = parse_arguments()
    logging.debug('args:')
    logging.debug(args)

    # run the pipeline
    snakemake.snakemake(
        snakefile=snakefile,
        config=args,
        cores=args['threads'],
        lock=False,
        printreason=True,
        printshellcmds=True)


if __name__ == '__main__':
    main()