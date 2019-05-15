#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io
import shutil
import snakemake
import subprocess
import sys
from pkg_resources import resource_filename


def main():
    # get the snakefile
    snakefile = resource_filename(__name__, 'Snakefile')

    # set args
    args = {
        'chunks': 1,
        'wait': 60,
        'fraction': 1,
        'reads': 'reads.fq',
        'assembly': 'assembly.fasta',
        'outdir': 'blerf',
        'output_filename': 'racon.fasta'
        }

    # store old stdout
    stdout = sys.stdout
    # call snakemake api and capture output
    sys.stdout = io.StringIO()
    snakemake.snakemake(
        snakefile,
        config=args,
        dryrun=True,
        printrulegraph=True,
        forceall=True)
    output = sys.stdout.getvalue()
    # restore sys.stdout
    sys.stdout = stdout
    # write output
    dag_prefix = 'dag'
    if shutil.which('dot'):
        svg_file = '{}.svg'.format(dag_prefix)
        # pipe the output to dot
        with open(svg_file, 'wb') as svg:
            dot_process = subprocess.Popen(
                ['dot', '-Tsvg'],
                stdin=subprocess.PIPE,
                stdout=svg)
        dot_process.communicate(input=output.encode())
    else:
        # just write the dag to file
        dag_file = '{}.dag'.format(dag_prefix)
        with open(dag_file, 'wt') as file:
            file.write(output)



    # # run the pipeline
    # snakemake.snakemake(
    #     snakefile=snakefile,
    #     config=args,
    #     cores=args['threads'],
    #     lock=False,
    #     printreason=True,
    #     printshellcmds=True,
    #     dryrun=True if args['dry_run'] else False)

if __name__ == '__main__':
    main()
