#!/usr/bin/env python3

import multiprocessing

###########
# GLOBALS #
###########

ways = 100

samtools = 'shub://TomHarrop/singularity-containers:samtools_1.9'
bbmap = 'shub://TomHarrop/singularity-containers:bbmap_38.00'

########
# MAIN #
########

all_chunks = [str(x) for x in range(1, ways + 1)]

#########
# RULES #
#########


# chunk the fasta file
rule partition:
    input:
        'data/flye_denovo_full.racon.fasta'
    output:
        expand('output/010_data/chunks/chunk_{chunk}.fasta',
               chunk=all_chunks)
    params:
        outfile = 'output/010_data/chunks/chunk_%.fasta',
        ways = ways
    log:
        'logs/010_data/partition.log'
    threads:
        1
    priority:
        10
    singularity:
        bbmap
    shell:
        'partition.sh '
        'in={input} '
        'out={params.outfile} '
        'ways={params.ways} '
        '2> {log}'

# sort and index the bwa aln.sam file
rule sort_sam:
    input:
        'data/aln.sam'
    output:
        bam = 'output/010_data/aln_sorted.bam',
        bai = 'output/010_data/aln_sorted.bam.bai'
    log:
        'logs/010_data/sort_sam.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        samtools
    shell:
        'samtools sort '
        '-l 0 '
        '-m 7G '
        '-O BAM '
        '-@ {threads} '
        '{input} '
        '> {output.bam} '
        '2> {log} '
        '&& '
        'samtools index {output.bam} {output.bai} '
        '2>> {log}'
