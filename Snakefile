#!/usr/bin/env python3

import multiprocessing

###########
# GLOBALS #
###########

n_chunks = 1000

samtools = 'shub://TomHarrop/singularity-containers:samtools_1.9'
bbmap = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
bwa = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'

########
# MAIN #
########

all_chunks = [str(x) for x in range(0, n_chunks)]

reads = 'data/pe_reads.fq'
assembly = 'data/flye_denovo_full.racon.fasta'
alignment = 'data/aln.sam'

#########
# RULES #
#########

rule target:
    input:
        expand('output/010_chunks/chunk_{chunk}.fasta',
               chunk=all_chunks),
        'output/020_alignment/aln_sorted.bam'

# chunk the fasta file
rule partition:
    input:
        assembly
    output:
        expand('output/010_chunks/chunk_{chunk}.fasta',
               chunk=all_chunks)
    params:
        outfile = 'output/010_chunks/chunk_%.fasta',
        ways = n_chunks
    log:
        'logs/010_chunks/partition.log'
    threads:
        1
    singularity:
        bbmap
    shell:
        'partition.sh '
        'in={input} '
        'out={params.outfile} '
        'ways={params.ways} '
        '2> {log}'

# sort and index the bwa aln.sam file
rule index_bam:
    input:
        bam = 'output/020_alignment/aln_sorted.bam',
    output:
        bai = 'output/020_alignment/aln_sorted.bam.bai'
    log:
        'logs/020_alignment/index_bam.log'
    threads:
        1
    singularity:
        samtools
    shell:
        'samtools index {input.bam} {output.bai} '
        '2> {log}'

rule sort_sam:
    input:
        'output/020_alignment/aln.sam'
    output:
        bam = 'output/020_alignment/aln_sorted.bam',
    log:
        'logs/020_alignment/sort_sam.log'
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

# map the reads to the whole fasta
rule map_reads:
    input:
        index = expand('output/020_alignment/index.{suffix}',
                       suffix=['amb', 'ann', 'bwt', 'pac', 'sa']),
        fq = reads
    output:
        temp('output/020_alignment/aln.sam')
    params:
        prefix = 'output/020_alignment/index'
    log:
        'logs/020_alignment/bwa-mem.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        bwa
    shell:
        'bwa mem '
        '-t {threads} '
        '-p -C '
        '{params.prefix} '
        '{input.fq} '
        '> {output} '
        '2> {log}'

rule index_assembly:
    input:
        fasta = assembly
    output:
        expand('output/020_alignment/index.{suffix}',
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix = 'output/020_alignment/index'
    log:
        'logs/020_alignment/bwa-index.log'
    threads:
        1
    singularity:
        bwa
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input.fasta} '
        '2> {log} '

