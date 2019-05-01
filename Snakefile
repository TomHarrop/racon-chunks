#!/usr/bin/env python3

import multiprocessing


#############
# FUNCTIONS #
#############

def read_contig_list(contig_list):
    with open(contig_list, 'rt') as f:
        contigs = [x.rstrip() for x in f.readlines()]
    return(' '.join(contigs))


###########
# GLOBALS #
###########

n_chunks = 1000
wait_mins = 30

racon_chunks = 'shub://TomHarrop/singularity-containers:racon-chunks'

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

singularity: racon_chunks

rule target:
    input:
        expand('output/050_racon/chunk_{chunk}.fasta',
               # chunk=all_chunks)
               chunk=['87']) # just test the pipeline

# run racon on the chunks
rule racon:
    input:
        fasta = 'output/010_chunks/chunk_{chunk}.fasta',
        aln = 'output/030_bam-chunks/chunk_{chunk}.sam',
        fq = 'output/040_read-chunks/chunk_{chunk}.fq'
    output:
        'output/050_racon/chunk_{chunk}.fasta'
    params:
        wait_mins = f'{wait_mins}m'
    log:
        'logs/050_racon/chunk_{chunk}.log'
    threads:
        multiprocessing.cpu_count()
    shell:
        'timeout {params.wait_mins} '
        'racon '
        '-t {threads} '
        '{input.fq} '
        '{input.aln} '
        '{input.fasta} '
        '> {output} '
        '2> {log}'

# retrieve the reads from the bam chunk
rule repair_reads:
    input:
        r1 = 'output/040_read-chunks/chunk_{chunk}_r1.fq',
        r2 = 'output/040_read-chunks/chunk_{chunk}_r2.fq'
    output:
        'output/040_read-chunks/chunk_{chunk}.fq'
    log:
        'logs/040_read-chunks/repair_reads_{chunk}.fq'
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'repair=t '
        'out={output} '
        '&> {log}'

rule retrieve_reads:
    input:
        sam = 'output/030_bam-chunks/chunk_{chunk}.sam',
        r1_idx = 'output/000_reads/r1.idx',
        r2_idx = 'output/000_reads/r2.idx'
    output:
        r1 = temp('output/040_read-chunks/chunk_{chunk}_r1.fq'),
        r2 = temp('output/040_read-chunks/chunk_{chunk}_r2.fq')
    log:
        'logs/040_read-chunks/retrieve_reads_{chunk}.log'
    benchmark:
        'benchmarks/040_read-chunks/retrieve_reads_{chunk}.log'
    script:
        'src/retrieve_reads.py'


# subset the BAM by the chunk list
rule chunk_bam:
    input:
        bam = 'output/020_alignment/aln_sorted.bam',
        bai = 'output/020_alignment/aln_sorted.bam.bai',
        contig_list = 'output/010_chunks/chunk_{chunk}_contigs.txt'
    output:
        'output/030_bam-chunks/chunk_{chunk}.sam'
    log:
        'logs/030_bam-chunks/view_{chunk}.log'
    threads:
        1
    shell:
        # the horrible sed command replaces the newlines with spaces in
        # contig_list. This allows the workflow to run before contig_list is
        # created. Adapted from https://stackoverflow.com/a/1252191
        'contigs="$(sed -e \':a\' -e \'N\' -e \'$!ba\' '
        '-e \'s/\\n/ /g\' {input.contig_list})" ; '
        'samtools view '
        '-h '
        '-F 256 '       # exclude secondary alignments
        '-O SAM '
        '{input.bam} '
        '${{contigs}} '
        '> {output} '
        '2> {log}'

# chunk the fasta file
rule list_contigs:
    input:
        'output/010_chunks/chunk_{chunk}.fasta'
    output:
        'output/010_chunks/chunk_{chunk}_contigs.txt'
    threads:
        1
    shell:
        'grep "^>" {input} | cut -d " " -f1 | sed -e \'s/>//g\' > {output}'

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
    shell:
        'samtools sort '
        '-l 0 '
        '-m 7G '
        '-O BAM '
        '--threads {threads} '
        '{input} '
        '> {output.bam} '
        '2> {log} '

# map the reads to the whole fasta
rule map_reads:
    input:
        index = expand('output/020_alignment/index.{suffix}',
                       suffix=['amb', 'ann', 'bwt', 'pac', 'sa']),
        r1 = 'output/000_reads/r1.fq',
        r2 = 'output/000_reads/r2.fq'
    output:
        temp('output/020_alignment/aln.sam')
    params:
        prefix = 'output/020_alignment/index'
    log:
        'logs/020_alignment/bwa-mem.log'
    threads:
        multiprocessing.cpu_count()
    shell:
        'bwa mem '
        '-t {threads} '
        '{params.prefix} '
        '{input.r1} '
        '{input.r2} '
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
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input.fasta} '
        '2> {log} '

# index the reads
rule index_reads:
    input:
        'output/000_reads/r{r}.fq'
    output:
        'output/000_reads/r{r}.idx'
    script:
        'src/index_reads.py'

rule split_reads:
    input:
        reads
    output:
        r1 = 'output/000_reads/r1.fq',
        r2 = 'output/000_reads/r2.fq'
    log:
        'logs/000_reads/split_reads.log'
    threads:
        1
    shell:
        'reformat.sh '
        'in={input} '
        'out=stdout.fq '
        'verifyinterleaved=t '
        'int=t '
        'trimreaddescription=t '
        '2> {log} '
        '| rename.sh '
        'in=stdin.fq '
        'out={output.r1} '
        'out2={output.r2} '
        'int=t '
        '2>> {log}'
