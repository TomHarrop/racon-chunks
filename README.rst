racon_chunks
==================

A python3 wrapper for the running racon on chunks, using snakemake_.

.. image:: https://raw.githubusercontent.com/TomHarrop/racon_chunks/master/dag.svg

Requirements
------------

If Singularity_ is available, download and run the racon_chunks_ container.

Otherwise, the following must be installed:

* ``racon``
* ``samtools``
* ``bwa``
* ``bbmap``

Installation
------------

``pip3 install git+git://github.com/tomharrop/racon_chunks.git``

Usage
-----

.. code::

  racon_chunks [-h] --reads READS --assembly ASSEMBLY --outdir OUTDIR
               [--output_filename OUTPUT_FILENAME] [--threads int]
               [--chunks int] [--wait_min int] [--fraction float]

  optional arguments:
    -h, --help            show this help message and exit
    --reads READS         Reads in [interleaved,] uncompressed fastq
    --assembly ASSEMBLY   Genome assembly in uncompressed fasta
    --outdir OUTDIR       Output directory
    --output_filename OUTPUT_FILENAME
                          Filename for the polished fasta file
    --threads int         Number of threads. Default: 1
    --chunks int          Number of chunks. Default: 2000
    --wait_min int        Number of minutes to wait for each racon job. Default:
                          60
    --fraction float      Fraction of fastq reads to use. Default: 1.000000


.. _Singularity: https://www.sylabs.io/singularity/
.. _snakemake: https://snakemake.readthedocs.io/en/stable/
.. _racon_chunks: https://www.singularity-hub.org/containers/8716