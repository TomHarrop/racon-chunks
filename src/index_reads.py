#!/usr/bin/env python3

import Bio
from Bio import SeqIO
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


read_file = snakemake.input[0]
db_file = snakemake.output[0]

try:
    read_index = SeqIO.index_db(db_file,
                                read_file,
                                'fastq')
except Exception as e:
    logging.exception('')
    raise e
