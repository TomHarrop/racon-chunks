#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools import find_packages


# load README.rst
def readme():
    with open('README.rst') as file:
        return file.read()


setup(
    name='racon_chunks',
    version='0.0.6',
    description='python3 wrapper for racon_chunks',
    long_description=readme(),
    url='https://github.com/tomaharrop/racon_chunks',
    author='Tom Harrop',
    author_email='twharrop@gmail.com',
    license='GPL-3',
    packages=find_packages(),
    install_requires=[
        'biopython>=1.73',
        'psutil>=5.6.2',
        'snakemake>=5.4.5'
    ],
    scripts=[
        'racon_chunks/src/retrieve_reads.py',
    ],
    entry_points={
        'console_scripts': [
            'racon_chunks = racon_chunks.__main__:main',
            'print_dag = racon_chunks.print_dag:main'
            ],
    },
    package_data={
        'racon_chunks': [
            'Snakefile',
            'README.rst'
        ],
    },
    zip_safe=False)
