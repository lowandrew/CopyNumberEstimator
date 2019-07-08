#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="copy_number_estimator",
    version="0.0.1",
    packages=find_packages(),
    entry_points={
       'console_scripts': [
            'copy_number_estimator = copy_number_estimator.copy_number_estimator:main',
            'copy_number_database_setup = copy_number_estimator.download_single_copy_genes:main',
            'copy_number_visualizer = copy_number_estimator.copy_visualizer:main'
       ],
    },
    author="Andrew Low",
    author_email="andrew.low@canada.ca",
    url="https://github.com/OLC-Bioinformatics/ConFindr",
    install_requires=['biopython',
                      'matplotlib',
                      'pysam',
                      'numpy',
                      'scipy',
                      'requests']
)
