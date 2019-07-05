#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import argparse
import csv


def create_plot(mean, stdev, max_copy, stdev_multiplier, gene_depth, gene_name):
    xticks = list()
    for i in range(1, max_copy + 1):
        xticks.append(mean * i)
        samples = np.random.normal(loc=mean * i, scale=stdev, size=1000)
        plt.hist(samples, alpha=0.5, label='{}X'.format(i), bins=20)
        stdev *= stdev_multiplier
    plt.plot([gene_depth, gene_depth], [0, 125], 'black', label=gene_name)
    # plt.legend()
    plt.xticks(xticks)
    plt.title('Copy Number for {}'.format(gene_name))
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--copy_number_report',
                        type=str,
                        required=True,
                        help='Report created by running copy_number_estimator.py')
    parser.add_argument('-g', '--gene',
                        type=str,
                        required=True,
                        help='Name of gene in report you want to create a visualization for.')
    parser.add_argument('-m', '--max_copy_number',
                        type=int,
                        default=5,
                        help='Copy number to go up to.')
    args = parser.parse_args()
    gene_found = False
    with open(args.copy_number_report) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['Gene'] == args.gene:
                create_plot(mean=float(row['SingleCopyDepth']),
                            stdev=float(row['SingleCopyStdev']),
                            max_copy=args.max_copy_number,
                            stdev_multiplier=float(row['StdevMultiplier']),
                            gene_depth=float(row['Depth']),
                            gene_name=row['Gene'])
                gene_found = True
    if gene_found is False:
        print('Could not find specified gene in report.')

