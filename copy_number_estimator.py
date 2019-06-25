#!/usr/bin/env python

from scipy.stats import norm
import multiprocessing
from Bio import SeqIO
import numpy as np
import tempfile
import argparse
import sqlite3
import logging
import shutil
import pysam
import glob
import os


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--forward_reads',
                        type=str,
                        required=True,
                        help='Path to forward reads.')
    parser.add_argument('-2', '--reverse_reads',
                        type=str,
                        required=True,
                        help='Path to forward reads.')
    parser.add_argument('-f', '--gene_fasta_file',
                        type=str,
                        required=True,
                        help='Path to FASTA-formatted file containing your gene (or genes?) of interest.')
    parser.add_argument('-p', '--p_value',
                        default=0.05,
                        type=float,
                        help='P-value to report on. NOT YET IMPLEMENTED')
    parser.add_argument('-t', '--threads',
                        default=multiprocessing.cpu_count(),
                        type=int,
                        help='Number of threads to run analysis on. NOT YET IMPLEMENTED.')
    parser.add_argument('--multiple_test_correction',
                        default=False,
                        action='store_true',
                        help='Add this flag to do multiple test correction if you have multiple genes. NOT YET IMPLEMENTED')
    # TODO: Store this db someplace better by default.
    parser.add_argument('--database',
                        default='test_db',
                        type=str,
                        help='Database to store some stuff in.')
    return parser.parse_args()


def find_variant_proportions(bamfile):
    # Not sure we'll ever actually end up using this...
    bamfile = pysam.AlignmentFile(bamfile, 'rb')
    proportions_list = list()
    for column in bamfile.pileup():
        base_list = [0, 0, 0, 0]
        for read in column.pileups:
            if not read.is_del and not read.is_refskip:
                if read.alignment.query_sequence[read.query_position] == 'A':
                    base_list[0] += 1
                if read.alignment.query_sequence[read.query_position] == 'C':
                    base_list[1] += 1
                if read.alignment.query_sequence[read.query_position] == 'T':
                    base_list[2] += 1
                if read.alignment.query_sequence[read.query_position] == 'G':
                    base_list[3] += 1
        proportion = max(base_list)/sum(base_list)
        if 0.1 <= proportion <= 0.9:
            proportions_list.append(proportion)
    bamfile.close()
    return np.array(proportions_list)


def get_copy_of_gene(forward_reads, reverse_reads, multigene_fasta):
    # Use KMA to figure out which allele of a gene is present in a set of reads.

    # Create kma database if it doesn't exist already.
    if not os.path.isfile(multigene_fasta.replace('.fasta', '_kma.name')):
        cmd = 'kma index -i {} -o {}'.format(multigene_fasta, multigene_fasta.replace('.fasta', '_kma'))
        os.system(cmd)
    # Run KMA, writing results to a tmpdir as we go.
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = 'kma -t 12 -ipe {} {} -t_db {} -o {}'.format(forward_reads, reverse_reads, multigene_fasta.replace('.fasta', '_kma'),
                                                           os.path.join(tmpdir, 'kma_results'))
        os.system(cmd)
        # Parse the result file to find out what version of the gene we actually have present.
        with open(os.path.join(tmpdir, 'kma_results.res')) as f:
            lines = f.readlines()
        gene_name = lines[1].split()[0]
        # Probably use SeqIO.index here, should be faster than parsing through the whole file...
        for s in SeqIO.parse(multigene_fasta, 'fasta'):
            if s.id == gene_name:
                return s


def main():
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    args = get_args()
    """
    Here's the plan for CopyNumberEstimator (name is subject to change) - probably wouldn't be too bad to adapt this to
    work with PacBio/Nanopore as well, but given that you get good read length with them it's much less of an issue.

    Big questions:
    How many single copy genes to include? More is probably better, but will eventually hit a point where we don't need
    any more to approximate the distribution of depths and we'll just be burning compute time.

    How to approximate distributions for multi-copy genes? My intuition is that the variance on distributions for genes
    with more copies will be bigger, but this should actually get backed up somehow if possible. Do some empirical
    testing to find out?

    What distribution to fit depth of coverage to? This is definitely very important and should be thought about lots.
    Extremely preliminary, but a normal distribution looks about right. Can then estimate mean and stdev, and then
    get probability fairly easily.

    Implementation:
    0) For all of our conserved single copy genes, figure out which allele is actually present in our reads. Use KMA
    for this for now, but may need to find something better at some point. KMA has been more buggy than I'd like it
    to be.
    1) Take a set of reads (going to assume Illumina, paired end) and align them back to a set of conserved, single
    copy genes.
    2) Also align that set of reads against whatever our target gene (or genes are) to figure out their read depths.
    3) Potentially do some correction of read depths based on GC percent, which literature says can bias things. This
    would likely get done as part of steps 1 and 2 and wouldn't really be a separate step. CNOGpro paper has an equation
    that might be steal-able
    4) Take the depths from our known single copy genes and fit them to some sort of distribution (read some papers to
    figure out what kind of distribution we'd expect - normal distribution might work?)
    5) Compare the depths for each of our genes of interest to the distribution and do some stats to figure out if
    we're significantly outside of the distribution. Could extrapolate distribution to other copy numbers (so if
    single copy genes average 50X depth, assume same distribution around 100X corresponds to two copies? Likely to end up
    with more variance as copy number gets higher, so that would have to be taken into consideration as well).
    6) Optionally, can also do a search for minor variants that correspond to multiple alleles of a gene being present
    in a genome.
    """
    # TODO: Implement GC correction on depth calcuations!
    if check_dependencies() is False:
        logging.error('Some necessary dependencies not found. Exiting..')
        quit(code=1)

    sample_name = os.path.split(args.forward_reads)[1].split('_R1')[0]
    db_check = check_for_sample_in_database(sample_name=sample_name,
                                            db_name=args.database)
    if db_check is False:
        ref_fastas = glob.glob('single_copy_genes/*.fasta')
        gene_names = list()
        sequences = list()
        for fasta in ref_fastas:
            s = get_copy_of_gene(args.forward_reads, args.reverse_reads, fasta)
            gene_names.append(s.id)
            sequences.append(s)
        SeqIO.write(sequences, 'genes_and_stuff.fasta', 'fasta')
        depths = get_depths(args.forward_reads, args.reverse_reads, 'genes_and_stuff.fasta', gene_names)
        depth_list = list()
        for gene_name in depths:
            depth_list.append(depths[gene_name])
            print('{},{}'.format(gene_name, depths[gene_name]))
        avg, std = norm.fit(np.array(depth_list))
        add_sample_to_database(sample_name=sample_name,
                               avgdepth=avg,
                               stdev=std,
                               db_name=args.database)

    else:
        avg, std = db_check
    # Now find out if our gene of interest is outside our distribution.
    gene_names = list()
    for s in SeqIO.parse(args.gene_fasta_file, 'fasta'):
        gene_names.append(s.id)
    depths = get_depths(args.forward_reads, args.reverse_reads, args.gene_fasta_file, gene_names)
    logging.info('Average:{}\nStdev:{}\n'.format(avg, std))
    print('Gene,Depth,CDF')
    for gene_name in depths:
        depth = depths[gene_name]
        print(find_most_likely_copy_number(single_copy_mean=avg,
                                           single_copy_stdev=std,
                                           max_copy_number=10,
                                           depth=depth))


def check_for_sample_in_database(sample_name, db_name):
    # Check our sqlite database to see if we've already examined this same sample before and can skip a bunch
    # of fairly intensive computation.
    sample_name = (sample_name, )  # sqlite expects this to be a tuple, not just a string.
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    # Step 1: Create our table if it doesn't already exist.
    c.execute('CREATE TABLE IF NOT EXISTS depths (strain text, avgdepth real, stdev real)')
    conn.commit()
    # See if our sample is there.
    c.execute('SELECT * FROM depths WHERE strain=?', sample_name)
    table_row = c.fetchone()
    conn.close()
    # If sample isn't present, return false
    if table_row is None:
        return False
    # If sample is present, return average depth and standard deviation
    else:
        name, avg, stdev = table_row
        return (avg, stdev)


def add_sample_to_database(sample_name, avgdepth, stdev, db_name):
    data_to_add = (sample_name, avgdepth, stdev)
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    # Step 1: Create our table if it doesn't already exist.
    c.execute('CREATE TABLE IF NOT EXISTS depths (strain text, avgdepth real, stdev real)')
    conn.commit()
    # Now add our sample and commit changes.
    c.execute('INSERT INTO depths VALUES (?, ?, ?)', data_to_add)
    conn.commit()
    conn.close()


def find_most_likely_copy_number(single_copy_mean, single_copy_stdev, max_copy_number, depth):
    copy_number_dict = dict()
    current_mean = single_copy_mean
    current_stdev = single_copy_stdev
    for i in range(max_copy_number):
        prob = norm.cdf(depth, current_mean, current_stdev)
        copy_number_dict[i + 1] = prob
        current_mean += single_copy_mean
        # Add 10 percent to current standard deviation so distributions get fatter as we get more copy number.
        # This is definitely my intuition on how things should work, but needs to actually be modeled at some point.
        current_stdev = current_stdev * 1.1
    return copy_number_dict


def get_depths(forward_reads, reverse_reads, gene_fasta, gene_names):
    depths = dict()
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpfasta = os.path.join(tmpdir, gene_fasta)
        shutil.copy(gene_fasta, tmpfasta)
        outbam = os.path.join(tmpdir, 'out.bam')
        # Make our bowtie2 index
        bt2_index = os.path.join(tmpdir, 'bowtie_db')
        cmd = 'bowtie2-build {} {}'.format(tmpfasta, bt2_index)
        os.system(cmd)
        # Now align our reads against created database. Use very sensitive local bowtie settings so we can be
        # relatively sure we've gotten everything.
        cmd = 'bowtie2 -p 12 --very-sensitive-local -x {} -1 {} -2 {} | samtools view -bS > {}'.format(bt2_index,
                                                                                                       forward_reads,
                                                                                                       reverse_reads,
                                                                                                       outbam)
        os.system(cmd)
        # Sort and index created bamfile so we can parse it
        sorted_bam = os.path.join(tmpdir, 'out_sorted.bam')
        cmd = 'samtools sort {} > {}'.format(outbam, sorted_bam)
        os.system(cmd)
        cmd = 'samtools index {}'.format(sorted_bam)
        os.system(cmd)
        os.system('samtools faidx {}'.format(gene_fasta))
        # Now parse through the bamfile and get average depth for each of the sequences.
        bamfile = pysam.AlignmentFile(sorted_bam, 'rb')
        for gene_name in gene_names:
            total_depth = 0
            bases = 0
            for column in bamfile.pileup(gene_name, stepper='samtools', ignore_orphans=False,
                                         fastafile=pysam.FastaFile(gene_fasta),
                                         min_base_quality=0):
                bases += 1
                total_depth += column.n
            try:
                depths[gene_name] = total_depth/bases
            except ZeroDivisionError:
                pass
        bamfile.close()
    return depths


def check_dependencies():
    logging.info('Checking for external dependencies.')
    all_dependencies_found = True
    dependencies = ['kma', 'bowtie2', 'bowtie2-build']
    for dependency in dependencies:
        if shutil.which(dependency) is None:
            all_dependencies_found = False
            logging.warning('WARNING: {} not found.'.format(dependency))
        else:
            logging.info('Checking for {}... Found!'.format(dependency))
    return all_dependencies_found


if __name__ == '__main__':
    main()
