#!/usr/bin/env python

from scipy.stats import norm
import multiprocessing
from Bio import SeqIO
import numpy as np
import subprocess
import tempfile
import argparse
import sqlite3
import logging
import shutil
import pysam
import json
import glob
import os


def run_cmd(cmd):
    """
    Runs a command using subprocess, and returns both the stdout and stderr from that command
    If exit code from command is non-zero, raises subproess.CalledProcessError
    :param cmd: command to run as a string, as it would be called on the command line
    :return: out, err: Strings that are the stdout and stderr from the command called.
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    out = out.decode('utf-8')
    err = err.decode('utf-8')
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, cmd=cmd)
    return out, err


# Methods for putting JSON into/out of our sqlite database.
# Adapted/taken from: https://chrisostrouchov.com/post/python_sqlite/

def adapt_json(data):
    return (json.dumps(data, sort_keys=True)).encode()


def convert_json(blob):
    return json.loads(blob.decode())


sqlite3.register_adapter(dict, adapt_json)
sqlite3.register_adapter(list, adapt_json)
sqlite3.register_adapter(tuple, adapt_json)
sqlite3.register_converter('JSON', convert_json)


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
    parser.add_argument('-a', '--assembly',
                        type=str,
                        help='If you have an assembly for your genome of interest, specify it here. If you do not have '
                             'an assembly, a quick and dirty assembly will be created with SKESA.')
    parser.add_argument('-p', '--p_value',
                        default=0.05,
                        type=float,
                        help='P-value to report on.')
    parser.add_argument('-t', '--threads',
                        default=multiprocessing.cpu_count(),
                        type=int,
                        help='Number of threads to run analysis on. NOT YET IMPLEMENTED.')
    # TODO: Store this db someplace better by default.
    parser.add_argument('--database',
                        default='test_db',
                        type=str,
                        help='Database to store some stuff in.')
    parser.add_argument('-s', '--single_copy_gene_folder',
                        type=str,
                        required=True,
                        help='Path to folder that contains single copy genes. Must have a .fasta extension.')
    parser.add_argument('-o', '--output_report',
                        type=str,
                        default='copy_number_report.csv',
                        help='CSV file to store your output in.')
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


def quick_and_dirty_assembly(forward_reads, reverse_reads, output_assembly):
    # Some explanation on this: All we want to do with our output assembly is run a quick BLAST to get the locations
    # of some (or hopefully all) of our universal single copy genes. Doesn't matter if our assembly is generally not
    # particularly good as long as those genes get assembled. Parameters for quick assembly from: https://github.com/ncbi/SKESA/issues/11
    # TODO: Investigate kmer size.
    cmd = 'skesa --fastq {forward_reads} --fastq {reverse_reads} --steps 1 --kmer 99 --vector_percent 1 ' \
          '--contigs_out {output_assembly}'.format(forward_reads=forward_reads,
                                                   reverse_reads=reverse_reads,
                                                   output_assembly=output_assembly)
    out, err = run_cmd(cmd)


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
    # Still have a ludicrous number of TODOs before this is anywhere near production ready. Rough list of things
    # is below:

    # Output formatting: Need to decide what output reports will actually look like and then do that instead of having
    # information put up as print statements here and there

    # Add more command line options/make the ones we have do something: threads and pvalue don't do anything right now,
    # the database location option should be set to a sensible default, and more options will almost certainly come up.
    # Should also probably get this set up as one main command with subparsers - at the very least one for running copy
    # number estimation and another for any database manipulations that one might have. Should also have one for
    # downloading/updating databases of single copy genes once I figure out what genes we're actually going to use.

    # Lots of code right now is hilariously inefficient - for example, bam files are getting created multiple times
    # instead of saving one copy and re-parsing it later. Also, option should get added to save the bam/other potential
    # intermediate files that get created to somewhere instead of having them vanish from the face of the earth forever
    # as soon as the leave the scope of the tempdir.

    # Alleles selected for each gene are currently written to hardcoded genes_and_stuff.fasta. Set this to a reasonable
    # default instead.

    # No way to call this from within python right now, only available as CLI. Make some sort of method so other scripts
    # can import and call from within python

    # Lots of places need lots more logging statements.

    # At some point, turn this into an actual package instead of just having it as a script.
    make_copy_number_report(forward_reads=args.forward_reads,
                            reverse_reads=args.reverse_reads,
                            genes_of_interest=args.gene_fasta_file,
                            single_copy_gene_folder=args.single_copy_gene_folder,
                            assembly=args.assembly,
                            pvalue=args.p_value,
                            database=args.database)


def make_copy_number_report(forward_reads, reverse_reads, genes_of_interest, single_copy_gene_folder,
                            assembly=None, pvalue=0.05, threads=1, database='test_db', stdev_multiplier=1.1,
                            max_copy_number=10, read_identifier='_R1', report_file='copy_number_report.csv'):
    if check_dependencies() is False:
        logging.error('Some necessary dependencies not found. Exiting..')
        quit(code=1)

    # Assign a name to the sample
    sample_name = os.path.split(forward_reads)[1].split(read_identifier)[0]
    db_check = check_for_sample_in_database(sample_name=sample_name, db_name=database)

    with tempfile.TemporaryDirectory() as tmpdir:
        if db_check is False:
            if assembly is None:
                logging.info('No assembly provided. Creating rough assembly...')
                assembly_to_use = os.path.join(tmpdir, 'assembly.fasta')
                quick_and_dirty_assembly(forward_reads, reverse_reads, assembly_to_use)
            else:
                shutil.copy(assembly, tmpdir)
                assembly_to_use = os.path.join(tmpdir, os.path.split(assembly)[1])
            seq_index = SeqIO.index(assembly_to_use, 'fasta')
            conserved_count = 0
            # Now time to find single copy genes in our assembly!
            logging.info('Finding single copy genes...')
            single_copy_gene_file = os.path.join(tmpdir, 'single_copy_genes.fasta')
            gene_names = list()
            for fasta_file in sorted(glob.glob(os.path.join(single_copy_gene_folder, '*.fasta'))):
                db_name = fasta_file.replace('.fasta', '')
                # Create diamond database if it doesn't exist already.
                if not os.path.isfile(db_name):
                    cmd = 'diamond makedb --db {db_name} --in {ref_fasta}'.format(db_name=db_name,
                                                                                  ref_fasta=fasta_file)
                    run_cmd(cmd)
                diamond_report = os.path.join(tmpdir, 'diamond.tsv')
                cmd = 'diamond blastx --query {assembly} --db {db_name} -o {outfile}'.format(assembly=assembly_to_use,
                                                                                             db_name=db_name,
                                                                                             outfile=diamond_report)
                run_cmd(cmd)
                with open(diamond_report) as f:
                    lines = f.readlines()
                with open(single_copy_gene_file, 'a+') as f:
                    if len(lines) > 0:
                        x = lines[0].split()
                        contig = x[0]
                        start_pos = int(x[6])
                        end_pos = int(x[7])
                        percent_id = float(x[2])
                        if percent_id > 95:
                            conserved_count += 1
                            f.write('>gene_{}\n'.format(conserved_count))
                            gene_names.append('gene_{}'.format(conserved_count))
                            if start_pos < end_pos:
                                f.write('{}\n'.format(seq_index[contig].seq[start_pos:end_pos]))
                            else:
                                f.write('{}\n'.format(seq_index[contig].seq[end_pos:start_pos]))
            bam = os.path.join(tmpdir, 'sorted_bam.bam')
            logging.info('Calculating depths for single copy genes...')
            depths, gc_dict = get_depths(forward_reads, reverse_reads, single_copy_gene_file, gene_names, out_bam=bam)
            logging.info('Correcting for GC content...')
            corrected_depths = get_corrected_depths(forward_reads, reverse_reads, single_copy_gene_file,
                                                    gene_names, gc_dict, depths, bam=bam)
            uncorr_avg, uncorr_std = norm.fit(np.array(depths))
            logging.info('UNCORRECTED STATS:')
            logging.info('Median depth: {}'.format(np.median(depths)))
            logging.info('Normal distribution mean: {}'.format(uncorr_avg))
            logging.info('Normal distribution stdev: {}'.format(uncorr_std))

            logging.info('\nGC CORRECTED STATS:')
            avg, std = norm.fit(np.array(corrected_depths))
            logging.info('Median depth: {}'.format(np.median(corrected_depths)))
            logging.info('Normal distribution mean: {}'.format(avg))
            logging.info('Normal distribution stdev: {}'.format(std))
            add_sample_to_database(sample_name=sample_name,
                                   avgdepth=avg,
                                   stdev=std,
                                   db_name=database,
                                   gc_dict=gc_dict)
        else:
            avg, std, gc_dict = db_check
            # SQLite returns a bytes string for gc_dict. Fix that.
            gc_dict = json.loads(gc_dict.decode('utf-8'))
            # Also, keys in the dictionary became strings, so fix that too.
            for percent in gc_dict:
                gc_dict[int(percent)] = gc_dict.pop(percent)

        gene_names = list()
        for s in SeqIO.parse(genes_of_interest, 'fasta'):
            gene_names.append(s.id)
        logging.info('Calculating depths for target genes...')
        depths = get_corrected_depths_per_gene(forward_reads, reverse_reads, genes_of_interest, gene_names, gc_dict, avg)
        with open(report_file, 'w') as f:
            f.write('Gene,Depth,PossibleCopyNumbers,SingleCopyDepth,SingleCopyStdev,StdevMultiplier\n')
            for gene_name in depths:
                depth = depths[gene_name]
                possible_numbers = find_possible_copy_numbers(single_copy_mean=avg,
                                                              single_copy_stdev=std,
                                                              stdev_multiplier=stdev_multiplier,
                                                              max_copy_number=max_copy_number,
                                                              gene_depth=depth,
                                                              pval_cutoff=pvalue)
                f.write('{},{},{},{},{},{}\n'.format(gene_name,
                                                     depth,
                                                     ';'.join(possible_numbers),
                                                     avg,
                                                     std,
                                                     stdev_multiplier))
        logging.info('Done! Thanks for using CopyNumberEstimator.')

def find_possible_copy_numbers(single_copy_mean, single_copy_stdev, stdev_multiplier, max_copy_number, gene_depth, pval_cutoff):
    possible_copy_numbers = list()
    current_stdev = single_copy_stdev
    for i in range(1, max_copy_number + 1):
        prob = norm.cdf(gene_depth, single_copy_mean * i, current_stdev)
        if prob <= 0.5:
            p = 2 * prob
        else:
            p = (1 - prob) * 2
        if p >= pval_cutoff:
            possible_copy_numbers.append(str(i))
        current_stdev = current_stdev * stdev_multiplier
    return possible_copy_numbers


def check_for_sample_in_database(sample_name, db_name):
    # Check our sqlite database to see if we've already examined this same sample before and can skip a bunch
    # of fairly intensive computation.
    sample_name = (sample_name, )  # sqlite expects this to be a tuple, not just a string.
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    # Step 1: Create our table if it doesn't already exist.
    c.execute('CREATE TABLE IF NOT EXISTS depths (strain text, avgdepth real, stdev real, gc json)')
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
        name, avg, stdev, gc = table_row
        return (avg, stdev, gc)


def add_sample_to_database(sample_name, avgdepth, stdev, db_name, gc_dict):
    data_to_add = (sample_name, avgdepth, stdev, gc_dict)
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    # Step 1: Create our table if it doesn't already exist.
    c.execute('CREATE TABLE IF NOT EXISTS depths (strain text, avgdepth real, stdev real, gc json)')
    conn.commit()
    # Now add our sample and commit changes.
    c.execute('INSERT INTO depths VALUES (?, ?, ?, ?)', data_to_add)
    conn.commit()
    conn.close()


def find_most_likely_copy_number(single_copy_mean, single_copy_stdev, max_copy_number, depth, stdev_multiplier=1.1):
    copy_number_dict = dict()
    current_mean = single_copy_mean
    current_stdev = single_copy_stdev
    for i in range(max_copy_number):
        prob = norm.cdf(depth, current_mean, current_stdev)
        if prob <= 0.5:
            pvalue = 2 * prob
        else:
            pvalue = (1 - prob) * 2
        copy_number_dict[i + 1] = pvalue
        current_mean += single_copy_mean
        # Add 10 percent to current standard deviation so distributions get fatter as we get more copy number.
        # This is definitely my intuition on how things should work, but needs to actually be modeled at some point.
        current_stdev = current_stdev * stdev_multiplier
    return copy_number_dict


def get_corrected_depths_per_gene(forward_reads, reverse_reads, gene_fasta, gene_names, gc_dict, corrected_median_depth):
    corrected_depths = dict()
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpfasta = os.path.join(tmpdir, gene_fasta)
        try:
            shutil.copy(gene_fasta, tmpfasta)
        except shutil.SameFileError:
            pass
        outbam = os.path.join(tmpdir, 'out.bam')
        # Make our bowtie2 index
        bt2_index = os.path.join(tmpdir, 'bowtie_db')
        cmd = 'bowtie2-build {} {}'.format(tmpfasta, bt2_index)
        run_cmd(cmd)
        # Now align our reads against created database. Use very sensitive local bowtie settings so we can be
        # relatively sure we've gotten everything.
        cmd = 'bowtie2 -p 12 --very-sensitive-local -x {} -1 {} -2 {} | samtools view -bS > {}'.format(bt2_index,
                                                                                                       forward_reads,
                                                                                                       reverse_reads,
                                                                                                       outbam)
        run_cmd(cmd)
        # Sort and index created bamfile so we can parse it
        sorted_bam = os.path.join(tmpdir, 'out_sorted.bam')
        cmd = 'samtools sort {} > {}'.format(outbam, sorted_bam)
        run_cmd(cmd)
        cmd = 'samtools index {}'.format(sorted_bam)
        run_cmd(cmd)
        run_cmd('samtools faidx {}'.format(gene_fasta))
        seq_index = SeqIO.index(tmpfasta, 'fasta')
        bamfile = pysam.AlignmentFile(sorted_bam, 'rb')
        for gene_name in gene_names:
            gene_window_depths = list()
            # Find out what depth is and put it into corrected_depths dictionary.
            gene_length = len(seq_index[gene_name].seq)
            windows = generate_windows(gene_length, window_size=100)
            for window in windows:
                total_depth = 0
                bases = 0
                gc_bases = 0
                start_base, end_base = window
                for column in bamfile.pileup(contig=gene_name, start=start_base, end=end_base,
                                             stepper='samtools', ignore_orphans=False,
                                             fastafile=pysam.FastaFile(gene_fasta),
                                             min_base_quality=0):
                    ref_base = seq_index[gene_name].seq[column.reference_pos]
                    if ref_base.upper() == 'G' or ref_base.upper() == 'C':
                        gc_bases += 1
                    bases += 1
                    total_depth += column.nsegments
                try:
                    depth = total_depth/bases
                    gc_content = round(100*gc_bases/bases)
                    # If we haven't ever seen the gc content of this segment before, can't correct for the gc content.
                    # Unlikely this will ever happen, but better this than crashing with a KeyError
                    if gc_content not in gc_dict:
                        gene_window_depths.append(depth)
                    # If we have seen this GC content before, apply the correction described in:
                    # https://genome.cshlp.org/content/19/9/1586.full and
                    # https://academic.oup.com/bioinformatics/article/31/11/1708/2365681
                    else:
                        corrected_depth = (depth * corrected_median_depth)/np.median(gc_dict[gc_content])
                        gene_window_depths.append(corrected_depth)
                except ZeroDivisionError:
                    pass
            corrected_depths[gene_name] = np.median(gene_window_depths)
        bamfile.close()
    return corrected_depths


def get_corrected_depths(forward_reads, reverse_reads, gene_fasta, gene_names, gc_dict, uncorrected_depths, bam=None):
    corrected_depths = list()
    uncorrected_median = np.median(uncorrected_depths)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpfasta = os.path.join(tmpdir, gene_fasta)
        try:
            shutil.copy(gene_fasta, tmpfasta)
        except shutil.SameFileError:
            pass
        seq_index = SeqIO.index(tmpfasta, 'fasta')
        if bam is None:
            outbam = os.path.join(tmpdir, 'out.bam')
            # Make our bowtie2 index
            bt2_index = os.path.join(tmpdir, 'bowtie_db')
            cmd = 'bowtie2-build {} {}'.format(tmpfasta, bt2_index)
            run_cmd(cmd)
            # Now align our reads against created database. Use very sensitive local bowtie settings so we can be
            # relatively sure we've gotten everything.
            cmd = 'bowtie2 -p 12 --very-sensitive-local -x {} -1 {} -2 {} | samtools view -bS > {}'.format(bt2_index,
                                                                                                           forward_reads,
                                                                                                           reverse_reads,
                                                                                                           outbam)
            run_cmd(cmd)
            # Sort and index created bamfile so we can parse it
            sorted_bam = os.path.join(tmpdir, 'out_sorted.bam')
            cmd = 'samtools sort {} > {}'.format(outbam, sorted_bam)
            run_cmd(cmd)
            cmd = 'samtools index {}'.format(sorted_bam)
            run_cmd(cmd)
            run_cmd('samtools faidx {}'.format(gene_fasta))
        else:
            sorted_bam = bam
            run_cmd('samtools faidx {}'.format(gene_fasta))
            cmd = 'samtools index {}'.format(sorted_bam)
            run_cmd(cmd)
        bamfile = pysam.AlignmentFile(sorted_bam, 'rb')
        for gene_name in gene_names:
            # Need to break our gene up into 100 bp windows. Unless the gene we're looking at happens to be a multiple
            # of 100 in terms of length, this won't quite work. To avoid super small sample sizes that could mess up
            # our GC correction, set our last window to be the last 100 bases of the gene.
            # Example: for a gene of length 450, windows go from base 1-100, 101-200, 201-300, 301-400, and 351-450
            # Does mean we end up oversampling some bases, but seems to me to be better than the alternative of having
            # really short windows with potentially extremely biased GC content/depths
            gene_length = len(seq_index[gene_name].seq)
            windows = generate_windows(gene_length, window_size=100)
            for window in windows:
                total_depth = 0
                bases = 0
                gc_bases = 0
                start_base, end_base = window
                for column in bamfile.pileup(contig=gene_name, start=start_base, end=end_base,
                                             stepper='samtools', ignore_orphans=False,
                                             fastafile=pysam.FastaFile(gene_fasta),
                                             min_base_quality=0):
                    ref_base = seq_index[gene_name].seq[column.reference_pos]
                    if ref_base.upper() == 'G' or ref_base.upper() == 'C':
                        gc_bases += 1
                    bases += 1
                    total_depth += column.nsegments
                try:
                    depth = total_depth/bases
                    gc_content = round(100*gc_bases/bases)
                    # If we haven't ever seen the gc content of this segment before, can't correct for the gc content.
                    # Unlikely this will ever happen, but better this than crashing with a KeyError
                    if gc_content not in gc_dict:
                        corrected_depths.append(depth)
                    # If we have seen this GC content before, apply the correction described in:
                    # https://genome.cshlp.org/content/19/9/1586.full and
                    # https://academic.oup.com/bioinformatics/article/31/11/1708/2365681
                    else:
                        corrected_depth = (depth * uncorrected_median)/np.median(gc_dict[gc_content])
                        corrected_depths.append(corrected_depth)
                except ZeroDivisionError:
                    pass
        bamfile.close()
    return corrected_depths


def get_depths(forward_reads, reverse_reads, gene_fasta, gene_names, out_bam=None):
    depths = list()
    gc_dict = dict()
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpfasta = os.path.join(tmpdir, gene_fasta)
        try:
            shutil.copy(gene_fasta, tmpfasta)
        except shutil.SameFileError:
            pass
        outbam = os.path.join(tmpdir, 'out.bam')
        # Make our bowtie2 index
        bt2_index = os.path.join(tmpdir, 'bowtie_db')
        cmd = 'bowtie2-build {} {}'.format(tmpfasta, bt2_index)
        run_cmd(cmd)
        # Now align our reads against created database. Use very sensitive local bowtie settings so we can be
        # relatively sure we've gotten everything.
        cmd = 'bowtie2 -p 12 --very-sensitive-local -x {} -1 {} -2 {} | samtools view -bS > {}'.format(bt2_index,
                                                                                                       forward_reads,
                                                                                                       reverse_reads,
                                                                                                       outbam)
        run_cmd(cmd)
        # Sort and index created bamfile so we can parse it
        sorted_bam = os.path.join(tmpdir, 'out_sorted.bam')
        cmd = 'samtools sort {} > {}'.format(outbam, sorted_bam)
        run_cmd(cmd)
        cmd = 'samtools index {}'.format(sorted_bam)
        run_cmd(cmd)
        run_cmd('samtools faidx {}'.format(gene_fasta))
        seq_index = SeqIO.index(tmpfasta, 'fasta')
        bamfile = pysam.AlignmentFile(sorted_bam, 'rb')
        for gene_name in gene_names:
            # Need to break our gene up into 100 bp windows. Unless the gene we're looking at happens to be a multiple
            # of 100 in terms of length, this won't quite work. To avoid super small sample sizes that could mess up
            # our GC correction, set our last window to be the last 100 bases of the gene.
            # Example: for a gene of length 450, windows go from base 1-100, 101-200, 201-300, 301-400, and 351-450
            # Does mean we end up oversampling some bases, but seems to me to be better than the alternative of having
            # really short windows with potentially extremely biased GC content/depths
            gene_length = len(seq_index[gene_name].seq)
            windows = generate_windows(gene_length, window_size=100)
            for window in windows:
                total_depth = 0
                bases = 0
                gc_bases = 0
                start_base, end_base = window
                for column in bamfile.pileup(contig=gene_name, start=start_base, end=end_base,
                                             stepper='samtools', ignore_orphans=False,
                                             fastafile=pysam.FastaFile(gene_fasta),
                                             min_base_quality=0):
                    ref_base = seq_index[gene_name].seq[column.reference_pos]
                    if ref_base.upper() == 'G' or ref_base.upper() == 'C':
                        gc_bases += 1
                    bases += 1
                    total_depth += column.nsegments
                try:
                    depth = total_depth/bases
                    gc_content = round(100*gc_bases/bases)
                    if gc_content in gc_dict:
                        gc_dict[gc_content].append(depth)
                    else:
                        gc_dict[gc_content] = [depth]
                    depths.append(depth)
                except ZeroDivisionError:
                    pass
        bamfile.close()
        if out_bam is not None:
            try:
                shutil.copy(sorted_bam, out_bam)
            except shutil.SameFileError:
                pass
    return depths, gc_dict


def generate_windows(gene_length, window_size=100):
    windows = list()  # This will be a list of tuples, where each tuple will be start coordinate and end coordinate
    if gene_length < window_size:
        return windows
    current_position = 0
    while current_position < gene_length:
        windows.append((current_position, current_position + window_size - 1))
        current_position += window_size
    windows.pop()
    # Pysam uses zero based coordinates, so need the -1 so we don't overshoot end of gene
    windows.append((gene_length - window_size - 1, gene_length - 1))
    return windows


def check_dependencies():
    logging.info('Checking for external dependencies.')
    all_dependencies_found = True
    dependencies = ['skesa', 'bowtie2', 'bowtie2-build', 'diamond', 'samtools']
    for dependency in dependencies:
        if shutil.which(dependency) is None:
            all_dependencies_found = False
            logging.warning('WARNING: {} not found.'.format(dependency))
        else:
            logging.info('Checking for {}... Found!'.format(dependency))
    return all_dependencies_found


if __name__ == '__main__':
    main()
