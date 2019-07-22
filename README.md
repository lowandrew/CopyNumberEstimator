# Copy Number Estimator

This is a still fairly crude attempt to estimate gene copy numbers in Illumina reads. Definitely not near production
ready yet.

### Installation

It's assumed you're using python 3.5+.

Clone this repository and install as a package using `setup.py` (you should do this in a virtualenv.)

You'll need the following external programs available somewhere in your $PATH. Other versions may work, I haven't tested
anything thoroughly:

- skesa v2.3.0
- bowtie2 v2.3.5
- diamond v0.9.24 
- samtools v1.9

### Usage

```
usage: copy_number_estimator [-h] -1 FORWARD_READS -2 REVERSE_READS -f
                             GENE_FASTA_FILE [-a ASSEMBLY] [-p P_VALUE]
                             [-t THREADS] [--database DATABASE]
                             [-s SINGLE_COPY_GENE_FOLDER] [-o OUTPUT_REPORT]

optional arguments:
  -h, --help            show this help message and exit
  -1 FORWARD_READS, --forward_reads FORWARD_READS
                        Path to forward reads.
  -2 REVERSE_READS, --reverse_reads REVERSE_READS
                        Path to forward reads.
  -f GENE_FASTA_FILE, --gene_fasta_file GENE_FASTA_FILE
                        Path to FASTA-formatted file containing your gene (or
                        genes?) of interest.
  -a ASSEMBLY, --assembly ASSEMBLY
                        If you have an assembly for your genome of interest,
                        specify it here. If you do not have an assembly, a
                        quick and dirty assembly will be created with SKESA.
  -p P_VALUE, --p_value P_VALUE
                        P-value to report on.
  -t THREADS, --threads THREADS
                        Number of threads to run analysis on. NOT YET
                        IMPLEMENTED.
  --database DATABASE   Database to store some stuff in.
  -s SINGLE_COPY_GENE_FOLDER, --single_copy_gene_folder SINGLE_COPY_GENE_FOLDER
                        Path to folder that contains single copy genes. Must
                        have a .fasta extension.
  -o OUTPUT_REPORT, --output_report OUTPUT_REPORT
                        CSV file to store your output in.

```
