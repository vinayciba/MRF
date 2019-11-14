# MRF: Missing Regions Finder

MRF compares two genomes, a query and a reference, and lists the missing coding sequences (complete and partial) of query genome with respect to the annotations of the reference genome.

To use a web version of this tool, please visit bioinfo.ciba.res.in/mrf

Prerequisites

    Bash
    GNU Awk
    Perl
    
Installation

	There is no need for installation.Clone the repo, which consists of one bash wrapper script run-mrf-wrapper.sh that calls individual scripts.

Overview

To compare and list the missing coding sequences in a query genome with respect to a reference genome.

	bash run-mrf-wrapper.sh -q query fasta -r reference fasta -f reference gff3

	The above command generates two output files, pre_mr.txt and pre_mcr.txt.

To compare and list the missing genomic regions in a query genome with respect to a reference genome.

	bash run-mrf-wrapper.sh -q query fasta -r reference fasta

	The above command generates a output file called pre_mr.txt.

To print usage

	bash run-mrf-wrapper.sh -h

	The above command prints help text as below.

	Usage:  run-mrf-wrapper.sh -q query fasta -r reference fasta -f reference gff3

	options: -m      mummer exact match length [default: 20]
		 -o      output prefix [default: pre]
		 -h      prints usage

Use Cases:

	Use case 1: To identify deleted coding sequences in highly similar genomes with significant genome length difference and inconsistent annotations (Ex: White Spot Syndrome Virus)

	Use case 2: To tabulate deletions in virulence-related genes among DNA virus genomes (Ex. African Swine Fever Virus)

	Use case 3: To tabulate deletions in virulence-related genes in highly mutable RNA virus genomes (Ex. Human Immunodeficiency Virus)

	Use case 4: To detect the presence of one or more cloned fragments in a genome of interest (Ex. Shrimp Hemocyte Iridescent Virus)
