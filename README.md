# MRF: Missing Regions Finder

MRF is a virus comparative genomics tool. It takes a query fasta, reference fasta, reference gff3 as inputs and lists the missing coding sequences (complete and partial) of query genome with respect to the annotations of the reference genome. MRF can be run in simple mode (one vs one) in command line and web interface or batch mode (one vs many) in command line.

To use a web version of this tool, please visit http://bioinfo.ciba.res.in/mrf

**Dependencies**

There are very minimal dependencies for running MRF-batch program. Verify if you have the following and install accordingly. Alternatively one can create a conda environment with requisite packages (scroll below)

* Perl v5
* Bash
* R > v3.6 and the following packages.
	* RColorBrewer_1.1-2 
	* tidyr_1.1.2        
	* plyr_1.8.4        
	* gplots_3.0.1.1
	* ggplot2_3.1.1      
	* dplyr_1.0.2

**Installation**

If you meet the above dependencies, clone the repository, cd into the directory and start running the programs.

**Overview**

**simple mode** (This usage is equivalent to the web version of MRF)
To compare and list the missing coding sequences in a query genome with respect to a reference genome. 

	bash run-mrf-wrapper.sh -q query fasta -r reference fasta -f reference gff3

	The above command generates several output files, out of which mrfOUT.mum.mcr contains the missing coding sequences and mrfOUT.mum.mr contains the missing regions.
	The program also creates a folder called circos_files in which the configuration and data files for generating a circos plot are present

To compare and list the missing genomic regions in a query genome with respect to a reference genome.

	bash run-mrf-wrapper.sh -q query fasta -r reference fasta

	The above command without the annotation file (GFF3), will just generate the missing regions.

To print usage

	bash run-mrf-wrapper.sh -h

	The above command prints help text as below.

	Usage:  run-mrf-wrapper.sh -q query fasta -r reference fasta -f reference gff3

	options: -m      mummer exact match length [default: 20]
         	 -l      False match length threshold.Values below this will be screened and confirmed if they are true matches [default: 15]
         	 -n      Negative offset. This option requires False match length as mandatory argument [default: 1]
         	 -p      Positive offset. This option requires False match length as mandatory argument [default: 1]
         	 -c      Partial coding sequences below this set threshold are not shown [default: 0]
         	 -o      output prefix [default: mrfOUT]
         	 -h      prints usage

**batch mode**

Use Cases:

	Use case 1: To identify deleted coding sequences in highly similar genomes with significant genome length difference and inconsistent annotations (Ex: White Spot Syndrome Virus)

	Use case 2: To tabulate deletions in virulence-related genes among DNA virus genomes (Ex. African Swine Fever Virus)

	Use case 3: To tabulate deletions in virulence-related genes in highly mutable RNA virus genomes (Ex. Human Immunodeficiency Virus)

	Use case 4: To detect the presence of one or more cloned fragments in a genome of interest (Ex. Shrimp Hemocyte Iridescent Virus)
