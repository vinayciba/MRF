# MRF: Missing Regions Finder

MRF is a virus comparative genomics tool. It takes a query fasta, reference fasta, reference gff3 as inputs and lists the missing coding sequences (complete and partial) of query genome with respect to the annotations of the reference genome. MRF can be run in single mode (one vs one) in command line and web interface or batch mode (one vs many) in command line.

To use a web version of this tool, please visit http://14.139.181.163/mrf

**Dependencies**

There are very minimal dependencies for running MRF-batch program. Verify if you have the following and install accordingly. Alternatively one can create a conda environment with requisite packages (scroll below)

* Perl v5
* Bash
* Mummer
* R > v3.6 and the following packages.
	* RColorBrewer_1.1-2 
	* tidyr_1.1.2        
	* plyr_1.8.4        
	* gplots_3.0.1.1
	* ggplot2_3.1.1      
	* dplyr_1.0.2

**Installation**

If you meet the above dependencies, clone the repository, cd into the directory and start running the programs.Otherwise you can intall through conda as below.

* **Through Conda**

Clone this repository and follow the steps below

    cd MRF
    conda env create -f MRF.conda_env.yml
    source activate MRF.env



**Overview**

**single mode** (This usage is equivalent to the web version of MRF)

To compare and list the missing coding sequences in a query genome with respect to a reference genome. 

	bash run-mrf-wrapper.sh -q query fasta -r reference fasta -f reference gff3

The above command generates several output files, out of which mrfOUT.mum.mcr contains the missing coding sequences and mrfOUT.mum.mr contains the missing regions.The program also creates a folder called circos_files in which the configuration and data files for generating a circos plot are present

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

MRF can be used in batch mode (command line usage only) to analyze multiple query genomes against a single reference genome.

To see the available options, run the following command.

    bash batch-run-mrf-wrapper.sh -h

```
Usage:    batch-run-mrf-wrapper.sh -d query_fasta_directory -r reference_fasta -f reference_gff3


options: -m      mummer exact match length [default: 20]
         -l      False match length threshold.Values below this will be screened and confirmed if they are true matches [default: 15]
         -n      Negative offset. This option requires False match length as mandatory argument [default: 1]
         -p      Positive offset. This option requires False match length as mandatory argument [default: 1]
         -c      Partial coding sequences below this set threshold are not shown [default: 0]
         -o      output prefix [default: mrfOUT]

Output parsing options:

         -I      Search by protein ids [provide comma separated]
         -N      Search by protein names [provide comma separated]
         -Y      Filter by missing coding sequences length [default: 0]
         -F      List first n proteins [default:10]
         -Z      List last n proteins [default:10]
         -R      List by range[m,n]
         -U      List only top n affected genomes [default: 10]
         -A      List all proteins

Plot options:

        -L      Number of missing coding sequences to show in heatmap [default: 15]
        -B      Number of genomes to show in barplot [default: 50]

        -h      prints usage

```
To analyze the genomes present in a directory called *query* against a reference genome, run the following command.

    bash batch-run-mrf-wrapper.sh -d query -r Reference.fasta -f Reference.gff3
    
The above command produces three files called summaryAll.txt, cdsHeatmap.txt and cdsAll.txt. Apart from these three files, four easy to interpret figures,  a Scatterplot, a Barplot and two Heatmaps (clustered and unclustered) are produced.

**For detailed explanation along with the use cases refer to the manual**

