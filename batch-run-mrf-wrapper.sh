 #!/bin/bash

###This is a wrapper script which calls other programs and passes inputs to them

MUM_LEN=20
CDS_PROP_THRESH=0
FALSE_MATCH_THRESH=15
NEGOFF=1
POSOFF=1
OUTPUT_PREFIX=mrfOUT
TOP=10
#FIRST=10
#LAST=10
#CDS_NUM=15
usage()

{
echo "Usage:    $0 -d query_fasta_directory -r reference_fasta -f reference_gff3


options: -m	 mummer exact match length [default: 20]
	 -l 	 False match length threshold.Values below this will be screened and confirmed if they are true matches [default: 15]
	 -n 	 Negative offset. This option requires False match length as mandatory argument [default: 1]
	 -p      Positive offset. This option requires False match length as mandatory argument [default: 1]
	 -c 	 Partial coding sequences below this set threshold are not shown [default: 0]
	 -o 	 output prefix [default: mrfOUT]

Output parsing options:

 	 -I      Search by protein ids [provide comma separated]
         -N      Search by protein names [provide comma separated]
         -Y      Filter by missing coding sequences length [default: 0]
         -F      List first n proteins [default:10]
         -Z      List last n proteins [default:10]
         -R      List by range[m,n]
	 -U	 List only top n affected genomes [default: 10]
         -A      List all proteins

Plot options:

	-L      Number of missing coding sequences to show in heatmap [default: 15]
	-B	Number of genomes to show in barplot [default: 50]

	-h 	prints usage
	 " 1>&2; exit 1;
}

while getopts "AF:Z:I:N:Y:L:B:R:U:m:g:q:d:r:f:l:n:p:c:o:" name;
   do
     case "${name}" in
        m) MUM_LEN=${OPTARG} ;;           	# mummer mum length
        g) GENOME_SIZE=${OPTARG} ;;       	# genome size
        q) QUERY ;;     	  	# query genome
	d) QDIR=${OPTARG} ;;                   # query genome directory
        r) REFERENCE=${OPTARG} ;;         	# reference genome
        f) GFF3_FILE=${OPTARG} ;;         	# gff3 file
	l) FALSE_MATCH_THRESH=${OPTARG} ;;	# False match length cutoff value
	n) NEGOFF=${OPTARG} ;; 			# Negative Offset
	p) POSOFF=${OPTARG} ;;			# Positive offset
	c) CDS_PROP_THRESH=${OPTARG} ;;		# CDS below this percentage will not be shown
	o) OUTPUT_PREFIX=${OPTARG} ;;	  	# prefix for the output file name
	I) P_ID=${OPTARG} ;;			# Search by protein id
	N) P_NAME=${OPTARG} ;;                  # Search by protein name
	Y) FILTER=${OPTARG} ;;			# Filter by missing coding sequence length
	L) CDS_NUM=${OPTARG} ;;		# Number of CDS to show in heatmap
	B) NUM_BARS=${OPTARG} ;;         # Number of genomes to show in barplot
	R) RANGE=${OPTARG} ;;			# List by range
	F) FIRST=${OPTARG} ;;				# List first n proteins
	Z) LAST=${OPTARG} ;;				# List last n proteins
	U) TOP=${OPTARG} ;;			# List top n affected genomes
	A) ALL=all;;				# List all proteins
        *) usage ;;                       	# display usage and exit
     esac
   done
shift $((OPTIND-1))

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" #get source location of this script


if [ -z "${QDIR}" ] || [ -z "${REFERENCE}" ] || [ -z "${GFF3_FILE}" ]; then
  usage
else

	for i in $QDIR/*;
		do
		if [[ $i =~ .*\.(fasta|fa|fna) ]];
			then
		QUERY=$i;
		Q=${QUERY##*/}
		Q=${Q%.*}
		R=${REFERENCE##*/}
		R=${R%.*}
		OUT_PREFIX=${OUTPUT_PREFIX}_${Q}_${R}

		queryName=`grep ">" $QUERY | awk -F '>' '{print $2}'`
		referenceName=`grep ">" $REFERENCE | awk -F '>' '{print $2}'`

		echo "MUM_LEN = ${MUM_LEN}"
#echo "GENOME_SIZE = ${GENOME_SIZE}"
		echo "QUERY = ${QUERY}"
		echo "REFERENCE = ${REFERENCE}"
		echo "GFF3_FILE = ${GFF3_FILE}"
		echo "FALSE_MATCH_THRESH = $[FALSE_MATCH_THRESH]"
		echo "NEGOFF = $[NEGOFF]"
		echo "POSOFF = $[POSOFF]"
		echo "CDS_PROP_THRESH = $[CDS_PROP_THRESH]"

		echo $'\n';
#		mkdir circos_files;
		mummer -mum -l $MUM_LEN $REFERENCE $QUERY > $OUT_PREFIX.mum;

		ref_gen_len=`awk 'NR>1{s+=length()}END{print s}' $REFERENCE`; #size of the reference genome

		echo $'\n';

		echo "Generating missing regions and missing coding regions..."
		echo $'\n';

		perl $DIR/MRF_batch.pl --inputfile  $OUT_PREFIX.mum --genome_size $ref_gen_len --gff3 $GFF3_FILE --cds_threshold $CDS_PROP_THRESH --fmatch_len $FALSE_MATCH_THRESH --negoff $NEGOFF --posoff $POSOFF --queryName $queryName --referenceName $referenceName;
		echo "done"
		else echo $i "Not a valid format";
		fi;
	done;
fi

echo "Generated missing ooding sequences for all.";
echo $'\n';


if [ ! -z "${FIRST}" ];
	then
echo "showing first n proteins ..";
	perl $DIR/parse_mrf_output.pl --top $TOP --head $FIRST;
elif [ ! -z "${LAST}" ];
        then
echo "showing  last n proteins ..";
        perl $DIR/parse_mrf_output.pl --top $TOP --tail $LAST;
elif [ ! -z "${RANGE}" ];
        then
echo "filtering by range ..";
        perl $DIR/parse_mrf_output.pl --top $TOP --range $RANGE;
elif [ ! -z "${P_ID}" ];
        then
echo "filtering by protein id ..";
        perl $DIR/parse_mrf_output.pl --top $TOP --prot_ID $P_ID;
elif [ ! -z "${P_NAME}" ];
        then
echo "filtering by protein name ..";
        perl $DIR/parse_mrf_output.pl --top $TOP --prot_name $P_NAME;
elif [ ! -z "${CDS_NUM}" ];
        then
echo "showing top n number of proteins ..";
        perl $DIR/parse_mrf_output.pl --top $TOP --number $CDS_NUM;
elif [ ! -z "${FILTER}" ];
        then
echo "filter by CDS length ..";
        perl $DIR/parse_mrf_output.pl --top $TOP --filter $FILTER;
elif [ ! -z "${ALL}" ];
        then
echo "showing all proteins ..";
        perl $DIR/parse_mrf_output.pl --all;
else
	perl $DIR/parse_mrf_output.pl --top $TOP;

fi;


if [ ! -z "${NUM_BARS}" ];
        then

Rscript $DIR/genPlots.R $NUM_BARS;

else

Rscript $DIR/genPlots.R

fi;
