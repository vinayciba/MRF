#!/bin/bash

MUM_LEN=20
CDS_PROP_THRESH=0
FALSE_MATCH_THRESH=15
NEGOFF=1
POSOFF=1
OUTPUT_PREFIX=mrfOUT
usage()

{
echo "Usage:	$0 -q query fasta -r reference fasta -f reference gff3

options: -m	 mummer exact match length [default: 20]
	 -l 	 False match length threshold.Values below this will be screened and confirmed if they are true matches [default: 15]
	 -n 	 Negative offset. This option requires False match length as mandatory argument [default: 1]
	 -p      Positive offset. This option requires False match length as mandatory argument [default: 1]
	 -c 	 Partial coding sequences below this set threshold are not shown [default: 0]
	 -o 	 output prefix [default: mrfOUT]
	 -h 	 prints usage
	 " 1>&2; exit 1;
}

while getopts ":m:g:q:r:f:l:n:p:c:o:" name;
   do
     case "${name}" in
        m) MUM_LEN=${OPTARG} ;;           	# mummer mum length
        g) GENOME_SIZE=${OPTARG} ;;       	# genome size
        q) QUERY=${OPTARG} ;;     	  	# query genome
        r) REFERENCE=${OPTARG} ;;         	# reference genome
        f) GFF3_FILE=${OPTARG} ;;         	# gff3 file
	l) FALSE_MATCH_THRESH=${OPTARG} ;;	# False match length cutoff value
	n) NEGOFF=${OPTARG} ;; 			# Negative Offset
	p) POSOFF=${OPTARG} ;;			# Positive offset
	c) CDS_PROP_THRESH=${OPTARG} ;;		# CDS below this percentage will not be shown
	o) OUTPUT_PREFIX=${OPTARG} ;;	  	# prefix for the output file name
        *) usage ;;                       	# display usage and exit
     esac
   done
shift $((OPTIND-1))

if [ -z "${QUERY}" ] || [ -z "${REFERENCE}" ] ; then
  usage
fi
if [ ! -z "${QUERY}" ] & [ ! -z "${REFERENCE}" ] & [ -z "${GFF3_FILE}" ]; then


echo "GFF3 file not provided, so proceeding to generate only missing genomic regions. Please provide reference GFF3 file to generate missing coding sequences."
#echo $'\n';
echo "MUM_LEN = ${MUM_LEN}"
echo "QUERY = ${QUERY}"
echo "REFERENCE = ${REFERENCE}"
echo "GFF3_FILE = ${GFF3_FILE}"
echo "FALSE_MATCH_THRESH = $[FALSE_MATCH_THRESH]"
echo "NEGOFF = $[NEGOFF]"
echo "POSOFF = $[POSOFF]"
echo "CDS_PROP_THRESH = $[CDS_PROP_THRESH]"
echo $'\n';

mummer -mum -l $MUM_LEN $REFERENCE $QUERY > $OUTPUT_PREFIX.mum;

ref_gen_len=`awk 'NR>1{s+=length()}END{print s}' $REFERENCE`;


echo "Generating missing genomic regions..."

perl missing_regions_finder.pl --inputfile $OUTPUT_PREFIX.mum --genome_size $ref_gen_len;

echo $'\n';
echo "done"

else

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
mkdir circos_files;
mummer -mum -l $MUM_LEN $REFERENCE $QUERY > $OUTPUT_PREFIX.mum;

ref_gen_len=`awk 'NR>1{s+=length()}END{print s}' $REFERENCE`;

echo $'\n';

echo "Generating missing regions and missing coding regions..."

echo $'\n';

perl missing_regions_finder.pl --inputfile  $OUTPUT_PREFIX.mum --genome_size $ref_gen_len --gff3 $GFF3_FILE --cds_threshold $CDS_PROP_THRESH --fmatch_len $FALSE_MATCH_THRESH --negoff $NEGOFF --posoff $POSOFF;
echo "done"

fi


