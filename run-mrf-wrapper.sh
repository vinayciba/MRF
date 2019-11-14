#!/bin/bash

MUM_LEN=20 
 OUTPUT_PREFIX=pre
usage() 

{
echo "Usage:	$0 -q query fasta -r reference fasta -f reference gff3

options: -m	 mummer exact match length [default: 20]
	 -o 	 output prefix [default: pre]
	 -h 	 prints usage
	 " 1>&2; exit 1;
}

while getopts ":m:g:q:r:f:o:" name;
   do
     case "${name}" in
        m) MUM_LEN=${OPTARG} ;;           # mummer mum length
        g) GENOME_SIZE=${OPTARG} ;;       # genome size
        q) QUERY=${OPTARG} ;;     	  # query genome
        r) REFERENCE=${OPTARG} ;;         # reference genome
        f) GFF3_FILE=${OPTARG} ;;         # gff3 file
	o) OUTPUT_PREFIX=${OPTARG} ;;	  # prefix for the output file name
        *) usage ;;                       # display usage and exit
     esac
   done
shift $((OPTIND-1))

if [ -z "${QUERY}" ] || [ -z "${REFERENCE}" ] ; then
  usage
fi
if [ ! -z "${QUERY}" ] & [ ! -z "${REFERENCE}" ] & [ -z "${GFF3_FILE}" ]; then

echo "GFF3 file not provided, so proceeding to generate only missing genomic regions. Please provide reference GFF3 file to generate missing coding sequences."
echo "MUM_LEN = ${MUM_LEN}"
echo "QUERY = ${QUERY}"
echo "REFERENCE = ${REFERENCE}"
echo "GFF3_FILE = ${GFF3_FILE}"

mummer -mum -l $MUM_LEN $REFERENCE $QUERY > mum_new_out;
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' mum_new_out |  awk 'NR > 1' > mum_new_formatted;

ref_gen_len=`awk 'NR>1{s+=length()}END{print s}' $REFERENCE`;
#echo $'\n';
echo "Generating missing genomic regions..."
perl generate_missing_genomic_regions.pl mum_new_formatted $ref_gen_len > $OUTPUT_PREFIX"_mr.txt";
echo "done"

else

echo "MUM_LEN = ${MUM_LEN}"
#echo "GENOME_SIZE = ${GENOME_SIZE}"
echo "QUERY = ${QUERY}"
echo "REFERENCE = ${REFERENCE}"
echo "GFF3_FILE = ${GFF3_FILE}"


mummer -mum -l $MUM_LEN $REFERENCE $QUERY > mum_new_out;

awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' mum_new_out |  awk 'NR > 1' > mum_new_formatted;

ref_gen_len=`awk 'NR>1{s+=length()}END{print s}' $REFERENCE`;

#echo $ref_gen_len
#echo $'\n';

echo "Generating missing regions and missing coding regions..."
perl generate_missing_genomic_regions.pl mum_new_formatted $ref_gen_len >  $OUTPUT_PREFIX"_mr.txt";
perl generate_missing_coding_sequences.pl mum_new_formatted $GFF3_FILE $ref_gen_len >  $OUTPUT_PREFIX"_mcr.txt";
echo "done"

fi


rm mum_new_out mum_new_formatted
rm featuretable.gff3
