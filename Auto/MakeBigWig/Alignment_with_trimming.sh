#!/bin/bash -l


#==============================================================
# export following parameters
#==============================================================
#ORGANISM
#NAME
#DIR_DATA
#DIR_LIB

TIME_STAMP=$(date +"%Y-%m-%d")
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${NAME}" ] && echo "Please specify NAME" && exit 1
[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1
[ ! -n "${DIR_LIB}" ] && echo "Please specify DIR_LIB" && exit 1

PROGRAM_CIGAR=${DIR_LIB}/CigarFilter.pl
[ ! -e $PROGRAM_CIGAR ] && echo "CigarFilter.pl is not found in $PROGRAM_CIGAR" && exit 1

case $ORGANISM in
	pombe) BOWTIE_TARGET=ASM294v2.19
	       export BOWTIE2_INDEXES=/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v1.18 ;;
	human) BOWTIE_TARGET=hg19
	       export BOWTIE2_INDEXES=/wistar/noma/Data/Human_seq/hg19 ;;
	mouse) BOWTIE_TARGET=mm10
	       export BOWTIE2_INDEXES=/wistar/noma/Data/Mouse_seq/mm10 ;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac


cd ${DIR_DATA}

echo "alignment of ${NAME}.fastq"
let READ_LENGTH=`head -n 2 ${NAME}.fastq | tail -n 1 | wc -m`
let MAX2_TRIM=$READ_LENGTH-25
let MAX_TRIM=$READ_LENGTH-20

echo "alignment of ${READ_LENGTH} read"
/applications/bowtie2/current/bowtie2 -x ${BOWTIE_TARGET} -U ${NAME}.fastq -q -p 12 --no-unal --un ${NAME}_unaligned.fastq > ${NAME}.sam
mv ${NAME}_unaligned.fastq ${NAME}_tmp.fastq

for LEN in `seq 5 5 $MAX2_TRIM`
do
	echo "trimming $LEN bp"
	/applications/bowtie2/current/bowtie2 -x ${BOWTIE_TARGET} -U ${NAME}_tmp.fastq -q -3 $LEN -p 12 --no-hd --no-unal --un ${NAME}_unaligned.fastq >> ${NAME}.sam
	mv ${NAME}_unaligned.fastq ${NAME}_tmp.fastq
	/applications/bowtie2/current/bowtie2 -x ${BOWTIE_TARGET} -U ${NAME}_tmp.fastq -q -5 $LEN -p 12 --no-hd --no-unal --un ${NAME}_unaligned.fastq >> ${NAME}.sam
	mv ${NAME}_unaligned.fastq ${NAME}_tmp.fastq
done

echo "trimming $MAX_TRIM bp"
/applications/bowtie2/current/bowtie2 -x ${BOWTIE_TARGET} -U ${NAME}_tmp.fastq -q -3 ${MAX_TRIM} -p 12 --no-hd --no-unal --un ${NAME}_unaligned.fastq >> ${NAME}.sam
mv ${NAME}_unaligned.fastq ${NAME}_tmp.fastq
/applications/bowtie2/current/bowtie2 -x ${BOWTIE_TARGET} -U ${NAME}_tmp.fastq -q -5 ${MAX_TRIM} -p 12 --no-hd >> ${NAME}.sam
rm ${NAME}_tmp.fastq

echo "alignment of partial alignment reads"
perl ${PROGRAM_CIGAR} -s ${NAME}.sam -o ${NAME}_fastqList.txt > ${NAME}_tmp.sam
mv ${NAME}_tmp.sam ${NAME}.sam
LIST_NUM=`cat ${NAME}_fastqList.txt | wc -l`
for i in `seq 1 ${LIST_NUM}`
do
	TRIM_OPTION=`head -n $i ${NAME}_fastqList.txt | tail -n 1 | cut -f1`
	FASTQ_REANALYZE=`head -n $i ${NAME}_fastqList.txt | tail -n 1 | cut -f2`
	echo "Command : /applications/bowtie2/current/bowtie2 -x ${BOWTIE_TARGET} -U ${FASTQ_REANALYZE} -q ${TRIM_OPTION} -p 12 --no-hd"
	/applications/bowtie2/current/bowtie2 -x ${BOWTIE_TARGET} -U ${FASTQ_REANALYZE} -q ${TRIM_OPTION} -p 12 --no-hd >> ${NAME}.sam
	rm ${FASTQ_REANALYZE}
done
rm ${NAME}_fastqList.txt

