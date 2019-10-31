#!/bin/bash
# RNA-seq pipeline


get_usage(){
	cat <<EOF

Usage : $0 [OPTION] [fastq files]

Description
	-h, --help
		show help

	-v, --version
		show version

	-x, --organism
		organism. {human, mouse, pombe, fly}

	-d, --directory
		data directory having target and control

	-n, --name
		sample name

	-L, --layout [single / pair]
		paired end or single (default : single)

	-o, --log
		log directory
	
	-t, --thread [thread number]
		specify thread number. default : 1
EOF
}

get_version(){
	echo "${0} version 1.0"
}


SHORT=hvx:d:n:L:o:t:
LONG=help,version,organism:,directory:,name:,layout:,log:,thread:
PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [[ $? -ne 0 ]]; then
    exit 2
fi
eval set -- "$PARSED"

while true; do
    case "$1" in
        -h|--help)
            get_usage
            exit 1
            ;;
        -v|--version)
            get_version
            exit 1
            ;;
        -x|--organism)
            ORGANISM="$2"
            shift 2
            ;;
        -d|--directory)
            DIR_DATA="$2"
            shift 2
            ;;
        -n|--name)
            NAME="$2"
            shift 2
            ;;
        -L|--layout)
            LAYOUT="$2"
            shift 2
            ;;
        -o|--log)
            DIR_LOG="$2"
            shift 2
            ;;
        -t|--thread)
            THREAD="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done

[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${NAME}" ] && echo "Please specify sample NAME" && exit 1
[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1
[ ! -n "${DIR_LOG}" ] && echo "Please specify log directory" && exit 1
LAYOUT=${PAIRED:-single}
THREAD=${THREAD:-1}

TIME_STAMP=$(date +"%Y-%m-%d")
DIR_LIB=$(dirname $0)
INPUT_FILES=$@

case $ORGANISM in
	fly)
		CHROM_SIZE=${HOME}/Genome/data/drosophila_melanogaster/BDGP6.22/LENGTH.txt
		STAR_DIR=${HOME}/Genome/data/drosophila_melanogaster/BDGP6.22/STAR
		RSEM_REF=${HOME}/Genome/data/drosophila_melanogaster/BDGP6.22/RSEM/dros_0
		;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac

module load racs-eb 2> /dev/null
module load samtools
module load deepTools 2> /dev/null
module load miniconda
conda activate rsem

[ ! -e ${DIR_DATA}/${NAME} ] && mkdir ${DIR_DATA}/${NAME}
cd ${DIR_DATA}/${NAME}

### STAR alignment
STAR --runThreadN ${THREAD} --genomeDir ${STAR_DIR} --readFilesIn ${INPUT_FILES} --outSAMattributes MD --outFileNamePrefix ${NAME}_ --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate

## make index file for bam
samtools index ${NAME}_Aligned.sortedByCoord.out.bam

## bamCoverage
bamCoverage -p ${THREAD} -b ${NAME}_Aligned.sortedByCoord.out.bam -o ${NAME}_Aligned.sortedByCoord.out.bw

## RSEM
if [ $LAYOUT = "single" ]; then
	RSEM_MODE=""
else
	RSEM_MODE="--paired-end"
fi
rsem-calculate-expression -p ${THREAD} $RSEM_MODE --bam ${NAME}_Aligned.toTranscriptome.out.bam ${RSEM_REF} ${NAME}

## Cleaning
rm -r ${NAME}.stat
rm ${NAME}*.bam ${NAME}*.bai ${NAME}*.wig ${NAME}*.tab


