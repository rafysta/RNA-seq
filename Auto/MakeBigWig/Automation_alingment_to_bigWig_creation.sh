#!/bin/bash
# trimming alignment, homer, bigWigの作成まで自動化で行うプログラム


get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-n, --name
		sample name

	-x, --organism
		organism. {human, mouse, pombe}

	-m, --method
		STAR or BOWTIE with trimming method

	-d, --directory
		data directory

	-o, --log
		log directory

	-q, --fastqc
		TRUE for doing fastqc analysis or FALSE for not doing (default TRUE)
EOF
}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvn:x:m:d:o:q:
LONG=help,version,name:,organism:,method:,directory:,log:,fastqc:
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
        -n|--name)
            NAME="$2"
            shift 2
            ;;
        -x|--organism)
            ORGANISM="$2"
            shift 2
            ;;
        -m|--method)
            METHOD="$2"
            shift 2
            ;;
        -d|--directory)
            DIR_DATA="$2"
            shift 2
            ;;
        -o|--log)
            DIR_LOG="$2"
            shift 2
            ;;
        -q|--fastqc)
            FLAG_fastqc="$2"
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
[ ! -n "${NAME}" ] && echo "Please specify NAME" && exit 1
[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1
[ ! -n "${DIR_LOG}" ] && echo "Please specify log directory" && exit 1
FLAG_fastqc=${FLAG_fastqc:-TRUE}

TIME_STAMP=$(date +"%Y-%m-%d")
DIR_LIB=$(dirname $0)


# fastqcディレクトリが存在しなかったら作成
[ "$FLAG_fastqc" = "TRUE" ] && [ ! -e "${DIR_DATA}/fastqc" ] && mkdir "${DIR_DATA}/fastqc"


case $ORGANISM in
	pombe) BOWTIE_TARGET=ASM294v2.19
		   DIR_contig=/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v2.30
		   CHROM_SIZE=/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v1.18/LENGTH.txt
	       BOWTIE2_INDEXES=/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v1.18 ;;
	human) BOWTIE_TARGET=hg19
		   DIR_contig=/wistar/noma/Data/Human_seq/hg19
		   CHROM_SIZE=/wistar/noma/Data/Human_seq/hg19/LENGTH.txt
	       BOWTIE2_INDEXES=/wistar/noma/Data/Human_seq/hg19 ;;
	mouse) BOWTIE_TARGET=mm10
		   DIR_contig=/wistar/noma/Data/Mouse_seq/mm10
		   CHROM_SIZE=/wistar/noma/Data/Mouse_seq/mm10/LENGTH.txt
	       BOWTIE2_INDEXES=/wistar/noma/Data/Mouse_seq/mm10 ;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac

DIR_STAR=${DIR_contig}/STAR


[ "$FLAG_fastqc" = "TRUE" ] && [ ! -e ${DIR_DATA}/fastqc/${NAME}_fastqc ] && sbatch -n 12 --job-name=fastqc_${NAME} -o "${DIR_LOG}/${TIME_STAMP}_fastqc_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; /applications/fastqc/current/fastqc -o fastqc/ --nogroup -t 12 ${NAME}.fastq"



###########
# alignment
###########
case $METHOD in
	STAR)
		sbatch -n 12 --job-name=Align_${NAME} -o "${DIR_LOG}/${TIME_STAMP}_Alignment_by_STAR_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; /applications/STAR/current/STAR --runThreadN 12  --genomeDir $DIR_STAR --readFilesIn ${NAME}.fastq --outFileNamePrefix ${NAME} && mv ${NAME}Aligned.out.sam ${NAME}.sam; mv ${NAME}Log.out ${DIR_LOG}_Alignment_STAR_${NAME}.Log.out; mv ${NAME}Log.progress.out ${DIR_LOG}/${TIME_STAMP}_Alignment_STAR_${NAME}.Log.progress.out; mv ${NAME}Log.final.out ${DIR_LOG}/${TIME_STAMP}_Alignment_STAR_${NAME}.Log.final.out; rm ${NAME}SJ.out.tab" ;;
	BOWTIE)
		sbatch -n 12 --export=NAME=${NAME},ORGANISM=${ORGANISM},DIR_DATA=${DIR_DATA},DIR_LIB=${DIR_LIB} --job-name=Align_${NAME} -o "${DIR_LOG}/${TIME_STAMP}_Alignment_by_trimming_${NAME}.log" --open-mode append ${DIR_LIB}/Alignment_with_trimming.sh ;;
	*)
		echo "Unknown method"
		exit 1;;
esac


###########
# homer make tag directory
###########
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep "Align_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=Hom_${NAME} $DEPEND -o "${DIR_LOG}/${TIME_STAMP}_homer_${NAME}.log"  --open-mode append --wrap="cd ${DIR_DATA}; /applications/homer/current/bin/makeTagDirectory ${NAME} ${NAME}.sam -format sam && /applications/homer/current/bin/makeUCSCfile ${NAME} -o auto -name ${NAME} -fragLength 150 && cd ${NAME} && gunzip ${NAME}.ucsc.bedGraph.gz && /applications/bedgraphtobigwig/current/bedGraphToBigWig ${NAME}.ucsc.bedGraph $CHROM_SIZE ${NAME}.bw && gzip ${NAME}.ucsc.bedGraph"




