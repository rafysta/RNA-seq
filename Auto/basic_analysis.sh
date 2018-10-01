#!/bin/bash
# RNA-seq pipeline


get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-x, --organism
		organism. {human, mouse, pombe}

	-d, --directory
		data directory having target and control

	-n, --name
		sample name

	-p, --paired
		paired end or not (default : FALSE)

	-o, --log
		log directory

EOF
}

get_version(){
	echo "${0} version 1.0"
}


SHORT=hvx:d:n:p:o:
LONG=help,version,organism:,directory:,name:,paired:,log:
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
        -p|--paired)
            PAIRED="$2"
            shift 2
            ;;
        -o|--log)
            DIR_LOG="$2"
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
PAIRED=${PAIRED:-FALSE}

TIME_STAMP=$(date +"%Y-%m-%d")
DIR_LIB=$(dirname $0)


case $ORGANISM in
	pombe) BOWTIE_TARGET=ASM294v2.19
		   CHROM_SIZE=/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v1.18/LENGTH.txt ;;
	human) BOWTIE_TARGET=hg19
		   STAR_REF=/wistar/noma/Data/Human_seq/hg19/RSEM_star/hg19
		   CHROM_SIZE=/wistar/noma/Data/Human_seq/hg19/LENGTH.txt ;;
	mouse) BOWTIE_TARGET=mm10
		   CHROM_SIZE=/wistar/noma/Data/Mouse_seq/mm10/LENGTH.txt ;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac



	sbatch -n 8  --job-name=rsem_${NAME} -o "${DIR_LOG}/${TIME_STAMP}_alignment_by_rsem_${NAME}.log" --open-mode append <<EOF
#!/bin/sh
export PATH=/cm/shared/wistar/RSEM/RSEM-1.2.31:$PATH
cd ${DIR_DATA}
[ ! -e ${NAME} ] && mkdir ${NAME}

# single read
[ $PAIRED = FALSE ] && rsem-calculate-expression -p 8 --star --star-path /applications/STAR/current --estimate-rspd --keep-intermediate-files --output-genome-bam  ${NAME}.fastq ${STAR_REF} ${NAME}/${NAME}

# STAR alignment paired-end
[ $PAIRED = TRUE ] && rsem-calculate-expression --paired-end -p 8 --star --star-path /applications/STAR/current --estimate-rspd --keep-intermediate-files --output-genome-bam ${NAME}_1.fastq ${NAME}_2.fastq ${DIR_contig}/RSEM_star/hg19 ${NAME}/${NAME}

EOF






