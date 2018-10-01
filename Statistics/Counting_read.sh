#!/bin/bash
# Counting read number from SAM file


get_usage(){
	cat <<EOF

Usage : $0 [OPTION] -i [sam file] -o [output file]

Description
	-h, --help
		show help

	-i, --input
		SAM file

	-o, --output
		output text file

EOF
}


SHORT=hi:o:
LONG=help,input:,output:
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
		-i|--input)
            FILE_SAM="$2"
            shift 2
            ;;
		-o|--output)
            FILE_OUT="$2"
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

if [ ! -e $FILE_SAM ]; then
	echo "$FILE_SAM is not exists"
	exit 1
fi


NO_ALIGN=$(samtools view -c -f 4 $FILE_SAM)
ALIGNED=$(samtools view -c -F 4 $FILE_SAM)
UNIQUE_ALIGN=$(samtools view -c -F 4 -q 10 $FILE_SAM)
let TOTAL=${NO_ALIGN}+${ALIGNED}
let LOW_QUALITY=${ALIGNED}-${UNIQUE_ALIGN}

NAME=$(basename ${FILE_SAM})
NAME=${NAME%.sam}

echo -e "name\ttotal\tno_aligned\tlow_quality\tunique_aligned" > ${FILE_OUT}
echo -e "${NAME}\t${TOTAL}\t${NO_ALIGN}\t${LOW_QUALITY}\t${UNIQUE_ALIGN}" >> ${FILE_OUT}



