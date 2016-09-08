#!/bin/bash

# Copyright 2015 Junko Tsuji

# Compute strand cross-correlation and fragment length
# with 'phantompeakqualtools'.

#### Software path
export PATH=$PATH:/home/tsujij/soft/samtools-1.0/
export PATH=$PATH:/home/tsujij/soft/bedtools2-2.20.1/bin/
RSCRIPT="/home/tsujij/soft/R-3.1.0/bin/Rscript"
RUN_SPP_NODUPS="/home/tsujij/soft/phantompeakqualtools/run_spp_nodups.R"

#### Usage
function usage {
cat <<EOF
Usage: $0 [options] -f <se | pe> -i <bam>

Compute strand cross-correlation and fragment length.

Required arguments:
  -f <se | pe>   Library type.
                   - 'se': single-end
                   - 'pe': paired-end
  -i <bam>       Input BAM file.

Options:
  -h             Show this help message and exit.
  -p <cpu-num>   Number of CPUs (default: 8 cores).
  -o <prefix>    Prefix for output files (default: prefix of <bam>).
EOF
}

#### Parameters
NTHREADS=8
NREADS=15000000
unset SEQ_TYPE
unset BAM_FILE
unset OFPREFIX

#### Read arguments and options
[[ $# -eq 0 ]] && usage && exit 1;
while getopts "hi:f:p:o:" OPT
do
  case ${OPT} in
    "i") BAM_FILE=${OPTARG} ;;
    "f") SEQ_TYPE=${OPTARG} ;;
    "p") NTHREADS=${OPTARG} ;;
    "o") OFPREFIX=${OPTARG} ;;
    "h") usage && exit 1    ;;
     * ) usage && exit 1    ;;
  esac
done
[[ $# -lt 4 ]] && echo "$0: need more arguments" && exit 1;

#### Check arguments
if [ "${SEQ_TYPE}" != "se" ] && [ "${SEQ_TYPE}" != "pe" ]
  then echo "$0: input 'se' or 'pe'" && exit 1; fi
if [ ! -f "${BAM_FILE}" ]
  then echo "$0: can't read BAM" && exit 1; fi
if [ "${SEQ_TYPE}" == "pe" ]
then
  READ_1=`samtools view -f64 ${BAM_FILE} | head -n 100 | wc -l`
  if [ ${READ_1} -eq 0 ]
    then echo "$0: input BAM is not paired-end BAM" && exit 1; fi
fi
if [ $( echo ${NTHREADS} | sed 's/^[+0-9][0-9]*//' | wc -c ) -ne 1 ]
  then echo "$0: can't interpret -p ${NTHREADS}" && exit 1; fi
if [ -z "${OFPREFIX}" ]
  then OFPREFIX=`basename ${BAM_FILE} | awk -F ".bam" '{print $1}'`; fi

#### Create subsample tagAlign file
SUBSAMPLE_TA="${OFPREFIX}.subsample.tagAlign.gz"
if [ "${SEQ_TYPE}" == "se" ]
then
  bedtools bamtobed -i ${BAM_FILE} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | \
  grep -v "chrM" | shuf -n ${NREADS} | gzip -c > ${SUBSAMPLE_TA}
else
  samtools view -f64 -b ${BAM_FILE} | bedtools bamtobed -i stdin | \
  awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | grep -v "chrM" | \
  shuf -n ${NREADS} | gzip -c > ${SUBSAMPLE_TA}
fi

#### Calculate strand cross-correlation and fragment length
CC_SCORES_FILE="${OFPREFIX}.QC.cc"
CC_PLOT_FILE="${OFPREFIX}.QC.cc.pdf"
${RSCRIPT} ${RUN_SPP_NODUPS} -c=${SUBSAMPLE_TA} -p=${NTHREADS} -filtchr=chrM \
                             -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE}
echo 1 | awk 'BEGIN{OFS="\t"}
              {print "#Filename","numReads","estFragLen","corr_estFragLen",
                     "PhantomPeak","corr_phantomPeak","argmin_corr","min_corr",
                     "phantomPeakCoef","relPhantomPeakCoef","QualityTag"}' > ${CC_SCORES_FILE}.temp
sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} | \
awk -F ".subsample.tagAlign.gz" '{print $1".bam"$2}' >> ${CC_SCORES_FILE}.temp
mv ${CC_SCORES_FILE}.temp ${CC_SCORES_FILE}
rm ${SUBSAMPLE_TA}
