from utils import job_runner

def get_mapstats():
    # mainly samtools flagstat
    flagstat_cmd = ""

def get_library_complexity():
#### Compute library complexity
# Tab-delimited output:
# [1] TotalReadPairss
# [2] DistinctReadPairs
# [3] OneReadPairs
# [4] TwoReadPairs
# [5] NRF=Distinct/Total
# [6] PBC1=OnePair/Distinct
# [7] PBC2=OnePair/TwoPair
PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"
samtools sort -n ${FILT_BAM_FILE} ${FILT_BAM_FILE}_tmp
mv ${FILT_BAM_FILE}_tmp.bam ${FILT_BAM_FILE}
echo 1 | awk '{print "#Total\tDistinct\tOne\tTwo\tNRF\tPBC1\tPBC2"}' > ${PBC_FILE_QC}
bedtools bamtobed -bedpe -i ${FILT_BAM_FILE} | \
awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | \
awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}
     ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1}
     END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' >> ${PBC_FILE_QC}


def clean():
    # remove unnecessary files
    rm_cmd = ""

    #### Clean up
    rm ${FILT_BAM}
    mv ${FINAL_BEDPE_FILE} "${OFPREFIX}.bedpe.gz"
    mv ${FINAL_BAM_FILE} "${OFPREFIX}.bam"
    mv ${FINAL_BAM_INDEX_FILE} "${OFPREFIX}.bai"
    mv ${PBC_FILE_QC} "${OFPREFIX}".QC.pbc
    mv ${DUP_FILE_QC} "${OFPREFIX}".QC.markdups
    mv ${RAW_BAM_FILE_MAPSTATS} "${OFPREFIX}".QC.raw.mapstats
    mv ${FINAL_BAM_FILE_MAPSTATS} "${OFPREFIX}".QC.final.mapstats

