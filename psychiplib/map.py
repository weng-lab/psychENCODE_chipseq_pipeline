#!/usr/bin/env python

import os.path
import subprocess
from argparse import ArgumentParser

from utils import job_runner

def map(fq1, fq2, bwa_index, prefix, cpu_num):
    """Genome mapping function for ChIP-seq reads.

    """
    job = JobRunner()

    params = locals() # get all arguments as a dictionary

    ## mapping reads to the genome
    raw_sam = "{}.raw.sam.gz".format(prefix)
    params.setdefault({"raw_sam": raw_sam})
    bwa_cmd = \
      "bwa mem -M -t {cpu_num} {bwa_index} {fq1} {fq2} | gzip -c > {raw_sam}" 
    job.append([[ bwa_cmd.format(**params) ]])
    job.run()


    ## find bad CIGAR strings
    badcigar_file = "{}.badcigar".format(prefix)
    find_badcigar_awk_cmd = \
        """zcat {0} | \
           awk 'BEGIN{FS="\t"; OFS="\t"}
             !/^@/ && $6!="*" {
                 cigar=$6; gsub("[0-9]+D","",cigar);
                 n=split(cigar, vals, "[A-Z]");
                 s=0; for(i=1;i<=n;i++) s=s+vals[i];
                 seqlen=length($10);
                 if(s!=seqlen) print $1"\t";}' | \
           sort | uniq > {1}""".format(raw_sam, badcigar_file)
    job.append([[ find_badcigar_awk_cmd ]])
    job.run()

    ## remove bad CIGAR read pairs
    with open(badcigar_file) as badcigar_fp:
        peek = [x for i,x in enumerate(fp) if i < 10 and x.rstrip()]

    raw_bam = "{}.raw.srt.bam".format(prefix)
    badcigar_run_cmd = "{} samtooves view -Su {} | samtools sort - {}"
    if len(peek):
        badcigar_flt_cmd = \
            "zcat {} | grep -v -F -f {} |".format(raw_sam, badcigar_file)
        badcigar_run_cmd.format(badcigar_flt_cmd, "-", raw_bam[:-4])
    else:
        badcigar_run_cmd.format(" ", raw_sam, raw_bam[:-4])
    job.append([[ badcigar_run_cmd ]])
    job.run()

    ## QC part: there should be the same command somewhere
    ## We should make the command simpler
    raw_bam_mapstats = "{}.flagstat.qc".format(raw_bam[:-4])    
    qc_cmd = \
        "samtools flagstat {} > {}".format(raw_bam, raw_bam_mapstats)
    job.append([[ qc_cmd ]])
    # removing intermediate files... which one should we keep?
    # job.append([[ "rm {badcigar_file} {raw_sam}" ]])

    ## filter unmapped, mate unmapped, not primary alignment,
    ## and low quality reads
FILT_BAM_PREFIX="${OFPREFIX}.filt.srt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
TMP_FILT_BAM_PREFIX="${FILT_BAM_PREFIX}.tmp.nmsrt"
TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam"
TMP_FILT_BAM_FILE_CLEAN="${TMP_FILT_BAM_PREFIX}.clean.bam"
samtools view -F 1804 -f 2 -q ${MAPQ} -u ${RAW_BAM_FILE} | \
samtools sort -n - ${TMP_FILT_BAM_PREFIX}


#### Filter orphan reads (pair was removed) and read pairs mapping to different chromosomes
samtools fixmate -r ${TMP_FILT_BAM_FILE} ${TMP_FILT_BAM_FILE_CLEAN}
samtools view -F 1804 -f 2 -u ${TMP_FILT_BAM_FILE_CLEAN} | \
samtools sort - ${FILT_BAM_PREFIX}
rm ${TMP_FILT_BAM_FILE} ${TMP_FILT_BAM_FILE_CLEAN}


#### Mark duplicates with picard
TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc"
java -Xmx4G -jar ${PICARD_JAR} MarkDuplicates INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} \
METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}
rm ${RAW_BAM_FILE}


#### Remove duplicates and sponges, and create final name sorted BAM
FINAL_BAM_PREFIX="${OFPREFIX}.filt.srt.nodup"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam"
FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai"
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc"
FINAL_NMSRT_BAM_PREFIX="${OFPREFIX}.filt.nmsrt.nodup"
FINAL_NMSRT_BAM_FILE="${FINAL_NMSRT_BAM_PREFIX}.bam"
if [ ! -z "${EXCLUDE}" ]
then
    samtools view -F 1804 -f 2 -h ${FILT_BAM_FILE} | \
    grep -v -F -f ${EXCLUDE} - | samtools view -Su - > ${FINAL_BAM_FILE}
else
    samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}
fi
samtools sort -n ${FINAL_BAM_FILE} ${FINAL_NMSRT_BAM_PREFIX}

#### Index final sorted BAM
samtools index ${FINAL_BAM_FILE}
mv ${FINAL_BAM_FILE}.bai ${FINAL_BAM_INDEX_FILE}
samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}


#### Compute library complexity
# Tab-delimited output:
# [1] TotalReadPairs
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
rm ${FILT_BAM_FILE}



    BEDPE 2 files
    bedpe3.gz
    bedpe7.gz

#### Convert final name sorted BAM to BEDPE
FINAL_BEDPE_FILE="${FINAL_NMSRT_BAM_PREFIX}.bedpe.gz"
bedtools bamtobed -bedpe -mate1 -i ${FINAL_NMSRT_BAM_FILE} | gzip -c > ${FINAL_BEDPE_FILE}
rm ${FINAL_NMSRT_BAM_FILE}

#### Clean up
mv ${FINAL_BEDPE_FILE} "${OFPREFIX}.bedpe.gz"
mv ${FINAL_BAM_FILE} "${OFPREFIX}.bam"
mv ${FINAL_BAM_INDEX_FILE} "${OFPREFIX}.bai"
mv ${PBC_FILE_QC} "${OFPREFIX}".QC.pbc
mv ${DUP_FILE_QC} "${OFPREFIX}".QC.markdups
mv ${RAW_BAM_FILE_MAPSTATS} "${OFPREFIX}".QC.raw.mapstats
mv ${FINAL_BAM_FILE_MAPSTATS} "${OFPREFIX}".QC.final.mapstats

