from utils import job_runner

def remove_badcigar(raw_sam, badcigar_file, cpu_num, job):
    """Returns filtered BAM file name
       Finds reads with bad CIGAR strings, filters those reads,
       and delete raw SAM file and the bad CIGAR list.

    """
    badcigar_file = raw_sam + ".badcigar"
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

    raw_bam = raw_sam.split(".sam")[0] + ".raw.srt"

    with open(badcigar_file) as badcigar_fp:
        peek = [x for i,x in enumerate(fp) if i < 10 and x.rstrip()]

    badcigar_run_cmd = \
        "{} samtools view -u {} | samtools sort -@{} - {}"
    if len(peek):
        badcigar_flt_cmd = \
            "zcat {} | grep -vF -f {} |".format(raw_sam, badcigar_file)
        badcigar_run_cmd.format(badcigar_flt_cmd, "-", cpu_num, raw_bam)
    else:
        badcigar_run_cmd.format(" ", raw_sam, cpu_num, raw_bam)

    job.append([[ badcigar_run_cmd ]])
    job.append([[ "rm {} {}".format(badcigar_file, raw_sam) ]])
    job.run()

    return raw_bam + ".bam"


def remove_aberrant_reads(raw_bam, cpu_num, job):
    """Returns filtered BAM file name
       Filters unmapped, mate unmapped, not primary alignment,
       low quality reads, orphan reads (pair was removed), and
       read pairs mapping to different chromosomes.

    """

    flt_bam = raw_bam.split(".raw.")[0] + ".flt.srt.bam"

    tmp_flt_bam = flt_bam[:-4] + ".tmp.bam"
    tmp_clean_flt_bam = flt_bam[:-4] + ".tmp_clean.bam"

    tmp_flt_bam = "{}.tmp.nmsrt.bam".format(flt_bam[:-4])
    tmp_flt_bam_clean = "{}.clean.bam".format(tmp_flt_bam[:-4])

    filter_unmapped_cmd = \
        "samtools view -F1804 -f2 -q20 -u {} | \
         samtools sort -@{} -n - {}".format(raw_bam, cpu_num,
                                            tmp_flt_bam[:-4])
    fix_mate_cmd = \
        "samtools fixmate -r {} {}".format(tmp_flt_bam,
                                           tmp_flt_bam_clean)
    filter_orphan_cmd = \
        "samtools view -F1804 -f2 -u {} | \
         samtools -@{} sort - {}".format(tmp_flt_bam_clean, cpu_num,
                                         flt_bam[:-4])
    job.append([[ filter_unmapped_cmd ]])
    job.append([[ fix_mate_cmd ]])
    job.append([[ filter_orphan_cmd ]])
    job.append([[ "rm {} {}".format(tmp_flt_bam, tmp_flt_bam_clean) ]])
    job.run()

    return flt_bam

def remove_dups_sponges(picard_jar, flt_bam, cpu_num, sponge, job)
    """Returns

    """
    final_bam = flt_bam[:-4] + ".nodup.bam"

    tmp_bam = flt_bam[:-4] + ".dupmark.bam"
    dup_qc = flt_bam[:-4] + ".dup.qc"
    picard_cmd = \
        """java -Xmx4G -jar {} MarkDuplicates INPUT={} \
            OUTPUT={} METRICS_FILE={} ASSUME_SORTED=true \
            VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false
        """.format(picard_jar, flt_bam, tmp_bam, dup_qc)
    job.append([[ picard_cmd ]])

    #mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}
    #rm ${RAW_BAM_FILE}

    if sponge:
    samtools view -F 1804 -f 2 -h ${FILT_BAM_FILE} | \
    grep -v -F -f ${EXCLUDE} - | samtools view -Su - > ${FINAL_BAM_FILE}
    else:
    samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}

#### Remove duplicates and sponges, and create final name sorted BAM
FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai"
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc"


def filter(raw_sam, cpu_num):
    """Genome mapping function for ChIP-seq reads.

    """
    job = JobRunner()

    ## Filter reads with bad CIGAR string
    raw_bam = remove_badcigar(raw_sam, badcigar_file, cpu_num, job)

    ## TO-DO
    ## QC: samtools flagstat {prefix.raw.srt.bam} > {prefix.raw.srt.flagstat.qc}

    ## Filter other aberrant reads
    flt_bam = remove_aberrant_reads(raw_bam, cpu_num, job)

    ## Filter PCR duplicates
    remove_dups_sponges(flt_bam, cpu_num, job)

FINAL_NMSRT_BAM_PREFIX="${OFPREFIX}.filt.nmsrt.nodup"
FINAL_NMSRT_BAM_FILE="${FINAL_NMSRT_BAM_PREFIX}.bam"
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

