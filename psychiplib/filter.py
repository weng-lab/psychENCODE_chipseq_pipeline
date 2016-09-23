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
    final_bam = flt_bam[:-4] + ".nodup.bai"

    tmp_bam = flt_bam[:-4] + ".dupmark.bam"
    dup_qc = flt_bam[:-4] + ".dup.qc"
    picard_cmd = \
        """java -Xmx4G -jar {} MarkDuplicates INPUT={} \
           OUTPUT={} METRICS_FILE={} ASSUME_SORTED=true \
           VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false
        """.format(picard_jar, flt_bam, tmp_bam, dup_qc)
    job.append([[ picard_cmd ]])
    job.append([[ "mv {} {}".format(tmp_bam, flt_bam) ]])

    #rm ${RAW_BAM_FILE} ## TODO

    if sponge:
        filter_dups_sponge_cmd = \
            """samtools view -F 1804 -f 2 -h {} | \
               grep -vF -f{} - | \
               samtools view -Su - > {}
            """.format(flt_bam, sponge, final_bam)
        job.append([[ filter_dups_sponge_cmd ]])
    else:
        filter_dups_cmd = \
            "samtools view -F1804 -f2 -b {} > {}".format(flt_bam,
                                                         final_bam)
        job.append([[ filter_dups_cmd ]])

    make_index_cmd = \
        "samtools index {} {}".foramt(final_bam, final_index)
    job.append([[ make_index_cmd ]])

    ## TODO
    ## FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc"
    ## QC: samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}"

    job.run()

    return final_bam


def generate_bedpe(final_bam, cpu_num, job):
    """Returns BEDPE file name
       Makes BEDPE file with the final cleaned BAM file.
    
    """
    bedpe3 = final_bam[:-4] + ".bedpe3.gz"
    bedpe7 = final_bam[:-4] + ".bedpe7.gz"
    nmsrt_bam = final_bam[:-4] + ".nmsrt.bam"

    name_sort_cmd = "samtools sort -n {} {}".format(final_bam,
                                                    nmsrt_bam[:-4])
    job.append([[ name_sort_cmd ]])

    bedtools_cmd = \
        """bedtools bamtobed -bedpe -mate1 -i {} | \
         gzip -c > {}""".format(nmsrt_bam, bedpe7)
    job.append([[ bedtools_cmd ]])

    bedpe3_cmd = \
        """zcat {} | \
           awk 'BEGIN{OFS="\t"; FS="\t"}
               {chrom=$1; beg=$2; end=$6;
               if($2>$5){beg=$5} if($3>$6){end=$3}
               print chrom,beg,end}' - > {}
        """.format(bedpe7, bedpe3)
    job.append([[ bedpe3_cmd ]])
    job.append([[ "rm {}".format(nmsrt_bam) ]])

    job.run()

    return (bedpe3, bedpe7)


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
    final_bam = remove_dups_sponges(flt_bam, cpu_num, job)

    ## Generate BEDPE file
    bedpe3, bedpe7 = generate_bedpe(final_bam, cpu_num, job)
