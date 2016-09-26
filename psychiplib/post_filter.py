########################################################################
# psychip - ChIP-seq pipeline for psychENCODE project
# Copyright (C) 2016  Weng's Lab @ UMASS medical school - Worcester (MA)
#
# Developers:
# Eugenio Mattei : eugenio.mattei@umassmed.edu
# Junko Tsuji : junko.tsuji@umassmed.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################


#!/usr/bin/env python

import os

import utils.job_runner as jr


class PostFilter():
    def __init__(self, output_dir):
        self.output_dir = output_dir

    def remove_badcigar(self, raw_sam, threads=1):
        job = jr.JobRunnter()

        prefix = os.path.basename(raw_sam).split(".raw.")[0]

        # Find reads with bad CIGAR strings
        badcigar = "%s/%s.badcigar" %(self.output_dir, prefix)
        command = """zcat %s | \
                     awk 'BEGIN{FS="\t"; OFS="\t"}
                          !/^@/ && $6!="*" {
                          cigar=$6; gsub("[0-9]+D","",cigar);
                          n=split(cigar, vals, "[A-Z]");
                          s=0; for(i=1;i<=n;i++) s=s+vals[i];
                          seqlen=length($10);
                          if(s!=seqlen) print $1"\t";}' | \
                     sort | uniq > %s""" %(raw_sam, badcigar)
        job.append([[ command ]])
        job.run()

        raw_bam = "%s/%s.raw.bam" %(self.output_dir, prefix)

        with open(badcigar) as fp:
            peek = [x for i,x in enumerate(fp) if i<10 and x.rstrip()]

        # Remove the reads with bad CIGAR strings
        if len(peek):
            command = "zcat %s | grep -vF -f %s | " %(raw_sam, badcigar)
                       + "samtools view -u - | "
                       + "samtools sort -@ %s -o %s -" %(threads, raw_bam)
        else:
            command = "samtools view -u %s | " %(raw_sam)
                       + "samtools sort -@ %s -o %s -" %(threads, raw_bam)
        job.append([[ command ]])
        job.append([[ "rm %s" %(badcigar) ]])
        job.run()

        return raw_bam


    def remove_aberrant_reads(self, raw_bam, threads=1):
        job = jr.JobRunner()

        prefix = os.path.basename(raw_bam).split(".raw.")
        flt_bam = "%s/%s.flt.bam" %(self.output_dir, prefix)
        tmp_bam = "%s/%s.tmp.bam" %(self.output_dir, prefix)
        tmp_clean_bam = "%s/%s.tmp.clean.bam" %(self.output_dir, prefix)

        # Filter unmapped and low quality (MAPQ<20) reads
        command = "samtools view -F1804 -f2 -q20 "
                   + "-u %s | " %(raw_bam)
                   + "samtools sort -@ %s -o %s -n -" %(threads, tmp_bam)
        job.append([[ command ]])

        # Filter orphan reads (pair was removed), and read pairs
        # mapping to different chromosomes
        command = "samtools fixmate -r %s %s" %(tmp_bam, tmp_clean_bam)
        job.append([[ command ]])
        command = "samtools view -F1804 -f2 -u %s | " %(tmp_clean_bam)
                   + "samtools sort -@ %s -o %s -" %(threads, flt_bam)
        job.append([[ command ]])

        job.append([[ "rm %s %s" %(tmp_bam, tmp_clean_bam) ]])
        job.run()

        return flt_bam


    def remove_artifacts(self, picard_jar, flt_bam, sponge=None)
        job = jr.JobRunner()

        prefix = os.path.basename(flt_bam).split(".flt.")
        final_bam = "%s/%s.final.bam" %(self.output_dir, prefix)
        final_bai = "%s/%s.final.bai" %(self.output_dir, prefix)
        tmp_bam = "%s/%s.tmp.bam" % (self.output_dir, prefix)

        qc_dups = "%s/%s_PCRdups.QC.txt"

        # Mark PCR duplicates and generate QC file with Picard
        command = "java -Xmx4G -jar %s MarkDuplicates " %(picard_jar)
                   + "INPUT=%s OUTPUT=%s " %(flt_bam, tmp_bam)
                   + "METRICS_FILE=%s ASSUME_SORTED=true " %(qc_dups)
                   + "VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false"
        job.append([[ command ]])
        job.append([[ "mv %s %s" %(tmp_bam, flt_bam) ]])

        if sponge:
            command = "samtools view -F1804 -f2 -h %s | " %(flt_bam)
                       + "grep -vF -f %s - | " %(sponge)
                       + "samtools view -o %s -Su -" %(final_bam)
            job.append([[ command ]])
        else:
            command = "samtools view -F1804 -f2 -b "
                       + "-o %s %s" %(final_bam, flt_bam)
            job.append([[ command ]])

        # Make BAM index file
        command = "samtools index %s %s" %(final_bam, final_bai)
        job.append([[ command ]])
        job.run()
   
        return final_bam


    def generate_bed(self, final_bam, threads=1):
        job = jr.JobRunner()

        prefix = os.path.basename(final_bam).split(".final.")[0]

        bed = "%s/%s.bed.gz" %(self.output_dir, prefix)
        bedpe = "%s/%s.bedpe.gz" %(self.output_dir, prefix)
        tmp_bam = "%s/%s.tmp.bam" %(self.output_dir, prefix)

        command = "samtools sort -n -o %s %s" %(tmp_bam, final_bam)
        job.append([[ command ]])

        command = "bamToBed -bedpe -mate1 -i %s | " %(tmp_bam)
                   + "gzip -c > %s" %(bedpe)
        job.append([[ command ]])

        command = """zcat %s | \
                     awk 'BEGIN{OFS="\t"; FS="\t"}
                          { chrom=$1; beg=$2; end=$6;
                            if($2>$5){beg=$5} if($3>$6){end=$3}
                            print chrom,beg,end
                          }' - > %s""" %(bedpe, bed)
        job.append([[ command ]])
        job.append([[ "rm %s" %(tmp_bam) ]])
        job.run()

        return (bed, bedpe)
