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
    def __init__(self, raw_sam, output_dir):
        self.output_dir = output_dir
        self.prefix = os.path.basename(raw_sam).split(".raw.")[0]

        if not os.path.basename(raw_sam):
            raise Exception("raw SAM file does not exist")
        self.raw_sam = raw_sam

        self.raw_bam = "%s/%s.raw.bam" %(output_dir, self.prefix)
        self.flt_bam = "%s/%s.flt.bam" %(output_dir, self.prefix)
        self.final_bam = "%s/%s.final.bam" %(output_dir, self.prefix)
        self.final_bai = "%s/%s.final.bai" %(output_dir, self.prefix)
        self.bed = "%s/%s.bed.gz" %(output_dir, self.prefix)
        self.bedpe = "%s/%s.bedpe.gz" %(output_dir, self.prefix)


    def remove_badcigar(self, threads=1):
        job = jr.JobRunner()

        # Find reads with bad CIGAR strings
        badcigar = "%s/%s.badcigar" %(self.output_dir, self.prefix)
        command = """zcat %s | \
                     awk 'BEGIN{FS="\t"; OFS="\t"}
                          !/^@/ && $6!="*" {
                          cigar=$6; gsub("[0-9]+D","",cigar);
                          n=split(cigar, vals, "[A-Z]");
                          s=0; for(i=1;i<=n;i++) s=s+vals[i];
                          seqlen=length($10);
                          if(s!=seqlen) print $1"\t";}' | \
                     sort | uniq > %s""" %(self.raw_sam, badcigar)
        job.append([[ command ]])
        job.run()

        with open(badcigar) as fp:
            peek = [x for i,x in enumerate(fp) if i<10 and x.rstrip()]

        # Remove the reads with bad CIGAR strings
        if len(peek):
            command = "zcat %s | grep -vF -f %s | " %(self.raw_sam, badcigar)\
                       + "samtools view -u - | " \
                       + "samtools sort -@ %s -o %s -" %(threads, self.raw_bam)
        else:
            command = "samtools view -u %s | " %(self.raw_sam)\
                       + "samtools sort -@ %s -o %s -" %(threads, self.raw_bam)
        job.append([[ command ]])
        job.run()

        job.append([[ "rm %s" %(badcigar) ]])
        job.run()


    def remove_aberrant_reads(self, threads=1):
        if not os.path.exists(self.raw_bam):
            raise Exception("raw BAM file does not exist")

        job = jr.JobRunner()

        tmp_bam = "%s/%s.tmp.bam" %(self.output_dir, self.prefix)
        tmp_clean_bam = "%s/%s.tmp.clean.bam" %(self.output_dir, self.prefix)

        # Filter unmapped and low quality (MAPQ<20) reads
        command = "samtools view -F1804 -f2 -q20 "\
                   + "-u %s | " %(self.raw_bam)\
                   + "samtools sort -@ %s -o %s -n -" %(threads, tmp_bam)
        job.append([[ command ]])
        job.run()

        # Filter orphan reads (pair was removed), and read pairs
        # mapping to different chromosomes
        command = "samtools fixmate -r %s %s" %(tmp_bam, tmp_clean_bam)
        job.append([[ command ]])
        job.run()
        command = "samtools view -F1804 -f2 -u %s | " %(tmp_clean_bam)\
                   + "samtools sort -@ %s -o %s -" %(threads, self.flt_bam)
        job.append([[ command ]])
        job.run()

        job.append([[ "rm %s %s" %(tmp_bam, tmp_clean_bam) ]])
        job.run()


    def remove_artifacts(self, picard_jar, sponge=None):
        if not os.path.exists(self.flt_bam):
            raise Exception("filtered BAM file does not exist!")

        job = jr.JobRunner()
        tmp_bam = "%s/%s.tmp.bam" %(self.output_dir, self.prefix)

        qc_dups = "%s/%s.pcrDups.QC.txt" %(self.output_dir, self.prefix)

        # Mark PCR duplicates and generate QC file with Picard
        command = "java -Xmx4G -jar %s MarkDuplicates " %(picard_jar)\
                   + "INPUT=%s OUTPUT=%s " %(self.flt_bam, tmp_bam)\
                   + "METRICS_FILE=%s ASSUME_SORTED=true " %(qc_dups)\
                   + "VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false"
        job.append([[ command ]])
        job.run()

        job.append([[ "mv %s %s" %(tmp_bam, self.flt_bam) ]])
        job.run()

        if sponge:
            command = "samtools view -F1804 -f2 -h %s | " %(self.flt_bam)\
                       + "grep -vF -f %s - | " %(sponge)\
                       + "samtools view -o %s -Su -" %(self.final_bam)
            job.append([[ command ]])
        else:
            command = "samtools view -F1804 -f2 -b "\
                       + "-o %s %s" %(self.final_bam, self.flt_bam)
            job.append([[ command ]])
        job.run()

        # Make BAM index file
        command = "samtools index %s %s" %(self.final_bam, self.final_bai)
        job.append([[ command ]])
        job.run()


    def generate_bed(self, threads=1):
        if not os.path.exists(self.final_bam):
            raise Exception("final BAM file does not exist")

        job = jr.JobRunner()

        tmp_bam = "%s/%s.tmp.bam" %(self.output_dir, self.prefix)

        command = "samtools sort -n -o %s %s" %(tmp_bam, self.final_bam)
        job.append([[ command ]])
        job.run()

        command = "bamToBed -bedpe -mate1 -i %s | " %(tmp_bam)\
                   + "gzip -c > %s" %(self.bedpe)
        job.append([[ command ]])
        job.run()

        command = """zcat %s | \
                     awk 'BEGIN{OFS="\t"; FS="\t"}
                          { chrom=$1; beg=$2; end=$6;
                            if($2>$5){beg=$5} if($3>$6){end=$3}
                            print chrom,beg,end
                          }' - > %s""" %(self.bedpe, self.bed)
        job.append([[ command ]])
        job.run()

        job.append([[ "rm %s" %(tmp_bam) ]])
        job.run()


    def clean(self):
        # Clean up intermediate files
        job = jr.JobRunner()
        job.append([[ "rm %s" %(self.raw_sam) ]])
        job.append([[ "rm %s" %(self.raw_bam) ]])
        job.append([[ "rm %s" %(self.flt_bam) ]])
        job.run()

        self.raw_sam = None
        self.raw_bam = None
        self.flt_bam = None

