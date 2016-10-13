########################################################################
# psychip - ChIP-seq pipeline for psychENCODE project
# Copyright (C) 2016  Weng's Lab @ UMASS medical school - Worcester (MA)
#
# Developers:
# Eugenio Mattei : eugenio.mattei@umassmed.edu
# Junko Tsuji : junko.tsuji@umassmed.edu
# 
# modified from https://github.com/ENCODE-DCC/chip-seq-pipeline
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

import os, subprocess, shlex, time
import utils.job_runner as jr

class Xcor():
    def __init__(self, input_bam, spp_file, cpus):
        self.input_bam_filename = input_bam
        self.input_bam_basename = self.input_bam_filename.rstrip('.bam')
        self.paired_end = False # hardcoded. to be changed in the future
        self.cpus = cpus
        self.spp = spp_file

    def process(self):
        job = jr.JobRunner()
        
	self.intermediate_TA_filename = self.input_bam_basename + ".tagAlign"
	if self.paired_end:
	  end_infix = 'PE2SE'
	else:
	  end_infix = 'SE'
	self.final_TA_filename = self.input_bam_basename + '.' + end_infix + '.tagAlign.gz'

	# ===================
	# Create tagAlign file
	# ===================

        job.append([[
        "bamToBed -i %s | " %(self.input_bam_filename),
        r"""awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'""",
        "| tee %s | " %(self.intermediate_TA_filename),
        "gzip -c > %s" %(self.final_TA_filename)]])
        job.run()
        
	# ================
	# Create BEDPE file
	# ================
        if self.paired_end:
            self.final_BEDPE_filename = self.input_bam_basename + ".bedpe.gz"
            #need namesorted bam to make BEDPE
            self.final_nmsrt_bam_prefix = self.input_bam_basename + ".nmsrt"
            self.final_nmsrt_bam_filename = self.final_nmsrt_bam_prefix + ".bam"
            job.append([["samtools sort -@ %s -n -o %s %s" \
			 %(self.cpus, self.final_nmsrt_bam_filename, self.input_bam_filename)]])
            job.run()
            job.append([[
            "bedtools bamtobed -bedpe -mate1 -i %s | " %(self.final_nmsrt_bam_filename),
            "gzip -c > %s" %(self.final_BEDPE_filename)]])
            job.run()
        # =================================
        # Subsample tagAlign file
        # ================================
        NREADS=15000000
        if self.paired_end:
            end_infix = 'MATE1'
        else:
            end_infix = 'SE'
        self.subsampled_TA_filename = self.input_bam_basename + ".filt.nodup.sample.%d.%s.tagAlign.gz" %(NREADS/1000000, end_infix)
        steps = 'grep -v "chrM" %s | ' %(self.intermediate_TA_filename) + \
            'shuf -n %d | ' %(NREADS)
        if self.paired_end:
            steps += r"""awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | """
        steps += 'gzip -c > %s' %(self.subsampled_TA_filename)
        print(steps)
        job.append([[ steps ]])
        job.run()
        # Calculate Cross-correlation QC scores
        self.CC_scores_filename = self.subsampled_TA_filename + ".cc.qc"
        self.CC_plot_filename = self.subsampled_TA_filename + ".cc.plot.pdf"

        # CC_SCORE FILE format
        # Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag
        print("Rscript %s -c=%s -p=%d -filtchr=chrM -savp=%s -out=%s" %(self.spp, self.subsampled_TA_filename, self.cpus, self.CC_plot_filename, self.CC_scores_filename))
        job.append([[
            "Rscript %s -c=%s -p=%d -filtchr=chrM -savp=%s -out=%s > /dev/null 2>&1" \
                %(self.spp, self.subsampled_TA_filename, self.cpus, self.CC_plot_filename, self.CC_scores_filename)]])
        job.run()
        
        job.append([[r"""sed -r  's/,[^\t]+//g' %s > %s """ %(self.CC_scores_filename, "temp")]])
        job.run()
        
        job.append([[
            "mv temp %s" %(self.CC_scores_filename)]])
        job.run()


