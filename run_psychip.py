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

from __future__ import print_function

import utils.job_runner
import os,sys
import psychiplib.mapping as mp
import psychiplib.post_filter as pf
import psychiplib.qc as qc
import 
import
import


inputs = sys.argv[1]
idx = int(sys.argv[2])


test_map = mp.BwaMapper(inputs['inputs'][idx][0], # fastq read 1
                        inputs['inputs'][idx][1], # fastq read 2
                        inputs['prefix_inputs'][idx], #prefix
                        inputs['sponges'], # bwa index
                        os.path.join(inputs['outdir'],\
                                     inputs['prefix_inputs'][idx]))
                                     # output directory

test_map.run(inputs['cpus']) # CPU for BWA
                 # Output: R1.raw.sam.gz
                 # Format: <prefix>.raw.sam.gz


### Post alignment filtering
test_filter = pf.PostFilter("%s/%s.raw.sam.gz" %(os.path.join(inputs['outdir'], inputs['prefix_inputs'][idx]))) #raw sam file
                            inputs['outdir']) # output directory

## Step 1: remove reads with bad CIGAR strings
test_filter.remove_badcigar(inputs['cpus']) # Threads for SAMtools
                                # Output: R1.raw.bam
                                # Format: <prefix>.raw.bam

## Step 2: remove unmapped, low quality, orphan, and supplemental reads
test_filter.remove_aberrant_reads(inputs['cpus']) # Threads for SAMtools
                                      # Output: R1.flt.bam
                                      # Format: <prefix>.flt.bam

## Step 3: remove PCR duplicates and sponge mapping reads
##         (currently, "java -Xmx4G". Please change it if needed)
test_filter.remove_artifacts("path/to/picard.jar",  # path to picard.jar
                             inputs['sponges']) # sponge name list
                                                    # Output: R1.final.[bam|bai], R1.pcrDups.QC.txt

## Step 4: generate BED and BEDPE files
##         (those files are *unsorted*)
test_filter.generate_bed(inputs['cpus']) # Threads for SAMtools
                             # Output: R1.bed.gz (for pseudorep)
                             # Output: R1.bedpe.gz (bedtool raw output)


### QC part
#test_qc = qc.QC("test_out") # output direcotry

## Step 1: Get map stats
#test_qc.get_mapstats("test_out/R1.raw.bam",   # Output: R1.raw.flagstat.QC.txt
#                     "test_out/R1.final.bam") # Output: R1.final.flagstat.QC.txt

## Step 2: Compute library complexity
#test_qc.get_library_complexity("test_out/R1.flt.bam", 16) # Threads for SAMtools
                                                          # Output: R1.libComplexity.QC.txt

## Step 3: *!NOT TESTED!* Cross correlation coefficient
#test_qc.cross_correlation("test_out/R1.bedpe.gz",                  # BEDPE file
#                          "phantompeakqualtools/run_spp_nodups.R", # phantompeakqualtools
#                          32)                                      # Number of CPU for phantompeakqualtools

## Clean up intermediate files (if needed)
## It removes R1.raw.sam.gz, R1.raw.bam, R1.flt.bam
#test_filter.clean()


### PSEUDOREP

### PEAK

### OVERLAP

