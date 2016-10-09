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

import utils.job_runner as jr
import os,sys,json,shutil
import psychiplib.mapping as mp
import psychiplib.post_filter as pf
import psychiplib.qc as qc
import psychiplib.pseudoreps as psr
import psychiplib.peak_calling as pc
import psychiplib.overlap as ovrlp
import traceback



def map_inputs(inputs,idx):
    ### ALIGN INPUT LIBRARIES
    run_map_input = mp.BwaMapper(inputs['inputs'][idx][0], # fastq read 1
                            inputs['inputs'][idx][1], # fastq read 2
                            inputs['prefix_inputs'][idx], #prefix
                            inputs['index'], # bwa index
                            os.path.join(os.path.join(inputs['outdir'],\
                                         inputs['prefix_inputs'][idx]),\
                                         "bwa_out"))
                                         # output directory

    run_map_input.run(inputs['cpus']) # CPU for BWA
                     # Output: R1.raw.sam.gz
                     # Format: <prefix>.raw.sam.gz


    ### Post alignment filtering
    run_filter_input = pf.PostFilter("%s/%s/bwa_out/%s.raw.sam.gz" %(inputs['outdir'], inputs['prefix_inputs'][idx],inputs['prefix_inputs'][idx]), \
                               "%s/%s/bwa_out" %(inputs['outdir'], inputs['prefix_inputs'][idx])) # output directory

    ## Step 1: remove reads with bad CIGAR strings
    run_filter_input.remove_badcigar(inputs['cpus']) # Threads for SAMtools
                                    # Output: R1.raw.bam
                                    # Format: <prefix>.raw.bam

    ## Step 2: remove unmapped, low quality, orphan, and supplemental reads
    run_filter_input.remove_aberrant_reads(inputs['cpus']) # Threads for SAMtools
                                          # Output: R1.flt.bam
                                          # Format: <prefix>.flt.bam

    ## Step 3: remove PCR duplicates and sponge mapping reads
    ##         (currently, "java -Xmx4G". Please change it if needed)
    run_filter_input.remove_artifacts(inputs['picard'],  # path to picard.jar
                                 inputs['threads'],
                                 inputs['sponges']) # sponge name list
                                                        # Output: R1.final.[bam|bai], R1.pcrDups.QC.txt

    ## Step 4: generate BED and BEDPE files
    ##         (those files are *unsorted*)
    run_filter_input.generate_bed(inputs['cpus']) # Threads for SAMtools
                                 # Output: R1.bed.gz (for pseudorep)
                                 # Output: R1.bedpe.gz (bedtool raw output)


    if not os.path.isdir("%s/%s/macs_input/" %(inputs['outdir'], inputs['prefix_inputs'][idx])):
        os.mkdir("%s/%s/macs_input/" %(inputs['outdir'], inputs['prefix_inputs'][idx]))
    if not os.path.isfile("%s/%s/macs_input/%s.bed.gz" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx])):
        shutil.move("%s/%s/bwa_out/%s.bed.gz" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx]),\
               "%s/%s/macs_input/%s.bed.gz" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx]))
    
    ### QC part
    test_qc = qc.QC("%s/%s/bwa_out/" %(inputs['outdir'], inputs['prefix_inputs'][idx])) # output direcotry

    ## Step 1: Get map stats
    test_qc.get_mapstats("%s/%s/bwa_out/%s.raw.bam" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx]),
                         "%s/%s/bwa_out/%s.final.bam" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx])) # Output: R1.final.flagstat.QC.txt

    ## Step 2: Compute library complexity
    test_qc.get_library_complexity("%s/%s/bwa_out/%s.flt.bam" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx]), inputs['cpus']) # Threads for SAMtools
                                                              # Output: R1.libComplexity.QC.txt
    ''' 
    ## Step 3: *!NOT TESTED!* Cross correlation coefficient
    test_qc.cross_correlation("%s/%s/bwa_out/%s.bedpe.gz" %(inputs['outdir'], inputs['prefix_controls'][idx], inputs['prefix_controls'][idx]), # BEDPE file
                              "%s/common/run_spp_nodups.R" %(os.path.dirname(os.path.abspath(__file__))), # phantompeakqualtools
                              inputs['cpus'])                                      # Number of CPU for phantompeakqualtools
    '''
    job=jr.JobRunner()
    insert = 150
    job.append([["java -jar %s CollectInsertSizeMetrics I=%s/%s/bwa_out/%s.final.bam\
                  O=%s/%s/bwa_out/InsertSizeMetrics.txt AS=F H=%s/%s/bwa_out/InsertSizeMetrics.histogram" \
                  %(inputs['picard'],inputs['outdir'],inputs['prefix_inputs'][idx],inputs['prefix_inputs'][idx],\
                    inputs['outdir'],inputs['prefix_inputs'][idx],\
                    inputs['outdir'],inputs['prefix_inputs'][idx])]])
    job.run()

    job.append([["grep -A2 \"## METRICS CLASS\" %s/%s/bwa_out/InsertSizeMetrics.txt | tail -1 | cut -f1" %(inputs['outdir'],inputs['prefix_inputs'][idx])]])
    insert=job.run()[0].rstrip()

    run_filter_input.clean()


    ### PSEUDOREP
    run_psr = psr.PseudorepsGenerator("%s/%s/macs_input/%s.bed.gz" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx]),\
                                      inputs['prefix_inputs'][idx],
                                      "%s/%s/macs_input/" %(inputs['outdir'], inputs['prefix_inputs'][idx]),
                                      inputs['threads'],
                                      inputs['cpus'])

    run_psr.run()


def map_controls(inputs, idx):
    ### ALIGN CONTROL LIBRARIES
    run_map_control = mp.BwaMapper(inputs['controls'][idx][0], # fastq read 1
                            inputs['controls'][idx][1], # fastq read 2
                            inputs['prefix_controls'][idx], #prefix
                            inputs['index'], # bwa index
                            os.path.join(os.path.join(inputs['outdir'],\
                                         inputs['prefix_controls'][idx]),\
                                         "bwa_out"))
                                         # output directory

    run_map_control.run(inputs['cpus']) # CPU for BWA
                     # Output: R1.raw.sam.gz
                     # Format: <prefix>.raw.sam.gz


    ### Post alignment filtering
    run_filter_control = pf.PostFilter("%s/%s/bwa_out/%s.raw.sam.gz" %(inputs['outdir'], inputs['prefix_controls'][idx],inputs['prefix_controls'][idx]),
                               "%s/%s/bwa_out" %(inputs['outdir'], inputs['prefix_controls'][idx])) # output directory

    ## Step 1: remove reads with bad CIGAR strings
    run_filter_control.remove_badcigar(inputs['cpus']) # Threads for SAMtools
                                    # Output: R1.raw.bam
                                    # Format: <prefix>.raw.bam

    ## Step 2: remove unmapped, low quality, orphan, and supplemental reads
    run_filter_control.remove_aberrant_reads(inputs['cpus']) # Threads for SAMtools
                                          # Output: R1.flt.bam
                                          # Format: <prefix>.flt.bam

    ## Step 3: remove PCR duplicates and sponge mapping reads
    ##         (currently, "java -Xmx4G". Please change it if needed)
    run_filter_control.remove_artifacts(inputs['picard'],  # path to picard.jar
                                 inputs['threads'],
                                 inputs['sponges']) # sponge name list
                                                        # Output: R1.final.[bam|bai], R1.pcrDups.QC.txt

    ## Step 4: generate BED and BEDPE files
    ##         (those files are *unsorted*)
    run_filter_control.generate_bed(inputs['cpus']) # Threads for SAMtools
                                 # Output: R1.bed.gz (for pseudorep)
                                 # Output: R1.bedpe.gz (bedtool raw output)


    if not os.path.isdir("%s/%s/macs_input/" %(inputs['outdir'], inputs['prefix_controls'][idx])):
        os.mkdir("%s/%s/macs_input/" %(inputs['outdir'], inputs['prefix_controls'][idx]))
    if not os.path.isfile("%s/%s/macs_input/%s.bed.gz" %(inputs['outdir'], inputs['prefix_controls'][idx], inputs['prefix_controls'][idx])):
        shutil.move("%s/%s/bwa_out/%s.bed.gz" %(inputs['outdir'], inputs['prefix_controls'][idx], inputs['prefix_controls'][idx]),\
               "%s/%s/macs_input/%s.bed.gz" %(inputs['outdir'], inputs['prefix_controls'][idx], inputs['prefix_controls'][idx]))
    run_filter_control.clean()



def pool_libraries(inputs):
    if not os.path.isdir("%s/oracle_peaks" %(inputs['outdir'])):
        os.mkdir("%s/oracle_peaks" %(inputs['outdir']))
    job=jr.JobRunner()
    cmd_inputs = "zcat "
    for prefix in inputs['prefix_inputs']:
        cmd_inputs = cmd_inputs + "%s/%s/macs_input/%s.bed.gz " %(inputs['outdir'], prefix, prefix)
    cmd_inputs = cmd_inputs + "| sort --parallel=%s -k1,1 -k2,2n | gzip -c > %s/oracle_peaks/pooled_inputs.bed.gz" %(inputs['cpus'], inputs['outdir'])
    print(cmd_inputs)
    cmd_controls = "zcat "
    for prefix in inputs['prefix_controls']:
        cmd_controls = cmd_controls + "%s/%s/macs_input/%s.bed.gz " %(inputs['outdir'], prefix, prefix) 
    cmd_controls = cmd_controls + "| sort --parallel=%s -k1,1 -k2,2n | gzip -c > %s/oracle_peaks/pooled_controls.bed.gz" %(inputs['cpus'], inputs['outdir'])    
    print(cmd_controls)
    job.append([[ cmd_controls ]])
    job.append([[ cmd_inputs  ]])
    job.run()


def peak(inputs, idx):
    print("SONOQUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII")
    job=jr.JobRunner()
    job.append([["grep -A2 \"## METRICS CLASS\" %s/%s/bwa_out/InsertSizeMetrics.txt | tail -1 | cut -f1" %(inputs['outdir'],inputs['prefix_inputs'][idx])]])
    insert=job.run()[0].rstrip()
    ### PEAK
    run_peak=pc.Macs2PeakCaller("%s/%s/macs_input/%s.bed.gz" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx]),
                         "%s/oracle_peaks/pooled_controls.bed.gz" %(inputs['outdir']),
                         insert,
                         "BEDPE",
                         "%s/%s/macs_output/%s" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx]),
                         inputs['genome'],"hs", inputs['macs2'],inputs['common'])
    run_peak.run()

    run_peak=pc.Macs2PeakCaller("%s/%s/macs_input/psr_%s.00.bedpe.gz" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx]),
                         "%s/oracle_peaks/pooled_controls.bed.gz" %(inputs['outdir']),
                         insert,
                         "BEDPE",
                         "%s/%s/macs_output/psr_%s.00" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx]),
                         inputs['genome'],
                         "hs", inputs['macs2'],inputs['common'])
    run_peak.run()

    run_peak=pc.Macs2PeakCaller("%s/%s/macs_input/psr_%s.01.bedpe.gz" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx]),
                         "%s/oracle_peaks/pooled_controls.bed.gz" %(inputs['outdir']),
                         insert,
                         "BEDPE",
                         "%s/%s/macs_output/psr_%s.01" %(inputs['outdir'], inputs['prefix_inputs'][idx], inputs['prefix_inputs'][idx]),
                         inputs['genome'],
                         "hs",inputs['macs2'],inputs['common'])
    run_peak.run()

    ### OVERLAP
    run_overlap=ovrlp.Overlap(inputs['prefix_inputs'][idx],inputs['outdir'])
    run_overlap.run()

def peak_oracle(inputs):
    insert = 0
    job=jr.JobRunner()
    for prefix in inputs['prefix_inputs']:
        job.append([["grep -A2 \"## METRICS CLASS\" %s/%s/bwa_out/InsertSizeMetrics.txt | tail -1 | cut -f1" %(inputs['outdir'], prefix)]])
        insert = insert + int(job.run()[0].rstrip())
    insert = insert/len(inputs['prefix_inputs'])

    ### PEAK
    run_peak=pc.Macs2PeakCaller("%s/oracle_peaks/pooled_inputs.bed.gz" %(inputs['outdir']),
                         "%s/oracle_peaks/pooled_controls.bed.gz" %(inputs['outdir']),
                         insert,
                         "BEDPE",
                         "%s/oracle_peaks/macs_output" %(inputs['outdir']),
                         inputs['genome'],"hs", inputs['macs2'],inputs['common'])
    run_peak.run()




try:
  
    inputs={}
    with open(sys.argv[1], 'r') as inputs_json:
        inputs=json.load(inputs_json)
    
    if sys.argv[2] == "map":
        if sys.argv[4] == "controls":
            map_controls(inputs,int(sys.argv[3]))
        else:
            map_inputs(inputs,int(sys.argv[3]))
    elif sys.argv[2] == "pooling":
        pool_libraries(inputs)
    elif sys.argv[2] == "peak":
        if int(sys.argv[3]) < len(inputs['inputs']):
            peak(inputs,int(sys.argv[3]))
        else:
            peak_oracle(inputs)

except Exception, err:
    try:
        exc_info = sys.exc_info()
        try:
            raise TypeError("Error in one of the module")
        except:
            pass
    finally:
        # Display the *original* exception
        traceback.print_exception(*exc_info)
        del exc_info
