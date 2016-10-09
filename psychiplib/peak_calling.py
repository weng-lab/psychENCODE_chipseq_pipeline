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

import os, sys, time,subprocess
import utils.common
import utils.job_runner as jr

class Macs2PeakCaller():
    def __init__(self, experiment, control, xcor_scores_input,
                     input_type, output_name, chrom_sizes, genomesize, macs2, common):
        # Initialize data object inputs on the platform

        self.experiment_name = experiment
        self.control_name = control
        self.xcor_scores_input_name = xcor_scores_input
        self.input_type = input_type
        self.output_name = output_name
        self.genomesize = genomesize
        self.narrowPeak_as_name = os.path.join(common,"narrowPeak.as")
        self.gappedPeak_as_name = os.path.join(common,"gappedPeak.as")
        self.broadPeak_as_name = os.path.join(common,"broadPeak.as")
        self.chrom_sizes_name = chrom_sizes
        self.macs2 = macs2
        self.common = common

    def run(self):
        job=jr.JobRunner()

        #Define the output filenames
        peaks_dirname = self.output_name
        if not os.path.exists(peaks_dirname):
            os.makedirs(peaks_dirname)
        prefix =self.experiment_name.split("/")[-1]
        if prefix.endswith('.bedpe.gz'):
            prefix = prefix[:-9]
        elif prefix.endswith('.bed.gz') :
            prefix = prefix[:-7]

        narrowPeak_fn    = "%s/%s.narrowPeak" %(peaks_dirname, prefix)
        gappedPeak_fn    = "%s/%s.gappedPeak" %(peaks_dirname, prefix)
        broadPeak_fn     = "%s/%s.broadPeak"  %(peaks_dirname, prefix)
        self.narrowPeak_gz_fn = narrowPeak_fn + ".gz"
        self.gappedPeak_gz_fn = gappedPeak_fn + ".gz"
        self.broadPeak_gz_fn  = broadPeak_fn  + ".gz"
        self.narrowPeak_bb_fn = "%s.bb" %(narrowPeak_fn)
        self.gappedPeak_bb_fn = "%s.bb" %(gappedPeak_fn)
        self.broadPeak_bb_fn  = "%s.bb" %(broadPeak_fn)
        self.fc_signal_fn = "%s/%s.fc_signal.bw" \
                             %(peaks_dirname, prefix)
        self.pvalue_signal_fn = "%s/%s.pvalue_signal.bw" \
                                 %(peaks_dirname, prefix)

        # Extract the fragment length estimate from \
        # column 3 of the cross-correlation scores file
#        if self.xcor_scores_input_name == "150":
        fraglen=self.xcor_scores_input_name
#        else:
#            with open(self.xcor_scores_input_name,'r') as fh:
#                firstline = fh.readline()
#                fraglen = firstline.split()[2] #third column
#                print "Fraglen %s" %(fraglen)
        #===========================================
        command = '%s callpeak '%(self.macs2) + \
                  '-t %s -c %s '\
                   %(self.experiment_name, self.control_name) + \
                  '-f %s -n %s/%s '\
                   %(self.input_type, peaks_dirname, prefix) + \
                  '-g %s -p 1e-2 --nomodel --shift 0 --extsize %s \
                  --keep-dup all -B --SPMR' %(self.genomesize, fraglen)
        job.append([["%s" %(command)]])
        job.run()
        #print command
        #returncode = common.block_on(command)
        #print "MACS2 exited with returncode %d" %(returncode)
        #assert returncode == 0, "MACS2 non-zero return"

        # Rescale Col5 scores to range 10-1000 to conform \
        # to narrowPeak.as format (score must be <1000)
        rescaled_narrowpeak_fn = utils.common.rescale_scores(
                                 '%s/%s_peaks.narrowPeak' \
                                  %(peaks_dirname, prefix), 
                                  scores_col=5)

        # Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank>
        command = "sort -k 8gr,8gr %s |\
                   awk 'BEGIN{OFS=\"\\t\"}{$4=\"Peak_\"NR ; \
                   print $0}' | tee %s | gzip -c > %s" \
                   %(rescaled_narrowpeak_fn, narrowPeak_fn,
                   self.narrowPeak_gz_fn)

        job.append([["%s" %(command)]])
        job.run()

        # remove additional files
        command ="rm -f %s/%s_peaks.xls %s/%s_peaks.bed %s_summits.bed"\
                  %(peaks_dirname, prefix,
                    peaks_dirname, prefix, prefix)
        job.append([["%s" %(command)]])
        job.run()
        #===========================================
        # Generate Broad and Gapped Peaks
        #============================================
        command = '%s callpeak ' %(self.macs2) + \
                  '-t %s -c %s ' \
                  %(self.experiment_name, self.control_name) + \
                  '-f %s -n %s/%s '%(self.input_type, 
                                    peaks_dirname, prefix) + \
                  '-g %s -p 1e-2 --broad --nomodel --shift 0 \
                  --extsize %s --keep-dup all' \
                  %(self.genomesize, fraglen)
        #print command
        job.append([["%s" %(command)]])
        job.run()
        # Rescale Col5 scores to range 10-1000 to conform to narrowPeak.as format (score must be <1000)
        rescaled_broadpeak_fn = utils.common.rescale_scores('%s/%s_peaks.broadPeak'
                                               %(peaks_dirname, prefix),
                                               scores_col=5)
        # Sort by Col8 (for broadPeak) or Col 14(for gappedPeak)  in descending order and replace long peak names in Column 4 with Peak_<peakRank>
        command = "sort -k 8gr,8gr %s | \
                   awk 'BEGIN{OFS=\"\\t\"}{$4=\"Peak_\"NR ; print $0}'|\
                   tee %s | gzip -c > %s"%(rescaled_broadpeak_fn, 
                                           broadPeak_fn,
                                           self.broadPeak_gz_fn)
        #print command

        job.append([["%s" %(command)]])
        job.run()

        # Rescale Col5 scores to range 10-1000 to conform to narrowPeak.as format (score must be <1000)
        rescaled_gappedpeak_fn = utils.common.rescale_scores('%s/%s_peaks.gappedPeak' 
                                       %(peaks_dirname, prefix),
                                       scores_col=5)

        command = "sort -k 14gr,14gr %s | \
                   awk 'BEGIN{OFS=\"\\t\"}{$4=\"Peak_\"NR ; print $0}'|\
                   tee %s | gzip -c > %s"%(rescaled_gappedpeak_fn,\
                                           gappedPeak_fn,
                                           self.gappedPeak_gz_fn)
        #print command

        job.append([["%s" %(command)]])
        job.run()

        # remove additional file
        job.append([["rm -f %s/%s_peaks.xls %s/%s_peaks.bed %s_summits.bed" %(peaks_dirname,prefix,peaks_dirname,prefix,prefix)]])
        job.run()

        #===========================================
        # For Fold enrichment signal tracks
        #============================================
        # This file is a tab delimited file with 2 columns Col1 (chromosome name), Col2 (chromosome size in bp).
        command = '%s bdgcmp ' %(self.macs2) + \
                  '-t %s/%s_treat_pileup.bdg ' %(peaks_dirname, prefix) + \
                  '-c %s/%s_control_lambda.bdg ' %(peaks_dirname, prefix) + \
                  '--outdir %s -o %s_FE.bdg ' %(peaks_dirname, prefix) + \
                  '-m FE'
        #print command

        job.append([["%s" %(command)]])
        job.run()
        # Remove coordinates outside chromosome sizes (stupid MACS2 bug)
        job.append([["bedtools slop -i %s/%s_FE.bdg -g %s -b 0 |\
                    %s/bedClip stdin %s %s/%s.fc.signal.bedgraph"
                    %(peaks_dirname, prefix, self.chrom_sizes_name, self.common, 
                      self.chrom_sizes_name, peaks_dirname, prefix)]])
        #print command
        job.run()

        #rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}_FE.bdg
        # Convert bedgraph to bigwig
        command = '%s/bedGraphToBigWig '%(self.common) + \
                  '%s/%s.fc.signal.bedgraph ' %(peaks_dirname, prefix)+\
                  '%s %s' %(self.chrom_sizes_name, self.fc_signal_fn)
        #print command
        job.append([["%s" %(command)]])
        job.run()

        #rm -f ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.bedgraph
        #===========================================
        # For -log10(p-value) signal tracks
        #===========================================
        # Compute sval = min(no. of reads in ChIP, no. of reads in control) / 1,000,000
        if (self.input_type=="BED" or self.input_type=="BEDPE") :
            job.append([['gzip -dc %s' %(self.experiment_name)+'|wc -l']])
            chipReads = int(job.run()[0])
            job.append([['gzip -dc %s' %(self.control_name)+'|wc -l']])
            controlReads = int(job.run()[0])
            sval=str(min(float(chipReads), float(controlReads))/1000000)
        else:
            job.append([["samtools idxstats %s | awk '{sum=sum+$3}END{print sum}'"%(self.experiment_name)]])
            chipReads = int(job.run()[0])
            job.append([["samtools idxstats %s | awk '{sum=sum+$3}END{print sum}'"%(self.control_name)]])
            controlReads = int(job.run()[0])
            sval=str(min(float(chipReads), float(controlReads))/1000000)
        #    print sval,chipReads,controlReads
        print "chipReads = %s, controlReads = %s, sval = %s" %(chipReads, controlReads, sval)

        job.append([['%s bdgcmp ' %(self.macs2) + \
                '-t %s/%s_treat_pileup.bdg ' %(peaks_dirname, prefix) + \
                '-c %s/%s_control_lambda.bdg ' %(peaks_dirname, prefix) + \
                '--outdir %s -o %s_ppois.bdg ' %(peaks_dirname, prefix) + \
                '-m ppois -S %s' %(sval)]])
        job.run()
        # Remove coordinates outside chromosome sizes (stupid MACS2 bug)
        job.append([['bedtools slop -i %s/%s_ppois.bdg -g %s -b 0' 
                   %(peaks_dirname, prefix, self.chrom_sizes_name)+ \
                '| %s/bedClip stdin %s %s/%s.pval.signal.bedgraph'
                 %(self.common, self.chrom_sizes_name,peaks_dirname, prefix)]])
        job.run()

        job.append([["rm -rf %s/%s_ppois.bdg" %(peaks_dirname,prefix)]])
        job.run()

        # Convert bedgraph to bigwig
        command = '%s/bedGraphToBigWig ' %(self.common) + \
        '%s/%s.pval.signal.bedgraph ' %(peaks_dirname, prefix) + \
        '%s %s' %(self.chrom_sizes_name, self.pvalue_signal_fn)

        job.append([["%s" %(command)]])
        job.run()
        job.append([["rm -f %s/%s.pval.signal.bedgraph" %(peaks_dirname,prefix)]])
        job.append([["rm -f %s/%s_treat_pileup.bdg %s_control_lambda.bdg" %(peaks_dirname,prefix,prefix)]])
        job.run()
        #===========================================
        # Generate bigWigs from beds to support trackhub visualization of peak files
        #============================================
        narrowPeak_bb_fname = utils.common.bed2bb('%s' %(narrowPeak_fn),
                                        '%s' %(self.chrom_sizes_name),
                                        '%s' %(self.narrowPeak_as_name),
                                         bed_type='bed6+4')
        gappedPeak_bb_fname = utils.common.bed2bb('%s' %(gappedPeak_fn),
                                        '%s' %(self.chrom_sizes_name),
                                        '%s' %(self.gappedPeak_as_name),
                                        bed_type='bed12+3')
        broadPeak_bb_fname =  utils.common.bed2bb('%s' %(broadPeak_fn),
                                        '%s' %(self.chrom_sizes_name),
                                        '%s' %(self.broadPeak_as_name),
                                        bed_type='bed6+3')

