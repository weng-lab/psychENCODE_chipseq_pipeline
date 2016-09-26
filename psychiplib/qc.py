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


import utils.job_runner as jr


class QC():
    def __init__(self, output_dir):
        self.output_dir = output_dir

    def get_mapstats(self, raw_bam, final_bam):
        job = jr.JobRunner()

        raw_bam_qc = raw_bam[:-3] + "flagstat.QC.txt"
        final_bam_qc = final_bam[:-3] + "flagstat.QC.txt"

        command = "samtools flagstat %s > %s"
        job.append([[ command %(raw_bam, raw_bam_qc) ]])
        job.append([[ command %(final_bam, final_bam_qc) ]])
        job.run()

        return (raw_bam_qc, final_bam_qc)


    def get_library_complexity(self, flt_bam, threads=1):
        job = jr.JobRunner()

        prefix = os.path.basename(flt_bam).split(".flt.")[0]

        tmp_bam = "%s/%s.tmp.bam" %(self.output_dir, prefix)
        qc_pbc = "%s/%s_libComplexity.QC.txt" %(self.output_dir, prefix)

        command = "samtools sort -n -@ %s " %(threads)
                   + "-o %s %s" %(tmp_bam, flt_bam)
        job.append([[ command ]])
        job.append([[ "mv %s %s" %(tmp_bam, flt_bam) ]])

        # Library Complexity
        # [1] TotalReadPairss
        # [2] DistinctReadPairs
        # [3] OneReadPairs
        # [4] TwoReadPairs
        # [5] NRF=Distinct/Total
        # [6] PBC1=OnePair/Distinct
        # [7] PBC2=OnePair/TwoPair
        command = """echo 1 | \
                     awk '{print "#Total\tDistinct\tOne\tTwo\tNRF\tPBC1\tPBC2"}' \
                     > %s""" %(qc_pbc)
        job.append([[ command ]])

        command = """bamToBed -bedpe -i %s | \
                     awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | \
                     grep -v 'chrM' | sort | uniq -c | \
                     awk 'BEGIN{mt=0;m0=0;m1=0;m2=0}
                               ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1}
                          END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'
                     >> %s""" %(flt_bam, qc_pbc)
        job.append([[ command ]])
        job.run()

        return qc_pbc


    def clean(self, final_bam):
        prefix = os.path.basename(final_bam).split(".final.")[0]
        raw_sam = "%s/%s.raw.sam.gz" %(self.output_dir, prefix)
        raw_bam = "%s/%s.raw.bam" %(self.output_dir, prefix)
        flt_bam = "%s/%s.flt.bam" %(self.output_dir, prefix)

        # Clean up intermediate files
        job = jr.JobRunner()
        job.append([[ "rm %s" %(raw_sam) ]])
        job.append([[ "rm %s" %(raw_bam) ]])
        job.append([[ "rm %s" %(flt_bam) ]])
        job.run()


