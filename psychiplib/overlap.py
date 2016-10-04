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

import os,sys

import utils.job_runner as jr


class Overlap:

    def __init__(self, arr_inputs, out_folder):
        self.input_file = arr_inputs[0]
        self.out = out_folder
        #self.names = self.find_names()

    def run(self):
        job=jr.JobRunner()
        # narrowPeak
        np_oracle = "%s/macs_%s/%s.narrowPeak.gz" %(self.out,self.input_file,self.input_file)
        np_psr1 = "%s/macs_%s_psr_00/psr_%s.00.narrowPeak.gz"  %(self.out,self.input_file,self.input_file)
        np_psr2 = "%s/macs_%s_psr_01/psr_%s.01.narrowPeak.gz"  %(self.out,self.input_file,self.input_file)

        # broadPeak
        bp_oracle = "%s/macs_%s/%s.broadPeak.gz"  %(self.out,self.input_file,self.input_file)
        bp_psr1 = "%s/macs_%s_psr_00/psr_%s.00.broadPeak.gz"  %(self.out,self.input_file,self.input_file)
        bp_psr2 = "%s/macs_%s_psr_01/psr_%s.01.broadPeak.gz"  %(self.out,self.input_file,self.input_file)

        # gappedPeak
        gp_oracle = "%s/macs_%s/%s.gappedPeak.gz"  %(self.out,self.input_file,self.input_file)
        gp_psr1 = "%s/macs_%s_psr_00/psr_%s.00.gappedPeak.gz"  %(self.out,self.input_file,self.input_file)
        gp_psr2 = "%s/macs_%s_psr_01/psr_%s.01.gappedPeak.gz"  %(self.out,self.input_file,self.input_file)

        job.append([["bedtools intersect \
                     -a %s -b %s -f 0.50 -F 0.50 -e |\
                     bedtools intersect \
                     -a stdin -b %s -f 0.50 -F 0.50 -e > \
                     %s/%s.final.narrowPeak.gz"\
                     %(np_oracle,np_psr1,np_psr2, \
                       self.out, self.input_file)]])

        job.append([["bedtools intersect \
                     -a %s -b %s -f 0.50 -F 0.50 -e |\
                     bedtools intersect \
                     -a stdin -b %s -f 0.50 -F 0.50 -e > \
                     %s/%s.final.broadPeak.gz"\
                     %(bp_oracle,bp_psr1,np_psr2, \
                       self.out, self.input_file)]])
        job.append([["bedtools intersect \
                     -a %s -b %s -f 0.50 -F 0.50 -e |\
                     bedtools intersect \
                     -a stdin -b %s -f 0.50 -F 0.50 -e > \
                     %s/%s.final.gappedPeak.gz"\
                     %(gp_oracle,gp_psr1,np_psr2, \
                       self.out, self.input_file)]])
        job.run()

    def find_names(self):
        #return [item for sublist in self.helper_names()
        #        for item in sublist]
        return sum(self.helper_names(),[])

