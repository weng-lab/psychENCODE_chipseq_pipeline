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


class PseudorepsGenerator:

    def __init__(self, filepath,  prefix, out, threads, cpus ):
        self.prefix = prefix
        self.threads = threads
        self.input_file = filepath
        self.out = out
        self.cpus= cpus

    def run(self):
        job = jr.JobRunner()

        # Shuffle input
        print("Shuffle input")
        shuf_out = os.path.join(self.out,
                                "shuf_%s.bedpe" %(self.prefix))
        job.append([["zcat %s | shuf > %s" %(self.input_file,shuf_out)]])
        job.run()

        # Split into two files
        print("split in two")
        psr_prefix = os.path.join(self.out,"psr_%s." %(self.prefix))
        job.append([["split -d -nl/2 --additional-suffix=\".bedpe\"\
                      %s %s" %(shuf_out,psr_prefix)]])
        job.run()

        # TODO core to be set according to input parameters
        print("sort")
        job = jr.JobRunner()
        job.append([["sort --parallel=%s -S 2G \
                     -k1,1 -k2,2n %s00.bedpe | gzip -c > %s00.bedpe.gz" 
                     %(self.cpus, psr_prefix, psr_prefix)]])
        job.append([["sort --parallel=%s -S 2G \
                     -k1,1 -k2,2n %s01.bedpe | gzip -c > %s01.bedpe.gz" 
                     %(self.cpus, psr_prefix, psr_prefix)]])
        #job.append([["zcat %s | sort --parallel=%s -S 2G \
        #             -k1,1 -k2,2n -o %s && gzip -c %s" 
        #             %(self.input_file, self.cpus, self.input_file, self.input_file )]])
        job.run()
        
        # Clean
        job.append([["rm %s00.bedpe %s01.bedpe %s" 
                     %(psr_prefix,psr_prefix,shuf_out)]])
        job.run() 


    def find_names(self):
        #return [item for sublist in self.helper_names() 
        #        for item in sublist]
        return sum(self.helper_names(),[])

    def helper_names(self):
         for dirs in self.inputs:
             yield [each for each
                    in os.listdir(os.path.join(self.out,dirs))
                    if each.endswith('.bed')]

