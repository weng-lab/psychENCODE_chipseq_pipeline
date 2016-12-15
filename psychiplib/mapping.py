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
import re
import utils.job_runner as jr

class BwaMapper():
    def __init__(self, prefix, genome_index, output_dir, fastq1, fastq2=""):
        self.fastq1 = fastq1
        self.fastq2 = fastq2 # empty string = single-end
        self.genome_index = genome_index
        self.prefix = prefix
        self.output_dir = output_dir


    def run(self, cpu_num):
        if cpu_num <= 0:
            raise Exeption("The number of CPU must be > 0!")

        job = jr.JobRunner()

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        sam_output = "%s/%s.raw.sam.gz" %(self.output_dir, self.prefix)

        command = "bwa mem -M -t %s %s " %(cpu_num, self.genome_index) \
                   + "%s %s " %(self.fastq1, self.fastq2) \
                   + "| gzip -c > %s" %(sam_output)

        job.append([[command]])
        job.run()
