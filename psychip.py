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

import os,argparse,sys,json
import utils.utils as u
import utils.job_runner as jr


def input_integrity(filename):
    try:
        inputs=[]
        controls=[]
        count = 1
        with open(filename,"r") as f:
            for line in f:
                fields = line.rstrip().split()
                if len(fields) %2 != 0:
                    raise Exception("The input file is not well formatted. Line %s is missing files"%(str(count)))
                else:
                    test = os.path.isfile(fields[0])
                    test = test + os.path.isfile(fields[1])
                    test = test + os.path.isfile(fields[2])
                    test = test + os.path.isfile(fields[3])
                    if test %2 != 0:
                        raise Exception("One or more files in line %s don't exist. Please check your input" %(str(count)))
                    
                    inputs.append([fields[0], fields[1]])
                    controls.append([fields[2], fields[3]])
                    count = count+1

    except Exception as e:
        raise(e)

    return inputs,controls

def check_requirements():
    checker=u.Utils()
    try:
        checker.which('bwa')
        checker.which('bedtools')
        checker.which('samtools')
        #check sort version
        #check samtools version
        #check ('java -jar picard.jar') #How to check this properly?
    except Exception as e:
        raise(e) 
    

def check_arguments(args):
    inputs = {}
    try:
        os.path.isfile(args['f'])
        inputs['input'],inputs['controls']=input_integrity(args['f'])
        print(inputs)
    except Exception as e:
        raise(e)


    return inputs


def main(inputs):
    print("welcome to the pipeline")
    #runpipe(json)


if __name__ == "__main__":
    usage = """psychip [options] -b <bowtie_index> -g <genome_chrInfo> -f <input_file> -o <output_folder>"""
    description = "Histone ChIP-seq processing pipeline based on \
                   the ENCODE(phase-3) with support for paired-end \
                   reads and multiple replicates."

    parser = argparse.ArgumentParser(usage = usage, description=description, prog="psychip")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-f', metavar='FILE', required = True, help = "Input file")
    parser.add_argument('-o', metavar='FOLDER', required = True, help = "Output file")
    parser.add_argument('-i', metavar='FILE', required = True, \
                        help = "BWA index file")
    parser.add_argument('-g', metavar='FILE', required = True,  help = "ChromInfo file for the genome")
    parser.add_argument('-t', metavar='INT', default=1, help = "Number of threads")
    parser.add_argument('-c', metavar='INT', default=1, \
                        help = """Number of cpus per thread. \
                               Total number of cpus used = \
                               #threads * #cpus.""")
    parser.add_argument('--sponges', metavar='FILE', help = "Sponges file used to \
                                             remove sequences")
    #parser.add_argument('--single', help = "If the library is \
    #                                        single-end")
    #parser.add_argument('--cluster', help = "#developping")
    #parser.add_argument('-d', help = "#developping \
    #                                  Delete large files not used")
    #parser.add_argument('--caller', help = "#developping \
    #                                        use your own peak caller")

    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    
    try:
        check_requirements()
    except Exception as e:
        print(e)
        sys.exit(1)

    inputs = {}

    try:
        inputs=check_arguments(vars(args))
    except Exception as e:
        print(e)
        sys.exit(1)

    sys.exit(main(inputs))
