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


def input_integrity(input_file):
    try:
        inputs = []
        prefix_inputs = []
        controls = []
        prefix_counts = []
        count = 1
        with open(input_file,"r") as f:
            for line in f:
                files = line.rstrip().split()
                if len(files) %2 != 0:
                    raise Exception("The input file is not well formatted. Line %s is missing files"%(str(count)))
                else:
                    test = 0
                    for filename in files:
                        test = test + os.path.isfile(filename)
                    
                    if test %2 != 0:
                        raise Exception("One or more files in line %s don't exist. Please check your input" %(str(count)))
                    
                    inputs.append([fields[0], fields[1]])
                    prefix_inputs.append([extract_prefix(fields[0]),extract_prefix(fields[1])])
                    
                    if len(fields) == 4 :
                        controls.append([fields[2], fields[3]])
                        prefix_controls.append([extract_prefix(fields[2]),extract_prefix(fields[3])])
                    count = count+1

    except Exception as e:
        raise(e)

    return inputs,controls,prefix_inputs,prefix_controls

def extract_prefix(filename):
    basename = os.path.basename(filename)
    ext = re.search(".f[ast]*q", basename)
    if ext:
        return basename[:ext.start()]
    else:
        return basename


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
        # Check if the input files exist and are well formatted
        os.path.isfile(args['f'])
#        inputs['input'], inputs['controls'], inputs['prefix_inputs'], inputs['prefix_controls'] = input_integrity(args['f'])

        # Check if the output folder exist
        if os.path.isdir(args['o']):
            raise Exception("mkdir: cannot create directory '%s': File exists." %(args['o']))
        else :
            inputs['outdir'] = args['o']

        # Check if Bowtie index exists
        if args['b'].endswith('.bwt') :
            prefix_index = os.path.splitext(args['b'])[0]
        else :
            prefix_index = args['b']
        
        for suffix in ['.bwt', '.amb', '.ann', '.pac', '.sa'] :
            if not os.path.isfile("%s%s" %(prefix_index,suffix)):
                raise Exception("The BWA index file '%s%s' does not exist." %(prefix_index,suffix))
        inputs['index'] = prefix_index

        # Check if genome file exists 
        if not os.path.isfile(args['g']):
            raise Exception("The genome chromInfo file '%s' does not exist." %(args['g']))
        inputs['genome'] = args['g']
        
        if args['sponges'] != None :
            if not os.path.isfile(args['sponges']):
                raise Exception("The sponges file '%s' does not exist." %(args['sponges']))
            inputs['sponges'] = args['sponges']
        else:
            inputs['sponges'] = None
        
        inputs['threads'] = args['t']
        inputs['cpus'] = args['c']
        
        
        
        # If not enough controls are supplied, use pooled controls.
        if len(inputs['inputs']) != len(inputs['controls']):
            inputs['pooled'] = True
        else:
            inputs['pooled'] = False

    except Exception as e:
        raise(e)
    


    return inputs


def main(inputs):
    job = jr.JobRunner(cpus = inputs['threads'])
    
    print("--- Welcome to the pipeline ---")
    if inputs['pooled']:
        print("Warning : The controls libraries will be pooled together.")
        for idx in len(inputs['inputs']):
            job.append([["python run_psychip.py map %s inputs" %(idx)]])
        for idx in len(inputs['controls']):
            job.append([["python run_psychip.py map %s controls" %(idx)]])
        job.run()

        job.append([["python run_psychip.py pooling"]])
        job.run()

        for idx in len(inputs['inputs']):
            job.append([["python run_psychip.py peak %s" %(idx)]])
        job.run()
    else:
        for idx in len(inputs['inputs']):
            job.append([["python run_psychip.py complete %s" %(n)]])
        job.run()


if __name__ == "__main__":
    usage = """psychip [options] -b <bowtie_index> -g <genome_chrInfo> -f <input_file> -o <output_folder>"""
    description = "Histone ChIP-seq processing pipeline based on \
                   the ENCODE(phase-3) with support for paired-end \
                   reads and multiple replicates."

    parser = argparse.ArgumentParser(usage = usage, description=description, prog="psychip")
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-f', metavar = 'FILE', required = True, help = "Input file")
    parser.add_argument('-o', metavar = 'FOLDER', required = True, help = "Output file")
    parser.add_argument('-b', metavar = 'FILE', required = True, \
                        help = "BWA index file")
    parser.add_argument('-g', metavar = 'FILE', required = True,  help = "ChromInfo file for the genome")
    parser.add_argument('-t', metavar = 'INT', type = int,  default=1, help = "Number of threads")
    parser.add_argument('-c', metavar = 'INT', type = int,  default=1, \
                        help = """Number of cpus per thread. \
                               Total number of cpus used = \
                               #threads * #cpus.""")
    parser.add_argument('--sponges', metavar = 'FILE', help = "Sponges file used to \
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
