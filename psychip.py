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
import re

def input_integrity(input_file):
    try:
        inputs = []
        prefix_inputs = []
        controls = []
        prefix_controls= []
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
                    
                    inputs.append([files[0], files[1]])
                    prefix_inputs.append(extract_prefix(files[0]))
                    
                    if len(files) == 4 :
                        controls.append([files[2], files[3]])
                        prefix_controls.append(extract_prefix(files[2]))
                    count = count+1
        print("Input Integrity Verification Completed")
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
        checker.which('R')
        #check sort version
        #check samtools version
    except Exception as e:
        raise(e) 
    

def check_arguments(args):
    checker=u.Utils()
    inputs = {}
    try:
        # Check if the input files exist and are well formatted
        os.path.isfile(args['f'])
        inputs['inputs'], inputs['controls'], inputs['prefix_inputs'], inputs['prefix_controls'] = input_integrity(args['f'])
        
        # Check if the output folder exist
        if os.path.isdir(args['o']):
            raise Exception("mkdir: cannot create directory '%s': File exists." %(args['o']))
        else :
            os.makedirs(args['o'])
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

        if args['macs2'] != None:
            if os.path.exists(args['macs2']):
                inputs['macs2'] = args['macs2']
            else:
                raise Exception("Macs2 not found in '%s'" %(args['macs2']))
        else:
            try:
                checker.which('macs2')
                inputs['macs2'] = "macs2"
            except Exception as e:
                raise(e)

#        if not os.path.exists(args['picard']):
#            raise Exception("Picard not found in '%s'" %(args['picard']))
#        else:
#            inputs['picard'] = args['picard']
        inputs['fraglen'] = -1
        if args['fraglen'] != None:
            if args['spp']:
                raise Exception("'--spp' and '--fraglen' are mutually exclusive. Please select only one.")
            else:
                inputs['fraglen'] = args['fraglen']
        else:
            if args['spp']:
                try:
                    j = jr.JobRunner()
                    j.append([['Rscript ./common/test.R']])
                    j.run()
                    inputs['spp'] = True
                except:
                   raise Exception("Spp is not currently installed or cannot be loaded.")
            else:
                inputs['spp'] = False
                inputs['fraglen'] = -1
                print('Picard tools will be used to estimate insert size.')
                
        
        inputs['threads'] = args['t']
        inputs['cpus'] = args['c']
        inputs['common'] = os.path.join(os.path.dirname(os.path.abspath(__file__)),"common") 
        inputs['picard'] = os.path.join(os.path.dirname(os.path.abspath(__file__)),"common/picard-tools-1.141/picard.jar") 
        
        # If not enough controls are supplied, use pooled controls.
        if len(inputs['inputs']) != len(inputs['controls']):
            inputs['pooled'] = True
        else:
            inputs['pooled'] = False
        
        print("Argument Verification completed")
    except Exception as e:
        raise(e)
    


    return inputs


def main(inputs):
    with open('%s/inputs.json' %(inputs['outdir']), 'w') as outfile:
        json.dump(inputs, outfile)

    job = jr.JobRunner(cpus = inputs['threads'])
    
    print("--- Welcome to the pipeline ---")
    if inputs['pooled']:
        print("Warning : The controls libraries will be pooled together.")
        
        for idx in range(len(inputs['inputs'])):
            job.append([["python run_psychip_pool_input.py %s/inputs.json map %s inputs" %(inputs['outdir'], idx)]])
        for idx in range(len(inputs['controls'])):
            job.append([["python run_psychip_pool_input.py %s/inputs.json map %s controls" %(inputs['outdir'], idx)]])
        job.run()
        
        job.append([["python run_psychip_pool_input.py %s/inputs.json pooling" %(inputs['outdir'])]])
        job.run()
        for idx in range(len(inputs['inputs']) + 1): # the plus one is for the pooled input
            job.append([["python run_psychip_pool_input.py %s/inputs.json peak %s" %(inputs['outdir'], idx)]])
        job.run()
    else:
        for idx in range(len(inputs['inputs'])):
            job.append([["python run_psychip.py %s/inputs.json %s" %(inputs['outdir'], idx)]])
        job.run()
        job.append([["python run_psychip_pool_input.py %s/inputs.json pooling" %(inputs['outdir'])]])
        job.run()
        job.append([["python run_psychip_pool_input.py %s/inputs.json peak %s" %(inputs['outdir'], len(inputs['inputs'])+1)]])
        job.run()


if __name__ == "__main__":
    usage = """psychip [options] -b <bowtie_index> -g <genome_chrInfo> -f <input_file> -o <output_folder>"""
    description = "Histone ChIP-seq processing pipeline based on \
                   the ENCODE with support for paired-end \
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
                        help = "Number of cpus per thread. \
                               Total number of cpus used = \
                               #threads * #cpus.")
    parser.add_argument('--sponges', metavar = 'FILE', help = "Sponges file used to \
                                             remove sequences")
    parser.add_argument('--spp', action='store_true', help = "Compute Metrics for input file using SPP package (requires SPP package)")
    parser.add_argument('--fraglen', metavar = 'INT', help = "Insert size for paired end reads")
    parser.add_argument('--macs2', metavar = 'FILE', help = "Location of macs2 executable e.g. (/home/user/macs2/bin/macs2)")
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
