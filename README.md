## Summary

psychip pipeline is an improved version of the ENCODE pipeline for histone marks developed for the psychENCODE project.

The orginal specs for the ENCODE pipeline can be found here : https://goo.gl/KqHjKH

## Dependencies

NOTE: The current package does not run on a MacOS or Windows.
* unix, gnu utilscore (>=8.20), R (>=3.1), awk, samtools (>1.0), java, bedtools (>=2.25), bwa, python (v2.7)
* modified version of macs2 working with BEDPE inputs. Available at https://github.com/weng-lab/MACS
Optional but highly recommended:
* R packages: spp (>=1.14) (https://github.com/kundajelab/phantompeakqualtools/blob/master/README.md)

## Introduction

psychip consists of a main script calling different module performing alignment, filtering, quality control, pseudo-replicates creations, peak calling and overlapping.



## General Usage

Usage: `python psychip.py -f test/inputfile -o outputdir -b /home/path/to/hg38 -g /home/path/to/hg38.chromInfo <options>`

1. Required:
    * -f inputfile : this file contains the path to the libraries you want to processing
    
    * -o outputdir : the directory where the output will be written
    * -b bwaIndex  : bwa index for the genome
    * -g chromInfo : genome file which simply lists the names of the chromosomes and their size (in basepairs). Genome files must be tab-delimited.
                     (Exmaple: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz)
    
2. Optional:
    * -t threads   : number of threads to use
    * -c cpus      : number of cpus per threads
    * --sponges    : names of mapping targets of sequences sponges file to remove the highly repetitive as presented here :http://nar.oxfordjournals.org/content/43/20/e133
                     (an example of the file can be found under : common/GRCh38_sponges.names)
    * --macs2      : path to the macs2 bin if not in the $PATH variable
    * --spp        : use spp package to compute library complexity, quality controls and estimate fragment size 
    * --fraglen    : instead of using spp package to estimare the insert fragment size, supply it manually.
                    If either spp nor fraglen are supplied, picard will be used to estimate the insert fragment size.


## Test

You can test the pipeline by entering the main project folder and run
```
$ python psychip.py -f test/inputfile2 -o .test/output -b your_bwa_index -g your_genome_file
```

## Input file format


1. **Inputfile**: The input file contains the path to all the librarie and controls . It contains 4 tab delimited columns.

   | ip.R1 | ip.R2 | input1.R1 | input2.R2 |
   |-----------|-----------|-------------|-------------|
   |  /home/path/to/syn00001.R1.fastq.gz | /home/path/to/syn00001.R2.fastq.gz | /home/path/to/syn00001.input.R1.fastq.gz | /home/path/to/syn00001.input.R1.fastq.gz |
   |  /home/path/to/syn00002.R1.fastq.gz | /home/path/to/syn00002.R2.fastq.gz | /home/path/to/syn00002.input.R1.fastq.gz | /home/path/to/syn00002.input.R1.fastq.gz |
   |  /home/path/to/syn00003.R1.fastq.gz | /home/path/to/syn00003.R2.fastq.gz | /home/path/to/syn00003.input.R1.fastq.gz | /home/path/to/syn00003.input.R1.fastq.gz |
   |  /home/path/to/syn00004.R1.fastq.gz | /home/path/to/syn00004.R2.fastq.gz | /home/path/to/syn00004.input.R1.fastq.gz | /home/path/to/syn00004.input.R1.fastq.gz |
   
In the above scenario each ip library has its own input that will be used during the peak calling procedure. In case controls are not available for all the ips it is possible to omit the last two columns. If so, the all the inputs will be pooled togehter and used for all the ips

    
   | ip.R1 | ip.R2 | input1.R1 | input2.R2 |
   |-----------|-----------|-------------|-------------|
   |  /home/path/to/syn00001.R1.fastq.gz | /home/path/to/syn00001.R2.fastq.gz | /home/path/to/syn00001.input.R1.fastq.gz | /home/path/to/syn00001.input.R1.fastq.gz |
   |  /home/path/to/syn00002.R1.fastq.gz | /home/path/to/syn00002.R2.fastq.gz | /home/path/to/syn00002.input.R1.fastq.gz | /home/path/to/syn00002.input.R1.fastq.gz |
   |  /home/path/to/syn00003.R1.fastq.gz | /home/path/to/syn00003.R2.fastq.gz | | |
   |  /home/path/to/syn00004.R1.fastq.gz | /home/path/to/syn00004.R2.fastq.gz | | |
   
In the above example input libraries are missing for syn00003 and syn00004. Therefore, syn00001.input and syn00002.control will be pooled together and used for all the ips.
   
   
## Output file formats

Considering as inputfile for psychip:

    | ip.R1 | ip.R2 | input1.R1 | input2.R2 |
    |-----------|-----------|-------------|-------------|
    |  /home/path/to/syn00001.R1.fastq.gz | /home/path/to/syn00001.R2.fastq.gz | /home/path/to/syn00001.input.R1.fastq.gz | /home/path/to/syn00001.input.R1.fastq.gz |
    |  /home/path/to/syn00002.R1.fastq.gz | /home/path/to/syn00002.R2.fastq.gz | /home/path/to/syn00002.input.R1.fastq.gz | /home/path/to/syn00002.input.R1.fastq.gz |


psychip will produce for each library 3 files containing narrow peaks, gapped peaks and broad peaks.

e.g.

/home/path/to/outdir/syn00001.R1.narrowPeak\
/home/path/to/outdir/syn00001.R1.broadPeak\
/home/path/to/outdir/syn00001.R1.gappedPeak\
/home/path/to/outdir/syn00002.R1.narrowPeak\
/home/path/to/outdir/syn00002.R1.broadPeak\
/home/path/to/outdir/syn00002.R1.gappedPeak


## Future improvement
* Add support for other organisms (right now works only with human)
* Add checkpoint to re-run just a part of the pipeline
* Use BAMPE input instead of BEDPE in order to use other peak caller
* Add check for sort and samtools versions

## References


## Contributors
* Eugenio Mattei - Postdoctoral Fellow, Weng Lab, UMASS Medical School, Worcester(MA)
* Junko Tsuji - Postdoctoral Fellow, Weng Lab, UMASS Medical School, Worcester (MA)
