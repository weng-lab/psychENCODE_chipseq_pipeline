from utils import job_runner

def map(fq1, fq2, bwa_index, prefix, cpu_num):
    """Returns SAM file name generated by mapping.

    """
    job = JobRunner()

    params = locals() # get all arguments as a dictionary

    raw_sam = "{}.raw.sam.gz".format(prefix)

    ## mapping reads to the genome
    params.setdefault({"raw_sam": raw_sam})
    bwa_cmd = \
      "bwa mem -M -t {cpu_num} {bwa_index} {fq1} {fq2} | gzip -c > {raw_sam}" 
    job.append([[ bwa_cmd.format(**params) ]])
    job.run()

    return raw_sam
