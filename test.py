import peak_calling as pc
import utils.job_runner as jr

p=pc.Macs2PeakCaller("./test/input1/input1_reduced.bedpe.gz",
                     "./test/control1/control1_reduced.bedpe.gz",
                     "150",
                     "BEDPE",
                     "./test/macs",
                     "./common/hg38_sponges.chromInfo",
                     "hs")
p.run()
