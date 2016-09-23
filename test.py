import peak_calling as pc
import utils.job_runner as jr

p=pc.Macs2PeakCaller("./test/input1/psr_input1.00.bedpe.gz",
                     "./test/control1/control1.bedpe.gz",
                     "150",
                     "BEDPE",
                     "./test/macs",
                     "./common/hg38_sponges.chromInfo",
                     "hs")
p.run()
