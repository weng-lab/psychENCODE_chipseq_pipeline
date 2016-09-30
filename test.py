import peak_calling as pc
import utils.job_runner as jr

p=pc.Macs2PeakCaller("./test/input1/input1.bedpe.gz",
                     "./test/control1/control1.bedpe.gz",
                     "150",
                     "BEDPE",
                     "./test/macs_rep",
                     "./common/hg38_sponges.chromInfo",
                     "hs")
p.run()

p=pc.Macs2PeakCaller("./test/input1/psr_input1.00.bedpe.gz",
                     "./test/control1/control1.bedpe.gz",
                     "150",
                     "BEDPE",
                     "./test/macs_psr00",
                     "./common/hg38_sponges.chromInfo",
                     "hs")
p.run()
p=pc.Macs2PeakCaller("./test/input1/psr_input1.01.bedpe.gz",
                     "./test/control1/control1.bedpe.gz",
                     "150",
                     "BEDPE",
                     "./test/macs_psr01",
                     "./common/hg38_sponges.chromInfo",
                     "hs")
p.run()
