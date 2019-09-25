import pprint
from pyVCFparallel.ParseVcfHeader import Parsers
from pyVCFparallel.ParallelReader import Reader

# headerInfo = Parsers.ParseHeader(filename="D:/NeoWork/PyVCF-P/pyVCFparallel/test/ALL_chr16_sample.vcf")
# reader = Reader(filename="D:/NeoWork/PyVCF-P/pyVCFparallel/test/ALL_chr16_sample.vcf", header=headerInfo);

# print(__name__+" => "+__file__)

if __name__ != "__mp_main__":
    filename = "E:/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"
#    filename = "D:/DATASETS/1000_genomes/data/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"
#    filename = "D:/NeoWork/PyVCF-P/pyVCFparallel/test/ALL_chr16_sample.vcf"
    reader = Reader(filename=filename)
    reader.run(4,300)  # 4 worker processes with 300 item communications buffers

    line = 0
    for record in reader:
        line += 1
#        print(str(line) + " >> " + str(record))

