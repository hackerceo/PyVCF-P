import time
import pprint
from pyVCFparallel.Parsers import Parsers
from pyVCFparallel.ParallelReader import Reader

# headerInfo = Parsers.ParseHeader(filename="D:/NeoWork/PyVCF-P/pyVCFparallel/test/ALL_chr16_sample.vcf")
# reader = Reader(filename="D:/NeoWork/PyVCF-P/pyVCFparallel/test/ALL_chr16_sample.vcf", header=headerInfo);

# print(__name__+" => "+__file__)

if __name__ != "__mp_main__":
#    filename = "E:/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"
#    filename = "D:/DATASETS/1000_genomes/data/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"
    filename = "D:/NeoWork/PyVCF-P/pyVCFparallel/test/ALL_chr16_sample.vcf"
#    filename = "D:/DATASETS/dbSNP/human_9606_b151_GRCh37p13/head.vcf"
    reader = Reader(filename=filename)
    reader.run(5,100)  # 4 worker processes with 30 item communications buffers

    line = 0
    time_start = time.time()
    time_last_tick = 0


# TODO STUFF HERE:
# handle ##PEDIGREE and ##pedigreeDB

    variant_db = {}

    for record in reader:
        line = line + 1
        if line % 1000 == 0:
            exec_time = time.time() - time_start
            print(str(line).ljust(10) + " records @ " + str( line/exec_time ) + " records/sec avg.  Instantaneous: " + str(1000 / (exec_time - time_last_tick)) + " records/sec")
            time_last_tick = exec_time
#        if line > 10: break
#        pp = pprint.PrettyPrinter(indent=4)
#        pp.pprint(record)

        # for key in record["alt"]:
        #     idx = record["ref"][0] + ">" + key
        #     if idx in variant_db:
        #         variant_db[idx] = variant_db[idx] + 1
        #     else:
        #         variant_db[idx] = 1

    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(variant_db)

#
# # =================== temporary testing of an idea jason gave me
# import random
# from bitarray import bitarray
# import time
#
# cntVariants = 360000
# cntPatients = 2504
#
# def generateRandPatientArray(total_patients):
#     ret = total_patients * bitarray('0')
#     # how many patients are affected?
#     pt_to_set = int(random.gammavariate(2.3,5))
#     # now set that many patients as having mutation
#     for pt in range(0, pt_to_set):
#         ret[random.randrange(0, total_patients)] = True
#     return ret
#
# # create the patient arrays for all variants
# variantArrays = [None] * cntVariants
# for idx in range(len(variantArrays)):
#     variantArrays[idx] = generateRandPatientArray(cntPatients)
#
#
# # pick the first variant and compute counts for the entire row intersection
# variant_counts = cntVariants * [0]
# variant_row_bitmap = variantArrays[0]
#
# #start timer
# time_start = time.clock()
#
# # compute the counts
# for idx in range(len(variantArrays)):
#     variant_counts[idx] = (variant_row_bitmap & variantArrays[idx]).count()
#
# # display run time in seconds
# print(time.clock() - time_start)
#
#
#
#
