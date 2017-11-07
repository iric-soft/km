import unittest

import os
from km.utils import MutationFinder as umf
from km.utils.Jellyfish import Jellyfish


class kmMuttaionTest(unittest.TestCase):
    def testNPM1(self):
        seq_f = "./data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa"
        base_name = os.path.basename(seq_f)
        (ref_name, ext) = os.path.splitext(base_name)
        jf = Jellyfish("./data/jf/02H025_NPM1.jf", cutoff=0.05, n_cutoff=5)

        ref_seq = []
        for line in open(seq_f, "r"):
            line = line.strip()
            if line[0] == '>':
                continue
            ref_seq.append(line)
        ref_seq = ''.join(ref_seq)

        finder = umf.MutationFinder(
            ref_name, ref_seq, jf,
            False, 500
        )

        paths_quant = finder.get_paths_quant()
        paths_quant = paths_quant[len(paths_quant)-1]
        self.assertEqual(str(paths_quant.get_variant_name()),
                         "Insertion\t45:/TCTG:45",
                         "Test fail: NPM1 -> variant name")
        self.assertEqual(paths_quant.get_sequence(),
                         "CGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAG",
                         "Test fail: NPM1 -> sequence")
        # AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA

    def testFLT3(self):
        seq_f = "./data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa"
        base_name = os.path.basename(seq_f)
        (ref_name, ext) = os.path.splitext(base_name)
        jf = Jellyfish("./data/jf/03H116_ITD.jf", cutoff=0.05, n_cutoff=5)

        ref_seq = []
        for line in open(seq_f, "r"):
            line = line.strip()
            if line[0] == '>':
                continue
            ref_seq.append(line)
        ref_seq = ''.join(ref_seq)

        finder = umf.MutationFinder(
            ref_name, ref_seq, jf,
            False, 500
        )

        paths_quant = finder.get_paths_quant()
        paths_quant = paths_quant[len(paths_quant)-1]
        self.assertEqual(str(paths_quant.get_variant_name()),
                         "ITD\t204:/AACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACC:204",
                         "Test fail: FLT3-ITD -> variant name")
        self.assertEqual(paths_quant.get_sequence(),
                         "TACCTTCCCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTA",
                         "Test fail: FLT3-ITD -> sequence")
def runTests():
    unittest.main()


if __name__ == "__main__":
    runTests()
