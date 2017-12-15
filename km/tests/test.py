# cmd to run test: coverage run -m unittest discover ./km/tests

import unittest

import os
import sys

from argparse import Namespace
from km.tools import find_mutation as fm
from km.utils.Jellyfish import Jellyfish
from km.utils import MutationFinder as umf

from contextlib import contextmanager
from StringIO import StringIO


@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class kmMuttaionTest(unittest.TestCase):
    def testNPM1(self):
        args = Namespace(
            count=5,
            graphical=False,
            jellyfish_fn='./data/jf/02H025_NPM1.jf',
            ratio=0.05,
            steps=500,
            target_fn=["./data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa"],
            verbose=False
        )

        with captured_output() as (out, err):
            fm.main_find_mut(args, None)

        output = out.getvalue().split("\n")
        output = output[16].split("\t")

        self.assertEqual(output[2],
                         "Insertion",
                         "Test fail: NPM1 -> type")
        self.assertEqual(output[3],
                         "45:/TCTG:45",
                         "Test fail: NPM1 -> variant")
        self.assertEqual(output[7],
                         "CGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAG",
                         "Test fail: NPM1 -> sequence")
        # AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA

    def testFLT3_ITD(self):
        args = Namespace(
            count=5,
            graphical=False,
            jellyfish_fn='./data/jf/03H116_ITD.jf',
            ratio=0.05,
            steps=500,
            target_fn=["./data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa"],
            verbose=False
        )

        with captured_output() as (out, err):
            fm.main_find_mut(args, None)

        output = out.getvalue().split("\n")
        output = output[16].split("\t")

        self.assertEqual(output[2],
                         "ITD",
                         "Test fail: FLT3-ITD -> type")
        self.assertEqual(output[3],
                         "204:/AACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACC:204",
                         "Test fail: FLT3-ITD -> variant")
        self.assertEqual(output[7],
                         "TACCTTCCCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTA",
                         "Test fail: FLT3-ITD -> sequence")

    def testFLT3_TKD(self):
        args = Namespace(
            count=5,
            graphical=False,
            jellyfish_fn='./data/jf/05H094_FLT3-TKD_del.jf',
            ratio=0.05,
            steps=500,
            target_fn=["./data/catalog/GRCh38/FLT3-TKD_exon_20.fa"],
            verbose=False
        )

        with captured_output() as (out, err):
            fm.main_find_mut(args, None)

        output = out.getvalue().split("\n")
        output = output[16].split("\t")

        self.assertEqual(output[2],
                         "Deletion",
                         "Test fail: FLT3-TKD -> type")
        self.assertEqual(output[3],
                         "32:gat/:35",
                         "Test fail: FLT3-TKD -> variant")
        self.assertEqual(output[7],
                         "TGCCCCTGACAACATAGTTGGAATCACTCATATCTCGAGCCAATCCAAAGTCACATATCTT",
                         "Test fail: FLT3-TKD -> sequence")

    def testDNMT3A(self):
        args = Namespace(
            count=5,
            graphical=False,
            jellyfish_fn="./data/jf/02H033_DNMT3A_sub.jf",
            ratio=0.05,
            steps=500,
            target_fn=["./data/catalog/GRCh38/DNMT3A_R882_exon_23.fa"],
            verbose=False
        )

        with captured_output() as (out, err):
            fm.main_find_mut(args, None)

        output = out.getvalue().split("\n")
        output = output[16].split("\t")

        self.assertEqual(output[2],
                         "Substitution",
                         "Test fail: DNMT3A -> type")
        self.assertEqual(output[3],
                         "33:c/T:34",
                         "Test fail: DNMT3A -> variant")
        self.assertEqual(output[7],
                         "TGACCGGCCCAGCAGTCTCTGCCTCGCCAAGTGGCTCATGTTGGAGACGTCAGTATAGTGGA",
                         "Test fail: DNMT3A -> sequence")

    def test_not_linear(self):
        ref_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        k_len = 31
        ref_name = "not_linear"
        jf = Jellyfish("./data/jf/02H033_DNMT3A_sub.jf")

        finder = umf.MutationFinder("", "", jf, False, 500)

        with self.assertRaises(ValueError):
            finder.get_ref_kmer(ref_seq, 31, ref_name)


def runTests():
    unittest.main()


if __name__ == "__main__":
    runTests()
