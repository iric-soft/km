#!/usr/bin/env python

# cmd to run test: pytest ./km --cov=km/

import unittest

import os
import sys

from argparse import Namespace
from km.tools import find_mutation as fm
from km.tools import find_report as fr
from km.tools import linear_kmin as lk
from km.tools import min_cov as mc

from km.utils.Jellyfish import Jellyfish
from km.utils import MutationFinder as umf
from km.utils import common as uc

from contextlib import contextmanager
from io import StringIO


@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class TestkmMuttaion(unittest.TestCase):
    def test_NPM1(self):
        target = "./data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa"
        args = Namespace(
            count=5,
            graphical=False,
            jellyfish_fn='./data/jf/02H025_NPM1.jf',
            ratio=0.05,
            steps=500,
            branchs=10,
            target_fn=[target],
            verbose=False
        )

        with captured_output() as (out, err):
            fm.main_find_mut(args, None)

        fm_output = out.getvalue()
        find_output = fm_output.split("\n")
        find_output = find_output[11].split("\t")

        self.assertEqual(find_output[2],
                         "Insertion",
                         "Test fail: NPM1 -> find type")
        self.assertEqual(find_output[3],
                         "45:/TCTG:45",
                         "Test fail: NPM1 -> find variant")
        self.assertEqual(find_output[8],
                         "CGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAG",
                         "Test fail: NPM1 -> find sequence")

        args = Namespace(
            target=target,
            infile=StringIO(fm_output),
            info="vs_ref",
            min_cov=1,
            exclu="",
            format=None
        )

        with captured_output() as (out, err):
            fr.main_find_report(args, None)

        output = out.getvalue()
        report_output = output.split("\n")
        report_output = report_output[2].split("\t")

        self.assertEqual(report_output[2],
                         "chr5:171410544",
                         "Test fail: NPM1 -> report pos")
        self.assertEqual(report_output[3],
                         "ITD",
                         "Test fail: NPM1 -> report type")
        self.assertEqual(report_output[4],
                         "0",
                         "Test fail: NPM1 -> report Removed")
        self.assertEqual(report_output[5],
                         "4 | 1",
                         "Test fail: NPM1 -> report Added")
        self.assertEqual(report_output[6],
                         "2870.6",
                         "Test fail: NPM1 -> report Abnormal")
        self.assertEqual(report_output[7],
                         "3055.2",
                         "Test fail: NPM1 -> report Normal")
        self.assertEqual(report_output[8],
                         "0.484",
                         "Test fail: NPM1 -> report rVAF")
        self.assertEqual(report_output[9],
                         "2428",
                         "Test fail: NPM1 -> report Min_coverage")
        self.assertEqual(report_output[11],
                         "/TCTG",
                         "Test fail: NPM1 -> report variant")
        self.assertEqual(report_output[14],
                         "AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA",
                         "Test fail: NPM1 -> report sequence")

        args = Namespace(
            target=target,
            infile=StringIO(fm_output),
            info="vs_ref",
            min_cov=1,
            exclu="",
            format="vcf"
        )

        with captured_output() as (out, err):
            fr.main_find_report(args, None)

        output = out.getvalue()
        report_output = [o for o in output.split("\n") if o and o[0] != "#"]
        report_output = report_output[0].split("\t")

        self.assertEqual(report_output[1],
                         "171410539",
                         "Test fail: NPM1 -> vcf report pos")
        self.assertEqual(report_output[3],
                         "CTCTGG",
                         "Test fail: NPM1 -> vcf report ref")
        self.assertEqual(report_output[4],
                         "CTCTGTCTGG",
                         "Test fail: NPM1 -> vcf report alt")

    def test_FLT3_IandI(self):
        target = "./data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa"
        args = Namespace(
            count=5,
            graphical=False,
            jellyfish_fn='./data/jf/03H112_IandI.jf',
            ratio=0.05,
            steps=500,
            branchs=10,
            target_fn=[target],
            verbose=False
        )

        with captured_output() as (out, err):
            fm.main_find_mut(args, None)

        fm_output = out.getvalue()
        find_output = fm_output.split("\n")
        find_output = find_output[11].split("\t")

        self.assertEqual(find_output[2],
                         "ITD",
                         "Test fail: FLT3-IandI -> find type")
        self.assertEqual(find_output[3],
                         "152:/TCTTGCGTTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCAAACTCTAAATTTTCTCTTGGAA" +
                         "ACTCCCATTTGAGATCATATTC:152",
                         "Test fail: FLT3-IandI -> find variant")
        self.assertEqual(find_output[8],
                         "TTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCGTTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTC" +
                         "CCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCTCTTGCGTTCATCACTTTTCCAAAAGCACC" +
                         "TGATCCTAGTACCTTCCCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAAT" +
                         "CAACGTAGAAGTACTC",
                         "Test fail: FLT3-IandI -> find sequence")

        args = Namespace(
            target=target,
            infile=StringIO(fm_output),
            info="vs_ref",
            min_cov=1,
            exclu="",
            format=None
        )

        with captured_output() as (out, err):
            fr.main_find_report(args, None)

        output = out.getvalue()
        report_output = output.split("\n")
        report_output = report_output[2].split("\t")

        self.assertEqual(report_output[2],
                         "chr13:28034128",
                         "Test fail: FLT3-IandI -> report pos")
        self.assertEqual(report_output[3],
                         "I&I",
                         "Test fail: FLT3-IandI -> report type")
        self.assertEqual(report_output[4],
                         "0",
                         "Test fail: FLT3-IandI -> report Removed")
        self.assertEqual(report_output[5],
                         "93 | 1",
                         "Test fail: FLT3-IandI -> report Added")
        self.assertEqual(report_output[6],
                         "399.1",
                         "Test fail: FLT3-IandI -> report Abnormal")
        self.assertEqual(report_output[7],
                         "398.5",
                         "Test fail: FLT3-IandI -> report Normal")
        self.assertEqual(report_output[8],
                         "0.500",
                         "Test fail: FLT3-IandI -> report rVAF")
        self.assertEqual(report_output[9],
                         "285",
                         "Test fail: FLT3-IandI -> report Min_coverage")
        self.assertEqual(report_output[11],
                         "/TCTTGCGTTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCAAACTCTAAATTTTCTCTTGGAAACTCC" +
                         "CATTTGAGATCATATTC",
                         "Test fail: FLT3-IandI -> report variant")
        self.assertEqual(report_output[14],
                         "CTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCGTTCATCACTTTT" +
                         "CCAAAAGCACCTGATCCTAGTACCTTCCCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCTCT" +
                         "TGCGTTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCAAACTCTAAATTTTCTCTTGGAAACTCCCATT" +
                         "TGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGT" +
                         "AGCTGGCTTTCATACCTAAATTGCTTTTTGTACTTGTGACAAATTAGCAGGGTTAAAACGACAATGAAGAGGAGACA" +
                         "AACACCAATTGTTGCATAGAATGAGATGTTGTCTTGGATGAAAGGGAAGGGGC",
                         "Test fail: FLT3-IandI -> report sequence")

        args = Namespace(
            target=target,
            infile=StringIO(fm_output),
            info="vs_ref",
            min_cov=1,
            exclu="",
            format="vcf"
        )

        with captured_output() as (out, err):
            fr.main_find_report(args, None)

        output = err.getvalue()
        find_output = [o for o in output.split("\n") if o[:5] == "NOTE:"]
        find_output = find_output[0]

        self.assertEqual(find_output,
                         "NOTE: Mutation overlaps 2 exons or more, VCF output is disabled ")

    def test_FLT3_ITD(self):
        target = "./data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa"
        args = Namespace(
            count=5,
            graphical=False,
            jellyfish_fn='./data/jf/03H116_ITD.jf',
            ratio=0.05,
            steps=500,
            branchs=10,
            target_fn=[target],
            verbose=False
        )

        with captured_output() as (out, err):
            fm.main_find_mut(args, None)

        fm_output = out.getvalue()
        find_output = fm_output.split("\n")
        find_output = find_output[11].split("\t")

        self.assertEqual(find_output[2],
                         "ITD",
                         "Test fail: FLT3-ITD -> find type")
        self.assertEqual(find_output[3],
                         "204:/AACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACC:204",
                         "Test fail: FLT3-ITD -> find variant")
        self.assertEqual(find_output[8],
                         "TACCTTCCCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGA" +
                         "AGTACTCATTATCTGAGGAGCCGGTCACCAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAG" +
                         "TACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTA",
                         "Test fail: FLT3-ITD -> find sequence")

        args = Namespace(
            target=target,
            infile=StringIO(fm_output),
            info="vs_ref",
            min_cov=1,
            exclu="",
            format=None
        )

        with captured_output() as (out, err):
            fr.main_find_report(args, None)

        output = out.getvalue()
        report_output = output.split("\n")
        report_output = report_output[2].split("\t")

        self.assertEqual(report_output[2],
                         "chr13:28034180",
                         "Test fail: FLT3-ITD -> report pos")
        self.assertEqual(report_output[3],
                         "ITD",
                         "Test fail: FLT3-ITD -> report type")
        self.assertEqual(report_output[4],
                         "0",
                         "Test fail: FLT3-ITD -> report Removed")
        self.assertEqual(report_output[5],
                         "75 | 1",
                         "Test fail: FLT3-ITD -> report Added")
        self.assertEqual(report_output[6],
                         "417.6",
                         "Test fail: FLT3-ITD -> report Abnormal")
        self.assertEqual(report_output[7],
                         "1096.7",
                         "Test fail: FLT3-ITD -> report Normal")
        self.assertEqual(report_output[8],
                         "0.276",
                         "Test fail: FLT3-ITD -> report rVAF")
        self.assertEqual(report_output[9],
                         "443",
                         "Test fail: FLT3-ITD -> report Min_coverage")
        self.assertEqual(report_output[11],
                         "/AACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACC",
                         "Test fail: FLT3-ITD -> report variant")
        self.assertEqual(report_output[14],
                         "CTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCGTTCATCACTTTT" +
                         "CCAAAAGCACCTGATCCTAGTACCTTCCCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATA" +
                         "TTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCAACTCCCATTTGAGATCATATTCATATT" +
                         "CTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTA" +
                         "AATTGCTTTTTGTACTTGTGACAAATTAGCAGGGTTAAAACGACAATGAAGAGGAGACAAACACCAATTGTTGCATA" +
                         "GAATGAGATGTTGTCTTGGATGAAAGGGAAGGGGC",
                         "Test fail: FLT3-ITD -> report sequence")


        args = Namespace(
            target=target,
            infile=StringIO(fm_output),
            info="vs_ref",
            min_cov=1,
            exclu="",
            format="vcf"
        )

        with captured_output() as (out, err):
            fr.main_find_report(args, None)

        output = out.getvalue()
        report_output = [o for o in output.split("\n") if o and o[0] != "#"]
        report_output = report_output[0].split("\t")

        self.assertEqual(report_output[1],
                         "28034104",
                         "Test fail: FLT3-ITD -> vcf report pos")
        self.assertEqual(report_output[3],
                         "AAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCT",
                         "Test fail: FLT3-ITD -> vcf report ref")
        self.assertEqual(report_output[4],
                         "AAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCA" +
                         "ACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCT",
                         "Test fail: FLT3-ITD -> vcf report alt")

    def test_FLT3_TKD(self):
        target = "./data/catalog/GRCh38/FLT3-TKD_exon_20.fa"
        args = Namespace(
            count=5,
            graphical=False,
            jellyfish_fn='./data/jf/05H094_FLT3-TKD_del.jf',
            ratio=0.05,
            steps=500,
            branchs=10,
            target_fn=[target],
            verbose=False
        )

        with captured_output() as (out, err):
            fm.main_find_mut(args, None)

        fm_output = out.getvalue()
        find_output = fm_output.split("\n")
        find_output = find_output[11].split("\t")

        self.assertEqual(find_output[2],
                         "Deletion",
                         "Test fail: FLT3-TKD -> find type")
        self.assertEqual(find_output[3],
                         "32:gat/:35",
                         "Test fail: FLT3-TKD -> find variant")
        self.assertEqual(find_output[8],
                         "TGCCCCTGACAACATAGTTGGAATCACTCATATCTCGAGCCAATCCAAAGTCACATATCTT",
                         "Test fail: FLT3-TKD -> find sequence")

        args = Namespace(
            target=target,
            infile=StringIO(fm_output),
            info="vs_ref",
            min_cov=1,
            exclu="",
            format=None
        )

        with captured_output() as (out, err):
            fr.main_find_report(args, None)

        output = out.getvalue()
        report_output = output.split("\n")
        report_output = report_output[2].split("\t")

        self.assertEqual(report_output[2],
                         "",
                         "Test fail: FLT3-TKD -> report pos")
        self.assertEqual(report_output[3],
                         "Deletion",
                         "Test fail: FLT3-TKD -> report type")
        self.assertEqual(report_output[11],
                         "gat/",
                         "Test fail: FLT3-TKD -> report variant")
        self.assertEqual(report_output[14],
                         "TGCCCCTGACAACATAGTTGGAATCACTCATATCTCGAGCCAATCCAAAGTCACATATCTTCACC",
                         "Test fail: FLT3-TKD -> report sequence")

        args = Namespace(
            target=target,
            infile=StringIO(fm_output),
            info="vs_ref",
            min_cov=1,
            exclu="",
            format="vcf"
        )

        with captured_output() as (out, err):
            fr.main_find_report(args, None)

        output = out.getvalue()
        report_output = [o for o in output.split("\n") if o and o[0] != "#"]
        report_output = report_output[0].split("\t")

        self.assertEqual(report_output[1],
                         "28018497",
                         "Test fail: FLT3-TKD -> vcf report pos")
        self.assertEqual(report_output[3],
                         "CATGATA",
                         "Test fail: FLT3-TKD -> vcf report ref")
        self.assertEqual(report_output[4],
                         "CATA",
                         "Test fail: FLT3-TKD -> vcf report alt")

    def test_DNMT3A(self):
        target = "./data/catalog/GRCh38/DNMT3A_R882_exon_23.fa"
        args = Namespace(
            count=5,
            graphical=False,
            jellyfish_fn="./data/jf/02H033_DNMT3A_sub.jf",
            ratio=0.05,
            steps=500,
            branchs=10,
            target_fn=[target],
            verbose=False
        )

        with captured_output() as (out, err):
            fm.main_find_mut(args, None)

        fm_output = out.getvalue()
        find_output = fm_output.split("\n")
        find_output = find_output[11].split("\t")

        self.assertEqual(find_output[2],
                         "Substitution",
                         "Test fail: DNMT3A -> find type")
        self.assertEqual(find_output[3],
                         "33:c/T:34",
                         "Test fail: DNMT3A -> find variant")
        self.assertEqual(find_output[8],
                         "TGACCGGCCCAGCAGTCTCTGCCTCGCCAAGTGGCTCATGTTGGAGACGTCAGTATAGTGGA",
                         "Test fail: DNMT3A -> find sequence")

        args = Namespace(
            target=target,
            infile=StringIO(fm_output),
            info="vs_ref",
            min_cov=1,
            exclu="",
            format=None
        )

        with captured_output() as (out, err):
            fr.main_find_report(args, None)

        output = out.getvalue()
        report_output = output.split("\n")
        report_output = report_output[2].split("\t")

        self.assertEqual(report_output[2],
                         "chr2:25234373",
                         "Test fail: DNMT3A -> report pos")
        self.assertEqual(report_output[3],
                         "Substitution",
                         "Test fail: DNMT3A -> report type")
        self.assertEqual(report_output[4],
                         "1",
                         "Test fail: DNMT3A -> report Removed")
        self.assertEqual(report_output[5],
                         "1",
                         "Test fail: DNMT3A -> report Added")
        self.assertEqual(report_output[6],
                         "33.7",
                         "Test fail: DNMT3A -> report Abnormal")
        self.assertEqual(report_output[7],
                         "48.7",
                         "Test fail: DNMT3A -> report Normal")
        self.assertEqual(report_output[8],
                         "0.409",
                         "Test fail: DNMT3A -> report rVAF")
        self.assertEqual(report_output[9],
                         "32",
                         "Test fail: DNMT3A -> report Min_coverage")
        self.assertEqual(report_output[11],
                         "c/T",
                         "Test fail: DNMT3A -> report variant")
        self.assertEqual(report_output[14],
                         "ATGACCGGCCCAGCAGTCTCTGCCTCGCCAAGTGGCTCATGTTGGAGACGTCAGTATAGTGGACT",
                         "Test fail: DNMT3A -> report sequence")

        args = Namespace(
            target=target,
            infile=StringIO(fm_output),
            info="vs_ref",
            min_cov=1,
            exclu="",
            format="vcf"
        )

        with captured_output() as (out, err):
            fr.main_find_report(args, None)

        output = out.getvalue()
        report_output = [o for o in output.split("\n") if o and o[0] != "#"]
        report_output = report_output[0].split("\t")

        self.assertEqual(report_output[1],
                         "25234373",
                         "Test fail: DNMT3A -> vcf report pos")
        self.assertEqual(report_output[3],
                         "C",
                         "Test fail: DNMT3A -> vcf report ref")
        self.assertEqual(report_output[4],
                         "T",
                         "Test fail: DNMT3A -> vcf report alt")

    def test_not_linear(self):
        ref_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        k_len = 31
        ref_name = "not_linear"

        with self.assertRaises(ValueError):
            uc.get_ref_kmer(ref_seq, k_len, ref_name)

    def test_linear_kmin(self):
        target = "./data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa"
        args = Namespace(
            start=5,
            target_fn=[target]
        )

        with captured_output() as (out, err):
            lk.main_linear_kmin(args, None)

        output = out.getvalue()
        report_output = output.split("\n")
        report_output = report_output[1].split("\t")

        self.assertEqual(report_output[1],
                         "10",
                         "Test fail: linear_kmin -> wrong kmin")

    def test_min_cov(self):
        target = "./data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa"
        args = Namespace(
            jellyfish_fn=["./data/jf/"],
            target_fn=target
        )

        with captured_output() as (out, err):
            mc.main_min_cov(args, None)

        output = out.getvalue()
        report_output = output.split("\n")

        for report_line in report_output:
            report_line = report_line.split("\t")
            if report_line[0] == "./data/jf/02H025_NPM1.jf":
                self.assertEqual(report_line[1],
                                 "0",
                                 "Test fail: min_cov1 -> wronng count")

                self.assertEqual(report_line[2],
                                 "345",
                                 "Test fail: min_cov1 -> wronng length")

                self.assertEqual(report_line[3],
                                 "0",
                                 "Test fail: min_cov1 -> wronng min")

                self.assertEqual(report_line[4],
                                 "0",
                                 "Test fail: min_cov1 -> wronng max")

                self.assertEqual(report_line[5],
                                 "0.00",
                                 "Test fail: min_cov1 -> wronng mean")

                self.assertEqual(report_line[6],
                                 "315",
                                 "Test fail: min_cov1 -> wronng kmer_nb")

                self.assertEqual(report_line[7],
                                 "315",
                                 "Test fail: min_cov1 -> wronng kmer_nb_0")

            if report_line[0] == "./data/jf/03H112_IandI.jf":
                self.assertEqual(report_line[1],
                                 "275596",
                                 "Test fail: min_cov2 -> wronng count")

                self.assertEqual(report_line[2],
                                 "345",
                                 "Test fail: min_cov2 -> wronng length")

                self.assertEqual(report_line[3],
                                 "618",
                                 "Test fail: min_cov2 -> wronng min")

                self.assertEqual(report_line[4],
                                 "1368",
                                 "Test fail: min_cov2 -> wronng max")

                self.assertEqual(report_line[5],
                                 "874.91",
                                 "Test fail: min_cov2 -> wronng mean")

                self.assertEqual(report_line[6],
                                 "315",
                                 "Test fail: min_cov2 -> wronng kmer_nb")

                self.assertEqual(report_line[7],
                                 "0",
                                 "Test fail: min_cov2 -> wronng kmer_nb_0")


def runTests():
    unittest.main()


if __name__ == "__main__":
    runTests()
