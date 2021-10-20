#                             -*- Mode: Python -*-
# Report.py
#

import sys
import re
from collections import namedtuple

from . import common as uc


Exon = namedtuple('Exon', ['chro', 'strand', 'nts'])


class Line:
    def __init__(self, tok, ref_exons, exclu):
        self.samp = tok[0]
        self.query = tok[1]
        self.ratio = tok[4]
        self.alt_exp = tok[5]
        self.ref_exp = tok[6]
        self.min_cov = tok[7]
        self.ref_min_cov = tok[8]
        self.alt_seq = tok[9]
        self.ref_seq = tok[10]
        self.info = tok[11]

        self.variant = (tok[2], tok[3])

        if 'offset' in self.info:
            self.start_off = self.info.split('offset=')[1]
        else:
            self.start_off = '0'

        if exclu != "" and alt_seq != "":
            res = uc.get_cov(exclu, alt_seq)
            self.min_exclu = str(res[2])
        else:
            self.min_exclu = ""

        ref_pair = self.query.split('::')
        assert len(ref_pair) <= 2

        chro = set([ref_exons[r].chro for r in ref_pair])
        assert len(chro) == 1  # TODO: Implement > 1 chromosome
        self.chro = chro.pop()

        strand = set([ref_exons[r].strand for r in ref_pair])
        assert len(strand) == 1  # TODO: Implement mixed strand
        self.strand = strand.pop()

        self.nts = [x for r in ref_pair for x in ref_exons[r].nts]

        # case: entries with no mutations
        if self.variant[0] == 'Reference':
            self.isref = True

            first, last = self.nts[0], self.nts[-1]
            if self.strand == "-":
                first, last = last, first

            self.region = "{}:{}-{}".format(self.chro, first, last)

            self.location = '-'
            self.insert_type = self.variant[0]
            self.removed = '0'
            self.added = '0'
            self.ref_exp = self.alt_exp
            self.alt_exp = '0.0'
            self.mod = '-'
            self.alt_seq = ''
            self.ref_seq = ''
        else:

            self.isref = False

            self.region = None
            self.location = None
            self.insert_type = None
            self.removed = None
            self.added = None
            self.mod = None


    def process_variant(self):
        # case: there is a mutation

        start, self.mod, stop = self.variant[1].split(":")
        delet, insert = self.mod.split("/")

        self.added = str(len(insert))
        self.removed = str(len(delet))

        self.insert_type = self.variant[0]

        # start and end positions in 0-based coordinates
        pos = int(start) - 1
        pos -= int(self.start_off)
        end = int(stop) - 2  # one to go back to last position, the other for 0-base
        end -= int(self.start_off)

        start_pos, end_pos = self.nts[pos], self.nts[end]
        if self.strand == "-":
            start_pos, end_pos = end_pos, start_pos
        self.region = "{}:{}-{}".format(self.chro, start_pos, end_pos + 1)

        if len(delet) == 0 and len(insert) != 0:
            start_pos, end_pos = self.nts[pos], self.nts[end + 1]  # insertions end at last position
            if self.strand == "-":
                start_pos, end_pos = end_pos, start_pos
            self.region = "{}:{}-{}".format(self.chro, start_pos, end_pos + 1)

            self.location = self.chro + ":" + str(end_pos)

            if self.insert_type == "ITD":
                self.added += " | " + str(end_pos - start_pos + 1)
            elif self.insert_type == "I&I":
                self.added += " | " + str(end_pos - start_pos + 1)

        elif self.variant[0] == 'Deletion':
            self.location = ""

        elif self.variant[0] == 'Substitution':
            self.location = self.chro + ":" + str(start_pos)

        elif self.variant[0] == 'Indel':
            self.location = self.chro + ":" + str(end_pos)

        else:
            sys.stderr.write("WARNING: This variant isn't taken account\n")
            sys.stderr.write(" - variant: " + str(variant[0]) + "\n")
            sys.stderr.write(" - line: " + line)
            sys.exit()

        self.delet = delet
        self.insert = insert
        self.pos = pos

    def __str__(self):
        return "\t".join([
                self.samp,
                self.region,
                self.location,
                self.insert_type,
                self.removed,
                self.added,
                self.alt_exp,
                self.ref_exp,
                self.ratio,
                self.min_cov,
                self.ref_min_cov,
                self.min_exclu,
                self.mod,
                self.query,
                self.info,
                self.alt_seq,
                self.ref_seq
            ])


class Report:
    header = [
            "Sample",
            "Region",
            "Location",
            "Type",
            "Removed",
            "Added",
            "Abnormal",
            "Normal",
            "rVAF",
            "Min_coverage",
            "Ref_min_coverage",
            "Exclu_min_cov",
            "Variant",
            "Target",
            "Info",
            "Variant_sequence",
            "Reference_sequence"
        ]

    def __init__(self, infile, info, min_cov, exclu):
        self.infile = infile
        self.info = info
        self.min_cov = min_cov
        self.exclu = exclu

    def print_header(self):
        sys.stdout.write("\t".join(self.header) + "\n")

    def parse(self):
        for line in self._parse():
            if not line.isref:
                line.process_variant()
            sys.stdout.write(str(line) + "\n")

    def _parse(self):
        ready = False
        parsed_exons = {}
        ref_exons = {}

        for line in self.infile:
            # filter header
            if line[0] == "#":
                if line.startswith('#target:'):
                    exon = line.strip().split(':')[1].replace('/', ':')
                    exon_id, exon_loc = exon.split('=')
                    parsed_exons[exon_id] = self.parse_exon(exon_loc)
                    ready = False
            else:
                if not ready:
                    ref_exons = parsed_exons
                    parsed_exons = {}
                    ready = True
                tok = line.strip("\n").split("\t")
                # filter on info column
                if re.search(self.info, line) and tok[0] != "Database" and len(tok) > 1:
                    line_obj = Line(tok, ref_exons, self.exclu)
                    if int(line_obj.min_cov) >= self.min_cov:
                        yield line_obj

    @staticmethod
    def parse_exon(exon):
        chro, loc, strand = exon.split(':')
        refstart, refstop = loc.split('-')
        refstart = int(refstart)
        refstop = int(refstop)
        nts = list(range(refstart, refstop + 1)) # BE 1-based
        if strand == "-":
            nts = nts[::-1]
        nts = tuple(nts)

        return Exon(chro, strand, nts)


class TableReport(Report):
    def print_header(self):
        return

    def parse(self):
        variants = {}
        samples = {}
        data = {}

        for line in self._parse():
            if not line.isref:
                line.process_variant()

            if "/" not in line.variant[0]:
                var_name = line.variant[0] + "/" + line.query
            else:
                var_name = line.variant[0]
            region_mod = line.region + ":" + line.mod if line.mod != "-" else line.region
            var = (var_name, region_mod)
            if var not in variants:
                variants[var] = 0
            variants[var] += 1

            samp = line.samp

            if samp not in samples:
                samples[samp] = set()
            samples[samp].add(var)

            if samp not in data:
                data[samp] = {}
            data[samp][var] = float(line.ratio)

        sorted_variants = sorted(variants, key=variants.get, reverse=True)

        sys.stdout.write("Sample")
        for v in sorted_variants:
            if v[0].split("/")[0] == "Reference":
                sys.stdout.write("\t" + v[0])
            else:
                sys.stdout.write("\t" + v[1])
        sys.stdout.write("\n")

        for s, sv in samples.items():
            sys.stdout.write(s)
            for v in sorted_variants:
                if v in sv:
                    if 'Reference' not in v[0] and (not data[s][v]):
                        sys.stdout.write("\t" + ".")
                    else:
                        sys.stdout.write("\t" + str(data[s][v]))
                else:
                    sys.stdout.write("\t" + ".")
            sys.stdout.write("\n")


class VCFReport(Report):
    header  = '##fileformat=VCFv4.1\n'
    header += '##INFO=<ID=TYPE,Number=A,Type=String,Description='
    header += '"The type of variant: Insertion, ITD, I&I, Deletion, Substitution or Indel.">\n'
    header += '##INFO=<ID=TARGET,Number=A,Type=String,Description='
    header += '"Name of the sequencing that contains the mutation.">\n'
    header += '##INFO=<ID=RATIO,Number=A,Type=String,Description="Mutation to reference ratio.">\n'
    header += '##INFO=<ID=MINCOV,Number=A,Type=String,Description='
    header += '"Minimum k-mer coverage of alternative allele.">\n'
    header += '##INFO=<ID=REMOVED,Number=A,Type=String,Description="Number of removed bases.">\n'
    header += '##INFO=<ID=ADDED,Number=A,Type=String,Description="Number of added bases.">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def print_header(self):
        sys.stdout.write(self.header)

    def parse(self):
        for line in self._parse():
            if line.isref:
                continue

            line.process_variant()

            res = self.canonicalize(line.delet, line.insert, line.pos, line.ref_seq)
            ref_var, alt_var, ibef, iaft = res

            if line.strand == "+":
                loc_var = line.nts[ibef]
                end_var = line.nts[ibef+len(ref_var)-1]
            elif line.strand == "-":
                loc_var = line.nts[iaft]
                end_var = line.nts[iaft-len(ref_var)+1]

            if loc_var + len(ref_var) - 1 != end_var:
                sys.stderr.write("NOTE: Mutation overlaps 2 exons or more, VCF output is disabled \n")
            else:
                if line.strand == '-':
                    complement = str.maketrans('ATGCU', 'TACGA')
                    ref_var = ref_var.translate(complement)[::-1]
                    alt_var = alt_var.translate(complement)[::-1]

                flags = "TYPE=%s;TARGET=%s;RATIO=%s;MINCOV=%s;REMOVED=%s;ADDED=%s" % (
                        line.insert_type,
                        line.query,
                        line.ratio,
                        line.min_cov,
                        line.removed,
                        line.added
                    )

                out_line = "\t".join([
                        line.chro,
                        str(loc_var),
                        ".",
                        ref_var,
                        alt_var,
                        ".",
                        ".",
                        flags
                    ])

                sys.stdout.write(out_line + "\n")


    @staticmethod
    def canonicalize(delet, insert, pos, ref_seq):
        sys.setrecursionlimit(10000)

        def get_extremities(va, p, rs):
            if p - 1 > 0 and rs[p - 1] == va[-1]:
                return get_extremities(rs[p - 1] + va[:-1], p - 1, rs)
            return p - 1

        delet = delet.upper()
        insert = insert.upper()
        ref_seq = ref_seq.upper()
        reflen = len(ref_seq)

        #for a, b in list(zip(delet, insert)):
        #    if a == b:
        #        delet = delet[1:]
        #        insert = insert[1:]
        #        continue
        #    break
        #for a, b in list(zip(delet[::-1], insert[::-1])):
        #    if a == b:
        #        delet = delet[:-1]
        #        insert = insert[:-1]
        #        continue
        #    break

        end = pos + len(delet)

        if delet and insert:
            if len(delet) == len(insert):
                ibef = pos
                iaft = end
                ref_var = delet
                alt_var = insert
            else:
                ibef = pos - 1
                iaft = end + 1
                ref_var = ref_seq[pos-1] + delet + ref_seq[end + 1]
                alt_var = ref_seq[pos-1] + insert + ref_seq[end + 1]
        elif not delet:
            ibef = get_extremities(insert, pos, ref_seq)  # include current position
            iaft = reflen - get_extremities(insert[::-1], reflen-pos, ref_seq[::-1]) - 1
            b = ref_seq[ibef:pos]  # before
            a = ref_seq[pos:iaft+1]  # after
            ref_var = b + a
            alt_var = b + insert + a
        elif not insert:
            ibef = get_extremities(delet, pos, ref_seq)
            iaft = reflen - get_extremities(delet[::-1], reflen-pos-len(delet), ref_seq[::-1]) - 1
            b = ref_seq[ibef:pos]
            a = ref_seq[pos+len(delet):iaft+1]
            ref_var = b + delet + a
            alt_var = b + a
        else:
            assert delet or insert  # will fail

        return ref_var, alt_var, ibef, iaft
