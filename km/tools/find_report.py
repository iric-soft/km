import sys
import re

from .. utils import common as uc


def print_line(sample, region, location, type_var, removed,
               added, abnormal, normal, ratio, min_cov, min_ref,
               variant, target, info, var_seq, ref_seq):
    if min_ref is None:
        line = "\t".join([sample, region, location, type_var, removed,
                          added, abnormal, normal, ratio, min_cov, "",
                          variant, target, info, var_seq, ref_seq])
    else:
        line = "\t".join([sample, region, location, type_var, removed,
                          added, abnormal, normal, ratio, min_cov, min_ref,
                          variant, target, info, var_seq, ref_seq])

    sys.stdout.write(line + "\n")


def init_ref_seq(arg_ref):
    # BE 1-based
    nts = [-1]
    ref_seq = []
    chro = None
    if arg_ref:
        # CODE for multiple >LOC lines:
        for line in open(arg_ref, "r"):
            line = line.strip()
            if line[0] == '>':
                exon = line.strip(">")
                chro, pos = exon.split(":")
                refstart, refstop = pos.split("-")
                for i in xrange(int(refstart), int(refstop) + 1):
                    nts += [i]
        ref_seq = "X" + ''.join(ref_seq)

    return(nts, ref_seq, chro)


def create_report(args):
    variants = {}
    samples = {}
    data = {}

    (nts, ref_seq, chro) = init_ref_seq(args.target)

    print_line("Sample", "Region", "Location", "Type", "Removed",
               "Added", "Abnormal", "Normal", "Ratio", "Min_coverage",
               "Min_ref_cov", "Variant", "Target", "Info", "Variant_sequence",
               "Reference_sequence")

    for line in args.infile:
        # filter header
        if line[0] == "#":
            # sys.stderr.write("Filtred: " + line)
            continue

        # filter on info column
        if not re.search(args.info, line):
            # sys.stderr.write("Filtred: " + line)
            continue

        tok = line.strip("\n").split("\t")
        if len(tok) > 1 and tok[0] != "Database":
            pathvals = {
                'project': '',
                'sample': tok[0]
            }
            samp = (pathvals['project'], pathvals['sample'])

            variant = (tok[2], tok[3])
            if variant not in variants:
                variants[variant] = 0
            variants[variant] += 1

            if samp not in samples:
                samples[samp] = set()
            samples[samp].add(variant)

            if samp not in data:
                data[samp] = {}

            data[samp][variant] = tok

            query = tok[1]
            ratio = tok[4]
            alt_exp = tok[5]
            ref_exp = tok[10]
            min_cov = tok[6]
            start_off = tok[7]
            alt_seq = tok[8]
            refSeq = tok[11]
            min_ref = None

            if args.ref != "":
                res = uc.get_cov(args.ref, alt_seq)
                min_ref = str(res[2])

            if int(min_cov) < args.min_cov:
                continue

            if variant[0] == 'Reference':
                print_line(samp[1], samp[0], '-', variant[0], '0', '0',
                           '0.0', alt_exp, tok[4], min_cov, min_ref, '-',
                           query, tok[-1], "", "")
                continue

            start, mod, stop = variant[1].split(":")
            delet, insert = mod.split("/")

            added = str(len(insert))

            # Reinterpret mutation for small ITD.
            # INSERTIONS
            if len(delet) == 0 and len(insert) != 0:
                pos = int(start)-1
                pos -= int(start_off)
                upstream = alt_seq[pos-len(insert):pos]
                match = 0
                # careful, going upstream may put us outside the reference.
                if pos-len(insert) >= 0:
                    for i in xrange(0, len(insert)):
                        if insert[i] == upstream[i]:
                            match += 1
                    match = float(match)/len(insert)

                insert_type = "Insertion"
                region = chro + ":" + str(nts[int(pos-1)]) + "-" + str(nts[int(stop)-1])
                if pos-len(insert) >= 0 and len(insert) >= 3 and insert == upstream:
                    insert_type = "ITD"
                    region = chro + ":" + str(nts[int(pos-len(insert)+1)]) + "-" + str(nts[int(stop)-1])
                    added += " | " + str(nts[int(stop)-1] - nts[int(pos-len(insert)+1)] + 1)
                elif pos-len(insert) >= 0 and len(insert) >= 3 and match > 0.5:
                    insert_type = "I&I"
                    region = chro + ":" + str(nts[int(pos-len(insert)+1)]) + "-" + str(nts[int(stop)-1])
                    added += " | " + str(nts[int(stop)-1] - nts[int(pos-len(insert)+1)] + 1)

                location = chro + ":" + str(nts[int(stop)])

            elif variant[0] == 'Deletion':
                region = chro + ":" + str(nts[int(start)-1]+1) + "-" + str(nts[int(stop)-1])
                location = ""
                insert_type = variant[0]
            elif variant[0] == 'Substitution':
                region = chro + ":" + str(nts[int(start)-1]+1) + "-" + str(nts[int(stop)-1])
                location = chro + ":" + str(nts[int(stop)-1])
                insert_type = variant[0]
            elif variant[0] == 'Indel':
                region = chro + ":" + str(nts[int(start)-1]+1) + "-" + str(nts[int(stop)-1])
                location = chro + ":" + str(nts[int(stop)])
                insert_type = variant[0]
            else:
                sys.stderr.write("WARNING: This variant isn't take account\n")
                sys.stderr.write(" - variant: " + str(variant[0]) + "\n")
                sys.stderr.write(" - line: " + line)
                sys.exit()

            print_line(samp[1], region, location, insert_type,
                       str(len(delet)), added, alt_exp, ref_exp, ratio,
                       min_cov, min_ref, mod, query, tok[-1],
                       alt_seq, refSeq)


def main_find_report(args, argparser):

    if args.infile.isatty() or args.target is None:
        argparser.print_help()
        sys.exit()

    create_report(args)
