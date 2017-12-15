import sys
import re


def print_line(sample, region, location, type_var, removed,
               added, abnormal, normal, ratio, min_cov,
               variant, target, info, var_seq, ref_seq):
    line = "\t".join([sample, region, location, type_var, removed,
                      added, abnormal, normal, ratio, min_cov,
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


def create_report(arg_ref, infile, args_info, args_min_cov):
    variants = {}
    samples = {}
    data = {}

    (nts, ref_seq, chro) = init_ref_seq(arg_ref)

    print_line("Sample", "Region", "Location", "Type", "Removed",
               "Added", "Abnormal", "Normal", "Ratio", "Min coverage",
               "Variant", "Target", "Info", "Variant sequence",
               "Reference sequence")

    for line in infile:
        # filter header
        if line[0] == "#":
            # sys.stderr.write("Filtred: " + line)
            continue

        # filter on info column
        if not re.search(args_info, line):
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
            alt_ratio = tok[5]
            ref_ratio = tok[9]
            min_cov = tok[6]
            alt_seq = tok[7]
            refSeq = tok[10]

            if int(min_cov) <= args_min_cov:
                continue

            if variant[0] == 'Reference':
                print_line(samp[1], samp[0], '-', variant[0], '0', '0',
                            '0.0', alt_ratio, tok[4], min_cov, '-', query,
                            tok[-1], "", "")
                continue

            start, mod, stop = variant[1].split(":")
            delet, insert = mod.split("/")

            added = str(len(insert))

            # Reinterpret mutation for small ITD.
            # INSERTIONS
            if len(delet) == 0 and len(insert) != 0:
                pos = int(start)-1
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
            # SNP
            elif variant[0] == 'Substitution':
                region = chro + ":" + str(nts[int(start)-1]+1) + "-" + str(nts[int(stop)-1])
                location = chro + ":" + str(nts[int(stop)-1])
                insert_type = variant[0]
            else:
                sys.stderr.write("WARNING: This case isn't take account")
                sys.exit()

            print_line(samp[1], region, location, insert_type,
                       str(len(delet)), added, alt_ratio,
                       ref_ratio, ratio, min_cov, mod, query, tok[-1],
                       alt_seq, refSeq)


def main_find_report(args, argparser):

    if args.infile.isatty() or args.target is None:
        argparser.print_help()
        sys.exit()

    create_report(args.target, args.infile, args.info, args.min_cov)
