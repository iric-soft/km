import sys
import re


def main_find_report(args, argparser):

    if args.infile.isatty() or args.target is None:
        argparser.print_help()
        sys.exit()

    arg_ref = args.target
    infile = args.infile

    variants = {}
    samples = {}
    data = {}

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

    print "\t".join(['Sample', 'Region', 'Location', 'Type', 'Removed',
                     'Added', 'Abnormal', 'Normal', 'Ratio', 'Min coverage',
                     'Variant', 'Target', 'Info', 'Variant sequence',
                     'Reference sequence'])

    for line in infile:
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
            alt_ratio = tok[5]
            ref_ratio = tok[9]
            min_cov = tok[6]
            alt_seq = tok[7]
            refSeq = tok[10]

            if int(min_cov) <= args.min_cov:
                continue

            if variant[0] == 'Reference':
                print "\t".join([samp[1], samp[0], '-', variant[0], '0', '0',
                                 '0.0', alt_ratio, tok[4], min_cov, '-', query,
                                 tok[-1], "", ""])
                continue

            start, mod, stop = variant[1].split(":")
            delet, insert = mod.split("/")
            # sys.stderr.write("start: " + str(start) + "\n")
            # sys.stderr.write("stop: " + str(stop) + "\n")
            # sys.stderr.write("delet: " + str(delet) + "\n")
            # sys.stderr.write("insert: " + str(insert) + "\n")

            if int(start) == len(nts) or int(stop) == len(nts):
                sys.stderr.write("WARNING: Mutation point outside reference range")
                sys.stderr.write("WARNING: " + "\t".join([samp[1], samp[0],
                                                          chro + ":" + str(nts[int(stop)-1] + 1),
                                                          variant[0],
                                                          str(len(delet)),
                                                          str(len(insert)),
                                                          alt_ratio,
                                                          ref_ratio, tok[4],
                                                          mod, query]))
                print "\t".join([samp[1], samp[0],
                                 chro + ":" + str(nts[int(stop)-1]+1),
                                 variant[0], str(len(delet)), str(len(insert)),
                                 alt_ratio, ref_ratio, ratio, min_cov, mod,
                                 query, tok[-1], tok[7], refSeq])
            else:
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
                    added = str(len(insert))
                    if pos-len(insert) >= 0 and len(insert) >= 3 and insert == upstream:
                        insert_type = "ITD"
                        region = chro + ":" + str(nts[int(pos-len(insert)+1)]) + "-" + str(nts[int(stop)-1])
                        added += " | " + str(nts[int(stop)-1] - nts[int(pos-len(insert)+1)] + 1)
                    elif pos-len(insert) >= 0 and len(insert) >= 3 and match > 0.5:
                        insert_type = "I&I"
                        region = chro + ":" + str(nts[int(pos-len(insert)+1)]) + "-" + str(nts[int(stop)-1])
                        added += " | " + str(nts[int(stop)-1] - nts[int(pos-len(insert)+1)] + 1)

                    print "\t".join([samp[1],
                                     region,
                                     chro + ":" + str(nts[int(stop)]),
                                     insert_type, str(len(delet)),
                                     added,
                                     alt_ratio, ref_ratio, ratio, min_cov, mod,
                                     query, tok[-1], alt_seq, refSeq])

                elif variant[0] == 'Deletion':
                    print "\t".join([samp[1],
                                     chro + ":" + str(nts[int(start)-1]+1) + "-" + str(nts[int(stop)-1]),
                                     "", variant[0], str(len(delet)),
                                     str(len(insert)), alt_ratio, ref_ratio,
                                     ratio, min_cov, mod, query, tok[-1],
                                     alt_seq, refSeq])
                # SNP
                elif variant[0] == 'Substitution':
                    print "\t".join([samp[1],
                                     chro + ":" + str(nts[int(start)-1]+1) + "-" + str(nts[int(stop)-1]),
                                     chro + ":" + str(nts[int(stop)-1]),
                                     variant[0], str(len(delet)),
                                     str(len(insert)), alt_ratio, ref_ratio,
                                     ratio, min_cov, mod, query, tok[-1],
                                     alt_seq, refSeq])
                else:
                    print "\t".join([samp[1],
                                     chro + ":" + str(nts[int(start)-1]+1) + "-" + str(nts[int(stop)-1]),
                                     chro + ":" + str(nts[int(stop)]),
                                     variant[0], str(len(delet)),
                                     str(len(insert)), alt_ratio, ref_ratio,
                                     ratio, min_cov, mod, query, tok[-1],
                                     alt_seq, refSeq])
