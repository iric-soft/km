import sys
import re

from .. utils import common as uc


def print_line(sample, region, location, type_var, removed,
               added, abnormal, normal, ratio, min_cov, min_exclu,
               variant, target, info, var_seq, ref_seq):
    line = "\t".join([sample, region, location, type_var, removed,
                      added, abnormal, normal, ratio, min_cov, min_exclu,
                      variant, target, info, var_seq, ref_seq])
    sys.stdout.write(line + "\n")


def print_vcf_header():
    header  = '##fileformat=VCFv4.1\n'
    header += '##INFO=<ID=TYPE,Number=A,Type=String,Description='
    header += '"The type of variant, either Insertion, ITD, I&I, Deletion, Substitution or Indel.">\n'
    header += '##INFO=<ID=TARGET,Number=A,Type=String,Description='
    header += '"Name of the sequencing that contains the mutation.">\n'
    header += '##INFO=<ID=RATIO,Number=A,Type=String,Description="Ratio of mutation to reference.">\n'
    header += '##INFO=<ID=MINCOV,Number=A,Type=String,Description='
    header += '"Minimum k-mer coverage of alternative allele.">\n'
    header += '##INFO=<ID=REMOVED,Number=A,Type=String,Description="Number of removed bases.">\n'
    header += '##INFO=<ID=ADDED,Number=A,Type=String,Description="Number of added bases.">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    sys.stdout.write(header)


def print_vcf_line(chro, loc, ref_var, alt_var, type_var, target, ratio, min_cov, rem, ad):
    line = "\t".join([chro, str(loc), ".", ref_var, alt_var, ".", ".",
                      "TYPE="+type_var+";TARGET="+target+";RATIO="+ratio+";MINCOV="+min_cov +
                      ";REMOVED="+str(rem)+";ADDED="+str(ad)])
    sys.stdout.write(line + "\n")


def init_ref_seq(arg_ref):
    if not arg_ref:
        sys.exit("ERROR: Target file is empty\n")

    # BE 1-based
    nts = []
    chro = None
    strand = None

    # CODE for multiple >LOC lines:
    for line in open(arg_ref, "r"):
        line = line.strip()
        nt = []
        # Parse attributes
        if line[0] == '>':
            # sanity check
            loc = line.split(" ")[0]
            if ":" not in loc or "-" not in loc:
                sys.exit('ERROR: Fasta entries do not contain a correctly ' +
                         'formatted location: {}\n'.format(loc))

            # look up attributes in fasta file
            line = line.replace(">", "location=", 1)
            attr = {x.split("=")[0].strip(): x.split("=")[1].strip() for x in line.split("|")}
            exon = attr["location"]
            chro, pos = exon.split(":")
            refstart, refstop = pos.split("-")

            # get nt coordinates on the genome
            if 'strand' not in list(attr.keys()):
                attr['strand'] = '+'
                sys.stderr.write("WARNING: Strand is assumed to be '+' \n")
            strand = attr["strand"]
            for i in range(int(refstart), int(refstop) + 1):
                nt += [i]
            nt = nt[::-1] if strand == "-" else nt
            nts.extend(nt)

    return nts, chro, strand


def create_report(args):

    # Find correct extremities of a mutation
    sys.setrecursionlimit(10000)

    def get_extremities(va, p, rs):
        if p - 1 > 0 and rs[p - 1] == va[-1]:
            return get_extremities(rs[p - 1] + va[:-1], p - 1, rs)
        return p - 1

    if args.format == "vcf" and args.info == "cluster":
        # Note: could salvage that option if we get the fill ref from vs_ref entries
        sys.exit("ERROR: -f vcf and -i cluster options are incompatible")

    variants = {}
    samples = {}
    data = {}
    vcf = True if args.format == 'vcf' else False
    table = True if args.format == 'table' else False

    (nts, chro, strand) = init_ref_seq(args.target)

    if vcf:
        print_vcf_header()
    elif not table:
        print_line("Sample", "Region", "Location", "Type", "Removed",
                   "Added", "Abnormal", "Normal", "rVAF", "Min_coverage",
                   "Exclu_min_cov", "Variant", "Target", "Info", "Variant_sequence",
                   "Reference_sequence")

    for line in args.infile:
        # filter header
        if line[0] == "#":
            # sys.stderr.write("Filtred: " + line)
            continue

        tok = line.strip("\n").split("\t")

        # filter on info column
        if not re.search(args.info, line) or tok[0] == "Database" or len(tok) <= 1:
            # sys.stderr.write("Filtered: " + line)
            continue

        samp = tok[0]
        query = tok[1]
        ratio = tok[4]
        alt_exp = tok[5]
        ref_exp = tok[9]
        min_cov = tok[6]
        start_off = tok[7]
        alt_seq = tok[8]
        refSeq = tok[10]
        info = tok[11]

        min_exclu = ""
        variant = (tok[2], tok[3])
        ref_seq = refSeq.upper()

        if args.exclu != "" and alt_seq != "":
            res = uc.get_cov(args.exclu, alt_seq)
            min_exclu = str(res[2])

        if int(min_cov) < args.min_cov:
            continue

        # case: entries with no mutations
        if variant[0] == 'Reference':
            mod = ""
            if strand == "-":
                region = "{}:{}-{}".format(chro, nts[-1], nts[0])
            else:
                region = "{}:{}-{}".format(chro, nts[0], nts[-1])
            if not vcf and not table:
                print_line(samp, region, '-', variant[0], '0', '0',
                           '0.0', alt_exp, tok[4], min_cov, min_exclu, '-',
                           query, tok[-1], "", "")
                continue
            elif vcf:
                continue

        # case: there is a mutation
        else:
            start, mod, stop = variant[1].split(":")
            delet, insert = mod.split("/")

            added = str(len(insert))
            removed = str(len(delet))

            # start and end positions in 0-based coordinates
            pos = int(start) - 1
            pos -= int(start_off)
            end = int(stop) - 2  # one to go back to last position, the other for 0-base
            end -= int(start_off)

            if strand == "+":
                start_pos = nts[pos]
                end_pos = nts[end]
            elif strand == "-":
                start_pos = nts[end]
                end_pos = nts[pos]

            region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)

            ref_var = delet.upper()
            alt_var = insert.upper()
            loc_var = start_pos
            end_var = end_pos

            if len(delet) == 0 and len(insert) != 0:
                if strand == "+":
                    start_pos = nts[pos]
                    end_pos = nts[end + 1]  # insertions end at last position
                elif strand == "-":
                    start_pos = nts[end + 1]
                    end_pos = nts[pos]
                region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)

                var = insert.upper()
                ibef = get_extremities(var, pos, ref_seq)  # include current position
                before = ref_seq[ibef:pos]
                iaft = get_extremities(var[::-1], len(ref_seq)-pos, ref_seq[::-1])
                after = ref_seq[::-1][iaft:len(ref_seq)-pos][::-1]
                iaft = len(ref_seq) - iaft - 1
                ref_var = before + after
                alt_var = before + var + after
                loc_var = nts[iaft] if strand == "-" else nts[ibef]
                end_var = nts[iaft-len(ref_var)+1] if strand == "-" else nts[ibef+len(ref_var)-1]

                if loc_var + len(ref_var) - 1 != end_var and vcf:
                    sys.stderr.write("NOTE: Mutation overlaps 2 exons or more, VCF output is disabled \n")
                    continue

                # Reinterpret mutations for small ITDs
                # careful, going upstream may put us outside the reference.
                upstream = alt_seq[pos-len(insert):pos]
                match = 0
                if pos-len(insert) >= 0:
                    for i in range(0, len(insert)):
                        if insert[i] == upstream[i]:
                            match += 1
                    match = float(match)/len(insert)

                insert_type = "Insertion"
                if pos-len(insert) >= 0 and len(insert) >= 3 and insert == upstream:
                    insert_type = "ITD"
                    added += " | " + str(end_pos - start_pos + 1)
                elif pos-len(insert) >= 0 and len(insert) >= 3 and match > 0.5:
                    insert_type = "I&I"
                    added += " | " + str(end_pos - start_pos + 1)

                location = chro + ":" + str(end_pos)

            elif variant[0] == 'Deletion':
                region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)
                location = ""
                insert_type = variant[0]

                var = delet.upper()
                ibef = get_extremities(var, pos, ref_seq)
                before = ref_seq[ibef:pos]
                iaft = get_extremities(var[::-1], len(ref_seq)-pos-1-len(var)+1, ref_seq[::-1])
                after = ref_seq[::-1][iaft:len(ref_seq)-pos-1-len(var)+1][::-1]
                iaft = len(ref_seq) - iaft - 1
                ref_var = before + var + after
                alt_var = before + after
                loc_var = nts[iaft] if strand == "-" else nts[ibef]
                end_var = nts[iaft-len(ref_var)+1] if strand == "-" else nts[ibef+len(ref_var)-1]

                if loc_var + len(ref_var) - 1 != end_var and vcf:
                    continue

            elif variant[0] == 'Substitution':
                location = chro + ":" + str(start_pos)
                insert_type = variant[0]

                if loc_var + len(ref_var) - 1 != end_var and vcf:
                    sys.stderr.write("NOTE: Mutation overlaps 2 exons or more, VCF output is disabled \n")
                    continue

            elif variant[0] == 'Indel':
                location = chro + ":" + str(end_pos)
                insert_type = variant[0]

                ref_var = ref_seq[pos-1] + delet.upper() + ref_seq[end + 1]
                alt_var = ref_seq[pos-1] + insert.upper() + ref_seq[end + 1]
                loc_var = start_pos - 1
                end_var = end_pos + 1

                if loc_var + len(ref_var) - 1 != end_var and vcf:
                    sys.stderr.write("NOTE: Mutation overlaps 2 exons or more, VCF output is disabled \n")
                    continue

            else:
                sys.stderr.write("WARNING: This variant isn't taken account\n")
                sys.stderr.write(" - variant: " + str(variant[0]) + "\n")
                sys.stderr.write(" - line: " + line)
                sys.exit()

        if not vcf and not table:
            print_line(samp, region, location, insert_type,
                       removed, added, alt_exp, ref_exp, ratio,
                       min_cov, min_exclu, mod, query, info,
                       alt_seq, refSeq)

        elif vcf:
            complement = str.maketrans('ATGCU', 'TACGA')
            ref_var = ref_var.translate(complement)[::-1] if strand == '-' else ref_var
            alt_var = alt_var.translate(complement)[::-1] if strand == '-' else alt_var
            print_vcf_line(chro, loc_var, ref_var, alt_var, insert_type,
                           query, ratio, min_cov, removed, added.replace(" ", ""))

        elif table:
            var_name = variant[0] + "/" + query if "/" not in variant[0] else variant[0]
            region_mod = region + ":" + mod if mod else region
            var = (var_name, region_mod)
            if var not in variants:
                variants[var] = 0
            variants[var] += 1

            if samp not in samples:
                samples[samp] = set()
            samples[samp].add(var)

            if samp not in data:
                data[samp] = {}
            data[samp][var] = float(ratio)

    if table:
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


def main_find_report(args, argparser):

    if args.infile.isatty() or args.target is None:
        argparser.print_help()
        sys.exit()

    create_report(args)
