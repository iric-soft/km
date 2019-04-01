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
    header='##fileformat=VCFv4.1\n'
    header+='##INFO=<ID=TYPE,Number=A,Type=String,Description=\
"The type of variant, either Insertion, ITD, I&I, Deletion, Substitution or Indel.">\n'
    header+='##INFO=<ID=TARGET,Number=A,Type=String,Description=\
"Name of the target sequencing containing mutation.">\n'
    header+='##INFO=<ID=RATIO,Number=A,Type=String,Description="Ration of mutation to reference.">\n'
    header+='##INFO=<ID=MINCOV,Number=A,Type=String,Description="Minimum k-mer coverage on \
mutated sequence.">\n'
    header+='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'

    sys.stdout.write(header)


def print_vcf_line(region, ref_var, alt_var, type_var, target, ratio, min_cov):
    chro, location = region.split(":")
    location = str(int(location.split("-")[0])-1)
    line = "\t".join([chro, location, ".", ref_var, alt_var, ".", ".",
                      "TYPE="+type_var+";TARGET="+target+";RATIO="+ratio+";MINCOV="+min_cov])
    sys.stdout.write(line + "\n")


def init_ref_seq(arg_ref):
    # BE 1-based
    chro = None
    strand = None
    ref_attributes = {}
    all_nts = []
    all_ref_seq = []
    exons = True
    if arg_ref:
        # CODE for multiple >LOC lines:
        for line in open(arg_ref, "r"):
            line = line.strip()
            if line[0] == '>':
                nts = []
                ref_seq = []
                line = line.replace(">", "location=", 1)
                attr = {x.split("=")[0].strip():x.split("=")[1].strip() for x in line.split("|")}
                exon = attr["location"]
                chro, pos = exon.split(":")
                refstart, refstop = pos.split("-")
                attr["refstart"], attr["refstop"] = refstart, refstop
                for i in xrange(int(refstart), int(refstop) + 1):
                    nts += [i]
                attr["nts"] = nts
                all_nts.extend(nts)  # mutation mode
            else:
                ref_seq += line.strip()
                attr["seq"] = ref_seq
                all_ref_seq += ref_seq  # mutation mode
            try:
                strand = attr["strand"]
            except KeyError:
                exons = False
                strand = None
            try:
                ref_attributes[attr["name"] + "e" + attr["n"]] = attr
                if strand is None:
                    sys.exit("ERROR: Strand is not specified for exon n=%s\n" % attr["n"])
            except KeyError:
                ""
    ref_attributes["all_nts"] = all_nts
    ref_attributes["all_ref"] = all_ref_seq
    
    return (chro, strand, ref_attributes)


def create_report(args):
    variants = {}
    samples = {}
    data = {}
    mode = "mutation"
    header = False
    
    if args.target:
        (chro, strand, attributes) = init_ref_seq(args.target)
    
    if args.annot:
        print_vcf_header()
        vcf = True
    elif args.table:
        vcf=False
    else:
        print_line("Sample", "Region", "Location", "Type", "Removed",
                   "Added", "Abnormal", "Normal", "Ratio", "Min_coverage",
                   "Exclu_min_cov", "Variant", "Target", "Info", "Variant_sequence",
                   "Reference_sequence")
        vcf=False
    
    for line in args.infile:
        # filter header
        if line[0] == "#":
            header = True
            if "mode" in line:
                mode = line.strip().split(":")[1]
            # sys.stderr.write("Filtred: " + line)
            continue
        
        if header and not args.target:
            sys.exit("ERROR: Target file not specified. Please use -t option.")
        
        # filter on info column
        if not re.search(args.info, line):
            # sys.stderr.write("Filtred: " + line)
            continue
        
        tok = line.strip("\n").split("\t")
        if len(tok) > 1 and tok[0] != "Database" and header:
            #pathvals = {
            #    'project': '',
            #    'sample': tok[0]
            #}
            #samp = (pathvals['project'], pathvals['sample'])
            samp = tok[0]
            
            query = tok[1]
            ratio = tok[4]
            alt_exp = tok[5]
            ref_exp = tok[10]
            min_cov = tok[6]
            start_off = tok[7]
            alt_seq = tok[8]
            refSeq = tok[11]
            info = tok[-1]
            min_exclu = "" 
            
            variant = (tok[2], tok[3])
            fusion = ""
            if mode == "fusion":
                variant_name, exon = variant[0].split("/")
                if "Fusion" in variant_name:
                    fusion = "Fusion-"
                    variant_name = variant_name.replace("Fusion-", "")
                #exon = [e.split("e")[-1] for e in exon.split("::")]
                #exons = [x for e in exon for x in e.split('-')]
                print exon
                exons = [e.split("e")[0] + "e" + ee for e in exon.split("::")
                                                    for ee in e.split("e")[1].split("-")]
                print exons
                #gene = attributes[exons[0]]["name"]
                #exon = "/" + "::".join([gene + "e" + e for e in exons])
                exon = "/" + exon
                if strand == "-":
                    nts = [n for e in exons for n in attributes[e]["nts"][::-1]]
                else:
                    nts = [n for e in exons for n in attributes[e]["nts"]]
            else:
                variant_name = variant[0]
                exon = ""
                nts = attributes["all_nts"]
                if strand is None:
                    nts = [-12]*len(nts)
                if strand == "-":
                    nts = nts[::-1]
            
            if int(min_cov) < args.min_cov:
                continue
            
            if variant_name == 'Reference' or variant_name == 'Fusion':
                if strand == "-":
                    region = "{}:{}-{}".format(chro, nts[-1], nts[0])
                else:
                    region = "{}:{}-{}".format(chro, nts[0], nts[-1])
                if not vcf and not args.table:
                    print_line(samp, region, '-', variant_name + exon, '0', '0',
                               '0.0', alt_exp, tok[4], min_cov, min_exclu, '-',
                               query, tok[-1], "", "")
                continue
            
            start, mod, stop = variant[1].split(":")
            delet, insert = mod.split("/")
            
            added = str(len(insert))
            removed = str(len(delet))
            
            pos = int(start) - 1
            pos -= int(start_off)
            end = int(stop) - 1 -1
            end -= int(start_off)
            
            if strand == "-":
                end_pos = nts[pos - 1] - 1
                start_pos = nts[end - 1] - 1
            else:
                start_pos = nts[pos - 1] + 1
                end_pos = nts[end - 1] + 1
            
            region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)
            
            if len(delet) == 0 and len(insert) != 0:
                insert_type = "Insertion"
                end += 1
                
                if strand == "-":
                    end_pos = nts[pos - 1] - 1
                    start_pos = nts[end - 1] - 1
                else:
                    start_pos = nts[pos - 1] + 1
                    end_pos = nts[end - 1] + 1
                
                region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)
                
                # Reinterpret mutations for small ITDs
                # careful, going upstream may put us outside the reference.
                # I&I discovered this way do not take into account insertions and deletions.
                if insert[:int(len(insert)/2)] == insert[int(len(insert)/2):]:
                    possible_itd = insert[:int(len(insert)/2)] 
                else:
                    possible_itd = insert
                upstream = alt_seq[pos-len(possible_itd):pos]
                match = 0
                if pos-len(possible_itd) >= 0:
                    for i in xrange(0, len(possible_itd)):
                        if possible_itd[i] == upstream[i]:
                            match += 1
                    match = float(match)/len(possible_itd)
                
                if pos-len(insert) >= 0 and len(insert) >= 3 and insert == upstream:
                    insert_type = "ITD"
                    added += " | " + str(end_pos - start_pos + 1)
                elif pos-len(insert) >= 0 and len(insert) >= 3 and match > 0.5:
                    insert_type = "I&I"
                    added += " | " + str(end_pos - start_pos + 1)
                insert_type += exon
                
                location = chro + ":" + str(end_pos)

            elif variant_name == 'Substitution':
                region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)
                location = chro + ":" + str(start_pos)
                insert_type = variant_name + exon
            elif variant_name == 'Deletion':
                region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)
                location = chro + ":" + str(start_pos)
                insert_type = variant_name + exon
            elif variant_name == 'Indel':
                region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)
                location = chro + ":" + str(start_pos)
                insert_type = variant_name + exon
            else:
                sys.stderr.write("WARNING: This variant isn't taken account\n")
                sys.stderr.write(" - variant: " + str(variant_name + exon) + "\n")
                sys.stderr.write(" - line: " + line)
                sys.exit()
            
            insert_type = fusion + insert_type
        
        elif not header:
            samp = tok[0]
            region = tok[1]
            location = tok[2]
            insert_type = tok[3]
            removed = tok[4]
            added = tok[5]
            alt_exp = tok[6]
            ref_exp = tok[7]
            ratio = tok[8]
            min_cov = tok[9]
            min_exclu = tok[10]
            mod = tok[11]
            query = tok[12]
            info = tok[13]
            alt_seq = tok[14]
            refSeq = tok[15]
            
            variant = (insert_type, )
        
        if vcf:
            ref_var =  refSeq[int(start)-1-1] + delet + refSeq[int(stop)-1]
            alt_var =  refSeq[int(start)-1-1] + insert + refSeq[int(stop)-1]
            print_vcf_line(region, ref_var, alt_var, insert_type, query, ratio, min_cov)
        elif args.table:
            ""
        else:
            print_line(samp, region, location, insert_type,
                       removed, added, alt_exp, ref_exp, ratio,
                       min_cov, min_exclu, mod, query, info,
                       alt_seq, refSeq)
       
        if "/" not in variant[0]:
            var = (variant[0] + "/" + query, region + ":" + mod)
        else:
            var = (variant[0], region + ":" + mod)
        
        if var not in variants:
            variants[var] = 0
        variants[var] += 1
        
        if samp not in samples:
            samples[samp] = set()
        samples[samp].add(var)
       
        if samp not in data:
            data[samp] = {}
        data[samp][var] = float(ratio)
        
        if args.exclu != "" and alt_seq != "":
            res = uc.get_cov(args.exclu, alt_seq)
            min_exclu = str(res[2])
    
    if not args.table:
        return
    
    sorted_variants = sorted(variants, key=variants.get, reverse=True)

    sys.stdout.write("Sample")
    for v in sorted_variants:
        if "Reference" in v[0]:
            sys.stdout.write("\t" + v[0])
        else:
            sys.stdout.write("\t" + v[1])
    sys.stdout.write("\n")

    for s, sv in samples.iteritems():
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

    if args.infile.isatty():
        argparser.print_help()
        sys.exit()

    create_report(args)
