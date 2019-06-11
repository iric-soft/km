import sys
import re
from string import maketrans  # won't be needed in python3 (use str.maketrans)

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
                      "TYPE="+type_var+";TARGET="+target+";RATIO="+ratio+";MINCOV="+min_cov+\
                      ";REMOVED="+str(rem)+";ADDED="+str(ad)])
    sys.stdout.write(line + "\n")


def init_ref_seq(arg_ref):
    if not arg_ref:
        sys.exit("Target file is empty\n")

    # BE 1-based
    ref_attributes = {}
    all_nts = []
    all_ref_seq = []
    strand_list = []
    chro_list = []
    
    # CODE for multiple >LOC lines (find_mutation) or for exons (find_fusion)
    attkeys = None
    for line in open(arg_ref, "r"):
        line = line.strip()
        
        # Parse attributes
        if line[0] == '>':
            nts = []
            ref_seq = []
            
            loc = line.split(" ")[0]
            if "chr" not in loc or ":" not in loc or "-" not in loc:
                sys.exit('ERROR: Fasta entries do not contain a correctly ' +\
                         'formatted location: {}\n'.format(loc))
            line = line.replace(">", "location=", 1)
            # look up attributes in fasta file
            attr = {x.split("=")[0].strip():x.split("=")[1].strip() for x in line.split("|")}
            attk = '|'.join(attr.keys())
            if attkeys is None:
                #sys.stderr.write("Found attributes: " + attk + '\n')
                attkeys = attk
            elif attkeys != attk:
                sys.exit("ERROR: Attributes do not match in multi-fasta file (-t)\n")
            
            exon = attr["location"]
            chro, pos = exon.split(":")
            refstart, refstop = pos.split("-")
            
            # Get nt coordinates on the genome
            if 'strand' not in attr.keys():
                attr['strand'] = '+'
                sys.stderr.write("WARNING: Strand is assumed to be '+' \n")
            strand = attr["strand"]
            if 'cigar' not in attr.keys():
                for i in xrange(int(refstart), int(refstop) + 1):
                    nts += [i]
            else:
                cigar = attr["cigar"]
                if strand == "-":
                    cigar = re.split("([^\d])", cigar)[:-1][::-1]
                    cigar = ''.join([x for i in range(0, len(cigar)-1, 2)
                                       for x in [cigar[i+1], cigar[i]]])
                matches = []
                for i, mm in enumerate(re.split("[^\d]", cigar)[:-1]):
                    for j in range(int(mm)):
                        if i % 2:
                            matches.append(False)
                        else:
                            matches.append(True)
                for ind, i in enumerate(xrange(int(refstart), int(refstop) + 1)):
                    if matches[ind]:
                        nts += [i]
            
            attr["refstart"], attr["refstop"], attr["nts"] = refstart, refstop, nts
            all_nts.extend(nts)  #  for mutation mode
            
        # Parse sequence
        else:
            ref_seq += line.strip()
            attr["seq"] = ref_seq
            all_ref_seq += ref_seq  # for mutation mode
        
        try:
            ref_attributes[attr["name"] + "e" + attr["n"]] = attr
        except KeyError:
            ''  # We're in mutation mode
        
        strand_list.append(strand)
        chro_list.append(chro)
    
    if not ref_attributes:  # we're in mutation mode
        ref_attributes["all_nts"] = all_nts
        ref_attributes["all_ref"] = all_ref_seq
    
    assert len(set(chro_list)) == 1
    chro = chro_list[0]
    assert len(set(strand_list)) == 1
    strand = strand_list[0]
    
    return (chro, strand, ref_attributes)


def create_report(args):
    # Find correct extremities of a mutation
    sys.setrecursionlimit(10000)
    def get_bong_bong(var, p, rs):  # p = pos, rs = ref_seq
        if p - 1 > 0 and rs[p - 1] == var[-1]:
            return get_bong_bong(rs[p - 1] + var[:-1], p - 1, rs)
        return p - 1
    
    if args.format == "vcf" and args.info == "cluster":
        sys.exit("ERROR: -f vcf and -i cluster options are incompatible")
    
    variants = {}
    samples = {}
    data = {}
    mode = "mutation"
    header = False
    vcf = False
    table = False
    
    if args.target:
        (chro, strand, attributes) = init_ref_seq(args.target)
    
    if args.format == 'vcf':
        print_vcf_header()
        vcf = True
    elif args.format == 'table':
        table = True
    else:
        print_line("Sample", "Region", "Location", "Type", "Removed",
                   "Added", "Abnormal", "Normal", "Ratio", "Min_coverage",
                   "Exclu_min_cov", "Variant", "Target", "Info", "Variant_sequence",
                   "Reference_sequence")
    
    for line in args.infile:
        # filter header and initialize mode
        if line[0] == "#":
            header = True
            if "mode" in line:
                mode = line.strip().split(":")[1]
            #sys.stderr.write("Filtered: " + line)
            continue
        
        # sanity check
        if header and not args.target:
            sys.exit("ERROR: Target file not specified for km mode '{}'. " +\
                     "Please use -t option.".format(mode))
        
        # filter on info column
        if not re.search(args.info, line):
            # sys.stderr.write("Filtered: " + line)
            continue
        
        tok = line.strip("\n").split("\t")
        
        if not header:
            # keep everything as is
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
        
        elif len(tok) > 1 and tok[0] != "Database" and header:
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
            info = tok[12]
            
            min_exclu = "" 
            variant = (tok[2], tok[3])
            ref_seq = refSeq.upper()
            
            if int(min_cov) < args.min_cov:
                continue
            
            # get nts and other values based on mode
            fusion = ""
            if mode == "fusion":
                variant_name, exon = variant[0].split("/")
                if "Fusion" in variant_name:
                    fusion = "Fusion-"
                    variant_name = variant_name.replace("Fusion-", "")
                all_exons = exon.split('|')
                exon = "/"
                for ex in all_exons:
                    # any value exons takes should return the same nts
                    exons = [e.split("e")[0] + "e" + ee for e in ex.split("::")
                                                        for ee in e.split("e")[1].split("-")]
                    exon += "::".join(exons) + "|"
                exon = exon.rstrip("|")
                if strand == "-":
                    nts = [n for e in exons for n in attributes[e]["nts"][::-1]]
                    nts = nts[int(start_off):int(start_off)+len(ref_seq)]
                else:
                    nts = [n for e in exons for n in attributes[e]["nts"]]
                    nts = nts[int(start_off):int(start_off)+len(ref_seq)]
                variant = (variant_name + exon, variant[1])
                
            elif mode == "mutation":
                variant_name = variant[0]
                exon = ""
                nts = attributes["all_nts"]
                if strand == "+":
                    nts = [-12]*len(nts)
                if strand == "-":
                    nts = nts[::-1]
            assert len(nts) == len(ref_seq)
            
            seen = {}
            duplicate_count = {}
            cntd = True
            broke = False
            start_pos = 0
            if len(nts) > len(set(nts)):
                start, mod, stop = variant[1].split(":")
                delet, insert = mod.split("/")
                pos = int(start) - 1
                pos -= int(start_off)
                end = int(stop) - 1 - 1
                end -= int(start_off)
                
                var = delet.upper()
                ibef = get_bong_bong(delet, pos, ref_seq)
                before = ref_seq[ibef:pos]
                iaft = get_bong_bong(var[::-1], len(ref_seq)-pos-1-len(var)+1, ref_seq[::-1])
                after = ref_seq[::-1][iaft:len(ref_seq)-pos-1-len(var)+1][::-1]
                iaft = len(ref_seq) - iaft - 1
                
                start = ibef + 1
                stop = iaft + 2
                delet = before + delet + after
                delet = delet.lower()
                insert = before + insert + after
                
                new_nts = []
                for nt, n in zip(nts, ref_seq):
                    if nt not in seen:
                        if broke:
                            cntd = False
                        seen[nt] = start_pos
                        new_nts.append(nt)
                    else:
                        duplicate_count[start_pos] = [2]  # will point to same list with seen[nt]
                        duplicate_count[seen[nt]] = duplicate_count[start_pos]
                        broke = True
                        if cntd == False:
                            sys.stderr.write("nts sequence doesn't make sense")
                            exit(1)
                    start_pos += 1
                
                false_del = ""
                for spos, n in zip(range(start - 1, stop - 1), delet.upper()):
                    if spos in duplicate_count and duplicate_count[spos] == [2]:
                        false_del += n
                        duplicate_count[spos][0] -= 1
                
                nts = new_nts
                old_delet = delet
                delet = delet.replace(false_del.lower(), "", 1)
                assert old_delet != delet
                for i, j in zip(delet.upper(), insert.upper()):
                    if i == j:
                        delet = delet[1:]
                        insert = insert[1:]
                        start += 1
                    else:
                        break
                for i, j in zip(delet[::-1].upper(), insert[::-1].upper()):
                    if i == j:
                        delet = delet[:-1]
                        insert = insert[:-1]
                        stop -= 1
                    else:
                        break
                reduct_len = len(false_del)
                reduct_offset = max(0, start - min(duplicate_count.keys()))
                ref_seq = ref_seq.replace(false_del, "", 1)
                if not delet and not insert:
                    variant_name = "Reference"
                    var = ""
                elif delet and insert:
                    variant_name = "Indel"
                    cur_start = start - reduct_offset
                    cur_end = stop - len(false_del) - reduct_offset
                    if len(delet) == len(insert):
                        variant_name = "Substitution"
                    var = '%d:%s/%s:%d' % (cur_start, delet, insert, cur_end)
                elif delet and not insert:
                    variant_name = "Deletion"
                    var = '%d:%s/%s:%d' % (start - reduct_offset, delet, insert,
                                           stop - len(false_del) - reduct_offset)
                elif not delet and insert:
                    variant_name = "Insertion"
                    var = '%d:%s/%s:%d' % (start - reduct_offset, delet, insert,
                                           stop - len(false_del) - reduct_offset)
                variant = (variant_name + exon, var)
            
            assert len(set(nts)) == len(nts)
            
            # case: entries with no mutations
            if variant_name == 'Reference' or variant_name == 'Fusion':
                mod = ""
                if strand == "-":
                    region = "{}:{}-{}".format(chro, nts[-1], nts[0])
                else:
                    region = "{}:{}-{}".format(chro, nts[0], nts[-1])
                if not vcf and not table:
                    print_line(samp, region, '-', variant_name + exon, '0', '0',
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
                end = int(stop) - 1 - 1
                end -= int(start_off)
                
                if strand == "-":
                    end_pos = nts[pos - 1] - 1
                    start_pos = nts[end - 1] - 1
                else:
                    start_pos = nts[pos - 1] + 1
                    end_pos = nts[end - 1] + 1
                
                ref_var =  delet.upper()
                alt_var =  insert.upper()
                loc_var = start_pos
                end_var = end_pos
                
                region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)
                
                if len(delet) == 0 and len(insert) != 0:
                    insert_type = "Insertion"
                    
                    var = insert.upper()
                    # include current position in before
                    ibef = get_bong_bong(var, pos + 1, ref_seq)
                    before = ref_seq[ibef:pos+1]
                    iaft = get_bong_bong(var[::-1], len(ref_seq)-pos-1, ref_seq[::-1])
                    after = ref_seq[::-1][iaft:len(ref_seq)-pos-1][::-1]
                    iaft = len(ref_seq) - iaft - 1
                    ref_var = before + after
                    alt_var = before + var + after
                    loc_var = nts[iaft] if strand == "-" else nts[ibef]
                    end_var = nts[iaft-len(ref_var)+1] if strand == "-" else nts[ibef+len(ref_var)-1]
                    
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
                    
                    if loc_var + len(ref_var) - 1 < end_var:
                        insert_type = 'PossibleIntron' + exon
                    
                    location = chro + ":" + str(end_pos)
                    
                elif variant_name == 'Substitution':
                    region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)
                    location = chro + ":" + str(start_pos)
                    insert_type = variant_name + exon
                    
                    # NOTE: Substitutions that end at junctions are considered mutations not altsplice
                    if loc_var + len(ref_var) - 1 < end_var:
                        insert_type = 'PossibleAltsplice' + exon
                    
                elif variant_name == 'Deletion':
                    region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)
                    location = chro + ":" + str(start_pos)
                    insert_type = variant_name + exon
                    
                    var = delet.upper()
                    ibef = get_bong_bong(var, pos, ref_seq)
                    before = ref_seq[ibef:pos]
                    iaft = get_bong_bong(var[::-1], len(ref_seq)-pos-1-len(var)+1, ref_seq[::-1])
                    after = ref_seq[::-1][iaft:len(ref_seq)-pos-1-len(var)+1][::-1]
                    iaft = len(ref_seq) - iaft - 1
                    ref_var = before + var + after
                    alt_var = before + after
                    loc_var = nts[iaft] if strand == "-" else nts[ibef]
                    end_var = nts[iaft-len(ref_var)+1] if strand == "-" else nts[ibef+len(ref_var)-1]
                    
                    if loc_var + len(ref_var) - 1 < end_var:
                        insert_type = 'Altsplice' + exon
                    
                elif variant_name == 'Indel':
                    region = "{}:{}-{}".format(chro, start_pos, end_pos + 1)
                    location = chro + ":" + str(start_pos)
                    insert_type = variant_name + exon
                    
                    ref_var =  ref_seq[pos-1] + delet.upper() + ref_seq[end + 1]
                    alt_var =  ref_seq[pos-1] + insert.upper() + ref_seq[end + 1]
                    loc_var = start_pos - 1
                    end_var = end_pos + 1
                    
                    if loc_var + len(ref_var) - 1 < end_var:
                        insert_type = 'PossibleAltspliceIntron' + exon
                    
                else:
                    sys.stderr.write("WARNING: This variant isn't taken account\n")
                    sys.stderr.write(" - variant: " + str(variant_name + exon) + "\n")
                    sys.stderr.write(" - line: " + line)
                    sys.exit()
                
                insert_type = fusion + insert_type
        
        if args.exclu != "" and alt_seq != "":
            res = uc.get_cov(args.exclu, alt_seq)
            min_exclu = str(res[2])
       
        if not vcf and not table:
            print_line(samp, region, location, insert_type,
                       removed, added, alt_exp, ref_exp, ratio,
                       min_cov, min_exclu, mod, query, info,
                       alt_seq, refSeq)
            
        elif vcf and header:
            complement = maketrans('ATGCU', 'TACGA')
            ref_var = ref_var.translate(complement)[::-1] if strand == '-' else ref_var
            alt_var = alt_var.translate(complement)[::-1] if strand == '-' else alt_var
            if '/' in insert_type:
                insert_type, query = insert_type.split('/')
            try:
                assert loc_var + len(ref_var) - 1 <= end_var
            except AssertionError:
                sys.stderr.write(
                        "WARNING: Throwing away an aberration on {} for {}\n".format(chro, query)
                        )
                continue
            if '::' not in query or args.junction:
                print_vcf_line(chro, loc_var, ref_var, alt_var, insert_type,
                               query, ratio, min_cov, removed, added)
            
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
            if v[0].split("/")[0] == "Reference" or v[0].split("/")[0] == "Fusion":
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
