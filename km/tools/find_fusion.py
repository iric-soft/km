# find mutation ---
#
#   Usage:  find_mutation <region_fasta or directory> <jellyfish_db>
import os
import sys
import time
import logging as log
from .. utils import MutationFinder as umf
from .. utils import common as uc
from .. utils.Jellyfish import Jellyfish


# ###########################################################################
# Main function
def main_find_fus(args, argparser):
    time_start = time.time()

    if args.verbose:
        log.basicConfig(level=log.DEBUG, format="VERBOSE: %(message)s")

    for k, v in vars(args).iteritems():
        sys.stdout.write("#" + str(k) + ':' + str(v) + "\n")

    jf = Jellyfish(args.jellyfish_fn, cutoff=args.ratio, n_cutoff=args.count)

    seq_files = uc.target_2_seqfiles(args.target_fn)

    umf.MutationFinder.output_header()

    ref_name= uc.file_2_fus_names(seq_files)

    ref_seq = uc.exons_2_fusion_seq(seq_files)

    first_gene_ref = uc.file_2_exon_list(seq_files[0])

    finder = umf.MutationFinder(
            ref_name, ref_seq, jf,
            args.graphical, first_gene_ref, args.steps, args.branchs
        )

    for path in finder.get_paths():
        sys.stdout.write(str(path) + "\n")

    sys.stdout.write("#Elapsed time:" + str(time.time() - time_start) + "\n")
