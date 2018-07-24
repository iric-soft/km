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

    new_seq = []

    for i in range(1):
        new_seq.append(ref_seq[11])
        new_seq.append(ref_seq[0])

    for i in range(len(target_seq)):
        print(ref_name[i])
        print(new_seq[i])

    finder = umf.MutationFinder(
            ref_name, new_seq, jf,
            args.graphical, args.steps, args.branchs
        )

    for path in finder.get_paths():
        sys.stdout.write(str(path) + "\n")

    sys.stdout.write("#Elapsed time:" + str(time.time() - time_start) + "\n")
