# find_mutation ---
#
#   Usage:  find_mutation <region_fasta or directory> <jellyfish_db>
import os
import sys
import time
import logging as log
from .. utils import MutationFinder as umf
from .. utils import common as uc
from .. utils import Sequence as us
from .. utils.Jellyfish import Jellyfish


# ###########################################################################
# Main function

def main_find_mut(args, argparser):
    time_start = time.time()

    if args.verbose:
        log.basicConfig(level=log.INFO, format="VERBOSE: %(message)s")

    if args.debug:
        log.basicConfig(level=log.DEBUG, format="VERBOSE: %(message)s")

    for k, v in vars(args).items():
        sys.stdout.write("#" + str(k) + ':' + str(v) + "\n")

    jf = Jellyfish(args.jellyfish_fn, cutoff=args.ratio, n_cutoff=args.count)

    seq_files = uc.target_2_seqfiles(args.target_fn)

    umf.MutationFinder.output_header()

    refpaths = []

    for seq_f in seq_files:

        (ref_name, ext) = os.path.splitext(os.path.basename(seq_f))

        ref_seqs, ref_attr = uc.file_2_seq(seq_f)

        refpath = us.RefSeq(''.join(ref_seqs), ref_name, jf.k)

        refpaths.append(refpath)

    for refpath in refpaths:

        finder = umf.MutationFinder(
            refpath, jf, args.steps, args.branchs, args.nodes
        )

        finder.graph_analysis()
        finder.quantify_paths(args.graphical)
        finder.quantify_clusters(args.graphical)

        for path in finder.get_paths(sort=True):
            sys.stdout.write(str(path) + "\n")

    sys.stdout.write("#Elapsed time:" + str(time.time() - time_start) + "\n")
