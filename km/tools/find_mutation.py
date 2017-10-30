# find_mutation ---
#
#   Usage:  find_mutation <region_fasta or directory> <jellyfish_db>
import os
import time
import utils.MutationFinder as umf
import logging as log
from utils.Jellyfish import Jellyfish


# ###########################################################################
# Main function
def main_find_mut(args, argv):
    time_start = time.time()

    if args.verbose:
        log.basicConfig(level=log.DEBUG, format="VERBOSE: %(message)s")

    print '#', " ".join(argv)
    for k, v in vars(args).iteritems():
        print '#', k, ':', v

    jf = Jellyfish(args.jellyfish_fn, cutoff=args.ratio, n_cutoff=args.count)

    # Gather file names for ref. sequences.
    if len(args.reference_fn) > 1:
        seq_files = args.reference_fn
    else:
        if os.path.isdir(args.reference_fn[0]):
            seq_files = map(
                lambda f: os.path.join(args.reference_fn[0], f),
                os.listdir(args.reference_fn[0]))
        else:
            seq_files = args.reference_fn

    umf.MutationFinder.output_header()

    for seq_f in seq_files:

        (ref_name, ext) = os.path.splitext(os.path.basename(seq_f))

        ref_seq = []
        for line in open(seq_f, "r"):
            line = line.strip()
            if line[0] == '>':
                continue
            ref_seq.append(line)
        ref_seq = ''.join(ref_seq)

        finder = umf.MutationFinder(
            ref_name, ref_seq, jf,
            args.graphical, args.steps
        )

        for path in finder.get_paths():
            print path

    print "#Elapsed time:", time.time() - time_start
