# find_mutation ---
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
def main_find_mut(args, argparser):
    time_start = time.time()

    if args.verbose:
        log.basicConfig(level=log.DEBUG, format="VERBOSE: %(message)s")

    for k, v in vars(args).iteritems():
        sys.stdout.write("#" + str(k) + ':' + str(v) + "\n")

    jf = Jellyfish(args.jellyfish_fn, cutoff=args.ratio, n_cutoff=args.count)

    seq_files = uc.target_2_seqfiles(args.target_fn)

    umf.MutationFinder.output_header()

    for seq_f in seq_files:

        #this will have to be switched to an agr parser 
        (ref_name, ext) = os.path.splitext(os.path.basename(seq_f))

        temp = ref_name
        ref_name = []

        for i in range(0, 5):
            ref_name.append(temp)

        ref_seq = []
        for i in range(0, 5):
            ref_seq.append(uc.file_2_seq(seq_f))

        
        finder = umf.MutationFinder(
                ref_name, ref_seq, jf,
                args.graphical, args.steps
            )

        
        for i in range(0, 2):
            for path in finder.get_paths(target_num = i):
                sys.stdout.write(str(path) + "\n")

    sys.stdout.write("#Elapsed time:" + str(time.time() - time_start) + "\n")
