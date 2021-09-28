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

    refpaths_list = []

    for seq_f in seq_files:

        (ref_name, ext) = os.path.splitext(os.path.basename(seq_f))

        ref_seqs, ref_attr = uc.file_2_seq(seq_f)

        refpaths = []
        for seq, att in zip(ref_seqs, ref_attr):
            if '_filename' in ref_attr and ref_name != ref_attr['_filename']:
                sys.stderr.write('Fasta header _filename field will be overwritten.\n')
            att['_filename'] = ref_name

            assert len(seq) >= jf.k
            refseq = us.RefSeq(seq, att, jf.k)
            refpaths.append(refseq)

            sys.stdout.write("#target:" + str(refseq) + '\n')

        refpaths_list.append(refpaths)

    umf.MutationFinder.output_header()

    for refpaths in refpaths_list:

        finder = umf.MutationFinder(
            refpaths, jf, args.steps, args.branchs, args.nodes
        )

        finder.graph_analysis()
        finder.quantify_paths(args.graphical)
        finder.quantify_clusters(args.graphical)

        for path in finder.get_paths(sort=True):
            sys.stdout.write(str(path) + "\n")

    sys.stdout.write("#Elapsed time:" + str(time.time() - time_start) + "\n")
