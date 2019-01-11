# find_mutation ---
#
#   Usage:  find_mutation <region_fasta or directory> <jellyfish_db>
import os
import sys
import time
import logging as log
#from collections import OrderedDict

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
    
    if len(seq_files) != 2:
        raise IndexError("find_fusion needs exactly two fasta files (given= %d)" % len(seq_files))
    
    umf.MutationFinder.output_header()
    
    # Get exon sequences and remove duplicates, keeping input order
    (name, ext) = os.path.splitext(os.path.basename(seq_files[0]))
    ref_left_seq = uc.file_2_seq(seq_files[0], individual=True)
    ref_left_name = [name + "_exon" + str(i) for i in range(len(ref_left_seq))]
    
    (name, ext) = os.path.splitext(os.path.basename(seq_files[1]))
    ref_right_seq = uc.file_2_seq(seq_files[1], individual=True)
    ref_right_name = [name + "_exon" + str(i) for i in range(len(ref_right_seq))]
    
    finder = umf.MutationFinder(
        [ref_left_name, ref_right_name], [ref_left_seq, ref_right_seq], jf,
        args.graphical, args.steps, args.branchs, mode='fusion'
    )

    finder.altsplice = args.altsplice

    finder.graph_analysis(args.graphical)
        
    for path in finder.get_paths():
        sys.stdout.write(str(path) + "\n")

    sys.stdout.write("#Elapsed time:" + str(time.time() - time_start) + "\n")
