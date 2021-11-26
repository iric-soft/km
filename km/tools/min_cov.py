# min_cov ---
#
import os
import sys
from .. utils import common as uc


# ###########################################################################
# Main function
def main_min_cov(args, argparser):
    lst_files = uc.args_2_list_files(args.jellyfish_fn)

    ref_seq = args.target_fn
    if os.path.isfile(args.target_fn):
        ref_seq, ref_attr = uc.file_2_seq(args.target_fn)
        ref_seq = ''.join(ref_seq)

    sys.stdout.write("DB\tcount\tlength\tmin\tmax\tmean\tkmer_nb\tkmer_nb_0\n")

    for jf_file in lst_files:
        res = uc.get_cov(jf_file, ref_seq)
        res_str = "%s\t%d\t%d\t%d\t%d\t%.2f\t%d\t%d" % (
                  jf_file, res[0], res[1], res[2], res[3],
                  res[4], res[5], res[6])
        sys.stdout.write(res_str + "\n")
