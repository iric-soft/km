import os
import sys

from .. utils import common as uc


def find_kmin(ref_name, ref_seq, start):
    # -1: This variant is incremented at the beginning of the loop
    k_len = start - 1
    uniq = False
    linear = False
    while (not uniq or not linear) and k_len < len(ref_seq):
        k_len += 1

        try:
            target_kmers = uc.get_target_kmers(ref_seq, k_len, ref_name)
            uniq = True
        except ValueError:
            continue

        cpt_kmer = 0
        for i in range(len(target_kmers)):
            cpt_forward = 0
            cpt_backward = 0
            for j in range(len(target_kmers)):
                if i != j:
                    kmer_i = target_kmers[i]
                    kmer_j = target_kmers[j]

                    if kmer_i[1:len(kmer_i)] == kmer_j[0:(len(kmer_j)-1)]:
                        cpt_forward += 1
                    if kmer_i[0:(len(kmer_i)-1)] == kmer_j[1:len(kmer_j)]:
                        cpt_backward += 1

                    if cpt_forward > 1 or cpt_backward > 1:
                        break

            if cpt_forward > 1 or cpt_backward > 1:
                break
            else:
                cpt_kmer += 1

        if cpt_kmer == len(target_kmers):
            linear = True

    sys.stdout.write(ref_name + "\t" + str(k_len) + "\n")


def main_linear_kmin(args, argparser):

    sys.stdout.write("target_name\tlinear_kmin\n")

    seq_files = uc.target_2_seqfiles(args.target_fn)

    for seq_f in seq_files:

        (ref_name, ext) = os.path.splitext(os.path.basename(seq_f))

        ref_seq = uc.file_2_seq(seq_f)
        find_kmin(ref_name, ref_seq, args.start)
