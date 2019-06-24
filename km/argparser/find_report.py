import sys
import argparse
from .common import *


def get_argparser_find_report(parser):
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument(
        "-t",
        dest="target",
        help="Filename of the target sequence file",
        type=lambda x: is_valid_file(parser, x))

    required.add_argument(
        'infile',
        nargs='?',
        type=argparse.FileType('r'),
        default=sys.stdin)

    optional.add_argument(
        "-i",
        dest="info",
        help="Filter on info column (Default: vs_ref)",
        default="vs_ref",
        type=str)

    optional.add_argument(
        "-m",
        dest="min_cov",
        help="Min coverage allowed (Default: 1)",
        default=1,
        type=int)

    optional.add_argument(
        "-e", "--exclu",
        dest="exclu",
        help="Filename of a jf database, containing k-mers which can create false positive variants (as, a jf build on the transcriptome)",
        default="",
        type=str)
    
    optional.add_argument(
        "-f", "--format",
        dest="format",
        help="Option 'vcf': Output variants in VCF-like file format -- Option 'table':" +\
             "Group variants by position and return per-sample ratio",
        choices=['vcf', 'table'])
