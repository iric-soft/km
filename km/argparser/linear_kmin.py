import sys
import argparse
from .common import *


def get_argparser_linear_kmin(parser):
    parser.add_argument("-s",
                        "--start",
                        help="starting length (default: -s 10)",
                        action="store",
                        nargs='?',
                        default=10,
                        type=int)
    parser.add_argument("target_fn",
                        help="Filename of the reference sequence file or directory.",
                        nargs='*')
