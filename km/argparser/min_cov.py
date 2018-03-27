from .common import *


def get_argparser_min_cov(parser):
    parser.add_argument(
        "target_fn",
        help="Filename of the target sequence file or directory."
        )

    parser.add_argument(
        "jellyfish_fn",
        help="Filename of the jellyfish database.",
        nargs='*')
