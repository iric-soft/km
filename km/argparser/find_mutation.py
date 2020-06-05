from .common import *


def get_argparser_find_mut(parser):
    parser.add_argument(
        "-c", "--count",
        help="Minimum occurence needed for exploration of alternative (default: -c 5)",
        action="store",
        nargs='?',
        default=5,
        type=int)
    parser.add_argument(
        "-p", "--ratio",
        help="Minimum occurence ratio needed for exploration of alternative (default: -p 0.05)",
        action="store",
        nargs='?',
        default=0.05,
        type=float)
    parser.add_argument(
        "-s", "--steps",
        help="Maximum steps to discover a new branch on a target sequence (default: -s 500)",
        action="store",
        nargs='?',
        default=500,
        type=int)
    parser.add_argument(
        "-b", "--branchs",
        help="Maximum branchs until getback to target sequence (default: -b 10)",
        action="store",
        nargs='?',
        default=10,
        type=int)
    parser.add_argument(
        "-g", "--graphical",
        help="Display coverage graph.",
        action="store_true")
    parser.add_argument(
        "-v", "--verbose",
        help="Get more information.",
        action="store_true")
    parser.add_argument(
        "target_fn",
        help="Filename of the target sequence file or directory.",
        nargs='*')
    parser.add_argument(
        "jellyfish_fn",
        help="Filename of the jellyfish database.")
