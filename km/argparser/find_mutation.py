from argparser.common import *


def get_parser_find_mut(argparser):
    argparser.add_argument(
        "-c", "--count",
        help="Minimum occurence needed for exploration of alternative (default: -c 5)",
        action="store",
        nargs='?',
        default=5,
        type=int)
    argparser.add_argument(
        "-p", "--ratio",
        help="Minimum occurence ratio needed for exploration of alternative (default: -p 0.05)",
        action="store",
        nargs='?',
        default=0.05,
        type=float)
    argparser.add_argument(
        "-s", "--steps",
        help="Maximum steps to discover a new branch on a target sequence (default: -s 500)",
        action="store",
        nargs='?',
        default=500,
        type=int)
    argparser.add_argument(
        "-g", "--graphical",
        help="Display coverage graph.",
        action="store_true")
    argparser.add_argument(
        "-v", "--verbose",
        help="Get more information.",
        action="store_true")
    argparser.add_argument(
        "reference_fn",
        help="Filename of the reference sequence file or directory.",
        nargs='*')
    argparser.add_argument(
        "jellyfish_fn",
        help="Filename of the jellyfish database.")
