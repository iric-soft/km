import argparse
import sys

from .argparser.find_mutation import get_argparser_find_mut
from .argparser.find_report import get_argparser_find_report
from .argparser.linear_kmin import get_argparser_linear_kmin
from .argparser.min_cov import get_argparser_min_cov

from .tools.find_mutation import main_find_mut
from .tools.find_report import main_find_report
from .tools.linear_kmin import main_linear_kmin
from .tools.min_cov import main_min_cov


# ###########################################################################
# Main function
def main():
    # sys.stderr.write("\n--------------------------------------------------------------\n")
    # sys.stderr.write("km.py: Tools for targeted variant detection.\n")
    # sys.stderr.write("This program was written by IRIC's bioinformatics platform\n")
    # sys.stderr.write("----------------------------------------------------------------\n\n")

    parser = argparse.ArgumentParser(prog='PROG')
    subparsers = parser.add_subparsers(help='sub-command help')

    # create the argparser for the "find_mutation" command
    find_mut = subparsers.add_parser(
        'find_mutation',
        help='Identify and quantify mutations from a target sequence and a k-mer database.'
    )
    find_mut.set_defaults(func=main_find_mut)
    get_argparser_find_mut(find_mut)

    # create the argparser for the "find_report" command
    find_report = subparsers.add_parser(
        'find_report',
        help='Parse find_mutation output and reformat it in a more user-friendly tabulated file.'
    )
    find_report.set_defaults(func=main_find_report)
    get_argparser_find_report(find_report)

    # create the argparser for the "linear_kmin" command
    linear_kmin = subparsers.add_parser(
        'linear_kmin',
        help='Find min k-length to decompose a target sequence in a linear graph.'
    )
    linear_kmin.set_defaults(func=main_linear_kmin)
    get_argparser_linear_kmin(linear_kmin)

    # create the argparser for the "min_cov" command
    min_cov = subparsers.add_parser(
        'min_cov',
        help='Compute coverage of target sequences.'
    )
    min_cov.set_defaults(func=main_min_cov)
    get_argparser_min_cov(min_cov)

    # display help if no arguments given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # recover arguments
    args = parser.parse_args()

    # execute the command
    args.func(args, parser)
