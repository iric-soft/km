import argparse

from .argparser.find_mutation import *
from .argparser.find_report import *
from .argparser.linear_kmin import *
from .argparser.min_cov import *

from .tools.find_mutation import main_find_mut
from .tools.find_report import main_find_report
from .tools.linear_kmin import main_linear_kmin
from .tools.min_cov import main_min_cov


# ###########################################################################
# Main function
def main():
    # print("\n--------------------------------------------------------------")
    # print("km.py: Tools for targeted variant detection.")
    # print("This program was written by IRIC's bioinformatic platform")
    # print("----------------------------------------------------------------\n")

    argparser = argparse.ArgumentParser(prog='PROG')
    subparsers = argparser.add_subparsers(help='sub-command help')

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
        help='Parse find_mutation output to reformat it in tabulated file more user friendly.'
    )
    find_report.set_defaults(func=main_find_report)
    get_argparser_find_report(find_report)

    # create the argparser for the "linear_kmin" command
    linear_kmin = subparsers.add_parser(
        'linear_kmin',
        help='Find min k length to decompose a target sequence in a linear graph.'
    )
    linear_kmin.set_defaults(func=main_linear_kmin)
    get_argparser_linear_kmin(linear_kmin)

    # create the argparser for the "linear_kmin" command
    min_cov = subparsers.add_parser(
        'min_cov',
        help='Compute coverage of target sequences.'
    )
    min_cov.set_defaults(func=main_min_cov)
    get_argparser_min_cov(min_cov)

    # recover arguments
    args = argparser.parse_args()

    # execute the command
    args.func(args, argparser)
