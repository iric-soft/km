import sys
import argparse

from argparser.find_mutation import *

from tools.find_mutation import main_find_mut


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
        help='Identify and quantify mutations from a k-mer database.'
    )
    find_mut.set_defaults(func=main_find_mut)
    get_parser_find_mut(find_mut)

    # recover arguments
    args = argparser.parse_args()

    # execute the command
    args.func(args, sys.argv)


main()
