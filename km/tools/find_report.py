import sys

from .. utils import Report as ur


def main_find_report(args, argparser):

    if args.infile.isatty():
        argparser.print_help()
        sys.exit()

    if args.target is not None:
        sys.stderr.write(
            'DEPRECATED: Target file (-t option) is ignored and ' +
            'will be removed from future versions of km.\n\n'
        )

    if args.format == "vcf" and args.info == "cluster":
        # Note: could salvage that option if we get the fill ref from vs_ref entries
        sys.exit("ERROR: -f vcf and -i cluster options are incompatible")

    if args.format == "vcf":
        cls = ur.VCFReport
    elif args.format == "table":
        cls = ur.TableReport
    else:
        cls = ur.Report

    report = cls(args.infile, args.info, args.min_cov, args.exclu)
    report.print_header()
    report.parse()
