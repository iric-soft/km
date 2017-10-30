import os


def is_valid_path(parser, p_file):
    p_file = os.path.abspath(p_file)
    if not os.path.exists(p_file):
        parser.error("The path %s does not exist!" % p_file)
    else:
        return p_file


def is_valid_file(parser, n_file):
    if not os.path.isfile(n_file):
        parser.error("The file %s does not exist!" % n_file)
    else:
        return n_file
