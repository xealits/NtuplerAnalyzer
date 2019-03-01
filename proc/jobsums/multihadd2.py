import argparse
import logging
import os
from itertools import product


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "produce hadd commands to run on groups of filenames, use with xargs",
    epilog = "Example:\n$ python multihadd2.py -n 10 out.root in1.root in2.root"
    )

parser.add_argument('-n', '--number-of-groups', type=int, help="number of groups")
parser.add_argument('-i', '--items-per-groups', type=int, help="number of items per group")
parser.add_argument('-d', '--debug',    action='store_true', help="DEBUG level of logging")

parser.add_argument('out_file',   type=str, help="output directory, the hadd results")
parser.add_argument('inp_files', nargs='+', help="the files to hadd")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

def chunks_of_n_items(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))

if args.number_of_groups:
    groups = split(args.inp_files, args.number_of_groups)
elif args.items_per_groups:
    groups = chunks_of_n_items(args.inp_files, args.items_per_groups)

for i, gr in enumerate(groups):
    logging.debug(repr(gr))
    command = 'hadd %s.%d %s' % (args.out_file, i, ' '.join(gr))
    print command

