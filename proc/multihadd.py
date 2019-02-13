import argparse
import logging
import os
from itertools import product


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "run hadds filenames with several strings",
    epilog = "Example:\n$ python multihadd.py quick-test/v36-test1/ quick-test/v36-test1_hadds/ MC2016_Summer16_TTJets_powheg-MC2016_Summer16_DYJetsToLL_50toInf_madgraph:el_sel,-el_sel_lj,-el_sel_ljout"
    )

parser.add_argument('inp_dir',    help="input directory, the filed to hadd")
parser.add_argument('out_dir',    help="output directory, the hadd results")
parser.add_argument('defs',       type=str, help="define hadd groups, def1[:def2..], each def = group1[:group2..], the hadd result is direct product of the definitions")
parser.add_argument('--hadd-first',    action='store_true', help="also hadd results by the first group (like by the dtag)")

parser.add_argument("-d", "--debug",    action='store_true', help="DEBUG level of logging")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

groups = product(*[d.split('-') for d in args.defs.split(':')])

#logging.debug(repr(list(groups)))

for gr in groups:
    logging.debug(repr(gr))
    inp_files = args.inp_dir + '/*%s*root' % '*'.join(gr)
    out_file  = args.out_dir + '/%s.root'  % '_'.join(gr)
    logging.debug(out_file)
    os.system('hadd %s %s' % (out_file, inp_files))

if args.hadd_first:
    for first_group in [d for d in args.defs.split(':')[0].split('-')]:
        grp_file  = args.out_dir + '/%s.root' % first_group
        if os.path.isfile(grp_file):
            logging.debug('first group file %s exists, skipping' % grp_file)
            continue
        grp_files = args.out_dir + '/%s*.root' % first_group
        logging.debug(grp_file)
        os.system('hadd %s %s' % (grp_file, grp_files))

