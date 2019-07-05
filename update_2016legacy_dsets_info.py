import argparse
import logging
import textwrap
from yaml import load, dump
import os


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "update dsets_info with legacy 2016 MC",
    epilog = "Example:\n$ python update_2016legacy_dsets_info.py `cat dsets_ext_wjets_dyjets_2016legacy dsets_all_usual_mc_2016legacy`"
    )

parser.add_argument('dsets', nargs='+', help="""the to update if they are not present in the --dsets-file""")
parser.add_argument("--dsets-file", type=str, default="dsets_info.yaml", help="the target file with dsets info to update")

parser.add_argument("--debug", action="store_true", help="debug logging")
#parser.add_argument("-o", "--output-dir", type=str, default='python/crab_cfgs/', help="directory where config files are created for this vertsion (in version/ subdir)")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

def dset_firstname(dset):
    return dset.split('/')[1]


# get the existing info file
with open(args.dsets_file) as f:
    dsets_info = load(f)


# the target info
update_info = {}

for dset in args.dsets:
    if dset in dsets_info:
        logging.info("the dset is present: %s" % dset)
        continue

    # update, based on the dset first name
    first_name = dset_firstname(dset)

    # find the corresponding dset in the info
    for existing_dset in dsets_info:
        e_dset_first_name = dset_firstname(existing_dset)

        if e_dset_first_name == first_name:
            # get the info
            update_info[dset] = dsets_info[existing_dset]
            dsets_info.pop(existing_dset)
            break

    # if no info is found report and continue
    if dset not in update_info:
        logging.warning("no existing info for the dataset: %s" % dset)


# just print the update info
print dump(update_info)

