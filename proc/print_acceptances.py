import argparse
import logging
import os
from os.path import isfile, basename
import ctypes
from array import array
from sys import exit




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "signal acceptance printout",
    epilog = """Example:
python print_acceptances.py temp/NtuplerAnalyzer_test_METfiltersOFF_TTJets2_signal_recordAll.root
"""
    )

#parser.add_argument("input_filename",  type=str, help="the path tot he file with acceptances")
parser.add_argument("--sys-weights", action='store_true', help="print weight-based systematics too")
parser.add_argument("--yields",      action='store_true', help="print number of events processed")
parser.add_argument("--ratio",       action='store_true', help="calculate ratio to the first value")
parser.add_argument("--all-procs",   action='store_true', help="print for all processes")
parser.add_argument("--debug",       action='store_true', help="DEBUG level of logging")

parser.add_argument('input_files', nargs="+", help="""files with acceptances""")


args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

logging.info("import ROOT")

import ROOT
from ROOT import TH1D, TFile, TTree, gROOT
from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
from ROOT import TLegend
#from ROOT import THStack

logging.info("done")

if args.all_procs:
    results = [
    ("mutau"    , []),
    ("eltau"    , []),
    ("muj"      , []),
    ("elj"      , []),
    ("taumutauh", []),
    ("taueltauh", []),
    ("taumuj"   , []),
    ("tauelj"   , []),
    ("other"    , [])]
else:
    results = [
    ("mutau"    , []),
    ("eltau"    , [])]

range_length = 50 if args.sys_weights else 2

for file_i, input_filename in enumerate(args.input_files):
    if not isfile(input_filename):
        logging.info("missing: " + input_filename)
        exit(1)

    logging.debug(input_filename)

    tfile = TFile(input_filename)

    for name, numbers in results:
        histo_all     = tfile.Get("all_ev_" + name)
        histo_cut     = tfile.Get("cut_ev_" + name)
        histo_cut.Divide(histo_all)

        if args.yields:
            numbers.append(histo_all.GetBinContent(1))

        for i in range(1,range_length):
            numbers.append(histo_cut.GetBinContent(i))


for name, numbers in results:
    #
    print "%20s " % name + ' '*5, ' '.join("%10.3f" % numbers[i] for i in range(len(numbers)))

if args.ratio:
    for name, numbers in results:
        print "%20s " % name + ' '*5, numbers[0] / (sum(numbers[1:])/(len(numbers) - 1))


