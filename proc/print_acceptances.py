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
parser.add_argument("--with-stat",   action='store_true', help="print stat error too")
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

range_length = 75 if args.sys_weights else 2

for file_i, input_filename in enumerate(args.input_files):
    if not isfile(input_filename):
        logging.info("missing: " + input_filename)
        exit(1)

    logging.debug(input_filename)

    tfile = TFile(input_filename)

    all_histos = [tfile.Get("all_ev_" + name) for name, numbers in results]
    cut_histos = [tfile.Get("cut_ev_" + name) for name, numbers in results]

    sum_all_histos = all_histos[0].Clone()
    sum_cut_histos = cut_histos[0].Clone()

    sum_all_histos.SetName("sum_all_histos")
    sum_cut_histos.SetName("sum_cut_histos")

    for histo_cut, histo_all in zip(cut_histos[1:], all_histos[1:]):
        sum_all_histos += histo_all
        sum_cut_histos += histo_cut


    for name, numbers in results:
        histo_all = tfile.Get("all_ev_" + name)
        histo_cut = tfile.Get("cut_ev_" + name)
        histo_cut.Divide(histo_all)

        if args.yields:
            numbers.append((histo_all.GetBinContent(1), histo_all.GetBinError(1)))

        for i in range(1,range_length):
            numbers.append((histo_cut.GetBinContent(i), histo_cut.GetBinError(i)))


for name, numbers in results:
    #
    if args.with_stat:
        print "%20s " % name + ' '*5, ' '.join("%10.4f +- %6.4f" % (numbers[i][0], numbers[i][1]) for i in range(len(numbers)))
    else:
        print "%20s " % name + ' '*5, ' '.join("%10.4f" % numbers[i][0] for i in range(len(numbers)))

sum_cut_histos.Divide(sum_all_histos)

sum_num = []
if args.yields:
    sum_num.append((sum_all_histos.GetBinContent(1), sum_all_histos.GetBinError(1)))

for i in range(1,range_length):
    sum_num.append((sum_cut_histos.GetBinContent(i), sum_cut_histos.GetBinError(i)))

if args.with_stat:
    print "%20s " % "sum" + ' '*5, ' '.join("%10.4f +- %6.4f" % (sum_num[i][0], sum_num[i][1]) for i in range(len(sum_num)))
else:
    print "%20s " % "sum" + ' '*5, ' '.join("%10.4f" % sum_num[i][0] for i in range(len(sum_num)))

if args.ratio:
    for name, numbers in results:
        print "%20s " % name + ' '*5, numbers[0][0] / (sum(numbers[1:][0])/(len(numbers) - 1))


