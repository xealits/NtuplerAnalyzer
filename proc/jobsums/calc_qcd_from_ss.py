import argparse
import logging
from os.path import isfile

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "copy selected distributions for all processes and systematics, calculate QCD from Same Signal region if needed and store in the output file",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument("input_file",  help="input filename")
parser.add_argument("output_file", help="output filename")
parser.add_argument("-c", "--channel",     type=str, default='mu_sel',
    help="selection channels of events to consider (can be a list of channels for shape evolution study)")
parser.add_argument("-d", "--distr-name",  type=str, default='Mt_lep_met',
    help='distribution to consider')

parser.add_argument("--qcd",   type=float, default=1., help="get QCD from corresponding _ss channel and transfer it with the given factor (try --qcd 1)")


args = parser.parse_args()

assert isfile(args.input_file)

logging.info("import ROOT")

import ROOT
from ROOT import TFile, THStack, TLegend, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
#import plotting_root


'''
for the channel loop over all processes and systematics
copying the required distribution
'''


