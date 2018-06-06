from os.path import isfile
import argparse


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "find the systematic variations of the yield in a channel for given process",
    epilog = """Example:
python yield_variations.py mc_v25v26pR5_updowns_with_data.root -c ctr_old_mu_sel_lj -p tt_taultauh,tt_lj,tt_other,tt_mutau,tt_taulj -s FSRUp,FSRDown,ISRUp,ISRDown,HDAMPUp,HDAMPDown,TuneCUETP8M2T4Up,TuneCUETP8M2T4Down"""
    )

parser.add_argument("data_file",    help="full Data+MC file name")
parser.add_argument("-c", "--channel",   type=str, help="the channel to look in")
parser.add_argument("-p", "--processes", type=str, default='tt_taultauh,tt_lj,tt_other,tt_mutau,tt_taulj', help="processes to consider")
parser.add_argument("-d", "--distr",     type=str, default='yield', help="the distribution (yield by default)")
parser.add_argument("-s", "--sys",       type=str, default='FSRUp,FSRDown,ISRUp,ISRDown,HDAMPUp,HDAMPDown,TuneCUETP8M2T4Up,TuneCUETP8M2T4Down', help="systematics to consider")
parser.add_argument("--debug",     action = "store_true", help="log debug output")


args = parser.parse_args()

import logging

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

logging.info("import ROOT")

import ROOT
from ROOT import TFile, THStack, TLegend, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
import plotting_root

assert isfile(args.data_file)

fdata = TFile(args.data_file)

for process in args.processes.split(','):
    nominal = fdata.Get("{chan}/{proc}/NOMINAL/{chan}_{proc}_NOMINAL_{distr}".format(chan=args.channel, proc=process, distr=args.distr))
    print process
    for sys in args.sys.split(','):
        systematic = fdata.Get("{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}".format(chan=args.channel, proc=process, sys=sys, distr=args.distr))
        print sys, systematic.Integral() / nominal.Integral()

