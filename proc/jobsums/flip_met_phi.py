import argparse
import logging
from os.path import isfile
from sys import exit

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "to compare met phi in MC and Data flip them",
    epilog = 'Example:\npython flip_met_phi.py mc_v25v26pR5_updowns_HISTOSEL_met_phi_ctr_old_mu_sel_NOMINAL.root '
    )

parser.add_argument('input_file',  type=str, help="inp")
parser.add_argument('output_file', type=str, help="out")

parser.add_argument("-c", "--channel", type=str, default='ctr_old_mu_sel', help="channel")
parser.add_argument('-s', '--systematic',  type=str, default='NOMINAL', help="systematic")
parser.add_argument('-d', '--distr',  type=str, default='met_phi', help="distr")

args = parser.parse_args()

assert isfile(args.input_file)

logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, gPad, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TLine, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
gROOT.SetBatch()
gStyle.SetOptStat(0)

input_file = TFile(args.input_file)

data = input_file.Get("{chan}/data/{sys}/{chan}_data_{sys}_{distr}".format(sys=args.systematic, chan=args.channel, distr=args.distr))
mc   = input_file.Get("{chan}/sums_NOMINAL/mc_sum2_NOMINAL".format(chan=args.channel))

flipped_mc = mc.Clone()
flipped_mc.SetName("flipped_mc")
flipped_mc.SetDirectory(0)
flipped_mc.Reset()

flipped_mc_2 = mc.Clone()
flipped_mc_2.SetName("flipped_mc_2")
flipped_mc_2.SetDirectory(0)
flipped_mc_2.Reset()

flipped_mc_3 = mc.Clone()
flipped_mc_3.SetName("flipped_mc_3")
flipped_mc_3.SetDirectory(0)
flipped_mc_3.Reset()

flipped_data = data.Clone()
flipped_data.SetName("flipped_data")
flipped_data.SetDirectory(0)
flipped_data.Reset()

zero_bin   = mc.FindBin(0)
pi_bin_pos = mc.FindBin(3.14)
pi_bin_neg = mc.FindBin(-3.14)

# the flipping part

for i, bini in enumerate(range(pi_bin_neg, zero_bin)):
    content = mc.GetBinContent(bini)
    error   = mc.GetBinError(bini)
    flipped_bin = zero_bin + i
    flipped_mc.SetBinContent (flipped_bin, content)
    flipped_mc.SetBinError   (flipped_bin, error)
    logging.debug("flipped %2d %2d %4f" % (bini, flipped_bin, content))

for i, bini in enumerate(range(zero_bin, pi_bin_pos+1)):
    content = mc.GetBinContent(bini)
    error   = mc.GetBinError(bini)
    flipped_bin = pi_bin_neg+i
    flipped_mc.SetBinContent (flipped_bin, content)
    flipped_mc.SetBinError   (flipped_bin, error)
    logging.debug("flipped %2d %2d %4f" % (bini, flipped_bin, content))


#for i, bini in enumerate(range(pi_bin_neg, zero_bin)):
for i, bini in enumerate(range(pi_bin_pos, zero_bin, -1)):
    content = mc.GetBinContent(bini)
    error   = mc.GetBinError(bini)
    flipped_bin = zero_bin + i
    flipped_mc_2.SetBinContent (flipped_bin, content)
    flipped_mc_2.SetBinError   (flipped_bin, error)
    logging.debug("flipped %2d %2d %4f" % (bini, flipped_bin, content))

for i, bini in enumerate(range(zero_bin, pi_bin_neg-1, -1)):
    content = mc.GetBinContent(bini)
    error   = mc.GetBinError(bini)
    flipped_bin = pi_bin_neg+i
    flipped_mc_2.SetBinContent (flipped_bin, content)
    flipped_mc_2.SetBinError   (flipped_bin, error)
    logging.debug("flipped %2d %2d %4f" % (bini, flipped_bin, content))


for i, bini in enumerate(range(pi_bin_pos, pi_bin_neg-1, -1)):
    content = mc.GetBinContent(bini)
    error   = mc.GetBinError  (bini)
    flipped_bin = pi_bin_neg + i
    flipped_mc_3.SetBinContent (flipped_bin, content)
    flipped_mc_3.SetBinError   (flipped_bin, error)
    logging.debug("flipped %2d %2d %4f" % (bini, flipped_bin, content))


for i, bini in enumerate(range(pi_bin_pos, pi_bin_neg-1, -1)):
    content = data.GetBinContent(bini)
    error   = data.GetBinError  (bini)
    flipped_bin = pi_bin_neg + i
    flipped_data.SetBinContent (flipped_bin, content)
    flipped_data.SetBinError   (flipped_bin, error)
    logging.debug("flipped %2d %2d %4f" % (bini, flipped_bin, content))

fout = TFile(args.output_file, 'CREATE')
fout.Write()
flipped_mc.Write()
flipped_mc_2.Write()
flipped_mc_3.Write()
flipped_data.Write()
mc.Write()
data.Write()

fout.Write()
fout.Close()

