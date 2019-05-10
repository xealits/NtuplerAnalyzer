import argparse
import logging
from os.path import isfile



parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "plot the pulls in the file",
    epilog = "Example:\npython pulls_plot.py pulls_of_nuisances_both_v37_test13_bunchFULLFIT2_noMCstat_shapeshort.root"
    )

parser.add_argument("pulls_file", help="the file with pulls")

parser.add_argument("--label-width", type=float, default=4.5, help="the width of the channel label")
parser.add_argument("--debug", action='store_true', help="DEBUG logging")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

assert isfile(args.pulls_file)

import ROOT
from ROOT import TFile, TPaveText, gStyle, gROOT

tfile = TFile(args.pulls_file)

gStyle.SetOptStat(0)
gROOT.SetBatch()

#c1 = TCanvas("c1","Pull Shorter Version",200,10,800,800);
##c1.SetGrid();
#c1.DrawFrame(0.85,0.75,1.15,1.25);


nuisancs = tfile.Get('nuisancs')
nuisancs.Draw()

pt = TPaveText(.5, 1.7, 0.5 + args.label_width, 2.7)
if   'pulls_of_nuisances_both' in args.pulls_file:
    pt.AddText("both channels")
elif 'pulls_of_nuisances_mu' in args.pulls_file:
    pt.AddText("muon-tau channel")
else:
    pt.AddText("electron-tau channel")

pt.Draw("same")

nuisancs.SaveAs(args.pulls_file.replace('.root', '.png'))


