import argparse
import logging
from os.path import isfile, basename

#root [7] TH2D* lj = (TH2D*) ctr_old_mu_sel->Get("tt_lj/NOMINAL/ctr_old_mu_sel_tt_lj_NOMINAL_tau_dalitzes")




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "sumup 1 hist in root files",
    epilog = """Example:\npython printout_tau_dalitzes.py file_w_merged_histos.root"""
    )

parser.add_argument('histo_file', type=str, help='file with merged proc histos')

parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")


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

from plotting_root import rgb

gROOT.SetBatch()



tfile  = TFile(args.histo_file)

lj    = tfile.Get("ctr_old_mu_sel/tt_lj/NOMINAL/ctr_old_mu_sel_tt_lj_NOMINAL_tau_dalitzes")
mutau = tfile.Get("ctr_old_mu_sel/tt_mutau/NOMINAL/ctr_old_mu_sel_tt_mutau_NOMINAL_tau_dalitzes")

print "x,y,mutau,lj"
for nbin_x in range(mutau.GetNbinsX()):
    for nbin_y in range(mutau.GetNbinsY()):
        x = mutau.GetXaxis().GetBinCenter( nbin_x )
        y = mutau.GetYaxis().GetBinCenter( nbin_y )
        print x, y, mutau.GetBinContent(nbin_x, nbin_y), lj.GetBinContent(nbin_x, nbin_y)

