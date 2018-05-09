import argparse
import logging
from os.path import isfile
from sys import exit

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "straign in each bin per data uncertainty",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument("datafile", type=str, help="file with distrs")
parser.add_argument("-c", "--channel", type=str, default='mu_sel', help="channel")
parser.add_argument("-p", "--process", type=str, default='tt_lj', help="process")
parser.add_argument('-s', '--systematic',  type=str, default='NOMINAL', help="systematic")
parser.add_argument('-d', '--distr',  type=str, default='Mt_lep_met_f', help="distr")
#parser.add_argument('--logy',  action='store_true', help="log Y")


args = parser.parse_args()

assert isfile(args.datafile)

logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, gPad, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TLine, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
gROOT.SetBatch()
gStyle.SetOptStat(0)

tfile = TFile(args.datafile)
# ctr_old_el_sel_lj->Get("sums_NOMINAL/mc_sum1_NOMINAL")
histo_mc_sum = tfile.Get('{chan}/sums_{sys}/mc_sum1_{sys}'.format(chan=args.channel, sys=args.systematic))
# ctr_old_el_sel_lj->Get("data/NOMINAL/ctr_old_el_sel_lj_data_NOMINAL_Mt_lep_met_f")->ls()
histo_data   = tfile.Get('{chan}/data/NOMINAL/{chan}_data_NOMINAL_{distr}'.format(chan=args.channel, distr=args.distr))

## don't
#set mc errors = 0
for bini in range(histo_mc_sum.GetSize()):
    histo_mc_sum.SetBinError(bini, 0.)


#from plotting_root import rgb, nick_colour

cst = TCanvas("cst","stacked hists",10,10,700,700)

pad1 = TPad("pad1","This is pad1", 0., 0.,  1., 1.)

pad1.Draw()

pad1.cd()

histo_data.Divide(histo_mc_sum)
histo_mc_sum.Divide(histo_mc_sum)

histo_mc_sum.Draw()
histo_data.Draw("same")

cst.SaveAs('./prefit-bin-significance_{chan}_{sys}_{distr}.png'.format(chan=args.channel, sys=args.systematic, distr=args.distr))

# calc chi2 in the normalized data
chi2 = 0
nbins = 0
for bini in range(histo_data.GetSize()):
    content = histo_data.GetBinContent (bini)
    error   = histo_data.GetBinError   (bini)
    if error == 0:
        print "0 error bin", bini, content
        continue
    supposed_error = content**0.5
    print ((error - supposed_error)/supposed_error)**2
    chi2 += ((content - 1.) / error)**2
    nbins += 1

print nbins, chi2, chi2 / nbins, (chi2/nbins - 1) * nbins**0.5

