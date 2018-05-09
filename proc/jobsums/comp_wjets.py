import argparse
import logging
from os.path import isfile
from sys import exit

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "wjets comp files",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument('--normalize',     action='store_true', help="normalize integrals to 1")

args = parser.parse_args()

'''
parser.add_argument("data_file", help="data file name")
parser.add_argument("-c", "--channel", type=str, default='mu_sel', help="selection channels of events to consider (can be a list of channels for shape evolution study)")
parser.add_argument("-d", "--distr", type=str, default='Mt_lep_met_f', help="distribution name")
parser.add_argument('-s', '--systematic',  type=str, default=None, help="systematics to include (SYS1,SYS2,...)")

parser.add_argument('--no-data',     action='store_true', help="don't draw data")
parser.add_argument('-r', '--ratio', action='store_true', help="add ratio plot")

parser.add_argument("-o", "--output-directory", type=str, default='./', help="optional output directory")

parser.add_argument("--ratio-range", type=float, default=0.5, help="range of ratio plot (1-range 1+range)")
parser.add_argument("--y-max",       type=float,              help="set the maximum on Y axis")


assert isfile(args.data_file)

if args.systematic:
    sys_names = args.systematic.split(',') 
else:
    exit(1)
'''

logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
gROOT.SetBatch()
gStyle.SetOptStat(0)

from plotting_root import rgb, nick_colour


f_am = TFile('v20v21p1_mu_mc_aMCatNLO-wjets.root')
f_mg = TFile('v20v21p1_mu_mc.root')

distr = "Mt_lep_met_f"

wjets_am = f_am.Get("ctr_mu_wjet/wjets/NOMINAL/ctr_mu_wjet_wjets_NOMINAL_%s" % distr)
wjets_mg = f_mg.Get("ctr_mu_wjet/wjets/NOMINAL/ctr_mu_wjet_wjets_NOMINAL_%s" % distr)

wjets_am.SetName("WJets aMCatNLO")
wjets_mg.SetName("WJets Madgraph")

wjets_mg.Scale(0.72)
wjets_am.Scale()

#wjets_mg.Scale(0.80)
#wjets_am.Scale(1.15)

normalize = args.normalize

if normalize:
    wjets_mg.Scale(1./wjets_mg.Integral())
    wjets_am.Scale(1./wjets_am.Integral())

wjets_mg.SetLineColor(nick_colour['wjets'])
wjets_am.SetLineColor(kRed)

wjets_mg.SetLineWidth(2)
wjets_am.SetLineWidth(2)



wjets_mg.SetTitle('normalized' if normalize else '')
wjets_am.SetTitle('normalized' if normalize else '')

wjets_mg.SetXTitle(distr)
wjets_am.SetXTitle(distr)

wjets_mg.SetYTitle('Evts.')
wjets_am.SetYTitle('Evts.')

wjets_mg.GetYaxis().SetTitleFont(63)
wjets_am.GetYaxis().SetTitleSize(20)

wjets_mg.GetYaxis().SetLabelSize(0.02)
wjets_am.GetXaxis().SetLabelSize(0.02)

wjets_mg.GetYaxis().SetTitleOffset(1.5) # place the title not overlapping with labels...
wjets_am.GetYaxis().SetTitleOffset(1.5)


leg = TLegend(0.6,0.7, 0.9,0.9)

leg.AddEntry(wjets_mg, "WJets Madgraph", "l")
leg.AddEntry(wjets_am, "WJets aMCatNLO", "l")

cst = TCanvas("cst","stacked hists",10,10,700,700)

wjets_mg.Draw()
wjets_am.Draw('same')

leg.Draw("same")

output_directory = './'
cst.SaveAs(output_directory + '/wjets-comp_%s%s.png' % (distr, '_normalized' if normalize else ''))



