import argparse
import logging
from os.path import isfile
from sys import exit

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "to compare shapes from different files",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument("-c", "--channel", type=str, default='mu_sel', help="channel")
parser.add_argument("-p", "--process", type=str, default='tt_lj', help="process")
parser.add_argument('-s', '--systematic',  type=str, default='NOMINAL', help="systematic")
parser.add_argument('-d', '--distr',  type=str, default='Mt_lep_met_f', help="distr")

parser.add_argument('input_files', nargs='+', help="""the files process, `histo_name:filename`""")

args = parser.parse_args()

logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, gPad, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TLine, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
gROOT.SetBatch()
gStyle.SetOptStat(0)

from plotting_root import rgb, nick_colour

leg = TLegend(0.6,0.7, 0.9,0.9)

histos = []

for i, (nick, filename) in enumerate([fileparameter.split(':') for fileparameter in args.input_files]):
    if not isfile(filename):
        logging.info("missing: " + filename)
        continue

    tfile = TFile(filename)
    histo = tfile.Get("{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}".format(chan=args.channel, proc=args.process, sys=args.systematic, distr=args.distr))

    logging.debug("histo integral: %f" % histo.Integral())
    histo.SetDirectory(0)
    histo.SetLineColor(1 + i)
    histo.Scale(1./histo.Integral())
    histo.SetName(nick)
    histos.append(histo)
    leg.AddEntry(histo, nick, "l")


cst = TCanvas("cst","stacked hists",10,10,700,700)

pad1 = TPad("pad1","This is pad1", 0., 0.,  1., 1.)

pad1.Draw()

pad1.SetLogy()

pad1.cd()

histos[0].Draw()
for h in histos[1:]:
    h.Draw("same")

leg.Draw("same")

cst.SaveAs('./compare_%s_%s_%s_%s.png' % (args.channel, args.process, args.systematic, args.distr))

