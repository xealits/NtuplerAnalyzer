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
parser.add_argument('--logy',  action='store_true', help="log Y")

parser.add_argument('--formula',  type=str, help="to plot h1 overlayed with 123*h2+982*h3: `h1 : 123 h2 , 982h3`")

parser.add_argument('input_files', nargs='+', help="""the files process, `histo_name:filename`""")

args = parser.parse_args()

logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, gPad, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TLine, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
gROOT.SetBatch()
gStyle.SetOptStat(0)

from plotting_root import rgb, nick_colour

leg = TLegend(0.6,0.7, 0.9,0.9)

histos = {}

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
    histos[nick] = histo


cst = TCanvas("cst","stacked hists",10,10,700,700)

pad1 = TPad("pad1","This is pad1", 0., 0.,  1., 1.)

pad1.Draw()

if args.logy:
    pad1.SetLogy()

pad1.cd()

if not args.formula:
    for i, (nick, h) in enumerate(histos.items()):
        if i == 0:
            h.Draw()
        else:
            h.Draw("same")
        leg.AddEntry(histo, nick, "l")

else:
    # parse the formula
    import re
    drawn = False
    for subform in args.formula.split(':'): # sub-formulas
        # ' number  nick ' pairs
        logging.debug(subform)
        nicks = [nick for nick in histos.keys() if nick in subform]
        logging.debug(nicks)
        nick_pairs = []
        for nick_pair in subform.split(','):
            for nick in nicks:
                if nick in nick_pair:
                    logging.debug(nick_pair)
                    number = nick_pair.replace(nick, '')
                    nick_pairs.append((float(number) if number.strip() else 1., nick))
                    break

        logging.debug(nick_pairs)

        # construct the histo of the subform
        factor, nick = nick_pairs[0]
        histo = histos[nick]
        histo.Scale(factor)
        for factor, nick in nick_pairs[1:]:
            h = histos[nick]
            h.Scale(factor)
            histo.Add(h)

        # normalize the histo
        histo.Scale(1./histo.Integral())
        # and draw
        if drawn:
            histo.Draw("same")
        else:
            histo.Draw()
            drawn = True

        leg.AddEntry(histo, subform, "l")

leg.Draw("same")

cst.SaveAs('./compare_%s_%s_%s_%s_%s%s.png' % (args.channel, args.process, args.systematic, args.distr, ','.join(histos.keys()), '_formula' if args.formula else ''))

