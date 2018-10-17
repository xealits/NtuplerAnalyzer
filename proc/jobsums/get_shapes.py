import argparse
import logging
from os.path import isfile
from sys import exit

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "to compare shapes from different files",
    epilog = 'Example:\npython jobsums/get_shapes.py --formula "h1 : h2 " -c ctr_old_mu_sel_lj -d Mt_lep_met_f -p tt_lj h1:jobsums/AN-v1_mu_file-orig.root h2:merge-sets/v23/p10/MC2016_Summer16_TTJets_powheg.root'
    )

parser.add_argument("-c", "--channel", type=str, default='mu_sel', help="channel")
parser.add_argument("-p", "--process", type=str, default='tt_lj', help="process")
parser.add_argument('-s', '--systematic',  type=str, default='NOMINAL', help="systematic")
parser.add_argument('-d', '--distr',  type=str, default='Mt_lep_met_f', help="distr")
parser.add_argument('--logy',     action='store_true', help="log Y")
parser.add_argument('--no-norm',  action='store_true', help="don't normalize each histo to 1")
parser.add_argument('--norm-binwidth',  action='store_true', help="normalize each bin to width")
parser.add_argument('--norm-formulas',     action='store_true', help="normalize the formulas of histos to 1")
parser.add_argument("--y-range",     type=str,      help="set Y range as `ymin,ymax`")
parser.add_argument("--x-title",     type=str,      help="title of X axis")

parser.add_argument('--formula',  type=str, help="to plot h1 overlayed with 123*h2+982*h3: `h1 : 123 h2 , 982h3`")

parser.add_argument('--set-file',  type=str, help=" histos come from this file, skip filename in the input_files definition, but : must be there")

parser.add_argument('input_files', nargs='+', help="""the files process, `histo_name:filename[:chan,proc,sys[,distr]]`""")

args = parser.parse_args()

if args.y_range:
    y_min, y_max = [float(x) for x in args.y_range.split(',')]
else:
    y_min = y_max = None


logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, gPad, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TLine, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
gROOT.SetBatch()
gStyle.SetOptStat(0)

from plotting_root import rgb, nick_colour

leg = TLegend(0.6,0.7, 0.9,0.9)

histos = {}

for i, fileparameter in enumerate(args.input_files):
    pars = fileparameter.split(':')
    assert len(pars) == 2 or len(pars) == 3
    if len(pars) == 3:
        nick, filename, opts = pars
    else:
        nick, filename = pars
        opts = None

    if args.set_file:
        filename = args.set_file

    if not isfile(filename):
        logging.info("missing: " + filename)
        continue
    logging.debug("%s" % filename)

    channel = process = systematic = distr = ''
    if opts and len(opts.split(',')) == 4:
        channel, process, systematic, distr = opts.split(',')
    elif opts and len(opts.split(',')) == 3:
        channel, process, systematic = opts.split(',')
    elif opts and len(opts.split(',')) == 2:
        channel, process = opts.split(',')

    if not channel:
        channel    = args.channel
    if not process:
        process    = args.process
    if not systematic:
        systematic = args.systematic
    if not distr:
        distr = args.distr

    histo_path = "{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}".format(chan=channel, proc=process, sys=systematic, distr=distr)

    tfile = TFile(filename)
    logging.debug("%s" % histo_path)
    histo = tfile.Get(histo_path)

    histo.Sumw2()

    logging.debug("histo    : %20s" % nick)
    logging.debug("integral : %f" % histo.Integral())
    histo.SetDirectory(0)
    histo.SetLineColor(1 + i)
    if not args.no_norm:
        histo.Scale(1./histo.Integral())
    histo.SetName(nick)
    if y_max is not None and y_min is not None:
        histo.SetMaximum(y_max)
        histo.SetMinimum(y_min)
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
            if y_max is not None and y_min is not None:
                h.SetMaximum(y_max)
                h.SetMinimum(y_min)
            if args.x_title:
                #nom_MC.SetTitle(args.process)
                h.SetXTitle(args.x_title)
            h.Draw()
        else:
            h.Draw("same")
        leg.AddEntry(histo, nick, "l")

else:
    # parse the formula
    import re
    drawn = False
    form_histos = []

    for form_i, subform in enumerate(args.formula.split(':')): # sub-formulas
        # ' number  nick ' pairs
        logging.debug("subform " + subform)
        nicks = [nick for nick in histos.keys() if nick in subform]
        logging.debug("nicks %s" % str(nicks))
        nick_pairs = []
        for nick_pair in subform.split('+'):
            for nick in nicks:
                if nick in nick_pair:
                    logging.debug("nick_pair " + nick_pair)
                    number = nick_pair.replace(nick, '')
                    logging.debug("number " + number)
                    factor = float(number) if number.strip() else 1.
                    logging.debug("factor %f" % factor)
                    nick_pairs.append((factor, nick))
                    break

        logging.debug(nick_pairs)

        ## normalize the histo
        #if not args.no_norm:
        #    histo.Scale(1./histo.Integral())

        # construct the histo of the subform
        factor, nick = nick_pairs[0]
        histo = histos[nick].Clone()
        histo.SetName("formula%d" % form_i)
        histo.SetDirectory(0)
        histo.Scale(factor)
        for factor, nick in nick_pairs[1:]:
            # not sure if pyroot can do h += h2*factor correctly, without mutation
            h = histos[nick].Clone()
            h.Scale(factor)
            histo.Add(h)

        if args.norm_formulas:
            histo.Scale(1./histo.Integral())

        logging.debug("histo %20s %f" % (nick, histo.Integral()))

        form_histos.append(histo)

        leg.AddEntry(histo, subform, "l")

    # and draw
    for histo in form_histos:
        # normalize bin width if needed
        if args.norm_binwidth:
            #
            for bini in range(histo.GetSize()):
                content = histo.GetBinContent(bini)
                error   = histo.GetBinError(bini)
                width   = histo.GetXaxis().GetBinUpEdge(bini) - histo.GetXaxis().GetBinLowEdge(bini)
                histo.SetBinContent(bini, content/width)
                histo.SetBinError(bini, error/width)

            # rescale to 1 just in case
            histo.Scale(1./histo.Integral())

        if drawn:
            if y_max is not None and y_min is not None:
                logging.debug("setting min-max")
                histo.SetMaximum(y_max)
                histo.SetMinimum(y_min)
            histo.Draw("same")
        else:
            if args.x_title:
                #nom_MC.SetTitle(args.process)
                histo.SetXTitle(args.x_title)
            if y_max is not None and y_min is not None:
                logging.debug("setting min-max")
                histo.SetMaximum(y_max)
                histo.SetMinimum(y_min)

            # draw the first one with line
            histo.SetFillColor(0)
            histo.Draw("hist e")
            drawn = True

leg.Draw("same")

cst.SaveAs('./compare_%s_%s_%s_%s_%s%s%s.png' % (args.channel, args.process, args.systematic, args.distr, ','.join(histos.keys()), '_formula' if args.formula else '', '_nonorm' if args.no_norm else ''))

