import argparse
import logging
from os.path import isfile
from sys import exit

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "plot data and different SYS of MC sums w.r. to NOMINAL MC",
    epilog = "Example:\npython mc_per_sys.py histosel-roots/v25p2_mu_noqcd/histosels_v25p2_mu_ljbins.root -p tt_lj -c ctr_old_mu_sel_lj -d Mt_lep_met_f -s PUUp,PUDown"
    )

parser.add_argument("data_file", help="data file name")
parser.add_argument("-c", "--channel",    type=str, default='mu_sel', help="selection channels of events to consider (can be a list of channels for shape evolution study)")
parser.add_argument("-d", "--distr",      type=str, default='Mt_lep_met_f', help="distribution name")
parser.add_argument('-s', '--systematic', type=str, default=None, help="systematics to include (SYS1,SYS2,...)")
parser.add_argument('-p', '--process',    type=str, default=None, help="process to target, if given only the process is taken instead of sum of MC")

parser.add_argument('--no-data',     action='store_true', help="don't draw data")
parser.add_argument('-r', '--ratio', action='store_true', help="add ratio plot")

parser.add_argument('--add-titles', action='store_true', help="the left and right titles")
parser.add_argument("--title-x", type=str, help="set title of X axis on the plot")

parser.add_argument("-o", "--output-directory", type=str, default='./', help="optional output directory")

parser.add_argument("--ratio-range", type=float, default=0.3, help="range of ratio plot (1-range 1+range)")
parser.add_argument("--y-max",       type=float,              help="set the maximum on Y axis")


args = parser.parse_args()

assert isfile(args.data_file)

if args.systematic:
    sys_names = args.systematic.split(',') 
else:
    exit(1)

logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, kBlue
gROOT.SetBatch()
gStyle.SetOptStat(0)

from plotting_root import rgb

fdata = TFile(args.data_file)

leg = TLegend(0.6,0.6, 1.,1.)

# get nominal MC
sys = 'NOMINAL'
if args.process:
    nom_MC = fdata.Get('{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}'.format(chan=args.channel, proc=args.process, sys=sys, distr=args.distr))
else:
    nom_MC = fdata.Get(args.channel + '/sums_%s/mc_sum1_%s' % (sys, sys))
nom_MC.Print()

nom_MC.SetFillStyle(3004);
nom_MC.SetFillColor(1);
nom_MC.SetMarkerStyle(1)
nom_MC.SetMarkerColor(0)

nom_MC.GetXaxis().SetLabelFont(63)
nom_MC.GetXaxis().SetLabelSize(14) # labels will be 14 pixels
nom_MC.GetYaxis().SetLabelFont(63)
nom_MC.GetYaxis().SetLabelSize(14) # labels will be 14 pixels
nom_MC.GetXaxis().SetTitleFont(63)
nom_MC.GetXaxis().SetTitleSize(20)
nom_MC.GetYaxis().SetTitleFont(63)
nom_MC.GetYaxis().SetTitleSize(20)

if args.title_x and not args.ratio:
    title_x = args.title_x
    nom_MC.SetXTitle(title_x)
    nom_MC.GetXaxis().SetTitleOffset(2)

#leg.AddEntry(nom_MC, "MC NOMINAL", "L")
leg.AddEntry(nom_MC, "NOMINAL", "f")

# and data if not specified otherwise
if not args.no_data:
    data = fdata.Get(args.channel + '/data/NOMINAL/%s_data_NOMINAL_%s' % (args.channel, args.distr))
    data.Print()
    leg.AddEntry(data, "Data", "lep")

sys_color = {'TOPPTUp': kBlue,
    'FragUp':   kBlue,
    'FragDown': kGreen,
    'PUUp':   kBlue,
    'PUDown': kGreen,
    'bSFUp':   kBlue,
    'bSFDown': kGreen + 1,
    'JERUp':   kBlue,
    'JERDown': kGreen + 1,
    'JESUp':   kBlue,
    'JESDown': kGreen + 1,
    'TESUp':   kBlue,
    'TESDown': kGreen + 1,
    'ISRUp'   : kOrange,
    'ISRDown' : kOrange + 1,
    'FSRUp'   : kOrange,
    'FSRDown' : kOrange + 1,
    'TuneCUETP8M2T4Up'   : kOrange,
    'TuneCUETP8M2T4Down' : kOrange + 1,
}

sys_MCs = []
# and systematic MC-s
for sys in sys_names:
    if args.process:
        sys_MC = fdata.Get('{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}'.format(chan=args.channel, proc=args.process, sys=sys, distr=args.distr))
    else:
        sys_MC = fdata.Get(args.channel + '/sums_%s/mc_sum1_%s' % (sys, sys))
    sys_MC.Print()
    if sys in sys_color:
        sys_MC.SetLineColor(sys_color[sys])
    elif 'Up' in sys:
        sys_MC.SetLineColor(kRed)
    elif 'Down' in sys:
        sys_MC.SetLineColor(kRed+1)
    leg.AddEntry(sys_MC, "%s" % sys, "L")
    sys_MCs.append(sys_MC)


cst = TCanvas("cst","stacked hists",10,10,700,700)

if args.ratio:
    pad1 = TPad("pad1","This is pad1", 0., 0.3,  1., 1.)
    pad2 = TPad("pad2","This is pad2", 0., 0.05,  1., 0.3)
    pad1.Draw()
    pad2.Draw() # these have to be before cd()
    # now set margins:
    pad1.cd()
    ROOT.gPad.SetBottomMargin(0.02)
    #ROOT.gPad.SetTopMargin(0.01)
    pad2.cd()
    #ROOT.gPad.SetBottomMargin(0.001)
    ROOT.gPad.SetTopMargin(0.02)
else:
    pad1 = TPad("pad1","This is pad1", 0., 0.,  1., 1.)
    pad1.Draw()

pad1.cd()
if args.y_max:
    data.SetMaximum(args.y_max)

if args.ratio:
    #data.GetXaxis().SetLabelOffset(999)
    #data.GetXaxis().SetLabelSize(0)
    nom_MC.GetXaxis().SetLabelOffset(999)
    nom_MC.GetXaxis().SetLabelSize(0)

nom_MC.Draw("e2")
for h in sys_MCs:
    h.Print()
    h.Draw("same")

#nom_MC.Draw("same e2")
if not args.no_data:
    data.Draw("same")

leg.Draw("same")

if args.add_titles:
    logging.debug("adding left titles")
    left_title = TPaveText(0.12, 0.82, 0.28, 0.88, "brNDC")
    left_title.AddText("CMS")
    left_title.SetTextFont(1)
    left_title.SetFillColor(0)
    left_title .Draw("same")

    left_title2 = TPaveText(0.12, 0.76, 0.28, 0.82, "brNDC")
    left_title2.AddText("simulation")
    left_title2.SetTextFont(1)
    left_title2.SetFillColor(0)
    left_title2 .Draw("same")

if args.ratio:
    pad2.cd()

    nom_MC_rel = nom_MC.Clone()
    if not args.no_data:
        data_rel = data.Clone()

    nom_MC_rel.GetXaxis().SetLabelFont(63)
    nom_MC_rel.GetXaxis().SetLabelSize(14) # labels will be 14 pixels
    nom_MC_rel.GetYaxis().SetLabelFont(63)
    nom_MC_rel.GetYaxis().SetLabelSize(14) # labels will be 14 pixels
    nom_MC_rel.GetXaxis().SetTitleFont(63)
    nom_MC_rel.GetXaxis().SetTitleSize(20)
    nom_MC_rel.GetYaxis().SetTitleFont(63)
    nom_MC_rel.GetYaxis().SetTitleSize(20)

    offset_x_title = 3

    nom_MC_rel.GetXaxis().SetLabelOffset(0.01)
    nom_MC_rel.GetXaxis().SetTitleOffset(offset_x_title)
    #nom_MC_rel.GetYaxis().SetTitleOffset(1.4)
    #nom_MC_rel.GetXaxis().SetLabelSize(14)

    if args.title_x:
        title_x = args.title_x
        nom_MC_rel.SetXTitle(title_x)
        #nom_MC_rel.GetXaxis().SetTitleOffset(2.)

    ratio_max = 1. + args.ratio_range
    ratio_min = 1. - args.ratio_range
    #data_rel.SetMaximum(ratio_max)
    #data_rel.SetMinimum(ratio_min)
    nom_MC_rel.SetMaximum(ratio_max)
    nom_MC_rel.SetMinimum(ratio_min)

    sys_MCs_rel = []
    for h in sys_MCs:
        h_rel = h.Clone()
        h_rel.SetName('rel_'+h.GetName())
        h_rel.SetDirectory(0)
        h_rel.Print()
        h_rel.Divide(nom_MC)

        #if args.title_x:
        #    title_x = args.title_x
        #    h_rel.SetXTitle(title_x)

        sys_MCs_rel.append(h_rel)

    nom_MC_rel.Divide(nom_MC)
    nom_MC_rel.Draw("H")
    nom_MC_rel.Draw("e2")
    for h_rel in sys_MCs_rel:
        h_rel.GetXaxis().SetTitleOffset(offset_x_title)
        h_rel.Draw("same")

    if not args.no_data:
        data_rel.Divide(nom_MC)
        data_rel.GetXaxis().SetTitleOffset(offset_x_title)
        data_rel.Draw("same")

    #right_title = TPaveText(0.65, 0.9, 0.9,  0.95, "brNDC")
    #right_title.AddText("%s fb^{-1} (13 TeV)" % (args.lumi / 1000. if args.lumi else args.lumi_label))
    #right_title.SetTextFont(132)
    #right_title.SetFillColor(0)
    #right_title.Draw("same")


cst.SaveAs(args.output_directory + '/MC-systs_%s_%s_%s_%s_%s.png' % (args.data_file, args.channel, args.process, args.distr, args.systematic))

