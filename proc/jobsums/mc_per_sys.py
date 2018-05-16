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

parser.add_argument("-o", "--output-directory", type=str, default='./', help="optional output directory")

parser.add_argument("--ratio-range", type=float, default=0.5, help="range of ratio plot (1-range 1+range)")
parser.add_argument("--y-max",       type=float,              help="set the maximum on Y axis")


args = parser.parse_args()

assert isfile(args.data_file)

if args.systematic:
    sys_names = args.systematic.split(',') 
else:
    exit(1)

logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
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

#leg.AddEntry(nom_MC, "MC NOMINAL", "L")
leg.AddEntry(nom_MC, "MC NOMINAL", "f")

# and data if not specified otherwise
if not args.no_data:
    data = fdata.Get(args.channel + '/data/NOMINAL/%s_data_NOMINAL_%s' % (args.channel, args.distr))
    data.Print()
    leg.AddEntry(data, "Data", "lep")

sys_color = {'TOPPTUp': kRed,
    'PUUp':   kGreen,
    'PUDown': kGreen + 1,
    'bSFUp':   kGreen,
    'bSFDown': kGreen + 1,
    'JERUp':   kGreen,
    'JERDown': kGreen + 1,
    'JESUp':   kGreen,
    'JESDown': kGreen + 1,
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
    sys_MC.SetLineColor(sys_color[sys])
    leg.AddEntry(sys_MC, "MC %s" % sys, "L")
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

if args.ratio:
    pad2.cd()

    nom_MC_rel = nom_MC.Clone()
    if not args.no_data:
        data_rel = data.Clone()

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
        sys_MCs_rel.append(h_rel)

    nom_MC_rel.Divide(nom_MC)
    nom_MC_rel.Draw("e2")
    for h_rel in sys_MCs_rel:
        h_rel.Draw("same")

    if not args.no_data:
        data_rel.Divide(nom_MC)
        data_rel.Draw("same")

cst.SaveAs(args.output_directory + '/MC-systs_%s_%s_%s_%s.png' % (args.data_file, args.channel, args.distr, args.systematic))

