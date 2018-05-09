import argparse
import logging
from os.path import isfile
from sys import exit

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "plot data and different SYS of MC sums w.r. to NOMINAL MC",
    epilog = '''Example:
python distr_ratios.py mc_v25p2_v26p2.root -c ctr_old_mu_sel:Mt_lep_met_f -p tt_lj --ratio-channel-distrs ctr_old_mu_sel:Mt_lep_met_t2,ctr_old_mu_sel:Mt_lep_met_t3,ctr_old_mu_sel:Mt_lep_met_t4,ctr_old_mu_sel:Mt_lep_met_t5,ctr_old_mu_sel:Mt_lep_met_t6
python distr_ratios.py mc_v25v26pF5.root -c ctr_old_mu_presel:tt_lj:met_lep_phi --ratio-channel-distrs ctr_old_mu_presel:tt_mutau:met_lep_phi,ctr_old_mu_presel:tt_taultauh:met_lep_phi,ctr_old_mu_presel:tt_other:met_lep_phi'''
    )

parser.add_argument("data_file", help="data file name")
parser.add_argument("-c", "--channel-distr", type=str, default='mu_sel:Mt_lep_met_f', help="the reference channel:distr (mu_sel:Mt_lep_met_f)")
#parser.add_argument("-d", "--distr", type=str, default='Mt_lep_met_f', help="distribution name")
parser.add_argument("-p", "--process", type=str, default='tt_lj', help="process")
parser.add_argument('-s', '--systematic',  type=str, default='NOMINAL', help="systematics to include (SYS1,SYS2,...)")

parser.add_argument("--ratio-channel-distrs", type=str, default='mu_sel:Mt_lep_met_t2', help="distribution names")

parser.add_argument('--no-data',     action='store_true', help="don't draw data")
parser.add_argument('-r', '--ratio', action='store_true', help="add ratio plot")

parser.add_argument('--logy', action='store_true', help="log scale Y axis")

parser.add_argument('--leg-chan', action='store_true', help="put channels in the legend")

parser.add_argument("-o", "--output-directory", type=str, default='./', help="optional output directory")
parser.add_argument("--custom-filename", type=str, help="custom filename")

parser.add_argument("--ratio-range", type=float, default=0.5, help="range of ratio plot (1-range 1+range)")
#parser.add_argument("--y-max",       type=float,              help="set the maximum on Y axis")
parser.add_argument("--y-range",     type=str,              help="set Y range as `ymin,ymax`")

parser.add_argument("--title",   type=str, default=None, help="title for the plot")
parser.add_argument("--title-x", type=str, default=None, help="title for the plot")

parser.add_argument("--add-y-line", type=float, help="add horizontal line at the given Y value")

args = parser.parse_args()

assert isfile(args.data_file)

if args.title is not None:
    title = args.title
else:
    title = args.process

refs = args.channel_distr.split(':')
if len(refs) == 3:
    ref_chan, ref_proc, ref_distr = refs
else:
    # assume the process is the same
    ref_chan, ref_distr = refs
    ref_proc = args.process

if args.title_x is not None:
    title_x = args.title_x
else:
    title_x = ref_distr

logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, gPad, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TLine, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
gROOT.SetBatch()
gStyle.SetOptStat(0)

from plotting_root import rgb, nick_colour

fdata = TFile(args.data_file)

leg = TLegend(0.6,0.7, 0.9,0.9)

distr_nick = {
    'Mt_lep_met_f':  'inclusive',
    'Mt_lep_met_t1': 'pt < 30',
    'Mt_lep_met_t2': 'pt < 40',
    'Mt_lep_met_t3': 'pt < 50',
    'Mt_lep_met_t4': 'pt < 70',
    'Mt_lep_met_t5': 'pt < 90',
    'Mt_lep_met_t6': 'pt > 90',
}

# get nominal MC
sys = args.systematic

nom_MC = fdata.Get(ref_chan + '/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}'.format(proc=ref_proc, chan=ref_chan, sys=args.systematic, distr=ref_distr))
nom_MC.Print()

if args.leg_chan:
    leg_entry = ref_chan
else:
    leg_entry = distr_nick.get(ref_distr, ref_distr)
leg.AddEntry(nom_MC, leg_entry, "L")

#leg.AddEntry(nom_MC, "MC NOMINAL", "L")
#leg.AddEntry(nom_MC, "MC NOMINAL", "f")
sys = args.systematic

## and data if not specified otherwise
#if not args.no_data:
#    data = fdata.Get(args.channel + '/data/{sys}/{chan}_data_{sys}_{distr}'.format(chan=args.channel, distr=args.distr, sys=sys))
#    data.Print()
#    leg.AddEntry(data, "Data", "lep")

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

histo_MCs = []
# and systematic MC-s
for i, chan_def in enumerate(chd.split(':') for chd in args.ratio_channel_distrs.split(',')):
    if len(chan_def) == 3:
        channel, proc, distr = chan_def
    else:
        channel, distr = chan_def
        proc = args.process
    sys_MC = fdata.Get(channel + '/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}'.format(proc=proc, chan=channel, sys=args.systematic, distr=distr))
    sys_MC.Print()
    sys_MC.SetLineColor   (nick_colour[proc]  + 1*i)
    sys_MC.SetMarkerColor (nick_colour[proc]  + 1*i)
    sys_MC.Divide(nom_MC)
    sys_MC.Print()
    #sys_MC.SetFillColor(0)

    if args.leg_chan:
        leg_entry = channel
    else:
        leg_entry = distr_nick.get(distr, distr)
    leg.AddEntry(sys_MC, leg_entry, "L")
    histo_MCs.append(sys_MC)


cst = TCanvas("cst","stacked hists",10,10,700,700)

pad1 = TPad("pad1","This is pad1", 0., 0.,  1., 1.)

pad1.Draw()

pad1.cd()
#if args.y_max:
#    data.SetMaximum(args.y_max)

if args.ratio:
    data.GetXaxis().SetLabelOffset(999)
    data.GetXaxis().SetLabelSize(0)

print 'title', title

nom_MC.SetTitle(title)
nom_MC.SetXTitle(title_x)

nom_MC.Divide(nom_MC)
nom_MC.Print()

if args.y_range:
    y_min, y_max = [float(x) for x in args.y_range.split(',')]
else:
    y_max = 2
    if args.logy:
        y_min = 0.001
    else:
        y_min = 0.

nom_MC.SetMaximum(y_max)
nom_MC.SetMinimum(y_min)

if args.logy:
    #cst.SetLogy()
    #pad1.SetLogy()
    gPad.SetLogy()

nom_MC.SetTitle(title)
nom_MC.SetLineWidth(2)
nom_MC.Draw("e")

for histo in histo_MCs:
    histo.Print()
    histo.SetTitle(title)
    histo.SetLineWidth(2)
    histo.Draw("same e")

if args.add_y_line:
    x_min = nom_MC.GetXaxis().GetXmin()
    x_max = nom_MC.GetXaxis().GetXmax()
    l = TLine(x_min, args.add_y_line, x_max, args.add_y_line)
    l.Draw("same")

leg.Draw("same")

if args.custom_filename:
    filename = args.custom_filename + '.png'
else:
    filename = 'ratios_%s_%s_%s_%s_%s.png' % (ref_chan, args.process, args.systematic, ref_distr, args.ratio_channel_distrs)

cst.SaveAs(args.output_directory + '/' + filename)

