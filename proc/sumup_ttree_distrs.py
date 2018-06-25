import argparse
import logging
from os.path import isfile, basename




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "sumup TTree Draw",
    epilog = """Example:\npython sumup_ttree_distrs.py "met_init.pt()" --histo-range 200,0,200 --histo-name data_met_init --output data_met_init.root gstore_outdirs/v28/SingleMuon/Ntupler_v28_Data13TeV_SingleMuon2016*/1*/*/*root"""
    )

parser.add_argument('draw_com', type=str, help='Draw("??")')
parser.add_argument('--cond-com', type=str, default="", help='Draw("", "??")')
parser.add_argument('--ttree',    type=str, default="ntupler/reduced_ttree", help='path to the TTree in the file')
parser.add_argument('--histo-name',  type=str, default="out_histo", help='the ROOTName for output')
parser.add_argument('--histo-range', type=str, default=None, help='optionally set the range')
parser.add_argument('--histo-color', type=str, default=None, help='optional rgb color, like `255,255,255`')

parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")
parser.add_argument("--output", type=str, default="output.root", help="filename for output")

parser.add_argument("--per-weight",  action='store_true', help="normalize by event weight of datasets")
parser.add_argument("--scan",        action='store_true', help="also scan events ant print out")
parser.add_argument("--get-maximum", action='store_true', help="just find maxima on the draw_com")

parser.add_argument('input_files', nargs='+', help="""the files to sum up, passed by shell, as:

/gstore/t3cms/store/user/otoldaie/v19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v19_MC2016_Summer16_TTJets_powheg/180226_022336/0000/MC2016_Summer16_TTJets_powheg_*.root""")


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



input_files = args.input_files


'''
out_histo  = TH1D("os_tauall", "", 44, -2, 20)
tau_b      = TH1D("os_taub", "", 44, -2, 20)
tau_w      = TH1D("os_tauw", "", 44, -2, 20)
tau_w_c    = TH1D("os_tauw_c", "", 44, -2, 20)
tau_w_notc = TH1D("os_tauw_notc", "", 44, -2, 20)

ss_tau_all    = TH1D("ss_tauall", "", 44, -2, 20)
ss_tau_b      = TH1D("ss_taub", "", 44, -2, 20)
ss_tau_w      = TH1D("ss_tauw", "", 44, -2, 20)
ss_tau_w_c    = TH1D("ss_tauw_c", "", 44, -2, 20)
ss_tau_w_notc = TH1D("ss_tauw_notc", "", 44, -2, 20)
'''


if args.histo_color:
    logging.debug("color: " + args.histo_color)
else:
    logging.debug("no color")

# Draw command
# 
if args.histo_range:
    draw_command = args.draw_com + ">>h(%s)" % args.histo_range
else:
    # TOFIX: without explicit range the histos won't add up
    draw_command = args.draw_com

logging.debug("draw: " + draw_command)
logging.debug("cond: " + args.cond_com)

out_histo = None
weight_counter = None

if args.get_maximum:
   maximum = -1111111.

for filename in input_files:
    if not isfile(filename):
        logging.info("missing: " + filename)
        continue

    logging.debug(filename)

    tfile = TFile(filename)
    ttree = tfile.Get(args.ttree)

    # Draw the file and sum up
    # 
    # TOFIX: without explicit range the histos won't add up
    if args.get_maximum:
        m = ttree.GetMaximum(args.draw_com)
        print "%30s %f" % (basename(filename), m)
        maximum = max(m, maximum)
    else:
        ttree.Draw(draw_command,  args.cond_com)
        histo = ttree.GetHistogram()
        histo.SetDirectory(0)

        if args.per_weight:
            wcounter = tfile.Get('ntupler/weight_counter')
            if not weight_counter:
                weight_counter = wcounter
                weight_counter.SetDirectory(0)
                weight_counter.SetName("sumup_weight_counter")
            else:
                weight_counter.Add(wcounter)

        if not out_histo:
            out_histo = histo
            out_histo.SetName(args.histo_name)
        else:
            out_histo.Add(histo)

    if args.scan:
        print filename
        ttree.Scan(args.draw_com, args.cond_com)

    tfile.Close()

if args.per_weight:
    #weight_counter = tfile.Get('ntupler/weight_counter')
    histo.Scale(1./weight_counter.GetBinContent(2))

if args.histo_color:
    out_histo.SetLineColor(rgb(*(int(i) for i in args.histo_color.split(','))))


if args.get_maximum:
   print "Max:", maximum

''' no legend for now
leg = TLegend(0.5, 0.5, 0.9, 0.9)
leg.SetName("legOS")

leg.AddEntry(tau_all    , "tt->lj, OS, all", 'L')
leg.AddEntry(tau_b      , "tt->lj, OS, b prod", 'L')
leg.AddEntry(tau_w      , "tt->lj, OS, w prod", 'L')
leg.AddEntry(tau_w_c    , "tt->lj, OS, w, C flav", 'L')
leg.AddEntry(tau_w_notc , "tt->lj, OS, w, not C", 'L')

leg2 = TLegend(0.5, 0.5, 0.9, 0.9)
leg2.SetName("legSS")

leg2.AddEntry(tau_all    , "tt->lj, SS, all", 'L')
leg2.AddEntry(tau_b      , "tt->lj, SS, b prod", 'L')
leg2.AddEntry(tau_w      , "tt->lj, SS, w prod", 'L')
leg2.AddEntry(tau_w_c    , "tt->lj, SS, w, C flav", 'L')
leg2.AddEntry(tau_w_notc , "tt->lj, SS, w, not C", 'L')
'''


fout = TFile(args.output, "RECREATE")
fout.Write()

out_histo.Write()

fout.Write()
fout.Close()

