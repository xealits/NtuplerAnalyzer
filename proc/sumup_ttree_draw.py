import argparse
import logging
from os.path import isfile, basename
import ctypes




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "sumup TTree Draw",
    epilog = """Example:\npython sumup_ttree_draw.py "met_init.pt()" --histo-range 200,0,200 --histo-name data_met_init --output data_met_init.root gstore_outdirs/v28/SingleMuon/Ntupler_v28_Data13TeV_SingleMuon2016*/1*/*/*root
python sumup_ttree_draw.py "event_leptons[0].pt()" --ttree ttree_out --cond-com "selection_stage == 5" --histo-range 50,0,200 --output test1_lep_pt.root --histo-name foo/bar/test1_lep_pt --save-weight MC2016_Summer16_TTJets_powheg_test1.root """
    )

parser.add_argument('draw_com', type=str, help='Draw("??")')
parser.add_argument('--cond-com', type=str, default="", help='Draw("", "??")')
parser.add_argument('--ttree',    type=str, default="ntupler/reduced_ttree", help='path to the TTree in the file')
parser.add_argument('--histo-name',  type=str, default="out_histo", help='the ROOTName for output')
parser.add_argument('--std-histo-name',  type=str, default="distr", help='construct standard name for histogram: take "path" part from --histo-name and "distr" part from here, attach as chan_proc_sys_distr')
parser.add_argument('--histo-color', type=str, default=None, help='optional rgb color, like `255,255,255`')
parser.add_argument('--histo-range',  type=str, default=None, help='optionally set the range')
parser.add_argument('--custom-range', type=str, default=None, help='optionally set the range in custom bins')
parser.add_argument("--try-xsec",  action='store_true', help="try to normalize to cross-section, if the output filename contains dtag")

parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")
parser.add_argument("--output", type=str, default="output.root", help="filename for output")

parser.add_argument("--save-weight", action='store_true', help="save the weight counters in the output file")
parser.add_argument("--per-weight",  action='store_true', help="normalize by event weight of datasets")
parser.add_argument("--scale", type=float, help="the value to scale histo")
parser.add_argument("--scan",        action='store_true', help="also scan events ant print out")
parser.add_argument("--get-maximum", action='store_true', help="just find maxima on the draw_com")

parser.add_argument('input_files', nargs='+', help="""the files to sum up, passed by shell, as:

/gstore/t3cms/store/user/otoldaie/v19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v19_MC2016_Summer16_TTJets_powheg/180226_022336/0000/MC2016_Summer16_TTJets_powheg_*.root""")


args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

# [/]histo_path/histo_name
histo_name = args.histo_name.split('/')[-1]
histo_path = args.histo_name.split('/')[:-1]
# skip empty parts in the path (sequences of ////)
histo_path = [part for part in histo_path if part]

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

'''
root [4] double pt_bins[] = {20, 30, 40, 50, 70, 90, 120, 150, 200, 300, 2000};
root [5] 
root [5] sizeof
ROOT_prompt_5:2:1: error: expected expression
;
^
root [6] sizeof(pt_bins)
(unsigned long) 88
root [7] 
root [7] sizeof(double)
(unsigned long) 8
root [8] 
root [8] int pt_bins_n = 10
(int) 10
root [9] 
root [9] 
root [9] TH1D* hpts = new TH1D("hpts", "", pt_bins_n, pt_bins)
(TH1D *) 0x4db6750
root [10] hpts
(TH1D *) 0x4db6750
root [11] 
root [11] 
root [11] ttree_out->Draw("event_leptons[0].pt()>>hpts")
'''

if args.try_xsec:
    from per_dtag_script import dtags

# Draw command
# 
temp_output_histo = None # histo-template for custom bins
if args.custom_range:
    bin_edges = [float(b) for b in args.custom_range.split(',')]
    n_bin_edges = len(bin_edges)
    n_bins = n_bin_edges - 1
    logging.debug("making %d bins from %s" % (n_bins, args.custom_range))

    # first method
    root_bin_edges = (ctypes.c_double * n_bin_edges)(*bin_edges)
    temp_output_histo = TH1D("histotemp", "", n_bins, root_bin_edges) # root commands can access it by the name

    draw_command = args.draw_com + '>>' + "histotemp"
elif args.histo_range:
    draw_command = args.draw_com + ">>h(%s)" % args.histo_range
else:
    # TOFIX: without explicit range the histos won't add up
    draw_command = args.draw_com

logging.debug("draw: " + draw_command)
logging.debug("cond: " + args.cond_com)
if temp_output_histo:
    logging.debug("temp: " + temp_output_histo.GetName())
    #temp_output_histo.SetDirectory(0)
    # and some preparation for root stuff
    ROOT.gROOT.ProcessLine('TFile* intfile;')
    ROOT.gROOT.ProcessLine('TTree* inttree;')

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
    if temp_output_histo:
        temp_output_histo.SetDirectory(tfile)
        # in ROOT the TTree.Draw command "sees" only histograms in current working directory
        # after opening new file ROOT moves working directory into that file
        # therefore move temp tehre too
        # probably it moving current directory to some 0 directory could work
        # and as probable is the possibility to bring more troubles that way

    # Draw the file and sum up
    # 
    # TOFIX: without explicit range the histos won't add up, need to pick up the range of the first histo and propagate to the following?
    if args.get_maximum:
        m = ttree.GetMaximum(args.draw_com)
        print "%30s %f" % (basename(filename), m)
        maximum = max(m, maximum)
    else:
        ttree.Draw(draw_command,  args.cond_com)

        if temp_output_histo:
            # if there is temp histo then the histo was written there
            logging.debug(temp_output_histo.Integral())
            logging.debug(ROOT.histotemp.Integral())
            histo = temp_output_histo
        else:
            histo = ttree.GetHistogram()
            histo.SetDirectory(0)

        # handle errors
        histo.Sumw2()

        if args.per_weight or args.save_weight:
            wcounter = tfile.Get('ntupler/weight_counter')
            if not wcounter:
                # try top level
                wcounter = tfile.Get('weight_counter')
            assert bool(wcounter)
            if not weight_counter:
                weight_counter = wcounter
                weight_counter.SetDirectory(0)
                weight_counter.SetName("sumup_weight_counter")
            else:
                weight_counter.Add(wcounter)

        if not out_histo:
            out_histo = histo.Clone()
            out_histo.SetName(histo_name)
            out_histo.SetDirectory(0)
            out_histo.Sumw2()
        else:
            out_histo.Add(histo)

    if args.scan:
        print filename
        ttree.Scan(args.draw_com, args.cond_com)

    # on closing the file all objects in it might be deleted
    # therefore move away the temp histo
    if temp_output_histo:
        temp_output_histo.SetDirectory(0)

    tfile.Close()

if args.per_weight:
    #weight_counter = tfile.Get('ntupler/weight_counter')
    out_histo.Scale(1./weight_counter.GetBinContent(2))

if args.try_xsec:
    matched_dtag = None
    for dtag in dtags:
        if dtag in args.output:
            matched_dtag = dtag
            break

    nickname, xsec = dtags.get(matched_dtag, 1.)
    out_histo.Scale(xsec)

if args.scale:
    out_histo.Scale(args.scale)

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

fout.cd()

# corresponding to old protocol everything is saved as
# channel/process/systematic/channel_process_systematic_distr
# -- the final part is the name of the histogram, which must be unique in root
#    the rest is for convenience
# in principle, root handles many same name histos in different directories, but can run into some weirdness on practice
# that's why this protocol is still kept

# if trying the xsec -- plug the nickname too and construct standard name if needed
if args.try_xsec and len(histo_path) > 1:
    histo_path[1] = nickname
    logging.debug("new path: " + '/'.join(histo_path))
    if args.std_histo_name:
        #
        distr = args.std_histo_name
        histo_name = "{path_part}_{distr}".format(path_part='_'.join(histo_path), distr=distr)
        logging.debug("new histo name: " + histo_name)
        out_histo.SetName(histo_name)

# check if this directory already exists in the file (a feature for future)
out_dir_name = ''.join(part + '/' for part in histo_path)
if out_dir_name and fout.Get(out_dir_name):
    logging.debug('found  ' + out_dir_name)
    out_dir = fout.Get(out_dir_name)
else:
    logging.debug('making ' + out_dir_name)
    # iteratively create each directory fout -> part0 -> part1 -> ...
    # somehow root did not work for creating them in 1 go
    out_dir = fout
    for directory in histo_path:
        nested_dir = out_dir.Get(directory) if out_dir.Get(directory) else out_dir.mkdir(directory)
        nested_dir.cd()
        out_dir = nested_dir

# save the histogram in this nested directory
#for histo in histos.values():
#    histo.SetDirectory(out_dir)
#    histo.Write()

out_histo.SetDirectory(out_dir)
out_histo.Write()

# save weight counters etc in the top directory in the file
if args.save_weight:
    fout.cd()
    weight_counter.SetName('weight_counter')
    weight_counter.Write()

fout.Write()
fout.Close()

