
'''
Attaching file fitDiagnosticsMuShapesPostFitDistr.root as _file0...
.l(TFile *) 0xc21770
root [1] .ls
TFile**		fitDiagnosticsMuShapesPostFitDistr.root	
 TFile*		fitDiagnosticsMuShapesPostFitDistr.root	
  KEY: TDirectoryFile	shapes_prefit;1	shapes_prefit
  KEY: RooArgSet	norm_prefit;1	Set of RooAbsArg objects
  KEY: TProcessID	ProcessID0;1	08e517dc-22df-11e8-9717-0864c10abeef
  KEY: RooFitResult	nuisances_prefit_res;1	
  KEY: RooArgSet	nuisances_prefit;1	Set of RooAbsArg objects
  KEY: RooFitResult	fit_s;1	
  KEY: TDirectoryFile	shapes_fit_s;1	shapes_fit_s
  KEY: RooArgSet	norm_fit_s;1	Set of RooAbsArg objects
  KEY: TTree	tree_fit_sb;1	tree_fit_sb
  KEY: TTree	tree_fit_b;1	tree_fit_b
  KEY: TTree	tree_prefit;1	tree_prefit

root [3] shapes_fit_s->Get("ctr_old_mu_sel_lj")->ls()
TDirectoryFile*		ctr_old_mu_sel_lj	ctr_old_mu_sel_lj
 KEY: TGraphAsymmErrors	data;1	dibosons in ctr_old_mu_sel_lj
 KEY: TH1F	dibosons;1	dibosons in ctr_old_mu_sel_lj
 KEY: TH1F	dy_other;1	dy_other in ctr_old_mu_sel_lj
 KEY: TH1F	dy_tautau;1	dy_tautau in ctr_old_mu_sel_lj
 KEY: TH1F	qcd;1	qcd in ctr_old_mu_sel_lj
 KEY: TH1F	s_top_lj;1	s_top_lj in ctr_old_mu_sel_lj
 KEY: TH1F	s_top_mutau;1	s_top_mutau in ctr_old_mu_sel_lj
 KEY: TH1F	s_top_other;1	s_top_other in ctr_old_mu_sel_lj
 KEY: TH1F	tt_lj;1	tt_lj in ctr_old_mu_sel_lj
 KEY: TH1F	tt_mutau;1	tt_mutau in ctr_old_mu_sel_lj
 KEY: TH1F	tt_other;1	tt_other in ctr_old_mu_sel_lj
 KEY: TH1F	tt_taultauh;1	tt_taultauh in ctr_old_mu_sel_lj
 KEY: TH1F	wjets;1	wjets in ctr_old_mu_sel_lj
 KEY: TH1F	total;1	Total signal+background in ctr_old_mu_sel_lj
 KEY: TH1F	total_signal;1	Total signal in ctr_old_mu_sel_lj
 KEY: TH1F	total_background;1	Total background in ctr_old_mu_sel_lj
 KEY: TH2F	total_covar;1	Covariance signal+background

'''

import argparse
import logging
from os.path import isfile
from sys import exit

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "stack the postfit distrs from higgsCombine",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument("data_file", help="data file name")
parser.add_argument("channel",   help="channel name")

parser.add_argument("-o", "--output-directory", type=str, default='./', help="optional output directory")
parser.add_argument("--ratio-range", type=float, default=0.5, help="range of ratio plot (1-range 1+range)")
parser.add_argument("-r", "--ratio", action='store_true', help="add ratio plot")
parser.add_argument("--y-max",       type=float,              help="set the maximum on Y axis")

parser.add_argument("--lumi",  type=float, default=35.8, help="lumi reported")
parser.add_argument("--title", type=str,   default='',   help="optional title")
parser.add_argument("--title-y", type=str,   default='Evts.',   help="optional title")
parser.add_argument("--title-x", type=str,   default='M_{T}',   help="optional title")

parser.add_argument("--mu", action='store_true', help="aim mu processes")

parser.add_argument("--prefit", action='store_true', help="plot prefit instead")

parser.add_argument("--no-data", action='store_true', help="don't plot data")

parser.add_argument("--coarse-bins", action='store_true', help="the 10 bins from 0 to 200 for Mt")



args = parser.parse_args()

assert isfile(args.data_file)

logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
gROOT.SetBatch()
gStyle.SetOptStat(0)

import plotting_root

gROOT.SetBatch()
gStyle.SetOptStat(0)

from plotting_root import nick_order, nick_colour, rgb

fdata = TFile(args.data_file)
chan_dir = fdata.Get(("shapes_prefit/" if args.prefit else "shapes_fit_s/") + args.channel)


mc_processes_mu = [
      "dibosons",
      "dy_other",
      "dy_tautau",
      "qcd",
      "s_top_lj",
      "s_top_mutau",
      "s_top_other",
      "tt_lj",
      "tt_mutau",
      "tt_other",
      "tt_taulj",
      "tt_taultauh",
      "wjets",
      "wjets_taul",
]

'''
 KEY: TH1F	qcd;1	qcd in ctr_old_mu_sel_lj
 KEY: TH1F	dibosons;1	dibosons in ctr_old_mu_sel_lj
 KEY: TH1F	dy_other;1	dy_other in ctr_old_mu_sel_lj
 KEY: TH1F	dy_tautau;1	dy_tautau in ctr_old_mu_sel_lj
 KEY: TH1F	s_top_lj;1	s_top_lj in ctr_old_mu_sel_lj
 KEY: TH1F	s_top_mutau;1	s_top_mutau in ctr_old_mu_sel_lj
 KEY: TH1F	s_top_other;1	s_top_other in ctr_old_mu_sel_lj
 KEY: TH1F	tt_lj;1	tt_lj in ctr_old_mu_sel_lj
 KEY: TH1F	tt_mutau;1	tt_mutau in ctr_old_mu_sel_lj
 KEY: TH1F	tt_other;1	tt_other in ctr_old_mu_sel_lj
 KEY: TH1F	tt_taulj;1	tt_taulj in ctr_old_mu_sel_lj
 KEY: TH1F	tt_taultauh;1	tt_taultauh in ctr_old_mu_sel_lj
 KEY: TH1F	wjets;1	wjets in ctr_old_mu_sel_lj
 KEY: TH1F	wjets_taul;1	wjets_taul in ctr_old_mu_sel_lj
'''

mc_processes_el = [
      "dibosons",
      "dy_other",
      "dy_tautau",
      "qcd",
      "s_top_lj",
      "s_top_eltau",
      "s_top_other",
      "tt_lj",
      "tt_eltau",
      "tt_other",
      "tt_taulj",
      "tt_taultauh",
      "wjets",
      "wjets_taul",
]

mc_processes = mc_processes_mu if args.mu else mc_processes_el


mc_processes.sort(key=lambda k: nick_order[k])

data_higComb   = chan_dir.Get("data")
mc_sum_higComb = chan_dir.Get("total")

## test tgraph/histo
#data_higComb.Divide(mc_sum_higComb)
# -- does not work simply, need something more

def th_postfit(name):
    if args.coarse_bins:
        return TH1D(name, '', 10, 0, 200)
    else:
        return TH1D(name, '', 20, 0, 250)

data   = th_postfit(data_higComb.GetName() + '_u')
mc_sum = th_postfit(mc_sum_higComb.GetName() + '_u')

# convert the higcomb TGraph data to histogram
for bini in range(data.GetSize()):
    tgraph_content = data_higComb.Eval(bini-0.5) # must fit the center of the bin in tgraph
    tgraph_error   = data_higComb.GetErrorY(bini) # just symmetrical error
    data.SetBinContent (bini, tgraph_content)
    data.SetBinError   (bini, tgraph_error)

#data.SetFillStyle(3004)
#data.SetFillColor(1)
data.SetMarkerStyle(1)
data.SetMarkerColor(0)

# set mc_sum histogram
for bini in range(mc_sum_higComb.GetSize()+1):
    mc_sum.SetBinContent(bini,  mc_sum_higComb.GetBinContent(bini))
    mc_sum.SetBinError  (bini,  mc_sum_higComb.GetBinError(bini))

mc_sum.SetFillStyle(3004)
mc_sum.SetFillColor(1)
mc_sum.SetMarkerStyle(1)
mc_sum.SetMarkerColor(0)


hs = THStack("mc_stack", "mc_stack")

mc_histos = []
for nick in mc_processes:
    col = nick_colour[nick]
    histo_higComb = chan_dir.Get(nick)

    histo = th_postfit(histo_higComb.GetName() +'_u')
    for bini in range(histo_higComb.GetSize()+1):
        histo.SetBinContent (bini, histo_higComb.GetBinContent(bini))
        histo.SetBinError   (bini, histo_higComb.GetBinError(bini))

    histo.SetFillColor( col );
    #histo.SetLineColor( col ); # it's needed for shapes
    histo.SetMarkerStyle(20);
    #histo.SetLineStyle(proc_ocurance);
    histo.SetMarkerColor(col);
    #used_histos.append(histo) # hopefully root wont screw this up
    hs.Add(histo, "HIST")
    mc_histos.append(histo)

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

# copy these before the drawing settings applied
data_rel = data.Clone()
mc_sum_rel = mc_sum.Clone()


hs.SetTitle(args.title)
mc_sum.SetTitle(args.title)

#data.GetYaxis().SetTitleFont(63)
#data.GetYaxis().SetTitleSize(20)
mc_sum.GetYaxis().SetTitleFont(63)
mc_sum.GetYaxis().SetTitleSize(20)
mc_sum.GetXaxis().SetTitleFont(63)
mc_sum.GetXaxis().SetTitleSize(20)

# axis labels
#data.GetYaxis().SetLabelSize(0.02)
#data.GetXaxis().SetLabelSize(0.02)
mc_sum.GetYaxis().SetLabelSize(0.02)
mc_sum.GetXaxis().SetLabelSize(0.02)

#hs             .GetYaxis().SetTitleOffset(1.4)
#data  .GetYaxis().SetTitleOffset(1.5) # place the title not overlapping with labels...
mc_sum.GetYaxis().SetTitleOffset(1.5)

#hs             .SetYTitle(title_y)
#data  .SetYTitle(args.title_y)
mc_sum.SetYTitle(args.title_y)
mc_sum.SetXTitle(args.title_x)

if args.y_max:
    print "SETING Y MAXIMUM"
    data.SetMaximum(args.y_max)
    data.GetYaxis().SetRangeUser(0, args.y_max)
    mc_sum.SetMaximum(args.y_max)
    #mc_sum.GetYaxis().SetRange(0, args.y_max)
    mc_sum.GetYaxis().SetRangeUser(0, args.y_max)
    mc_sum.SetAxisRange(0, args.y_max, "Y")

data.GetXaxis().SetLabelFont(63)
data.GetXaxis().SetLabelSize(14) # labels will be 14 pixels
mc_sum.GetXaxis().SetLabelFont(63)
mc_sum.GetXaxis().SetLabelSize(14) # labels will be 14 pixels

data.GetYaxis().SetLabelFont(63)
data.GetYaxis().SetLabelSize(14) # labels will be 14 pixels
mc_sum.GetYaxis().SetLabelFont(63)
mc_sum.GetYaxis().SetLabelSize(14) # labels will be 14 pixels

if args.ratio:
    mc_sum .GetXaxis().SetLabelOffset(999)
    mc_sum .GetXaxis().SetLabelSize(0)
    data   .GetXaxis().SetLabelOffset(999)
    data   .GetXaxis().SetLabelSize(0)

mc_sum.Draw("e2")
hs.Draw("same")
mc_sum.Draw("same e2")

if not args.no_data:
    data.SetLineWidth(3)
    data.Draw("same P")

left_title = TPaveText(0.1, 0.9, 0.4, 0.94, "brNDC")
left_title.AddText("CMS preliminary at 13 TeV")
left_title.SetTextFont(1)

right_title = TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC")
right_title.AddText("L = %s fb^{-1}" % args.lumi)
right_title.SetTextFont(132)
right_title.SetFillColor(0)

left_title .Draw("same")
right_title.Draw("same")

''' TGraph crap
Y_title = TPaveText(0.01, 0.88, 0.1, 0.93, "brNDC")
Y_title.AddText(args.title_y)
Y_title.SetTextFont(132)
Y_title.SetFillColor(0)

X_title = TPaveText(0.7, 0.01, 0.9, 0.06, "brNDC")
X_title.AddText(args.title_x)
X_title.SetTextFont(132)
X_title.SetFillColor(0)

Y_title.Draw("same")
X_title.Draw("same")
'''


'''
if args.ratio:
    data.GetXaxis().SetLabelOffset(999)
    data.GetXaxis().SetLabelSize(0)
    data.Divide(mc
'''

if args.ratio:
    pad2.cd()

    # converting TGraph to THist
    # TGraph also does not have convenient get X of a point
    # due to this crap called "root" I cannot have the usual ratio plot in reasonable amount of time
    # thus I'm skipping it for now and will do it later
    #for bini in range(data.GetN()):
    #    
    #    if qcd_hist.GetBinContent(bini) < 0:
    #        qcd_hist.SetBinContent(bini, 0)

    #data_rel = data.Clone()
    #mc_sum_rel = mc_sum.Clone()

    #data_rel.GetXaxis().SetLabelOffset(999)
    #data_rel.GetXaxis().SetLabelSize(0)
    data_rel.GetXaxis().SetLabelFont(63)
    data_rel.GetXaxis().SetLabelSize(14) # labels will be 14 pixels
    mc_sum_rel.GetXaxis().SetLabelFont(63)
    mc_sum_rel.GetXaxis().SetLabelSize(14) # labels will be 14 pixels

    data_rel.GetYaxis().SetLabelFont(63)
    data_rel.GetYaxis().SetLabelSize(14) # labels will be 14 pixels
    mc_sum_rel.GetYaxis().SetLabelFont(63)
    mc_sum_rel.GetYaxis().SetLabelSize(14) # labels will be 14 pixels

    data_rel   .GetXaxis().SetTitleOffset(4.) # place the title not overlapping with labels...
    mc_sum_rel .GetXaxis().SetTitleOffset(4.)

    mc_sum_rel .GetYaxis().SetTitleFont(63)
    mc_sum_rel .GetYaxis().SetTitleSize(14)
    data_rel   .GetYaxis().SetTitleFont(63)
    data_rel   .GetYaxis().SetTitleSize(14)

    mc_sum_rel .GetXaxis().SetTitleFont(63)
    mc_sum_rel .GetXaxis().SetTitleSize(14)
    data_rel   .GetXaxis().SetTitleFont(63)
    data_rel   .GetXaxis().SetTitleSize(14)

    mc_sum_rel .GetYaxis().SetTitleOffset(1.5)
    data_rel   .GetYaxis().SetTitleOffset(1.5)

    mc_sum_rel .GetXaxis().SetTitleOffset(4.0)
    data_rel   .GetXaxis().SetTitleOffset(4.0)

    #mc_sum_rel .GetXaxis().SetLabelOffset(4.)
    #mc_sum_rel .GetXaxis().SetLabelSize(14)
    #data_rel   .GetXaxis().SetLabelOffset(4.)
    #data_rel   .GetXaxis().SetLabelSize(14)

    data_rel   .SetYTitle("Data/MC")
    mc_sum_rel .SetYTitle("Data/MC")

    data_rel   .SetXTitle(args.title_x)
    mc_sum_rel .SetXTitle(args.title_x)

    ratio_max = 1. + args.ratio_range
    ratio_min = 1. - args.ratio_range
    data_rel.SetMaximum(ratio_max)
    data_rel.SetMinimum(ratio_min)
    mc_sum_rel.SetMaximum(ratio_max)
    mc_sum_rel.SetMinimum(ratio_min)

    data_rel.Divide(mc_sum)
    #     data_rel.Divide(mc_sum)
    #TypeError: void TGraphAsymmErrors::Divide(const TH1* pass, const TH1* total, const char* opt = "cp") =>
    # a really degenerate system..
    data_rel.Draw()

    mc_sum_rel.Divide(mc_sum)
    mc_sum_rel.Draw("e2")
    data_rel.Draw("same")


out_file = args.output_directory + '/postfit-distr_%s_%s' % (args.data_file.split('.')[0], args.channel)
if args.prefit:
    out_file += '_prefit'

out_file += '.png'

cst.SaveAs(out_file)



