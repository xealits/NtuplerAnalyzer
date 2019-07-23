import argparse, logging
from pdb import set_trace

import ctypes
from ctypes import c_double

# Allocate array of double*
n = 3
#x     = (ctypes.POINTER(ctypes.c_double) * n)()
#y     = (ctypes.POINTER(ctypes.c_double) * n)()
#ymin  = (ctypes.POINTER(ctypes.c_double) * n)()
#ymax  = (ctypes.POINTER(ctypes.c_double) * n)()
#ymin2 = (ctypes.POINTER(ctypes.c_double) * n)()
#ymax2 = (ctypes.POINTER(ctypes.c_double) * n)()

x     = (ctypes.c_double * n)()
y     = (ctypes.c_double * n)()
ymin  = (ctypes.c_double * n)()
ymax  = (ctypes.c_double * n)()
ymin2 = (ctypes.c_double * n)()
ymax2 = (ctypes.c_double * n)()


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "plot the NLL scans for the given fitting release",
    epilog = "Example:\n    $ python plot_uncertainty_scans.py v37_test13_bunchFULLFIT2_2"
    )

parser.add_argument("fit_release",    help="the name tag of the fit")
parser.add_argument("--nll-limit",  type=float, default=10., help="limit the maximum of NLL on the plot (deafult 10)")
parser.add_argument("--xsec-scale", type=float, default=143.23, help="the name tag of the fit (default 143.23)")
parser.add_argument("--font-size",   type=int, default=40,  help="main font size")
parser.add_argument("--canvas-size", type=int, default=800, help="the size of the canvas -- adjust to fit the relative font size")
parser.add_argument("--title-x",    type=str,   default='\\text{visible cross section [pb]}', help="the name tag of the fit")
parser.add_argument("--offset-x",   type=float,  default=0.8, help='')
parser.add_argument("--offset-y",   type=float,  default=0.8, help='')
parser.add_argument("--output-type",      type=str,   default='png', help="the output format (png by default)")
parser.add_argument("--horizontal-lines", type=str,   help="add dashed grey horizontal lines at the Y axis given as <y1>[,<y2>]...")
parser.add_argument("--debug", action="store_true", help="debug logging")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

logging.info("importing ROOT")
import ROOT
from ROOT import TCanvas, TGraph, TLine, gStyle, gROOT, gPad, TFile, TTree, TPaveText, TLegend, kBlack, TColor, kGray
#include "TGraph.h"
#include "TAxis.h"
#include "TH1.h"
logging.info("done")

gStyle.SetTitleFontSize(args.font_size)
gROOT.ForceStyle()

gStyle.SetTitleW(0.7) # //per cent of the pad width
gStyle.SetTitleH(0.5) # //per cent of the pad height
# none of these work


def rgb(r, g, b):
    '''rgb(r, g, b):

    from:
    http://webhome.phy.duke.edu/~dmb60/the-guide/
    TColor* color = gROOT->GetColor(TColor::GetColor(red,green,blue));//Use ints from 0 to 255 
    color->SetAlpha(0.5);//0 is fully transparent, 1 fully opaque
    hist->SetFillColor(color->GetNumber());
    '''
    return TColor.GetColor(r, g, b)

colors = {'el': 43, 'mu': 45, 'both': kBlack}
colors = {'el': rgb(254,204,92),  'mu': rgb(255,255,178), 'both': kBlack}
colors = {'el': rgb(240,59,32),   'mu': rgb(254,204,92),  'both': kBlack}
colors = {'el': rgb(153,153,255), 'mu': rgb(146,0,255),   'both': kBlack}

names = {'el': 'e#tau_{h}', 'mu': '#mu#tau_{h}', 'both': 'combined'}

# accept el, bu, both: 140.318  146.057  143.229
draw_command = "2*deltaNLL:r*%f" % args.xsec_scale # visible cross section in both channels # full space "2*deltaNLL:r*831.76" # "2*deltaNLL:r"

'''
higgsCombineMuFullUncertainty.MultiDimFit.mH120.root higgsCombineMuNoTau.MultiDimFit.mH120.root higgsCombineMuNoSysButTOPPT.MultiDimFit.mH120.root

TTree* ttree_full = (TTree*) _file0->Get("limit")
n = ttree_full->Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
TGraph *g2 = new TGraph(430, ttree_full->GetV2(), ttree_full->GetV1())
'''


#the_files = {
#'el': ["higgsCombineElFullUncertainty.MultiDimFit.mH120.root", "higgsCombineElNoTau.MultiDimFit.mH120.root", "higgsCombineElNoSys.MultiDimFit.mH120.root"],
#'mu': ["higgsCombineMuFullUncertainty.MultiDimFit.mH120.root", "higgsCombineMuNoTau.MultiDimFit.mH120.root", "higgsCombineMuNoSys.MultiDimFit.mH120.root"],
#'both': ["higgsCombineBothFullUncertainty.MultiDimFit.mH120.root", "higgsCombineBothNoTau.MultiDimFit.mH120.root", "higgsCombineBothNoSys.MultiDimFit.mH120.root"],
#}
#the_files_expected = {
#'el': ["higgsCombineExpectedElFullUncertainty.MultiDimFit.mH120.root", "higgsCombineExpectedElNoTau.MultiDimFit.mH120.root", "higgsCombineExpectedElNoSys.MultiDimFit.mH120.root"],
#'mu': ["higgsCombineExpectedMuFullUncertainty.MultiDimFit.mH120.root", "higgsCombineExpectedMuNoTau.MultiDimFit.mH120.root", "higgsCombineExpectedMuNoSys.MultiDimFit.mH120.root"],
#'both': ["higgsCombineExpectedBothFullUncertainty.MultiDimFit.mH120.root", "higgsCombineExpectedBothNoTau.MultiDimFit.mH120.root", "higgsCombineExpectedBothNoSys.MultiDimFit.mH120.root"],
#}

nn = args.fit_release

nn_el   = nn + '_el'
nn_mu   = nn + '_mu'
nn_both = nn + '_both'

the_files = {
'el':   [ "higgsCombine%sFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn_el,   "higgsCombine%sNoTau.MultiDimFit.mH120.root" % nn_el,   "higgsCombine%sNoSys.MultiDimFit.mH120.root" % nn_el],
'mu':   [ "higgsCombine%sFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn_mu,   "higgsCombine%sNoTau.MultiDimFit.mH120.root" % nn_mu,   "higgsCombine%sNoSys.MultiDimFit.mH120.root" % nn_mu],
'both': [ "higgsCombine%sFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn_both, "higgsCombine%sNoTau.MultiDimFit.mH120.root" % nn_both, "higgsCombine%sNoSys.MultiDimFit.mH120.root" % nn_both],
}

the_files_expected = {
'el':   [ "higgsCombine%s_ExpectedFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn_el,   "higgsCombine%s_ExpectedNoTau.MultiDimFit.mH120.root" % nn_el,   "higgsCombine%s_ExpectedNoSys.MultiDimFit.mH120.root" % nn_el],
'mu':   [ "higgsCombine%s_ExpectedFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn_mu,   "higgsCombine%s_ExpectedNoTau.MultiDimFit.mH120.root" % nn_mu,   "higgsCombine%s_ExpectedNoSys.MultiDimFit.mH120.root" % nn_mu],
'both': [ "higgsCombine%s_ExpectedFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn_both, "higgsCombine%s_ExpectedNoTau.MultiDimFit.mH120.root" % nn_both, "higgsCombine%s_ExpectedNoSys.MultiDimFit.mH120.root" % nn_both],
}



gROOT.SetBatch();
gStyle.SetOptStat(0);

def plot(chan, plot_expected, plot_data, report_lumi=True):
   #c1.SetGrid();

   #c1 = TCanvas("c1","Expected Bias Band Graph",200,10,800,800);
   c1 = TCanvas("c1","Expected Bias Band Graph",200,10, args.canvas_size, args.canvas_size);
   c1.DrawFrame(0.85,0.75,1.15,1.25);

   gStyle.SetTitleW(0.7) # //per cent of the pad width
   gStyle.SetTitleH(0.5) # //per cent of the pad height
   # none of these work

   # prepare data files
   if plot_data:
           fn_full, fn_notau, fn_stat = the_files[chan]

           file_full  = TFile(fn_full)
           file_notau = TFile(fn_notau)
           file_stat  = TFile(fn_stat)

           ttree_full = file_full.Get("limit")
           n = ttree_full.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           g_full = TGraph(n, ttree_full.GetV2(), ttree_full.GetV1())

           ttree_notau = file_notau.Get("limit")
           n = ttree_notau.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           g_notau = TGraph(n, ttree_notau.GetV2(), ttree_notau.GetV1())

           ttree_stat = file_stat.Get("limit")
           n = ttree_stat.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           g_stat = TGraph(n, ttree_stat.GetV2(), ttree_stat.GetV1())

           g_full.SetLineWidth(4)
           g_notau.SetLineWidth(2)
           g_stat.SetLineWidth(2)
           g_stat.SetLineStyle(7)

           # removing the title
           #g_full .SetTitle(";\\text{fitted } #hat{r};") # ROOT latex cannot put a hat on a letter
           #g_full .SetTitle(";\\text{fitted } r;-2\\Delta log L")
           #g_full .SetTitle(";\\text{fitted signal strength};-2\\Delta log L")

           #g_full .SetTitle(";%s;-2#Delta log (L)" % args.title_x)
           g_full .SetTitle(";;")
           g_notau.SetTitle(";;")
           g_stat .SetTitle(";;")

           #g_full .GetXaxis().SetRange(0.75, 1.35)
           #g_notau.GetXaxis().SetRange(0.75, 1.35)
           #g_stat .GetXaxis().SetRange(0.75, 1.35)
           g_full .GetXaxis().SetRangeUser(0.75, 1.35)
           g_notau.GetXaxis().SetRangeUser(0.75, 1.35)
           g_stat .GetXaxis().SetRangeUser(0.75, 1.35)

           # font
           #set_trace()
           g_full .GetXaxis().SetTitleFont(133)
           g_notau.GetXaxis().SetTitleFont(133)
           g_stat .GetXaxis().SetTitleFont(133)

           g_full .GetXaxis().SetTitleSize(args.font_size)
           g_notau.GetXaxis().SetTitleSize(args.font_size)
           g_stat .GetXaxis().SetTitleSize(args.font_size)

           g_full .GetYaxis().SetTitleFont(133)
           g_notau.GetYaxis().SetTitleFont(133)
           g_stat .GetYaxis().SetTitleFont(133)

           g_full .GetYaxis().SetTitleSize(args.font_size)
           g_notau.GetYaxis().SetTitleSize(args.font_size)
           g_stat .GetYaxis().SetTitleSize(args.font_size)

           print "set up data plots", g_full

   # prepare expected files
   if plot_expected:
           fn_full, fn_notau, fn_stat = the_files_expected[chan]

           exp_file_full  = TFile(fn_full)
           exp_file_notau = TFile(fn_notau)
           exp_file_stat  = TFile(fn_stat)

           print fn_full, fn_notau, fn_stat

           exp_ttree_full = exp_file_full.Get("limit")
           n = exp_ttree_full.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           exp_g_full = TGraph(n, exp_ttree_full.GetV2(), exp_ttree_full.GetV1())

           exp_ttree_notau = exp_file_notau.Get("limit")
           n = exp_ttree_notau.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           exp_g_notau = TGraph(n, exp_ttree_notau.GetV2(), exp_ttree_notau.GetV1())

           exp_ttree_stat = exp_file_stat.Get("limit")
           n = exp_ttree_stat.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           exp_g_stat = TGraph(n, exp_ttree_stat.GetV2(), exp_ttree_stat.GetV1())

           exp_g_full.SetLineWidth(2)
           exp_g_notau.SetLineWidth(2)
           exp_g_stat.SetLineWidth(2)
           exp_g_stat.SetLineStyle(7)

           # removing the title
           #g_full .SetTitle(";\\text{fitted } #hat{r};") # ROOT latex cannot put a hat on a letter
           #exp_g_full .SetTitle(";\\text{fitted } r;-2\\Delta log L")
           #exp_g_full .SetTitle(";\\text{fitted signal strength};-2\\Delta log L")
           exp_g_full .SetTitle(";%s;-2 #Delta log (L)" % args.title_x)
           exp_g_notau.SetTitle(";;")
           exp_g_stat .SetTitle(";;")

           #exp_g_full .GetXaxis().SetRange(0.75, 1.35)
           #exp_g_notau.GetXaxis().SetRange(0.75, 1.35)
           #exp_g_stat .GetXaxis().SetRange(0.75, 1.35)
           exp_g_full .GetXaxis().SetRangeUser(0.75, 1.35)
           exp_g_notau.GetXaxis().SetRangeUser(0.75, 1.35)
           exp_g_stat .GetXaxis().SetRangeUser(0.75, 1.35)

           # fonts
           exp_g_full .GetXaxis().SetTitleFont(133)
           exp_g_notau.GetXaxis().SetTitleFont(133)
           exp_g_stat .GetXaxis().SetTitleFont(133)

           exp_g_full .GetXaxis().SetTitleSize(args.font_size)
           exp_g_notau.GetXaxis().SetTitleSize(args.font_size)
           exp_g_stat .GetXaxis().SetTitleSize(args.font_size)

           exp_g_full .GetYaxis().SetTitleFont(133)
           exp_g_notau.GetYaxis().SetTitleFont(133)
           exp_g_stat .GetYaxis().SetTitleFont(133)

           exp_g_full .GetYaxis().SetTitleSize(args.font_size)
           exp_g_notau.GetYaxis().SetTitleSize(args.font_size)
           exp_g_stat .GetYaxis().SetTitleSize(args.font_size)

           exp_g_full .GetXaxis().SetTitleOffset(args.offset_x)
           exp_g_notau.GetXaxis().SetTitleOffset(args.offset_x)
           exp_g_stat .GetXaxis().SetTitleOffset(args.offset_x)

           exp_g_full .GetYaxis().SetTitleOffset(args.offset_y)
           exp_g_notau.GetYaxis().SetTitleOffset(args.offset_y)
           exp_g_stat .GetYaxis().SetTitleOffset(args.offset_y)

           print "set up expectation plots", exp_g_full


   ## add X axis title
   #g_full .SetXTitle("fitted \\mu")
   #g_notau.SetXTitle("fitted \\mu")
   #g_stat .SetXTitle("fitted \\mu")
   leg = TLegend(0.6, 0.7, 0.89, 0.89)
   leg.SetTextFont(133)
   leg.SetTextSize(args.font_size)

   # draw horizontal lines first
   drawn = False
   if args.horizontal_lines:
       for y in [float(s) for s in args.horizontal_lines.split(',')]:
           l = TLine(100, y, 1000, y)
           l.SetLineColor(kGray)
           if not drawn: # root is awesome, isn't it?
               drawn = True
               l.Draw()
           else:
               l.Draw('same')
           print "drew line", y

   if plot_expected and plot_data:
      print "plotting both"
      exp_g_full.SetLineStyle(7)

      if not drawn:
          exp_g_full .Draw("ap") # this cast makes the following work
          exp_g_full .Draw("L")
      else:
          exp_g_full .Draw("same ap")
          exp_g_full .Draw("same L")

      if args.horizontal_lines:
         for y in [float(s) for s in args.horizontal_lines.split(',')]:
             l = TLine(0, y, 1000, y)
             l.SetLineColor(kGray)
             l.Draw('same L')
             print "drew line", y

      exp_g_full .Draw("same ap")
      exp_g_full .Draw("same L")

      #exp_g_full .Draw("ap") # this cast makes the following work
      #exp_g_full .Draw("L")
      #exp_g_notau.Draw("L same")
      #exp_g_stat .Draw("L same")

      leg.AddEntry(exp_g_full, "expected", "L")
      #leg.AddEntry(exp_g_full, "exp. full unc.", "L")
      #leg.AddEntry(exp_g_stat, "exp. stat unc.", "L")

      g_full .Draw("L same")
      #g_notau.Draw("L same")
      #g_stat .Draw("L same")

      leg.AddEntry(g_full, "observed", "L")
      #leg.AddEntry(g_full, "fitted full unc.", "L")
      #leg.AddEntry(g_stat, "fitted stat unc.", "L")

   elif plot_data:
      print "plotting data"
      g_full .Draw("ap") # this cast makes the following work
      g_full .Draw("L")
      g_notau.Draw("L same")
      g_stat .Draw("L same")

      leg.AddEntry(g_full,  "fitted full unc.", "L")
      leg.AddEntry(g_notau, "#splitline{fitted unc.}{w.o. tau ID}", "L")
      leg.AddEntry(g_stat,  "fitted stat unc.", "L")

   elif plot_expected:
      print "plotting expected"

      exp_g_full .SetLineColor(43)
      exp_g_notau.SetLineColor(43)
      exp_g_stat .SetLineColor(43)

      exp_g_full .Draw("ap")
      exp_g_full .Draw("L")
      exp_g_notau.Draw("L same")
      exp_g_stat .Draw("L same")

      leg.AddEntry(exp_g_full,  "exp. full unc.", "L")
      leg.AddEntry(exp_g_notau, "#splitline{exp. unc.}{w.o. tau ID}", "L")
      leg.AddEntry(exp_g_stat,  "exp. stat unc.", "L")

   #left_title = TPaveText(0.1, 0.9, 0.4, 0.94, "brNDC")
   #left_title.AddText("CMS preliminary at 13 TeV")
   #left_title.SetTextFont(1)

   #left_title = TPaveText(0.15, 0.82, 0.3, 0.88, "brNDC")
   #left_title.AddText("CMS")
   #left_title.SetTextFont(1)
   #left_title.SetFillColor(0)

   left_title = TPaveText(0.1, 0.92, 0.5, 0.99, "brNDC")
   left_title.SetTextAlign(13)
   left_title.SetMargin(0)
   left_title.AddText("CMS")
   left_title.SetFillColor(0)

   left_title.Draw("same")

   #right_title = TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC")
   #right_title.AddText("L = %s fb^{-1}" % (31.3 if chan == 'el' else 35.8))
   #right_title.SetTextFont(132)
   #right_title.SetFillColor(0)

   #right_title = TPaveText(0.5, 0.9, 0.9, 0.95, "brNDC")
   right_title = TPaveText(0.5, 0.92, 0.9, 0.99, "brNDC")
   both = True
   if report_lumi:
       right_title.AddText("%s fb^{-1} (13 TeV)" % 35.8)
   elif both:
       right_title.AddText("%s fb^{-1} (13 TeV)" % '35.8')
   else:
       right_title.AddText("(13 TeV)")

   right_title.SetMargin(0)
   right_title.SetTextAlign(33)
   right_title.SetTextFont(132)
   right_title.SetFillColor(0)
   right_title.Draw("same")

   leg.SetBorderSize(0)
   leg.Draw("same")

   plotted = ''
   plotted += '_exp' if plot_expected else ''
   plotted += '_obs' if plot_data else ''
   c1.SaveAs("uncertainty_scans_%s_%s%s.%s" % (args.fit_release, chan, plotted, args.output_type))


def plot_all(report_lumi=True):
    # prepare all plots
    plots = {}
    for chan in ('el', 'mu', 'both'):
           fn_full, fn_notau, fn_stat = the_files[chan]

           file_full  = TFile(fn_full)
           file_notau = TFile(fn_notau)
           file_stat  = TFile(fn_stat)

           ttree_full = file_full.Get("limit")
           n = ttree_full.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           g_full = TGraph(n, ttree_full.GetV2(), ttree_full.GetV1())

           ttree_notau = file_notau.Get("limit")
           n = ttree_notau.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           g_notau = TGraph(n, ttree_notau.GetV2(), ttree_notau.GetV1())

           ttree_stat = file_stat.Get("limit")
           n = ttree_stat.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           g_stat = TGraph(n, ttree_stat.GetV2(), ttree_stat.GetV1())

           if chan == 'both':
               g_full.SetLineWidth(4)
           else:
               g_full.SetLineWidth(2)
           g_notau.SetLineWidth(2)
           g_stat.SetLineWidth(2)
           g_stat.SetLineStyle(7)

           # removing the title
           #g_full .SetTitle(";\\text{fitted } #hat{r};") # ROOT latex cannot put a hat on a letter
           #g_full .SetTitle(";\\text{fitted } r;-2\\Delta log L")
           #g_full .SetTitle(";\\text{fitted signal strength};-2\\Delta log L")
           g_full .SetTitle(";%s;-2#Delta log (L)" % args.title_x)
           g_notau.SetTitle(";;")
           g_stat .SetTitle(";;")

           #g_full .GetXaxis().SetRange(0.75, 1.35)
           #g_notau.GetXaxis().SetRange(0.75, 1.35)
           #g_stat .GetXaxis().SetRange(0.75, 1.35)
           g_full .GetXaxis().SetRangeUser(0.75, 1.35)
           g_notau.GetXaxis().SetRangeUser(0.75, 1.35)
           g_stat .GetXaxis().SetRangeUser(0.75, 1.35)

           plots[chan] = g_full, g_notau, g_stat
           print "set up observed plots"

    plots_expected = {}
    # prepare expected plots
    for chan in ('el', 'mu', 'both'):
           fn_full, fn_notau, fn_stat = the_files_expected[chan]

           exp_file_full  = TFile(fn_full)
           exp_file_notau = TFile(fn_notau)
           exp_file_stat  = TFile(fn_stat)

           print fn_full, fn_notau, fn_stat

           exp_ttree_full = exp_file_full.Get("limit")
           n = exp_ttree_full.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           exp_g_full = TGraph(n, exp_ttree_full.GetV2(), exp_ttree_full.GetV1())

           exp_ttree_notau = exp_file_notau.Get("limit")
           n = exp_ttree_notau.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           exp_g_notau = TGraph(n, exp_ttree_notau.GetV2(), exp_ttree_notau.GetV1())

           exp_ttree_stat = exp_file_stat.Get("limit")
           n = exp_ttree_stat.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %f" % args.nll_limit, "L")
           exp_g_stat = TGraph(n, exp_ttree_stat.GetV2(), exp_ttree_stat.GetV1())

           if chan == 'both':
               exp_g_full.SetLineWidth(4)
           else:
               exp_g_full.SetLineWidth(2)
           exp_g_notau.SetLineWidth(2)
           exp_g_stat.SetLineWidth(2)
           exp_g_stat.SetLineStyle(7)

           # removing the title
           #g_full .SetTitle(";\\text{fitted } #hat{r};") # ROOT latex cannot put a hat on a letter
           #exp_g_full .SetTitle(";\\text{fitted } r;-2\\Delta log L")
           #exp_g_full .SetTitle(";\\text{fitted signal strength};-2\\Delta log L")
           exp_g_full .SetTitle(";%s;-2#Delta log (L)" % args.title_x)
           exp_g_notau.SetTitle(";;")
           exp_g_stat .SetTitle(";;")

           #exp_g_full .GetXaxis().SetRange(0.75, 1.35)
           #exp_g_notau.GetXaxis().SetRange(0.75, 1.35)
           #exp_g_stat .GetXaxis().SetRange(0.75, 1.35)
           exp_g_full .GetXaxis().SetRangeUser(0.75, 1.35)
           exp_g_notau.GetXaxis().SetRangeUser(0.75, 1.35)
           exp_g_stat .GetXaxis().SetRangeUser(0.75, 1.35)

           plots_expected[chan] = exp_g_full, exp_g_notau, exp_g_stat
           print "set up expectation plots"

    leg = TLegend(0.6, 0.7, 0.89, 0.89)

    # plotting
    # first plots to set up the canvas

    # draw horizontal lines first
    drawn = False
    if args.horizontal_lines:
        for y in [float(s) for s in args.horizontal_lines.split(',')]:
            l = TLine(0.1, y, 0.9, y)
            l.SetLineColor(kGray)
            if not drawn: # root is awesome, isn't it?
                drawn = True
                l.Draw()
            else:
                l.Draw('same')

    chan = 'el'
    print 'plotting %s' % chan
    exp_g_full = plots_expected[chan][0]
    g_full     = plots[chan][0]

    exp_g_full.SetLineColor(colors[chan])
    g_full    .SetLineColor(colors[chan])

    if not drawn:
        g_full .Draw("ap") # this cast makes the following work
        g_full .Draw("L")
    else:
        g_full .Draw("same ap")
        g_full .Draw("same L")
    #g_notau.Draw("L same")
    #g_stat .Draw("L same")

    name = names[chan]
    leg.AddEntry(g_full, "%s" % name, "L")
    #leg.AddEntry(g_full, "%s observed" % name, "L")
    #leg.AddEntry(g_full, "fitted full unc.", "L")
    #leg.AddEntry(g_stat, "fitted stat unc.", "L")

    #exp_g_full.SetLineStyle(7)

    ##exp_g_full .Draw("ap")
    #exp_g_full .Draw("L same")
    ##exp_g_notau.Draw("L same")
    ##exp_g_stat .Draw("L same")

    #leg.AddEntry(exp_g_full, "%s expected" % chan, "L")
    ##leg.AddEntry(exp_g_full, "exp. full unc.", "L")
    ##leg.AddEntry(exp_g_stat, "exp. stat unc.", "L")

    # add el and mu
    for chan in ('mu', 'both'):
        exp_g_full = plots_expected[chan][0]
        g_full     = plots[chan][0]

        exp_g_full.SetLineColor(colors[chan])
        g_full    .SetLineColor(colors[chan])

        name = names[chan]
        g_full .Draw("L same")
        #leg.AddEntry(g_full, "%s observed" % name, "L")
        leg.AddEntry(g_full, "%s" % name, "L")

        #exp_g_full.SetLineStyle(7)
        ##exp_g_full .Draw("same ap") # this cast makes the following work
        #exp_g_full .Draw("same L")
        #leg.AddEntry(exp_g_full, "%s expected" % name, "L")


    #left_title = TPaveText(0.15, 0.82, 0.3, 0.88, "brNDC")
    #left_title.AddText("CMS")
    #left_title.SetTextFont(1)
    #left_title.SetFillColor(0)

    left_title = TPaveText(0.1, 0.92, 0.5, 0.99, "brNDC")
    left_title.SetTextAlign(13)
    left_title.SetMargin(0)
    left_title.AddText("CMS")
    left_title.SetFillColor(0)

    left_title.Draw("same")

    right_title = TPaveText(0.5, 0.9, 0.9, 0.99, "brNDC")
    both = True

    if report_lumi:
        right_title.AddText("%s fb^{-1} (13 TeV)" % (35.8 if chan == 'el' else 35.8))
    elif both:
        right_title.AddText("%s fb^{-1} (13 TeV)" % '35.8')
    else:
        right_title.AddText("(13 TeV)")

    right_title.SetTextFont(132)
    right_title.SetFillColor(0)
    right_title.Draw("same")

    leg.SetBorderSize(0)
    leg.Draw("same")

    plotted = ''
    plotted += '_exp'
    plotted += '_obs'
    c1.SaveAs("uncertainty_scans_%s_%s%s.%s" % (args.fit_release, 'all', plotted, args.output_type))

#plot('mu', True, True)
#c1.SaveAs("uncertainty_scans_%s_%s%s.png" % (args.fit_release, 'mu', '_exp_obs'))
#c1.Clear()
#plot('mu', True, False)
#c1.SaveAs("uncertainty_scans_%s_%s%s.png" % (args.fit_release, 'mu', '_exp'))
#c1.Clear()

#plot('el', True, True)
#plot('el', True, False)
#c1.Clear()
#
##plot('mu', False, True)
##plot('el', False, True)
#
## 
plot('el', True, True , True)
plot('el', True, False, True)

plot('mu', True, True , True)
plot('mu', True, False, True)

plot('both', True, True , True)
plot('both', True, False, True)
#c1.Clear()

#plot_all()

