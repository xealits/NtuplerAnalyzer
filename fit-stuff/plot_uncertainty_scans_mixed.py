import argparse, logging

import ctypes
from ctypes import c_double

import pdb

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
parser.add_argument("--title-x",    type=str,   default='\\text{visible cross section [pb]}', help="the name tag of the fit")
parser.add_argument("--output-type",      type=str,   default='png', help="the output format (png by default)")
parser.add_argument("--horizontal-lines", type=str,   help="add dashed grey horizontal lines at the Y axis given as <y1>[,<y2>]...")
parser.add_argument("--only-both", action="store_true", help="draw only the fit in both channels")
parser.add_argument("--inter", action="store_true", help="no batch mode")
parser.add_argument("--debug", action="store_true", help="debug logging")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

import ROOT
from ROOT import TCanvas, TGraph, TLine, gStyle, gROOT, gPad, TFile, TTree, TPaveText, TLegend, TLatex
from ROOT import TColor, kRed, kBlue, kBlack, kGreen

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
colors = {'el': kRed, 'mu': kBlue, 'both': kBlack}


#include "TGraph.h"
#include "TAxis.h"
#include "TH1.h"

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


'''
higgsCombineMuFullUncertainty.MultiDimFit.mH120.root higgsCombineMuNoTau.MultiDimFit.mH120.root higgsCombineMuNoSysButTOPPT.MultiDimFit.mH120.root

TTree* ttree_full = (TTree*) _file0->Get("limit")
n = ttree_full->Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
TGraph *g2 = new TGraph(430, ttree_full->GetV2(), ttree_full->GetV1())
'''

draw_command = "2*deltaNLL:r*%f" % args.xsec_scale # visible cross section in both channels # full space "2*deltaNLL:r*831.76" # draw_command

nn = args.fit_release

the_files = {
'el':   [  "higgsCombineel_%sFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn,   "higgsCombineel_%sNoTau.MultiDimFit.mH120.root" % nn,   "higgsCombineel_%sNoSys.MultiDimFit.mH120.root" % nn],
'mu':   [  "higgsCombinemu_%sFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn,   "higgsCombinemu_%sNoTau.MultiDimFit.mH120.root" % nn,   "higgsCombinemu_%sNoSys.MultiDimFit.mH120.root" % nn],
'both': ["higgsCombineboth_%sFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn, "higgsCombineboth_%sNoTau.MultiDimFit.mH120.root" % nn, "higgsCombineboth_%sNoSys.MultiDimFit.mH120.root" % nn],
}

the_files_expected = {
'el':   [  "higgsCombineel_%s_ExpectedFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn,   "higgsCombineel_%s_ExpectedNoTau.MultiDimFit.mH120.root" % nn,   "higgsCombineel_%s_ExpectedNoSys.MultiDimFit.mH120.root" % nn],
'mu':   [  "higgsCombinemu_%s_ExpectedFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn,   "higgsCombinemu_%s_ExpectedNoTau.MultiDimFit.mH120.root" % nn,   "higgsCombinemu_%s_ExpectedNoSys.MultiDimFit.mH120.root" % nn],
'both': ["higgsCombineboth_%s_ExpectedFullUncertaintyNoLumi.MultiDimFit.mH120.root" % nn, "higgsCombineboth_%s_ExpectedNoTau.MultiDimFit.mH120.root" % nn, "higgsCombineboth_%s_ExpectedNoSys.MultiDimFit.mH120.root" % nn],
}




gStyle.SetOptStat(0);
if not args.inter:
    gROOT.SetBatch()

c1 = TCanvas("c1","Expected Bias Band Graph",200,10,800,800);

#c1.SetGrid();
c1.DrawFrame(0.85,0.75,1.15,1.25);
c1.cd()

leg = TLegend(0.6, 0.7, 0.89, 0.89)

file_full_el_exp   = TFile("higgsCombineel_%s_ExpectedFullUncertainty.MultiDimFit.mH120.root" % nn)
file_full_mu_exp   = TFile("higgsCombinemu_%s_ExpectedFullUncertainty.MultiDimFit.mH120.root" % nn)
file_full_both_exp = TFile("higgsCombineboth_%s_ExpectedFullUncertainty.MultiDimFit.mH120.root" % nn)

file_full_el_obs   = TFile("higgsCombineel_%sFullUncertainty.MultiDimFit.mH120.root" % nn)
file_full_mu_obs   = TFile("higgsCombinemu_%sFullUncertainty.MultiDimFit.mH120.root" % nn)
file_full_both_obs = TFile("higgsCombineboth_%sFullUncertainty.MultiDimFit.mH120.root" % nn)


ttree_full_el_obs   = file_full_el_obs.Get("limit")
ttree_full_mu_obs   = file_full_mu_obs.Get("limit")
ttree_full_both_obs = file_full_both_obs.Get("limit")

ttree_full_el_exp   = file_full_el_exp.Get("limit")
ttree_full_mu_exp   = file_full_mu_exp.Get("limit")
ttree_full_both_exp = file_full_both_exp.Get("limit")


n             = ttree_full_el_obs.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %s" % args.nll_limit, "L")
g_full_el_obs = TGraph(n, ttree_full_el_obs.GetV2(), ttree_full_el_obs.GetV1())
n             = ttree_full_mu_obs.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %s" % args.nll_limit, "L")
g_full_mu_obs = TGraph(n, ttree_full_mu_obs.GetV2(), ttree_full_mu_obs.GetV1())
n             = ttree_full_both_obs.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %s" % args.nll_limit, "L")
g_full_both_obs = TGraph(n, ttree_full_both_obs.GetV2(), ttree_full_both_obs.GetV1())

n             = ttree_full_el_exp.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %s" % args.nll_limit, "L")
g_full_el_exp = TGraph(n, ttree_full_el_exp.GetV2(), ttree_full_el_exp.GetV1())
n             = ttree_full_mu_exp.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %s" % args.nll_limit, "L")
g_full_mu_exp = TGraph(n, ttree_full_mu_exp.GetV2(), ttree_full_mu_exp.GetV1())
n               = ttree_full_both_exp.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %s" % args.nll_limit, "L")
g_full_both_exp = TGraph(n, ttree_full_both_exp.GetV2(), ttree_full_both_exp.GetV1())

# the 1sigma and 2sigma graphs
sigma = 1.
n                      = ttree_full_both_exp.Draw(draw_command, "2*deltaNLL>0 && 2*deltaNLL< %s" % str(sigma), "L")
#pdb.set_trace()

# pyROOT does not implement Python's interface "len" for a TGraph -- it shows maximum entries for the extrapolation
# find the values
# to extend them for the 1-sgima area

sigma1_cros, sigma1_nll = [], []
try_next = True
i = 0
while try_next:
    cros, nll = ttree_full_both_exp.GetV2()[i], ttree_full_both_exp.GetV1()[i]
    i += 1
    if nll < sigma:
        sigma1_cros .append(cros)
        sigma1_nll  .append(nll)
    else:
        break

g_full_both_exp_1sigma = TGraph(n, ttree_full_both_exp.GetV2(), ttree_full_both_exp.GetV1())
# and add the points to close the graph
g_full_both_exp_1sigma.SetPoint(g_full_both_exp_1sigma.GetN(), sigma1_cros[-1], 0.)
g_full_both_exp_1sigma.SetPoint(g_full_both_exp_1sigma.GetN(), sigma1_cros[ 0], 0.)

g_full_both_exp_1sigma .SetFillColor(kGreen);
#g_full_both_exp_1sigma .SetFillColor(2);
g_full_both_exp_1sigma .SetFillStyle(3004);
g_full_both_exp_1sigma .SetLineColor(8);
#g_full_both_exp_1sigma .Draw("AF")

# https://github.com/cms-analysis/HiggsAnalysis-HiggsToTauTau/blob/90d9f39732a5b9f477558fcc39e43254656c6102/src/plottingTanb.cxx
# and kGreen for 1sigma 3004 not sure where it is used
#minus2sigma[i]->SetFillColor(TColor::GetColor(252,241,15));

# more examples:
## g_full_both_exp
#g_full_both_exp.SetFillColor(2);
#g_full_both_exp.SetFillStyle(3002);
#g_full_both_exp.SetLineColor(8);
#g_full_both_exp.SetLineWidth(-602);
##g_full_both_exp.Draw("AL");

#pdb.set_trace()

#g_full_both_exp_1sigma = TGraph(n, ttree_full_both_exp.GetV2(), ttree_full_both_exp.GetV1())
#g_full_both_exp_1sigma = ttree_full_both_exp.GetHistogram()

c1.Clear()

g_full_el_obs  .SetLineWidth(2)
g_full_el_exp  .SetLineWidth(2)
g_full_mu_obs  .SetLineWidth(2)
g_full_mu_exp  .SetLineWidth(2)

g_full_both_obs.SetLineWidth(4)
g_full_both_exp.SetLineWidth(4)

#g_full_el_obs  .SetLineColor(kRed)
#g_full_el_exp  .SetLineColor(kRed)
#g_full_mu_obs  .SetLineColor(kBlue)
#g_full_mu_exp  .SetLineColor(kBlue)
#g_full_both_obs  .SetLineColor(kBlack)
#g_full_both_exp  .SetLineColor(kBlack)

g_full_el_obs    .SetLineColor(colors['el'])
g_full_el_exp    .SetLineColor(colors['el'])
g_full_mu_obs    .SetLineColor(colors['mu'])
g_full_mu_exp    .SetLineColor(colors['mu'])
g_full_both_obs  .SetLineColor(colors['both'])
g_full_both_exp  .SetLineColor(colors['both'])

g_full_el_exp  .SetLineStyle(7)
g_full_mu_exp  .SetLineStyle(7)
g_full_both_exp.SetLineStyle(7)

#g_full_el_obs .SetTitle(";\\text{fitted signal strength};-2\\Delta ln L")
#g_full_el_exp .SetTitle(";\\text{fitted signal strength};-2\\Delta ln L")
g_full_el_obs .SetTitle(";%s;-2#Delta log (L);" % args.title_x)
g_full_el_exp .SetTitle(";%s;-2#Delta log (L);" % args.title_x)

g_full_mu_obs .SetTitle(";%s;-2#Delta log (L);" % args.title_x)
g_full_mu_exp .SetTitle(";%s;-2#Delta log (L);" % args.title_x)

g_full_both_obs .SetTitle(";%s;-2#Delta log (L);" % args.title_x)
g_full_both_exp .SetTitle(";%s;-2#Delta log (L);" % args.title_x)

g_full_both_exp .GetXaxis().SetRangeUser(600,1000) # works weirdly for tgraphs...

'''
g_full_el_obs.GetYaxis().SetTitle("-2\\Delta ln L")
g_full_el_exp.GetYaxis().SetTitle("-2\\Delta ln L")

g_full_mu_obs.GetYaxis().SetTitle("-2\\Delta ln L")
g_full_mu_exp.GetYaxis().SetTitle("-2\\Delta ln L")

g_full_both_obs.GetYaxis().SetTitle("-2\\Delta ln L")
g_full_both_exp.GetYaxis().SetTitle("-2\\Delta ln L")

g_full_both_obs.GetHistogram().GetYaxis().SetTitle("-2\\Delta ln L")
g_full_both_exp.GetHistogram().GetYaxis().SetTitle("-2\\Delta ln L")
'''

if args.only_both:
    #g_full_both_exp.Draw()
    ##g_full_both_exp.Draw('AL*')
    g_full_both_exp.Draw('ap')
    g_full_both_exp.Draw('L')
    #g_full_both_exp.Draw('same')
    g_full_both_obs.Draw('same')

    # 1sigma area
    #g_full_both_exp_1sigma.Draw("same")

    leg.AddEntry(g_full_both_obs, "observed", "L")
    leg.AddEntry(g_full_both_exp, "expected", "L")

    channels_drawn = ['both']

else:
    g_full_el_exp.Draw()
    g_full_el_obs.Draw('same')
    g_full_mu_exp.Draw('same')
    g_full_mu_obs.Draw('same')
    g_full_both_exp.Draw('same')
    g_full_both_obs.Draw('same')

    leg.AddEntry(g_full_el_obs, "observed, e#tau", "L")
    leg.AddEntry(g_full_el_exp, "expected, e#tau", "L")
    leg.AddEntry(g_full_mu_obs, "observed, #mu#tau", "L")
    leg.AddEntry(g_full_mu_exp, "expected, #mu#tau", "L")
    leg.AddEntry(g_full_both_obs, "observed, both", "L")
    leg.AddEntry(g_full_both_exp, "expected, both", "L")

    channels_drawn = ['el', 'mu', 'both']


## Y  axis title, which disappears from TGraph y axis
#y_title = TPaveText(0.0, 0.52, 0.2, 0.99, "brNDC")
#y_title.SetTextFont(132)
##y_title.SetTextAlign(13)
##y_title.SetTextAlign(22)
#y_title.SetTextAlign(12)
#y_title.SetMargin(0)
#y_title_text = y_title.AddText("-2 Delta ln L")
#y_title_text.SetTextAngle(90.)
#y_title_text.SetTextFont(132)
## -- the font size is wrong, it calculates the font size before rotation
#y_title.SetFillColor(0)
#y_title.SetTextFont(132)
#y_title.Draw("same")

#y_title = TLatex()
#y_title.SetTextSize(0.1)
#y_title.DrawLatex(0.07,0.5, "-2 Delta ln L")

# -- all this does not work,
#    to get the plot out of ROOT in pdf format
#    turn off the batch mode and draw interactively, adjusting the shape of the window if needed

#left_title = TPaveText(0.15, 0.82, 0.3, 0.88, "brNDC")
##left_title = TPaveText(0.12, 0.8, 0.2, 0.88, "brNDC")
##left_title.AddText("CMS preliminary at 13 TeV")
#left_title.AddText("CMS")
#left_title.SetTextFont(1)
#left_title.SetFillColor(0)
#left_title .Draw("same")

left_title = TPaveText(0.1, 0.92, 0.5, 0.99, "brNDC")
left_title.SetTextAlign(13)
left_title.SetMargin(0)
left_title.AddText("CMS")
left_title.SetFillColor(0)

left_title.Draw("same")

plot_lumi = False
#if plot_lumi:
#    #right_title = TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC")
#    right_title = TPaveText(0.65, 0.9, 0.9, 0.95, "brNDC")
#    right_title.AddText("%s fb^{-1} (13 TeV)" % (31.3 if chan == 'el' else 35.8))
#    right_title.SetTextFont(132)
#    right_title.SetFillColor(0)
#    right_title.Draw("same")

#right_title = TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC")
#right_title.AddText("L = %s fb^{-1}" % (31.3 if chan == 'el' else 35.8))
#right_title.SetTextFont(132)
#right_title.SetFillColor(0)

#right_title = TPaveText(0.5, 0.9, 0.9, 0.95, "brNDC")
right_title = TPaveText(0.5, 0.92, 0.9, 0.99, "brNDC")
if plot_lumi:
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

if args.inter:
    pdb.set_trace()

print "SAVING"
c1.SaveAs("uncertainty_scans_%s.%s" % ('_'.join(channels_drawn), args.output_type))
