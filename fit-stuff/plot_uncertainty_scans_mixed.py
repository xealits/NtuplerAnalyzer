import ROOT
from ROOT import TCanvas, TGraph, TLine, gStyle, gROOT, gPad, TFile, TTree, TPaveText, TLegend
from ROOT import kRed, kBlue, kBlack

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
n = ttree_full->Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
TGraph *g2 = new TGraph(430, ttree_full->GetV2(), ttree_full->GetV1())
'''


the_files = {
'el': ["higgsCombineElFullUncertainty.MultiDimFit.mH120.root", "higgsCombineElNoTau.MultiDimFit.mH120.root", "higgsCombineElNoSys.MultiDimFit.mH120.root"],
'mu': ["higgsCombineMuFullUncertainty.MultiDimFit.mH120.root", "higgsCombineMuNoTau.MultiDimFit.mH120.root", "higgsCombineMuNoSys.MultiDimFit.mH120.root"],
'both': ["higgsCombineBothFullUncertainty.MultiDimFit.mH120.root", "higgsCombineBothNoTau.MultiDimFit.mH120.root", "higgsCombineBothNoSys.MultiDimFit.mH120.root"],
}

the_files_expected = {
'el': ["higgsCombineExpectedElFullUncertainty.MultiDimFit.mH120.root", "higgsCombineExpectedElNoTau.MultiDimFit.mH120.root", "higgsCombineExpectedElNoSys.MultiDimFit.mH120.root"],
'mu': ["higgsCombineExpectedMuFullUncertainty.MultiDimFit.mH120.root", "higgsCombineExpectedMuNoTau.MultiDimFit.mH120.root", "higgsCombineExpectedMuNoSys.MultiDimFit.mH120.root"],
'both': ["higgsCombineExpectedBothFullUncertainty.MultiDimFit.mH120.root", "higgsCombineExpectedBothNoTau.MultiDimFit.mH120.root", "higgsCombineExpectedBothNoSys.MultiDimFit.mH120.root"],
}



gROOT.SetBatch();
gStyle.SetOptStat(0);

c1 = TCanvas("c1","Expected Bias Band Graph",200,10,800,800);

#c1.SetGrid();
c1.DrawFrame(0.85,0.75,1.15,1.25);
c1.cd()

leg = TLegend(0.6, 0.7, 0.89, 0.89)

file_full_el_exp   = TFile("higgsCombineExpectedElFullUncertainty.MultiDimFit.mH120.root")
file_full_mu_exp   = TFile("higgsCombineExpectedMuFullUncertainty.MultiDimFit.mH120.root")
file_full_both_exp = TFile("higgsCombineExpectedBothFullUncertainty.MultiDimFit.mH120.root")

file_full_el_obs   = TFile("higgsCombineElFullUncertainty.MultiDimFit.mH120.root")
file_full_mu_obs   = TFile("higgsCombineMuFullUncertainty.MultiDimFit.mH120.root")
file_full_both_obs = TFile("higgsCombineBothFullUncertainty.MultiDimFit.mH120.root")


ttree_full_el_obs   = file_full_el_obs.Get("limit")
ttree_full_mu_obs   = file_full_mu_obs.Get("limit")
ttree_full_both_obs = file_full_both_obs.Get("limit")

ttree_full_el_exp   = file_full_el_exp.Get("limit")
ttree_full_mu_exp   = file_full_mu_exp.Get("limit")
ttree_full_both_exp = file_full_both_exp.Get("limit")


n             = ttree_full_el_obs.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
g_full_el_obs = TGraph(n, ttree_full_el_obs.GetV2(), ttree_full_el_obs.GetV1())
n             = ttree_full_mu_obs.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
g_full_mu_obs = TGraph(n, ttree_full_mu_obs.GetV2(), ttree_full_mu_obs.GetV1())
n               = ttree_full_both_obs.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
g_full_both_obs = TGraph(n, ttree_full_both_obs.GetV2(), ttree_full_both_obs.GetV1())

n             = ttree_full_el_exp.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
g_full_el_exp = TGraph(n, ttree_full_el_exp.GetV2(), ttree_full_el_exp.GetV1())
n             = ttree_full_mu_exp.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
g_full_mu_exp = TGraph(n, ttree_full_mu_exp.GetV2(), ttree_full_mu_exp.GetV1())
n               = ttree_full_both_exp.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
g_full_both_exp = TGraph(n, ttree_full_both_exp.GetV2(), ttree_full_both_exp.GetV1())


c1.Clear()

g_full_el_obs  .SetLineWidth(2)
g_full_el_exp  .SetLineWidth(2)
g_full_mu_obs  .SetLineWidth(2)
g_full_mu_exp  .SetLineWidth(2)

g_full_both_obs.SetLineWidth(4)
g_full_both_exp.SetLineWidth(4)

g_full_el_obs  .SetLineColor(kRed)
g_full_el_exp  .SetLineColor(kRed)
g_full_mu_obs  .SetLineColor(kBlue)
g_full_mu_exp  .SetLineColor(kBlue)
g_full_both_obs  .SetLineColor(kBlack)
g_full_both_exp  .SetLineColor(kBlack)

g_full_el_exp  .SetLineStyle(7)
g_full_mu_exp  .SetLineStyle(7)
g_full_both_exp.SetLineStyle(7)

g_full_el_obs .SetTitle(";\\text{fitted signal strength};-2\\Delta ln L")
g_full_el_exp .SetTitle(";\\text{fitted signal strength};-2\\Delta ln L")

g_full_mu_obs .SetTitle(";;")
g_full_mu_exp .SetTitle(";;")

g_full_both_obs .SetTitle(";;")
g_full_both_exp .SetTitle(";;")



leg.AddEntry(g_full_el_obs, "observed, e#tau", "L")
leg.AddEntry(g_full_el_exp, "expected, e#tau", "L")

leg.AddEntry(g_full_mu_obs, "observed, #mu#tau", "L")
leg.AddEntry(g_full_mu_exp, "expected, #mu#tau", "L")

leg.AddEntry(g_full_both_obs, "observed, both", "L")
leg.AddEntry(g_full_both_exp, "expected, both", "L")

g_full_el_exp.Draw()
g_full_el_obs.Draw('same')
g_full_mu_exp.Draw('same')
g_full_mu_obs.Draw('same')
g_full_both_exp.Draw('same')
g_full_both_obs.Draw('same')


left_title = TPaveText(0.15, 0.82, 0.3, 0.88, "brNDC")
#left_title = TPaveText(0.12, 0.8, 0.2, 0.88, "brNDC")
#left_title.AddText("CMS preliminary at 13 TeV")
left_title.AddText("CMS")
left_title.SetTextFont(1)
left_title.SetFillColor(0)
left_title .Draw("same")

plot_lumi = False
if plot_lumi:
    #right_title = TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC")
    right_title = TPaveText(0.65, 0.9, 0.9, 0.95, "brNDC")
    right_title.AddText("%s fb^{-1} (13 TeV)" % (31.3 if chan == 'el' else 35.8))
    right_title.SetTextFont(132)
    right_title.SetFillColor(0)
    right_title.Draw("same")

leg.Draw("same")

print "SAVING"
c1.SaveAs("uncertainty_scans_%s.png" % '_'.join(['el', 'mu', 'both']))

