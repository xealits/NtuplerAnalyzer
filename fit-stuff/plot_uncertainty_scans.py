import ROOT
from ROOT import TCanvas, TGraph, TLine, gStyle, gROOT, gPad, TFile, TTree, TPaveText, TLegend
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




def plot(chan, plot_expected, plot_data, report_lumi=True):
   gROOT.SetBatch();
   gStyle.SetOptStat(0);

   c1 = TCanvas("c1","Expected Bias Band Graph",200,10,800,800);

   #c1.SetGrid();
   c1.DrawFrame(0.85,0.75,1.15,1.25);

   # prepare data files
   if plot_data:
           fn_full, fn_notau, fn_stat = the_files[chan]

           file_full  = TFile(fn_full)
           file_notau = TFile(fn_notau)
           file_stat  = TFile(fn_stat)

           ttree_full = file_full.Get("limit")
           n = ttree_full.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
           g_full = TGraph(n, ttree_full.GetV2(), ttree_full.GetV1())

           ttree_notau = file_notau.Get("limit")
           n = ttree_notau.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
           g_notau = TGraph(n, ttree_notau.GetV2(), ttree_notau.GetV1())

           ttree_stat = file_stat.Get("limit")
           n = ttree_stat.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
           g_stat = TGraph(n, ttree_stat.GetV2(), ttree_stat.GetV1())

           g_full.SetLineWidth(3)
           g_notau.SetLineWidth(2)
           g_stat.SetLineWidth(2)
           g_stat.SetLineStyle(7)

           # removing the title
           #g_full .SetTitle(";\\text{fitted } #hat{r};") # ROOT latex cannot put a hat on a letter
           #g_full .SetTitle(";\\text{fitted } r;-2\\Delta ln L")
           g_full .SetTitle(";\\text{fitted signal strength};-2\\Delta ln L")
           g_notau.SetTitle(";;")
           g_stat .SetTitle(";;")

           #g_full .GetXaxis().SetRange(0.75, 1.35)
           #g_notau.GetXaxis().SetRange(0.75, 1.35)
           #g_stat .GetXaxis().SetRange(0.75, 1.35)
           g_full .GetXaxis().SetRangeUser(0.75, 1.35)
           g_notau.GetXaxis().SetRangeUser(0.75, 1.35)
           g_stat .GetXaxis().SetRangeUser(0.75, 1.35)

           print "set up data plots", g_full

   # prepare expected files
   if plot_expected:
           fn_full, fn_notau, fn_stat = the_files_expected[chan]

           exp_file_full  = TFile(fn_full)
           exp_file_notau = TFile(fn_notau)
           exp_file_stat  = TFile(fn_stat)

           exp_ttree_full = exp_file_full.Get("limit")
           n = exp_ttree_full.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
           exp_g_full = TGraph(n, exp_ttree_full.GetV2(), exp_ttree_full.GetV1())

           exp_ttree_notau = exp_file_notau.Get("limit")
           n = exp_ttree_notau.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
           exp_g_notau = TGraph(n, exp_ttree_notau.GetV2(), exp_ttree_notau.GetV1())

           exp_ttree_stat = exp_file_stat.Get("limit")
           n = exp_ttree_stat.Draw("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL< 10", "L")
           exp_g_stat = TGraph(n, exp_ttree_stat.GetV2(), exp_ttree_stat.GetV1())

           exp_g_full.SetLineWidth(3)
           exp_g_full.SetLineStyle(7)
           exp_g_notau.SetLineWidth(2)
           exp_g_stat.SetLineWidth(2)
           exp_g_stat.SetLineStyle(7)

           #exp_g_full .SetLineColor(43)
           exp_g_notau.SetLineColor(43)
           exp_g_stat .SetLineColor(43)

           # removing the title
           #g_full .SetTitle(";\\text{fitted } #hat{r};") # ROOT latex cannot put a hat on a letter
           #exp_g_full .SetTitle(";\\text{fitted } r;-2\\Delta ln L")
           exp_g_full .SetTitle(";\\text{fitted signal strength};-2\\Delta ln L")
           exp_g_notau.SetTitle(";;")
           exp_g_stat .SetTitle(";;")

           #exp_g_full .GetXaxis().SetRange(0.75, 1.35)
           #exp_g_notau.GetXaxis().SetRange(0.75, 1.35)
           #exp_g_stat .GetXaxis().SetRange(0.75, 1.35)
           exp_g_full .GetXaxis().SetRangeUser(0.75, 1.35)
           exp_g_notau.GetXaxis().SetRangeUser(0.75, 1.35)
           exp_g_stat .GetXaxis().SetRangeUser(0.75, 1.35)

           print "set up expectation plots", exp_g_full


   ## add X axis title
   #g_full .SetXTitle("fitted \\mu")
   #g_notau.SetXTitle("fitted \\mu")
   #g_stat .SetXTitle("fitted \\mu")
   leg = TLegend(0.6, 0.7, 0.89, 0.89)

   if plot_expected and plot_data:
      print "plotting both"
      exp_g_full .Draw("ap") # this cast makes the following work
      exp_g_full .Draw("L")
      #exp_g_notau.Draw("L same")
      #exp_g_stat .Draw("L same")

      leg.AddEntry(exp_g_full, "exp. full unc.", "L")
      #leg.AddEntry(exp_g_stat, "exp. stat unc.", "L")

      g_full .Draw("L same")
      #g_notau.Draw("L same")
      #g_stat .Draw("L same")

      leg.AddEntry(g_full, "fitted full unc.", "L")
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

   left_title = TPaveText(0.15, 0.82, 0.3, 0.88, "brNDC")
   left_title.AddText("CMS")
   left_title.SetTextFont(1)
   left_title.SetFillColor(0)
   left_title.Draw("same")

   #right_title = TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC")
   #right_title.AddText("L = %s fb^{-1}" % (31.3 if chan == 'el' else 35.8))
   #right_title.SetTextFont(132)
   #right_title.SetFillColor(0)

   right_title = TPaveText(0.5, 0.9, 0.9, 0.95, "brNDC")
   if report_lumi:
       right_title.AddText("%s fb^{-1} (13 TeV)" % (31.3 if chan == 'el' else 35.8))
   else:
       right_title.AddText("(13 TeV)")
   right_title.SetTextFont(132)
   right_title.SetFillColor(0)
   right_title.Draw("same")

   leg.Draw("same")

   plotted = ''
   plotted += '_exp' if plot_expected else ''
   plotted += '_obs' if plot_data else ''
   c1.SaveAs("uncertainty_scans_%s%s.png" % (chan, plotted))

plot('mu', True, True)
plot('el', True, True)

plot('mu', True, False)
plot('el', True, False)

plot('mu', False, True)
plot('el', False, True)


# 
plot('both', True, True , False)
plot('both', True, False, False)

