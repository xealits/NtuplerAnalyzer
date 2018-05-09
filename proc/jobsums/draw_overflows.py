import logging
import argparse
import ctypes

logging.basicConfig(level=logging.INFO)


logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, gPad, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TLine, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
gROOT.SetBatch()
gStyle.SetOptStat(0)


'''
TH1F *DrawOverflow(TH1F* h)
{
//function to paint the histogram h with an extra bin for overflows
UInt_t nx = h->GetNbinsX()+1; Double_t *xbins= new Double_t[nx+1];
for (UInt_t i=0;i<nx;i++)
        xbins[i]=h->GetBinLowEdge(i+1);
xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
//book a temporary histogram having extra bins for overflows
TH1F *htmp = new TH1F(h->GetName(), h->GetTitle(), nx, xbins);
htmp->Sumw2();
//fill the new histogram including the overflows
for (UInt_t i=1; i<=nx; i++)
        {
        htmp->SetBinContent(htmp->FindBin(htmp->GetBinCenter(i)),h->GetBinContent(i)); htmp->SetBinError(htmp->FindBin(htmp->GetBinCenter(i)),h->GetBinError(i));
        }
htmp->SetBinContent(htmp->FindBin(h->GetBinLowEdge(1)-1), h->GetBinContent(0)); htmp->SetBinError(htmp->FindBin(h->GetBinLowEdge(1)-1), h->GetBinError(0));
// Restore the number of entries
htmp->SetEntries(h->GetEffectiveEntries());
return htmp;
}
'''

def DrawOverflow(h, overflow_bin_width = 10.):
    # function to paint the histogram h with an extra bin for overflows
    nx = h.GetNbinsX()+1
    #Double_t *xbins= new Double_t[nx+1];
    #for (UInt_t i=0;i<nx;i++)
    #        xbins[i]=h->GetBinLowEdge(i+1);
    #xbins[nx]=xbins[nx-1] + h->GetBinWidth(nx);
    edges = [h.GetBinLowEdge(i+1) for i in range(nx)]
    edges.append(edges[-1] + overflow_bin_width)
    xbins = (ctypes.c_double * len(edges))(* edges)
    # book a temporary histogram having extra bins for overflows
    htmp = TH1D(h.GetName(), h.GetTitle(), nx, xbins)
    htmp.Sumw2()
    # fill the new histogram including the overflows
    for i in range(nx+1):
        htmp.SetBinContent(htmp.FindBin(htmp.GetBinCenter(i)), h.GetBinContent(i))
        htmp.SetBinError  (htmp.FindBin(htmp.GetBinCenter(i)), h.GetBinError(i))
    # underflow ??
    #htmp->SetBinContent(htmp->FindBin(h->GetBinLowEdge(1)-1), h->GetBinContent(0));
    #htmp->SetBinError  (htmp->FindBin(h->GetBinLowEdge(1)-1), h->GetBinError(0));
    # Restore the number of entries
    htmp.SetEntries(h.GetEffectiveEntries())
    return htmp


if __name__ == '__main__':
    infile = TFile("merge-sets/v21/p4_ljbins/MC2016_Summer16_TTJets_powheg.root")
    histo = infile.Get("ctr_old_mu_presel/tt_lj/NOMINAL/ctr_old_mu_presel_tt_lj_NOMINAL_dijet_trijet_mass")
    cst = TCanvas("cst","stacked hists",10,10,700,700)
    h_out = DrawOverflow(histo)
    h_out.Draw()
    cst.SaveAs("overflows.png")

