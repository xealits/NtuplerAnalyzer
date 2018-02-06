import ROOT
from ROOT import TCanvas, TGraph, TLine, gStyle, gROOT, gPad, TPaveText, TPad
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


def grshade(channel_mu):
   gROOT.SetBatch();
   gStyle.SetOptStat(0);

   c1 = TCanvas("c1","Expected Bias Band Graph",10,10,800,800);

   #c1.SetGrid();
   c1.DrawFrame(0.85,0.75,1.15,1.25) # for some reason it draws the graph outside the canvas without this optin
   # I've never had to use this option when drawing histograms or stacks of histograms

   #pad1 = TPad("pad1","This is pad1", 0., 0.,  1., 1.)
   #pad1.Draw()
   #pad1.cd()


   '''
   const Int_t n = 3;
   Double_t x[n], y[n],
       ymin[n], ymax[n],
       ymin2[n], ymax2[n];
   Int_t i;
   '''

   '''
   n = 3
   x = [0, 0, 0]
   y = [0, 0, 0]
   ymin = [0, 0, 0]
   ymax = [0, 0, 0]
   ymin2 = [0, 0, 0]
   ymax2 = [0, 0, 0]
   '''

   '''
   for (i=0;i<n;i++) {
     x[i] = 0.1+i*0.1;
     ymax[i] = 10*sin(x[i]+0.2);
     ymin[i] = 8*sin(x[i]+0.1);
     y[i] = 9*sin(x[i]+0.15);
   }
   '''

   if channel_mu:
      x    [0] = c_double(0.9);
      y    [0] = c_double(0.8999); # median
      ymin [0] = c_double((0.8395 + 0.8405) * 0.5); # sigma dev
      ymax [0] = c_double(0.9645);
      ymin2[0] = c_double((0.7835 + 0.7845) * 0.5); # 2 sigma dev
      ymax2[0] = c_double((2*1.0345 + 4*1.0335) / (2+4));

      x    [1] = c_double(1.0);
      y    [1] = c_double(1.005);
      ymin [1] = c_double(0.934);
      ymax [1] = c_double(1.07);
      ymin2[1] = c_double(0.8735);
      ymax2[1] = c_double((1*1.1465 + 4*1.1455) / (1+4));

      x    [2] = c_double(1.1);
      y    [2] = c_double(1.0999);
      ymin [2] = c_double((1.0285 + 1.0295) * 0.5);
      ymax [2] = c_double(1.1764);
      ymin2[2] = c_double(0.9625);
      ymax2[2] = c_double((2*1.2585 + 4*1.2575) / (2+4));

   else:
      x    [0] = c_double(0.9);
      y    [0] = c_double((0.8995 + 0.9) * 0.5);
      ymin [0] = c_double((0.8375 + 0.8385) * 0.5);
      ymax [0] = c_double((0.9665 + 5*0.9675) / 6.);
      ymin2[0] = c_double(0.780);
      ymax2[0] = c_double(1.038);

      x    [1] = c_double(1.0);
      y    [1] = c_double(1.0);
      ymin [1] = c_double(0.9325);
      ymax [1] = c_double(1.0725);
      ymin2[1] = c_double(0.869);
      ymax2[1] = c_double(1.150);

      x    [2] = c_double(1.1);
      y    [2] = c_double((1.0995 + 1.1) * 0.5);
      ymin [2] = c_double((1.0265 + 2*1.0275) / 3.);
      ymax [2] = c_double(1.1786);
      ymin2[2] = c_double(0.96);
      ymax2[2] = c_double(1.263);


   grmin = TGraph(n,x,ymin);
   grmax = TGraph(n,x,ymax);
   gr    = TGraph(n,x,y);
   grshade  = TGraph(2*n);
   grshade2 = TGraph(2*n);

   for i in range(n): #(i=0;i<n;i++) {
      grshade.SetPoint(i,x[i],ymax[i]);
      grshade.SetPoint(n+i,x[n-i-1],ymin[n-i-1]);
      grshade2.SetPoint(i,x[i],ymax2[i]);
      grshade2.SetPoint(n+i,x[n-i-1],ymin2[n-i-1]);

   #gr.SetTitle("median;generated \\mu;fitted \\mu");
   gr      .SetTitle("median;generated;fitted");
   grmin   .SetTitle("median;generated;fitted");
   grmax   .SetTitle("median;generated;fitted");
   grshade .SetTitle("median;generated;fitted");
   grshade2.SetTitle("median;generated;fitted");

   gr      .GetYaxis().SetTitle("fitted \\mu");
   gr      .GetXaxis().SetTitle("generated \\mu");
   grmin   .GetYaxis().SetTitle("fitted \\mu");
   grmin   .GetXaxis().SetTitle("generated \\mu");
   grmax   .GetYaxis().SetTitle("fitted \\mu");
   grmax   .GetXaxis().SetTitle("generated \\mu");
   grshade .GetYaxis().SetTitle("fitted \\mu");
   grshade .GetXaxis().SetTitle("generated \\mu");
   grshade2.GetYaxis().SetTitle("fitted \\mu");
   grshade2.GetXaxis().SetTitle("generated \\mu");

   grshade2.GetHistogram().GetYaxis().SetTitle("fitted \\mu");
   grshade2.GetHistogram().GetXaxis().SetTitle("generated \\mu");

   ##gPad.Modified();

   #grshade.SetFillStyle(3013);
   #grshade2.SetFillColor(20);
   #grshade2.SetFillColor(18);
   grshade2.SetFillColor(5); # brazil plot is about 5 yellow 2sigma 3 green 1sigma
   grshade2.Draw("f");
   #grshade.SetFillColor(25);
   grshade.SetFillColor(3);
   grshade.Draw("same f");

   #grmin.Draw("l");
   #grmax.Draw("l");

   gr.SetLineWidth(0);
   #gr.SetLineStyle(10);
   gr.SetMarkerColor(4);
   gr.SetMarkerStyle(21);
   gr.Draw("CP");

   l = TLine(0.86, 0.86, 1.14, 1.14);
   l.SetLineWidth(1);
   #l.SetLineStyle(4);
   #l.SetLineColor(kRed);
   #l.GetXaxis.SetTitle("generated");
   #l.GetYaxis.SetTitle("fitted");
   l.Draw("same");

   left_title = TPaveText(0.1, 0.9, 0.4, 0.94, "brNDC")
   left_title.AddText("CMS preliminary at 13 TeV")
   left_title.SetTextFont(1)

   right_title = TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC")
   right_title.AddText("L = %s fb^{-1}" % (35.8 if channel_mu else 31.3))
   right_title.SetTextFont(132)
   right_title.SetFillColor(0)

   left_title .Draw("same")
   right_title.Draw("same")

   #gPad.Modified();

   # root has some absolutely different protocol
   # hence the axis titles are done manually
   title_x = TPaveText(0.7,0.12,0.9,.18, "brNDC")
   text_x = title_x.AddText("generated #mu")
   text_x.SetTextAlign(22)
   text_x.SetTextFont(43)
   text_x.SetTextSize(30)
   title_x.SetFillStyle(0)
   title_x.SetBorderSize(0)
   title_x.Draw('same')

   title_y = TPaveText(0.12,0.7,0.18,.9, "brNDC")
   text_y = title_y.AddText("fitted #mu")
   text_y.SetTextAlign(22)
   text_y.SetTextFont(43)
   text_y.SetTextSize(30)
   text_y.SetTextAngle(90.)
   title_y.SetFillStyle(0)
   title_y.SetBorderSize(0)
   title_y.Draw('same')

   c1.Update();
   gr.GetHistogram().GetXaxis().SetTitle("layer No");
   c1.Modified();

   if channel_mu:
      c1.SaveAs("expected-limits-band_mu.png");
   else:
      c1.SaveAs("expected-limits-band_el.png");


grshade(True)
grshade(False)

