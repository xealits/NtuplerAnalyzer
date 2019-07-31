from ROOT import *
from my_combined_histofit import *

# https://root.cern.ch/root/html/tutorials/roofit/rf605_profilell.C.html
#
# // C o n s t r u c t   p l a i n   l i k e l i h o o d
#  // ---------------------------------------------------
#
#  // Construct unbinned likelihood
#  RooAbsReal* nll = model.createNLL(*data,NumCPU(2)) ;
#
#  // Minimize likelihood w.r.t all parameters before making plots
#  RooMinuit(*nll).migrad() ;
#
#  // Plot likelihood scan frac 
#  RooPlot* frame1 = frac.frame(Bins(10),Range(0.01,0.95),Title("LL and profileLL in frac")) ;
#  nll->plotOn(frame1,ShiftToZero()) ;
#
#  // Plot likelihood scan in sigma_g2
#  RooPlot* frame2 = sigma_g2.frame(Bins(10),Range(3.3,5.0),Title("LL and profileLL in sigma_g2")) ;
#  nll->plotOn(frame2,ShiftToZero())

nll = models[0][1].createNLL(data_rh[0], RooFit.NumCPU(2))

#RooMinuit(nll).migrad() # maybe it's minimized already?

frame1 = tt_factor.frame(RooFit.Bins(10), RooFit.Range(0.01, 1.95), RooFit.Title("LL"))
nll.plotOn(frame1, RooFit.ShiftToZero())

#// The profile likelihood estimator on nll for frac will minimize nll w.r.t
#  // all floating parameters except frac for each evaluation
#
#  RooAbsReal* pll_frac = nll->createProfile(frac) ;
#
#  // Plot the profile likelihood in frac
#  pll_frac->plotOn(frame1,LineColor(kRed)) ;
#
#  // Adjust frame maximum for visual clarity
#  frame1->SetMinimum(0) ;
#  frame1->SetMaximum(3) 

pll_tt = nll.createProfile(RooArgSet(tt_factor))

#// Plot the profile likelihood in frac
pll_tt.plotOn(frame1, RooFit.LineColor(kRed))

#// Adjust frame maximum for visual clarity
frame1.SetMinimum(0) ;
frame1.SetMaximum(3) 



dataCanvas = TCanvas("dataCanvas")
dataCanvas.Divide(1, 1)
dataCanvas.SetLogy(0)
dataCanvas.cd(1)
frame1.Draw("LL")

dataCanvas.SaveAs("my_combined_histofit_MCstudyLL.png")



