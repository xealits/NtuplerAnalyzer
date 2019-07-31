from ROOT import *
from my_combined_histofit import *

# https://root.cern.ch/roofit-20-minutes

#ifitResult = simModel.fitTo(data_hist_combined, RooFit.ExternalConstraints(contraints_set))
#nll = models[0][1].createNLL(data_rh[0], RooFit.NumCPU(8), RooFit.ExternalConstraints(contraints_set))
nll = simModel.createNLL(data_hist_combined, RooFit.NumCPU(8), RooFit.ExternalConstraints(contraints_set))
#nll = simModel.createNLL(data_hist_combined, RooFit.NumCPU(8))
pll = nll.createProfile(RooArgSet(tt_factor))

frame = tt_factor.frame()
frame.SetTitle(categories_name)
nll.plotOn(frame)
pll.plotOn(frame, RooFit.LineColor(kRed))

frame.GetYaxis().SetRange(0,5)
frame.GetYaxis().SetRangeUser(0,5)

#Draw the results
c1 = TCanvas()
c1.Divide(1,1)

c1.cd(1)
frame.Draw()

c1.SaveAs("my_combined_histo_NLL.png")


