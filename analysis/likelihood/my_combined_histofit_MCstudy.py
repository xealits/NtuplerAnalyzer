from ROOT import *
from my_combined_histofit import *
from my_combined_histofit import distr

print distr

# mc study (NO NLL scan plot!)

#Get the observable and PDF out of the Workspace
#mass = ws.var("mass")     # distr
#totPDF = ws.pdf("totPDF") # models [0][1]
#mc_study = RooMCStudy(models[0][1], RooArgSet(distr), RooFit.Extended(1), RooFit.FitOptions(RooFit.Save(1)) )

# distr is out observable parameter
# sample is Index of event categories
#mc_study = RooMCStudy(simModel, RooArgSet(distr, sample), RooFit.Extended(1), RooFit.FitOptions(RooFit.Save(1)), RooFit.ExternalConstraints(contraints_set) )
mc_study = RooMCStudy(simModel, RooArgSet(distr, sample), RooFit.Extended(1), RooFit.FitOptions(RooFit.Save(1)))
mc_study.generateAndFit(1)

frame_cross_par = mc_study.plotParam(tt_factor, RooFit.Bins(40))
frame_cross_err = mc_study.plotError(tt_factor, RooFit.Bins(40), RooFit.FrameRange(0.,3.) )
frame_cross_pul = mc_study.plotPull (tt_factor, RooFit.Bins(40), RooFit.FitGauss(1) )
#Also, let's see the distribution of the NLL for all the fits
frame_nll = mc_study.plotNLL(RooFit.Bins(40))
#Now plot the whole thing
gStyle.SetOptStat(0)

mcstudy_Canvas = TCanvas("mcstudy_Canvas")
mcstudy_Canvas.Divide(2,2)

mcstudy_Canvas.cd(1)
frame_cross_par.Draw()

mcstudy_Canvas.cd(2)
frame_cross_err.Draw()

mcstudy_Canvas.cd(3)
frame_cross_pul.Draw()

mcstudy_Canvas.cd(4)
frame_nll.Draw()

mcstudy_Canvas.SaveAs("my_combined_histofit_MCStudy.png")



