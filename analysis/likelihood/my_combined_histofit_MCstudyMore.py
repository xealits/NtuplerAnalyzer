from ROOT import *
from my_combined_histofit import *

# mc study (NO NLL scan plot)
ws = RooWorkspace("ws") # now, this crap just depends on workspace to be set... holly cows on root, one thought finally we got rid of c++ -- nope

# POI is Parameter Of Interest
# it will be scaned
poi = RooArgSet(tt_factor)

#Set the RooModelConfig and let it know what the content of the workspace is about

model = RooStats.ModelConfig()
model.SetWorkspace(ws)
model.SetPdf(models[0][1])
model.SetParametersOfInterest(poi)

#Here we explicitly set the value of the parameters for the null hypothesis
#We want no signal contribution, so cross_psi = 0

#cross_psi = ws.var("cross_psi")
nullParams = poi.snapshot()
#nullParams.setRealValue(tt_factor, 0.)
nullParams.setRealValue("tt_factor", 1.)

#Build the profile likelihood calculator
plc = RooStats.ProfileLikelihoodCalculator()

plc.SetData(data_rh[0])
plc.SetModel(model)
plc.SetParameters(poi)
plc.SetNullParameters(nullParams)

#We get a HypoTestResult out of the calculator, and we can query it.
htr = plc.GetHypoTest()

print "-------------------------------------------------"
print "The p-value for the null is ", htr.NullPValue()
print "Corresponding to a signifcance of ", htr.Significance()
print "-------------------------------------------------"

#PyROOT sometimes fails cleaning memory, this helps
del plc

'''
#Set confidence level
confidenceLevel = 0.68

plc.SetModel(model)
plc.SetParameters(poi)
plc.SetConfidenceLevel(confidenceLevel)

#Get the interval
pl_Interval = plc.GetInterval()

#Let's make a plot
dataCanvas = TCanvas("dataCanvas")
#dataCanvas.Divide(2,1)
#dataCanvas.Divide(1,1)

dataCanvas.cd(1)
plot_Interval = RooStats.LikelihoodIntervalPlot(pl_Interval)
plot_Interval.SetTitle("Profile Likelihood Ratio")
plot_Interval.SetMaximum(3.)
plot_Interval.Draw()

dataCanvas.SaveAs("my_combined_histofit_MCstudy_NLLscan.png")


'''

