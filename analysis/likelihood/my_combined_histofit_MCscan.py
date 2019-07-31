from ROOT import *
from my_combined_histofit import *

# mc study (NLL scan plot finally? NO)
ws = RooWorkspace("ws") # now, this crap just depends on workspace to be set... holly cows on root, one thought finally we got rid of c++ -- nope

# POI is Parameter Of Interest
# it will be scaned
poi = RooArgSet(tt_factor)
poi.find("tt_factor").setRange(0.,2.)  #this is mostly for plotting

model = RooStats.ModelConfig()
model.SetWorkspace(ws)
model.SetPdf(models[0][1])
model.SetName("Model")
model.SetParametersOfInterest(poi)

bModel = model.Clone()
bModel.SetPdf(models[0][1])
bModel.SetName(model.GetName() + "_with_poi_1")
poi.find("tt_factor").setVal(1)
bModel.SetSnapshot(poi)

#First example is with a frequentist approach
fc = RooStats.FrequentistCalculator(data_rh[0], bModel, model)
fc.SetToys(2500,1500)

#Create hypotest inverter passing the desired calculator 
calc = RooStats.HypoTestInverter(fc)

#set confidence level (e.g. 95% upper limits)
calc.SetConfidenceLevel(0.95)

#use CLs
calc.UseCLs(1)

#reduce the noise
calc.SetVerbose(0)

#Configure ToyMC Samler
toymcs = calc.GetHypoTestCalculator().GetTestStatSampler()

#Use profile likelihood as test statistics 
profll = RooStats.ProfileLikelihoodTestStat(model.GetPdf())

#for CLs (bounded intervals) use one-sided profile likelihood
profll.SetOneSided(1)

#set the test statistic to use for toys
toymcs.SetTestStatistic(profll)

npoints = 8 #Number of points to scan
# min and max for the scan (better to choose smaller intervals)
poimin = poi.find("tt_factor").getMin()
poimax = poi.find("tt_factor").getMax()

print "Doing a fixed scan  in interval : ", poimin, " , ", poimax
calc.SetFixedScan(npoints,poimin,poimax);

result = calc.GetInterval() #This is a HypoTestInveter class object
upperLimit = result.UpperLimit()


#Now let's print the result of the two methods
#First the CLs
print "################"
print "The observed CLs upper limit is: ", upperLimit

#Compute expected limit
print "Expected upper limits, using the B (alternate) model : "
print " expected limit (median) ", result.GetExpectedUpperLimit(0)
print " expected limit (-1 sig) ", result.GetExpectedUpperLimit(-1)
print " expected limit (+1 sig) ", result.GetExpectedUpperLimit(1)
print "################"


#Plot now the result of the scan 

#First the CLs
freq_plot = RooStats.HypoTestInverterPlot("HTI_Result_Plot","Frequentist scan result for psi xsec",result)
#Then the Bayesian posterior
#bc_plot = bc.GetPosteriorPlot()

#Plot in a new canvas with style
dataCanvas = TCanvas("dataCanvas")
#dataCanvas.Divide(2,1)
dataCanvas.Divide(1, 1)
dataCanvas.SetLogy(0)
dataCanvas.cd(1)
freq_plot.Draw("2CL")
#dataCanvas.cd(2)
#bc_plot.Draw()
dataCanvas.SaveAs("my_combined_histofit_MCstudyCL.png")


