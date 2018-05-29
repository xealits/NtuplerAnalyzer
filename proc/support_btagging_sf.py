from os import environ
import logging
import ROOT
from ROOT import gSystem, TFile,  TH1D



'''
b-tagging SF

itak

standalone even if you get to it and  run .L
colapses with
error: constructor for 'BTagCalibration' must explicitly
      initialize the member 'data_'

-- same error as when the wrapper was made...

it should work like:
import ROOT

ROOT.gROOT.ProcessLine('.L BTagCalibrationStandalone.cc+') 

calib = ROOT.BTagCalibration("csv", "CSV.csv")
reader = ROOT.BTagCalibrationReader(0, "central")  # 0 is for loose op
reader.load(calib, 0, "comb")  # 0 is for b flavour, "comb" is the measurement type

# in your event loop
reader.eval(0, 1.2, 30.)  # jet flavor, eta, pt

-- maybe loading ttbar lib the calibrator class should be there?

yup:
>>> import ROOT
>>> ROOT.gROOT.Reset()
>>> ROOT.gROOT.ProcessLine(".L pu_weight.C+") # this is also needed for stuf to run
0L
>>> ROOT.gSystem.Load("libUserCodettbar-leptons-80X.so")
0
>>> 
>>> ROOT.BTagCalibration
<class 'ROOT.BTagCalibration'>

-- test tomorrow
'''

# SUPPOSE THE ttbar library is loaded already
#ROOT.gSystem.Load("libUserCodettbar-leptons-80X.so")


logging.info("loading b-tagging SF stuff")

#ROOT.gROOT.ProcessLine(".L pu_weight.C+") # this is also needed for stuf to run
# not sure I use it right now
#ROOT.gROOT.ProcessLine("set_bSF_calibrators()")


#calib = ROOT.BTagCalibration("csv", "CSV.csv")
#reader = ROOT.BTagCalibrationReader(0, "central")  # 0 is for loose op
#reader.load(calib, 0, "comb")  # 0 is for b flavour, "comb" is the measurement type

logging.info("loading b-tagging SF callibration")
#bCalib_filename = "${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/data/weights/CSVv2_Moriond17_B_H.csv"
bCalib_filename = environ['CMSSW_BASE'] + '/src/UserCode/ttbar-leptons-80X/data/weights/CSVv2_Moriond17_B_H.csv'
#gSystem.ExpandPathName(bCalib_filename)

logging.info("btag SFs from " + bCalib_filename)
btagCalib = ROOT.BTagCalibration("CSVv2", bCalib_filename)
#BTagCalibration* btagCalib; // ("CSVv2", string(std::getenv("CMSSW_BASE"))+"/src/UserCode/ttbar-leptons-80X/data/weights/CSVv2_Moriond17_B_H.csv");

logging.info("loading Reader")
btagCal_sys = ROOT.vector('string')()
btagCal_sys.push_back("up")
btagCal_sys.push_back("down")
btagCal = ROOT.BTagCalibrationReader(ROOT.BTagEntry.OP_MEDIUM,  # operating point
# BTagCalibrationReader btagCal(BTagEntry::OP_TIGHT,  // operating point
		 "central",            # central sys type
		 btagCal_sys)      # other sys types
# the type is:
# const std::vector<std::string> & otherSysTypes={}


logging.info("...flavours")
btagCal.load(btagCalib,        # calibration instance
    ROOT.BTagEntry.FLAV_B,      # btag flavour
#            "comb");          # they say "comb" is better precision, but mujets are independent from ttbar dilepton channels
    "mujets")              #
btagCal.load(btagCalib,        # calibration instance
    ROOT.BTagEntry.FLAV_C,      # btag flavour
    "mujets")              # measurement type
btagCal.load(btagCalib,        # calibration instance
    ROOT.BTagEntry.FLAV_UDSG,   # btag flavour
    "incl")                # measurement type

logging.info("loaded b-tagging callibration")

logging.info("loading b-tagging efficiencies")
#bTaggingEfficiencies_filename = std::string(std::getenv("CMSSW_BASE")) + "/src/UserCode/ttbar-leptons-80X/jobing/configs/b-tagging-efficiencies.root";
bTaggingEfficiencies_filename = '${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/analysis/b-tagging/v9.38-for-b-effs/beff_histos.root'
gSystem.ExpandPathName(bTaggingEfficiencies_filename)
bTaggingEfficiencies_file = TFile(bTaggingEfficiencies_filename)

bEff_histo_b = bTaggingEfficiencies_file.Get('MC2016_Summer16_TTJets_powheg_btag_b_hadronFlavour_candidates_tagged')
bEff_histo_c = bTaggingEfficiencies_file.Get('MC2016_Summer16_TTJets_powheg_btag_c_hadronFlavour_candidates_tagged')
bEff_histo_udsg = bTaggingEfficiencies_file.Get('MC2016_Summer16_TTJets_powheg_btag_udsg_hadronFlavour_candidates_tagged')

logging.info("loaded b-tagging efficiencies")

h_control_btag_eff_b    = TH1D("control_btag_eff_b",    "", 150, 0, 2)
h_control_btag_eff_c    = TH1D("control_btag_eff_c",    "", 150, 0, 2)
h_control_btag_eff_udsg = TH1D("control_btag_eff_udsg", "", 150, 0, 2)

h_control_btag_weight_b    = TH1D("control_btag_weight_b",    "", 150, 0, 2)
h_control_btag_weight_c    = TH1D("control_btag_weight_c",    "", 150, 0, 2)
h_control_btag_weight_udsg = TH1D("control_btag_weight_udsg", "", 150, 0, 2)

h_control_btag_weight_notag_b    = TH1D("control_btag_weight_notag_b",    "", 150, 0, 2)
h_control_btag_weight_notag_c    = TH1D("control_btag_weight_notag_c",    "", 150, 0, 2)
h_control_btag_weight_notag_udsg = TH1D("control_btag_weight_notag_udsg", "", 150, 0, 2)


#def calc_btag_sf_weight(hasCSVtag: bool, flavId: int, pt: float, eta: float) -> float:
def calc_btag_sf_weight(hasCSVtag, flavId, pt, eta, sys="central"):
    # int flavId=jet.partonFlavour();
    #int flavId=jet.hadronFlavour();
    # also: patJet->genParton().pdgId()
    # fill_btag_eff(string("mc_all_b_tagging_candidate_jets_pt_eta"), jet.pt(), eta, weight);

    sf, eff = 1.0, 1.0
    # If the jet is tagged -- weight *= SF of the jet
    # if not weight *= (1 - eff*SF)/(1 - eff)
    # 

    if abs(flavId) == 5:
        # get SF for the jet
        sf = btagCal.eval_auto_bounds(sys, ROOT.BTagEntry.FLAV_B, eta, pt, 0.)
        # get eff for the jet
        #eff = bTagging_b_jet_efficiency(pt, eta)
        eff = bEff_histo_b.GetBinContent(bEff_histo_b.FindBin(pt, eta))
        h_control_btag_eff_b.Fill(eff)
    elif abs(flavId) == 4:
        sf = btagCal.eval_auto_bounds(sys, ROOT.BTagEntry.FLAV_C, eta, pt, 0.)
        #eff = bTagging_c_jet_efficiency(pt, eta)
        eff = bEff_histo_c.GetBinContent(bEff_histo_c.FindBin(pt, eta))
        h_control_btag_eff_c.Fill(eff)
    else:
        sf = btagCal.eval_auto_bounds(sys, ROOT.BTagEntry.FLAV_UDSG, eta, pt, 0.)
        #eff = bTagging_udsg_jet_efficiency(pt, eta)
        eff = bEff_histo_udsg.GetBinContent(bEff_histo_udsg.FindBin(pt, eta))
        h_control_btag_eff_udsg.Fill(eff)

    jet_weight_factor = 1.
    if hasCSVtag: # a tagged jet
        jet_weight_factor = sf
        if abs(flavId) == 5:
            h_control_btag_weight_b.Fill(jet_weight_factor)
        elif abs(flavId) == 4:
            h_control_btag_weight_c.Fill(jet_weight_factor)
        else:
            h_control_btag_weight_udsg.Fill(jet_weight_factor)
    else: # not tagged
        # truncate efficiency to 0 and 0.99
        #eff = 0. if eff < 0. else (0.999 if eff > 0.999 else eff)
        if 1 - eff > 0:
            jet_weight_factor = (1 - sf*eff) / (1 - eff)
        else:
            jet_weight_factor = 1
        if abs(flavId) == 5:
            h_control_btag_weight_notag_b.Fill(jet_weight_factor)
        elif abs(flavId) == 4:
            h_control_btag_weight_notag_c.Fill(jet_weight_factor)
        else:
            h_control_btag_weight_notag_udsg.Fill(jet_weight_factor)

    return jet_weight_factor, sf, eff


