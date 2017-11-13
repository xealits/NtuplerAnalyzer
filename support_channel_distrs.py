import os
from os import environ
from array import array
from collections import OrderedDict
import cProfile
import logging as Logging


logging = Logging.getLogger("common")

logging.info('importing ROOT')
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TLorentzVector, TVector3, gROOT, gSystem, TCanvas, TGraphAsymmErrors, TMath, TString

# the lib is needed for BTagCalibrator and Recoil corrections
# TODO: somehow these 2 CMSSW classes should be acceptable from pyROOT on its' own without the whole lib
ROOT.gROOT.Reset()
ROOT.gSystem.Load("libUserCodettbar-leptons-80X.so")

pileup_ratio = array('d', [0, 0.360609416811339, 0.910848525427002, 1.20629960507795, 0.965997726573782, 1.10708082813183, 1.14843491548622, 0.786526251164482, 0.490577792661333, 0.740680941110478,
0.884048630953726, 0.964813189764159, 1.07045369167689, 1.12497267309738, 1.17367530613108, 1.20239808206413, 1.20815108390021, 1.20049333094509, 1.18284686347315, 1.14408796655615,
1.0962284704313, 1.06549162803223, 1.05151011089581, 1.05159666626121, 1.05064452078328, 1.0491726301522, 1.05772537082991, 1.07279673875566, 1.0837536468865, 1.09536667397119,
1.10934472980173, 1.09375894592864, 1.08263679568271, 1.04289345879947, 0.9851490341672, 0.909983816540809, 0.821346330143864, 0.71704523475871, 0.609800913869359, 0.502935245638477,
0.405579825620816, 0.309696044611377, 0.228191137503131, 0.163380359253309, 0.113368437957202, 0.0772279997453792, 0.0508111733313502, 0.0319007262683943, 0.0200879459309245, 0.0122753366005436,
0.00739933885813127, 0.00437426967257811, 0.00260473545284139, 0.00157047254226743, 0.000969500595715493, 0.000733193118123283, 0.000669817107713128, 0.000728548958604492, 0.000934559691182011, 0.00133719688378802,
0.00186652283903214, 0.00314422244976771, 0.00406954793369611, 0.00467888840511915, 0.00505224284441512, 0.00562827194936864, 0.0055889504870752, 0.00522867039470319, 0.00450752163476433, 0.00395300774604375,
0.00330577167682956, 0.00308353042577215, 0.00277846504893301, 0.00223943190687725, 0.00196650068765464, 0.00184742734258922,])
# 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

pileup_ratio_up = array('d', [0, 0.351377216124927, 0.717199649125846, 1.14121536968772, 0.84885826611733, 1.00700929402897, 1.03428595270903, 0.717444379696992, 0.344078389355127, 0.499570875027422,
0.606614916257104, 0.632584599390169, 0.731450949466174, 0.827511723989754, 0.910682115553867, 0.960170981598162, 0.988896170761361, 1.02468865580207, 1.05296667126403, 1.05112033565679,
1.0269129153969, 1.00548641752714, 0.998316130432865, 1.01492587998551, 1.03753749807849, 1.05742218946485, 1.08503978097083, 1.12134132247053, 1.15585339474274, 1.19214399856171,
1.23308400947467, 1.24528804633732, 1.26786364716917, 1.26101551498967, 1.23297806722714, 1.18042533075471, 1.10534683838101, 1.00275591661645, 0.889094305531985, 0.768791254270252,
0.655054015673457, 0.533361034358457, 0.423095146361996, 0.329177839117034, 0.250352385505809, 0.188377378855567, 0.137852651411779, 0.0968577167707531, 0.0686240187247059, 0.0473889635126706,
0.0323695027438475, 0.0216752397914914, 0.0145119352923332, 0.00961177893634792, 0.00615582219138384, 0.00430085627914427, 0.00305735512896403, 0.00223567790438986, 0.00189369737638594, 0.00199585978316291,
0.00236236592656064, 0.00372472999463276, 0.00474687312579969, 0.00549508151576102, 0.00603023110946686, 0.0068545111910253, 0.00695838760530896, 0.00666224781277046, 0.00588243140681038, 0.00528714370892014,
0.00453424615273565, 0.00433985030329723, 0.00401493171035719, 0.00332436608713241, 0.00300063798808221, 0.00289925128977536,])
#0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

pileup_ratio_down = array('d', [0, 0.37361294640242, 1.1627791004568, 1.26890787896295, 1.10266790442705, 1.23456697093644, 1.26278991594152, 0.909648777562084, 0.759569490571151, 1.09035651921682,
1.34530547603283, 1.48713160105, 1.52535976889483, 1.49730550773404, 1.49792998045778, 1.49767851097519, 1.44431045398336, 1.3681909492045, 1.29912252494785, 1.2274279217797,
1.16525969099909, 1.12531044676724, 1.09094501417685, 1.06405434433422, 1.03997120824565, 1.0185716022098, 1.00560949501652, 0.997570939806059, 0.985543761409897, 0.972557804582185,
0.957832827239337, 0.9139572640153, 0.872252387173971, 0.808388185417578, 0.733817960498049, 0.650440963845892, 0.561688505024782, 0.466564380334112, 0.374428618658619, 0.28845274688129,
0.214909665968644, 0.149991974352384, 0.100014138338029, 0.0642260884603397, 0.0396553405911344, 0.0238687936736627, 0.0137921542898078, 0.00756854010632403, 0.00415483516246187, 0.00221776872027937,
0.00118249725637452, 0.000641889697310868, 0.000383647166012176, 0.000273637590071334, 0.000242902582071058, 0.000291239677209452, 0.000394091114279828, 0.000542541231466254, 0.000771067920964491, 0.00113596447675107,
0.00158061353194779, 0.00261959852500539, 0.00331800452823827, 0.00372426930370732, 0.00392086545082614, 0.00425479965493548, 0.00411256966391362, 0.00374240422174387, 0.00313603438166934, 0.00267155793176928,
0.00216878198028599, 0.00196249821290853, 0.00171433839159669, 0.00133866519755926, 0.00113810604240254, 0.00103447940224886,])
#0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0])









logging.info("leptonic SFs")

muon_effs_dirname = "${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/analysis/muon-effs/"
gSystem.ExpandPathName(muon_effs_dirname    )
#TString muon_effs_dirname = "analysis/muon-effs/";
logging.info("muon SFs from " + muon_effs_dirname)

muon_effs_tracking_BCDEF_file  = TFile(muon_effs_dirname + "/2016_23Sep_tracking_more_BCDEF_fits.root" )
muon_effs_tracking_GH_file     = TFile(muon_effs_dirname + "/2016_23Sep_tracking_more_GH_fits.root" )
muon_effs_tracking_BCDEF_graph = muon_effs_tracking_BCDEF_file.Get("ratio_eff_aeta_dr030e030_corr")
muon_effs_tracking_GH_graph    = muon_effs_tracking_GH_file.Get("ratio_eff_aeta_dr030e030_corr")

print type(muon_effs_tracking_BCDEF_graph)

muon_effs_id_BCDEF_file  = TFile(muon_effs_dirname + "/2016_23Sep_MuonID_EfficienciesAndSF_BCDEF.root" )
muon_effs_id_GH_file     = TFile(muon_effs_dirname + "/2016_23Sep_MuonID_EfficienciesAndSF_GH.root" )
muon_effs_id_BCDEF_histo = muon_effs_id_BCDEF_file.Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta").Get("pt_abseta_ratio")
muon_effs_id_GH_histo    = muon_effs_id_GH_file   .Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta").Get("pt_abseta_ratio")

print type(muon_effs_id_BCDEF_histo)

muon_effs_iso_BCDEF_file  = TFile(muon_effs_dirname + "/2016_23Sep_MuonISO_EfficienciesAndSF_BCDEF.root" )
muon_effs_iso_GH_file     = TFile(muon_effs_dirname + "/2016_23Sep_MuonISO_EfficienciesAndSF_GH.root" )
muon_effs_iso_BCDEF_histo = muon_effs_iso_BCDEF_file.Get("TightISO_TightID_pt_eta").Get("pt_abseta_ratio")
muon_effs_iso_GH_histo    = muon_effs_iso_GH_file   .Get("TightISO_TightID_pt_eta").Get("pt_abseta_ratio")

print type(muon_effs_iso_BCDEF_histo)

# --- yep, everywhere here Tight ID and ISO is used, since that's the leptons I use

# ------------------------------------------------------------------
# Muon trigger SF
# 

muon_effs_trg_BCDEF_file  = TFile(muon_effs_dirname + "/2016_23Sep_SingleMuonTrigger_EfficienciesAndSF_RunBtoF.root" )
muon_effs_trg_GH_file     = TFile(muon_effs_dirname + "/2016_23Sep_SingleMuonTrigger_EfficienciesAndSF_Period4.root" )
muon_effs_trg_BCDEF_histo = muon_effs_trg_BCDEF_file.Get("IsoMu24_OR_IsoTkMu24_PtEtaBins").Get("pt_abseta_ratio")
muon_effs_trg_GH_histo    = muon_effs_trg_GH_file   .Get("IsoMu24_OR_IsoTkMu24_PtEtaBins").Get("pt_abseta_ratio")

print type(muon_effs_trg_BCDEF_file)

muon_effs_id_BCDEF_histo_max_x = muon_effs_id_BCDEF_histo.GetXaxis().GetXmax()
muon_effs_id_BCDEF_histo_max_y = muon_effs_id_BCDEF_histo.GetYaxis().GetXmax()
muon_effs_iso_BCDEF_histo_max_x = muon_effs_iso_BCDEF_histo.GetXaxis().GetXmax()
muon_effs_iso_BCDEF_histo_max_y = muon_effs_iso_BCDEF_histo.GetYaxis().GetXmax()

print muon_effs_id_BCDEF_histo_max_x, muon_effs_id_BCDEF_histo_max_y, muon_effs_iso_BCDEF_histo_max_x, muon_effs_iso_BCDEF_histo_max_y

muon_effs_id_GH_histo_max_x = muon_effs_id_GH_histo.GetXaxis().GetXmax()
muon_effs_id_GH_histo_max_y = muon_effs_id_GH_histo.GetYaxis().GetXmax()
muon_effs_iso_GH_histo_max_x = muon_effs_iso_GH_histo.GetXaxis().GetXmax()
muon_effs_iso_GH_histo_max_y = muon_effs_iso_GH_histo.GetYaxis().GetXmax()

print muon_effs_id_GH_histo_max_x, muon_effs_id_GH_histo_max_y, muon_effs_iso_GH_histo_max_x, muon_effs_iso_GH_histo_max_y

h_weight_mu_trk_bcdef_pt = TH1D("weight_mu_trk_bcdef_pt", "", 50, 0, 200)
h_weight_mu_idd_bcdef_pt = TH1D("weight_mu_idd_bcdef_pt", "", 50, 0, 200)
h_weight_mu_iso_bcdef_pt = TH1D("weight_mu_iso_bcdef_pt", "", 50, 0, 200)
h_weight_mu_trg_bcdef_pt = TH1D("weight_mu_trg_bcdef_pt", "", 50, 0, 200)

h_weight_mu_trk_bcdef_eta = TH1D("weight_mu_trk_bcdef_eta", "", 50, 0, 3)
h_weight_mu_idd_bcdef_eta = TH1D("weight_mu_idd_bcdef_eta", "", 50, 0, 3)
h_weight_mu_iso_bcdef_eta = TH1D("weight_mu_iso_bcdef_eta", "", 50, 0, 3)
h_weight_mu_trg_bcdef_eta = TH1D("weight_mu_trg_bcdef_eta", "", 50, 0, 3)

#
#h_weight_mu_idd_bcdef_nbins = TH2D("weight_mu_idd_bcdef_nbins", "", 100, 0, muon_effs_id_BCDEF_histo.GetSize(), 100, 0, muon_effs_id_BCDEF_histo.GetSize())

def lepton_muon_SF(abs_eta, pt): #, SingleMuon_data_bcdef_fraction, SingleMuon_data_gh_fraction):
    #weight *= 0.98; // average mu trig SF
    #weight *= 0.97; // average mu ID
    weight_muon_sfs = 1

    # muon ID SFs:
    #double weight_muon_sfs = 1;
    bcdef_weight = 1
    gh_weight = 1
    #double SingleMuon_data_bcdef_fraction = 19716.274 / (19716.274 + 15931.028);
    #double SingleMuon_data_gh_fraction    = 15931.028 / (19716.274 + 15931.028);

    #double abs_eta = abs(lep0_eta);
    #double pt = lep0_pt;

    # for control (TODO: remove this later)
    #double bcdef_weight_trk, bcdef_weight_id, bcdef_weight_iso,
    #        gh_weight_trk, gh_weight_id, gh_weight_iso;

    # hopefully tracking won't overflow in eta range:
    bcdef_weight_trk = muon_effs_tracking_BCDEF_graph.Eval(abs_eta)
    #h_weight_mu_trk_bcdef_pt = TH1D("weight_mu_trk_bcdef_pt", "", 50, 0, 200),
    h_weight_mu_trk_bcdef_eta.Fill(abs_eta)
    # the id-s totally can overflow:
    bin_x = pt      if pt < muon_effs_id_BCDEF_histo_max_x      else muon_effs_id_BCDEF_histo_max_x - 1
    bin_y = abs_eta if abs_eta < muon_effs_id_BCDEF_histo_max_y else muon_effs_id_BCDEF_histo_max_y - 0.01 # checked. the binnint is about 0.2 there
    bcdef_weight_id = muon_effs_id_BCDEF_histo.GetBinContent(muon_effs_id_BCDEF_histo.FindBin(bin_x, bin_y))
    h_weight_mu_idd_bcdef_pt .Fill(bin_x)
    h_weight_mu_idd_bcdef_eta.Fill(bin_y)

    # these too:
    bin_x = pt      if pt      < muon_effs_iso_BCDEF_histo_max_x else muon_effs_iso_BCDEF_histo_max_x - 1
    bin_y = abs_eta if abs_eta < muon_effs_iso_BCDEF_histo_max_y else muon_effs_iso_BCDEF_histo_max_y - 0.01
    bcdef_weight_iso = muon_effs_iso_BCDEF_histo.GetBinContent(muon_effs_iso_BCDEF_histo.FindBin(bin_x, bin_y))
    h_weight_mu_iso_bcdef_pt .Fill(bin_x)
    h_weight_mu_iso_bcdef_eta.Fill(bin_y)
    #bin_x = (pt < muon_effs_iso_BCDEF_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_iso_BCDEF_histo->GetXaxis()->GetXmax() - 1);
    #bin_y = (abs_eta < muon_effs_iso_BCDEF_histo->GetYaxis()->GetXmax() ? abs_eta : muon_effs_iso_BCDEF_histo->GetYaxis()->GetXmax() - 1);
    #bcdef_weight_iso = muon_effs_iso_BCDEF_histo->GetBinContent (muon_effs_iso_BCDEF_histo->FindBin(bin_x, bin_y));

    #fill_1d(string("weight_muon_effs_BCDEF_trk"),  200, 0., 1.1,   bcdef_weight_trk, 1);
    #fill_1d(string("weight_muon_effs_BCDEF_id"),   200, 0., 1.1,   bcdef_weight_id,  1);
    #fill_1d(string("weight_muon_effs_BCDEF_iso"),  200, 0., 1.1,   bcdef_weight_iso, 1);
    #bcdef_weight *= bcdef_weight_trk * bcdef_weight_id * bcdef_weight_iso;

    gh_weight_trk = muon_effs_tracking_GH_graph.Eval(abs_eta)
    #h_weight_mu_trk_gh_eta.Fill(abs_eta)
    # the id-s totally can overflow:
    bin_x = pt      if pt < muon_effs_id_GH_histo_max_x      else muon_effs_id_GH_histo_max_x - 1
    bin_y = abs_eta if abs_eta < muon_effs_id_GH_histo_max_y else muon_effs_id_GH_histo_max_y - 0.01
    gh_weight_id = muon_effs_id_GH_histo.GetBinContent(muon_effs_id_GH_histo.FindBin(bin_x, bin_y))
    #h_weight_mu_idd_gh_pt .Fill(bin_x)
    #h_weight_mu_idd_gh_eta.Fill(bin_y)
    # these too:
    bin_x = pt      if pt < muon_effs_iso_GH_histo_max_x      else muon_effs_iso_GH_histo_max_x - 1
    bin_y = abs_eta if abs_eta < muon_effs_iso_GH_histo_max_y else muon_effs_iso_GH_histo_max_y - 0.01
    gh_weight_iso = muon_effs_iso_GH_histo.GetBinContent(muon_effs_iso_GH_histo.FindBin(bin_x, bin_y))
    #h_weight_mu_iso_gh_pt .Fill(bin_x)
    #h_weight_mu_iso_gh_eta.Fill(bin_y)

    return bcdef_weight_trk, bcdef_weight_id, bcdef_weight_iso, gh_weight_trk, gh_weight_id, gh_weight_iso


muon_effs_trg_BCDEF_histo_max_x = muon_effs_trg_BCDEF_histo.GetXaxis().GetXmax()
muon_effs_trg_BCDEF_histo_max_y = muon_effs_trg_BCDEF_histo.GetYaxis().GetXmax()

print muon_effs_trg_BCDEF_histo_max_x, muon_effs_trg_BCDEF_histo_max_y

# MUON Trigger
#double lepton_muon_trigger_SF ()
def lepton_muon_trigger_SF(abs_eta, pt): #, double SingleMuon_data_bcdef_fraction, double SingleMuon_data_gh_fraction)
    no_mu_trig = 1;
    #double mu_trig_weight = 1;
    # calculate it the inverse-probbility way
    #double abs_eta = abs(NT_lep_eta_0);
    #double pt = NT_lep_pt_0;
    bin_x = pt      if pt      < muon_effs_trg_BCDEF_histo_max_x else  muon_effs_trg_BCDEF_histo_max_x - 1
    bin_y = abs_eta if abs_eta < muon_effs_trg_BCDEF_histo_max_y else  muon_effs_trg_BCDEF_histo_max_y - 0.01
    #no_mu_trig *= SingleMuon_data_bcdef_fraction * (1 - muon_effs_trg_BCDEF_histo->GetBinContent( muon_effs_trg_BCDEF_histo->FindBin(bin_x, bin_y) )) +
    #        SingleMuon_data_gh_fraction * (1 - muon_effs_trg_GH_histo->GetBinContent( muon_effs_trg_GH_histo->FindBin(bin_x, bin_y) ));
    #mu_trig_weight = 1 - no_trig; // so for 1 muon it will = to the SF, for 2 there will be a mix
    #fill_1d(string("weight_trigger_no_muon"),  200, 0., 1.1,   no_mu_trig, 1);

    h_weight_mu_trg_bcdef_pt .Fill(bin_x)
    h_weight_mu_trg_bcdef_eta.Fill(bin_y)

    #weight_muon_trig = 1 - no_mu_trig;
    return muon_effs_trg_BCDEF_histo.GetBinContent(muon_effs_trg_BCDEF_histo.FindBin(bin_x, bin_y)), muon_effs_trg_GH_histo.GetBinContent(muon_effs_trg_GH_histo.FindBin(bin_x, bin_y))


'''
/* ---------------------------------------------------
 * now, electrons have
 *      track(reconstruction) efficiency, which is recommended per eta of muon now (however there should be something about N vertices too..
 *      and ID sf
 *      also trigger
 *
 * the trig eff for dilepton case is: apply negative of it for both leptons
 */
'''
logging.info("unpacking electron eff SFs")

electron_effs_dirname = "${CMSSW_BASE}/src/UserCode/ttbar-leptons-80X/analysis/electron-effs"
gSystem.ExpandPathName(electron_effs_dirname)

electron_effs_tracking_all_file  = TFile(electron_effs_dirname + "/2016_Sept23_ElectronReconstructionSF_egammaEffi.txt_EGM2D.root")
electron_effs_tracking_all_histo = electron_effs_tracking_all_file.Get("EGamma_SF2D")
logging.info("Y tracking (reconstruction)")

# for the selected electrons, Tight ID
# not for Veto
electron_effs_id_all_file  = TFile(electron_effs_dirname + "/2016_Sept23_ElectronID_TightCutBased_egammaEffi.txt_EGM2D.root")
electron_effs_id_all_histo = electron_effs_id_all_file.Get("EGamma_SF2D")
logging.info("Y id")

#analysis/electron-effs/2016_03Feb_TriggerSF_Run2016All_v1.root
electron_effs_trg_all_file  = TFile(electron_effs_dirname + "/2016_03Feb_TriggerSF_Run2016All_v1.root")
electron_effs_trg_all_histo = electron_effs_trg_all_file.Get("Ele27_WPTight_Gsf")
logging.info("Y trigger")

# --- these SFs will be applied to the selected leptons independently

electron_effs_tracking_all_histo_max_y = electron_effs_tracking_all_histo.GetYaxis().GetXmax()
electron_effs_id_all_histo_max_y = electron_effs_id_all_histo.GetYaxis().GetXmax()

h_weight_el_trk_pt = TH1D("weight_mu_trk_bcdef_pt", "", 50, 0, 200)
h_weight_el_idd_pt = TH1D("weight_mu_idd_bcdef_pt", "", 50, 0, 200)
h_weight_mu_iso_bcdef_pt = TH1D("weight_mu_iso_bcdef_pt", "", 50, 0, 200)
h_weight_mu_trg_bcdef_pt = TH1D("weight_mu_trg_bcdef_pt", "", 50, 0, 200)

def lepton_electron_SF(eta, pt):
    #double weight_reco, weight_id;

    # here X axis is eta, Y axis is pt
    # X is from -2.5 to 2.5 -- our eta is up to 2.4 (2.5 in wide ntuples), should be ok
    #double bin_x = (pt < electron_effs_tracking_all_histo->GetXaxis()->GetXmax()      ? pt : muon_effs_id_BCDEF_histo->GetXaxis()->GetXmax() - 1);
    bin_x = eta;
    bin_y = pt if pt < electron_effs_tracking_all_histo_max_y else electron_effs_tracking_all_histo_max_y - 1
    sf_reco = electron_effs_tracking_all_histo.GetBinContent (electron_effs_tracking_all_histo.FindBin(bin_x, bin_y))

    #bin_x = eta;
    bin_y = pt if pt < electron_effs_id_all_histo_max_y else electron_effs_id_all_histo_max_y - 1
    sf_id = electron_effs_id_all_histo.GetBinContent (electron_effs_id_all_histo.FindBin(bin_x, bin_y))

    return sf_reco, sf_id

electron_effs_trg_all_histo_max_x = electron_effs_trg_all_histo.GetXaxis().GetXmax()
electron_effs_trg_all_histo_max_y = electron_effs_trg_all_histo.GetYaxis().GetXmax()

def lepton_electron_trigger_SF(eta, pt):
    #double no_ele_trig = 1;
    # (calculate it the inverse-probbility way)
    # no, just return the SF, assume 1-lepton case
    #pat::Electron& el = selElectrons[i];
    #double eta = el.superCluster()->position().eta();
    #double pt = el.pt();
    # here X axis is pt, Y axis is eta (from -2.5 to 2.5)
    bin_x = pt  if pt  < electron_effs_trg_all_histo_max_x else electron_effs_trg_all_histo_max_x - 1

    bin_y = eta
    if eta > electron_effs_trg_all_histo_max_y:
        bin_y = electron_effs_trg_all_histo_max_y - 0.01
    elif eta < - electron_effs_trg_all_histo_max_y:
        bin_y = - electron_effs_trg_all_histo_max_y + 0.01

    return electron_effs_trg_all_histo.GetBinContent(electron_effs_trg_all_histo.FindBin(bin_x, bin_y))
    #el_trig_weight = 1 - no_trig; // so for 1 lepton it will = to the SF, for 2 there will be a mix
    #fill_1d(string("weight_trigger_no_electron"),  200, 0., 1.1,   no_ele_trig, 1);




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


with_bSF = True

if with_bSF:
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

    return jet_weight_factor



def top_pT_SF(x):
    # the SF function is SF(x)=exp(a+bx)
    # where x is pT of the top quark (at generation?)
    # sqrt(s)         channel             a             b
    # 7 TeV         all combined         0.199         -0.00166
    # 7 TeV         l+jets              0.174         -0.00137
    # 7 TeV         dilepton            0.222         -0.00197
    # 8 TeV         all combined         0.156         -0.00137
    # 8 TeV         l+jets               0.159         -0.00141
    # 8 TeV         dilepton             0.148         -0.00129
    # 13 TeV        all combined        0.0615        -0.0005
    # -- taking all combined 13 TeV
    a = 0.0615
    b = -0.0005
    return TMath.Exp(a + b*x)

def ttbar_pT_SF(t_pt, tbar_pt):
    return TMath.Sqrt(top_pT_SF(t_pt) * top_pT_SF(tbar_pt))

def transverse_mass(v1, v2):
    return TMath.Sqrt(2*v1.pt()*v2.pt()*(1 - TMath.Cos(v1.phi() - v2.phi())))

def transverse_mass_pts(v1_x, v1_y, v2_x, v2_y):
    v1v2 = TMath.Sqrt((v1_x*v1_x + v1_y*v1_y)*(v2_x*v2_x + v2_y*v2_y))
    return TMath.Sqrt(2*(v1v2 - (v1_x*v2_x + v1_y*v2_y)))


#lj_var = calc_lj_var(jets, jets_b)
def calc_lj_var(light_jets, b_jets):
    if len(b_jets) == 0 or len(light_jets) < 2: return 5000., 5000., 5000.
    # (not doing it now) loop over light jets -- check their mass
    # loop over light jets, check mass of their pairs to be close to 80
    # find the closest vector
    # then loop over b jets finding the combination with mass closest to 173

    # light jets close to W
    dist_W = 99999.
    ## closest_vector
    #for j, mult in light_jets:
    #    new_dist = abs(j.mass() * mult - 80) # DANGER: the LorentzVector is used, .M() is spacial magnitude
    #    if new_dist < dist_W:
    #        dist_W = new_dist
    #        closest_to_W = j * mult

    # pairs of light jets
    for i in range(len(light_jets)):
      for u in range(i):
        ji, multi = light_jets[i]
        ju, multu = light_jets[u]
        pair = ji * multi + ju * multu
        new_dist = abs(pair.mass() - 80)
        if new_dist < dist_W:
            dist_W = new_dist
            closest_to_W = pair

    # closest to 173
    dist_t = 99999.
    for j, mult in b_jets:
        pair = j * mult + closest_to_W
        new_dist = abs(pair.mass() - 173)
        if new_dist < dist_t:
            dist_t = new_dist
            closest_to_t = pair

    return TMath.Sqrt(dist_W*dist_W + dist_t*dist_t), closest_to_W.mass(), closest_to_t.mass()


# no data types/protocols in ROOT -- looping has to be done manually
def full_loop(tree, dtag, lumi_bcdef, lumi_gh, range_min, range_max, logger):
    '''full_loop(tree, dtag)

    TTree tree
    dtag string
    '''

    ratio_bcdef = lumi_bcdef / (lumi_bcdef + lumi_gh)
    ratio_gh    = lumi_gh / (lumi_bcdef + lumi_gh)

    logger.write('dtag=%s\n' % dtag)
    logger.write('lumi_bcdef=%f lumi_gh=%f\n' % (lumi_bcdef, lumi_gh))

    isMC = 'MC' in dtag
    #save_control = False
    save_weights = True and isMC
    aMCatNLO = 'amcatnlo' in dtag
    isTT = 'TT' in dtag
    isSTop = 'SingleT' in dtag or 'tchannel' in dtag or 'schannel' in dtag
    isSTopTSchannels = 'tchannel' in dtag or 'schannel' in dtag
    isDY = 'DY' in dtag
    isWJets = 'WJet' in dtag or 'W1Jet' in dtag or 'W2Jet' in dtag or 'W3Jet' in dtag or 'W4Jet' in dtag
    isQCD = 'QCD' in dtag
    isDibosons = 'WW' in dtag or 'ZZ' in dtag or 'WZ' in dtag

    logger.write(' '.join(name + '=' + str(setting) for name, setting in
        [('isMC', isMC), ('save_weights', save_weights), ('aMCatNLO', aMCatNLO), ('isTT', isTT), ('isSTop', isSTop), ('isSTopTSchannels', isSTopTSchannels), ('isDY', isDY), ('isWJets', isWJets), ('isQCD', isQCD), ('isDibosons', isDibosons)]) + '\n')


    logging.info("load root_funcs")
    ROOT.gROOT.ProcessLine(".L /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/root_funcs.C+") # TODO: change absolute path to relative

    # Recoil corrections
    doRecoilCorrections = isWJets or isDY # TODO: check if it is needed for DY
    if doRecoilCorrections:
        logging.info("will use Recoil Corrections")
        ROOT.gROOT.ProcessLine(".L /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/recoil_corrections.C+") # TODO: change absolute path to relative
        #recoil_corrections_data_file = TString("/HTT-utilities/RecoilCorrections/data/TypeIPFMET_2016BCD.root")
        #recoilPFMetCorrector = ROOT.RecoilCorrector(recoil_corrections_data_file);
        #ROOT.BTagCalibrationReader


    # Z pt mass weight
    # std::string zPtMassWeights_filename = std::string(std::getenv("CMSSW_BASE")) + "/src/UserCode/zpt_weights_2016.root";
    # TFile* zPtMassWeights_file  = TFile::Open(zPtMassWeights_filename.c_str());
    # TH2D*  zPtMassWeights_histo = (TH2D*) zPtMassWeights_file->Get("zptmass_histo");
    # TH2D*  zPtMassWeights_histo_err = (TH2D*) zPtMassWeights_file->Get("zptmass_histo_err");
    # 
    # float zPtMass_weight(float genMass, float genPt)
    #         {
    #         return zPtMassWeights_histo->GetBinContent(zPtMassWeights_histo->GetXaxis()->FindBin(genMass), zPtMassWeights_histo->GetYaxis()->FindBin(genPt));
    #         }

    if isDY:
        logging.info("loading Z pt mass weights")
        zPtMass_filename = environ['CMSSW_BASE'] + '/src/UserCode/zpt_weights_2016.root'
        zPtMassWeights_file  = TFile(zPtMass_filename)
        zPtMassWeights_file.Print()
        zPtMassWeights_histo     = zPtMassWeights_file.Get("zptmass_histo")
        zPtMassWeights_histo_err = zPtMassWeights_file.Get("zptmass_histo_err")

        def zPtMass_weight(genMass, genPt):
            return zPtMassWeights_histo.GetBinContent(zPtMassWeights_histo.GetXaxis().FindBin(genMass), zPtMassWeights_histo.GetYaxis().FindBin(genPt))


    #set_bSF_effs_for_dtag(dtag)
    if with_bSF: logger.write(' '.join(str(id(h)) for h in (bEff_histo_b, bEff_histo_c, bEff_histo_udsg)) + '\n')
    #print b_alljet, b_tagged, c_alljet, c_tagged, udsg_alljet, udsg_tagged
    #global bTagging_b_jet_efficiency, bTagging_c_jet_efficiency, bTagging_udsg_jet_efficiency
    #bTagging_b_jet_efficiency = bTagging_X_jet_efficiency(b_alljet, b_tagged)
    #bTagging_c_jet_efficiency = bTagging_X_jet_efficiency(c_alljet, c_tagged)
    #bTagging_udsg_jet_efficiency = bTagging_X_jet_efficiency(udsg_alljet, udsg_tagged)

    control_hs = OrderedDict([
    ('weight_pu', TH1D("weight_pu", "", 50, 0, 2)),
    ('weight_pu_up', TH1D("weight_pu_up", "", 50, 0, 2)),
    ('weight_pu_dn', TH1D("weight_pu_dn", "", 50, 0, 2)),
    ('weight_top_pt', TH1D("weight_top_pt", "", 50, 0, 2)),

    ('weight_z_mass_pt', TH1D("weight_z_mass_pt", "", 50, 0, 2)),
    ('weight_bSF',       TH1D("weight_bSF", "", 50, 0, 2)),

    ('weight_mu_trk_bcdef', TH1D("weight_mu_trk_bcdef", "", 50, 0, 2)),
    ('weight_mu_id_bcdef' , TH1D("weight_mu_id_bcdef", "", 50, 0, 2)),
    ('weight_mu_iso_bcdef', TH1D("weight_mu_iso_bcdef", "", 50, 0, 2)),
    ('weight_mu_trg_bcdef', TH1D("weight_mu_trg_bcdef", "", 50, 0, 2)),
    ('weight_mu_all_bcdef', TH1D("weight_mu_all_bcdef", "", 50, 0, 2)),

    ('weight_mu_trk_gh', TH1D("weight_mu_trk_gh", "", 50, 0, 2)),
    ('weight_mu_id_gh' , TH1D("weight_mu_id_gh", "", 50, 0, 2)),
    ('weight_mu_iso_gh', TH1D("weight_mu_iso_gh", "", 50, 0, 2)),
    ('weight_mu_trg_gh', TH1D("weight_mu_trg_gh", "", 50, 0, 2)),
    ('weight_mu_all_gh', TH1D("weight_mu_all_gh", "", 50, 0, 2)),

    ('weight_mu_bSF', TH1D("weight_mu_bSF", "", 50, 0, 2)),
    ('weight_mu_bSF_up', TH1D("weight_mu_bSF_up", "", 50, 0, 2)),
    ('weight_mu_bSF_down', TH1D("weight_mu_bSF_down", "", 50, 0, 2)),

    ('weight_el_trk', TH1D("weight_el_trk", "", 50, 0, 2)),
    ('weight_el_idd', TH1D("weight_el_idd", "", 50, 0, 2)),
    ('weight_el_trg', TH1D("weight_el_trg", "", 50, 0, 2)),
    ('weight_el_all', TH1D("weight_el_all", "", 50, 0, 2)),

    ('weight_el_bSF', TH1D("weight_el_bSF", "", 50, 0, 2)),
    ('weight_el_bSF_up', TH1D("weight_el_bSF_up", "", 50, 0, 2)),
    ('weight_el_bSF_down', TH1D("weight_el_bSF_down", "", 50, 0, 2)),
    ])

    # strange, getting PyROOT_NoneObjects from these after output
    for _, h in control_hs.items():
        h.SetDirectory(0)

    #channels = {'el': ['foo'], 'mu': ['foo'], 'mujets': ['foo'],
    #'mujets_b': ['foo'], 'taumutauh': ['foo'], 'taumutauh_antiMt_pretau_allb': ['foo'], 'taumutauh_antiMt_pretau': ['foo'], 'taumutauh_antiMt': ['foo'], 'taumutauh_antiMt_OS': ['foo']}
    #channels = {'presel': ['foo'], 'el': ['foo'], 'mu': ['foo']}
    # TODO: I actually need:
    #    same set of selections for el and mu
    #    which is preselection (no tau requirement), selection of tau, bins of inside/outside lj and also outside lj 2 or 1 b-tagged jets
    # do simplest, separate things
    #    Data or MC -> process (for whole event)
    #    for each systematic -> list of passed channels
    #    process+passed channel -> processes to store (TODO: or tt_other)
    # again, do simple: for each kind of dtags make {channels: ([procs], default)
    # or use a defaultdict somewhere?
    # --- I'll just do a special care for TT
    # the others have 1 process per whole dtag for now
    # (split DY and single-top later

    tt_el_procs = ['tt_eltau', 'tt_lj', 'tt_taultauh', 'tt_other']
    tt_mu_procs = ['tt_mutau', 'tt_lj', 'tt_taultauh', 'tt_other']

    if isMC:
        systematic_names_all = ['NOMINAL',
                'JESUp'     ,
                'JESDown'   ,
                'JERUp'     ,
                'JERDown'   ,
                'TauESUp'   ,
                'TauESDown' ,
                'bSFUp'     ,
                'bSFDown'   ,
                'PUUp'      ,
                'PUDown'    ,
                ]
        systematic_names_pu = ['NOMINAL', 'PUUp', 'PUDown']
        systematic_names_pu_toppt = ['NOMINAL', 'PUUp', 'PUDown']
    else:
        systematic_names_all = ['NOMINAL']
        systematic_names_pu  = ['NOMINAL']
        systematic_names_pu_toppt  = ['NOMINAL']
    systematic_names_toppt = ['NOMINAL']

    if isTT:
        systematic_names_all.append('TOPPTUp')
        systematic_names_all.append('TOPPTDown')
        systematic_names_toppt = ['NOMINAL', 'TOPPTUp']
        systematic_names_pu_toppt.append('TOPPTUp')

    if isMC:
        if isTT:
            tt_procs_el =  (['tt_eltau', 'tt_lj', 'tt_taultauh', 'tt_other'], 'tt_other')
            tt_procs_mu =  (['tt_mutau', 'tt_lj', 'tt_taultauh', 'tt_other'], 'tt_other')
            channels = {'el_presel': (tt_procs_el, systematic_names_all),
                   'el_sel':         (tt_procs_el, systematic_names_all),
                   'el_lj':          (tt_procs_el, systematic_names_all),
                   'el_lj_out':      (tt_procs_el, systematic_names_all),
                   'mu_presel':      (tt_procs_mu, systematic_names_all),
                   'mu_sel':         (tt_procs_mu, systematic_names_all),
                   'mu_lj':          (tt_procs_mu, systematic_names_all),
                   'mu_lj_out':      (tt_procs_mu, systematic_names_all),
                   # same sign for some QCD control
                   'el_sel_ss':      (tt_procs_el, ['NOMINAL']),
                   'mu_sel_ss':      (tt_procs_mu, ['NOMINAL']),
                   'el_presel_ss':   (tt_procs_el, ['NOMINAL']),
                   'mu_presel_ss':   (tt_procs_mu, ['NOMINAL']),
                   # with tau POG selection
                   'pog_mu_presel':  (tt_procs_mu, systematic_names_toppt),
                   'pog_mu_pass':    (tt_procs_mu, systematic_names_toppt),
                   'pog_mu_pass_ss': (tt_procs_mu, systematic_names_toppt),
                   'pog_mu_fail':    (tt_procs_mu, systematic_names_toppt),
                   # with addition of no DY mass, no match to b-tag (could add a cut on small MT)
                   'adv_el_sel':     (tt_procs_el, systematic_names_pu),
                   'adv_mu_sel':     (tt_procs_mu, systematic_names_pu),
                   'sel_mu_min':     (tt_procs_mu, systematic_names_pu_toppt), # minumum muon/el thresholds, loose b, loose tau
                   'sel_mu_min_ss':  (tt_procs_mu, systematic_names_pu_toppt), # minumum muon/el thresholds, loose b, loose tau
                   'sel_mu_min_medtau':     (tt_procs_mu, systematic_names_pu_toppt), # minimum selection with Medium taus -- hopefully it will reduce QCD
                   # control selections: WJets, DY mumu and tautau, tt elmu
                   'ctr_mu_wjet':    (tt_procs_mu, ['NOMINAL']),
                   'ctr_mu_wjet_old':(tt_procs_mu, ['NOMINAL']),
                   'ctr_mu_dy_mumu': (tt_procs_mu, ['NOMINAL']),
                   'ctr_mu_dy_tt':   (tt_procs_mu, ['NOMINAL']),
                   'ctr_mu_dy_tt_ss':(tt_procs_mu, ['NOMINAL']),
                   'ctr_mu_dy_SV_tt':   (tt_procs_mu, systematic_names_all),
                   'ctr_mu_dy_SV_tt_ss':(tt_procs_mu, systematic_names_all),
                   'ctr_mu_tt_em':   (tt_procs_mu, systematic_names_toppt),
                   }
            usual_process = 'tt_other'

        if isWJets:
            wjets_procs = (['wjets'], 'wjets')
            channels = {'el_presel': (wjets_procs, systematic_names_all),
                   'el_sel':      (wjets_procs, systematic_names_all),
                   'el_lj':       (wjets_procs, systematic_names_all),
                   'el_lj_out':   (wjets_procs, systematic_names_all),
                   'mu_presel':   (wjets_procs, systematic_names_all),
                   'mu_sel':      (wjets_procs, systematic_names_all),
                   'mu_lj':       (wjets_procs, systematic_names_all),
                   'mu_lj_out':   (wjets_procs, systematic_names_all),
                   # selection steps to control shape of DY and WJets
                   # current definitions are:
                   #pass_single_lep_presel = large_met and has_3jets and has_bjets #and os_lep_med_tau
                   #pass_single_lep_sel = pass_single_lep_presel and os_lep_med_tau
                   'mu_sel_nomet':         (wjets_procs, ['NOMINAL']),
                   #'el_sel_nomet':         (wjets_procs, ['NOMINAL']),
                   'mu_sel_onlymet':       (wjets_procs, ['NOMINAL']),
                   #'el_sel_onlymet':       (wjets_procs, ['NOMINAL']),
                   'mu_sel_nobtag':        (wjets_procs, ['NOMINAL']),
                   #'el_sel_nobtag':        (wjets_procs, ['NOMINAL']),
                   'mu_sel_onlybtag':      (wjets_procs, ['NOMINAL']),
                   #'el_sel_onlybtag':      (wjets_procs, ['NOMINAL']),
                   'mu_sel_no3jets':       (wjets_procs, ['NOMINAL']),
                   #'el_sel_no3jets':       (wjets_procs, ['NOMINAL']),
                   'mu_sel_only3jets':     (wjets_procs, ['NOMINAL']),
                   #'el_sel_only3jets':     (wjets_procs, ['NOMINAL']),
                   # same sign for some QCD control
                   'el_sel_ss':      (wjets_procs, ['NOMINAL']),
                   'mu_sel_ss':      (wjets_procs, ['NOMINAL']),
                   'el_presel_ss':   (wjets_procs, ['NOMINAL']),
                   'mu_presel_ss':   (wjets_procs, ['NOMINAL']),
                   # with tau POG selection
                   'pog_mu_presel':  (wjets_procs, systematic_names_toppt),
                   'pog_mu_pass':    (wjets_procs, systematic_names_toppt),
                   'pog_mu_pass_ss': (wjets_procs, systematic_names_toppt),
                   'pog_mu_fail':    (wjets_procs, systematic_names_toppt),
                   # with addition of no DY mass, no match to b-tag (could add a cut on small MT)
                   'adv_el_sel':     (wjets_procs, systematic_names_pu),
                   'adv_mu_sel':     (wjets_procs, systematic_names_pu),
                   'sel_mu_min':     (wjets_procs, systematic_names_pu_toppt),
                   'sel_mu_min_ss':  (wjets_procs, systematic_names_pu_toppt),
                   'sel_mu_min_medtau':     (wjets_procs, systematic_names_pu_toppt),
                   # control selections: WJets, DY mumu and tautau, tt elmu
                   'ctr_mu_wjet':    (wjets_procs, ['NOMINAL']),
                   'ctr_mu_wjet_old':(wjets_procs, ['NOMINAL']),
                   'ctr_mu_dy_mumu': (wjets_procs, ['NOMINAL']),
                   'ctr_mu_dy_tt':   (wjets_procs, ['NOMINAL']),
                   'ctr_mu_dy_tt_ss':(wjets_procs, ['NOMINAL']),
                   'ctr_mu_dy_SV_tt':   (wjets_procs, systematic_names_all),
                   'ctr_mu_dy_SV_tt_ss':(wjets_procs, systematic_names_all),
                   'ctr_mu_tt_em':   (wjets_procs, systematic_names_toppt),
                   }
            usual_process = 'wjets'

        if isDY:
            dy_procs = (['dy_tautau', 'dy_other'], 'dy_other')
            channels = {'el_presel': (dy_procs, systematic_names_all),
                   'el_sel':         (dy_procs, systematic_names_all),
                   'el_lj':          (dy_procs, systematic_names_all),
                   'el_lj_out':      (dy_procs, systematic_names_all),
                   'mu_presel':      (dy_procs, systematic_names_all),
                   'mu_sel':         (dy_procs, systematic_names_all),
                   'mu_lj':          (dy_procs, systematic_names_all),
                   'mu_lj_out':      (dy_procs, systematic_names_all),
                   # selection steps to control shape of DY and WJets
                   # current definitions are:
                   #pass_single_lep_presel = large_met and has_3jets and has_bjets #and os_lep_med_tau
                   #pass_single_lep_sel = pass_single_lep_presel and os_lep_med_tau
                   'mu_sel_nomet':         (dy_procs, ['NOMINAL']),
                   #'el_sel_nomet':         (dy_procs, ['NOMINAL']),
                   'mu_sel_onlymet':       (dy_procs, ['NOMINAL']),
                   #'el_sel_onlymet':       (dy_procs, ['NOMINAL']),
                   'mu_sel_nobtag':        (dy_procs, ['NOMINAL']),
                   #'el_sel_nobtag':        (dy_procs, ['NOMINAL']),
                   'mu_sel_onlybtag':      (dy_procs, ['NOMINAL']),
                   #'el_sel_onlybtag':      (dy_procs, ['NOMINAL']),
                   'mu_sel_no3jets':       (dy_procs, ['NOMINAL']),
                   #'el_sel_no3jets':       (dy_procs, ['NOMINAL']),
                   'mu_sel_only3jets':     (dy_procs, ['NOMINAL']),
                   #'el_sel_only3jets':     (dy_procs, ['NOMINAL']),
                   # same sign for some QCD control
                   'el_sel_ss':      (dy_procs, ['NOMINAL']),
                   'mu_sel_ss':      (dy_procs, ['NOMINAL']),
                   'el_presel_ss':   (dy_procs, ['NOMINAL']),
                   'mu_presel_ss':   (dy_procs, ['NOMINAL']),
                   # with tau POG selection
                   'pog_mu_presel':  (dy_procs, systematic_names_toppt),
                   'pog_mu_pass':    (dy_procs, systematic_names_toppt),
                   'pog_mu_pass_ss': (dy_procs, systematic_names_toppt),
                   'pog_mu_fail':    (dy_procs, systematic_names_toppt),
                   # with addition of no DY mass, no match to b-tag (could add a cut on small MT)
                   'adv_el_sel':     (dy_procs, systematic_names_pu),
                   'adv_mu_sel':     (dy_procs, systematic_names_pu),
                   'sel_mu_min':     (dy_procs, systematic_names_pu_toppt),
                   'sel_mu_min_ss':  (dy_procs, systematic_names_pu_toppt),
                   'sel_mu_min_medtau':     (dy_procs, systematic_names_pu_toppt),
                   # control selections: WJets, DY mumu and tautau, tt elmu
                   'ctr_mu_wjet':    (dy_procs, ['NOMINAL']),
                   'ctr_mu_wjet_old':(dy_procs, ['NOMINAL']),
                   'ctr_mu_dy_mumu': (dy_procs, ['NOMINAL']),
                   'ctr_mu_dy_tt':   (dy_procs, ['NOMINAL']),
                   'ctr_mu_dy_tt_ss':(dy_procs, ['NOMINAL']),
                   'ctr_mu_dy_SV_tt':   (dy_procs, systematic_names_all),
                   'ctr_mu_dy_SV_tt_ss':(dy_procs, systematic_names_all),
                   'ctr_mu_tt_em':   (dy_procs, systematic_names_toppt),
                   }
            usual_process = 'dy_other'

        if isSTop:
            s_top_procs_el = (['s_top_eltau', 's_top_lj', 's_top_other'], 's_top_other')
            s_top_procs_mu = (['s_top_mutau', 's_top_lj', 's_top_other'], 's_top_other')
            channels = {'el_presel': (s_top_procs_el, systematic_names_all),
                   'el_sel':         (s_top_procs_el, systematic_names_all),
                   'el_lj':          (s_top_procs_el, systematic_names_all),
                   'el_lj_out':      (s_top_procs_el, systematic_names_all),
                   'mu_presel':      (s_top_procs_mu, systematic_names_all),
                   'mu_sel':         (s_top_procs_mu, systematic_names_all),
                   'mu_lj':          (s_top_procs_mu, systematic_names_all),
                   'mu_lj_out':      (s_top_procs_mu, systematic_names_all),
                   # same sign for some QCD control
                   'el_sel_ss':      (s_top_procs_el, ['NOMINAL']),
                   'mu_sel_ss':      (s_top_procs_mu, ['NOMINAL']),
                   'el_presel_ss':   (s_top_procs_el, ['NOMINAL']),
                   'mu_presel_ss':   (s_top_procs_mu, ['NOMINAL']),
                   # with tau POG selection
                   'pog_mu_presel':  (s_top_procs_mu, systematic_names_toppt),
                   'pog_mu_pass_ss': (s_top_procs_mu, systematic_names_toppt),
                   'pog_mu_fail':    (s_top_procs_mu, systematic_names_toppt),
                   # with addition of no DY mass, no match to b-tag (could add a cut on small MT)
                   'adv_el_sel':     (s_top_procs_el, systematic_names_pu),
                   'adv_mu_sel':     (s_top_procs_mu, systematic_names_pu),
                   'sel_mu_min':     (s_top_procs_mu, systematic_names_pu_toppt),
                   'sel_mu_min_ss':  (s_top_procs_mu, systematic_names_pu_toppt),
                   'sel_mu_min_medtau':     (s_top_procs_mu, systematic_names_pu_toppt),
                   # control selections: WJets, DY mumu and tautau, tt elmu
                   'ctr_mu_wjet':    (s_top_procs_mu, ['NOMINAL']),
                   'ctr_mu_wjet_old':(s_top_procs_mu, ['NOMINAL']),
                   'ctr_mu_dy_mumu': (s_top_procs_mu, ['NOMINAL']),
                   'ctr_mu_dy_tt':   (s_top_procs_mu, ['NOMINAL']),
                   'ctr_mu_dy_tt_ss':(s_top_procs_mu, ['NOMINAL']),
                   'ctr_mu_dy_SV_tt':   (s_top_procs_mu, systematic_names_all),
                   'ctr_mu_dy_SV_tt_ss':(s_top_procs_mu, systematic_names_all),
                   'ctr_mu_tt_em':   (s_top_procs_mu, systematic_names_toppt),
                   }
            usual_process = 's_top_other'

        if isQCD:
            qcd_procs = (['qcd'], 'qcd')
            channels = {'el_presel': (qcd_procs, systematic_names_all),
                   'el_sel':         (qcd_procs, systematic_names_all),
                   'el_lj':          (qcd_procs, systematic_names_all),
                   'el_lj_out':      (qcd_procs, systematic_names_all),
                   'mu_presel':      (qcd_procs, systematic_names_all),
                   'mu_sel':         (qcd_procs, systematic_names_all),
                   'mu_lj':          (qcd_procs, systematic_names_all),
                   'mu_lj_out':      (qcd_procs, systematic_names_all),
                   # same sign for some QCD control
                   'el_sel_ss':      (qcd_procs, ['NOMINAL']),
                   'mu_sel_ss':      (qcd_procs, ['NOMINAL']),
                   'el_presel_ss':   (qcd_procs, ['NOMINAL']),
                   'mu_presel_ss':   (qcd_procs, ['NOMINAL']),
                   # with tau POG selection
                   'pog_mu_presel':  (qcd_procs, systematic_names_toppt),
                   'pog_mu_pass':    (qcd_procs, systematic_names_toppt),
                   'pog_mu_pass_ss': (qcd_procs, systematic_names_toppt),
                   'pog_mu_fail':    (qcd_procs, systematic_names_toppt),
                   # with addition of no DY mass, no match to b-tag (could add a cut on small MT)
                   'adv_el_sel':     (qcd_procs, systematic_names_pu),
                   'adv_mu_sel':     (qcd_procs, systematic_names_pu),
                   'sel_mu_min':     (qcd_procs, systematic_names_pu_toppt),
                   'sel_mu_min_ss':  (qcd_procs, systematic_names_pu_toppt),
                   'sel_mu_min_medtau':     (qcd_procs, systematic_names_pu_toppt),
                   # control selections: WJets, DY mumu and tautau, tt elmu
                   'ctr_mu_wjet':    (qcd_procs, ['NOMINAL']),
                   'ctr_mu_wjet_old':(qcd_procs, ['NOMINAL']),
                   'ctr_mu_dy_mumu': (qcd_procs, ['NOMINAL']),
                   'ctr_mu_dy_tt':   (qcd_procs, ['NOMINAL']),
                   'ctr_mu_dy_tt_ss':(qcd_procs, ['NOMINAL']),
                   'ctr_mu_dy_SV_tt':   (qcd_procs, systematic_names_all),
                   'ctr_mu_dy_SV_tt_ss':(qcd_procs, systematic_names_all),
                   'ctr_mu_tt_em':   (qcd_procs, systematic_names_toppt),
                   }
            usual_process = 'qcd'

        if isDibosons:
            dibosons_procs = (['dibosons'], 'dibosons')
            channels = {'el_presel': (dibosons_procs, systematic_names_all),
                   'el_sel':         (dibosons_procs, systematic_names_all),
                   'el_lj':          (dibosons_procs, systematic_names_all),
                   'el_lj_out':      (dibosons_procs, systematic_names_all),
                   'mu_presel':      (dibosons_procs, systematic_names_all),
                   'mu_sel':         (dibosons_procs, systematic_names_all),
                   'mu_lj':          (dibosons_procs, systematic_names_all),
                   'mu_lj_out':      (dibosons_procs, systematic_names_all),
                   # same sign for some QCD control
                   'el_sel_ss':      (dibosons_procs, ['NOMINAL']),
                   'mu_sel_ss':      (dibosons_procs, ['NOMINAL']),
                   'el_presel_ss':   (dibosons_procs, ['NOMINAL']),
                   'mu_presel_ss':   (dibosons_procs, ['NOMINAL']),
                   # with tau POG selection
                   'pog_mu_presel':  (dibosons_procs, systematic_names_toppt),
                   'pog_mu_pass':    (dibosons_procs, systematic_names_toppt),
                   'pog_mu_pass_ss': (dibosons_procs, systematic_names_toppt),
                   'pog_mu_fail':    (dibosons_procs, systematic_names_toppt),
                   # with addition of no DY mass, no match to b-tag (could add a cut on small MT)
                   'adv_el_sel':     (dibosons_procs, systematic_names_pu),
                   'adv_mu_sel':     (dibosons_procs, systematic_names_pu),
                   'sel_mu_min':     (dibosons_procs, systematic_names_pu_toppt),
                   'sel_mu_min_ss':  (dibosons_procs, systematic_names_pu_toppt),
                   'sel_mu_min_medtau':     (dibosons_procs, systematic_names_pu_toppt),
                   # control selections: WJets, DY mumu and tautau, tt elmu
                   'ctr_mu_wjet':    (dibosons_procs, ['NOMINAL']),
                   'ctr_mu_wjet_old':(dibosons_procs, ['NOMINAL']),
                   'ctr_mu_dy_mumu': (dibosons_procs, ['NOMINAL']),
                   'ctr_mu_dy_tt':   (dibosons_procs, ['NOMINAL']),
                   'ctr_mu_dy_tt_ss':(dibosons_procs, ['NOMINAL']),
                   'ctr_mu_dy_SV_tt':   (dibosons_procs, systematic_names_all),
                   'ctr_mu_dy_SV_tt_ss':(dibosons_procs, systematic_names_all),
                   'ctr_mu_tt_em':   (dibosons_procs, systematic_names_toppt),
                   }
            usual_process = 'dibosons'

    else:
        data_procs = (['data'], 'data')
        channels = {'el_presel': (data_procs, systematic_names_all),
                   'el_sel':    (data_procs, systematic_names_all),
                   'el_lj':     (data_procs, systematic_names_all),
                   'el_lj_out': (data_procs, systematic_names_all),
                   'mu_presel': (data_procs, systematic_names_all),
                   'mu_sel':    (data_procs, systematic_names_all),
                   'mu_lj':     (data_procs, systematic_names_all),
                   'mu_lj_out': (data_procs, systematic_names_all),
                   # same sign for some QCD control
                   'el_sel_ss':      (data_procs, ['NOMINAL']),
                   'mu_sel_ss':      (data_procs, ['NOMINAL']),
                   'el_presel_ss':   (data_procs, ['NOMINAL']),
                   'mu_presel_ss':   (data_procs, ['NOMINAL']),
                   # with tau POG selection
                   'pog_mu_presel':  (data_procs, systematic_names_toppt),
                   'pog_mu_pass':    (data_procs, systematic_names_toppt),
                   'pog_mu_pass_ss': (data_procs, systematic_names_toppt),
                   'pog_mu_fail':    (data_procs, systematic_names_toppt),
                   # with addition of no DY mass, no match to b-tag (could add a cut on small MT)
                   'adv_el_sel':     (data_procs, systematic_names_pu),
                   'adv_mu_sel':     (data_procs, systematic_names_pu),
                   'sel_mu_min':     (data_procs, systematic_names_pu_toppt),
                   'sel_mu_min_ss':  (data_procs, systematic_names_pu_toppt),
                   'sel_mu_min_medtau':     (data_procs, systematic_names_pu_toppt),
                   # control selections: WJets, DY mumu and tautau, tt elmu
                   'ctr_mu_wjet':    (data_procs, ['NOMINAL']),
                   'ctr_mu_wjet_old':(data_procs, ['NOMINAL']),
                   'ctr_mu_dy_mumu': (data_procs, ['NOMINAL']),
                   'ctr_mu_dy_tt':   (data_procs, ['NOMINAL']),
                   'ctr_mu_dy_tt_ss':(data_procs, ['NOMINAL']),
                   'ctr_mu_dy_SV_tt':   (data_procs, systematic_names_all),
                   'ctr_mu_dy_SV_tt_ss':(data_procs, systematic_names_all),
                   'ctr_mu_tt_em':   (data_procs, systematic_names_toppt),
                    }
        usual_process = 'data'

    proc = usual_process

    # test
    logger.write('%s\n' % '\n'.join('%s_%s_%s' % (chan, proc, sys[0]) for chan, ((procs, _), sys) in channels.items() for proc in procs))

    # channel -- reco selection
    # proc    -- MC gen info, like inclusive tt VS tt->mutau and others,
    #            choice of subprocesses depends on channel (sadly),
    #            (find most precise subprocess ones, then store accordingly for channels)
    # sys     -- shape systematics
    out_hs = OrderedDict([((chan, proc, sys), {'met':        TH1D('%s_%s_%s_met'        % (chan, proc, sys), '', 30, 0, 300),
                                               'init_met':   TH1D('%s_%s_%s_init_met'   % (chan, proc, sys), '', 30, 0, 300),
                                               'lep_pt':     TH1D('%s_%s_%s_lep_pt'     % (chan, proc, sys), '', 20, 0, 200),
                                               'lep_eta':    TH1D('%s_%s_%s_lep_eta'    % (chan, proc, sys), '', 20, -2.5, 2.5),
                                               'tau_pt':     TH1D('%s_%s_%s_tau_pt'     % (chan, proc, sys), '', 30, 0, 150),
                                               'tau_eta':    TH1D('%s_%s_%s_tau_eta'    % (chan, proc, sys), '', 20, -2.5, 2.5),
                                               'jet_pt':     TH1D('%s_%s_%s_jet_pt'     % (chan, proc, sys), '', 30, 0, 300),
                                               'jet_eta':    TH1D('%s_%s_%s_jet_eta'    % (chan, proc, sys), '', 20, -2.5, 2.5),
                                               'bjet_pt':    TH1D('%s_%s_%s_bjet_pt'    % (chan, proc, sys), '', 30, 0, 300),
                                               'bjet_eta':   TH1D('%s_%s_%s_bjet_eta'   % (chan, proc, sys), '', 20, -2.5, 2.5),
                                               'b_discr_all_jets':  TH1D('%s_%s_%s_b_discr_all_jets' % (chan, proc, sys), '', 30, 0., 1.),
                                               'b_discr_lead_jet':  TH1D('%s_%s_%s_b_discr_lead_jet' % (chan, proc, sys), '', 30, 0., 1.),
                                               #'bdiscr_max':       TH1D('%s_%s_%s_b_discr_max'      % (chan, proc, sys), '', 30, 0, 300),
                                               'dphi_lep_met': TH1D('%s_%s_%s_dphi_lep_met' % (chan, proc, sys), '', 20, -3.2, 3.2),
                                               'cos_dphi_lep_met': TH1D('%s_%s_%s_cos_dphi_lep_met' % (chan, proc, sys), '', 20, -1.1, 1.1),
                                               #'Mt_lep_met_f_mth':   TH1D('%s_%s_%s_Mt_lep_met_f_mth'   % (chan, proc, sys), '', 20, 0, 250),
                                               #'Mt_lep_met_f_cos':   TH1D('%s_%s_%s_Mt_lep_met_f_cos'   % (chan, proc, sys), '', 20, 0, 250),
                                               #'Mt_lep_met_f_cos_c': TH1D('%s_%s_%s_Mt_lep_met_f_cos_c' % (chan, proc, sys), '', 20, 0, 250),
                                               #'Mt_lep_met_f_c':     TH1D('%s_%s_%s_Mt_lep_met_f_c'     % (chan, proc, sys), '', 20, 0, 250),
                                               #'Mt_lep_met_f_test':  TH1D('%s_%s_%s_Mt_lep_met_f_test'  % (chan, proc, sys), '', 20, 0, 250),
                                               'Mt_lep_met_f':       TH1D('%s_%s_%s_Mt_lep_met_f'       % (chan, proc, sys), '', 20, 0, 250),
                                               # control distrs for effect of different weights
                                               'Mt_lep_met_f_init':       TH1D('%s_%s_%s_Mt_lep_met_f_init'       % (chan, proc, sys), '', 20, 0, 250),
                                               #sys_weight = weight * weight_bSF * weight_PU * weight_top_pt
                                               'Mt_lep_met_f_w_in':       TH1D('%s_%s_%s_Mt_lep_met_f_w_in'       % (chan, proc, sys), '', 20, 0, 250),
                                               'Mt_lep_met_f_w_pu':       TH1D('%s_%s_%s_Mt_lep_met_f_w_pu'       % (chan, proc, sys), '', 20, 0, 250),
                                               'Mt_lep_met_f_w_bf':       TH1D('%s_%s_%s_Mt_lep_met_f_w_bf'       % (chan, proc, sys), '', 20, 0, 250),
                                               'Mt_lep_met_f_w_bf_min':   TH1D('%s_%s_%s_Mt_lep_met_f_w_bf_min'   % (chan, proc, sys), '', 20, 0, 250),
                                               'Mt_lep_met_f_w_tp':       TH1D('%s_%s_%s_Mt_lep_met_f_w_tp'       % (chan, proc, sys), '', 20, 0, 250),
                                               'Mt_lep_met':  TH1D('%s_%s_%s_Mt_lep_met' % (chan, proc, sys), '', 10, 0, 200),
                                               'Mt_tau_met':  TH1D('%s_%s_%s_Mt_tau_met' % (chan, proc, sys), '', 20, 0, 200),
                                               # for dileptons, it is practically the same as lep+tau, but for simplicity keeping them separate
                                               'M_lep_lep':   TH1D('%s_%s_%s_M_lep_lep'  % (chan, proc, sys), '', 20, 0, 150),
                                               'M_lep_tau':   TH1D('%s_%s_%s_M_lep_tau'  % (chan, proc, sys), '', 20, 0, 200),
                                               'tau_SV_sign': TH1D('%s_%s_%s_tau_SV_sign'% (chan, proc, sys), '', 21, -1, 20),
                                               'tau_SV_leng': TH1D('%s_%s_%s_tau_SV_leng'% (chan, proc, sys), '', 21, -0.1, 1.),
                                               'tau_jet_bdiscr': TH1D('%s_%s_%s_tau_jet_bdiscr'  % (chan, proc, sys), '', 20, -0.1, 1.1),
                                               'tau_sign_bdiscr':TH2D('%s_%s_%s_tau_sign_bdiscr' % (chan, proc, sys), '', 21, -1, 20, 20, 0., 1.),
                                               'tau_leng_bdiscr':TH2D('%s_%s_%s_tau_leng_bdiscr' % (chan, proc, sys), '', 21, -0.1, 1., 20, 0., 1.),
                                               'tau_sign_energy':TH2D('%s_%s_%s_tau_sign_energy' % (chan, proc, sys), '', 21, -1,  20., 20, 20., 150.),
                                               'tau_leng_energy':TH2D('%s_%s_%s_tau_leng_energy' % (chan, proc, sys), '', 21, -0.1, 1., 20, 20., 150.),
                                               'nvtx':        TH1D('%s_%s_%s_nvtx'       % (chan, proc, sys), '', 50, 0, 50),
                                               'nvtx_gen':    TH1D('%s_%s_%s_nvtx_gen'   % (chan, proc, sys), '', 100, 0, 100),
                                               # TODO: add rho to ntuples
                                               'rho':         TH1D('%s_%s_%s_rho'        % (chan, proc, sys), '', 50, 0, 50),
                                               'njets':       TH1D('%s_%s_%s_njets'      % (chan, proc, sys), '', 10, 0, 10),
                                               'nbjets':      TH1D('%s_%s_%s_nbjets'     % (chan, proc, sys), '', 5, 0, 5),
                                               'dijet_mass':  TH1D('%s_%s_%s_dijet_mass'  % (chan, proc, sys), '', 20, 0, 200),
                                               'trijet_mass': TH1D('%s_%s_%s_trijet_mass' % (chan, proc, sys), '', 20, 0, 400),
                                               '2D_dijet_trijet':   TH2D('%s_%s_%s_2D_dijet_trijet'   % (chan, proc, sys), '', 20, 0, 200, 20, 0, 300),
                                               'dijet_trijet_mass': TH1D('%s_%s_%s_dijet_trijet_mass' % (chan, proc, sys), '', 20, 0, 400) })
                for chan, ((procs, _), systs) in channels.items() for proc in procs for sys in systs])

    # strange, getting PyROOT_NoneObjects from these after output
    for _, histos in out_hs.items():
        for h in histos.values():
            h.SetDirectory(0)

    '''
    for d, histos in out_hs.items():
        for name, histo in histos.items():
            print(d, name)
            histo.Print()
    '''

    logger.write("N entries: %d\n" % tree.GetEntries())
    if not range_max or range_max > tree.GetEntries():
        range_mas = tree.GetEntries()

    profile = cProfile.Profile()
    profile.enable()
    #for iev in range(range_min, range_max): # root sucks
    for iev, ev in enumerate(tree):
        '''
        HLT_el && abs(leps_ID) == 11 && abs(lep0_p4.eta()) < 2.4 && lep0_dxy <0.01 && lep0_dz<0.02 && njets > 2 && met_corrected.pt() > 40 && nbjets > 0
        '''

        '''
        scheme:

        dtag -> process [subprocesses]
        for event:
            for each SYS -> jet_pts, tau_pts, weight
                which (reco) [channels] it passes
                find (gen) subprocess <------! different subprocess def-s for diff channels
                record distr-s for each
        '''

        if iev <  range_min: continue
        if iev >= range_max: break
        #ev = tree.GetEntry(iev)

        # the lepton requirements for all 1-lepton channels:
        # TODO: is there relIso cut now?

        pass_mu = abs(ev.leps_ID) == 13 and ev.HLT_mu and ev.lep_matched_HLT[0] and ev.lep_p4[0].pt() > 27 and abs(ev.lep_p4[0].eta()) < 2.4 and ev.lep_dxy[0] < 0.01 and ev.lep_dz[0] < 0.02 # lep_relIso[0]
        pass_el = abs(ev.leps_ID) == 11 and ev.HLT_el and ev.lep_matched_HLT[0] and ev.lep_p4[0].pt() > 30 and abs(ev.lep_p4[0].eta()) < 2.4 and ev.lep_dxy[0] < 0.01 and ev.lep_dz[0] < 0.02 # lep_relIso[0]
        pass_elmu = ev.leps_ID == -11*13 and ev.HLT_mu and \
            (ev.lep_matched_HLT[0] if abs(ev.lep_id[0]) == 13 else ev.lep_matched_HLT[1]) and \
            (ev.lep_p4[0].pt() > 30 and abs(ev.lep_p4[0].eta()) < 2.4 and ev.lep_dxy[0] < 0.01 and ev.lep_dz[0] < 0.02) and \
            (ev.lep_p4[1].pt() > 30 and abs(ev.lep_p4[1].eta()) < 2.4 and ev.lep_dxy[1] < 0.01 and ev.lep_dz[1] < 0.02)
        pass_mumu = ev.leps_ID == -13*13 and (ev.HLT_mu) and (ev.lep_matched_HLT[0] or ev.lep_matched_HLT[1]) and \
            (ev.lep_p4[0].pt() > 30 and abs(ev.lep_p4[0].eta()) < 2.4 and ev.lep_dxy[0] < 0.01 and ev.lep_dz[0] < 0.02) and \
            (ev.lep_p4[1].pt() > 30 and abs(ev.lep_p4[1].eta()) < 2.4 and ev.lep_dxy[1] < 0.01 and ev.lep_dz[1] < 0.02)

        # minimum possible pt threshold -- 24 GeV, = to HLT and recorded in Ntupler
        pass_min_mu = abs(ev.leps_ID) == 13 and ev.HLT_mu and ev.lep_matched_HLT[0] and ev.lep_p4[0].pt() > 24. and abs(ev.lep_p4[0].eta()) < 2.4 and ev.lep_dxy[0] < 0.01 and ev.lep_dz[0] < 0.02 # lep_relIso[0]

        # 1-lep channels, 2mu DY and el-mu ttbar
        if not (pass_min_mu or pass_mu or pass_el or pass_mumu or pass_elmu): continue

        # also at least some kind of tau in single-el:
        if (pass_mu or pass_el) and (not ev.tau_p4.size() > 0): continue # this is the only thing reduces computing

        # expensive calls and they don't depend on systematics now
        if doRecoilCorrections:
            #def transverse_mass_pts(v1_x, v1_y, v2_x, v2_y):
            #met_x = ev.pfmetcorr_ex
            #met_y = ev.pfmetcorr_ey
            # no, recalculated them
            met_x = ROOT.met_pt_recoilcor_x(
                ev.met_corrected.Px(), # uncorrected type I pf met px (float)
                ev.met_corrected.Py(), # uncorrected type I pf met py (float)
                ev.gen_genPx, # generator Z/W/Higgs px (float)
                ev.gen_genPy, # generator Z/W/Higgs py (float)
                ev.gen_visPx, # generator visible Z/W/Higgs px (float)
                ev.gen_visPy, # generator visible Z/W/Higgs py (float)
                0  # max recoil (I checked that -- plot in scrap/recoil-corrections-study)
                #ev.nalljets  # number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
                )
            met_y = ROOT.met_pt_recoilcor_y(
                ev.met_corrected.Px(), # uncorrected type I pf met px (float)
                ev.met_corrected.Py(), # uncorrected type I pf met py (float)
                ev.gen_genPx, # generator Z/W/Higgs px (float)
                ev.gen_genPy, # generator Z/W/Higgs py (float)
                ev.gen_visPx, # generator visible Z/W/Higgs px (float)
                ev.gen_visPy, # generator visible Z/W/Higgs py (float)
                0
                #ev.nalljets  # number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
                )
            Mt_lep_met_c   = ROOT.MTlep_met_pt_recoilcor(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), ev.met_corrected.Px(), ev.met_corrected.Py(), ev.gen_genPx, ev.gen_genPy, ev.gen_visPx, ev.gen_visPy, 0)

            #Mt_lep_met = transverse_mass_pts(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), ev.pfmetcorr_ex, ev.pfmetcorr_ey)
            #Mt_tau_met_nominal = transverse_mass_pts(ev.tau_p4[0].Px(), ev.tau_p4[0].Py(), ev.pfmetcorr_ex, ev.pfmetcorr_ey)
            # TODO: tau is corrected with systematic ES
        else:
            met_x = ev.met_corrected.Px()
            met_y = ev.met_corrected.Py()
            Mt_lep_met_c   = ROOT.MTlep_c(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), ev.met_corrected.Px(), ev.met_corrected.Py())
            #Mt_lep_met_test = transverse_mass_pts(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), met_x, met_y)
            #Mt_lep_met = transverse_mass(ev.lep_p4[0], ev.met_corrected)
            #Mt_tau_met = transverse_mass(ev.tau_p4[0], ev.met_corrected)
            #Mt_tau_met_nominal = transverse_mass_pts(ev.tau_p4[0].Px(), ev.tau_p4[0].Py(), ev.met_corrected.Px(), ev.met_corrected.Py())
        # also
        #Mt_lep_met_d = (ev.lep_p4[0] + ev.met_corrected).Mt()


        # TODO: pass jet systematics to met?
        Mt_lep_met = transverse_mass_pts(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), met_x, met_y)
        #met_vector = TLorentzVector(met_x, met_y, 0, TMath.Sqrt(met_x*met_x + met_y*met_y))
        #met_vector = ROOT.Math.LorentzVector(met_x, met_y, 0, TMath.Sqrt(met_x*met_x + met_y*met_y))
        #met_vector = ROOT.Math.LorentzVector(ROOT.Math.PxPyPzE4D(float))(met_x, met_y, 0, TMath.Sqrt(met_x*met_x + met_y*met_y))
        met_vector = ROOT.Math.LorentzVector(ROOT.Math.PxPyPzE4D('double'))(met_x, met_y, 0, TMath.Sqrt(met_x*met_x + met_y*met_y))
        d_lep_met = ev.lep_p4[0] - met_vector
        dphi_lep_met = d_lep_met.Phi() # ev.lep_p4[0].Phi() - met_vector.Phi()
        cos_dphi_lep_met = TMath.Cos(dphi_lep_met)
        met_pt = TMath.Sqrt(met_x*met_x + met_y*met_y)

        #float MTlep_met_pt_recoilcor(float lep_px, float lep_py,
        #        float met_px, float met_py,
        #        float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
        #        int njets)
        Mt_lep_met_cos = TMath.Sqrt(2*ev.lep_p4[0].pt() * met_vector.pt() * (1 - cos_dphi_lep_met))
        Mt_lep_met_mth = abs((ev.lep_p4[0] - met_vector).Mt())
        Mt_lep_met_cos_c = ROOT.transverse_mass(ev.lep_p4[0], met_vector)

        weight = 1. # common weight of event (1. for data)
        if isMC:
            try:
                weight_pu    = pileup_ratio[ev.nvtx_gen]
                weight_pu_up = pileup_ratio_up[ev.nvtx_gen]
                weight_pu_dn = pileup_ratio_down[ev.nvtx_gen]
                control_hs['weight_pu']   .Fill(weight_pu)
                control_hs['weight_pu_up'].Fill(weight_pu_up)
                control_hs['weight_pu_dn'].Fill(weight_pu_dn)
            except:
                #print i, ev.nvtx
                continue

            if aMCatNLO and ev.aMCatNLO_weight < 0:
                weight *= -1

            weight_z_mass_pt = 1.
            if isDY:
                # float zPtMass_weight(float genMass, float genPt)
                weight_z_mass_pt *= zPtMass_weight(ev.genMass, ev.genPt)
                weight *= weight_z_mass_pt
                # DY has this trick, which I had to solve cleanly in ntuples
                # but I just dumped all info downstream..
                # it might have actual Z particle in decay chain
                # or not (and no gamma too) -- then you judge from "prompt" leptons
                # the mutually exclusive case never appears
                if ev.gen_N_zdecays > 0:
                    lep1_id = abs(ev.gen_zdecays_IDs[0])
                    lep2_id = abs(ev.gen_zdecays_IDs[1])
                else:
                    # check prompt leptns
                    # if no Z decay the leptons are there
                    lep1_id = abs(ev.gen_pythia8_prompt_leptons_IDs[0])
                    lep2_id = abs(ev.gen_pythia8_prompt_leptons_IDs[1])
                # TODO: actually track tau decays fro DY? -- no need, it's a small background
                if lep1_id >= 15 and lep2_id >= 15:
                    proc = 'dy_tautau'
                else:
                    proc = 'dy_other'

            weight_top_pt = 1.
            # "Only top quarks in SM ttbar events must be reweighted, not single tops or tops from BSM production mechanisms."
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
            if isTT:
                weight_top_pt = ttbar_pT_SF(ev.gen_t_pt, ev.gen_tb_pt)
                #weight *= weight_top_pt # to sys
                control_hs['weight_top_pt']   .Fill(weight_top_pt)
                # basically only difference is eltau/mutau
                t_wid  = abs(ev.gen_t_w_decay_id)
                tb_wid = abs(ev.gen_tb_w_decay_id)
                if (t_wid > 15*15 and tb_wid == 13) or (t_wid == 13 and tb_wid > 15*15): # lt
                    proc = 'tt_mutau'
                elif (t_wid > 15*15 and tb_wid == 11) or (t_wid == 11 and tb_wid > 15*15): # lt
                    proc = 'tt_eltau'
                elif t_wid * tb_wid == 13 or t_wid * tb_wid == 11: # lj
                    proc = 'tt_lj'
                elif (t_wid > 15*15 and (tb_wid == 11*15 or tb_wid == 13*15)) or \
                     ((t_wid == 11*15 or t_wid == 13*15) and tb_wid > 15*15): # taul tauh
                    proc = 'tt_taultauh'
                else:
                    proc = 'tt_other'

            if isSTop:
                # basically only difference is eltau/mutau
                w1_id = abs(ev.gen_wdecays_IDs[0])
                w2_id = 1
                if not isSTopTSchannels:
                    w2_id = abs(ev.gen_wdecays_IDs[1])
                # t/s channels emit top and a quark -- top ID + jets
                if (w1_id > 15*15 and w2_id == 13) or (w1_id == 13 and w2_id > 15*15): # lt
                    proc = 's_top_mutau'
                elif (w1_id > 15*15 and w2_id == 11) or (w1_id == 11 and w2_id > 15*15): # lt
                    proc = 's_top_eltau'
                elif w1_id * w2_id == 13 or w1_id * w2_id == 11: # lj
                    proc = 's_top_lj'
                #elif (w1_id > 15*15 and (w2_id == 11*15 or w2_id == 13*15)) or
                #     ((w1_id == 11*15 or w1_id == 13*15) and w2_id > 15*15): # taul tauh
                #    proc = 'tt_taultauh'
                else:
                    proc = 's_top_other'

            if (pass_mu or pass_elmu or pass_mumu) and isMC:
                mu_sfs = lepton_muon_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                mu_trg_sf = lepton_muon_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                # bcdef gh eras
                weight *= ratio_bcdef * mu_trg_sf[0] * mu_sfs[0] * mu_sfs[1] * mu_sfs[2] + ratio_gh * mu_trg_sf[1] * mu_sfs[3] * mu_sfs[4] * mu_sfs[5]

            if pass_el and isMC:
                el_sfs = lepton_electron_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                el_trg_sf = lepton_electron_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                weight *= el_trg_sf * el_sfs[0] * el_sfs[1]

        #JETS

        all_jets_b_discrs = []
        lead_jet_b_discr = -1.

        if len(ev.jet_b_discr) > 0:
            lead_jet_b_discr = ev.jet_b_discr[0]

        b_tag_wp = 0.8484

        # lists of "light jets", not passing b-tag
        jets_nominal = [] # nominal jet pts
        jets_jer_up, jets_jer_down = [], []
        jets_jes_up, jets_jes_down = [], []
        # lists of "heavy jets", passing b-tag
        jets_b_nominal = [] # nominal jet pts
        jets_b_jer_up, jets_b_jer_down = [], []
        jets_b_jes_up, jets_b_jes_down = [], []

        #nbjets_nominal = 0
        #nbjets_jer_up, nbjets_jer_down = 0, 0
        #nbjets_jes_up, nbjets_jes_down = 0, 0

        weight_bSF = 1.
        weight_bSF_up, weight_bSF_down = 1., 1.
        weight_bSF_jer_up, weight_bSF_jer_down = 1., 1.
        weight_bSF_jes_up, weight_bSF_jes_down = 1., 1.

        jets_nominal_min = [] # nominal jet pts
        jets_b_nominal_min = [] # nominal jet pts
        weight_bSF_min = 1.

        for i in xrange(ev.jet_p4.size()):
            pfid, p4 = ev.jet_PFID[i], ev.jet_p4[i]
            if pfid < 1 or abs(p4.eta()) > 2.5: continue # Loose PFID and eta

            jet_b_discr = ev.jet_b_discr[i]
            all_jets_b_discrs.append(jet_b_discr)
            b_tagged = jet_b_discr > b_tag_wp
            flavId = ev.jet_hadronFlavour[i]

            if p4.pt() > 30: # nominal jet
                if b_tagged:
                    #nbjets_nominal += 1
                    jets_b_nominal.append((p4, 1))
                else:
                    jets_nominal.append((p4, 1))
                if isMC and with_bSF: # do the SF weights
                    weight_bSF      *= calc_btag_sf_weight(b_tagged, flavId, p4.pt(), p4.eta())
                    weight_bSF_up   *= calc_btag_sf_weight(b_tagged, flavId, p4.pt(), p4.eta(), "up")
                    weight_bSF_down *= calc_btag_sf_weight(b_tagged, flavId, p4.pt(), p4.eta(), "down")

            # also save "minimal" jets -- with minimal pt threshold
            # make them only at NOMINAL systematics for now
            if p4.pt() > 20: # nominal jet
                if b_tagged:
                    #nbjets_nominal += 1
                    jets_b_nominal_min.append((p4, 1))
                else:
                    jets_nominal_min.append((p4, 1))
                if isMC and with_bSF: # do the SF weights
                    weight_bSF_min  *= calc_btag_sf_weight(b_tagged, flavId, p4.pt(), p4.eta())
                    #weight_bSF_up   *= calc_btag_sf_weight(b_tagged, flavId, p4.pt(), p4.eta(), "up")
                    #weight_bSF_down *= calc_btag_sf_weight(b_tagged, flavId, p4.pt(), p4.eta(), "down")

            # TODO: probably these weights are better to be computed at channel recording?

            if isMC:
                # vary the jets with systematics, save b-tagged and other
                factor, up, down = ev.jet_jer_factor[i], ev.jet_jer_factor_up[i], ev.jet_jer_factor_down[i]
                #jes_shift = ev.jet_jes_correction_relShift[i]
                jes_shift = ev.jet_jes_uncertainty[i]
                # 
                mult = (up/factor) if factor > 0 else 0
                if p4.pt() * mult > 30: # jer up
                    if b_tagged:
                        #nbjets_jer_up += 1
                        jets_b_jer_up.append((p4, mult))
                    else:
                        jets_jer_up.append((p4, mult))
                    if isMC and with_bSF:
                        weight_bSF_jer_up *= calc_btag_sf_weight(b_tagged, flavId, p4.pt() * mult, p4.eta())

                mult = (down/factor) if factor > 0 else 0
                if p4.pt() * mult > 30: # jer down
                    if b_tagged:
                        #nbjets_jer_down += 1
                        jets_b_jer_down.append((p4, mult))
                    else:
                        jets_jer_down.append((p4, mult))
                    if isMC and with_bSF:
                        weight_bSF_jer_down *= calc_btag_sf_weight(b_tagged, flavId, p4.pt() * mult, p4.eta())

                mult = 1 + jes_shift
                if p4.pt() * mult > 30: # jer up
                    if b_tagged:
                        #nbjets_jes_up += 1
                        jets_b_jes_up.append((p4, mult))
                    else:
                        jets_jes_up.append((p4, mult))
                    if isMC and with_bSF:
                        weight_bSF_jes_up *= calc_btag_sf_weight(b_tagged, flavId, p4.pt() * mult, p4.eta())

                mult = 1 - jes_shift
                if p4.pt() * mult > 30: # jer down
                    if b_tagged:
                        #nbjets_jes_down += 1
                        jets_b_jes_down.append((p4, mult))
                    else:
                        jets_jes_down.append((p4, mult))
                    if isMC and with_bSF:
                        weight_bSF_jes_down *= calc_btag_sf_weight(b_tagged, flavId, p4.pt() * mult, p4.eta())

            # save jets and factors for lj_var
            # count b jets
            # how to calc b SF weight? <---- just calc them, only for nominal jets
            #

        # tau pt-s
        # ES correction
        # modes have different correction but the same uncertainty = +- 1.2% = 0.012
        # uncertainties are not correlated, but I'll do correlated UP/DOWN -- all modes UP or all modes DOWN
        #tau_pts_corrected = []
        #tau_pts_corrected_up = []
        #tau_pts_corrected_down = []
        #for i, (p4, DM, IDlev) in enumerate(zip(ev.tau_p4, ev.tau_decayMode, ev.tau_IDlev)):
        #Mt_tau_met_nominal, Mt_tau_met_up, Mt_tau_met_down = None, None, None
        Mt_tau_met_nominal, Mt_tau_met_up, Mt_tau_met_down = 0, 0, 0

        # so, actually what I need from taus is
        # whether there is a medium tau with pt 30, eta 2.4
        # and if it is OS with the lepton
        # usually it is the tau on first position (0)
        # should I really loop?

        #for i in range(ev.tau_p4.size()):
        #    # it should work like Python does and not copy these objects! (cast)
        #    p4, DM, IDlev = ev.tau_p4[i], ev.tau_decayMode[i], ev.tau_IDlev[i]
        #    if IDlev < 3 or abs(p4.eta()) < 2.4: continue # only Medium taus
        #    if DM == 0:
        #      factor = 0.995
        #      if isMC:
        #        factor_up   = 0.995 + 0.012
        #        factor_down = 0.995 - 0.012
        #    elif DM < 10:
        #      factor = 1.011
        #      if isMC:
        #        factor_up   = 1.011 + 0.012
        #        factor_down = 1.011 - 0.012
        #    else:
        #      factor = 1.006
        #      if isMC:
        #        factor_up   = 1.006 + 0.012
        #        factor_down = 1.006 - 0.012
        #    tau_pts_corrected.append(p4.pt() * factor)
        #    if not Mt_tau_met_nominal: Mt_tau_met_nominal = transverse_mass_pts(ev.tau_p4[0].Px()*factor, ev.tau_p4[0].Py()*factor, met_x, met_y)
        #    if isMC:
        #        tau_pts_corrected_up.append(p4.pt() * factor_up)
        #        tau_pts_corrected_down.append(p4.pt() * factor_down)
        #        if not Mt_tau_met_up:
        #            Mt_tau_met_up   = transverse_mass_pts(ev.tau_p4[0].Px()*factor_up, ev.tau_p4[0].Py()*factor_down, met_x, met_y)
        #            Mt_tau_met_down = transverse_mass_pts(ev.tau_p4[0].Px()*factor_up, ev.tau_p4[0].Py()*factor_down, met_x, met_y)

        #has_tau = False
        #has_tau_es_up = False
        #has_tau_es_down = False

        # p4 momenta of taus and the Energy Scale factor
        taus_nominal_min = []
        taus_nominal = []
        taus_es_up   = []
        taus_es_down = []

        # only top pt tau is treated but that's fine
        '''IDlev
        1 VLoose
        2 Loose
        3 Medium
        4 Tight
        5 VTight
        '''
        if ev.tau_p4.size() > 0 and ev.tau_IDlev[0] > 1 and abs(ev.tau_p4[0].eta()) < 2.4:
            # it should work like Python does and not copy these objects! (cast)
            #p4, DM, IDlev = ev.tau_p4[0], ev.tau_decayMode[0], ev.tau_IDlev[0]
            #if IDlev < 3 or abs(p4.eta()) < 2.4: continue # only Medium taus

            # check dR to lepton
            tlep_p4 = TLorentzVector(ev.lep_p4[0].X(), ev.lep_p4[0].Y(), ev.lep_p4[0].Z(), ev.lep_p4[0].T())
            ttau_p4 = TLorentzVector(ev.tau_p4[0].X(), ev.tau_p4[0].Y(), ev.tau_p4[0].Z(), ev.tau_p4[0].T())
            if tlep_p4.DeltaR(ttau_p4) > 0.3:
                if not isMC:
                    if ev.tau_p4[0].pt() > 20.:
                        taus_nominal_min.append((ev.tau_p4[0], 1.))
                    if ev.tau_p4[0].pt() > 30. and ev.tau_IDlev[0] > 2:
                        taus_nominal.append((ev.tau_p4[0], 1.))
                else:
                    DM = ev.tau_decayMode[0]
                    if DM == 0:
                        factor = 0.995
                        factor_up   = 0.995 + 0.012
                        factor_down = 0.995 - 0.012
                    elif DM < 10:
                        factor = 1.011
                        factor_up   = 1.011 + 0.012
                        factor_down = 1.011 - 0.012
                    else:
                        factor = 1.006
                        factor_up   = 1.006 + 0.012
                        factor_down = 1.006 - 0.012

                    # only nominal for min taus
                    if ev.tau_p4[0].pt() * factor > 20.:
                        taus_nominal_min.append((ev.tau_p4[0], factor))
                    ## calculate it later, inplace of record
                    #if not Mt_tau_met_nominal:
                    #    Mt_tau_met_nominal = transverse_mass_pts(ev.tau_p4[0].Px()*factor, ev.tau_p4[0].Py()*factor, met_x, met_y)
                    #has_tau_es_up   = (ev.tau_p4[0].pt() * factor_up  ) > 30
                    #has_tau_es_down = (ev.tau_p4[0].pt() * factor_down) > 30
                    if ev.tau_IDlev[0] > 2:
                        if ev.tau_p4[0].pt() * factor > 30.:
                            taus_nominal.append((ev.tau_p4[0], factor))
                        if ev.tau_p4[0].pt() * factor_up > 30.:
                            taus_es_up.append((ev.tau_p4[0], factor_up))
                        if ev.tau_p4[0].pt() * factor_down > 30.:
                            taus_es_down.append((ev.tau_p4[0], factor_down))
                    #if not Mt_tau_met_up:
                    #    Mt_tau_met_up   = transverse_mass_pts(ev.tau_p4[0].Px()*factor_up, ev.tau_p4[0].Py()*factor_up, met_x, met_y)
                    #    Mt_tau_met_down = transverse_mass_pts(ev.tau_p4[0].Px()*factor_down, ev.tau_p4[0].Py()*factor_down, met_x, met_y)

        #has_medium_tau = any(IDlev > 2 and p4.pt() > 30 for IDlev, p4 in zip(ev.tau_IDlev, ev.tau_p4))
        #has_medium_tau = ev.tau_IDlev.size() > 0 and ev.tau_IDlev[0] > 2 and ev.tau_p4[0].pt() > 30
        #has_medium_tau = bool(tau_pts_corrected)
        #TODO: propagate TES to MET?


        # shape systematics are:
        # corrected jet pts and tau pts
        # weight variations
        #nominal systematics
        # jet pts, tau pts, b weight (=1 for data), pu weight (=1 for data)

        if isMC:
            systematics = {'NOMINAL'   : [jets_nominal,  jets_b_nominal,  taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, weight_bSF,          weight_pu,    1],
                           'JESUp'     : [jets_jes_up,   jets_b_jes_up,   taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, weight_bSF_jes_up,   weight_pu,    1],
                           'JESDown'   : [jets_jes_down, jets_b_jes_down, taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, weight_bSF_jes_down, weight_pu,    1],
                           'JERUp'     : [jets_jer_up,   jets_b_jer_up,   taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, weight_bSF_jer_up,   weight_pu,    1],
                           'JERDown'   : [jets_jer_down, jets_b_jer_down, taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, weight_bSF_jer_down, weight_pu,    1],
                           'TauESUp'   : [jets_nominal,  jets_b_nominal,  taus_es_up  , jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min,  weight_bSF,          weight_pu,    1],
                           'TauESDown' : [jets_nominal,  jets_b_nominal,  taus_es_down, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min,  weight_bSF,          weight_pu,    1],
                           'bSFUp'     : [jets_nominal,  jets_b_nominal,  taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, weight_bSF_up  ,     weight_pu,    1],
                           'bSFDown'   : [jets_nominal,  jets_b_nominal,  taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, weight_bSF_down,     weight_pu,    1],
                           'PUUp'      : [jets_nominal,  jets_b_nominal,  taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, weight_bSF,          weight_pu_up, 1],
                           'PUDown'    : [jets_nominal,  jets_b_nominal,  taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, weight_bSF,          weight_pu_dn, 1],
                           }
            if isTT:
                systematics['TOPPTUp']   = [jets_nominal, jets_b_nominal, taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, weight_bSF, weight_pu, weight_top_pt]
                systematics['TOPPTDown'] = [jets_nominal, jets_b_nominal, taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, weight_bSF, weight_pu, 1.]
        else:
            systematics = {'NOMINAL': [jets_nominal, jets_b_nominal, taus_nominal, jets_nominal_min, jets_b_nominal_min, taus_nominal_min, weight_bSF_min, 1., 1., 1.]}

        # for each systematic
        # pass one of the reco selections
        # check the subprocess
        # store distr

        for sys_name, (jets, jets_b, taus, jets_min, jets_b_min, taus_min, weight_bSF_min, weight_bSF, weight_PU, weight_top_pt) in systematics.items():
            sys_weight = weight * weight_bSF * weight_PU * weight_top_pt
            sys_weight_min = weight * weight_bSF_min * weight_PU * weight_top_pt
            # pass reco selections

            Mt_tau_met = 0
            if taus:# there are taus
                # Mt_tau_met_nominal = transverse_mass_pts(ev.tau_p4[0].Px()*factor, ev.tau_p4[0].Py()*factor, met_x, met_y)
                tau_p4, factor = taus[0]
                Mt_tau_met = transverse_mass_pts(tau_p4.Px()*factor, tau_p4.Py()*factor, met_x, met_y)

            # from channel_distrs
            # --- 'mu': {'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0]  && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4',
            # met_corrected.pt() > 40 && nbjets > 0 && tau_IDlev[0] > 2 && tau_id[0]*lep_id[0] < 0

            # --- 'mujets':
            #{'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4
            #&& tau_IDlev[0] > 0 && transverse_mass(lep_p4[0].pt(), met_corrected.pt(), lep_p4[0].phi(), met_corrected.phi()) > 50',
            # --- 'mujets_b':
            #{'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4
            # nbjets > 0 
            # transverse_mass(lep_p4[0].pt(), met_corrected.pt(), lep_p4[0].phi(), met_corrected.phi()) > 50',

            # --- 'taumutauh':
            #{'common': 'HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] 
            # nbjets == 0 
            # tau_IDlev[0] > 2 && tau_p4[0].pt() > 30 
            # lep_p4[0].pt() > 30 && fabs(lep_p4[0].eta()) < 2.4 
            # sqrt(2*lep_p4[0].pt()*met_corrected.pt() * (1 - cos(lep_p4[0].phi() - met_corrected.phi()))) < 40 
            # divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) > 45 
            # divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) < 85',
            # --- taumutauh_antiMt <-- anti Mt and ant DY dilep mass
            # --- taumutauh_antiMt_pretau
            # --- taumutauh_antiMt_pretau_allb

            # piecemeal
            # njets > 2 && met_corrected.pt() > 40 && nbjets > 0
            #met_pt = ev.met_corrected.pt()
            large_met = met_pt > 40
            #large_Mt_lep = Mt_lep_met > 40
            nbjets = len(jets_b)

            tau_matched_dR_bjet = False
            if taus: # there are taus, check if b-jets overlap with it
                tau_p4 = taus[0][0]
                for j_p4, _ in jets_b:
                    # these are LorentzVectors -- they don't have DeltaR,
                    # TLorentzVectors do, but they don't have "mass" method...
                    tj_p4   = TLorentzVector(j_p4.X(), j_p4.Y(), j_p4.Z(), j_p4.T())
                    ttau_p4 = TLorentzVector(tau_p4.X(), tau_p4.Y(), tau_p4.Z(), tau_p4.T())
                    if tj_p4.DeltaR(ttau_p4) < 0.3:
                        tau_matched_dR_bjet = True
                        break

            njets = nbjets + len(jets)
            nbjets_no_tau_match = nbjets - tau_matched_dR_bjet
            njets_no_tau_match  = nbjets + njets - tau_matched_dR_bjet

            has_bjets = nbjets > 0
            has_bjets_no_tau_match = nbjets_no_tau_match > 0
            has_3jets = njets > 2
            has_3jets_no_tau_match = njets_no_tau_match > 2

            # the jets with minimum pt cut
            nbjets_min = len(jets_b_min)
            njets_min = nbjets_min + len(jets_min)
            has_3jets_min = njets_min > 2
            has_bjets_min = nbjets_min > 0

            #has_pre_tau = len(tau_pts) > 0
            #has_medium_tau = has_pre_tau and tau_pts[0] > 30 and ev.tau_IDlev[0] > 2
            has_medium_tau = len(taus) > 0
            has_loose_tau = len(taus_min) > 0

            # these are for single-lepton
            os_lep_med_tau = has_medium_tau and ev.tau_id[0]*ev.lep_id[0] < 0
            ss_lep_med_tau = has_medium_tau and ev.tau_id[0]*ev.lep_id[0] > 0
            # dilep_mass is not TES corrected
            # it's used only in control regions and shouldn't be a big deal
            #dy_dilep_mass = has_pre_tau and (45 < (ev.lep_p4[0] + ev.tau_p4[0]).mass() < 85)

            os_lep_loose_tau = has_loose_tau and ev.tau_id[0]*ev.lep_id[0] < 0
            ss_lep_loose_tau = has_loose_tau and ev.tau_id[0]*ev.lep_id[0] > 0

            if has_medium_tau:
            #if len(ev.tau_p4) > 0:
                lep_tau = ev.lep_p4[0] + taus[0][0] * taus[0][1]
                lep_tau_mass = lep_tau.mass()

            pass_single_lep_presel = large_met and has_3jets and has_bjets and (pass_el or pass_mu) #and os_lep_med_tau
            pass_single_lep_sel = pass_single_lep_presel and os_lep_med_tau


            # all the channel selections follow

            passed_channels = []
            if pass_el and pass_single_lep_presel:
                passed_channels.append('el_presel')
            if pass_mu and pass_single_lep_presel:
                passed_channels.append('mu_presel')

            if pass_el and pass_single_lep_sel:
                passed_channels.append('el_sel')
            if pass_mu and pass_single_lep_sel:
                passed_channels.append('mu_sel')

            lj_var, w_mass, t_mass = -11., -11., -11.
            lj_cut = 50
            if pass_single_lep_sel:
                # calc lj_var
                lj_var, w_mass, t_mass = calc_lj_var(jets, jets_b)
                large_lj = lj_var > lj_cut

                # all these channels are exclusive
                fit_category_selection = ('el' if pass_el else 'mu') + '_' + ('lj_out' if large_lj else 'lj')
                passed_channels.append(fit_category_selection)

            # steps for DY and WJets shape control
            if pass_mu and os_lep_med_tau and has_3jets and has_bjets:
                passed_channels.append('mu_sel_nomet')
            if pass_mu and os_lep_med_tau and large_met:
                passed_channels.append('mu_sel_onlymet')
            if pass_mu and os_lep_med_tau and large_met and has_3jets:
                passed_channels.append('mu_sel_nobtag')
            if pass_mu and os_lep_med_tau and has_bjets:
                passed_channels.append('mu_sel_onlybtag')
            if pass_mu and os_lep_med_tau and has_bjets and large_met:
                passed_channels.append('mu_sel_no3jets')
            if pass_mu and os_lep_med_tau and has_3jets:
                passed_channels.append('mu_sel_only3jets')

            # same sign for some QCD control
            if pass_el and pass_single_lep_presel and ss_lep_med_tau:
                passed_channels.append('el_sel_ss')
            if pass_mu and pass_single_lep_presel and ss_lep_med_tau:
                passed_channels.append('mu_sel_ss')
            if pass_el and pass_single_lep_presel and ss_lep_med_tau:
                passed_channels.append('el_presel_ss')
            if pass_mu and pass_single_lep_presel and ss_lep_med_tau:
                passed_channels.append('mu_presel_ss')

            # with tau POG selection
            # mu cand, tau cand, di-mu veto, el veto, mu veto, Mt < 40 and p-zeta > 40 -- DY -> tautau
            # in tt no Mt and p-zeta -- instead at least 1 b
            # and then -- pass or not the medium ID
            # TODO: I don't have di-muon veto + need to add the veto leptons to output somehow -- their counters for example
            if pass_mu and has_bjets:
                passed_channels.append('pog_mu_presel')
            #if pass_mu and has_bjets and has_medium_tau:
            if pass_mu and has_bjets and os_lep_med_tau:
                passed_channels.append('pog_mu_pass')
            if pass_mu and has_bjets and ss_lep_med_tau:
                passed_channels.append('pog_mu_pass_ss')
            if pass_mu and has_bjets and not has_medium_tau:
                passed_channels.append('pog_mu_fail')

            # with addition of no DY mass, no tau match to b-tag (could add a cut on small MT)
            if pass_el and large_met and has_3jets and has_bjets_no_tau_match and os_lep_med_tau: # and (lep_tau_mass < 45 or lep_tau_mass > 85):
                passed_channels.append('adv_el_sel')
            if pass_mu and large_met and has_3jets and has_bjets_no_tau_match and os_lep_med_tau: # and (lep_tau_mass < 45 or lep_tau_mass > 85):
                passed_channels.append('adv_mu_sel')
            # check the mass cut in distribution of lep+tau mass in usual selection

            # the max efficiency requirements of our selection
            if pass_min_mu and large_met and has_3jets_min and has_bjets_min and os_lep_loose_tau: # and (lep_tau_mass < 45 or lep_tau_mass > 85):
                passed_channels.append('sel_mu_min')
            if pass_min_mu and large_met and has_3jets_min and has_bjets_min and ss_lep_loose_tau: # and (lep_tau_mass < 45 or lep_tau_mass > 85):
                passed_channels.append('sel_mu_min_ss')
            # the minimum selection with Medium taus -- hopefully it reduces the QCD
            if pass_min_mu and large_met and has_3jets_min and has_bjets_min and os_lep_med_tau: # and (lep_tau_mass < 45 or lep_tau_mass > 85):
                passed_channels.append('sel_mu_min_medtau')

            # control selections: WJets, DY mumu and tautau, tt elmu
            # WJets control
            #HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && abs(lep_p4[0].eta()) < 2.4 && lep_dxy[0] <0.01 && lep_dz[0]<0.02
            # && tau_IDlev[0] > -1
            # && nbjets == 0
            # && (divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) < 45 || divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) > 85)
            # also may add MT > 50 or met > X -- cuts QCD

            if pass_mu and nbjets == 0 and (Mt_lep_met > 50 or met_pt > 40): # skipped lep+tau mass -- hopefuly DY will be small (-- it became larger than tt)
                passed_channels.append('ctr_mu_wjet')
            # and some tau candidate
            # skipped lep+tau mass -- hopefuly DY will be small (-- it became larger than tt)
            # OS/SS for QCD control
            if pass_mu and nbjets == 0 and len(ev.tau_p4) > 0:
                passed_channels.append('ctr_mu_wjet_old')
            if pass_mu and nbjets == 0 and len(ev.tau_p4) > 0:
                if ev.tau_id[0]*ev.lep_id[0] < 0:
                    passed_channels.append('ctr_mu_wjet_old_os')
                elif ev.tau_id[0]*ev.lep_id[0] > 0:
                    passed_channels.append('ctr_mu_wjet_old_ss')

            # DY -> taumu tauh
            # HLT_mu && abs(leps_ID) == 13 && lep_matched_HLT[0] && abs(lep_p4[0].eta()) < 2.4 && lep_dxy[0] <0.01 && lep_dz[0]<0.02
            # && tau_IDlev[0] > 2 && lep_id[0]*tau_id[0] < 0
            # && divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) > 45 && divector_mass(lep_p4[0].x(), lep_p4[0].y(), lep_p4[0].z(), lep_p4[0].t(), tau_p4[0].x(), tau_p4[0].y(), tau_p4[0].z(), tau_p4[0].t()) < 85
            # && transverse_mass(lep_p4[0].pt(), met_corrected.pt(), lep_p4[0].phi(), met_corrected.phi()) < 40
            tau_SV_cut = 4.
            if pass_mu and os_lep_med_tau and lep_tau_mass > 25 and lep_tau_mass < 105 and Mt_lep_met < 40:
                passed_channels.append('ctr_mu_dy_tt')
                if ev.tau_refited_index[0] > -1. and ev.tau_SV_geom_flightLenSign[ev.tau_refited_index[0]] > tau_SV_cut:
                    passed_channels.append('ctr_mu_dy_SV_tt')
            if pass_mu and ss_lep_med_tau and lep_tau_mass > 25 and lep_tau_mass < 105 and Mt_lep_met < 40:
                passed_channels.append('ctr_mu_dy_tt_ss')
                if ev.tau_refited_index[0] > -1. and ev.tau_SV_geom_flightLenSign[ev.tau_refited_index[0]] > tau_SV_cut:
                    passed_channels.append('ctr_mu_dy_SV_tt')

            if pass_mumu:
                passed_channels.append('ctr_mu_dy_mumu')
            if pass_elmu:
                passed_channels.append('ctr_mu_tt_em')


            # keeping this as legacy?
            if not passed_channels:
                continue

            # this should not be needed now:
            #    proc is set to usual_process,
            #    if the MC is split more -- it is done
            #    the splitting is the same for all channels, it is per-MC not per-MC-in-channel as was supposed at the beginning
            ## now the process variety is restricted to 2 separate channels and TT
            #if isTT and pass_el:
            #    proc = proc if proc in tt_el_procs else 'tt_other'
            #elif isTT and pass_mu:
            #    proc = proc if proc in tt_mu_procs else 'tt_other'

            ## save distrs
            #for chan, proc in chan_subproc_pairs:
            #    out_hs[(chan, proc, sys_name)]['met'].Fill(met_pt, weight)
            #    out_hs[(chan, proc, sys_name)]['Mt_lep_met'].Fill(Mt_lep_met, weight)
            #    out_hs[(chan, proc, sys_name)]['Mt_lep_met_d'].Fill(Mt_lep_met_d, weight)
            #    out_hs[(chan, proc, sys_name)]['dijet_trijet_mass'].Fill(25, weight)
            for chan in passed_channels:
                # some channels might not have only inclusive processes or minimal systematics
                if (chan, proc, sys_name) not in out_hs:
                    continue # TODO: so it doesn't change amount of computing, systematics are per event, not per channel
                    # but it does reduce the amount of output -- no geom progression

                record_weight = sys_weight if chan not in ('sel_mu_min', 'sel_mu_min_ss', 'sel_mu_min_medtau') else sys_weight_min

                out_hs[(chan, proc, sys_name)]['met'].Fill(met_pt, record_weight)
                out_hs[(chan, proc, sys_name)]['init_met'].Fill(ev.met_corrected.pt(), record_weight) # for control
                out_hs[(chan, proc, sys_name)]['lep_pt']  .Fill(ev.lep_p4[0].pt(),  record_weight)
                out_hs[(chan, proc, sys_name)]['lep_eta'] .Fill(ev.lep_p4[0].eta(), record_weight)
                # not tagged jets
                if jets:
                    out_hs[(chan, proc, sys_name)]['jet_pt']  .Fill(jets[0][0].pt() * jets[0][1],  record_weight)
                    out_hs[(chan, proc, sys_name)]['jet_eta'] .Fill(jets[0][0].eta(), record_weight)
                # tagged jets
                if jets_b:
                    out_hs[(chan, proc, sys_name)]['bjet_pt']  .Fill(jets_b[0][0].pt() * jets_b[0][1],  record_weight)
                    out_hs[(chan, proc, sys_name)]['bjet_eta'] .Fill(jets_b[0][0].eta(), record_weight)

                # b-discr control
                out_hs[(chan, proc, sys_name)]['b_discr_lead_jet']  .Fill(lead_jet_b_discr,  record_weight)
                for discr in all_jets_b_discrs:
                    out_hs[(chan, proc, sys_name)]['b_discr_all_jets']  .Fill(discr,  record_weight)

                #out_hs[(chan, proc, sys_name)]['Mt_lep_met_f_c']    .Fill(Mt_lep_met_c, record_weight)
                #out_hs[(chan, proc, sys_name)]['Mt_lep_met_f_test'] .Fill(Mt_lep_met_test, record_weight)
                #out_hs[(chan, proc, sys_name)]['Mt_lep_met_f_cos']  .Fill(Mt_lep_met_cos, record_weight)
                #out_hs[(chan, proc, sys_name)]['Mt_lep_met_f_cos_c'].Fill(Mt_lep_met_cos_c, record_weight)
                #out_hs[(chan, proc, sys_name)]['Mt_lep_met_f_mth']  .Fill(Mt_lep_met_mth, record_weight)
                out_hs[(chan, proc, sys_name)]['Mt_lep_met_f']      .Fill(Mt_lep_met, record_weight)
                # controls for effect from weights
                out_hs[(chan, proc, sys_name)]['Mt_lep_met_f_init']      .Fill(Mt_lep_met)
                out_hs[(chan, proc, sys_name)]['Mt_lep_met_f_w_in']      .Fill(Mt_lep_met, weight)
                out_hs[(chan, proc, sys_name)]['Mt_lep_met_f_w_bf']      .Fill(Mt_lep_met, weight_bSF)
                out_hs[(chan, proc, sys_name)]['Mt_lep_met_f_w_bf_min']  .Fill(Mt_lep_met, weight_bSF_min)
                out_hs[(chan, proc, sys_name)]['Mt_lep_met_f_w_pu']      .Fill(Mt_lep_met, weight_PU)
                out_hs[(chan, proc, sys_name)]['Mt_lep_met_f_w_tp']      .Fill(Mt_lep_met, weight_top_pt)
                out_hs[(chan, proc, sys_name)]['Mt_lep_met']    .Fill(Mt_lep_met, record_weight)
                out_hs[(chan, proc, sys_name)]['dphi_lep_met']     .Fill(dphi_lep_met, record_weight)
                out_hs[(chan, proc, sys_name)]['cos_dphi_lep_met'] .Fill(cos_dphi_lep_met, record_weight)

                if pass_mumu or pass_elmu:
                    lep_lep = ev.lep_p4[0] + ev.lep_p4[1]
                    lep_lep_mass = lep_lep.mass()
                    out_hs[(chan, proc, sys_name)]['M_lep_lep']  .Fill(lep_lep_mass, record_weight)

                if has_medium_tau:
                #if len(ev.tau_p4) > 0:
                    #lep_tau = ev.lep_p4[0] + taus[0][0] # done above
                    #out_hs[(chan, proc, sys_name)]['M_lep_tau']  .Fill(lep_tau.mass(), record_weight)
                    out_hs[(chan, proc, sys_name)]['tau_pt']  .Fill(taus[0][0].pt() * taus[0][1],  record_weight)
                    out_hs[(chan, proc, sys_name)]['tau_eta'] .Fill(taus[0][0].eta(), record_weight)
                    out_hs[(chan, proc, sys_name)]['M_lep_tau']  .Fill(lep_tau_mass, record_weight)
                    out_hs[(chan, proc, sys_name)]['Mt_tau_met'] .Fill(Mt_tau_met, record_weight)
                    # we only work with 0th tau (highest pt)
                    tau_refit_index = ev.tau_refited_index[0]
                    tau_jet_index   = ev.tau_dR_matched_jet[0]
                    # require refit and dR quality of refit
                    refitted = tau_refit_index > -1 and ev.tau_SV_fit_track_OS_matched_track_dR[ev.tau_refited_index[0]] + ev.tau_SV_fit_track_SS1_matched_track_dR[ev.tau_refited_index[0]] + ev.tau_SV_fit_track_SS2_matched_track_dR[ev.tau_refited_index[0]] < 0.002
                    if refitted:
                        tau_SV_sign    = ev.tau_SV_geom_flightLenSign[tau_refit_index]
                        tau_SV_leng    = ev.tau_SV_geom_flightLen[tau_refit_index]
                        out_hs[(chan, proc, sys_name)]['tau_SV_sign'] .Fill(tau_SV_sign, record_weight)
                        out_hs[(chan, proc, sys_name)]['tau_SV_leng'] .Fill(ev.tau_SV_geom_flightLen[tau_refit_index], record_weight)
                    if tau_jet_index > -1:
                        tau_jet_bdiscr = ev.jet_b_discr[tau_jet_index]
                        out_hs[(chan, proc, sys_name)]['tau_jet_bdiscr'] .Fill(tau_jet_bdiscr, record_weight)
                    if refitted and tau_jet_index > -1:
                        out_hs[(chan, proc, sys_name)]['tau_sign_bdiscr'].Fill(tau_SV_sign, tau_jet_bdiscr, record_weight)
                        out_hs[(chan, proc, sys_name)]['tau_leng_bdiscr'].Fill(tau_SV_leng, tau_jet_bdiscr, record_weight)
                        tau_energy = taus[0][0].energy() * taus[0][1]
                        out_hs[(chan, proc, sys_name)]['tau_sign_energy'].Fill(tau_SV_sign, tau_energy, record_weight)
                        out_hs[(chan, proc, sys_name)]['tau_leng_energy'].Fill(tau_SV_leng, tau_energy, record_weight)

                #out_hs[(chan, proc, sys_name)]['Mt_lep_met_d'].Fill(Mt_lep_met_d, record_weight)
                out_hs[(chan, proc, sys_name)]['dijet_mass']  .Fill(w_mass, record_weight)
                out_hs[(chan, proc, sys_name)]['trijet_mass'] .Fill(t_mass, record_weight)
                out_hs[(chan, proc, sys_name)]['2D_dijet_trijet'] .Fill(w_mass, t_mass, record_weight)
                out_hs[(chan, proc, sys_name)]['dijet_trijet_mass'].Fill(lj_var, record_weight)

                out_hs[(chan, proc, sys_name)]['njets'].Fill(njets, record_weight)
                out_hs[(chan, proc, sys_name)]['nbjets'].Fill(nbjets, record_weight)
                out_hs[(chan, proc, sys_name)]['nvtx'].Fill(ev.nvtx, record_weight)
                if isMC:
                    out_hs[(chan, proc, sys_name)]['nvtx_gen'].Fill(ev.nvtx_gen, record_weight)

        if save_weights:
          #weight_bSF = 1.
          #weight_bSF_up, weight_bSF_down = 1., 1.
          #weight_bSF_jer_up, weight_bSF_jer_down = 1., 1.
          #weight_bSF_jes_up, weight_bSF_jes_down = 1., 1.
          control_hs['weight_z_mass_pt'] .Fill(weight_z_mass_pt)
          control_hs['weight_bSF']    .Fill(weight_bSF)
          #control_hs['weight_top_pt'] .Fill(weight_top_pt) # done above

          if pass_mu or pass_elmu or pass_mumu:
            # bcdef_weight_trk, bcdef_weight_id, bcdef_weight_iso, gh_weight_trk, gh_weight_id, gh_weight_iso
            #mu_sfs = lepton_muon_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
            #mu_trg_sf = lepton_muon_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())

            control_hs['weight_mu_trk_bcdef'].Fill(mu_sfs[0])
            control_hs['weight_mu_id_bcdef'] .Fill(mu_sfs[1])
            control_hs['weight_mu_iso_bcdef'].Fill(mu_sfs[2])
            control_hs['weight_mu_trg_bcdef'].Fill(mu_trg_sf[0])
            control_hs['weight_mu_all_bcdef'].Fill(mu_trg_sf[0] * mu_sfs[0] * mu_sfs[1] * mu_sfs[2])

            control_hs['weight_mu_trk_gh'].Fill(mu_sfs[3])
            control_hs['weight_mu_id_gh'] .Fill(mu_sfs[4])
            control_hs['weight_mu_iso_gh'].Fill(mu_sfs[5])
            control_hs['weight_mu_trg_gh'].Fill(mu_trg_sf[1])
            control_hs['weight_mu_all_gh'].Fill(mu_trg_sf[1] * mu_sfs[3] * mu_sfs[4] * mu_sfs[5])

            control_hs['weight_mu_bSF'].Fill(weight_bSF)
            control_hs['weight_mu_bSF_up']  .Fill(weight_bSF_up)
            control_hs['weight_mu_bSF_down'].Fill(weight_bSF_down)

          elif pass_el:
            #el_sfs = lepton_electron_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
            #el_trg_sf = lepton_electron_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())

            control_hs['weight_el_trk'].Fill(el_sfs[0])
            control_hs['weight_el_idd'].Fill(el_sfs[1])
            control_hs['weight_el_trg'].Fill(el_trg_sf)
            control_hs['weight_el_all'].Fill(el_trg_sf * el_sfs[0] * el_sfs[1])

            control_hs['weight_el_bSF'].Fill(weight_bSF)
            control_hs['weight_el_bSF_up']  .Fill(weight_bSF_up)
            control_hs['weight_el_bSF_down'].Fill(weight_bSF_down)

    profile.disable()
    #profile.print_stats()
    # there is no returning string
    #profile.dump_stats()

    return out_hs, control_hs, profile





#def main(input_dir, dtag, outdir, range_min, range_max):
def main(input_filename, outdir, range_min, range_max, lumi_bcdef=20263.3, lumi_gh=16518.588):
    '''main(input_filename, outdir, range_min, range_max, lumi_bcdef=20263.3, lumi_gh=16518.588)

    lumi defaults are from _full_ golden json for muon
    -- the bad lumis don't reduce it that much, should not affect the ratio too much
    '''


    f = TFile(input_filename)
    tree = f.Get('ntupler/reduced_ttree')
    if not range_max: range_max = tree.GetEntries()

    fout_name = input_filename.split('/')[-1].split('.root')[0] + "_%d-%d.root" % (range_min, range_max)
    logger_file = outdir + '/logs/' + fout_name.split('.root')[0] + '.log'

    # this dir should be made at spawning the threads
    if not os.path.exists(outdir + '/logs/'):
        os.makedirs(outdir + '/logs/')

    ## it doesn't deal with threads
    ##logger = logging.getLogger('job_processing_%d' % hash(logger_file))
    #logger = logging.getLogger()
    #hdlr = logging.FileHandler(logger_file)
    #formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    #hdlr.setFormatter(formatter)
    #logger.addHandler(hdlr) 
    ##logger.setLevel(logging.WARNING)
    #logger.setLevel(logging.INFO)
    logger = file(logger_file, 'w')

    #if '-w' in argv:
    #    input_filename, nick = '/eos/user/o/otoldaie/ttbar-leptons-80X_data/v12.7/MC2016_Summer16_WJets_amcatnlo.root', 'wjets'
    #    dtag = "MC2016_Summer16_WJets_madgraph"
    #    #init_bSF_call = 'set_bSF_effs_for_dtag("' + dtag + '");'
    #    #logger.write("init b SFs with: " + init_bSF_call)
    #    #gROOT.ProcessLine(init_bSF_call)
    #else:
    #    input_filename, nick = '/eos/user/o/otoldaie/ttbar-leptons-80X_data/v12.7/MC2016_Summer16_TTJets_powheg_1.root', 'tt'
    #    dtag = "MC2016_Summer16_TTJets_powheg"
    #    #init_bSF_call = 'set_bSF_effs_for_dtag("' + dtag + '");'
    #    #logger.write("init b SFs with: " + init_bSF_call)
    #    #gROOT.ProcessLine(init_bSF_call)

    #input_filename = input_dir + '/' + dtag + '.root'

    #dtag = input_filename.split('/')[-1].split('.')[0]
    #logger.write("dtag = " + dtag)
    logger.write("input file = %s\n" % input_filename)
    logger.write("output dir = %s\n" % outdir)
    #f = TFile('outdir/v12.3/merged-sets/MC2016_Summer16_TTJets_powheg.root')

    logger.write("N entries = %s\n" % tree.GetEntries())

    logger.write("range = %d, %d\n" % (range_min, range_max))

    logger.write("lumi BCDEF GH = %f %f\n" % (lumi_bcdef, lumi_gh))
    out_hs, c_hs, perf_profile = full_loop(tree, input_filename, lumi_bcdef, lumi_gh, range_min, range_max, logger)

    perf_profile.dump_stats(logger_file.split('.log')[0] + '.cprof')

    events_counter = f.Get('ntupler/events_counter')
    weight_counter = f.Get('ntupler/weight_counter')
    events_counter.SetDirectory(0)
    weight_counter.SetDirectory(0)

    f.Close()

    for name, h in c_hs.items():
        try:
            logger.write("%20s %9.5f\n" % (name, h.GetMean()))
        except Exception as e:
            #logger.error("%s\n%s\n%s" % (e.__class__, e.__doc__, e.message))
            logger.write("%s\n%s\n%s\n" % (e.__class__, e.__doc__, e.message))
            continue

    if with_bSF:
        logger.write("eff_b             %f\n" % h_control_btag_eff_b             .GetMean())
        logger.write("eff_c             %f\n" % h_control_btag_eff_c             .GetMean())
        logger.write("eff_udsg          %f\n" % h_control_btag_eff_udsg          .GetMean())
        logger.write("weight_b          %f\n" % h_control_btag_weight_b          .GetMean())
        logger.write("weight_c          %f\n" % h_control_btag_weight_c          .GetMean())
        logger.write("weight_udsg       %f\n" % h_control_btag_weight_udsg       .GetMean())
        logger.write("weight_notag_b    %f\n" % h_control_btag_weight_notag_b    .GetMean())
        logger.write("weight_notag_c    %f\n" % h_control_btag_weight_notag_c    .GetMean())
        logger.write("weight_notag_udsg %f\n" % h_control_btag_weight_notag_udsg .GetMean())

    '''
    for d, histos in out_hs.items():
        for name, histo in histos.items():
            print(d, name)
            histo.Print()
    '''

    # 
    #out_hs = OrderedDict([((chan, proc, sys), {'met': TH1D('met'+chan+proc+sys, '', 40, 0, 400),
    #                                           'Mt_lep_met': TH1D('Mt_lep_met'+chan+proc+sys, '', 20, 0, 200),
    #                                           #'Mt_lep_met_d': TH1D('Mt_lep_met_d'+chan+proc+sys, '', 20, 0, 200), # calculate with method of objects
    #                                           'Mt_tau_met': TH1D('Mt_tau_met'+chan+proc+sys, '', 20, 0, 200),
    #                                           'njets':  TH1D('njets'+chan+proc+sys,  '', 10, 0, 10),
    #                                           'nbjets': TH1D('nbjets'+chan+proc+sys, '', 5, 0, 5),
    #                                           'dijet_trijet_mass': TH1D('dijet_trijet_mass'+chan+proc+sys, '', 20, 0, 400) })
    #            for chan, (procs, _) in channels.items() for proc in procs for sys in systematic_names])

    #fout = TFile("lets_test.root", "RECREATE")
    #fout = TFile(outdir + '/' + dtag + "_%d-%d.root" % (range_min, range_max), "RECREATE")
    #fout_name = outdir + '/' + 
    fout = TFile(outdir + '/' + fout_name, "RECREATE")
    fout.Write()

    for (chan, proc, sys), histos in out_hs.items():
        fout.cd()
        out_dir_name = '%s/%s/%s/' % (chan, proc, sys)
        if fout.Get(out_dir_name):
            #logger.debug('found ' + out_dir_name)
            out_dir = fout.Get(out_dir_name)
        else:
            #logger.debug('made  ' + out_dir_name)
            out_dir_c = fout.Get(chan) if fout.Get(chan) else fout.mkdir(chan)
            out_dir_c.cd()
            out_dir_p = out_dir_c.Get(proc) if out_dir_c.Get(proc) else out_dir_c.mkdir(proc)
            out_dir_p.cd()
            out_dir   = out_dir_p.mkdir(sys)
        out_dir.cd()

        #out_dir.Print()
        for histo in histos.values():
            #histo.Print()
            histo.SetDirectory(out_dir)
            histo.Write()
        #out_dir.Print()

    #events_counter.SetDirectory()
    #weight_counter.SetDirectory()
    fout.cd()
    events_counter.Write() # hopefully these go to the root of the tfile
    weight_counter.Write()

    fout.Write()


