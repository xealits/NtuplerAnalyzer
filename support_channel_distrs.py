import os
from os import environ
from array import array
from collections import OrderedDict
import cProfile
import logging as Logging


logging = Logging.getLogger("common")

logging.info('importing ROOT')
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, gROOT, gSystem, TCanvas, TGraphAsymmErrors, TMath, TString

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
    # sqrt(s)         channel     	a     	b
    # 7 TeV         all combined 	0.199 	-0.00166
    # 7 TeV         l+jets      	0.174 	-0.00137
    # 7 TeV         dilepton    	0.222 	-0.00197
    # 8 TeV         all combined 	0.156 	-0.00137
    # 8 TeV         l+jets       	0.159 	-0.00141
    # 8 TeV         dilepton     	0.148 	-0.00129
    # 13 TeV        all combined	0.0615	-0.0005
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
    if len(b_jets) == 0 or len(light_jets) == 0: return 2000, 2000, 2000
    # loop over light jets -- check their mass and mass of their pairs to be close to 80
    # find the closest vector
    # then loop over b jets finding the combination with mass closest to 173

    # light jets close to W
    dist_W = 99999
    # closest_vector
    for j, mult in light_jets:
        new_dist = abs(j.mass() * mult - 80) # DANGER: the LorentzVector is used, .M() is spacial magnitude
        if new_dist < dist_W:
	    dist_W = new_dist
	    closest_to_W = j * mult

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
    dist_t = 99999
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

    save_control = False
    isMC = 'MC' in dtag
    aMCatNLO = 'amcatnlo' in dtag
    isTT = 'TT' in dtag
    isSTop = 'SingleT' in dtag or 'tchannel' in dtag or 'schannel' in dtag
    isDY = 'DY' in dtag
    isWJets = 'WJet' in dtag or 'W1Jet' in dtag or 'W2Jet' in dtag or 'W3Jet' in dtag or 'W4Jet' in dtag
    isQCD = 'QCD' in dtag
    isDibosons = 'WW' in dtag or 'ZZ' in dtag or 'WZ' in dtag

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

    tt_el_procs = ['tt_eltau', 'tt_lj', 'tt_other']
    tt_mu_procs = ['tt_mutau', 'tt_lj', 'tt_other']

    if isMC:
        if isTT:
            channels = {'el_presel': (['tt_eltau', 'tt_lj', 'tt_other'], 'tt_other'),
                   'el_sel':        (['tt_eltau', 'tt_lj', 'tt_other'], 'tt_other'),
                   'el_lj':     (['tt_eltau', 'tt_lj', 'tt_other'], 'tt_other'),
                   'el_lj_out': (['tt_eltau', 'tt_lj', 'tt_other'], 'tt_other'),
                   'mu_presel': (['tt_mutau', 'tt_lj', 'tt_other'], 'tt_other'),
                   'mu_sel':        (['tt_mutau', 'tt_lj', 'tt_other'], 'tt_other'),
                   'mu_lj':     (['tt_mutau', 'tt_lj', 'tt_other'], 'tt_other'),
                   'mu_lj_out': (['tt_mutau', 'tt_lj', 'tt_other'], 'tt_other')}
            process = 'tt_other'

        if isWJets:
            channels = {'el_presel': (['wjets'], 'wjets'),
                      'el_sel':        (['wjets'], 'wjets'),
                      'el_lj':     (['wjets'], 'wjets'),
                      'el_lj_out': (['wjets'], 'wjets'),
                      'mu_presel': (['wjets'], 'wjets'),
                      'mu_sel':        (['wjets'], 'wjets'),
                      'mu_lj':     (['wjets'], 'wjets'),
                      'mu_lj_out': (['wjets'], 'wjets')}
            process = 'wjets'

        if isDY:
            channels = {'el_presel': (['dy'], 'dy'),
                   'el_sel':        (['dy'], 'dy'),
                   'el_lj':     (['dy'], 'dy'),
                   'el_lj_out': (['dy'], 'dy'),
                   'mu_presel': (['dy'], 'dy'),
                   'mu_sel':        (['dy'], 'dy'),
                   'mu_lj':     (['dy'], 'dy'),
                   'mu_lj_out': (['dy'], 'dy')}
            process = 'dy'

        if isSTop:
            channels = {'el_presel': (['singletop'], 'singletop'),
                   'el_sel':          (['singletop'], 'singletop'),
                   'el_lj':       (['singletop'], 'singletop'),
                   'el_lj_out':   (['singletop'], 'singletop'),
                   'mu_presel':   (['singletop'], 'singletop'),
                   'mu_sel':          (['singletop'], 'singletop'),
                   'mu_lj':       (['singletop'], 'singletop'),
                   'mu_lj_out':   (['singletop'], 'singletop')}
            process = 'singletop'

        if isQCD:
            channels = {'el_presel': (['qcd'], 'qcd'),
                   'el_sel':          (['qcd'], 'qcd'),
                   'el_lj':       (['qcd'], 'qcd'),
                   'el_lj_out':   (['qcd'], 'qcd'),
                   'mu_presel':   (['qcd'], 'qcd'),
                   'mu_sel':          (['qcd'], 'qcd'),
                   'mu_lj':       (['qcd'], 'qcd'),
                   'mu_lj_out':   (['qcd'], 'qcd')}
            process = 'qcd'

        if isDibosons:
            channels = {'el_presel': (['dibosons'], 'dibosons'),
                   'el_sel':          (['dibosons'], 'dibosons'),
                   'el_lj':       (['dibosons'], 'dibosons'),
                   'el_lj_out':   (['dibosons'], 'dibosons'),
                   'mu_presel':   (['dibosons'], 'dibosons'),
                   'mu_sel':          (['dibosons'], 'dibosons'),
                   'mu_lj':       (['dibosons'], 'dibosons'),
                   'mu_lj_out':   (['dibosons'], 'dibosons')}
            process = 'dibosons'

    else:
        channels = {'el_presel': (['data'], 'data'),
                      'el_sel':        (['data'], 'data'),
                      'el_lj':     (['data'], 'data'),
                      'el_lj_out': (['data'], 'data'),
                      'mu_presel': (['data'], 'data'),
                      'mu_sel':        (['data'], 'data'),
                      'mu_lj':     (['data'], 'data'),
                      'mu_lj_out': (['data'], 'data')}
        process = 'data'

    proc = process

    if isMC:
        systematic_names = ['NOMINAL',
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
    else:
        systematic_names = ['NOMINAL']

    if isTT:
        systematic_names.append('TOPPTUp')
        systematic_names.append('TOPPTDown')

    # channel -- reco selection
    # proc    -- MC gen info, like inclusive tt VS tt->mutau and others,
    #            choice of subprocesses depends on channel (sadly),
    #            (find most precise subprocess ones, then store accordingly for channels)
    # sys     -- shape systematics
    out_hs = OrderedDict([((chan, proc, sys), {'met':        TH1D('%s_%s_%s_met' % (chan, sys, proc), '', 40, 0, 400),
                                               'Mt_lep_met': TH1D('%s_%s_%s_Mt_lep_met' % (chan, sys, proc), '', 20, 0, 200),
                                               #'Mt_lep_met_d': TH1D('Mt_lep_met_d'+chan+proc+sys, '', 20, 0, 200), # calculate with method of objects
                                               'Mt_tau_met': TH1D('%s_%s_%s_Mt_tau_met' % (chan, sys, proc), '', 20, 0, 200),
                                               'njets':      TH1D('%s_%s_%s_njets' % (chan, sys, proc),  '', 10, 0, 10),
                                               'nbjets':     TH1D('%s_%s_%s_nbjets' % (chan, sys, proc), '', 5, 0, 5),
                                               'dijet_trijet_mass': TH1D('%s_%s_%s_dijet_trijet_mass' % (chan, sys, proc), '', 20, 0, 400) })
                for chan, (procs, _) in channels.items() for proc in procs for sys in systematic_names])

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

        # only 1-lep channels
        if not (pass_mu or pass_el): continue

        # also at least some kind of tau:
        if not ev.tau_p4.size() > 0: continue

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
                ev.nalljets  # number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
                )
            met_y = ROOT.met_pt_recoilcor_y(
                ev.met_corrected.Px(), # uncorrected type I pf met px (float)
                ev.met_corrected.Py(), # uncorrected type I pf met py (float)
                ev.gen_genPx, # generator Z/W/Higgs px (float)
                ev.gen_genPy, # generator Z/W/Higgs py (float)
                ev.gen_visPx, # generator visible Z/W/Higgs px (float)
                ev.gen_visPy, # generator visible Z/W/Higgs py (float)
                ev.nalljets  # number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
                )

            #Mt_lep_met = transverse_mass_pts(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), ev.pfmetcorr_ex, ev.pfmetcorr_ey)
            #Mt_tau_met_nominal = transverse_mass_pts(ev.tau_p4[0].Px(), ev.tau_p4[0].Py(), ev.pfmetcorr_ex, ev.pfmetcorr_ey)
            # TODO: tau is corrected with systematic ES
        else:
            met_x = ev.met_corrected.Px()
            met_y = ev.met_corrected.Py()
            #Mt_lep_met = transverse_mass(ev.lep_p4[0], ev.met_corrected)
            #Mt_tau_met = transverse_mass(ev.tau_p4[0], ev.met_corrected)
            #Mt_tau_met_nominal = transverse_mass_pts(ev.tau_p4[0].Px(), ev.tau_p4[0].Py(), ev.met_corrected.Px(), ev.met_corrected.Py())
        # also
        #Mt_lep_met_d = (ev.lep_p4[0] + ev.met_corrected).Mt()


        # TODO: pass jet systematics to met?
        Mt_lep_met = transverse_mass_pts(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), met_x, met_y)
        met_pt = TMath.Sqrt(met_x*met_x + met_y*met_y)

        weight = 1. # common weight of event (1. for data)
        if isMC:
            try:
                weight_pu    = pileup_ratio[ev.nvtx]
                weight_pu_up = pileup_ratio_up[ev.nvtx]
                weight_pu_dn = pileup_ratio_down[ev.nvtx]
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

            weight_top_pt = 1.
            if isTT:
                weight_top_pt = ttbar_pT_SF(ev.gen_t_pt, ev.gen_tb_pt)
                #weight *= weight_top_pt # to sys
                control_hs['weight_top_pt']   .Fill(weight_top_pt)
                # basically only difference is eltau/mutau
                if (abs(ev.gen_t_w_decay_id) > 15*15 and abs(ev.gen_tb_w_decay_id) == 13) or (abs(ev.gen_t_w_decay_id) == 13 and abs(ev.gen_tb_w_decay_id) > 15*15): # lt
                    proc = 'tt_mutau'
                elif (abs(ev.gen_t_w_decay_id) > 15*15 and abs(ev.gen_tb_w_decay_id) == 11) or (abs(ev.gen_t_w_decay_id) == 11 and abs(ev.gen_tb_w_decay_id) > 15*15): # lt
                    proc = 'tt_eltau'
                elif abs(ev.gen_t_w_decay_id * ev.gen_tb_w_decay_id) == 13 or abs(ev.gen_t_w_decay_id * ev.gen_tb_w_decay_id) == 11: # lj
                    proc = 'tt_lj'
                else:
                    proc = 'tt_other'

            if pass_mu and isMC:
                mu_sfs = lepton_muon_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                mu_trg_sf = lepton_muon_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                # bcdef gh eras
                weight *= ratio_bcdef * mu_trg_sf[0] * mu_sfs[0] * mu_sfs[1] * mu_sfs[2] + ratio_gh * mu_trg_sf[1] * mu_sfs[3] * mu_sfs[4] * mu_sfs[5]

            if pass_el and isMC:
                el_sfs = lepton_electron_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                el_trg_sf = lepton_electron_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
                weight *= el_trg_sf * el_sfs[0] * el_sfs[1]

        #    #JES corrections
        #    # TODO: save only 30 pt 2.5 eta jets
        #    #[((1-f)*j.pt(), (1+f)*j.pt()) for f,j in zip(ev.jet_jes_correction_relShift, ev.jet_p4)]
        #    jet_pts_jes_up, jet_pts_jes_down = [], []
        #    #for f,j in zip(ev.jet_jes_correction_relShift, ev.jet_p4):
        #    #    jet_pts_jes_up  .append(j.pt()*(1-f))
        #    #    jet_pts_jes_down.append(j.pt()*(1+f))
        #    #JER
        #    #[(j.pt()*(down/factor), j.pt()*(up/factor)) for j, factor, up, down in zip(ev.jet_p4, ev.jet_jer_factor, ev.jet_jer_factor_up, ev.jet_jer_factor_down)]
        #    jet_pts_jer_up, jet_pts_jer_down = [], []
        #    #for j, factor, up, down in zip(ev.jet_p4, ev.jet_jer_factor, ev.jet_jer_factor_up, ev.jet_jer_factor_down):
        #    for i in range(ev.jet_p4.size()):
        #        j, factor, up, down = ev.jet_p4[i], ev.jet_jer_factor[i], ev.jet_jer_factor_up[i], ev.jet_jer_factor_down[i]
        #        jet_pts_jer_up  .append(j.pt()*(up/factor))
        #        jet_pts_jer_down.append(j.pt()*(down/factor))
        #        jes_shift = ev.jet_jes_correction_relShift[i]
        #        jet_pts_jes_up  .append(j.pt()*(1-jes_shift))
        #        jet_pts_jes_down.append(j.pt()*(1+jes_shift))
        ## implement lj_var?

        #jet_pts = [] # nominal jet pts
        #for j in ev.jet_p4:
        #    jet_pts.append(j.pt())

        ## number of b-tagged jets is not varied, it's the same for given jets
        ## TODO: however, the N b jets should be calculated for JER/JES varied jets
        ##       also the variation should be propagated to MET
        #N_bjets = 0
        #weight_bSF, weight_bSF_up, weight_bSF_down = 1., 1., 1.
        ##for i, (p4, flavId, b_discr) in enumerate(zip(ev.jet_p4, ev.jet_hadronFlavour, ev.jet_b_discr)):
        #for i in range(ev.jet_p4.size()):
        #    p4, flavId, b_discr = ev.jet_p4[i], ev.jet_hadronFlavour[i], ev.jet_b_discr[i]
        #    # calc_btag_sf_weight(hasCSVtag: bool, flavId: int, pt: float, eta: float) -> float:
        #    if p4.pt() < 30: continue
        #    N_bjets += b_discr > 0.8484
        #    if isMC and with_bSF:
        #        weight_bSF      *= calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta())
        #        weight_bSF_up   *= calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "up")
        #        weight_bSF_down *= calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "down")

        # I need from jets:
        #    cut [pt?], eta and PFID == Loose
        #    save jes/jer factors (if MC) <---- needed later only for lj_var
        #    find b tags for these jets
        #    and b SF factors (if MC)
        #
        #    get N jets passing selection for diff syst!
        #    N corresp b jets
        #    and b SF weights
        #    build lj_var for every variation of N jets (N b jets doesn't change)
        #
        # -- does pyROOT store pointers for everything, like Python, or it copies stuff?
        # test with TLorentzVector shows that it copies, as it should
        # therefore append stuff
        #
        # for each jet syst make ([p4-s], [sys_factors]) for jets passing eta, ID, pt
        # find N b-jets and SF weight
        # and for nominal collection of jets -- b SF up/down weight
        # build lj only for events passing the selection -- it's expensive

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

        for i in xrange(ev.jet_p4.size()):
            pfid, p4 = ev.jet_PFID[i], ev.jet_p4[i]
            if pfid < 1 or abs(p4.eta()) > 2.5: continue # Loose PFID and eta

            flavId = ev.jet_hadronFlavour[i]
            b_tagged = ev.jet_b_discr[i] > b_tag_wp

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

            if isMC:
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

        has_tau = False
        has_tau_es_up = False
        has_tau_es_down = False

        if ev.tau_p4.size() > 0 and ev.tau_IDlev[0] > 2 and abs(ev.tau_p4[0].eta()) < 2.4:
            # it should work like Python does and not copy these objects! (cast)
            #p4, DM, IDlev = ev.tau_p4[0], ev.tau_decayMode[0], ev.tau_IDlev[0]
            #if IDlev < 3 or abs(p4.eta()) < 2.4: continue # only Medium taus
            DM = ev.tau_decayMode[0]
            if DM == 0:
              factor = 0.995
              if isMC:
                factor_up   = 0.995 + 0.012
                factor_down = 0.995 - 0.012
            elif DM < 10:
              factor = 1.011
              if isMC:
                factor_up   = 1.011 + 0.012
                factor_down = 1.011 - 0.012
            else:
              factor = 1.006
              if isMC:
                factor_up   = 1.006 + 0.012
                factor_down = 1.006 - 0.012

            has_tau = (ev.tau_p4[0].pt() * factor) > 30
            if not Mt_tau_met_nominal:
	        Mt_tau_met_nominal = transverse_mass_pts(ev.tau_p4[0].Px()*factor, ev.tau_p4[0].Py()*factor, met_x, met_y)
            if isMC:
                has_tau_es_up   = (ev.tau_p4[0].pt() * factor_up  ) > 30
                has_tau_es_down = (ev.tau_p4[0].pt() * factor_down) > 30
                if not Mt_tau_met_up:
                    Mt_tau_met_up   = transverse_mass_pts(ev.tau_p4[0].Px()*factor_up, ev.tau_p4[0].Py()*factor_up, met_x, met_y)
                    Mt_tau_met_down = transverse_mass_pts(ev.tau_p4[0].Px()*factor_down, ev.tau_p4[0].Py()*factor_down, met_x, met_y)

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
            systematics = {'NOMINAL'   : [jets_nominal,  jets_b_nominal,  has_tau, weight_bSF,          weight_pu,    1, Mt_tau_met_nominal],
                           'JESUp'     : [jets_jes_up,   jets_b_jes_up,   has_tau, weight_bSF_jes_up,   weight_pu,    1, Mt_tau_met_nominal],
                           'JESDown'   : [jets_jes_down, jets_b_jes_down, has_tau, weight_bSF_jes_down, weight_pu,    1, Mt_tau_met_nominal],
                           'JERUp'     : [jets_jer_up,   jets_b_jer_up,   has_tau, weight_bSF_jer_up,   weight_pu,    1, Mt_tau_met_nominal],
                           'JERDown'   : [jets_jer_down, jets_b_jer_down, has_tau, weight_bSF_jer_down, weight_pu,    1, Mt_tau_met_nominal],
                           'TauESUp'   : [jets_nominal,  jets_b_nominal,  has_tau_es_up  , weight_bSF,     weight_pu,    1, Mt_tau_met_up],
                           'TauESDown' : [jets_nominal,  jets_b_nominal,  has_tau_es_down, weight_bSF,     weight_pu,    1, Mt_tau_met_down],
                           'bSFUp'     : [jets_nominal,  jets_b_nominal,  has_tau, weight_bSF_up  ,     weight_pu,    1, Mt_tau_met_nominal],
                           'bSFDown'   : [jets_nominal,  jets_b_nominal,  has_tau, weight_bSF_down,     weight_pu,    1, Mt_tau_met_nominal],
                           'PUUp'      : [jets_nominal,  jets_b_nominal,  has_tau, weight_bSF,          weight_pu_up, 1, Mt_tau_met_nominal],
                           'PUDown'    : [jets_nominal,  jets_b_nominal,  has_tau, weight_bSF,          weight_pu_dn, 1, Mt_tau_met_nominal],
                           }
            if isTT:
                systematics['TOPPTUp'] = [jets_nominal, jets_b_nominal, has_tau, weight_bSF, weight_pu, weight_top_pt, Mt_tau_met_nominal]
                systematics['TOPPTDown'] = [jets_nominal, jets_b_nominal, has_tau, weight_bSF, weight_pu, 1., Mt_tau_met_nominal]
        else:
            systematics = {'NOMINAL': [jets_nominal, jets_b_nominal, has_tau, 1., 1., 1., Mt_tau_met_nominal]}

        # for each systematic
        # pass one of the reco selections
        # check the subprocess
        # store distr

        for sys_name, (jets, jets_b, has_medium_tau, weight_bSF, weight_PU, weight_top_pt, Mt_tau_met) in systematics.items():
            sys_weight = weight * weight_bSF * weight_PU * weight_top_pt
            # pass reco selections

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
            njets =  nbjets + len(jets)
            has_bjets = len(jets_b) > 0
            has_3jets = (len(jets) + len(jets_b)) > 2
            #has_pre_tau = len(tau_pts) > 0
            #has_medium_tau = has_pre_tau and tau_pts[0] > 30 and ev.tau_IDlev[0] > 2
            #has_medium_tau = tau_pts[0] > 30 and ev.tau_IDlev[0] > 2
            os_lep_med_tau = has_medium_tau and ev.tau_id[0]*ev.lep_id[0] < 0
            # dilep_mass is not TES corrected
            # it's used only in control regions and shouldn't be a big deal
            #dy_dilep_mass = has_pre_tau and (45 < (ev.lep_p4[0] + ev.tau_p4[0]).mass() < 85)

            pass_single_lep_presel = large_met and has_3jets and has_bjets #and os_lep_med_tau
            pass_single_lep_sel = pass_single_lep_presel and os_lep_med_tau
            #pass_mujets    = pass_mu and has_pre_tau and large_Mt_lep
            #pass_mujets_b  = pass_mu and has_pre_tau and large_Mt_lep and has_bjets
            #pass_taumutauh = pass_mu and os_lep_med_tau and dy_dilep_mass and not large_Mt_lep and not has_bjets
            #pass_antiMt_taumutauh_pretau_allb = pass_mu and not dy_dilep_mass and large_Mt_lep
            #pass_antiMt_taumutauh_pretau      = pass_mu and not dy_dilep_mass and large_Mt_lep and not has_bjets
            #pass_antiMt_taumutauh             = pass_mu and not dy_dilep_mass and large_Mt_lep and not has_bjets and has_medium_tau
            #pass_antiMt_taumutauh_OS          = pass_mu and not dy_dilep_mass and large_Mt_lep and not has_bjets and os_lep_med_tau

            lj_cut = 50
            passed_channels = []
            if pass_el:
                if pass_single_lep_presel:
		    passed_channels.append('el_presel')
		    # calc lj_var
		    lj_var, w_mass, t_mass = calc_lj_var(jets, jets_b)
                    large_lj = lj_var > lj_cut # TODO: implement lj_var
                if pass_single_lep_sel:
		    passed_channels.append('el_sel')
                    if not large_lj: passed_channels.append('el_lj')
                    if large_lj: passed_channels.append('el_lj_out')
            else:
                if pass_single_lep_presel:
		    passed_channels.append('mu_presel')
		    lj_var, w_mass, t_mass = calc_lj_var(jets, jets_b)
                    large_lj = lj_var > lj_cut # TODO: implement lj_var
                if pass_single_lep_sel:
		    passed_channels.append('mu_sel')
                    if not large_lj: passed_channels.append('mu_lj')
                    if large_lj: passed_channels.append('mu_lj_out')

            # now the process variety is restricted to 2 separate channels and TT
            if isTT and pass_el:
                proc = proc if proc in tt_el_procs else 'tt_other'
            elif isTT and pass_mu:
                proc = proc if proc in tt_mu_procs else 'tt_other'

            ## save distrs
            #for chan, proc in chan_subproc_pairs:
            #    out_hs[(chan, proc, sys_name)]['met'].Fill(met_pt, weight)
            #    out_hs[(chan, proc, sys_name)]['Mt_lep_met'].Fill(Mt_lep_met, weight)
            #    out_hs[(chan, proc, sys_name)]['Mt_lep_met_d'].Fill(Mt_lep_met_d, weight)
            #    out_hs[(chan, proc, sys_name)]['dijet_trijet_mass'].Fill(25, weight)
            for chan in passed_channels:
                out_hs[(chan, proc, sys_name)]['met'].Fill(met_pt, sys_weight)
                out_hs[(chan, proc, sys_name)]['Mt_tau_met'].Fill(Mt_tau_met, sys_weight)
                out_hs[(chan, proc, sys_name)]['Mt_lep_met'].Fill(Mt_lep_met, sys_weight)
                #out_hs[(chan, proc, sys_name)]['Mt_lep_met_d'].Fill(Mt_lep_met_d, sys_weight)
                out_hs[(chan, proc, sys_name)]['dijet_trijet_mass'].Fill(lj_var, sys_weight)
                out_hs[(chan, proc, sys_name)]['njets'].Fill(njets, sys_weight)
                out_hs[(chan, proc, sys_name)]['nbjets'].Fill(nbjets, sys_weight)

        if save_control:
          if pass_mu:
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

            for i, (p4, flavId, b_discr) in enumerate(zip(ev.jet_p4, ev.jet_hadronFlavour, ev.jet_b_discr)):
                # calc_btag_sf_weight(hasCSVtag: bool, flavId: int, pt: float, eta: float) -> float:
                control_hs['weight_mu_bSF'].Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta()))
                control_hs['weight_mu_bSF_up']  .Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "up"))
                control_hs['weight_mu_bSF_down'].Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "down"))
            # I do need a working standalone b-tag calibrator instead of this jogling
            #b_taggin_SF(jet0_p4.pt(), jet0_p4.eta(), jet0_b_discr, jet0_hadronFlavour, 0.8484)

          elif pass_el:
            #el_sfs = lepton_electron_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
            #el_trg_sf = lepton_electron_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())

            control_hs['weight_el_trk'].Fill(el_sfs[0])
            control_hs['weight_el_idd'].Fill(el_sfs[1])
            control_hs['weight_el_trg'].Fill(el_trg_sf)
            control_hs['weight_el_all'].Fill(el_trg_sf * el_sfs[0] * el_sfs[1])

            for i, (p4, flavId, b_discr) in enumerate(zip(ev.jet_p4, ev.jet_hadronFlavour, ev.jet_b_discr)):
                # calc_btag_sf_weight(hasCSVtag: bool, flavId: int, pt: float, eta: float) -> float:
                control_hs['weight_el_bSF'].Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta()))
                control_hs['weight_el_bSF_up']  .Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "up"))
                control_hs['weight_el_bSF_down'].Fill(calc_btag_sf_weight(b_discr > 0.8484, flavId, p4.pt(), p4.eta(), "down"))

    profile.disable()
    #profile.print_stats()
    # there is no returning string
    #profile.dump_stats()

    return out_hs, control_hs, profile





#def main(input_dir, dtag, outdir, range_min, range_max):
def main(input_filename, outdir, range_min, range_max):
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
    #f = TFile('outdir/v12.3/merged-sets/MC2016_Summer16_TTJets_powheg.root')

    logger.write("N entries = %s\n" % tree.GetEntries())

    logger.write("range = %d, %d\n" % (range_min, range_max))
    out_hs, c_hs, perf_profile = full_loop(tree, input_filename, 0, 6175, range_min, range_max, logger)

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


