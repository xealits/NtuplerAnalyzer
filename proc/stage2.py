from datetime import datetime
import os
from os import environ
from os.path import isfile
from array import array
from collections import OrderedDict, namedtuple, defaultdict
import cProfile
import logging as Logging
import ctypes
from ctypes import CDLL, cdll, c_double, c_int, POINTER
from random import random as u_random
from gen_proc_defs import *
import math

#import pdb

OLD_MINIAOD_JETS = False
DO_W_STITCHING = False
ALL_JETS = False
with_bSF = True
WITH_RECOIL_CORRECTIONS = True
N_RECOIL_JETS = 0

ISO_LEPS    = True
JETS_PT_CUT = 30. # 21. # 
TAUS_PT_CUT = 30. # 21. # # 20GeV for DY->tautau selection
TAUS_ID_CUT_Tight  = 3  # Tight = 4-5, VTight = 5
TAUS_ID_CUT_Medium = 2
TAUS_ID_CUT_VLoose = 0
TAUS_ID_CUT = TAUS_ID_CUT_VLoose # TAUS_ID_CUT_Medium # Vloose cut for the shape
ONLY_3PI_TAUS = False
SV_SIGN_CUT = 2.5

METMuEGClean = True # the met for data
PROP_TAU     = True
REMOVE_LEPJET     = False
PROP_UNCORLEPJET  = False
PROP_LEPJET  = False
PROP_LEPJET_UNCLUSTER = False
PROP_JETS    = True

REQUIRE_MLB = False

def PASSES_FUNC(pass_mu, pass_elmu, pass_elmu_el, pass_mumu, pass_elel, pass_el, pass_mu_all, pass_el_all):
    #return pass_mu or pass_el or pass_elmu or pass_mu_all or pass_el_all or pass_mumu
    #return pass_elmu
    return pass_mu or pass_el or pass_elmu_el or pass_mumu or pass_elel


def passes_wjets_control_selection(passed_triggers, leps, N_jets, taus, proc_met):
    channel_stage = 0
    pass_mu, pass_elmu, pass_elmu_el, pass_mumu, pass_elel, pass_el, pass_mu_all, pass_el_all = passed_triggers

    #old_jet_sel = len(jets.medium) > 0 and (len(jets.taumatched[0]) + len(jets.taumatched[1]) + len(jets.medium) + len(jets.loose) + len(jets.rest)) > 2
    # N b-tagged, N all jets
    old_jet_sel = N_jets[0] == 0 and len(leps[0]) == 1 and (pass_mu or pass_el)
    old_jet_sel_tau_cand    = old_jet_sel and len(taus) > 0
    old_jet_sel_tau_cand_os = old_jet_sel_tau_cand and leps[4][0] * taus[0][2] < 0

    if   pass_mu and old_jet_sel_tau_cand_os:
        channel_stage = 7
    elif pass_mu and old_jet_sel_tau_cand:
        channel_stage = 6

    elif pass_el and old_jet_sel_tau_cand_os:
        channel_stage = 17
    elif pass_el and old_jet_sel_tau_cand:
        channel_stage = 16

    return channel_stage

# N_jets 
def passes_elmu_selection_stages(passed_triggers, leps, N_jets, taus, proc_met):
    channel_stage = 0
    pass_mu, pass_elmu, pass_elmu_el, pass_mumu, pass_elel, pass_el, pass_mu_all, pass_el_all = passed_triggers
    main_tt_elmu_trig = pass_elmu

    pass_elmu_leps = main_tt_elmu_trig and len(leps[4]) == 2 and leps[4][0] * leps[4][1] == -11*13
    pass_elmu_leps_jets = pass_elmu_leps and N_jets[1] > 1
    pass_elmu_leps_jets_bjet = pass_elmu_leps_jets and N_jets[0] > 0

    if   pass_elmu_leps_jets_bjet:
        channel_stage = 205
    elif pass_elmu_leps_jets:
        channel_stage = 204
    elif pass_elmu_leps:
        channel_stage = 203
    elif main_tt_elmu_trig:
        channel_stage = 202

    return channel_stage

def passes_tt_preselection_stages(passed_triggers, leps, N_jets, taus, proc_met):
    n_b_notau, n_b_all, n_all_jets = N_jets
    channel_stage = 0
    pass_mu, pass_elmu, pass_elmu_el, pass_mumu, pass_elel, pass_el, pass_mu_all, pass_el_all = passed_triggers

    #old_jet_sel = len(jets.medium) > 0 and (len(jets.taumatched[0]) + len(jets.taumatched[1]) + len(jets.medium) + len(jets.loose) + len(jets.rest)) > 2
    # N b-tagged, N all jets
    old_jet_presel      = n_b_all   > 0 and n_all_jets > 2 and len(leps[0]) == 1
    has_tau_cand = len(taus) > 0
    old_jet_presel_cand = n_b_notau > 0 and n_all_jets > 2 and len(leps[0]) == 1 and has_tau_cand
    tau_cand_os = has_tau_cand and leps[4][0] * taus[0][2] < 0

    if   pass_mu and old_jet_presel_cand and tau_cand_os:
        channel_stage = 9
    elif pass_mu and old_jet_presel_cand:
        channel_stage = 8

    elif pass_mu and old_jet_presel and tau_cand_os:
        channel_stage = 7
    elif pass_mu and old_jet_presel and has_tau_cand:
        channel_stage = 6
    elif pass_mu and old_jet_presel:
        channel_stage = 5

    elif pass_el and old_jet_presel_cand and tau_cand_os:
        channel_stage = 19
    elif pass_el and old_jet_presel_cand:
        channel_stage = 18

    elif pass_el and old_jet_presel and tau_cand_os:
        channel_stage = 17
    elif pass_el and old_jet_presel and has_tau_cand:
        channel_stage = 16
    elif pass_el and old_jet_presel:
        channel_stage = 15

    return channel_stage


def passes_tt_selection_stages(passed_triggers, leps, N_jets, taus, proc_met):
    channel_stage = 0
    pass_mu, pass_elmu, pass_elmu_el, pass_mumu, pass_elel, pass_el, pass_mu_all, pass_el_all = passed_triggers

    #old_jet_sel = len(jets.medium) > 0 and (len(jets.taumatched[0]) + len(jets.taumatched[1]) + len(jets.medium) + len(jets.loose) + len(jets.rest)) > 2
    # N b-tagged, N all jets
    old_jet_sel = N_jets[0] > 0 and N_jets[1] > 2 and len(leps[0]) == 1
    old_jet_sel_tau_cand    = old_jet_sel and len(taus) > 0
    old_jet_sel_tau_cand_os = old_jet_sel_tau_cand and leps[4][0] * taus[0][2] < 0

    #pass_old_presel    = old_jet_sel and len(taus_candidates) > 0
    #pass_old_presel_os = pass_old_presel and leps[4][0] * taus_candidates[0][2] < 0
    #taus_main.append((p4, (TES_factor, TES_factor_up, TES_factor_dn), tau_pdgID, i, tau_ID, jetmatched))
    tau_IDlev = 0
    if len(taus) > 0:
        tau_IDlev = taus[0][4]

    pass_Tight_sel     = old_jet_sel    and tau_IDlev > TAUS_ID_CUT_Tight
    pass_Tight_sel_os  = pass_Tight_sel and leps[4][0] * taus[0][2] < 0
    pass_old_sel       = old_jet_sel and tau_IDlev > TAUS_ID_CUT_Medium
    pass_old_sel_os    = pass_old_sel and leps[4][0] * taus[0][2] < 0
    pass_old_selVLoose       = old_jet_sel and tau_IDlev > TAUS_ID_CUT_VLoose
    pass_old_selVLoose_os    = pass_old_selVLoose and leps[4][0] * taus[0][2] < 0
    #lep_p4, lep_relIso, lep_matching_gen, lep_matching_gen_dR, lep_id = leps
    # taus_nom.medium.append((p4, TES_factor, tau_pdgID, i, jetmatched))

    if   pass_mu and pass_Tight_sel_os:
        channel_stage = 9
    elif pass_mu and pass_Tight_sel:
        channel_stage = 8
    elif pass_mu and pass_old_sel_os:
        channel_stage = 7
    elif pass_mu and pass_old_sel:
        channel_stage = 6
    elif pass_mu and pass_old_selVLoose_os:
        channel_stage = 5
    elif pass_mu and pass_old_selVLoose:
        channel_stage = 4
    elif pass_mu and old_jet_sel_tau_cand_os:
        channel_stage = 3
    elif pass_mu and old_jet_sel_tau_cand:
        channel_stage = 2
    #elif pass_mu and pass_old_presel_os:
    #    channel_stage = 3
    #elif pass_mu and pass_old_presel:
    #    channel_stage = 2

    elif pass_el and pass_Tight_sel_os:
        channel_stage = 19
    elif pass_el and pass_Tight_sel:
        channel_stage = 18
    elif pass_el and pass_old_sel_os:
        channel_stage = 17
    elif pass_el and pass_old_sel:
        channel_stage = 16
    elif pass_el and pass_old_selVLoose_os:
        channel_stage = 15
    elif pass_el and pass_old_selVLoose:
        channel_stage = 14
    elif pass_el and old_jet_sel_tau_cand_os:
        channel_stage = 13
    elif pass_el and old_jet_sel_tau_cand:
        channel_stage = 12
    #elif pass_el and pass_old_presel_os:
    #    channel_stage = 13
    #elif pass_el and pass_old_presel:
    #    channel_stage = 12

    return channel_stage

def passes_tt_selection_stages_alliso(passed_triggers, leps, N_jets, taus, proc_met):
    '''
    pass alliso leps and tau candidates, the N of jets cross-dR with alliso leptons
    '''
    channel_stage = 0
    pass_mu, pass_elmu, pass_elmu_el, pass_mumu, pass_elel, pass_el, pass_mu_all, pass_el_all = passed_triggers

    #old_jet_sel = len(jets.medium) > 0 and (len(jets.taumatched[0]) + len(jets.taumatched[1]) + len(jets.medium) + len(jets.loose) + len(jets.rest)) > 2
    # N b-tagged, N all jets
    n_bjets, n_jets = N_jets
    old_jet_sel = len(leps[0]) == 1 and n_bjets > 0 and n_jets > 2
    #pass_old_presel    = old_jet_sel and len(taus_candidates) > 0
    #pass_old_presel_os = pass_old_presel and leps[4][0] * taus_candidates[0][2] < 0
    pass_old_sel       = old_jet_sel and len(taus) > 0
    pass_old_sel_os    = pass_old_sel and leps[4][0] * taus[0][2] < 0
    #lep_p4, lep_relIso, lep_matching_gen, lep_matching_gen_dR, lep_id = leps
    # taus_nom.medium.append((p4, TES_factor, tau_pdgID, i, jetmatched))

    if   pass_mu_all and pass_old_sel_os:
        channel_stage = 503
    elif pass_mu_all and pass_old_sel:
        channel_stage = 502
    #elif pass_mu_all:
    #    channel_stage = 501
    #elif pass_mu and pass_old_presel_os:
    #    channel_stage = 3
    #elif pass_mu and pass_old_presel:
    #    channel_stage = 2

    elif pass_el_all and pass_old_sel_os:
        channel_stage = 513
    elif pass_el_all and pass_old_sel:
        channel_stage = 512
    #elif pass_el_all:
    #    channel_stage = 511
    #elif pass_el and pass_old_presel_os:
    #    channel_stage = 13
    #elif pass_el and pass_old_presel:
    #    channel_stage = 12

    return channel_stage


def passes_dy_tautau_selection_stages(passed_triggers, leps, N_jets, taus, proc_met):
    channel_stage = 0
    pass_mu, pass_elmu, pass_elmu_el, pass_mumu, pass_elel, pass_el, pass_mu_all, pass_el_all = passed_triggers

    # muon and OS tau and no b
    pass_dy_objects_mu = pass_mu and len(leps[0]) == 1 and len(taus) > 0 and N_jets[0] == 0
    pass_dy_objects_el = pass_el and len(leps[0]) == 1 and len(taus) > 0 and N_jets[0] == 0
    opposite_sign = (pass_dy_objects_mu or pass_dy_objects_el) and leps[4][0] * taus[0][2] < 0

    pass_dy_mass = False
    if pass_dy_objects_mu or pass_dy_objects_el:
        nom_tau = taus[0][0]*taus[0][1][0]
        pair    = leps[0][0] + nom_tau
        pair_mass = pair.mass()
        pass_dy_mass = 75. < pair_mass < 105.

    if   pass_dy_objects_mu and opposite_sign and pass_dy_mass:
        channel_stage = 103
    elif pass_dy_objects_mu and opposite_sign:
        channel_stage = 102

    elif pass_dy_objects_mu and pass_dy_mass:
        channel_stage = 203
    elif pass_dy_objects_mu:
        channel_stage = 202

    elif pass_dy_objects_el and opposite_sign and pass_dy_mass:
        channel_stage = 113
    elif pass_dy_objects_el and opposite_sign:
        channel_stage = 112

    elif pass_dy_objects_el and pass_dy_mass:
        channel_stage = 213
    elif pass_dy_objects_el:
        channel_stage = 212

    return channel_stage

def passes_dy_mumu_selection_stages(passed_triggers, leps, N_jets, taus, proc_met):
    channel_stage = 0
    pass_mu, pass_elmu, pass_elmu_el, pass_mumu, pass_elel, pass_el, pass_mu_all, pass_el_all = passed_triggers

    # muon and OS tau and no b
    pass_dy_objects_mu = pass_mumu and len(leps[0]) == 2
    pass_dy_objects_el = pass_elel and len(leps[0]) == 2
    no_b_jets = N_jets[0] == 0
    pass_dy_mass = False
    if pass_dy_objects_mu or pass_dy_objects_el:
        pair    = leps[0][0] + leps[0][1]
        pair_mass = pair.mass()
        pass_dy_mass = 75. < pair_mass < 105.

    if   pass_dy_objects_mu and no_b_jets and pass_dy_mass:
        channel_stage = 105
    elif pass_dy_objects_mu and no_b_jets:
        channel_stage = 103
    elif pass_dy_objects_mu:
        channel_stage = 102
    elif pass_dy_objects_el and no_b_jets and pass_dy_mass:
        channel_stage = 115
    elif pass_dy_objects_el and no_b_jets:
        channel_stage = 113
    elif pass_dy_objects_el:
        channel_stage = 112

    return channel_stage

logging = Logging.getLogger("common")

logging.info('importing ROOT')
import ROOT

ROOT.gROOT.Reset()


## maybe roccors are in tt lib now
#print '''loading Rochester Corrections'''
##libpath = "roccor_wrapper_cc"
##libpath = "RoccoR_cc"
#libpath = "RoccoR_cc.so"
###cdll.LoadLibrary( libpath ) # TODO: if error -- break
##roccor_lib = CDLL( libpath )
##roccor_lib.wrapper_kScaleDT.restype         = c_double
##roccor_lib.wrapper_kScaleFromGenMC.restype  = c_double
##roccor_lib.wrapper_kScaleAndSmearMC.restype = c_double
##ROOT.gROOT.ProcessLine(".L RoccoR.cc++")
#ROOT.gSystem.Load(libpath)
## it crashes on exit for some root reason
#print '''Rochester Corrections loaded'''

#ROOT.gSystem.Unload("RoccoR_cc")
#print "test:2", roccors.kScaleDT(1, 40., 0.1, 0.3, 0, 0)

# the lib is needed for BTagCalibrator and Recoil corrections
# TODO: somehow these 2 CMSSW classes should be acceptable from pyROOT on its' own without the whole lib
#ROOT.gSystem.Load("libUserCodettbar-leptons-80X.so")

print '''loaded all libraries'''



##ROOT.gSystem.Load("RoccoR_cc")
## -- trying to add it into ttbar lib
##roccors = ROOT.roccor_wrapper # wrapper does not work somehow
#roccors = ROOT.RoccoR("rcdata.2016.v3")
### works as roccors.kScaleDT(1, 40., 0.1, 0.3, 0, 0)
#print "rochcor test:", roccors.kScaleDT(1, 40., 0.1, 0.3, 0, 0)



from ROOT import TFile, TTree, TH1D, TH2D, TLorentzVector, TVector3, gROOT, gSystem, TCanvas, TGraphAsymmErrors, TMath, TString
from ROOT.Math import LorentzVector

# the names for the positions of renorm refact weights
MUf_nom_MUr_nom    = 0
MUf_up_MUr_nom     = 1
MUf_down_MUr_nom   = 2
MUf_nom_MUr_up     = 3
MUf_up_MUr_up      = 4
MUf_down_MUr_up    = 5
MUf_nom_MUr_down   = 6
MUf_up_MUr_down    = 7
MUf_down_MUr_down  = 8

pileup_ratio_h2 = array('d',
[       0., 5.2880843860098, 1.9057428051882, 0.8279489845062, 1.0183017803649, 0.8886546873859, 0.4586617471559, 0.4516527021066,
 0.6156808676490, 0.1704361822954, 0.3242275168011, 0.5576231467578, 0.6113064296958, 0.6618982312954, 0.7238796620111, 0.7496036493902, 0.7527610638253,
 0.7469665061269, 0.7399984493388, 0.7465103664963, 0.7983735170584, 0.8590931105313, 0.9034039953181, 0.9130886924812, 0.9080208719211, 0.9319379163103,
 0.9829435013036, 1.0322224123557, 1.1140955233618, 1.2058248250755, 1.2824965710179, 1.4313360174361, 1.5303677831147, 1.6810938978350, 1.8448654681967,
 1.9861677547885, 2.1231190233093, 2.2767850875912, 2.3717805455775, 2.4839160504037, 2.5332424087905, 2.4972378057180, 2.5549035019904, 2.5646400146497,
 2.5977594311101, 2.5410596186181, 2.4380019788525, 2.4582931794916, 2.4394333130376, 2.4135104510863, 2.3435484251709, 2.3026709406522, 2.2808207974466,
 2.1519700764652, 2.0425342255892, 1.5440178702333, 0.8925931999325, 0.4286677365298, 0.2057304468250, 0.0757408709610, 0.0356223199915, 0.0133823865731,
 0.0044643816760, 0.0021957442680, 0.0006438366822, 0.0002663753765, 0.0001247116374, 0.0000712745971, 0.0000307860041, 0.0000294100598, 0.0000103247078,
 0.0000055513056, 0.0000033501883, 0.0000035220198, 0.0000009221242])

pileup_ratio_b = array('d',
[       0.,  1.1931990788,  1.1075569611,  1.0661905165,  1.6281619271,  1.5759491716,  1.1654388423,  1.9665149072,  6.4816307673,  4.6824471835,
  3.5057573726,  2.7367422100,  2.0880922033,  1.8309627522,  1.7165427917,  1.6325677435,  1.5679104982,  1.5226132395,  1.4786172115,  1.4187556193,  1.3584199443,
  1.2535044794,  1.1332899913,  0.9984288822,  0.8518945520,  0.7154925510,  0.5827665258,  0.4457781053,  0.3302472894,  0.2310463295,  0.1498402182,  0.0965442096,
  0.0567106928,  0.0327699725,  0.0182363853,  0.0096751297,  0.0050124552,  0.0026185972,  0.0013947431,  0.0008483666,  0.0006315445,  0.0005978998,  0.0007339426,
  0.0010207996,  0.0015749231,  0.0025230767,  0.0042221245,  0.0078659845,  0.0152193676,  0.0308807634,  0.0644886301,  0.1423880343,  0.3304052121,  0.7589089584,
  1.8143503956,  3.5544006370,  5.4397543208,  7.0034023639,  9.0337430924,  8.8675115583, 10.9147825813, 10.4253507840,  8.5070481287,  9.7647211573,  6.3458905766,
  5.5327398483,  5.2282218781,  5.8413138671,  4.8338755929,  8.7575484032,  5.8155542511,  5.9275685169,  6.8122980083, 13.7528845721,  6.9540203430,  ])

pileup_ratio_sum = array('d',
[       0., 2.7084899, 1.2387248, 0.8567540, 1.1319902, 1.0625620, 0.6788884, 0.6404915, 1.5813397, 1.2265823, 1.1225145, 1.1346217, 1.0804378, 1.0719238,
 1.0481945, 1.0210914, 1.0035636, 0.9938051, 0.9824541, 0.9712776, 0.9866336, 0.9979909, 1.0111960, 1.0147509, 1.0078850, 1.0150851, 1.0273812, 1.0179873, 1.0242233,
 1.0233836, 0.9961096, 1.0091713, 0.9729596, 0.9599715, 0.9457656, 0.9168015, 0.8878596, 0.8702436, 0.8374816, 0.8195456, 0.7897366, 0.7430977, 0.7321259, 0.7130387,
 0.7051140, 0.6768290, 0.6399811, 0.6383632, 0.6289537, 0.6205741, 0.6051458, 0.6053456, 0.6291184, 0.6658740, 0.8105650, 0.9697623, 1.1145130, 1.2538530, 1.5301977,
 1.4704201, 1.7955171, 1.7098770, 1.3936616, 1.5989832, 1.0389566, 0.9057558, 0.8558735, 0.9562222, 0.7912980, 1.4335911, 0.9519911, 0.9703263, 1.1151534, 2.2513064,
 1.1383523])



pileup_ratio = array('d', [
0.360609416811339, 0.910848525427002, 1.20629960507795, 0.965997726573782, 1.10708082813183, 1.14843491548622, 0.786526251164482, 0.490577792661333, 0.740680941110478,
0.884048630953726,
0.964813189764159, 1.07045369167689, 1.12497267309738, 1.17367530613108, 1.20239808206413, 1.20815108390021, 1.20049333094509, 1.18284686347315, 1.14408796655615,
1.0962284704313, 1.06549162803223, 1.05151011089581, 1.05159666626121, 1.05064452078328, 1.0491726301522, 1.05772537082991, 1.07279673875566, 1.0837536468865, 1.09536667397119,
1.10934472980173, 1.09375894592864, 1.08263679568271, 1.04289345879947, 0.9851490341672, 0.909983816540809, 0.821346330143864, 0.71704523475871, 0.609800913869359, 0.502935245638477,
0.405579825620816, 0.309696044611377, 0.228191137503131, 0.163380359253309, 0.113368437957202, 0.0772279997453792, 0.0508111733313502, 0.0319007262683943, 0.0200879459309245, 0.0122753366005436,
0.00739933885813127, 0.00437426967257811, 0.00260473545284139, 0.00157047254226743, 0.000969500595715493, 0.000733193118123283, 0.000669817107713128, 0.000728548958604492, 0.000934559691182011, 0.00133719688378802,
0.00186652283903214, 0.00314422244976771, 0.00406954793369611, 0.00467888840511915, 0.00505224284441512, 0.00562827194936864, 0.0055889504870752, 0.00522867039470319, 0.00450752163476433, 0.00395300774604375,
0.00330577167682956, 0.00308353042577215, 0.00277846504893301, 0.00223943190687725, 0.00196650068765464, 0.00184742734258922,])

## recalc from ntuple, basically the same
#pileup_ratio = array('d', [
#0.344887, 0.869219, 1.12271, 0.968119, 1.04695, 1.10773, 0.740941, 0.477038, 0.719406, 0.857546,
#0.938079, 1.04656, 1.09508, 1.14765, 1.17614, 1.18508, 1.17563, 1.16215, 1.12725, 1.08148,
#1.05394, 1.04048, 1.04329, 1.04534, 1.04798, 1.05541, 1.07543, 1.08567, 1.10311, 1.11874,
#1.10373, 1.09798, 1.05651, 1.00496, 0.931962, 0.842275, 0.736296, 0.626743, 0.520142, 0.419807,
#0.323103, 0.239115, 0.170803, 0.119757, 0.0821237, 0.0541623, 0.0338, 0.0215088, 0.0131677, 0.00806879,
#0.00477981, 0.00283188, 0.00173487, 0.00107583, 0.000821263, 0.000757651, 0.000829695, 0.00102146, 0.00152372, 0.00205202,
#0.00372623, 0.00458362, 0.00632652, 0.00643546, 0.00646819, 0.00677855, 0.00618487, 0.00519075, 0.00488287, 0.00433706,
#0.00335919, 0.00359263, 0.00266798, 0.00270629, 0.00261116,])


pileup_ratio_up = array('d', [0.351377216124927, 0.717199649125846, 1.14121536968772, 0.84885826611733, 1.00700929402897, 1.03428595270903, 0.717444379696992, 0.344078389355127, 0.499570875027422,
0.606614916257104, 0.632584599390169, 0.731450949466174, 0.827511723989754, 0.910682115553867, 0.960170981598162, 0.988896170761361, 1.02468865580207, 1.05296667126403, 1.05112033565679,
1.0269129153969, 1.00548641752714, 0.998316130432865, 1.01492587998551, 1.03753749807849, 1.05742218946485, 1.08503978097083, 1.12134132247053, 1.15585339474274, 1.19214399856171,
1.23308400947467, 1.24528804633732, 1.26786364716917, 1.26101551498967, 1.23297806722714, 1.18042533075471, 1.10534683838101, 1.00275591661645, 0.889094305531985, 0.768791254270252,
0.655054015673457, 0.533361034358457, 0.423095146361996, 0.329177839117034, 0.250352385505809, 0.188377378855567, 0.137852651411779, 0.0968577167707531, 0.0686240187247059, 0.0473889635126706,
0.0323695027438475, 0.0216752397914914, 0.0145119352923332, 0.00961177893634792, 0.00615582219138384, 0.00430085627914427, 0.00305735512896403, 0.00223567790438986, 0.00189369737638594, 0.00199585978316291,
0.00236236592656064, 0.00372472999463276, 0.00474687312579969, 0.00549508151576102, 0.00603023110946686, 0.0068545111910253, 0.00695838760530896, 0.00666224781277046, 0.00588243140681038, 0.00528714370892014,
0.00453424615273565, 0.00433985030329723, 0.00401493171035719, 0.00332436608713241, 0.00300063798808221, 0.00289925128977536,])

pileup_ratio_down = array('d', [0.37361294640242, 1.1627791004568, 1.26890787896295, 1.10266790442705, 1.23456697093644, 1.26278991594152, 0.909648777562084, 0.759569490571151, 1.09035651921682,
1.34530547603283, 1.48713160105, 1.52535976889483, 1.49730550773404, 1.49792998045778, 1.49767851097519, 1.44431045398336, 1.3681909492045, 1.29912252494785, 1.2274279217797,
1.16525969099909, 1.12531044676724, 1.09094501417685, 1.06405434433422, 1.03997120824565, 1.0185716022098, 1.00560949501652, 0.997570939806059, 0.985543761409897, 0.972557804582185,
0.957832827239337, 0.9139572640153, 0.872252387173971, 0.808388185417578, 0.733817960498049, 0.650440963845892, 0.561688505024782, 0.466564380334112, 0.374428618658619, 0.28845274688129,
0.214909665968644, 0.149991974352384, 0.100014138338029, 0.0642260884603397, 0.0396553405911344, 0.0238687936736627, 0.0137921542898078, 0.00756854010632403, 0.00415483516246187, 0.00221776872027937,
0.00118249725637452, 0.000641889697310868, 0.000383647166012176, 0.000273637590071334, 0.000242902582071058, 0.000291239677209452, 0.000394091114279828, 0.000542541231466254, 0.000771067920964491, 0.00113596447675107,
0.00158061353194779, 0.00261959852500539, 0.00331800452823827, 0.00372426930370732, 0.00392086545082614, 0.00425479965493548, 0.00411256966391362, 0.00374240422174387, 0.00313603438166934, 0.00267155793176928,
0.00216878198028599, 0.00196249821290853, 0.00171433839159669, 0.00133866519755926, 0.00113810604240254, 0.00103447940224886,])

















pileup_ratio_ele = array('d', [
   0.413231   ,    1.01701    ,    1.19502     ,   0.883906   ,    1.05852    ,    1.11823    ,    0.789439   ,    0.515477   ,    0.81338    ,    0.990148   ,
   1.0919     ,    1.21784    ,    1.28268     ,   1.33936    ,    1.37267    ,    1.38001    ,    1.37224    ,    1.35253    ,    1.30805    ,    1.25303    ,
   1.21761    ,    1.20085    ,    1.1987      ,   1.19257    ,    1.1807     ,    1.17079    ,    1.15238    ,    1.10667    ,    1.03375    ,    0.935086   ,
   0.793376   ,    0.65125    ,    0.502727    ,   0.369298   ,    0.25859    ,    0.173207   ,    0.110361   ,    0.0677957  ,    0.0403186  ,    0.0236369  ,
   0.0133546  ,    0.00746494 ,    0.00417626  ,   0.00233773 ,    0.0013288  ,    0.000757718,    0.000432788,    0.000266239,    0.000177605,    0.000137241,
   0.000125696,    0.000137018,    0.000167806 ,   0.000215108,    0.000313214,    0.000464376,    0.000669855,    0.000981399,    0.00148275 ,    0.00211313 ,
   0.0035872  ,    0.00465614 ,    0.005359    ,   0.00578897 ,    0.00645001 ,    0.00640537 ,    0.00599263 ,    0.00516618 ,    0.00453067 ,    0.00378886 ,
   0.00353415 ,    0.00318451 ,    0.0025667   ,   0.00225388 ,    0.00211741 ,])

pileup_ratio_up_ele = array('d', [
   0.402668    ,   0.803377   ,    1.15963    ,    0.764147   ,    0.966328    ,   0.995159   ,    0.71563    ,    0.354304   ,   0.541943   ,    0.674778   ,
   0.713035    ,   0.830366   ,    0.942616   ,    1.03882    ,    1.09589     ,   1.12909    ,    1.17068    ,    1.20376    ,   1.20191    ,    1.17404    ,
   1.14928     ,   1.14083    ,    1.15906    ,    1.18279    ,    1.20082     ,   1.22277    ,    1.24559    ,    1.25129    ,   1.23607    ,    1.19534    ,
   1.09539     ,   0.978694   ,    0.82546    ,    0.662451   ,    0.50547     ,   0.367764   ,    0.25369    ,    0.168007   ,   0.10706    ,    0.0667404  ,
   0.0397874   ,   0.0233291  ,    0.0136533  ,    0.00799737 ,    0.00476279  ,   0.00284044 ,    0.00167744 ,    0.00103389 ,   0.000648432,    0.000427764,
   0.000303899 ,   0.000247672,    0.000236803,    0.000258026,    0.000345092 ,   0.000494341,    0.000708128,    0.00104444 ,   0.00159927 ,    0.00231779 ,
   0.00400894  ,   0.00530831 ,    0.00623822 ,    0.00688571 ,    0.00784455  ,   0.00797042 ,    0.00763388 ,    0.00674129 ,   0.00605947 ,    0.00519674 ,
   0.00497399  ,   0.00460162 ,    0.00381017 ,    0.00343914 ,    0.00332295  ,])

pileup_ratio_down_ele = array('d', [
   0.428107   ,   1.29388    ,    1.22078   ,    1.02596    ,    1.1769     ,    1.24377    ,    0.921862   ,    0.814769   ,    1.20901   ,    1.51527   ,
   1.68838    ,   1.73792    ,    1.70826   ,    1.70984    ,    1.71038    ,    1.65067    ,    1.56442    ,    1.48535    ,    1.40303   ,    1.33164   ,
   1.28514    ,   1.24342    ,    1.20714   ,    1.16839    ,    1.12262    ,    1.06993    ,    0.999693   ,    0.900043   ,    0.778486  ,    0.644942  ,
   0.497564   ,   0.37052    ,    0.259917  ,    0.174109   ,    0.111585   ,    0.0687061  ,    0.0404941  ,    0.0232033  ,    0.0129877 ,    0.00721863,
   0.00388086 ,   0.00206418 ,    0.00109735,    0.000585006,    0.000321242,    0.000184087,    0.000114353,    8.57172e-5 ,    7.68135e-5,    8.09537e-5,
   9.42381e-5 ,   0.000118387,    0.0001548 ,    0.000202628,    0.000294603,    0.000431519,    0.00061181 ,    0.000878668,    0.00129927,    0.0018102 ,
   0.00300153 ,   0.00380243 ,    0.0042683 ,    0.00449375 ,    0.00487654 ,    0.00471355 ,    0.0042893  ,    0.00359432 ,    0.00306198,    0.00248572,
   0.0022493  ,   0.00196487 ,    0.0015343 ,    0.00130443 ,    0.00118566 ,])


pu_per_epoch_BCDEF = array('d', [
 0.0095886208, 0.5776865397, 1.0233004203, 1.0083366096, 1.2158510340, 1.2887481725, 0.9705138670, 0.6573177335, 1.1635728135, 1.4857636927,
 1.5009482334, 1.4514586332, 1.4035598214, 1.4171230048, 1.4353069202, 1.4278634811, 1.4044837262, 1.3768186333, 1.3330031729, 1.2744359409, 1.2163075523,
 1.1616073030, 1.1240002373, 1.0941364455, 1.0660695295, 1.0416984593, 1.0130564673, 0.9700712357, 0.9203311560, 0.8692631761, 0.7965866756, 0.7314671223,
 0.6523463831, 0.5687054084, 0.4823132277, 0.3965994070, 0.3121499417, 0.2364186594, 0.1714799246, 0.1202111023, 0.0790221128, 0.0497517154, 0.0302802556,
 0.0178087370, 0.0102758135, 0.0057386316, 0.0030774403, 0.0016794149, 0.0009188412, 0.0005328167, 0.0003492515, 0.0002863425, 0.0002956655, 0.0003516098,
 0.0004988223, 0.0007342755, 0.0010578854, 0.0015503122, 0.0023433569, 0.0033407839, 0.0056724855, 0.0073637764, 0.0084760046, 0.0091564392, 0.0102022625,
 0.0101317755, 0.0094789852, 0.0081717643, 0.0071665366, 0.0059931677, 0.0055902693, 0.0050372087, 0.0040599721, 0.0035651639, 0.0033492907,])

pu_per_epoch_BCDEF_Down = array('d', [
 0.0153667969, 0.6891809486, 1.1601810247, 1.1655664961, 1.3539345418, 1.4539256296, 1.1448982617, 1.0782047058, 1.7930463380, 2.2164832902,
 2.1506625765, 1.9592786349, 1.8240756356, 1.7938228325, 1.7756989722, 1.6936020451, 1.5933227876, 1.5131795726, 1.4280715808, 1.3322524118, 1.2436532283,
 1.1648594314, 1.1058876595, 1.0528664466, 0.9963592864, 0.9389492536, 0.8781488618, 0.8099296851, 0.7417754407, 0.6760367200, 0.5957590110, 0.5236926486,
 0.4450245773, 0.3677395027, 0.2936508402, 0.2255192416, 0.1643016798, 0.1141439873, 0.0752804710, 0.0476059133, 0.0280347470, 0.0157204810, 0.0084836333,
 0.0044127581, 0.0022545087, 0.0011258103, 0.0005567695, 0.0003027543, 0.0001920254, 0.0001562929, 0.0001593515, 0.0001904067, 0.0002453016, 0.0003200674,
 0.0004653120, 0.0006818820, 0.0009671683, 0.0013893910, 0.0020547781, 0.0028630657, 0.0047475144, 0.0060144576, 0.0067514434, 0.0071080944, 0.0077136086,
 0.0074558179, 0.0067847592, 0.0056854570, 0.0048433925, 0.0039318873, 0.0035579083, 0.0031080075, 0.0024269317, 0.0020633280, 0.0018754582,])

pu_per_epoch_BCDEF_Up = array('d', [
 0.0059598412, 0.4901115533, 0.8955066803, 0.8795211517, 1.1062394600, 1.1404684587, 0.8650218669, 0.4437416134, 0.7393319145, 1.0097269213,
 1.0376404815, 1.0609548074, 1.0721267343, 1.1161700434, 1.1533323089, 1.1770118678, 1.2074397521, 1.2296721371, 1.2232723959, 1.1967643924, 1.1679136238,
 1.1379537690, 1.1207650375, 1.1100257565, 1.1032100074, 1.1060674610, 1.1108500776, 1.1022045800, 1.0827578262, 1.0565626732, 0.9997407140, 0.9499532715,
 0.8799800763, 0.8000519083, 0.7105812384, 0.6150169143, 0.5127635583, 0.4144487179, 0.3233486501, 0.2457555729, 0.1764673760, 0.1221868808, 0.0822719095,
 0.0537970753, 0.0346456708, 0.0216435771, 0.0129743697, 0.0078566879, 0.0046618186, 0.0027687511, 0.0016541260, 0.0010443045, 0.0007284114, 0.0005887501,
 0.0006412291, 0.0008277866, 0.0011377166, 0.0016555899, 0.0025265285, 0.0036605553, 0.0063339072, 0.0083903011, 0.0098630619, 0.0108888325, 0.0124065091,
 0.0126064207, 0.0120745864, 0.0106630150, 0.0095846516, 0.0082200638, 0.0078677807, 0.0072787942, 0.0060268729, 0.0054399860, 0.0052561842,])



pu_per_epoch_GH = array('d', [
 0.7923966165550341, 1.3206676450689692, 1.4314050681219410, 0.9139173915540429, 0.9732839512546894, 0.9758371835396235, 0.5602049252864559,
 0.2854729322811546, 0.2204862425770137, 0.1438852294551190, 0.3053192173605273, 0.6017836345742870, 0.7822855124164086, 0.8742128227661290, 0.9159002630545462,
 0.9378846334114083, 0.9495663163686077, 0.9442440878690979, 0.9117058293335381, 0.8770170417648185, 0.8799748963086315, 0.9160809496039715, 0.9625337586731879,
 0.9971456032972257, 1.0283877191853710, 1.0774400368137698, 1.1462819788479417, 1.2235925788952713, 1.3106768396716892, 1.4046659432469826, 1.4593082605398517,
 1.5146069014110246, 1.5233012690983123, 1.4974120340854682, 1.4360570232383967, 1.3438230304232581, 1.2151026963636524, 1.0690945958524307, 0.9106548859389182,
 0.7566089147801717, 0.5934459067043246, 0.4476876554380586, 0.3271054251238090, 0.2309154193601109, 0.1595851549791273, 0.1062544378149940, 0.0673559532806101,
 0.0427320898021242, 0.0262448490604951, 0.0158457770655627, 0.0093254025228037, 0.0054565650542544, 0.0031385995670769, 0.0017295617570310, 0.0010214898388716,
 0.0005905275011105, 0.0003234355927401, 0.0001771290570776, 0.0000995297289998, 0.0000530499437056, 0.0000342332435665, 0.0000173515193266, 0.0000080957416536,
 0.0000037146197696, 0.0000018526016169, 0.0000008674969736, 0.0000004002877156, 0.0000001760093781, 0.0000000804321989, 0.0000000354304911, 0.0000000174636467,
 0.0000000083012601, 0.0000000035114725, 0.0000000016057138, 0.0000000007810826,])

pu_per_epoch_GH_Down = array('d', [
 0.81428793844913028, 1.74534663553406610, 1.40265194859107756, 1.02529687327769592, 1.08773425686398384, 1.02767620586297292,
 0.62027093556673207, 0.36761988051584205, 0.22598535457795466, 0.27367834795556040, 0.67092979664080921, 0.99160106467545983, 1.09534930584835966,
 1.13395571674522722, 1.15568772724638946, 1.13765907121676846, 1.09125877593570952, 1.03581290492279665, 0.98061787954372870, 0.95984386651741038,
 0.97973858379080214, 1.00002335463852643, 1.01259558055105270, 1.02410890369605023, 1.04589437240884764, 1.08760823580738997, 1.14447117354317163,
 1.20156537301837574, 1.25644020730890071, 1.30446669755839539, 1.30536984313194937, 1.30101217348902454, 1.25535815482144364, 1.18412721357627415,
 1.08932481083630783, 0.97520676871619461, 0.83837446081404476, 0.69460195322691609, 0.55067398626496733, 0.42070831195011738, 0.30001033241545638,
 0.20370291959251538, 0.13279431867943567, 0.08300686750850418, 0.05045630187280000, 0.02937289581248214, 0.01619365100892409, 0.00889324082405334,
 0.00470961560362523, 0.00244481996169108, 0.00123545510112472, 0.00062135021495797, 0.00030849336596121, 0.00014798303586979, 0.00007711517872165,
 0.00004008260745848, 0.00002021197012521, 0.00001047534980043, 0.00000574164742812, 0.00000308308610710, 0.00000207130589446, 0.00000112578232819,
 0.00000057412526796, 0.00000028848668407, 0.00000015491212658, 0.00000007575632291, 0.00000003524537034, 0.00000001513400039, 0.00000000658427975,
 0.00000000271431289, 0.00000000123416400, 0.00000000053813129, 0.00000000020937973, 0.00000000008474601, 0.00000000003929975,])

pu_per_epoch_GH_Up = array('d', [
 0.77627151099057, 0.99653850182123, 1.44345865370983, 0.81114020684281, 0.88494710058229, 0.90367194860043, 0.53591083091285, 0.22148369164314,
 0.20464354023011, 0.11075115335543, 0.13432936165235, 0.32613193593120, 0.52661329157577, 0.65791369521069, 0.72256523581772, 0.75749739033748, 0.79988803112911,
 0.83560336934485, 0.83935773399918, 0.81798066008565, 0.80568660632538, 0.82654961208942, 0.88473425193572, 0.94837055783462, 1.00109871057968, 1.05917334508425,
 1.13424597667433, 1.22184675399576, 1.32669915771627, 1.45022086535295, 1.54733266600517, 1.65892159206283, 1.72972336154445, 1.76551607252308, 1.75837588885673,
 1.70849655365464, 1.60549032181041, 1.47295093869894, 1.31672565975546, 1.15852763206484, 0.97237258182634, 0.79323910186310, 0.63289428451977, 0.49213307046849,
 0.37748120225896, 0.28080020413459, 0.20004177233696, 0.14337334750013, 0.09994717770056, 0.06878108485635, 0.04630300598811, 0.03107832921615, 0.02053911732033,
 0.01300381784358, 0.00880252563489, 0.00579992419207, 0.00358626995251, 0.00218659280931, 0.00134309001403, 0.00076547631898, 0.00051520781680, 0.00026513104410,
 0.00012207373190, 0.00005371714956, 0.00002505645467, 0.00001080316276, 0.00000458637390, 0.00000188383115, 0.00000082589010, 0.00000035999690, 0.00000018068523,
 0.00000008947486, 0.00000004009622, 0.00000001966999, 0.00000001031619,])


def calc_pu_per_runs(pu_ele, pu_gld, ele_B = 19128.209 / (19128.209 + 12210.15), ele_H = 12210.15  / (19128.209 + 12210.15), gld_B = 19713.888 / (19713.888 + 16146.177), gld_H = 16146.177 / (19713.888 + 16146.177)):
    pu_H = (pu_ele    / ele_B - pu_gld / gld_B) / (ele_H / ele_B - gld_H / gld_B)
    pu_B = (pu_ele    / ele_H - pu_gld / gld_H) / (ele_B / ele_H - gld_B / gld_H)
    return pu_B, pu_H













from module_leptons import lepton_muon_SF, lepton_muon_trigger_SF, lepton_electron_SF, lepton_electron_trigger_SF, dilepton_or_sfs, dilepton_and_sfs

if with_bSF:
    import support_btagging_sf
    from support_btagging_sf import calc_btag_sf_weight, bEff_histo_b, bEff_histo_c, bEff_histo_udsg
    from support_btagging_sf import h_control_btag_eff_b, h_control_btag_eff_c, h_control_btag_eff_udsg, h_control_btag_weight_b, h_control_btag_weight_c, h_control_btag_weight_udsg, h_control_btag_weight_notag_b, h_control_btag_weight_notag_c, h_control_btag_weight_notag_udsg

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

def transverse_cos(v1, v2):
    t_cos = (v1.Px()*v2.Px() + v1.Py()*v2.Py()) / TMath.Sqrt((v1.Px()*v1.Px() + v1.Py()*v1.Py())* (v2.Px()*v2.Px() + v2.Py()*v2.Py()))
    return t_cos

m_w_boson = 80.385
m_t_quark = 172.5

def calc_lj_var(ev, light_jets, b_jets, save_all_permutations=False, isMC=False):
    closest_pair_gens = (0, 0)
    closest_b_gen = 0
    if len(b_jets) == 0 or len(light_jets) < 2: return 5000., 5000., 5000., (closest_pair_gens, closest_b_gen), []
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

    # save [(mW, mt)] of all permutations if requested
    all_masses = []

    # pairs of light jets
    light_jet_pairs = []
    for i in range(len(light_jets)):
      for u in range(i):
        ji, multi, _, _, _, _, _ = light_jets[i]
        ju, multu, _, _, _, _, _ = light_jets[u]
        pair = ji * multi[0] + ju * multu[0]
        if save_all_permutations:
            light_jet_pairs.append(pair)
        new_dist = abs(pair.mass() - m_w_boson)
        if new_dist < dist_W:
            dist_W = new_dist
            closest_to_W = pair
            if isMC:
                closest_pair_gens = (ev.jet_matching_gen[light_jets[i][6]], ev.jet_matching_gen[light_jets[u][6]])

    # closest to 173
    dist_t = 99999.
    for j, mult, _, _, _, _, jet_index in b_jets:
        b_jet = j * mult[0]
        pair = b_jet + closest_to_W
        if save_all_permutations:
            for cand_W in light_jet_pairs:
                cand_t = cand_W + b_jet
                all_masses.append((cand_W.mass(), cand_t.mass()))
        new_dist = abs(pair.mass() - m_t_quark)
        if new_dist < dist_t:
            dist_t = new_dist
            closest_to_t = pair
            if isMC:
                closest_b_gen = ev.jet_matching_gen[jet_index]

    return TMath.Sqrt(dist_W*dist_W + dist_t*dist_t), closest_to_W.mass(), closest_to_t.mass(), (closest_pair_gens, closest_b_gen), all_masses

def solve_quadratic(a, b, c):

    # calculate the discriminant
    d = (b**2) - (4*a*c)
    if d < 0:
        return False, (0, 0)

    # find two solutions
    sol1 = (-b-math.sqrt(d))/(2*a)
    sol2 = (-b+math.sqrt(d))/(2*a)

    return True, (sol1.real, sol2.real)

def solve_met_z_coef(met_pt, l_pt, l_z, l_e):
    # met*lep = (m_w^2 - m_l^2)/2
    # z = k*e + b

    b = (- met_pt * l_pt - m_w_boson**2 / 2) / l_z
    k = (l_e) / l_z

    return k, b

def solve_met_e_z(met_pt, l_pt, l_z, l_e):

    k, b = solve_met_z_coef(met_pt, l_pt, l_z, l_e)

    solved, (e1, e2) = solve_quadratic((k**2 - 1), 2*k*b, (b**2 + met_pt**2))
    if not solved:
        print "met E could not be found"

    if (e1 < 0) and (e2 < 0):
        print "both met E < 0"
    #elif (e1 > 0) and (e2 > 0):
    #    print "both met E > 0"

    e = max(e1, e2)

    z = k*e + b
    return e, z


#def PFTau_FlightLength_significance(TVector3 pv,TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov ){
def PFTau_FlightLength_significance(pv,  PVcov, sv, SVcov):
   run_test = False # option to print contents of some of the structures during calculation

   SVPV = sv - pv
   FD = ROOT.TVectorF()
   FD.ResizeTo(3);
   #FD(0) = SVPV.X();
   #FD(1) = SVPV.Y();
   #FD(2) = SVPV.Z();
   FD.SetElements(array('f', [SVPV.X(), SVPV.Y(), SVPV.Z()]))

   # input covs are
   # ROOT.Math.SMatrix(float, 3, 3, ROOT.Math.MatRepSym(float, 3) )

   #TMatrixT<double> PVcv;
   PVcv = ROOT.TMatrixT(float)()
   PVcv.ResizeTo(3,3);
   #TMatrixT<double> SVcv;
   SVcv = ROOT.TMatrixT(float)()
   SVcv.ResizeTo(3,3);
   #for(int nr =0; nr<PVcov.GetNrows(); nr++){
   #  for(int nc =0; nc<PVcov.GetNcols(); nc++){
   #    PVcv(nr,nc) = PVcov(nr,nc);
   #  }
   #}
   #for(int nr =0; nr<SVcov.GetNrows(); nr++){
   #  for(int nc =0; nc<SVcov.GetNcols(); nc++){
   #    SVcv(nr,nc) = SVcov(nr,nc);
   #  }
   #}
   # assume rows -- first index, coloumns -- second
   for nr in range(3):
      for nc in range(3):
         PVcv[nr][nc] = PVcov(nr, nc)
         SVcv[nr][nc] = SVcov(nr, nc)

   #TMatrixT<double> SVPVMatrix(3,1);
   #SVPVMatrix = ROOT.TMatrixT(float)(3,1) # Error in <operator*=(const TMatrixT &)>: source matrix has wrong shape
   SVPVMatrix = ROOT.TMatrixT(float)(3,3)
   #for(int i=0; i<SVPVMatrix.GetNrows();i++){
   #  SVPVMatrix(i,0)=FD(i);
   #}
   for i in range(SVPVMatrix.GetNrows()):
       SVPVMatrix[i][0] = FD(i)
       SVPVMatrix[i][1] = 0
       SVPVMatrix[i][2] = 0

   #TMatrixT<double> SVPVMatrixT=SVPVMatrix;
   SVPVMatrixT = SVPVMatrix.Clone()
   SVPVMatrixT.T()

   if run_test:
      SVcv.Print()
      PVcv.Print()
   SVcv += PVcv
   PVSVcv = SVcv
   if run_test:
      SVcv.Print()
      PVSVcv.Print()

   if run_test:
       SVPVMatrixT.Print()
       PVSVcv     .Print()
       SVPVMatrix .Print()

   #TMatrixT<double> lambda2 = SVPVMatrixT*(SVcv + PVcv)*SVPVMatrix;
   #lambda2 = SVPVMatrixT*(SVcv + PVcv)*SVPVMatrix
   #lambda2 = SVPVMatrixT*PVSVcv*SVPVMatrix
   # doc says "Compute target = target * source inplace"
   # https://root.cern.ch/doc/v608/classTMatrixT.html
   lambda2 = SVPVMatrixT
   lambda2 *= PVSVcv
   lambda2 *= SVPVMatrix
   if run_test:
       lambda2.Print()

   sigmaabs = TMath.Sqrt(lambda2(0,0))/SVPV.Mag()
   sign = SVPV.Mag()/sigmaabs

   return SVPV.Mag(), sigmaabs, sign


control_counters = TH1D("control_counters", "", 500, 0, 500)

# no data types/protocols in ROOT -- looping has to be done manually
def full_loop(tree, ttree_out, dtag, lumi_bcdef, lumi_gh, logger, channels_to_select):
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

    # fix the naming in dtag-dset info files
    tt_systematic_datasets = {'fsrup': 'FSRUp', 'fsrdown': 'FSRDown',
        'TuneCUETP8M2T4down': 'TuneCUETP8M2T4Down', 'TuneCUETP8M2T4up': 'TuneCUETP8M2T4Up',
        'isrup': 'ISRUp', 'isrdown': 'ISRDown',
        'hdampUP': 'HDAMPUp', 'hdampDOWN': 'HDAMPDown',
        'GluonMoveCRTune': 'GluonMoveCRTuneUp', 'QCDbasedCRTune': 'QCDbasedCRTuneUp'}

    def which_sys(dtag, systematics=tt_systematic_datasets):
        # if the dtag name is already fixed return it as is
        for sys_name in systematics.values():
            if sys_name in dtag:
                return sys_name
        # else try to translate the old name
        for sys_name in systematics.keys():
            if sys_name in dtag:
                return systematics[sys_name]
        return None

    isTT_systematic = isTT and which_sys(dtag)

    isSTop = 'SingleT' in dtag or 'tchannel' in dtag or 'schannel' in dtag
    isSTopTSchannels = 'tchannel' in dtag or 'schannel' in dtag
    isDY    = isMC and ('DY' in dtag)
    isWJets = isMC and ('WJet' in dtag or 'W1Jet' in dtag or 'W2Jet' in dtag or 'W3Jet' in dtag or 'W4Jet' in dtag)
    isWJetsInclusive = 'WJet' in dtag
    isQCD = 'QCD' in dtag
    isDibosons = 'WW' in dtag or 'ZZ' in dtag or 'WZ' in dtag

    logger.write(' '.join(name + '=' + str(setting) for name, setting in
        [('isMC', isMC), ('save_weights', save_weights), ('aMCatNLO', aMCatNLO), ('isTT', isTT), ('isTT_systematic', isTT_systematic), ('isSTop', isSTop), ('isSTopTSchannels', isSTopTSchannels), ('isDY', isDY), ('isWJets', isWJets), ('isQCD', isQCD), ('isDibosons', isDibosons)]) + '\n')


    logging.info("load root_funcs")
    ROOT.gROOT.ProcessLine(".L /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/root_funcs.C+") # TODO: change absolute path to relative

    # Recoil corrections
    doRecoilCorrections = WITH_RECOIL_CORRECTIONS and (isWJets or isDY) # TODO: check if it is needed for DY
    if doRecoilCorrections:
        logging.info("will use Recoil Corrections")
        #ROOT.gROOT.ProcessLine(".L /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/lib/slc6_amd64_gcc530/libHTT-utilitiesRecoilCorrections.so")
        ROOT.gROOT.ProcessLine(".L /lstore/cms/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/libHTT-utilitiesRecoilCorrections.so")

        #ROOT.gROOT.ProcessLine(".L /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/recoil_corrections.C+") # TODO: change absolute path to relative
        ROOT.gROOT.ProcessLine(".L /lstore/cms/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/recoil_corrections.C+")
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


    control_hs = OrderedDict([
    ('weight_pu',     TH1D("weight_pu", "", 50, 0, 2)),
    ('weight_pu_el', TH1D("weight_pu_el", "", 50, 0, 2)),
    ('weight_pu_mu', TH1D("weight_pu_mu", "", 50, 0, 2)),

    ('weight_pu_mu_up',     TH1D("weight_pu_up", "", 50, 0, 2)),
    ('weight_pu_mu_dn',     TH1D("weight_pu_dn", "", 50, 0, 2)),
    ('weight_pu_el_up', TH1D("weight_pu_ele_up", "", 50, 0, 2)),
    ('weight_pu_el_dn', TH1D("weight_pu_ele_dn", "", 50, 0, 2)),

    ('weight_pu_sum', TH1D("weight_pu_sum", "", 50, 0, 2)),
    ('weight_pu_b',   TH1D("weight_pu_b", "", 50, 0, 2)),
    ('weight_pu_h2',  TH1D("weight_pu_h2", "", 50, 0, 2)),
    ('weight_top_pt', TH1D("weight_top_pt", "", 50, 0, 2)),

    #('roccor_factor', TH1D("roccor_factor", "", 50, 0, 2)),

    ('weights_gen_weight_too',         TH1D("weights_gen_weight_too",         "", 50, 0, 2)),
    ('weights_gen_weight_norm',        TH1D("weights_gen_weight_norm",        "", 50, 0, 2)),
    ('weights_gen_weight_average',     TH1D("weights_gen_weight_average",     "", 50, 0, 2)),
    ('weights_gen_weight_alphasUp',    TH1D("weights_gen_weight_alphasUp",    "", 50, 0, 2)),
    ('weights_gen_weight_alphasDown',  TH1D("weights_gen_weight_alphasDown",  "", 50, 0, 2)),
    ('weights_gen_weight_centralFrag', TH1D("weights_gen_weight_centralFrag", "", 50, 0, 2)),
    ('weights_gen_weight_Peterson',    TH1D("weights_gen_weight_Peterson",    "", 50, 0, 2)),
    ('weights_gen_weight_FragUp',      TH1D("weights_FragUp",   "", 50, 0, 2)),
    ('weights_gen_weight_FragDown',    TH1D("weights_FragDown", "", 50, 0, 2)),
    ('weights_gen_weight_semilepbrUp',      TH1D("weights_semilepbrUp",   "", 50, 0, 2)),
    ('weights_gen_weight_semilepbrDown',    TH1D("weights_semilepbrDown", "", 50, 0, 2)),

    ('AlphaSUp',       TH1D("weight_AlphaSUp",       "", 50, 0, 2)),
    ('AlphaSDown',     TH1D("weight_AlphaSDown",     "", 50, 0, 2)),
    ('FragUp',         TH1D("weight_FragUp",         "", 50, 0, 2)),
    ('FragDown',       TH1D("weight_FragDown",       "", 50, 0, 2)),
    ('SemilepBRUp',    TH1D("weight_SemilepBRUp",         "", 50, 0, 2)),
    ('SemilepBRDown',  TH1D("weight_SemilepBRDown",       "", 50, 0, 2)),
    ('PetersonUp',     TH1D("weight_PetersonUp",         "", 50, 0, 2)),
    ('PetersonDown',   TH1D("weight_PetersonDown",       "", 50, 0, 2)),
    ('MrUp',      TH1D("weight_MrUp",         "", 50, 0, 2)),
    ('MrDown',    TH1D("weight_MrDown",       "", 50, 0, 2)),
    ('MfUp',      TH1D("weight_MfUp",         "", 50, 0, 2)),
    ('MfDown',    TH1D("weight_MfDown",       "", 50, 0, 2)),
    ('MfrUp',     TH1D("weight_MfrUp",        "", 50, 0, 2)),
    ('MfrDown',   TH1D("weight_MfrDown",      "", 50, 0, 2)),

    ('weights_gen_weight_nom'   , TH1D('weights_gen_weight_nom'  , "", 50, 0, 2)),
    ('weights_gen_weight_f_rUp' , TH1D('weights_gen_weight_f_rUp', "", 50, 0, 2)),
    ('weights_gen_weight_f_rDn' , TH1D('weights_gen_weight_f_rDn', "", 50, 0, 2)),
    ('weights_gen_weight_fUp_r' , TH1D('weights_gen_weight_fUp_r', "", 50, 0, 2)),
    ('weights_gen_weight_fDn_r' , TH1D('weights_gen_weight_fDn_r', "", 50, 0, 2)),
    ('weights_gen_weight_frUp'  , TH1D('weights_gen_weight_frUp' , "", 50, 0, 2)),
    ('weights_gen_weight_frDn'  , TH1D('weights_gen_weight_frDn' , "", 50, 0, 2)),

    ('weight_z_mass_pt', TH1D("weight_z_mass_pt", "", 50, 0, 2)),
    ('weight_bSF',       TH1D("weight_bSF", "", 50, 0, 2)),

    ('weight_mu_trk_bcdef', TH1D("weight_mu_trk_bcdef", "", 50, 0, 2)),
    ('weight_mu_trk_bcdef_vtx', TH1D("weight_mu_trk_bcdef_vtx", "", 50, 0, 2)),
    ('weight_mu_trk_bcdef_vtx_gen', TH1D("weight_mu_trk_bcdef_vtx_gen", "", 50, 0, 2)),
    ('weight_mu_id_bcdef' , TH1D("weight_mu_id_bcdef", "", 50, 0, 2)),
    ('weight_mu_iso_bcdef', TH1D("weight_mu_iso_bcdef", "", 50, 0, 2)),
    ('weight_mu_trg_bcdef', TH1D("weight_mu_trg_bcdef", "", 50, 0, 2)),
    ('weight_mu_all_bcdef', TH1D("weight_mu_all_bcdef", "", 50, 0, 2)),

    ('weight_mu_trk_gh', TH1D("weight_mu_trk_gh", "", 50, 0, 2)),
    ('weight_mu_trk_gh_vtx', TH1D("weight_mu_trk_gh_vtx", "", 50, 0, 2)),
    ('weight_mu_trk_gh_vtx_gen', TH1D("weight_mu_trk_gh_vtx_gen", "", 50, 0, 2)),
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

    if isTT:
        control_hs.update(dict(('PDFCT14n%dUp' % i, TH1D("weight_PDFCT14n%dUp" % i, "", 50, 0, 2)) for i in range(58) ))

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

    tt_el_procs = ['tt_eltau', 'tt_lj', 'tt_taultauh', 'tt_taulj', 'tt_other']
    tt_mu_procs = ['tt_mutau', 'tt_lj', 'tt_taultauh', 'tt_taulj', 'tt_other']

    systematic_names_nominal = ['NOMINAL']

    all_systematic_objects = ('NOMINAL', 'JESUp', 'JESDown', 'JERUp', 'JERDown', 'TauESUp', 'TauESDown', 'bSFUp', 'bSFDown')

    if isMC:
        systematic_names_all = ['NOMINAL',
                'LEPUp'     ,
                'LEPDown'   ,
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

    systematic_names_all_with_th = systematic_names_all[:]
    systematic_names_th_renorm_refact = ['NOMINAL']
    if isTT:
        systematic_names_all.extend(('TOPPTUp', 'TOPPTDown'))
        systematic_names_toppt = ['NOMINAL', 'TOPPTUp']
        systematic_names_pu_toppt.append('TOPPTUp')

        systematic_names_all_with_th.extend(('TOPPTUp', 'TOPPTDown'))
        systematic_names_all_with_th.extend(('AlphaSUp', 'AlphaSDown', 'FragUp', 'FragDown', 'SemilepBRUp', 'SemilepBRDown', 'PetersonUp', 'PetersonDown'))
        systematic_names_all_with_th.extend(('MrUp', 'MrDown', 'MfUp', 'MfDown', 'MfrUp', 'MfrDown'))
        systematic_names_all_with_th.extend(('PDF_TRIGGER',))

        systematic_names_th_renorm_refact.extend(['MrUp', 'MrDown', 'MfUp', 'MfDown', 'MfrUp', 'MfrDown'])

    if isMC:
        if isTT:
            # TODO: probably can clarify this
            if isTT_systematic:
                systematic_names_toppt    = [isTT_systematic]
                systematic_names_pu_toppt = [isTT_systematic]
                systematic_names_pu       = [isTT_systematic]
                systematic_names_all      = [isTT_systematic]
                systematic_names_nominal  = [isTT_systematic]
                systematic_names_all_with_th = [isTT_systematic]
                systematic_names_th_renorm_refact = [isTT_systematic]

            procs_el     = tt_procs_el     =  (['tt_eltau', 'tt_lj', 'tt_taultauh', 'tt_taulj', 'tt_other'], 'tt_other')
            procs_mu     = tt_procs_mu     =  (['tt_mutau', 'tt_lj', 'tt_taultauh', 'tt_taulj', 'tt_other'], 'tt_other')
            procs_el_3ch = tt_procs_el_3ch =  (['tt_eltau3ch', 'tt_eltau', 'tt_lj', 'tt_taultauh', 'tt_taulj', 'tt_other'], 'tt_other')
            procs_mu_3ch = tt_procs_mu_3ch =  (['tt_mutau3ch', 'tt_mutau', 'tt_lj', 'tt_taultauh', 'tt_taulj', 'tt_other'], 'tt_other')
            procs_el_3ch_fbw = tt_procs_el_3ch_fbw =  (['tt_eltau3ch', 'tt_eltau', 'tt_ljb', 'tt_ljw', 'tt_ljo', 'tt_lj', 'tt_taultauh', 'tt_taulj', 'tt_other'], 'tt_other')
            procs_mu_3ch_fbw = tt_procs_mu_3ch_fbw =  (['tt_mutau3ch', 'tt_mutau', 'tt_ljb', 'tt_ljw', 'tt_ljo', 'tt_lj', 'tt_taultauh', 'tt_taulj', 'tt_other'], 'tt_other')
            #procs_elmu   = tt_procs_elmu   =  (['tt_elmu', 'tt_taueltaumu', 'tt_other'], 'tt_other')
            procs_elmu   = tt_procs_elmu   =  (['tt_elmu', 'tt_ltaul', 'tt_taueltaumu', 'tt_other'], 'tt_other')
            usual_process   = 'tt_other'

        if isWJets:
            wjets_procs = (['wjets', 'wjets_tauh', 'wjets_taul'],  'wjets')
            procs_el_3ch = procs_el = procs_mu_3ch = procs_mu = procs_elmu = procs_mu_3ch_fbw = procs_el_3ch_fbw = wjets_procs
            usual_process = 'wjets'

        if isDY:
            dy_procs = (['dy_tautau', 'dy_other'], 'dy_other')
            procs_el_3ch = procs_el = procs_mu_3ch = procs_mu = procs_elmu = procs_mu_3ch_fbw = procs_el_3ch_fbw = dy_procs
            usual_process = 'dy_other'

        if isSTop:
            s_top_procs_el = (['s_top_eltau', 's_top_lj', 's_top_other'], 's_top_other')
            s_top_procs_mu = (['s_top_mutau', 's_top_lj', 's_top_other'], 's_top_other')
            s_top_procs_elmu = (['s_top_elmu', 's_top_other'], 's_top_other')
            procs_el_3ch = procs_el = procs_el_3ch_fbw = s_top_procs_el
            procs_mu_3ch = procs_mu = procs_mu_3ch_fbw = s_top_procs_mu
            procs_elmu = s_top_procs_elmu
            usual_process = 's_top_other'

        if isQCD:
            qcd_procs = (['qcd'], 'qcd')
            procs_el_3ch = procs_el = procs_mu_3ch = procs_mu = procs_elmu = procs_mu_3ch_fbw = procs_el_3ch_fbw = qcd_procs
            usual_process = 'qcd'

        if isDibosons:
            dibosons_procs = (['dibosons'], 'dibosons')
            procs_el_3ch = procs_el = procs_mu_3ch = procs_mu = procs_elmu = procs_mu_3ch_fbw = procs_el_3ch_fbw = dibosons_procs
            usual_process = 'dibosons'

    else:
        data_procs = (['data'], 'data')
        procs_el_3ch = procs_el = procs_mu_3ch = procs_mu = procs_elmu = procs_mu_3ch_fbw = procs_el_3ch_fbw = data_procs
        usual_process = 'data'

    channels_usual = {'el_presel':      (procs_el_3ch, systematic_names_pu_toppt), #systematic_names_all),
                'el_sel':         (procs_el_3ch, systematic_names_pu_toppt), #systematic_names_all),
                'el_lj':          (procs_el_3ch, systematic_names_all), #systematic_names_all),
                'el_lj_out':      (procs_el_3ch, systematic_names_all), #systematic_names_all),
                #'mu_prepresel':   (procs_mu_3ch, systematic_names_pu_toppt), #systematic_names_all),
                'mu_presel':      (procs_mu_3ch, systematic_names_pu_toppt), #systematic_names_all),
                'mu_sel':         (procs_mu_3ch, systematic_names_pu_toppt), #systematic_names_all),
                'mu_lj':          (procs_mu_3ch, systematic_names_all),
                'mu_lj_out':      (procs_mu_3ch, systematic_names_all),
                'mu_lj_ss':       (procs_mu_3ch, systematic_names_pu_toppt),
                'mu_lj_out_ss':   (procs_mu_3ch, systematic_names_pu_toppt),
                # same sign for some QCD control
                'el_sel_ss':      (procs_el, systematic_names_nominal),
                'mu_sel_ss':      (procs_mu, systematic_names_nominal),
                #'el_presel_ss':   (procs_el, systematic_names_nominal),
                #'mu_presel_ss':   (procs_mu, systematic_names_nominal),
                # with tau POG selection
                #'pog_mu_presel':  (procs_mu, systematic_names_nominal),
                #'pog_mu_pass':    (procs_mu, systematic_names_nominal),
                #'pog_mu_pass_ss': (procs_mu, systematic_names_nominal),
                #'pog_mu_fail':    (procs_mu, systematic_names_nominal),
                # with addition of no DY mass, no match to b-tag (could add a cut on small MT)
                #'adv_el_sel':       (procs_el_3ch, systematic_names_toppt), #systematic_names_pu),
                #'adv_el_sel_Sign4': (procs_el_3ch, systematic_names_toppt), #systematic_names_pu), # this is done with a hack in the following, watch closely

                # in adv selection jets cross-cleaned from tau
                'adv_mu_sel_Tight':          (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_Tight_ss':       (procs_mu_3ch, systematic_names_toppt),

                # preselection with all IDs...
                #'adv_mu_sel_preLoose':       (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_Loose':          (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_Loose_ss':       (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_Loose_lj':       (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_Loose_lj_ss':    (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_Loose_ljout':    (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_Loose_ljout_ss': (procs_mu_3ch, systematic_names_toppt),
                #'sel_mu_min':                (procs_mu,     systematic_names_all),
                #'sel_mu_min_ss':             (procs_mu,     systematic_names_all),
                #'sel_mu_min_lj':             (procs_mu,     systematic_names_all),
                #'sel_mu_min_lj_ss':          (procs_mu,     systematic_names_all),
                #'sel_mu_min_ljout':          (procs_mu,     systematic_names_all),
                #'sel_mu_min_ljout_ss':       (procs_mu,     systematic_names_all),
                # the advanced-advanced selection -- scanning the optimization
                # b-jets are Loose, tau is Medium
                # TODO: I need that uniform selection and distributions for this.....
                #'adv_mu_sel_preLoose2':       (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_LooseM':          (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_LooseM_ss':       (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_LooseM_lj':       (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_LooseM_lj_ss':    (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_LooseM_ljout':    (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_LooseM_ljout_ss': (procs_mu_3ch, systematic_names_toppt),

                # data-driven QCD only on "loose" selections
                'adv_mu_sel_Loose_alliso':          (procs_mu_3ch, systematic_names_toppt),
                'adv_mu_sel_Loose_alliso_ss':       (procs_mu_3ch, systematic_names_toppt),
                #'adv_mu_sel_Loose_alliso_lj':       (procs_mu_3ch, systematic_names_toppt),
                #'adv_mu_sel_Loose_alliso_lj_ss':    (procs_mu_3ch, systematic_names_toppt),
                #'adv_mu_sel_Loose_alliso_ljout':    (procs_mu_3ch, systematic_names_toppt),
                #'adv_mu_sel_Loose_alliso_ljout_ss': (procs_mu_3ch, systematic_names_toppt),
                #'sel_mu_min_alliso':                (procs_mu,     systematic_names_toppt),
                #'sel_mu_min_alliso_ss':             (procs_mu,     systematic_names_toppt),
                #'sel_mu_min_alliso_lj':             (procs_mu,     systematic_names_toppt),
                #'sel_mu_min_alliso_lj_ss':          (procs_mu,     systematic_names_toppt),
                #'sel_mu_min_alliso_ljout':          (procs_mu,     systematic_names_toppt),
                #'sel_mu_min_alliso_ljout_ss':       (procs_mu,     systematic_names_toppt),

                #'sel_mu_min_medtau':  (procs_mu, systematic_names_nominal), #systematic_names_pu_toppt), # minimum selection with Medium taus -- hopefully it will reduce QCD
                # control selections: WJets, DY mumu and tautau, tt elmu
                'ctr_mu_wjet':        (procs_mu, systematic_names_pu),
                #'ctr_mu_wjet_ss':    (procs_mu, systematic_names_nominal),
                #'ctr_mu_wjet_old':   (procs_mu, systematic_names_nominal),
                'ctr_mu_dy_mumu':     (procs_mu, systematic_names_pu),
                #'ctr_mu_dy_mumu_ss':  (procs_mu, systematic_names_nominal),
                'ctr_mu_dy_tt':       (procs_mu, systematic_names_pu),
                'ctr_mu_dy_tt_ss':    (procs_mu, systematic_names_pu),
                #'ctr_mu_dy_SV_tt':    (procs_mu, systematic_names_nominal), #systematic_names_all),
                #'ctr_mu_dy_SV_tt_ss': (procs_mu, systematic_names_nominal), #systematic_names_all),
                'ctr_mu_tt_em':       (procs_elmu, systematic_names_pu_toppt),
                }

    # for the special shape study of backgrounds
    if isWJets or isDY:
        channels_usual.update({
                   'mu_sel_nomet':         (procs_mu, ['NOMINAL']),
                   'mu_sel_onlymet':       (procs_mu, ['NOMINAL']),
                   'mu_sel_nobtag':        (procs_mu, ['NOMINAL']),
                   'mu_sel_onlybtag':      (procs_mu, ['NOMINAL']),
                   'mu_sel_no3jets':       (procs_mu, ['NOMINAL']),
                   'mu_sel_only3jets':     (procs_mu, ['NOMINAL'])
                   })

    channels_mu_only = {
                'mu_presel':      (procs_mu_3ch, systematic_names_nominal),
                'mu_sel':         (procs_mu_3ch, systematic_names_nominal),
                'mu_lj':          (procs_mu_3ch, systematic_names_nominal),
                'mu_lj_out':      (procs_mu_3ch, systematic_names_nominal),
                'mu_lj_ss':       (procs_mu_3ch, systematic_names_nominal),
                'mu_lj_out_ss':   (procs_mu_3ch, systematic_names_nominal),
                }

    # OPTIMIZATION
    # optimization scan over tau and b ID selections in mu-tau
    # tau ID = presel, loose, medium, tight
    # N b jets = 1L0M, 2L0M, 1L1M, 2L1M, 2L2M (2 means 2 or more)
    # = 20 selections
    # jets are cross-dR from tau
    # tau is OS to lep -- in principle need SS for QCD estimation... but works very well now (shape-wise) so skipping it here
    channels_mutau_optimization_scan = {
                'optmu_presel':         (procs_mu, systematic_names_nominal),
                'optmu_presel_ss':      (procs_mu, systematic_names_nominal),

                'optmu_presel_1L0M':      (procs_mu, systematic_names_nominal),
                'optmu_presel_2L0M':      (procs_mu, systematic_names_nominal),
                'optmu_presel_1L1M':      (procs_mu, systematic_names_nominal),
                'optmu_presel_2L1M':      (procs_mu, systematic_names_nominal),
                'optmu_presel_2L2M':      (procs_mu, systematic_names_nominal),

                'optmu_loose_1L0M':      (procs_mu, systematic_names_nominal),
                'optmu_loose_2L0M':      (procs_mu, systematic_names_nominal),
                'optmu_loose_1L1M':      (procs_mu, systematic_names_nominal),
                'optmu_loose_2L1M':      (procs_mu, systematic_names_nominal),
                'optmu_loose_2L2M':      (procs_mu, systematic_names_nominal),

                'optmu_medium_1L0M':      (procs_mu, systematic_names_nominal),
                'optmu_medium_2L0M':      (procs_mu, systematic_names_nominal),
                'optmu_medium_1L1M':      (procs_mu, systematic_names_nominal),
                'optmu_medium_2L1M':      (procs_mu, systematic_names_nominal),
                'optmu_medium_2L2M':      (procs_mu, systematic_names_nominal),

                'optmu_tight_1L0M':      (procs_mu, systematic_names_nominal),
                'optmu_tight_2L0M':      (procs_mu, systematic_names_nominal),
                'optmu_tight_1L1M':      (procs_mu, systematic_names_nominal),
                'optmu_tight_2L1M':      (procs_mu, systematic_names_nominal),
                'optmu_tight_2L2M':      (procs_mu, systematic_names_nominal),

                'optmu_presel_1L0M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_presel_2L0M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_presel_1L1M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_presel_2L1M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_presel_2L2M_ss':      (procs_mu, systematic_names_nominal),

                'optmu_loose_1L0M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_loose_2L0M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_loose_1L1M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_loose_2L1M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_loose_2L2M_ss':      (procs_mu, systematic_names_nominal),

                'optmu_medium_1L0M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_medium_2L0M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_medium_1L1M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_medium_2L1M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_medium_2L2M_ss':      (procs_mu, systematic_names_nominal),

                'optmu_tight_1L0M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_tight_2L0M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_tight_1L1M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_tight_2L1M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_tight_2L2M_ss':      (procs_mu, systematic_names_nominal),
               }

    '''
    optimization shows 2L1M as best choice for ttbar
    Tight tau ID gives about the same as old event yield ratios
    but about x2 more events
       14.5k signal VS 8.6k old signal
       27.4k all events VS 14.7k old events
       of which 8% is QCD, which is due to low pt cuts

    loose and medium WP don't give much more signal:
       18.8k with loose, 16.6k with medium

    but increase various backgrounds:
       50.5k all events with loose, 36.6k with medium

    2M1L has a bit more involved b-tagging correction, which I won't do right now [TODO]

    run a test to show the improved performance in the fit
    use 2L1M+tight for that, for mutau and eltau,
    add 2L1M+tight+pt cuts selection [should reduce QCD without cutting the signal too much]
    just in case add 2L1M+loose/tight for mutau only -- main computation is lj parameter, which is done only for events passing main selection (OS and SS tau)
    hence recomputing should not be that big

    also I need to run WJets and DY controls with PU systematics -- for WJets normalization and QCD OS/SS estimation
    and I need anti-iso region for QCD estimation
    -- I'll add a complete all_iso selection for each target selection,
       hopefully there are not so many anti-iso events (it seems so from the tests)

    thus there are:
        2L1M presel + (loose + medium) + tight * (1+2lj + 1cut) = 7 selections for mu-tau
        2L1M presel + tight * (1+2lj)                    = 4 for el-tau
                                        * 2ss * 2 all iso
                                                         = 44 selections
        WJets * 2ss * 2all iso + DY control              = 5
                                                         = 49 selections

    -- I do * 2ss * 2 all iso where I want to check QCD         (all main selections and WJets for QCD-WJets estimation)
       1 + 2lj only for the target selections for fit -- tight
       also only they will have all systematics (lj bins and all events for fit without lj)
       el has only the target                (for full fit)
       mutau has selection with cuts         (without lj, to see event yield posibilities)
       WJets and DY have PU systematics      (for WJets-DY estimation -- to cover the issue with bend at high Mt)

    I had 42 selections for optimization
    (with nominal only systematics and no all iso)
    let's hope it'll fly

    only tight selections have all systematics
    the rest are nominal,
    except DY and WJets os iso -- with PU Up/Down
    '''

    selection_definitions = {
        'channels_optimized_alliso' : {
                'optmu_alliso_presel_2L1M':         (procs_mu, systematic_names_nominal),
                'optmu_alliso_presel_2L1M_ss':      (procs_mu, systematic_names_nominal),
                #'optmu_alliso_loose_2L1M':          (procs_mu, systematic_names_nominal),
                #'optmu_alliso_loose_2L1M_ss':       (procs_mu, systematic_names_nominal),
                #'optmu_alliso_medium_2L1M':         (procs_mu, systematic_names_nominal),
                #'optmu_alliso_medium_2L1M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_alliso_tight_2L1M_cuts':     (procs_mu, systematic_names_nominal),
                'optmu_alliso_tight_2L1M_cuts_ss':  (procs_mu, systematic_names_nominal),
                'optmu_alliso_tight_2L1M':          (procs_mu, systematic_names_nominal),
                'optmu_alliso_tight_2L1M_ss':       (procs_mu, systematic_names_nominal),
                'optmu_alliso_tight_2L1M_lj':       (procs_mu, systematic_names_nominal),
                'optmu_alliso_tight_2L1M_lj_ss':    (procs_mu, systematic_names_nominal),
                'optmu_alliso_tight_2L1M_ljout':    (procs_mu, systematic_names_nominal),
                'optmu_alliso_tight_2L1M_ljout_ss': (procs_mu, systematic_names_nominal),

                'optel_alliso_presel_2L1M':         (procs_el, systematic_names_nominal),
                'optel_alliso_presel_2L1M_ss':      (procs_el, systematic_names_nominal),
                'optel_alliso_tight_2L1M':          (procs_el, systematic_names_nominal),
                'optel_alliso_tight_2L1M_ss':       (procs_el, systematic_names_nominal),
                'optel_alliso_tight_2L1M_lj':       (procs_el, systematic_names_nominal),
                'optel_alliso_tight_2L1M_lj_ss':    (procs_el, systematic_names_nominal),
                'optel_alliso_tight_2L1M_ljout':    (procs_el, systematic_names_nominal),
                'optel_alliso_tight_2L1M_ljout_ss': (procs_el, systematic_names_nominal),
        },

        'channels_optimized' : {
                'optmu_presel_2L1M':         (procs_mu, systematic_names_nominal),
                'optmu_presel_2L1M_ss':      (procs_mu, systematic_names_nominal),
                #'optmu_loose_2L1M':          (procs_mu, systematic_names_nominal),
                #'optmu_loose_2L1M_ss':       (procs_mu, systematic_names_nominal),
                #'optmu_medium_2L1M':         (procs_mu, systematic_names_nominal),
                #'optmu_medium_2L1M_ss':      (procs_mu, systematic_names_nominal),
                'optmu_tight_2L1M_cuts':     (procs_mu, systematic_names_nominal),
                'optmu_tight_2L1M_cuts_ss':  (procs_mu, systematic_names_nominal),

                'optmu_tight_2L1M':          (procs_mu, systematic_names_all),
                'optmu_tight_2L1M_ss':       (procs_mu, systematic_names_all),
                'optmu_tight_2L1M_lj':       (procs_mu, systematic_names_all),
                'optmu_tight_2L1M_lj_ss':    (procs_mu, systematic_names_all),
                'optmu_tight_2L1M_ljout':    (procs_mu, systematic_names_all),
                'optmu_tight_2L1M_ljout_ss': (procs_mu, systematic_names_all),

                'optmu_loose_2L1M':          (procs_mu, systematic_names_nominal),
                'optmu_loose_2L1M_ss':       (procs_mu, systematic_names_nominal),
                'optmu_loose_2L1M_lj':       (procs_mu, systematic_names_nominal),
                'optmu_loose_2L1M_lj_ss':    (procs_mu, systematic_names_nominal),
                'optmu_loose_2L1M_ljout':    (procs_mu, systematic_names_nominal),
                'optmu_loose_2L1M_ljout_ss': (procs_mu, systematic_names_nominal),

                'optel_presel_2L1M':         (procs_el, systematic_names_nominal),
                'optel_presel_2L1M_ss':      (procs_el, systematic_names_nominal),
                'optel_tight_2L1M':          (procs_el, systematic_names_all),
                'optel_tight_2L1M_ss':       (procs_el, systematic_names_all),
                'optel_tight_2L1M_lj':       (procs_el, systematic_names_all),
                'optel_tight_2L1M_lj_ss':    (procs_el, systematic_names_all),
                'optel_tight_2L1M_ljout':    (procs_el, systematic_names_all),
                'optel_tight_2L1M_ljout_ss': (procs_el, systematic_names_all),

                'ctr_mu_wjet':              (procs_mu, systematic_names_pu),
                'ctr_mu_wjet_ss':           (procs_mu, systematic_names_pu),
                'ctr_alliso_mu_wjet':       (procs_mu, systematic_names_nominal),
                'ctr_alliso_mu_wjet_ss':    (procs_mu, systematic_names_nominal),
                'ctr_mu_dy_mumu':           (procs_mu, systematic_names_pu),

                'ctr_mu_tt_em':             (procs_elmu, systematic_names_pu_toppt),
                'ctr_old_mu_sel':           (procs_mu, systematic_names_nominal),     # testing issue with event yield advantage
                'ctr_old_mu_sel_ss':        (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_lj':        (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_lj_ss':     (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ljout':     (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ljout_ss':  (procs_mu, systematic_names_nominal),
        },

        'channels_optimized_old_full_sys' : {
                'ctr_mu_wjet':              (procs_mu, systematic_names_pu),
                'ctr_mu_wjet_ss':           (procs_mu, systematic_names_pu),
                'ctr_el_wjet':              (procs_el, systematic_names_pu),
                'ctr_el_wjet_ss':           (procs_el, systematic_names_pu),
                'ctr_alliso_mu_wjet':       (procs_mu, systematic_names_nominal),
                'ctr_alliso_mu_wjet_ss':    (procs_mu, systematic_names_nominal),

                'ctr_mu_dy_mumu':           (procs_mu, systematic_names_pu),
                'ctr_mu_dy_elel':           (procs_el, systematic_names_pu),
                #'ctr_mu_tt_em':             (procs_elmu, systematic_names_pu_toppt),
                'ctr_mu_tt_em_close':       (procs_elmu, systematic_names_all_with_th),


                #'ctr_old_mu_presel':        (procs_mu_3ch_fbw, systematic_names_pu_toppt),     # testing issue with event yield advantage
                #'ctr_old_mu_presel_ss':     (procs_mu_3ch_fbw, systematic_names_pu_toppt),

                'ctr_old_mu_presel':        (procs_mu, systematic_names_pu_toppt),     # testing issue with event yield advantage
                'ctr_old_mu_presel_ss':     (procs_mu, systematic_names_pu_toppt),

                #'ctr_old_mu_presel_alliso':        (procs_mu, systematic_names_nominal),
                #'ctr_old_mu_presel_alliso_ss':     (procs_mu, systematic_names_nominal),
                #'ctr_old_el_presel_alliso':        (procs_el, systematic_names_nominal),
                #'ctr_old_el_presel_alliso_ss':     (procs_el, systematic_names_nominal),

                'ctr_old_mu_selVloose_alliso':      (procs_mu, systematic_names_nominal),
                'ctr_old_mu_selVloose_alliso_ss':   (procs_mu, systematic_names_nominal),
                'ctr_old_el_selVloose_alliso':      (procs_el, systematic_names_nominal),
                'ctr_old_el_selVloose_alliso_ss':   (procs_el, systematic_names_nominal),

                # testing issue with event yield advantage
                'ctr_old_mu_selVloose':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_ss':  (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel':           (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ss':        (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_lj':        (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_lj_ss':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ljout':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ljout_ss':  (procs_mu, systematic_names_all_with_th),

                'ctr_old_mu_sel_tauSign3':           (procs_mu_3ch_fbw, systematic_names_nominal),
                'ctr_old_mu_sel_tauSign3_ss':        (procs_mu_3ch_fbw, systematic_names_nominal),

                'ctr_old_el_presel':        (procs_el, systematic_names_pu_toppt),     # testing issue with event yield advantage
                'ctr_old_el_presel_ss':     (procs_el, systematic_names_pu_toppt),
                'ctr_old_el_selVloose':     (procs_el, systematic_names_all_with_th),
                'ctr_old_el_selVloose_ss':  (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel':           (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_ss':        (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_lj':        (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_lj_ss':     (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_ljout':     (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_ljout_ss':  (procs_el, systematic_names_all_with_th),
        },

        'channels_full_sys_electron_selections' : {
                'ctr_el_wjet':              (procs_el, systematic_names_pu),
                'ctr_el_wjet_ss':           (procs_el, systematic_names_pu),
                'ctr_old_el_presel':        (procs_el, systematic_names_pu_toppt),     # testing issue with event yield advantage
                'ctr_old_el_presel_ss':     (procs_el, systematic_names_pu_toppt),
                'ctr_old_el_selVloose':     (procs_el, systematic_names_all_with_th),
                'ctr_old_el_selVloose_ss':  (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel':           (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_ss':        (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_lj':        (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_lj_ss':     (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_ljout':     (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_ljout_ss':  (procs_el, systematic_names_all_with_th),
        },

        'channels_full_sys_muon_selections' : {
                'ctr_mu_wjet':                    (procs_mu, systematic_names_pu),
                'ctr_mu_wjet_ss':                 (procs_mu, systematic_names_pu),
                'ctr_old_mu_presel':              (procs_mu, systematic_names_pu_toppt),     # testing issue with event yield advantage
                'ctr_old_mu_presel_ss':           (procs_mu, systematic_names_pu_toppt),
                'ctr_old_mu_selVloose':           (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_ss':        (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_lj':        (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_lj_ss':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_ljout':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_ljout_ss':  (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel':                 (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ss':              (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_lj':              (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_lj_ss':           (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ljout':           (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ljout_ss':        (procs_mu, systematic_names_all_with_th),
        },

        'channels_full_sys_muon_sel' : {
                'ctr_old_mu_selVloose':           (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_ss':        (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_lj':        (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_lj_ss':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_ljout':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_ljout_ss':  (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel':                 (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ss':              (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_lj':              (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_lj_ss':           (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ljout':           (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ljout_ss':        (procs_mu, systematic_names_all_with_th),
        },

        'channels_nom_sys_muon_sel' : {
                'ctr_old_mu_selVloose':           (procs_mu, systematic_names_nominal),
                'ctr_old_mu_selVloose_ss':        (procs_mu, systematic_names_nominal),
                'ctr_old_mu_selVloose_lj':        (procs_mu, systematic_names_nominal),
                'ctr_old_mu_selVloose_lj_ss':     (procs_mu, systematic_names_nominal),
                'ctr_old_mu_selVloose_ljout':     (procs_mu, systematic_names_nominal),
                'ctr_old_mu_selVloose_ljout_ss':  (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel':                 (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ss':              (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_lj':              (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_lj_ss':           (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ljout':           (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ljout_ss':        (procs_mu, systematic_names_nominal),
        },

        'channels_nom_sys_muon_sel_main' : {
                'ctr_old_mu_sel':                 (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ss':              (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_lj':              (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_lj_ss':           (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ljout':           (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ljout_ss':        (procs_mu, systematic_names_nominal),
        },

        'channels_full_sys_lep_selections' : {
                #'ctr_el_wjet':              (procs_el, systematic_names_pu),
                #'ctr_el_wjet_ss':           (procs_el, systematic_names_pu),
                #'ctr_old_el_presel':        (procs_el, systematic_names_pu_toppt),     # testing issue with event yield advantage
                #'ctr_old_el_presel_ss':     (procs_el, systematic_names_pu_toppt),
                'ctr_old_el_selVloose':     (procs_el, systematic_names_all_with_th),
                'ctr_old_el_selVloose_ss':  (procs_el, systematic_names_all_with_th),
                'ctr_old_el_selVloose_lj':     (procs_el, systematic_names_all_with_th),
                'ctr_old_el_selVloose_lj_ss':  (procs_el, systematic_names_all_with_th),
                'ctr_old_el_selVloose_ljout':     (procs_el, systematic_names_all_with_th),
                'ctr_old_el_selVloose_ljout_ss':  (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel':           (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_ss':        (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_lj':        (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_lj_ss':     (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_ljout':     (procs_el, systematic_names_all_with_th),
                'ctr_old_el_sel_ljout_ss':  (procs_el, systematic_names_all_with_th),

                #'ctr_mu_wjet':              (procs_mu, systematic_names_pu),
                #'ctr_mu_wjet_ss':           (procs_mu, systematic_names_pu),
                #'ctr_old_mu_presel':        (procs_mu, systematic_names_pu_toppt),     # testing issue with event yield advantage
                #'ctr_old_mu_presel_ss':     (procs_mu, systematic_names_pu_toppt),
                'ctr_old_mu_selVloose':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_ss':  (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_lj':        (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_lj_ss':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_ljout':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_selVloose_ljout_ss':  (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel':           (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ss':        (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_lj':        (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_lj_ss':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ljout':     (procs_mu, systematic_names_all_with_th),
                'ctr_old_mu_sel_ljout_ss':  (procs_mu, systematic_names_all_with_th),
        },

        'channels_nominal_full_lep_selections' : {
                'ctr_old_el_presel':        (procs_el, systematic_names_nominal),     # testing issue with event yield advantage
                'ctr_old_el_presel_ss':     (procs_el, systematic_names_nominal),
                'ctr_old_el_sel':           (procs_el, systematic_names_nominal),
                'ctr_old_el_sel_ss':        (procs_el, systematic_names_nominal),
                'ctr_old_el_sel_lj':        (procs_el, systematic_names_nominal),
                'ctr_old_el_sel_lj_ss':     (procs_el, systematic_names_nominal),
                'ctr_old_el_sel_ljout':     (procs_el, systematic_names_nominal),
                'ctr_old_el_sel_ljout_ss':  (procs_el, systematic_names_nominal),

                'ctr_old_mu_presel':        (procs_mu, systematic_names_nominal),
                'ctr_old_mu_presel_ss':     (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel':           (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ss':        (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_lj':        (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_lj_ss':     (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ljout':     (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ljout_ss':  (procs_mu, systematic_names_nominal),
        },

        'channels_nominal_full_mu_selections' : {
                'ctr_old_mu_presel':        (procs_mu, systematic_names_nominal),
                'ctr_old_mu_presel_ss':     (procs_mu, systematic_names_nominal),
                'ctr_old_mu_presel2b':      (procs_mu, systematic_names_nominal),
                'ctr_old_mu_presel2b_ss':   (procs_mu, systematic_names_nominal),
                'ctr_old_mu_selVloose':     (procs_mu, systematic_names_nominal),
                'ctr_old_mu_selVloose_ss':  (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel':           (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ss':        (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_lj':        (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_lj_ss':     (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ljout':     (procs_mu, systematic_names_nominal),
                'ctr_old_mu_sel_ljout_ss':  (procs_mu, systematic_names_nominal),
        },

        'channels_nominal_full_el_selections' : {
                'ctr_old_el_presel':        (procs_el, systematic_names_nominal),
                'ctr_old_el_presel_ss':     (procs_el, systematic_names_nominal),
                'ctr_old_el_presel2b':      (procs_el, systematic_names_nominal),
                'ctr_old_el_presel2b_ss':   (procs_el, systematic_names_nominal),
                'ctr_old_el_selVloose':     (procs_el, systematic_names_nominal),
                'ctr_old_el_selVloose_ss':  (procs_el, systematic_names_nominal),
                'ctr_old_el_sel':           (procs_el, systematic_names_nominal),
                'ctr_old_el_sel_ss':        (procs_el, systematic_names_nominal),
                'ctr_old_el_sel_lj':        (procs_el, systematic_names_nominal),
                'ctr_old_el_sel_lj_ss':     (procs_el, systematic_names_nominal),
                'ctr_old_el_sel_ljout':     (procs_el, systematic_names_nominal),
                'ctr_old_el_sel_ljout_ss':  (procs_el, systematic_names_nominal),
        },

        'channels_refact_sys_lep_selections' : {
                #'ctr_el_wjet':              (procs_el, systematic_names_pu),
                #'ctr_el_wjet_ss':           (procs_el, systematic_names_pu),
                #'ctr_old_el_presel':        (procs_el, systematic_names_pu_toppt),     # testing issue with event yield advantage
                #'ctr_old_el_presel_ss':     (procs_el, systematic_names_pu_toppt),
                'ctr_old_el_selVloose':     (procs_el, systematic_names_th_renorm_refact),
                'ctr_old_el_selVloose_ss':  (procs_el, systematic_names_th_renorm_refact),
                'ctr_old_el_selVloose_lj':     (procs_el, systematic_names_th_renorm_refact),
                'ctr_old_el_selVloose_lj_ss':  (procs_el, systematic_names_th_renorm_refact),
                'ctr_old_el_selVloose_ljout':     (procs_el, systematic_names_th_renorm_refact),
                'ctr_old_el_selVloose_ljout_ss':  (procs_el, systematic_names_th_renorm_refact),
                'ctr_old_el_sel':           (procs_el, systematic_names_th_renorm_refact),
                'ctr_old_el_sel_ss':        (procs_el, systematic_names_th_renorm_refact),
                'ctr_old_el_sel_lj':        (procs_el, systematic_names_th_renorm_refact),
                'ctr_old_el_sel_lj_ss':     (procs_el, systematic_names_th_renorm_refact),
                'ctr_old_el_sel_ljout':     (procs_el, systematic_names_th_renorm_refact),
                'ctr_old_el_sel_ljout_ss':  (procs_el, systematic_names_th_renorm_refact),

                #'ctr_mu_wjet':              (procs_mu, systematic_names_pu),
                #'ctr_mu_wjet_ss':           (procs_mu, systematic_names_pu),
                #'ctr_old_mu_presel':        (procs_mu, systematic_names_pu_toppt),     # testing issue with event yield advantage
                #'ctr_old_mu_presel_ss':     (procs_mu, systematic_names_pu_toppt),
                'ctr_old_mu_selVloose':     (procs_mu, systematic_names_th_renorm_refact),
                'ctr_old_mu_selVloose_ss':  (procs_mu, systematic_names_th_renorm_refact),
                'ctr_old_mu_selVloose_lj':        (procs_mu, systematic_names_th_renorm_refact),
                'ctr_old_mu_selVloose_lj_ss':     (procs_mu, systematic_names_th_renorm_refact),
                'ctr_old_mu_selVloose_ljout':     (procs_mu, systematic_names_th_renorm_refact),
                'ctr_old_mu_selVloose_ljout_ss':  (procs_mu, systematic_names_th_renorm_refact),
                'ctr_old_mu_sel':           (procs_mu, systematic_names_th_renorm_refact),
                'ctr_old_mu_sel_ss':        (procs_mu, systematic_names_th_renorm_refact),
                'ctr_old_mu_sel_lj':        (procs_mu, systematic_names_th_renorm_refact),
                'ctr_old_mu_sel_lj_ss':     (procs_mu, systematic_names_th_renorm_refact),
                'ctr_old_mu_sel_ljout':     (procs_mu, systematic_names_th_renorm_refact),
                'ctr_old_mu_sel_ljout_ss':  (procs_mu, systematic_names_th_renorm_refact),
        },

        'channels_control_regions' : {
                'ctr_mu_wjet':              (procs_mu, systematic_names_pu),
                'ctr_mu_wjet_ss':           (procs_mu, systematic_names_pu),
                'ctr_el_wjet':              (procs_el, systematic_names_pu),
                'ctr_el_wjet_ss':           (procs_el, systematic_names_pu),
                'ctr_alliso_mu_wjet':       (procs_mu, systematic_names_nominal),
                'ctr_alliso_mu_wjet_ss':    (procs_mu, systematic_names_nominal),

                'ctr_mu_dy_mumu':           (procs_mu, systematic_names_pu),
                'ctr_mu_dy_elel':           (procs_el, systematic_names_pu),
                'ctr_mu_tt_em':             (procs_elmu, systematic_names_pu_toppt),
                #'ctr_mu_tt_em_close':       (procs_elmu, systematic_names_all_with_th), # for the ratio

                'ctr_old_mu_presel_alliso':        (procs_mu, systematic_names_nominal),
                'ctr_old_mu_presel_alliso_ss':     (procs_mu, systematic_names_nominal),
                'ctr_old_el_presel_alliso':        (procs_el, systematic_names_nominal),
                'ctr_old_el_presel_alliso_ss':     (procs_el, systematic_names_nominal),

                'ctr_old_mu_selVloose_alliso':      (procs_mu, systematic_names_nominal),
                'ctr_old_mu_selVloose_alliso_ss':   (procs_mu, systematic_names_nominal),
                'ctr_old_el_selVloose_alliso':      (procs_el, systematic_names_nominal),
                'ctr_old_el_selVloose_alliso_ss':   (procs_el, systematic_names_nominal),
        },

        'channels_control_regions_wjets' : {
                'ctr_mu_wjet':              (procs_mu, systematic_names_pu),
                'ctr_mu_wjet_ss':           (procs_mu, systematic_names_pu),
                'ctr_el_wjet':              (procs_el, systematic_names_pu),
                'ctr_el_wjet_ss':           (procs_el, systematic_names_pu),
                'ctr_alliso_mu_wjet':       (procs_mu, systematic_names_nominal),
                'ctr_alliso_mu_wjet_ss':    (procs_mu, systematic_names_nominal),
        },

        'channels_presels' : {
                'ctr_old_el_presel':        (procs_el, systematic_names_pu_toppt),     # testing issue with event yield advantage
                'ctr_old_el_presel_ss':     (procs_el, systematic_names_pu_toppt),
                'ctr_old_mu_presel':        (procs_mu, systematic_names_pu_toppt),     # testing issue with event yield advantage
                'ctr_old_mu_presel_ss':     (procs_mu, systematic_names_pu_toppt),
        },

        'channels_control_regions_elmu' : {
                'ctr_mu_tt_em':             (procs_elmu, systematic_names_pu_toppt),
                'ctr_mu_tt_em_close':       (procs_elmu, systematic_names_pu_toppt),
                'ctr_el_tt_em':             (procs_elmu, systematic_names_pu_toppt),
                'ctr_el_tt_em_close':       (procs_elmu, systematic_names_pu_toppt),
        },

        'channels_control_regions_elmu_allsys' : {
                'ctr_mu_tt_em':             (procs_elmu, systematic_names_all_with_th),
                'ctr_mu_tt_em_close':       (procs_elmu, systematic_names_all_with_th),
                'ctr_el_tt_em':             (procs_elmu, systematic_names_all_with_th),
                'ctr_el_tt_em_close':       (procs_elmu, systematic_names_all_with_th),
        },

        'channels_control_regions_dy' : {
                'ctr_mu_dy_mumu':           (procs_mu, systematic_names_pu),
                'ctr_mu_dy_elel':           (procs_el, systematic_names_pu),
        }
    }

    print systematic_names_all_with_th
    print systematic_names_pu_toppt
    print systematic_names_pu
    print systematic_names_nominal

    #channels = channels_mu_only
    #channels = channels_usual
    #channels, with_bSF = channels_mutau_optimization_scan, False
    #selected_channels = channels_optimized_old_full_sys
    #selected_channels = channels_full_sys_electron_selections
    #selected_channels = selection_definitions[channels_to_select]

    # find all requested systematics
    #requested_systematics  = set()
    #requested_objects      = set()
    #for k, (ch_name, systematic_names) in selected_channels.items():
    #  for systematic_name in systematic_names:
    #    requested_systematics.add(systematic_name)
    #    if systematic_name in all_systematic_objects:
    #        requested_objects.add(systematic_name)

    with_JER_sys   = with_JER   = isMC and True and not isTT_systematic
    with_JES_sys   = with_JES   = isMC and True and not isTT_systematic
    with_TauES_sys = with_TauES = isMC and True and not isTT_systematic

    with_bSF_sys = with_bSF # and ('bSFUp' in requested_systematics or 'bSFDown' in requested_systematics)

    with_AlphaS_sys  = True and isTT and not isTT_systematic # and ('AlphaSUp' in requested_systematics or 'AlphaSDown' in requested_systematics)
    with_Frag_sys    = True and isTT and not isTT_systematic # and ('FragUp'   in requested_systematics or 'FragDown'   in requested_systematics)
    with_MEscale_sys = True and isTT and not isTT_systematic # and ('MfUp'     in requested_systematics or 'MfDown'   in requested_systematics)
    with_PDF_sys     = True and isTT and not isTT_systematic # and ('PDF_TRIGGER'      in requested_systematics)

    #SystematicJets = namedtuple('Jets', 'nom sys_JERUp sys_JERDown sys_JESUp sys_JESDown sys_bUp sys_bDown')
    # all requested jet cuts
    JetBSplit = namedtuple('BJets', 'medium loose rest taumatched lepmatched')
    TauSplit  = namedtuple('Taus', 'candidates medium')
    #JetCutsPerSystematic = namedtuple('Jets', 'lowest cuts old') # TODO add jets loose
    #JetCutsPerSystematic = namedtuple('Jets', 'old old_alliso') # TODO add jets loose
    #TauCutsPerSystematic = namedtuple('Taus', 'lowest loose cuts old oldVloose presel_alliso presel')

    # get the corrections to JES for PS systematic datasets:
    correction_jes_endcap_b = None
    correction_jes_barrel_b = None
    correction_jes_endcap_c = None
    correction_jes_barrel_c = None
    correction_jes_endcap_udsg = None
    correction_jes_barrel_udsg = None

    # the usual dtag
    dtag_usual = '_'.join(dtag.split('/')[-1].split('.root')[0].split('_')[:-1])
    if isTT_systematic:
        logging.info("test TT sys %s" % dtag_usual)

        the_jes_corrections_fname = "studies/jes_bSF_corrections/corrections_to_nominal_jes.root"
        the_jes_corrections_file  = TFile(the_jes_corrections_fname)

        # Tune Up Down might not be present
        jes_cor_barrel_b = 'correction_' + dtag_usual + '_jes_cor_1d_gen_reco_barrel_b'
        if the_jes_corrections_file.Get(jes_cor_barrel_b):
            logging.info("test TT sys loading corrections to JES")

            correction_jes_endcap_b    = the_jes_corrections_file.Get('correction_' + dtag_usual + '_jes_cor_1d_gen_reco_endcap_b')
            correction_jes_barrel_b    = the_jes_corrections_file.Get('correction_' + dtag_usual + '_jes_cor_1d_gen_reco_barrel_b')
            correction_jes_endcap_c    = the_jes_corrections_file.Get('correction_' + dtag_usual + '_jes_cor_1d_gen_reco_endcap_c')
            correction_jes_barrel_c    = the_jes_corrections_file.Get('correction_' + dtag_usual + '_jes_cor_1d_gen_reco_barrel_c')
            correction_jes_endcap_udsg = the_jes_corrections_file.Get('correction_' + dtag_usual + '_jes_cor_1d_gen_reco_endcap_udsg')
            correction_jes_barrel_udsg = the_jes_corrections_file.Get('correction_' + dtag_usual + '_jes_cor_1d_gen_reco_barrel_udsg')


    #set_bSF_effs_for_dtag_usual(dtag_usual)
    if with_bSF:
        logger.write(' '.join(str(id(h)) for h in (bEff_histo_b, bEff_histo_c, bEff_histo_udsg)) + '\n')
        logging.info("test b SF weight with eff histo: %s" % bEff_histo_b.GetName())

        # test the calculator
        jet_pt  = 40.
        jet_eta = 1.1
        flavId = 5
        jet_weight_bSF_nom_tagged, sf_t, eff_t     = calc_btag_sf_weight(True,  flavId, jet_pt, jet_eta)
        jet_weight_bSF_nom_not_tagged, sf_n, eff_n = calc_btag_sf_weight(False, flavId, jet_pt, jet_eta)
        logging.info("test b SF weight for %d: %10.6f %10.6f,  %10.6f %10.6f" % (flavId, jet_weight_bSF_nom_not_tagged, jet_weight_bSF_nom_not_tagged, sf_t, eff_t))
        flavId = 1
        jet_weight_bSF_nom_tagged, sf_t, eff_t     = calc_btag_sf_weight(True,  flavId, jet_pt, jet_eta)
        jet_weight_bSF_nom_not_tagged, sf_n, eff_n = calc_btag_sf_weight(False, flavId, jet_pt, jet_eta)
        logging.info("test b SF weight for %d: %10.6f %10.6f,  %10.6f %10.6f" % (flavId, jet_weight_bSF_nom_not_tagged, jet_weight_bSF_nom_not_tagged, sf_t, eff_t))

        the_updown_bSF_efficiencies_filename = "studies/jes_bSF_corrections/%s.root" % dtag_usual
        if isTT_systematic and isfile(the_updown_bSF_efficiencies_filename):
            logging.info("test TT sys %s" % dtag_usual)

            # change the efficiency histograms
            the_updown_bSF_efficiencies = TFile(the_updown_bSF_efficiencies_filename)
            # the numerator
            bEff_histo_b_num    = the_updown_bSF_efficiencies.Get("btag_b_hadronFlavour_candidates_tagged")
            bEff_histo_c_num    = the_updown_bSF_efficiencies.Get("btag_c_hadronFlavour_candidates_tagged")
            bEff_histo_udsg_num = the_updown_bSF_efficiencies.Get("btag_udsg_hadronFlavour_candidates_tagged")
            # the denominator
            bEff_histo_b_denom    = the_updown_bSF_efficiencies.Get("btag_b_hadronFlavour_candidates")
            bEff_histo_c_denom    = the_updown_bSF_efficiencies.Get("btag_c_hadronFlavour_candidates")
            bEff_histo_udsg_denom = the_updown_bSF_efficiencies.Get("btag_udsg_hadronFlavour_candidates")

            ## for tests, reverse
            #bEff_histo_b_denom    = the_updown_bSF_efficiencies.Get("btag_b_hadronFlavour_candidates_tagged")
            #bEff_histo_c_denom    = the_updown_bSF_efficiencies.Get("btag_c_hadronFlavour_candidates_tagged")
            #bEff_histo_udsg_denom = the_updown_bSF_efficiencies.Get("btag_udsg_hadronFlavour_candidates_tagged")
            ## the denominator
            #bEff_histo_b_num    = the_updown_bSF_efficiencies.Get("btag_b_hadronFlavour_candidates")
            #bEff_histo_c_num    = the_updown_bSF_efficiencies.Get("btag_c_hadronFlavour_candidates")
            #bEff_histo_udsg_num = the_updown_bSF_efficiencies.Get("btag_udsg_hadronFlavour_candidates")

            bEff_histo_b_num    .Divide(bEff_histo_b_denom   )
            bEff_histo_c_num    .Divide(bEff_histo_c_denom   )
            bEff_histo_udsg_num .Divide(bEff_histo_udsg_denom)

            # set them in the calculator
            support_btagging_sf.bEff_histo_b    = bEff_histo_b_num
            support_btagging_sf.bEff_histo_c    = bEff_histo_c_num
            support_btagging_sf.bEff_histo_udsg = bEff_histo_udsg_num

            logging.info("test b SF weight with eff histo: %s" % bEff_histo_b.GetName())
            logging.info("test b SF weight with eff histo: %s" % support_btagging_sf.bEff_histo_b.GetName())

            # test the change
            flavId = 5
            jet_weight_bSF_nom_tagged, sf_t, eff_t     = calc_btag_sf_weight(True,  flavId, jet_pt, jet_eta)
            jet_weight_bSF_nom_not_tagged, sf_n, eff_n = calc_btag_sf_weight(False, flavId, jet_pt, jet_eta)
            logging.info("test TT sys %s b SF weight for %d: %10.6f %10.6f,  %10.6f %10.6f" % (dtag_usual, flavId, jet_weight_bSF_nom_not_tagged, jet_weight_bSF_nom_not_tagged, sf_t, eff_t))
            flavId = 1
            jet_weight_bSF_nom_tagged, sf_t, eff_t     = calc_btag_sf_weight(True,  flavId, jet_pt, jet_eta)
            jet_weight_bSF_nom_not_tagged, sf_n, eff_n = calc_btag_sf_weight(False, flavId, jet_pt, jet_eta)
            logging.info("test TT sys %s b SF weight for %d: %10.6f %10.6f,  %10.6f %10.6f" % (dtag_usual, flavId, jet_weight_bSF_nom_not_tagged, jet_weight_bSF_nom_not_tagged, sf_t, eff_t))

    #print b_alljet, b_tagged, c_alljet, c_tagged, udsg_alljet, udsg_tagged
    #global bTagging_b_jet_efficiency, bTagging_c_jet_efficiency, bTagging_udsg_jet_efficiency
    #bTagging_b_jet_efficiency = bTagging_X_jet_efficiency(b_alljet, b_tagged)
    #bTagging_c_jet_efficiency = bTagging_X_jet_efficiency(c_alljet, c_tagged)
    #bTagging_udsg_jet_efficiency = bTagging_X_jet_efficiency(udsg_alljet, udsg_tagged)

    class ZeroOutDistr:
        def Fill(self, *args):
            pass
        def SetDirectory(self, *args):
            pass
        def Write(self, *args):
            pass
    #zero_out_distr = ZeroOutDistr()

    # helper codes for phi distrs
    pdgId_codes = {11: 'elP_phi', -11: 'elN_phi', 13: 'muP_phi', -13: 'muN_phi'}


    # >>>>>>>>>>>>>> make output tree of the stage2 selection

    # lor vector for selected objects (particles and met)
    LorentzVector_Class = LorentzVector('ROOT::Math::PxPyPzE4D<double>')
    # vector of these for particles
    ROOT.gROOT.ProcessLine("typedef std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >> LorentzVectorS;")
    ROOT.gROOT.ProcessLine("typedef std::vector<int> IntVector;")
    ROOT.gROOT.ProcessLine("typedef std::vector<double> DoubleVector;")

    all_vector_branches = [] # for resetting them on each event

    #ttree_out = TTree( 'ttree_out', 'tree with stage2 selection' ) # INPUT NOW
    indexevents  = array( 'L', [ 0 ] )
    ttree_out.Branch( 'indexevents', indexevents, 'indexevents/l' )
    runNumber  = array( 'i', [ 0 ] )
    ttree_out.Branch( 'runNumber', runNumber, 'runNumber/I' )
    lumi  = array( 'i', [ 0 ] )
    ttree_out.Branch( 'lumi', lumi, 'lumi/I' )

    # ID of the process in this particular MC
    gen_proc_id = array( 'i', [ 0 ] )   # taken from the standard example, array is needed to pass the pointer to branch
    ttree_out.Branch( 'gen_proc_id', gen_proc_id, 'gen_proc_id/I' )

    # an event can pass multiple selections
    # I need to save some code which selections are passed
    # here I use the simplest appropriate method
    # usually I save linear sequences of selctions, tighter at each stage
    # therefore selection stage would be apropriate
    selection_stage = array( 'i', [ 0 ])
    ttree_out.Branch('selection_stage', selection_stage, 'selection_stage/I')

    selection_stage_presel = array( 'i', [ 0 ])
    ttree_out.Branch('selection_stage_presel', selection_stage_presel, 'selection_stage_presel/I')

    selection_stage_TESUp = array( 'i', [ 0 ])
    ttree_out.Branch('selection_stage_TESUp', selection_stage_TESUp, 'selection_stage_TESUp/I')
    selection_stage_TESDown = array( 'i', [ 0 ])
    ttree_out.Branch('selection_stage_TESDown', selection_stage_TESDown, 'selection_stage_TESDown/I')

    selection_stage_JESUp = array( 'i', [ 0 ])
    ttree_out.Branch('selection_stage_JESUp', selection_stage_JESUp, 'selection_stage_JESUp/I')
    selection_stage_JESDown = array( 'i', [ 0 ])
    ttree_out.Branch('selection_stage_JESDown', selection_stage_JESDown, 'selection_stage_JESDown/I')

    selection_stage_JERUp = array( 'i', [ 0 ])
    ttree_out.Branch('selection_stage_JERUp', selection_stage_JERUp, 'selection_stage_JERUp/I')
    selection_stage_JERDown = array( 'i', [ 0 ])
    ttree_out.Branch('selection_stage_JERDown', selection_stage_JERDown, 'selection_stage_JERDown/I')

    selection_stage_tt_alliso = array( 'i', [ 0 ])
    ttree_out.Branch('selection_stage_tt_alliso', selection_stage_tt_alliso, 'selection_stage_tt_alliso/I')

    selection_stage_wjets = array( 'i', [ 0 ])
    ttree_out.Branch('selection_stage_wjets', selection_stage_wjets, 'selection_stage_wjets/I')

    selection_stage_dy         = array('i', [0])
    selection_stage_dy_TESUp   = array('i', [0])
    selection_stage_dy_TESDown = array('i', [0])
    selection_stage_dy_JERUp   = array('i', [0])
    selection_stage_dy_JERDown = array('i', [0])
    selection_stage_dy_JESUp   = array('i', [0])
    selection_stage_dy_JESDown = array('i', [0])
    ttree_out.Branch('selection_stage_dy',         selection_stage_dy,         'selection_stage_dy/I' )
    ttree_out.Branch('selection_stage_dy_TESUp',   selection_stage_dy_TESUp  , 'selection_stage_dy_TESUp/I' )
    ttree_out.Branch('selection_stage_dy_TESDown', selection_stage_dy_TESDown, 'selection_stage_dy_TESDown/I' )
    ttree_out.Branch('selection_stage_dy_JERUp',   selection_stage_dy_JERUp  , 'selection_stage_dy_JERUp/I' )
    ttree_out.Branch('selection_stage_dy_JERDown', selection_stage_dy_JERDown, 'selection_stage_dy_JERDown/I' )
    ttree_out.Branch('selection_stage_dy_JESUp',   selection_stage_dy_JESUp  , 'selection_stage_dy_JESUp/I' )
    ttree_out.Branch('selection_stage_dy_JESDown', selection_stage_dy_JESDown, 'selection_stage_dy_JESDown/I' )

    selection_stage_dy_mumu         = array('i', [0])
    selection_stage_dy_mumu_TESUp   = array('i', [0])
    selection_stage_dy_mumu_TESDown = array('i', [0])
    selection_stage_dy_mumu_JERUp   = array('i', [0])
    selection_stage_dy_mumu_JERDown = array('i', [0])
    selection_stage_dy_mumu_JESUp   = array('i', [0])
    selection_stage_dy_mumu_JESDown = array('i', [0])
    ttree_out.Branch('selection_stage_dy_mumu',         selection_stage_dy_mumu,         'selection_stage_dy_mumu/I' )
    ttree_out.Branch('selection_stage_dy_mumu_TESUp',   selection_stage_dy_mumu_TESUp  , 'selection_stage_dy_mumu_TESUp/I' )
    ttree_out.Branch('selection_stage_dy_mumu_TESDown', selection_stage_dy_mumu_TESDown, 'selection_stage_dy_mumu_TESDown/I' )
    ttree_out.Branch('selection_stage_dy_mumu_JERUp',   selection_stage_dy_mumu_JERUp  , 'selection_stage_dy_mumu_JERUp/I' )
    ttree_out.Branch('selection_stage_dy_mumu_JERDown', selection_stage_dy_mumu_JERDown, 'selection_stage_dy_mumu_JERDown/I' )
    ttree_out.Branch('selection_stage_dy_mumu_JESUp',   selection_stage_dy_mumu_JESUp  , 'selection_stage_dy_mumu_JESUp/I' )
    ttree_out.Branch('selection_stage_dy_mumu_JESDown', selection_stage_dy_mumu_JESDown, 'selection_stage_dy_mumu_JESDown/I' )

    selection_stage_em         = array('i', [0])
    selection_stage_em_TESUp   = array('i', [0])
    selection_stage_em_TESDown = array('i', [0])
    selection_stage_em_JERUp   = array('i', [0])
    selection_stage_em_JERDown = array('i', [0])
    selection_stage_em_JESUp   = array('i', [0])
    selection_stage_em_JESDown = array('i', [0])
    ttree_out.Branch( 'selection_stage_em',         selection_stage_em,         'selection_stage_em/I' )
    ttree_out.Branch( 'selection_stage_em_TESUp',   selection_stage_em_TESUp  , 'selection_stage_em_TESUp/I' )
    ttree_out.Branch( 'selection_stage_em_TESDown', selection_stage_em_TESDown, 'selection_stage_em_TESDown/I' )
    ttree_out.Branch( 'selection_stage_em_JERUp',   selection_stage_em_JERUp  , 'selection_stage_em_JERUp/I' )
    ttree_out.Branch( 'selection_stage_em_JERDown', selection_stage_em_JERDown, 'selection_stage_em_JERDown/I' )
    ttree_out.Branch( 'selection_stage_em_JESUp',   selection_stage_em_JESUp  , 'selection_stage_em_JESUp/I' )
    ttree_out.Branch( 'selection_stage_em_JESDown', selection_stage_em_JESDown, 'selection_stage_em_JESDown/I' )

    nup = array( 'i', [ 0 ] )
    ttree_out.Branch( 'nup', nup, 'nup/I' )

    # lj_var of the jets in the event
    event_jets_lj_var = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_jets_lj_var', event_jets_lj_var, 'event_jets_lj_var/f' )
    event_jets_input_has = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_jets_input_has', event_jets_input_has, 'event_jets_input_has/f' )
    event_jets_found_has = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_jets_found_has', event_jets_found_has, 'event_jets_found_has/f' )

    event_jets_n_jets = array( 'i', [ 0 ] )
    ttree_out.Branch( 'event_jets_n_jets', event_jets_n_jets, 'event_jets_n_jets/I' )
    event_jets_n_bjets = array( 'i', [ 0 ] )
    ttree_out.Branch( 'event_jets_n_bjets', event_jets_n_bjets, 'event_jets_n_bjets/I' )
    event_jets_n_jets_lepmatched = array( 'i', [ 0 ] )
    ttree_out.Branch( 'event_jets_n_jets_lepmatched', event_jets_n_jets_lepmatched, 'event_jets_n_jets_lepmatched/I' )
    #jets_nom.lepmatched.append((p4, en_factors, bSF_weights, jet_b_discr, HF, PF, jet_index))

    event_jets_alliso_n_jets = array( 'i', [ 0 ] )
    ttree_out.Branch( 'event_jets_alliso_n_jets', event_jets_alliso_n_jets, 'event_jets_alliso_n_jets/I' )
    event_jets_alliso_n_bjets = array( 'i', [ 0 ] )
    ttree_out.Branch( 'event_jets_alliso_n_bjets', event_jets_alliso_n_bjets, 'event_jets_alliso_n_bjets/I' )

    # nominal weight of the event
    event_weight = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight', event_weight, 'event_weight/f' )

    event_weight_init = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_init', event_weight_init, 'event_weight_init/f' )

    event_weight_PU = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PU', event_weight_PU, 'event_weight_PU/f' )
    event_weight_PU_el = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PU_el', event_weight_PU_el, 'event_weight_PU_el/f' )
    event_weight_PU_mu = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PU_mu', event_weight_PU_mu, 'event_weight_PU_mu/f' )
    event_weight_PU_bcdef = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PU_bcdef', event_weight_PU_bcdef, 'event_weight_PU_bcdef/f' )
    event_weight_PU_gh = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PU_gh', event_weight_PU_gh, 'event_weight_PU_gh/f' )
    event_weight_PU_per_epoch = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PU_per_epoch', event_weight_PU_per_epoch, 'event_weight_PU_per_epoch/f' )

    event_weight_th = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_th', event_weight_th, 'event_weight_th/f' )
    event_weight_bSF = array( 'f', [ 0 ] )

    ttree_out.Branch( 'event_weight_bSF', event_weight_bSF, 'event_weight_bSF/f' )
    event_weight_bSFUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_bSFUp', event_weight_bSFUp, 'event_weight_bSFUp/f' )
    event_weight_bSFDown = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_bSFDown', event_weight_bSFDown, 'event_weight_bSFDown/f' )
    event_weight_bSF_JERUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_bSF_JERUp', event_weight_bSF_JERUp, 'event_weight_bSF_JERUp/f' )
    event_weight_bSF_JERDown = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_bSF_JERDown', event_weight_bSF_JERDown, 'event_weight_bSF_JERDown/f' )
    event_weight_bSF_JESUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_bSF_JESUp', event_weight_bSF_JESUp, 'event_weight_bSF_JESUp/f' )
    event_weight_bSF_JESDown = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_bSF_JESDown', event_weight_bSF_JESDown, 'event_weight_bSF_JESDown/f' )

    event_weight_toppt = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_toppt', event_weight_toppt, 'event_weight_toppt/f' )

    event_weight_LEPall = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPall', event_weight_LEPall, 'event_weight_LEPall/f' )
    event_weight_LEP_PU = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEP_PU', event_weight_LEP_PU, 'event_weight_LEP_PU/f' )

    event_weight_LEPmuBID = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmuBID', event_weight_LEPmuBID, 'event_weight_LEPmuBID/f' )
    event_weight_LEPmuBTRG = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmuBTRG', event_weight_LEPmuBTRG, 'event_weight_LEPmuBTRG/f' )

    event_weight_LEPmuHID = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmuHID', event_weight_LEPmuHID, 'event_weight_LEPmuHID/f' )
    event_weight_LEPmuHTRG = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmuHTRG', event_weight_LEPmuHTRG, 'event_weight_LEPmuHTRG/f' )

    event_weight_LEPmu0ID = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmu0ID', event_weight_LEPmu0ID, 'event_weight_LEPmu0ID/f' )
    event_weight_LEPmu0TRG = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmu0TRG', event_weight_LEPmu0TRG, 'event_weight_LEPmu0TRG/f' )

    event_weight_LEPmu0TRGUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmu0TRGUp', event_weight_LEPmu0TRGUp, 'event_weight_LEPmu0TRGUp/f' )
    event_weight_LEPmuTRGctrlUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmuTRGctrlUp', event_weight_LEPmuTRGctrlUp, 'event_weight_LEPmuTRGctrlUp/f' )
    event_weight_LEPmuTRGctrlDown = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmuTRGctrlDown', event_weight_LEPmuTRGctrlDown, 'event_weight_LEPmuTRGctrlDown/f' )

    event_weight_LEPelID = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPelID', event_weight_LEPelID, 'event_weight_LEPelID/f' )
    event_weight_LEPelTRG = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPelTRG', event_weight_LEPelTRG, 'event_weight_LEPelTRG/f' )

    event_weight_LEPmuIDUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmuIDUp', event_weight_LEPmuIDUp, 'event_weight_LEPmuIDUp/f' )
    event_weight_LEPmuIDDown = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmuIDDown', event_weight_LEPmuIDDown, 'event_weight_LEPmuIDDown/f' )
    event_weight_LEPmuTRGUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmuTRGUp', event_weight_LEPmuTRGUp, 'event_weight_LEPmuTRGUp/f' )
    event_weight_LEPmuTRGDown = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPmuTRGDown', event_weight_LEPmuTRGDown, 'event_weight_LEPmuTRGDown/f' )

    event_weight_LEPelIDUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPelIDUp', event_weight_LEPelIDUp, 'event_weight_LEPelIDUp/f' )
    event_weight_LEPelIDDown = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPelIDDown', event_weight_LEPelIDDown, 'event_weight_LEPelIDDown/f' )
    event_weight_LEPelTRGUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPelTRGUp', event_weight_LEPelTRGUp, 'event_weight_LEPelTRGUp/f' )
    event_weight_LEPelTRGDown = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_LEPelTRGDown', event_weight_LEPelTRGDown, 'event_weight_LEPelTRGDown/f' )

    event_weight_PUUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PUUp', event_weight_PUUp, 'event_weight_PUUp/f' )
    event_weight_PUDown = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PUDown', event_weight_PUDown, 'event_weight_PUDown/f' )

    event_weight_PUUp_el = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PUUp_el', event_weight_PUUp_el, 'event_weight_PUUp_el/f' )
    event_weight_PUDown_el = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PUDown_el', event_weight_PUDown_el, 'event_weight_PUDown_el/f' )

    event_weight_PUUp_mu = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PUUp_mu', event_weight_PUUp_mu, 'event_weight_PUUp_mu/f' )
    event_weight_PUDown_mu = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PUDown_mu', event_weight_PUDown_mu, 'event_weight_PUDown_mu/f' )

    event_weight_PUUp_bcdef = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PUUp_bcdef', event_weight_PUUp_bcdef, 'event_weight_PUUp_bcdef/f' )
    event_weight_PUDown_bcdef = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PUDown_bcdef', event_weight_PUDown_bcdef, 'event_weight_PUDown_bcdef/f' )

    event_weight_PUUp_gh = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PUUp_gh', event_weight_PUUp_gh, 'event_weight_PUUp_gh/f' )
    event_weight_PUDown_gh = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PUDown_gh', event_weight_PUDown_gh, 'event_weight_PUDown_gh/f' )

    event_weight_PetersonUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_PetersonUp', event_weight_PetersonUp, 'event_weight_PetersonUp/f' )
    event_weight_FragUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_FragUp', event_weight_FragUp, 'event_weight_FragUp/f' )
    event_weight_FragDown = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_FragDown', event_weight_FragDown, 'event_weight_FragDown/f' )

    event_weight_SemilepBRUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_SemilepBRUp', event_weight_SemilepBRUp, 'event_weight_SemilepBRUp/f' )
    event_weight_SemilepBRDown = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_SemilepBRDown', event_weight_SemilepBRDown, 'event_weight_SemilepBRDown/f' )


    event_weight_me_f_rUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_me_f_rUp', event_weight_me_f_rUp, 'event_weight_me_f_rUp/f' )
    event_weight_me_f_rDn = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_me_f_rDn', event_weight_me_f_rDn, 'event_weight_me_f_rDn/f' )

    event_weight_me_fUp_r = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_me_fUp_r', event_weight_me_fUp_r, 'event_weight_me_fUp_r/f' )
    event_weight_me_fDn_r = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_me_fDn_r', event_weight_me_fDn_r, 'event_weight_me_fDn_r/f' )

    event_weight_me_frUp = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_me_frUp', event_weight_me_frUp, 'event_weight_me_frUp/f' )
    event_weight_me_frDn = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_me_frDn', event_weight_me_frDn, 'event_weight_me_frDn/f' )


    event_weight_AlphaS_up = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_AlphaS_up', event_weight_AlphaS_up, 'event_weight_AlphaS_up/f' )
    event_weight_AlphaS_dn = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_weight_AlphaS_dn', event_weight_AlphaS_dn, 'event_weight_AlphaS_dn/f' )


    event_weight_pdf = ROOT.DoubleVector()
    ttree_out.Branch("event_weight_pdf", event_weight_pdf)
    all_vector_branches.append(event_weight_pdf)


    amcatnlo_w = array( 'f', [ 0 ] )
    ttree_out.Branch( 'amcatnlo_w', amcatnlo_w, 'amcatnlo_w/f' )

    # met lorentzvector
    event_met = LorentzVector_Class(0., 0., 0., 0.)
    ttree_out.Branch("event_met", event_met)
    event_met_JERUp = LorentzVector_Class(0., 0., 0., 0.)
    ttree_out.Branch("event_met_JERUp", event_met_JERUp)
    event_met_JERDown = LorentzVector_Class(0., 0., 0., 0.)
    ttree_out.Branch("event_met_JERDown", event_met_JERDown)
    event_met_JESUp = LorentzVector_Class(0., 0., 0., 0.)
    ttree_out.Branch("event_met_JESUp", event_met_JESUp)
    event_met_JESDown = LorentzVector_Class(0., 0., 0., 0.)
    ttree_out.Branch("event_met_JESDown", event_met_JESDown)
    event_met_TESUp = LorentzVector_Class(0., 0., 0., 0.)
    ttree_out.Branch("event_met_TESUp", event_met_TESUp)
    event_met_TESDown = LorentzVector_Class(0., 0., 0., 0.)
    ttree_out.Branch("event_met_TESDown", event_met_TESDown)

    event_met_init  = LorentzVector_Class(0., 0., 0., 0.)
    ttree_out.Branch("event_met_init", event_met_init)
    event_met_init2 = LorentzVector_Class(0., 0., 0., 0.)
    ttree_out.Branch("event_met_init2", event_met_init2)
    event_met_lep_mt_init  = array('f', [0])
    ttree_out.Branch('event_met_lep_mt_init', event_met_lep_mt_init, 'event_met_lep_mt_init/f')
    event_met_lep_mt_init2 = array('f', [0])
    ttree_out.Branch('event_met_lep_mt_init2', event_met_lep_mt_init2, 'event_met_lep_mt_init2/f')

    event_met_lep_mt = array('f', [0])
    ttree_out.Branch('event_met_lep_mt', event_met_lep_mt, 'event_met_lep_mt/f')
    event_met_lep_mt_JERUp   = array('f', [0])
    ttree_out.Branch('event_met_lep_mt_JERUp',   event_met_lep_mt_JERUp,   'event_met_lep_mt_JERUp/f')
    event_met_lep_mt_JERDown = array('f', [0])
    ttree_out.Branch('event_met_lep_mt_JERDown', event_met_lep_mt_JERDown, 'event_met_lep_mt_JERDown/f')
    event_met_lep_mt_JESUp   = array('f', [0])
    ttree_out.Branch('event_met_lep_mt_JESUp',   event_met_lep_mt_JESUp,   'event_met_lep_mt_JESUp/f')
    event_met_lep_mt_JESDown = array('f', [0])
    ttree_out.Branch('event_met_lep_mt_JESDown', event_met_lep_mt_JESDown, 'event_met_lep_mt_JESDown/f')
    event_met_lep_mt_TESUp   = array('f', [0])
    ttree_out.Branch('event_met_lep_mt_TESUp',   event_met_lep_mt_TESUp,   'event_met_lep_mt_TESUp/f')
    event_met_lep_mt_TESDown = array('f', [0])
    ttree_out.Branch('event_met_lep_mt_TESDown', event_met_lep_mt_TESDown, 'event_met_lep_mt_TESDown/f')

    event_dilep_mass = array('f', [0])
    ttree_out.Branch('event_dilep_mass', event_dilep_mass, 'event_dilep_mass/f')

    event_top_masses_medium = ROOT.DoubleVector()
    ttree_out.Branch("event_top_masses_medium", event_top_masses_medium)
    all_vector_branches.append(event_top_masses_medium)
    event_top_masses_loose = ROOT.DoubleVector()
    ttree_out.Branch("event_top_masses_loose", event_top_masses_loose)
    all_vector_branches.append(event_top_masses_loose)

    event_top_masses_closest = array( 'f', [ 0 ] )
    ttree_out.Branch( 'event_top_masses_closest', event_top_masses_closest, 'event_top_masses_closest/f' )

    event_leptons = ROOT.LorentzVectorS()
    ttree_out.Branch("event_leptons", event_leptons)
    all_vector_branches.append(event_leptons)
    event_leptons_ids = ROOT.IntVector()
    ttree_out.Branch("event_leptons_ids", event_leptons_ids)
    all_vector_branches.append(event_leptons_ids)

    # lep_dB -- 3D impact, there is also dxy and dz, from tau SV experience dz is very not precise and needs smart corrections
    # therefore store dxy
    event_leptons_dxy = ROOT.DoubleVector()
    ttree_out.Branch("event_leptons_dxy", event_leptons_dxy)
    all_vector_branches.append(event_leptons_dxy)
    # lep_matching_gen
    event_leptons_genmatch = ROOT.IntVector()
    ttree_out.Branch("event_leptons_genmatch", event_leptons_genmatch)
    all_vector_branches.append(event_leptons_genmatch)

    event_leptons_alliso_reliso = ROOT.DoubleVector()
    ttree_out.Branch("event_leptons_alliso_reliso", event_leptons_alliso_reliso)
    all_vector_branches.append(event_leptons_alliso_reliso)
    event_leptons_alliso_p4 = ROOT.LorentzVectorS()
    ttree_out.Branch("event_leptons_alliso_p4", event_leptons_alliso_p4)
    all_vector_branches.append(event_leptons_alliso_p4)
    event_leptons_alliso_pdgId = ROOT.IntVector()
    ttree_out.Branch("event_leptons_alliso_pdgId", event_leptons_alliso_pdgId)
    all_vector_branches.append(event_leptons_alliso_pdgId)

    event_taus_alliso_p4 = ROOT.LorentzVectorS()
    ttree_out.Branch("event_taus_alliso_p4", event_taus_alliso_p4)
    all_vector_branches.append(event_taus_alliso_p4)
    event_taus_alliso_pdgId = ROOT.IntVector()
    ttree_out.Branch("event_taus_alliso_pdgId", event_taus_alliso_pdgId)
    all_vector_branches.append(event_taus_alliso_pdgId)
    event_taus_alliso_IDlev = ROOT.IntVector()
    ttree_out.Branch("event_taus_alliso_IDlev", event_taus_alliso_IDlev)
    all_vector_branches.append(event_taus_alliso_IDlev)

    event_candidate_taus = ROOT.LorentzVectorS()
    ttree_out.Branch("event_candidate_taus", event_candidate_taus)
    all_vector_branches.append(event_candidate_taus)
    event_candidate_taus_ids = ROOT.IntVector()
    ttree_out.Branch("event_candidate_taus_ids", event_candidate_taus_ids)
    all_vector_branches.append(event_candidate_taus_ids)

    event_taus = ROOT.LorentzVectorS()
    ttree_out.Branch("event_taus", event_taus)
    all_vector_branches.append(event_taus)
    event_taus_ids = ROOT.IntVector()
    ttree_out.Branch("event_taus_ids", event_taus_ids)
    all_vector_branches.append(event_taus_ids)
    event_taus_IDlev = ROOT.IntVector()
    ttree_out.Branch("event_taus_IDlev", event_taus_IDlev)
    all_vector_branches.append(event_taus_IDlev)

    event_taus_TES_up = ROOT.DoubleVector()
    ttree_out.Branch("event_taus_TES_up", event_taus_TES_up)
    all_vector_branches.append(event_taus_TES_up)
    event_taus_TES_down = ROOT.DoubleVector()
    ttree_out.Branch("event_taus_TES_down", event_taus_TES_down)
    all_vector_branches.append(event_taus_TES_down)

    # tau_matching_gen
    event_taus_genmatch = ROOT.IntVector()
    ttree_out.Branch("event_taus_genmatch", event_taus_genmatch)
    all_vector_branches.append(event_taus_genmatch)

    # 3ch info
    event_taus_pat_sv_sign = ROOT.DoubleVector()
    ttree_out.Branch("event_taus_pat_sv_sign", event_taus_pat_sv_sign)
    all_vector_branches.append(event_taus_pat_sv_sign)
    # sign
    event_taus_sv_sign = ROOT.DoubleVector()
    ttree_out.Branch("event_taus_sv_sign", event_taus_sv_sign)
    all_vector_branches.append(event_taus_sv_sign)
    # length
    event_taus_sv_leng = ROOT.DoubleVector()
    ttree_out.Branch("event_taus_sv_leng", event_taus_sv_leng)
    all_vector_branches.append(event_taus_sv_leng)
    # dalitz m1
    event_taus_sv_dalitz_m1 = ROOT.DoubleVector()
    ttree_out.Branch("event_taus_sv_dalitz_m1", event_taus_sv_dalitz_m1)
    all_vector_branches.append(event_taus_sv_dalitz_m1)
    # dalitz m2
    event_taus_sv_dalitz_m2 = ROOT.DoubleVector()
    ttree_out.Branch("event_taus_sv_dalitz_m2", event_taus_sv_dalitz_m2)
    all_vector_branches.append(event_taus_sv_dalitz_m2)
    # sum of track energies
    event_taus_track_energy = ROOT.DoubleVector()
    ttree_out.Branch("event_taus_track_energy", event_taus_track_energy)
    all_vector_branches.append(event_taus_track_energy)
    # b discriminant of matching jet
    event_taus_jet_bdiscr = ROOT.DoubleVector()
    ttree_out.Branch("event_taus_jet_bdiscr", event_taus_jet_bdiscr)
    all_vector_branches.append(event_taus_jet_bdiscr)

    event_taus_l = ROOT.LorentzVectorS()
    ttree_out.Branch("event_taus_l", event_taus_l)
    all_vector_branches.append(event_taus_l)


    event_jets_b = ROOT.LorentzVectorS()
    ttree_out.Branch("event_jets_b", event_jets_b)
    all_vector_branches.append(event_jets_b)
    # b discriminant the jet
    event_jets_b_bdiscr = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_b_bdiscr", event_jets_b_bdiscr)
    all_vector_branches.append(event_jets_b_bdiscr)

    event_jets_JERUp = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_JERUp", event_jets_JERUp)
    all_vector_branches.append(event_jets_JERUp)
    event_jets_JERDown = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_JERDown", event_jets_JERDown)
    all_vector_branches.append(event_jets_JERDown)

    event_jets_JESUp = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_JESUp", event_jets_JESUp)
    all_vector_branches.append(event_jets_JESUp)
    event_jets_JESDown = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_JESDown", event_jets_JESDown)
    all_vector_branches.append(event_jets_JESDown)

    # genmatch
    event_jets_b_genmatch = ROOT.IntVector()
    ttree_out.Branch("event_jets_b_genmatch", event_jets_b_genmatch)
    all_vector_branches.append(event_jets_b_genmatch)
    # b SF weight
    event_jets_b_bSFweight = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_b_bSFweight", event_jets_b_bSFweight)
    all_vector_branches.append(event_jets_b_bSFweight)

    event_jets_r = ROOT.LorentzVectorS()
    ttree_out.Branch("event_jets_r", event_jets_r)
    all_vector_branches.append(event_jets_r)
    # b discriminant the jet
    event_jets_r_bdiscr = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_r_bdiscr", event_jets_r_bdiscr)
    all_vector_branches.append(event_jets_r_bdiscr)
    # genmatch
    event_jets_r_genmatch = ROOT.IntVector()
    ttree_out.Branch("event_jets_r_genmatch", event_jets_r_genmatch)
    all_vector_branches.append(event_jets_r_genmatch)
    # b SF weight
    event_jets_r_bSFweight = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_r_bSFweight", event_jets_r_bSFweight)
    all_vector_branches.append(event_jets_r_bSFweight)

    event_jets_l = ROOT.LorentzVectorS()
    ttree_out.Branch("event_jets_l", event_jets_l)
    all_vector_branches.append(event_jets_l)
    # b discriminant the jet
    event_jets_l_bdiscr = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_l_bdiscr", event_jets_l_bdiscr)
    all_vector_branches.append(event_jets_l_bdiscr)
    # genmatch
    event_jets_l_genmatch = ROOT.IntVector()
    ttree_out.Branch("event_jets_l_genmatch", event_jets_l_genmatch)
    all_vector_branches.append(event_jets_l_genmatch)
    # b SF weight
    event_jets_l_bSFweight = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_l_bSFweight", event_jets_l_bSFweight)
    all_vector_branches.append(event_jets_l_bSFweight)

    event_jets_t = ROOT.LorentzVectorS()
    ttree_out.Branch("event_jets_t", event_jets_t)
    all_vector_branches.append(event_jets_t)
    # b discriminant the jet
    event_jets_t_bdiscr = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_t_bdiscr", event_jets_t_bdiscr)
    all_vector_branches.append(event_jets_t_bdiscr)
    # genmatch
    event_jets_t_genmatch = ROOT.IntVector()
    ttree_out.Branch("event_jets_t_genmatch", event_jets_t_genmatch)
    all_vector_branches.append(event_jets_t_genmatch)
    # b SF weight
    event_jets_t_bSFweight = ROOT.DoubleVector()
    ttree_out.Branch("event_jets_t_bSFweight", event_jets_t_bSFweight)
    all_vector_branches.append(event_jets_t_bSFweight)

    # <<<<<<<<<<<<<<

    logger.write("N entries: %d\n" % tree.GetEntries())
    #if not range_max or range_max > tree.GetEntries():
        #range_mas = tree.GetEntries()

    #pdb.set_trace()

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

        #if iev > 10000: break
        control_counters.Fill(0)

        for vector_branch in all_vector_branches:
            vector_branch.clear()

        # NUP stitching for WJets
        #if isWJetsInclusive and DO_W_STITCHING:
        #    if ev.gen_NUP > 5: continue
        nup[0] = ev.gen_NUP

        #if iev <  range_min: continue
        #if iev >= 10: break
        #ev = tree.GetEntry(iev)

        '''
        # plug roccor right here
        #roccor_factor = 1.
        for i in range(len(ev.lep_id)):
            if abs(ev.lep_id[i]) == 13:
                #Q = c_int(ev.lep_id > 0)
                #pt = c_double(ev.lep_p4[i].pt())
                #eta = c_double(ev.lep_p4[i].eta())
                #phi = c_double(ev.lep_p4[i].phi())
                ##nl = c_int(ev.lep_N_trackerLayersWithMeasurement[i])
                #nl = c_int(14)
                #genPt = pt
                #u1 = c_double(0.5)
                #u2 = c_double(0.5)
                #s = c_int(0)
                #m = c_int(0)

                #
                Q = int(ev.lep_id > 0)
                pt  = ev.lep_p4[i].pt()
                eta = ev.lep_p4[i].eta()
                phi = ev.lep_p4[i].phi()
                #nl = ev.lep_N_trackerLayersWithMeasurement[i]
                nl = 12 # average in data and mc
                genPt = pt # TODO add this
                u1 = u_random() #0.5
                u2 = u_random() #0.5
                s = 0
                m = 0

                #if isMC and ev.lep_matching_gen[i]:
                #    #roccor_factor = roccors.wrapper_kScaleFromGenMC(Q, pt, eta, phi, nl, genPt, u1, s, m)
                #    roccor_factor = roccors.kScaleFromGenMC(Q, pt, eta, phi, nl, genPt, u1, s, m)
                #elif isMC and not ev.lep_matching_gen[i]:
                #    #roccor_factor = roccors.wrapper_kScaleAndSmearMC(Q, pt, eta, phi, nl, u1, u2, s, m)
                #    roccor_factor = roccors.kScaleAndSmearMC(Q, pt, eta, phi, nl, u1, u2, s, m)
                #else:
                #    #roccor_factor = roccors.wrapper_kScaleDT(Q, pt, eta, phi, s, m)
                #    roccor_factor = roccors.kScaleDT(Q, pt, eta, phi, s, m)

                ##ev.lep_p4[i] *= roccor_factor
                #control_hs['roccor_factor'].Fill(roccor_factor)
        '''

        if not isMC:
            if not (ev.pass_basic_METfilters and ev.METfilterbadChCand and ev.METfilterbadPFMuon):
                continue

        # the lepton requirements for all 1-lepton channels:
        # do relIso on 1/8 = 0.125, and "all iso" for QCD anti-iso factor

        # I'll make the iso distribution and get the factor over whole range
        pass_mu_id = abs(ev.leps_ID) == 13 and ev.HLT_mu and ev.lep_matched_HLT[0] and ev.no_iso_veto_leps and abs(ev.lep_dxy[0]) < 0.01 and abs(ev.lep_dz[0]) < 0.02
        pass_el_id = abs(ev.leps_ID) == 11 and ev.HLT_el and ev.lep_matched_HLT[0] and ev.no_iso_veto_leps

        if abs(ev.leps_ID) == 13:    control_counters.Fill(1)
        if abs(ev.leps_ID) == 13 and ev.HLT_mu:  control_counters.Fill(2)
        if abs(ev.leps_ID) == 13 and ev.HLT_mu and ev.lep_matched_HLT[0]:  control_counters.Fill(3)
        if abs(ev.leps_ID) == 13 and ev.HLT_mu and ev.lep_matched_HLT[0] and ev.no_iso_veto_leps:  control_counters.Fill(4)

        if abs(ev.leps_ID) == 11:    control_counters.Fill(6)
        if abs(ev.leps_ID) == 11 and ev.HLT_el:  control_counters.Fill(7)
        if abs(ev.leps_ID) == 11 and ev.HLT_el and ev.lep_matched_HLT[0]:  control_counters.Fill(8)
        if abs(ev.leps_ID) == 11 and ev.HLT_el and ev.lep_matched_HLT[0] and ev.no_iso_veto_leps:  control_counters.Fill(9)

        # impacts are embedded into ID procedure, no need to repeate them here, TODO: re-check
        # for mu suggested dB is 5mm dz and 2mm dxy
        # for el suggested dB 0.02 cm
        #pass_mu_impact = ev.lep_dz[0] < 0.5 and ev.lep_dxy[0] < 0.2
        #pass_el_impact = ev.lep_dB[0] < 0.02

        # for el the suggested relIso 0.0588 for barrel, 0.0571 for endcaps
        # data shows it does not matter
        pass_mu_iso = pass_mu_id and ev.lep_relIso[0] < 0.15  
        pass_el_iso = pass_el_id and ev.lep_relIso[0] < 0.0588

        # 1GeV above HLT pt
        # ele eta gap
        # suggested minimal offline
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#2016_Data
        # pt > 25, eta < 2.4
        pass_mu_kino = pass_mu_id and ev.lep_p4[0].pt() > 26. and abs(ev.lep_p4[0].eta()) < 2.4
        #pass_mu_kino = pass_mu_id and ev.lep_p4[0].pt() * roccor_factor > 26. and abs(ev.lep_p4[0].eta()) < 2.4
        pass_el_kino = pass_el_id and ev.lep_p4[0].pt() > 30. and abs(ev.lep_p4[0].eta()) < 2.4 and (abs(ev.lep_p4[0].eta()) < 1.4442 or abs(ev.lep_p4[0].eta()) > 1.5660)

        # (did) for optimization testing minimum pt cut --- review it after test results
        #    -- significant discrepancy at lowest lep pt bin -> UP 1 GeV from HLT and added a detailed distr of the trun on curve
        #    -- it seems the discrepancy were comming from trigger SF going down to 26 only -- fixed that, testing now

        pass_mu     = pass_mu_id and pass_mu_kino and pass_mu_iso # and pass_mu_impact
        pass_el     = pass_el_id and pass_el_kino and pass_el_iso

        if pass_el_id:    control_counters.Fill(11)
        if pass_el_iso:   control_counters.Fill(12)
        if pass_el_kino:  control_counters.Fill(13)
        if pass_el:       control_counters.Fill(14)

        if pass_mu_id:    control_counters.Fill(16)
        if pass_mu_iso:   control_counters.Fill(17)
        if pass_mu_kino:  control_counters.Fill(18)
        if pass_mu:       control_counters.Fill(19)

        # do this in place for the required selection
        #pass_mu_cuts = pass_mu and ev.lep_p4[0].pt() > 25. # so, it's smaller than the old cut 27 GeV

        #pass_mu_all = pass_mu_id_kino and ev.lep_relIso[0] >= 0.125

        pass_mu_id_all = abs(ev.leps_ID_allIso) == 13 and ev.HLT_mu and ev.lep_alliso_matched_HLT[0] and ev.nleps_veto_mu_all == 0 and ev.nleps_veto_el_all == 0
        pass_el_id_all = abs(ev.leps_ID_allIso) == 11 and ev.HLT_el and ev.lep_alliso_matched_HLT[0] and ev.nleps_veto_el_all == 0 and ev.nleps_veto_mu_all == 0

        pass_mu_kino_all = pass_mu_id_all and ev.lep_alliso_p4[0].pt() > 26. and abs(ev.lep_alliso_p4[0].eta()) < 2.4
        pass_el_kino_all = pass_el_id_all and ev.lep_alliso_p4[0].pt() > 30. and abs(ev.lep_alliso_p4[0].eta()) < 2.4 and (abs(ev.lep_alliso_p4[0].eta()) < 1.4442 or abs(ev.lep_alliso_p4[0].eta()) > 1.5660)

        pass_mu_all = pass_mu_id_all and pass_mu_kino_all
        pass_el_all = pass_el_id_all and pass_el_kino_all

        #pass_mu = pass_mu_id_kino and ev.lep_relIso[0] < 0.125
        #pass_el_all = pass_el_id and ev.lep_p4[0].pt() > 27. and abs(ev.lep_p4[0].eta()) < 2.4
        #pass_el     = pass_el_id and ev.lep_p4[0].pt() > 27. and abs(ev.lep_p4[0].eta()) < 2.4 and ev.lep_relIso[0] < 0.125

        # TODO: need to check the trigger SF-s then......
        # for now I'm just trying to get rid of el-mu mix in MC
        #pass_elmu = ev.leps_ID == -11*13 and ev.HLT_mu and not ev.HLT_el and ev.no_iso_veto_leps and \
        pass_elmu_id = ev.leps_ID == -11*13 and ev.HLT_mu and ev.no_iso_veto_leps and \
            ((ev.lep_matched_HLT[0] and abs(ev.lep_dxy[0]) < 0.01 and abs(ev.lep_dz[0]) < 0.02) if abs(ev.lep_id[0]) == 13 else (ev.lep_matched_HLT[1] and abs(ev.lep_dxy[1]) < 0.01 and abs(ev.lep_dz[1]) < 0.02)) and \
            (ev.lep_p4[0].pt() > 30 and abs(ev.lep_p4[0].eta()) < 2.4) and \
            (ev.lep_p4[1].pt() > 30 and abs(ev.lep_p4[1].eta()) < 2.4)

        #   (ev.lep_matched_HLT[0] if abs(ev.lep_id[0]) == 13 else ev.lep_matched_HLT[1]) and \
        #pass_elmu = pass_elmu_id and ((ev.lep_relIso[0] < 0.15 and ev.lep_relIso[1] < 0.0588) if abs(ev.lep_id[0]) == 13 else (ev.lep_relIso[1] < 0.15 and ev.lep_relIso[0] < 0.0588))

        # these are done in NT
        pass_elmu = pass_elmu_id
        # both HLTs
        #pass_elmu = pass_elmu_id and ev.lep_matched_HLT[0] and ev.HLT_el and ev.lep_matched_HLT[1]

        pass_elmu_el = ev.leps_ID == -11*13 and ev.HLT_el and ev.no_iso_veto_leps and \
            (ev.lep_matched_HLT[0] if abs(ev.lep_id[0]) == 11 else ev.lep_matched_HLT[1]) and \
            abs(ev.lep_p4[0].eta()) < 2.4 and abs(ev.lep_p4[1].eta()) < 2.4 and \
            (ev.lep_p4[0].pt() > 32. and ev.lep_p4[1].pt() > 32.)

        #((ev.lep_p4[0].pt() > 30 and ev.lep_p4[1].pt() > 26) if abs(ev.lep_id[0]) == 11 else (ev.lep_p4[1].pt() > 30 and ev.lep_p4[0].pt() > 26))
        # no need to check eta crack -- it is selected in the Ntupler
        #((ev.lep_p4[0].pt() > 1.5 ) if abs(ev.lep_id[0]) == 11 else (ev.lep_p4[1].pt() > 30 and ev.lep_p4[0].pt() > 26))

        pass_mumu = ev.leps_ID == -13*13 and ev.HLT_mu and (ev.lep_matched_HLT[0] or ev.lep_matched_HLT[1]) and ev.no_iso_veto_leps and \
            (ev.lep_p4[0].pt() > 30 and abs(ev.lep_p4[0].eta()) < 2.4) and \
            (ev.lep_p4[1].pt() > 30 and abs(ev.lep_p4[1].eta()) < 2.4) and \
            abs(ev.lep_dxy[0]) < 0.01 and abs(ev.lep_dz[0]) < 0.02 and \
            abs(ev.lep_dxy[1]) < 0.01 and abs(ev.lep_dz[1]) < 0.02

        pass_elel = ev.leps_ID == -11*11 and ev.HLT_el and (ev.lep_matched_HLT[0] or ev.lep_matched_HLT[1]) and ev.no_iso_veto_leps and \
            (ev.lep_p4[0].pt() > 30 and abs(ev.lep_p4[0].eta()) < 2.4) and \
            (ev.lep_p4[1].pt() > 30 and abs(ev.lep_p4[1].eta()) < 2.4)

        pass_mumu_ss = ev.leps_ID == 13*13 and ev.HLT_mu and (ev.lep_matched_HLT[0] or ev.lep_matched_HLT[1]) and ev.no_iso_veto_leps and \
            (ev.lep_p4[0].pt() > 30 and abs(ev.lep_p4[0].eta()) < 2.4) and \
            (ev.lep_p4[1].pt() > 30 and abs(ev.lep_p4[1].eta()) < 2.4) and \
            abs(ev.lep_dxy[0]) < 0.01 and abs(ev.lep_dz[0]) < 0.02 and \
            abs(ev.lep_dxy[1]) < 0.01 and abs(ev.lep_dz[1]) < 0.02

        if pass_elmu:    control_counters.Fill(21)
        if pass_mumu:    control_counters.Fill(22)
        if pass_elel:    control_counters.Fill(23)
        if pass_el_all:  control_counters.Fill(24)
        if pass_mu_all:  control_counters.Fill(25)

        # minimum possible pt threshold -- 24 GeV, = to HLT and recorded in Ntupler
        pass_mu_id_iso = pass_mu_id and pass_mu_iso
        pass_min_mu     = pass_mu_id_iso and abs(ev.lep_p4[0].eta()) < 2.4 and ev.lep_relIso[0] < 0.125
        pass_min_mu_all = pass_mu_id_iso and abs(ev.lep_p4[0].eta()) < 2.4 and ev.lep_relIso[0] >= 0.125

        # 1-lep channels, 2mu DY and el-mu ttbar, and anti-iso
        #if not (pass_min_mu or pass_min_mu_all or pass_mu_all or pass_mu or pass_el_all or pass_el or pass_mumu or pass_mumu_ss or pass_elmu): continue
        #passes = pass_min_mu or pass_min_mu_all or pass_mu_all or pass_mu or pass_el or pass_mumu or pass_mumu_ss or pass_elmu or pass_mu_id_iso
        # OPTIMIZATION tests are done only on pass_mu
        #passes_optimized = pass_mu_all or pass_el_all or pass_mumu or pass_elmu
        #passes_optimized = pass_mu or pass_el or pass_mumu or pass_elmu or pass_mu_all or pass_el_all or pass_elel
        #event_passes = pass_mu # pass_elmu # or pass_elmu_el # pass_mumu or pass_elel # pass_el # or pass_elmu_el # pass_mu or pass_el # passes_optimized #
        passed_triggers = (pass_mu, pass_elmu, pass_elmu_el, pass_mumu, pass_elel, pass_el, pass_mu_all, pass_el_all)
        event_passes = PASSES_FUNC(*passed_triggers)

        #if pass_elmu or pass_elmu_el or pass_mumu or pass_elel:
        #    # check dR of leptons
        #    tlep0_p4 = TLorentzVector(ev.lep_p4[0].X(), ev.lep_p4[0].Y(), ev.lep_p4[0].Z(), ev.lep_p4[0].T())
        #    tlep1_p4 = TLorentzVector(ev.lep_p4[1].X(), ev.lep_p4[1].Y(), ev.lep_p4[1].Z(), ev.lep_p4[1].T())
        #    if not tlep0_p4.DeltaR(tlep1_p4) > 0.3:
        #        continue

        if not event_passes: continue
        control_counters.Fill(51)

        pass_mus = pass_mu_all or pass_mu or pass_elmu or pass_mumu # or pass_mumu_ss
        # also at least some kind of tau in single-el:
        #if (pass_mu or pass_el) and (not ev.tau_p4.size() > 0): continue # this is the only thing reduces computing
        passed_ele_event = pass_el or pass_el_all or pass_elel or pass_elmu_el


        # expensive calls and they don't depend on systematics now
        if doRecoilCorrections:
            #def transverse_mass_pts(v1_x, v1_y, v2_x, v2_y):
            #met_x = ev.pfmetcorr_ex
            #met_y = ev.pfmetcorr_ey
            # no, recalculated them

            met_x = ROOT.met_pt_recoilcor_x(
                #ev.met_corrected.Px(),
                #ev.met_corrected.Py(),
                # RECORRECTED
                ev.met_init.Px(), # uncorrected type I pf met px (float)
                ev.met_init.Py(), # uncorrected type I pf met py (float)
                ev.gen_genPx, # generator Z/W/Higgs px (float)
                ev.gen_genPy, # generator Z/W/Higgs py (float)
                ev.gen_visPx, # generator visible Z/W/Higgs px (float)
                ev.gen_visPy, # generator visible Z/W/Higgs py (float)
                N_RECOIL_JETS  # max recoil (I checked that -- plot in scrap/recoil-corrections-study)
                #ev.nalljets  # number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
                )

            met_y = ROOT.met_pt_recoilcor_y(
                #ev.met_corrected.Px(),
                #ev.met_corrected.Py(),
                ev.met_init.Px(), # uncorrected type I pf met px (float)
                ev.met_init.Py(), # uncorrected type I pf met py (float)
                ev.gen_genPx, # generator Z/W/Higgs px (float)
                ev.gen_genPy, # generator Z/W/Higgs py (float)
                ev.gen_visPx, # generator visible Z/W/Higgs px (float)
                ev.gen_visPy, # generator visible Z/W/Higgs py (float)
                N_RECOIL_JETS
                #ev.nalljets  # number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
                )

            #Mt_lep_met_c   = ROOT.MTlep_met_pt_recoilcor(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), ev.met_corrected.Px(), ev.met_corrected.Py(), ev.gen_genPx, ev.gen_genPy, ev.gen_visPx, ev.gen_visPy, 0)
            #Mt_lep_met = transverse_mass_pts(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), ev.pfmetcorr_ex, ev.pfmetcorr_ey)
            #Mt_tau_met_nominal = transverse_mass_pts(ev.tau_p4[0].Px(), ev.tau_p4[0].Py(), ev.pfmetcorr_ex, ev.pfmetcorr_ey)
            # TODO: tau is corrected with systematic ES

        elif not isMC and METMuEGClean:
            met_x = ev.met_slimmedMETsMuEGClean.Px()
            met_y = ev.met_slimmedMETsMuEGClean.Py()
        # to test: is it = met_init

        else:
            #met_x = ev.met_corrected.Px()
            #met_y = ev.met_corrected.Py()
            # RECORRECTED
            # use miniaod met and jets, reapply corrections of the passed jets to met
            met_x = ev.met_init.Px()
            met_y = ev.met_init.Py()
            # NOMINAL NTUPLE CORRECTED
            #met_x = ev.met_corrected.Px()
            #met_y = ev.met_corrected.Py()

            #Mt_lep_met_c   = ROOT.MTlep_c(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), ev.met_corrected.Px(), ev.met_corrected.Py())
            #Mt_lep_met_test = transverse_mass_pts(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), met_x, met_y)
            #Mt_lep_met = transverse_mass(ev.lep_p4[0], ev.met_corrected)
            #Mt_tau_met = transverse_mass(ev.tau_p4[0], ev.met_corrected)
            #Mt_tau_met_nominal = transverse_mass_pts(ev.tau_p4[0].Px(), ev.tau_p4[0].Py(), ev.met_corrected.Px(), ev.met_corrected.Py())
        # also
        #Mt_lep_met_d = (ev.lep_p4[0] + ev.met_corrected).Mt()

        event_met_init.SetPx(ev.met_init.Px())
        event_met_init.SetPy(ev.met_init.Py())
        event_met_init2.SetPx(met_x)
        event_met_init2.SetPy(met_y)
        event_met_lep_mt_init[0]  = transverse_mass_pts(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), ev.met_init.Px(), ev.met_init.Py())
        event_met_lep_mt_init2[0] = transverse_mass_pts(ev.lep_p4[0].Px(), ev.lep_p4[0].Py(), met_x, met_y)

        proc_met         = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(met_x, met_y, 0., 0.)
        proc_met_JERUp   = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(met_x, met_y, 0., 0.)
        proc_met_JERDown = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(met_x, met_y, 0., 0.)
        proc_met_JESUp   = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(met_x, met_y, 0., 0.)
        proc_met_JESDown = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(met_x, met_y, 0., 0.)
        proc_met_TESUp   = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(met_x, met_y, 0., 0.)
        proc_met_TESDown = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(met_x, met_y, 0., 0.)
        proc_met_lepsub  = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(met_x, met_y, 0., 0.)

        # from PU tests, leaving it here for now
        if not isMC:
            weight_pu_sum = 1.
            weight_pu_b = 1.
            weight_pu_h2 = 1.

        weight_init = 1. # common weight of event (1. for data), which also includes the weights averaged per datataking epoch scale factors (PU and muon LEP SFs)
        weight_pu = 1.
        weight_pu_el = 1.
        weight_pu_mu = 1.
        weight_pu_bcdef = 1.
        weight_pu_gh    = 1.
        #weights_th = namedtuple('th_weights', 'AlphaSUp AlphaSDown FragUp FragDown')
        #weights_th = (1., 1., 1., 1.)
        weights_gen_weight_alphas = (1., 1.)
        weights_Frag   = (1., 1.)
        weights_gen_weight_centralFrag   = 1.

        gen_proc_id[0] = 0
        if isMC:
            try:
                weight_pu_el    = pileup_ratio_ele     [ev.nvtx_gen] # reduced json
                weight_pu_el_up = pileup_ratio_up_ele  [ev.nvtx_gen]
                weight_pu_el_dn = pileup_ratio_down_ele[ev.nvtx_gen]
                weight_pu_mu    = pileup_ratio[ev.nvtx_gen] # golden json
                weight_pu_mu_up = pileup_ratio_up[ev.nvtx_gen]
                weight_pu_mu_dn = pileup_ratio_down[ev.nvtx_gen]

                # per run-epoch PU
                #weight_pu_bcdef,    weight_pu_gh    = calc_pu_per_runs(weight_pu_el, weight_pu_mu)
                #weight_pu_bcdef_up, weight_pu_gh_up = calc_pu_per_runs(weight_pu_el_up, weight_pu_mu_up)
                #weight_pu_bcdef_dn, weight_pu_gh_dn = calc_pu_per_runs(weight_pu_el_dn, weight_pu_mu_dn)
                weight_pu_bcdef,    weight_pu_gh    = pu_per_epoch_BCDEF     [ev.nvtx_gen], pu_per_epoch_GH      [ev.nvtx_gen]
                weight_pu_bcdef_up, weight_pu_gh_up = pu_per_epoch_BCDEF_Up  [ev.nvtx_gen], pu_per_epoch_GH_Up   [ev.nvtx_gen]
                weight_pu_bcdef_dn, weight_pu_gh_dn = pu_per_epoch_BCDEF_Down[ev.nvtx_gen], pu_per_epoch_GH_Down [ev.nvtx_gen]

                if False and passed_ele_event: # v37+ run Ele w golden json
                    weight_pu    = weight_pu_el   
                    weight_pu_up = weight_pu_el_up
                    weight_pu_dn = weight_pu_el_dn
                    control_hs['weight_pu_el']   .Fill(weight_pu)
                    control_hs['weight_pu_el_up'].Fill(weight_pu_up)
                    control_hs['weight_pu_el_dn'].Fill(weight_pu_dn)
                else:
                    weight_pu    = weight_pu_mu   
                    weight_pu_up = weight_pu_mu_up
                    weight_pu_dn = weight_pu_mu_dn
                    control_hs['weight_pu_mu']   .Fill(weight_pu)
                    control_hs['weight_pu_mu_up'].Fill(weight_pu_up)
                    control_hs['weight_pu_mu_dn'].Fill(weight_pu_dn)
                # and the new PU-s
                weight_pu_sum  = pileup_ratio_sum[ev.nvtx_gen]
                weight_pu_b    = pileup_ratio_b[ev.nvtx_gen]
                weight_pu_h2   = pileup_ratio_h2[ev.nvtx_gen]
                control_hs['weight_pu_sum'] .Fill(weight_pu_sum)
                control_hs['weight_pu_b']   .Fill(weight_pu_b)
                control_hs['weight_pu_h2']  .Fill(weight_pu_h2)
            except:
                #print i, ev.nvtx
                continue

            if aMCatNLO:
                amcatnlo_w[0] = ev.aMCatNLO_weight
                if ev.aMCatNLO_weight < 0:
                    weight_init *= -1

            weight_z_mass_pt = 1.
            if isDY:
                # float zPtMass_weight(float genMass, float genPt)
                weight_z_mass_pt *= zPtMass_weight(ev.genMass, ev.genPt)
                weight_init *= weight_z_mass_pt
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
                    gen_proc_id[0] = genproc_dy_tautau
                else:
                    gen_proc_id[0] = genproc_dy_other

            if isWJets:
                if ev.gen_N_wdecays > 0:
                    lep1_id = abs(ev.gen_wdecays_IDs[0])
                #W does not need these prompt leptons, but I keep them just in case
                else:
                    # check prompt leptns
                    lep1_id = abs(ev.gen_pythia8_prompt_leptons_IDs[0])
                if lep1_id >= 15*15:
                    gen_proc_id[0] = genproc_wjets_tauh
                elif lep1_id >= 15: # 15*11 and 15*13
                    gen_proc_id[0] = genproc_wjets_taul
                # we use only WToLNu now, W->j does not happen
                #elif lep1_id == 1:
                #    proc = 'wjets_j'
                else:
                    gen_proc_id[0] = genproc_wjets

            weight_top_pt = 1.
            # "Only top quarks in SM ttbar events must be reweighted, not single tops or tops from BSM production mechanisms."
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
            if isTT:
                # th weights if needed
                #if with_AlphaS_sys:
                if with_PDF_sys:
                    #weights_gen_weight_norm = ev.gen_weight_norm
                    # weight_norm = 1 always
                    weights_gen_weight_alphas = (ev.gen_weight_alphas_1, ev.gen_weight_alphas_2)
                    # norm to average
                    #weights_gen_weight_norm = (weights_gen_weight_alphas[0] + weights_gen_weight_alphas[1]) / 2
                    weights_gen_weight_norm = ev.gen_weights_pdf_hessians[0]
                    # norm is the nominal PDF
                    if 0. <= weights_gen_weight_norm < 0.00001:
                        weights_gen_weight_norm = 0.00001
                    elif 0. > weights_gen_weight_norm > -0.00001:
                        weights_gen_weight_norm = -0.00001

                    #control_hs['weights_gen_weight_norm']   .Fill(ev.gen_weight_too)
                    control_hs['weights_gen_weight_too']        .Fill(ev.gen_weight_too)
                    control_hs['weights_gen_weight_norm']       .Fill(ev.gen_weight_norm)
                    control_hs['weights_gen_weight_average']    .Fill(weights_gen_weight_norm)
                    control_hs['weights_gen_weight_alphasUp']   .Fill(weights_gen_weight_alphas[0])
                    control_hs['weights_gen_weight_alphasDown'] .Fill(weights_gen_weight_alphas[1])

                    weights_gen_weight_alphas_up = weights_gen_weight_alphas[0] / weights_gen_weight_norm
                    weights_gen_weight_alphas_dn = weights_gen_weight_alphas[1] / weights_gen_weight_norm

                    # sane high weight in tt->taumu tauh as in pdfs
                    if weights_gen_weight_alphas_up > 2.:
                        weights_gen_weight_alphas_up = 1.
                    if weights_gen_weight_alphas_dn > 2.:
                        weights_gen_weight_alphas_dn = 1.

                if with_Frag_sys:
                    weights_gen_weight_centralFrag = ev.gen_weight_centralFrag if ev.gen_weight_centralFrag > 0. else 0.00001
                    weights_gen_weight_Frag = (ev.gen_weight_FragUp, ev.gen_weight_FragDown)
                    weights_gen_weight_semilepbr = (ev.gen_weight_semilepbrUp, ev.gen_weight_semilepbrDown)
                    weights_gen_weight_Peterson = ev.gen_weight_PetersonFrag
                    # sub to this naming in the following ntuple runs:
                    #weights_gen_weight_Frag = (ev.gen_weight_FragUp, ev.gen_weight_FragDown)
                    control_hs['weights_gen_weight_centralFrag']   .Fill(weights_gen_weight_centralFrag)
                    control_hs['weights_gen_weight_FragUp']   .Fill(weights_gen_weight_Frag[0])
                    control_hs['weights_gen_weight_FragDown'] .Fill(weights_gen_weight_Frag[1])
                    control_hs['weights_gen_weight_Peterson'] .Fill(weights_gen_weight_Peterson)
                    control_hs['weights_gen_weight_semilepbrUp']   .Fill(weights_gen_weight_semilepbr[0])
                    control_hs['weights_gen_weight_semilepbrDown'] .Fill(weights_gen_weight_semilepbr[1])

                if with_MEscale_sys:
                    """ these are calculated "with respect to the sum of nominal weights":
                      > Depending on the case, one may want to normalize the per-event weights
                      > to the sum of weights corresponding to the default scale choice...
                        i.e. they are not per-event normalized to the nominal weight
                        but they are normalized afterwards to the overall weight, i.e. to the average of the nominal
                    """
                    weights_gen_weight_nom   = ev.gen_weights_renorm_fact[MUf_nom_MUr_nom] if ev.gen_weights_renorm_fact[MUf_nom_MUr_nom] > 0. else 0.00001
                    weights_gen_weight_f_rUp = ev.gen_weights_renorm_fact[MUf_nom_MUr_up]
                    weights_gen_weight_f_rDn = ev.gen_weights_renorm_fact[MUf_nom_MUr_down]
                    weights_gen_weight_fUp_r = ev.gen_weights_renorm_fact[MUf_up_MUr_nom]
                    weights_gen_weight_fDn_r = ev.gen_weights_renorm_fact[MUf_down_MUr_nom]
                    weights_gen_weight_frUp  = ev.gen_weights_renorm_fact[MUf_up_MUr_up]
                    weights_gen_weight_frDn  = ev.gen_weights_renorm_fact[MUf_down_MUr_down]

                    # sub to this naming in the following ntuple runs:
                    #weights_gen_weight_Frag = (ev.gen_weight_FragUp, ev.gen_weight_FragDown)
                    control_hs['weights_gen_weight_nom']   .Fill(weights_gen_weight_nom  )
                    control_hs['weights_gen_weight_f_rUp'] .Fill(weights_gen_weight_f_rUp)
                    control_hs['weights_gen_weight_f_rDn'] .Fill(weights_gen_weight_f_rDn)
                    control_hs['weights_gen_weight_fUp_r'] .Fill(weights_gen_weight_fUp_r)
                    control_hs['weights_gen_weight_fDn_r'] .Fill(weights_gen_weight_fDn_r)
                    control_hs['weights_gen_weight_frUp']  .Fill(weights_gen_weight_frUp )
                    control_hs['weights_gen_weight_frDn']  .Fill(weights_gen_weight_frDn )

                weight_top_pt = ttbar_pT_SF(ev.gen_t_pt, ev.gen_tb_pt)
                #weight *= weight_top_pt # to sys
                control_hs['weight_top_pt']   .Fill(weight_top_pt)
                # basically only difference is eltau/mutau
                t_wid  = abs(ev.gen_t_w_decay_id)
                tb_wid = abs(ev.gen_tb_w_decay_id)
                ## save detailed proc control for elmu "other" contribution
                #tt_ids = [0, 0]
                #tt_id1 = 0
                #for i, t_w_id in enumerate((t_wid, tb_wid)):
                #    if t_w_id == 1:
                #        tt_ids[i] = 1
                #    elif t_w_id == 11:
                #        tt_ids[i] = 2
                #    elif t_w_id == 13:
                #        tt_ids[i] = 3
                #    elif t_w_id == 15*11:
                #        tt_ids[i] = 4
                #    elif t_w_id == 15*13:
                #        tt_ids[i] = 5
                #    elif 15*15 < t_w_id < 15*30:
                #        tt_ids[i] = 6
                #    elif t_w_id >= 15*30:
                #        tt_ids[i] = 7
                #    else:
                #        tt_ids[i] = 8

                if (t_wid > 15*15 and tb_wid == 13) or (t_wid == 13 and tb_wid > 15*15): # lt
                    # check if tau decayed in 3 charged particles
                    if (t_wid >= 15*30 and tb_wid == 13) or (t_wid == 13 and tb_wid >= 15*30): # lt
                        gen_proc_id[0] = genproc_tt_mutau3ch
                    else:
                        gen_proc_id[0] = genproc_tt_mutau
                elif (t_wid > 15*15 and tb_wid == 11) or (t_wid == 11 and tb_wid > 15*15): # lt
                    if (t_wid >= 15*30 and tb_wid == 11) or (t_wid == 11 and tb_wid >= 15*30): # lt
                        gen_proc_id[0] = genproc_tt_eltau3ch
                    else:
                        gen_proc_id[0] = genproc_tt_eltau
                elif ev.gen_t_w_decay_id * ev.gen_tb_w_decay_id == -13*11:
                    gen_proc_id[0] = genproc_tt_elmu
                #elif t_wid * tb_wid == 15*13*15*11: # this should work, but:
                elif (t_wid == 15*13 and tb_wid == 15*11) or (t_wid == 15*11 and tb_wid == 15*13):
                    gen_proc_id[0] = genproc_tt_taueltaumu
                elif (t_wid  == 13 and tb_wid == 15*11) or (t_wid  == 11 and tb_wid == 15*13) or \
                     (tb_wid == 13 and t_wid  == 15*11) or (tb_wid == 11 and  t_wid == 15*13):
                    gen_proc_id[0] = genproc_tt_ltaul # opposite leptons -- el-taumu etc
                elif t_wid * tb_wid == 13 or t_wid * tb_wid == 11: # lj
                    gen_proc_id[0] = genproc_tt_lj
                elif (t_wid > 15*15 and (tb_wid == 11*15 or tb_wid == 13*15)) or \
                     ((t_wid == 11*15 or t_wid == 13*15) and tb_wid > 15*15): # taul tauh
                    gen_proc_id[0] = genproc_tt_taultauh
                elif (t_wid == 1     and tb_wid == 13*15) or \
                     (t_wid == 13*15 and tb_wid == 1)     or \
                     (t_wid == 1     and tb_wid == 11*15) or \
                     (t_wid == 11*15 and tb_wid == 1):
                    gen_proc_id[0] = genproc_tt_taulj
                else:
                    gen_proc_id[0] = genproc_tt_other

            if isSTop:
                # basically only difference is eltau/mutau
                #w1_id = abs(ev.gen_wdecays_IDs[0])
                # check the top decay insted
                w1_id = abs(ev.gen_t_w_decay_id) if ev.gen_t_w_decay_id != -111 else abs(ev.gen_tb_w_decay_id)
                w2_id = 1
                if not isSTopTSchannels:
                    # in tW there is an additional W produced with top
                    w2_id = abs(ev.gen_wdecays_IDs[0])
                    #w2_id = abs(ev.gen_t_w_decay_id) if ev.gen_t_w_decay_id != -111 else abs(ev.gen_tb_w_decay_id)
                # t/s channels emit top and a quark -- top ID + jets
                if (w1_id > 15*15 and w2_id == 13) or (w1_id == 13 and w2_id > 15*15): # lt
                    gen_proc_id[0] = genproc_stop_mu
                elif (w1_id > 15*15 and w2_id == 11) or (w1_id == 11 and w2_id > 15*15): # lt
                    gen_proc_id[0] = genproc_stop_el
                elif (w1_id == 11 and w2_id == 13) or (w1_id == 13 and w2_id == 11): # is it faster than comparing to product?
                    gen_proc_id[0] = genproc_stop_elmu
                elif w1_id * w2_id == 13 or w1_id * w2_id == 11: # lj
                    gen_proc_id[0] = genproc_stop_lj
                #elif (w1_id > 15*15 and (w2_id == 11*15 or w2_id == 13*15)) or
                #     ((w1_id == 11*15 or w1_id == 13*15) and w2_id > 15*15): # taul tauh
                #    proc = 'tt_taultauh'
                else:
                    gen_proc_id[0] = genproc_stop_other

            # LEPTON SFs
            # and PU SFs averaged together
            weight_lep_pu    = 1.

            weight_lep_pu_elIDUp     = 1.
            weight_lep_pu_elIDDown   = 1.
            weight_lep_pu_elTRGUp    = 1.
            weight_lep_pu_elTRGDown  = 1.

            weight_lep_pu_PUUp    = 1.
            weight_lep_pu_PUDown  = 1.

            weight_lep_pu_muIDUp     = 1.
            weight_lep_pu_muIDDown   = 1.
            weight_lep_pu_muTRGUp    = 1.
            weight_lep_pu_muTRGDown  = 1.

            # also with no systematics TODO: add the lepton and PU systematics here
            weight_lepall_pu = 1. # DONE

            # separate lep SF weights for control (possibility to reconstruct the whole weight)
            weight_lep_EL_id    = 1.
            weight_lep_EL_trg   = 1.

            #weight_lepMU0_id     = 1.
            #weight_lepMU0_trg    = 1.
            #weight_lepMU0_trg_Up = 1.

            #weight_lep_MU_id    = 1.
            #weight_lep_MU_trg   = 1.
            #weight_lep_MU_id_Up    = 1.
            #weight_lep_MU_id_Down  = 1.
            #weight_lep_MU_trg_Up   = 1.
            #weight_lep_MU_trg_Down = 1.

            # 
            weight_lep_B_MU_id    = 1.
            weight_lep_B_MU_trg   = 1.
            weight_lep_H_MU_id    = 1.
            weight_lep_H_MU_trg   = 1.

            # no systematics in controls for now
            #weight_lep_EL_id_Up    = 1.
            #weight_lep_EL_id_Down  = 1.
            #weight_lep_EL_trg_Up   = 1.
            #weight_lep_EL_trg_Down = 1.
            #weight_lep_B_MU_id_Up    = 1.
            #weight_lep_B_MU_id_Down  = 1.
            #weight_lep_B_MU_trg_Up   = 1.
            #weight_lep_B_MU_trg_Down = 1.
            #weight_lep_H_MU_id_Up    = 1.
            #weight_lep_H_MU_id_Down  = 1.
            #weight_lep_H_MU_trg_Up   = 1.
            #weight_lep_H_MU_trg_Down = 1.

            # PU and LEP weights are per-epoch BCDEF and GH (only muon LEP SF, electron is common for whole year)
            # therefore they are averaged with per-epoch ratios together in the basic MC weight

            if isMC and pass_mu_all:
                #mu_sfs = lepton_muon_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt()) # old
                # 0, 1, 2 -- trk
                # 3, 4    -- id, iso
                mu_sfs_b, mu_sfs_h = lepton_muon_SF(ev.lep_alliso_p4[0].eta(), ev.lep_alliso_p4[0].pt(), ev.nvtx, ev.nvtx_gen)
                (mu_trg_sf_b, trg_b_unc), (mu_trg_sf_h, trg_h_unc) = lepton_muon_trigger_SF(ev.lep_alliso_p4[0].eta(), ev.lep_alliso_p4[0].pt())
                # bcdef gh eras
                # with gen_nvtx for tracking eff
                mu_b_trk, _ = mu_sfs_b[0:2] # unc in mu_sfs_b[1]
                mu_h_trk, _ = mu_sfs_h[0:2] # unc in mu_sfs_h[1]

                # test how much lepton SFs affect
                weight_lepall_pu = ratio_bcdef * mu_trg_sf_b * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * weight_pu_bcdef + \
                                      ratio_gh * mu_trg_sf_h * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * weight_pu_gh

                #weight_lep_pu_MU_id = ratio_bcdef * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * weight_pu_bcd + \
                #                         ratio_gh * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * weight_pu_gh
                #weight_lep_pu_MU_trg = ratio_bcdef * mu_trg_sf_b * weight_pu_bcdef + ratio_gh * mu_trg_sf_h * weight_pu_gh


            if isMC and pass_mu:
                #mu_sfs = lepton_muon_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt()) # old
                # 0, 1, 2 -- trk
                # 3, 4    -- id, iso
                mu_sfs_b, mu_sfs_h = lepton_muon_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt(), ev.nvtx, ev.nvtx_gen) # running tracking SF on reco nvtx and output on nvtx_gen for control
                (mu_trg_sf_b, trg_b_unc), (mu_trg_sf_h, trg_h_unc) = lepton_muon_trigger_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt())
                # bcdef gh eras
                #weight *= ratio_bcdef * mu_trg_sf[0] * mu_sfs_b[0] * mu_sfs_b[1] * mu_sfs_b[2] * mu_sfs_b[3] + ratio_gh * mu_trg_sf[1] * mu_sfs_h[0] * mu_sfs_h[1] * mu_sfs_h[2] * mu_sfs_h[3]
                # with gen_nvtx for tracking eff
                mu_b_trk, mu_b_trk_u = mu_sfs_b[0:2] # unc in mu_sfs_b[1]
                mu_h_trk, mu_h_trk_u = mu_sfs_h[0:2] # unc in mu_sfs_h[1]

                weight_lep_pu_muIDUp   = ratio_bcdef * mu_trg_sf_b * (mu_b_trk + mu_b_trk_u) * (mu_sfs_b[3][0] + mu_sfs_b[3][1]) * (mu_sfs_b[4][0] + mu_sfs_b[4][1]) * weight_pu_bcdef + \
                                            ratio_gh * mu_trg_sf_h * (mu_h_trk + mu_h_trk_u) * (mu_sfs_h[3][0] + mu_sfs_h[3][1]) * (mu_sfs_h[4][0] + mu_sfs_h[4][1]) * weight_pu_gh
                weight_lep_pu_muIDDown = ratio_bcdef * mu_trg_sf_b * (mu_b_trk - mu_b_trk_u) * (mu_sfs_b[3][0] - mu_sfs_b[3][1]) * (mu_sfs_b[4][0] - mu_sfs_b[4][1]) * weight_pu_bcdef + \
                                            ratio_gh * mu_trg_sf_h * (mu_h_trk - mu_h_trk_u) * (mu_sfs_h[3][0] - mu_sfs_h[3][1]) * (mu_sfs_h[4][0] - mu_sfs_h[4][1]) * weight_pu_gh
                weight_lep_pu_muTRGUp   = ratio_bcdef * (mu_trg_sf_b + trg_b_unc) * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * weight_pu_bcdef + \
                                             ratio_gh * (mu_trg_sf_h + trg_h_unc) * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * weight_pu_gh
                weight_lep_pu_muTRGDown = ratio_bcdef * (mu_trg_sf_b - trg_b_unc) * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * weight_pu_bcdef + \
                                             ratio_gh * (mu_trg_sf_h - trg_h_unc) * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * weight_pu_gh

                weight_lep_pu_PUUp = ratio_bcdef * mu_trg_sf_b * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * weight_pu_bcdef_up + \
                                        ratio_gh * mu_trg_sf_h * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * weight_pu_gh_up
                weight_lep_pu_PUDown = ratio_bcdef * mu_trg_sf_b * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * weight_pu_bcdef_dn + \
                                          ratio_gh * mu_trg_sf_h * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * weight_pu_gh_dn

                # test how much lepton SFs affect
                weight_lep_pu = ratio_bcdef * mu_trg_sf_b * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * weight_pu_bcdef + \
                                   ratio_gh * mu_trg_sf_h * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * weight_pu_gh

                # separate weights for control and reconstruction of weights
                weight_lep_B_MU_id = mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0]
                weight_lep_H_MU_id = mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0]

                weight_lep_B_MU_trg = mu_trg_sf_b
                weight_lep_H_MU_trg = mu_trg_sf_h

            if isMC and (pass_mumu or pass_mumu_ss):
                # our single-muon HLT demands an OR operation on the 2 muons
                #(mu_trg_sf_b, trg_b_unc), (mu_trg_sf_h, trg_h_unc) = lepton_muon_trigger_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt())
                (mu_trg_sf_b1, trg_b_unc1), (mu_trg_sf_h1, trg_h_unc1) = lepton_muon_trigger_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt())
                (mu_trg_sf_b2, trg_b_unc2), (mu_trg_sf_h2, trg_h_unc2) = lepton_muon_trigger_SF(ev.lep_p4[1].eta(), ev.lep_p4[1].pt())
                mu_trg_sf_b, trg_b_unc = dilepton_or_sfs(mu_trg_sf_b1, trg_b_unc1, mu_trg_sf_b2, trg_b_unc2)
                mu_trg_sf_h, trg_h_unc = dilepton_or_sfs(mu_trg_sf_h1, trg_h_unc1, mu_trg_sf_h2, trg_h_unc2)

                #mu_sfs = lepton_muon_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt()) # old
                # 0, 1, 2 -- trk
                # 3, 4    -- id, iso
                mu_sfs_b1, mu_sfs_h1 = lepton_muon_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt(), ev.nvtx, ev.nvtx_gen) # running tracking SF on reco nvtx and output on nvtx_gen for control
                mu_sfs_b2, mu_sfs_h2 = lepton_muon_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt(), ev.nvtx, ev.nvtx_gen)
                # the ID SFs are an AND operation on both muons

                # bcdef gh eras
                #weight *= ratio_bcdef * mu_trg_sf[0] * mu_sfs_b[0] * mu_sfs_b[1] * mu_sfs_b[2] * mu_sfs_b[3] + ratio_gh * mu_trg_sf[1] * mu_sfs_h[0] * mu_sfs_h[1] * mu_sfs_h[2] * mu_sfs_h[3]
                # with gen_nvtx for tracking eff
                mu_b_trk, mu_b_trk_u = dilepton_and_sfs(mu_sfs_b1[0], mu_sfs_b1[1], mu_sfs_b2[0], mu_sfs_b2[1])
                mu_h_trk, mu_h_trk_u = dilepton_and_sfs(mu_sfs_h1[0], mu_sfs_h1[1], mu_sfs_h2[0], mu_sfs_h2[1])

                mu_id_b,  mu_id_b_unc  = dilepton_and_sfs(mu_sfs_b1[3][0], mu_sfs_b1[3][1], mu_sfs_b2[3][0], mu_sfs_b2[3][1])
                mu_iso_b, mu_iso_b_unc = dilepton_and_sfs(mu_sfs_b1[4][0], mu_sfs_b1[4][1], mu_sfs_b2[4][0], mu_sfs_b2[4][1])

                mu_id_h,  mu_id_h_unc  = dilepton_and_sfs(mu_sfs_h1[3][0], mu_sfs_h1[3][1], mu_sfs_h2[3][0], mu_sfs_h2[3][1])
                mu_iso_h, mu_iso_h_unc = dilepton_and_sfs(mu_sfs_h1[4][0], mu_sfs_h1[4][1], mu_sfs_h2[4][0], mu_sfs_h2[4][1])

                weight_lep_pu_muIDUp   = ratio_bcdef * mu_trg_sf_b * (mu_b_trk + mu_b_trk_u) * (mu_id_b + mu_id_b_unc) * (mu_iso_b + mu_iso_b_unc) * weight_pu_bcdef + \
                                            ratio_gh * mu_trg_sf_h * (mu_h_trk + mu_h_trk_u) * (mu_id_h + mu_id_h_unc) * (mu_iso_h + mu_iso_h_unc) * weight_pu_gh
                weight_lep_pu_muIDDown = ratio_bcdef * mu_trg_sf_b * (mu_b_trk - mu_b_trk_u) * (mu_id_b - mu_id_b_unc) * (mu_iso_b - mu_iso_b_unc) * weight_pu_bcdef + \
                                            ratio_gh * mu_trg_sf_h * (mu_h_trk - mu_h_trk_u) * (mu_id_h - mu_id_h_unc) * (mu_iso_h - mu_iso_h_unc) * weight_pu_gh

                weight_lep_pu_muTRGUp = ratio_bcdef * (mu_trg_sf_b + trg_b_unc) * mu_b_trk * mu_id_b * mu_iso_b  * weight_pu_bcdef + \
                                           ratio_gh * (mu_trg_sf_h + trg_h_unc) * mu_h_trk * mu_id_h * mu_iso_h * weight_pu_gh
                weight_lep_pu_muTRGUp = ratio_bcdef * (mu_trg_sf_b - trg_b_unc) * mu_b_trk * mu_id_b * mu_iso_b  * weight_pu_bcdef + \
                                           ratio_gh * (mu_trg_sf_h - trg_h_unc) * mu_h_trk * mu_id_h * mu_iso_h * weight_pu_gh

                # tests and control
                #weight_lepMU0_trg_Up   = (ratio_bcdef * (mu_trg_sf_b1 + trg_b_unc1) * mu_sfs_b1[0] * mu_sfs_b1[3][0] * mu_sfs_b1[4][0] * weight_pu_bcdef + \
                #                             ratio_gh * (mu_trg_sf_h1 + trg_h_unc1) * mu_sfs_h1[0] * mu_sfs_h1[3][0] * mu_sfs_h1[4][0] * weight_pu_gh   )

                # test how much lepton SFs affect
                weight_lep_pu = ratio_bcdef * mu_trg_sf_b * mu_b_trk * mu_id_b * mu_iso_b * weight_pu_bcdef + \
                                   ratio_gh * mu_trg_sf_h * mu_h_trk * mu_id_h * mu_iso_h * weight_pu_gh

                weight_lep_pu_PUUp = ratio_bcdef * mu_trg_sf_b * mu_b_trk * mu_id_b * mu_iso_b * weight_pu_bcdef_up + \
                                        ratio_gh * mu_trg_sf_h * mu_h_trk * mu_id_h * mu_iso_h * weight_pu_gh_up
                weight_lep_pu_PUDown = ratio_bcdef * mu_trg_sf_b * mu_b_trk * mu_id_b * mu_iso_b * weight_pu_bcdef_dn + \
                                          ratio_gh * mu_trg_sf_h * mu_h_trk * mu_id_h * mu_iso_h * weight_pu_gh_dn

                weight_lep_B_MU_id = mu_b_trk * mu_id_b * mu_iso_b
                weight_lep_H_MU_id = mu_h_trk * mu_id_h * mu_iso_h
                weight_lep_B_MU_trg = mu_trg_sf_b
                weight_lep_H_MU_trg = mu_trg_sf_h


            # for now let's just use 1 elmu and it is electron-triggered one
            if isMC and pass_elmu:
                # find which lepton is mu and which is el
                mu_n, el_n = (0, 1) if abs(ev.lep_id[0]) == 13 else (1, 0)
                mu_sfs_b, mu_sfs_h     = lepton_muon_SF(ev.lep_p4[mu_n].eta(), ev.lep_p4[mu_n].pt(), ev.nvtx, ev.nvtx_gen)
                #(mu_trg_sf_b, trg_b_unc), (mu_trg_sf_h, trg_h_unc) = lepton_muon_trigger_SF(ev.lep_p4[mu_n].eta(), ev.lep_p4[mu_n].pt())
                el_sfs_reco, el_sfs_id = lepton_electron_SF(ev.lep_p4[el_n].eta(), ev.lep_p4[el_n].pt())
                #el_trg_sf, el_trg_unc  = lepton_electron_trigger_SF(ev.lep_p4[el_n].eta(), ev.lep_p4[el_n].pt())

                (mu_trg_sf_b, trg_b_unc), (mu_trg_sf_h, trg_h_unc) = lepton_muon_trigger_SF(ev.lep_p4[mu_n].eta(), ev.lep_p4[0].pt())

                #weight_lep_Up   *= weight * (el_trg_sf[0] + el_trg_sf[1]) * (el_sfs_reco[0] + el_sfs_reco[1]) * (el_sfs_id[0] + el_sfs_id[1])
                #weight_lep_Down *= weight * (el_trg_sf[0] - el_trg_sf[1]) * (el_sfs_reco[0] - el_sfs_reco[1]) * (el_sfs_id[0] - el_sfs_id[1])
                #weight *= el_trg_sf[0] * el_sfs_reco[0] * el_sfs_id[0]
                # w.o. el trig
                mu_b_trk, mu_b_trk_u = mu_sfs_b[0:2] # unc in mu_sfs_b[1]
                mu_h_trk, mu_h_trk_u = mu_sfs_h[0:2] # unc in mu_sfs_h[1]

                weight_lep_el = el_sfs_reco[0] * el_sfs_id[0]

                #weight_lep_mu = ratio_bcdef * mu_trg_sf_b * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] + \
                #                ratio_gh    * mu_trg_sf_h * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0]
                weight_lep_mu = ratio_bcdef * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * mu_trg_sf_b * weight_pu_bcdef + \
                                ratio_gh    * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * mu_trg_sf_h * weight_pu_gh

                weight_lep_mu_PUUp = ratio_bcdef * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * mu_trg_sf_b * weight_pu_bcdef_up + \
                                     ratio_gh    * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * mu_trg_sf_h * weight_pu_gh_up

                weight_lep_mu_PUDown = ratio_bcdef * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * mu_trg_sf_b * weight_pu_bcdef_dn + \
                                       ratio_gh    * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * mu_trg_sf_h * weight_pu_gh_dn

                weight_lep_pu_muTRGUp   = ratio_bcdef * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * (mu_trg_sf_b + trg_b_unc) * weight_pu_bcdef + \
                                          ratio_gh    * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * (mu_trg_sf_h + trg_h_unc) * weight_pu_gh

                weight_lep_pu_muTRGDown = ratio_bcdef * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * (mu_trg_sf_b - trg_b_unc) * weight_pu_bcdef + \
                                          ratio_gh    * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * (mu_trg_sf_h - trg_h_unc) * weight_pu_gh

                weight_lep_pu_muTRGUp   = weight_lep_pu
                weight_lep_pu_muTRGDown = weight_lep_pu

                # controls
                #weight_lep_EL_trg = el_trg_sf
                weight_lep_EL_id  = weight_lep_el
                weight_lep_B_MU_id  = mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0]
                weight_lep_H_MU_id  = mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0]
                # no muon trigger here

                # the overall weight and systematic
                weight_lep_pu = weight_lep_el * weight_lep_mu

                weight_lep_pu_PUUp   = weight_lep_el * weight_lep_mu_PUUp
                weight_lep_pu_PUDown = weight_lep_el * weight_lep_mu_PUDown

                weight_lep_pu_muIDUp    = weight_lep_el * \
                                         (ratio_bcdef * mu_trg_sf_b * (mu_b_trk + mu_b_trk_u) * (mu_sfs_b[3][0] + mu_sfs_b[3][1]) * (mu_sfs_b[4][0] + mu_sfs_b[4][1]) * weight_pu_bcdef + \
                                             ratio_gh * mu_trg_sf_h * (mu_h_trk + mu_h_trk_u) * (mu_sfs_h[3][0] + mu_sfs_h[3][1]) * (mu_sfs_h[4][0] + mu_sfs_h[4][1]) * weight_pu_gh   )
                weight_lep_pu_muIDDown  = weight_lep_el * \
                                         (ratio_bcdef * mu_trg_sf_b * (mu_b_trk - mu_b_trk_u) * (mu_sfs_b[3][0] - mu_sfs_b[3][1]) * (mu_sfs_b[4][0] - mu_sfs_b[4][1]) * weight_pu_bcdef + \
                                             ratio_gh * mu_trg_sf_h * (mu_h_trk - mu_h_trk_u) * (mu_sfs_h[3][0] - mu_sfs_h[3][1]) * (mu_sfs_h[4][0] - mu_sfs_h[4][1]) * weight_pu_gh   )

                weight_lep_pu_elIDUp    = weight_lep_mu * (el_sfs_reco[0] + el_sfs_reco[1]) * (el_sfs_id[0] + el_sfs_id[1])
                weight_lep_pu_elIDDown  = weight_lep_mu * (el_sfs_reco[0] - el_sfs_reco[1]) * (el_sfs_id[0] - el_sfs_id[1])

                weight_lep_pu_elTRGUp   = weight_lep_pu
                weight_lep_pu_elTRGDown = weight_lep_pu


            if False and isMC and pass_elmu_el:
                # find which lepton is mu and which is el
                mu_n, el_n = (0, 1) if abs(ev.lep_id[0]) == 13 else (1, 0)
                mu_sfs_b, mu_sfs_h     = lepton_muon_SF(ev.lep_p4[mu_n].eta(), ev.lep_p4[mu_n].pt(), ev.nvtx, ev.nvtx_gen)
                #(mu_trg_sf_b, trg_b_unc), (mu_trg_sf_h, trg_h_unc) = lepton_muon_trigger_SF(ev.lep_p4[mu_n].eta(), ev.lep_p4[mu_n].pt())
                el_sfs_reco, el_sfs_id = lepton_electron_SF(ev.lep_p4[el_n].eta(), ev.lep_p4[el_n].pt())
                el_trg_sf, el_trg_unc  = lepton_electron_trigger_SF(ev.lep_p4[el_n].eta(), ev.lep_p4[el_n].pt())
                # TODO !!!!!!!!!!!!!!!!!!!!!!!!! ELE TRIGGER IS ASYMETRIC IN ETA! CHECK ID TOO!

                #weight_lep_Up   *= weight * (el_trg_sf[0] + el_trg_sf[1]) * (el_sfs_reco[0] + el_sfs_reco[1]) * (el_sfs_id[0] + el_sfs_id[1])
                #weight_lep_Down *= weight * (el_trg_sf[0] - el_trg_sf[1]) * (el_sfs_reco[0] - el_sfs_reco[1]) * (el_sfs_id[0] - el_sfs_id[1])
                #weight *= el_trg_sf[0] * el_sfs_reco[0] * el_sfs_id[0]
                # w.o. el trig
                mu_b_trk, mu_b_trk_u = mu_sfs_b[0:2] # unc in mu_sfs_b[1]
                mu_h_trk, mu_h_trk_u = mu_sfs_h[0:2] # unc in mu_sfs_h[1]

                weight_lep_el = el_sfs_reco[0] * el_sfs_id[0]

                #weight_lep_mu = ratio_bcdef * mu_trg_sf_b * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] + \
                #                ratio_gh    * mu_trg_sf_h * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0]
                weight_lep_mu = ratio_bcdef * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * weight_pu_bcdef + \
                                ratio_gh    * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * weight_pu_gh

                weight_lep_mu_PUUp = ratio_bcdef * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * weight_pu_bcdef_up + \
                                     ratio_gh    * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * weight_pu_gh_up

                weight_lep_mu_PUDown = ratio_bcdef * mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0] * weight_pu_bcdef_dn + \
                                       ratio_gh    * mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0] * weight_pu_gh_dn

                # controls
                weight_lep_EL_trg = el_trg_sf
                weight_lep_EL_id  = weight_lep_el
                weight_lep_B_MU_id  = mu_b_trk * mu_sfs_b[3][0] * mu_sfs_b[4][0]
                weight_lep_H_MU_id  = mu_h_trk * mu_sfs_h[3][0] * mu_sfs_h[4][0]
                # no muon trigger here

                # the overall weight and systematic
                weight_lep_pu = weight_lep_el * weight_lep_mu * el_trg_sf

                weight_lep_pu_PUUp   = weight_lep_el * weight_lep_mu_PUUp   * el_trg_sf
                weight_lep_pu_PUDown = weight_lep_el * weight_lep_mu_PUDown * el_trg_sf

                weight_lep_pu_muIDUp    = weight_lep_el * \
                                         (ratio_bcdef * (mu_b_trk + mu_b_trk_u) * (mu_sfs_b[3][0] + mu_sfs_b[3][1]) * (mu_sfs_b[4][0] + mu_sfs_b[4][1]) * weight_pu_bcdef + \
                                             ratio_gh * (mu_h_trk + mu_h_trk_u) * (mu_sfs_h[3][0] + mu_sfs_h[3][1]) * (mu_sfs_h[4][0] + mu_sfs_h[4][1]) * weight_pu_gh   )
                weight_lep_pu_muIDDown  = weight_lep_el * \
                                         (ratio_bcdef * (mu_b_trk - mu_b_trk_u) * (mu_sfs_b[3][0] - mu_sfs_b[3][1]) * (mu_sfs_b[4][0] - mu_sfs_b[4][1]) * weight_pu_bcdef + \
                                             ratio_gh * (mu_h_trk - mu_h_trk_u) * (mu_sfs_h[3][0] - mu_sfs_h[3][1]) * (mu_sfs_h[4][0] - mu_sfs_h[4][1]) * weight_pu_gh   )

                weight_lep_pu_muTRGUp   = weight_lep_pu
                weight_lep_pu_muTRGDown = weight_lep_pu

                weight_lep_pu_elIDUp    = weight_lep_mu * (el_sfs_reco[0] + el_sfs_reco[1]) * (el_sfs_id[0] + el_sfs_id[1]) * el_trg_sf
                weight_lep_pu_elIDDown  = weight_lep_mu * (el_sfs_reco[0] - el_sfs_reco[1]) * (el_sfs_id[0] - el_sfs_id[1]) * el_trg_sf

                weight_lep_pu_elTRGUp   = weight_lep_mu * el_sfs_reco[0] * el_sfs_id[0] * (el_trg_sf + el_trg_unc)
                weight_lep_pu_elTRGDown = weight_lep_mu * el_sfs_reco[0] * el_sfs_id[0] * (el_trg_sf - el_trg_unc)



            if isMC and pass_el_all:
                el_sfs_reco, el_sfs_id = lepton_electron_SF(ev.lep_alliso_p4[0].eta(), ev.lep_alliso_p4[0].pt())
                el_trg_sf, el_trg_unc = lepton_electron_trigger_SF(ev.lep_alliso_p4[0].eta(), ev.lep_alliso_p4[0].pt())
                # on 0 position is the value, on 1 is uncertainty

                # test how much leptons SFs affect
                weight_lepall_pu = el_trg_sf * el_sfs_reco[0] * el_sfs_id[0] * weight_pu

                #weight_lepEL_trg = el_trg_sf
                #weight_lepEL_id  = el_sfs_reco[0] * el_sfs_id[0]


            if isMC and pass_el:
                el_sfs_reco, el_sfs_id = lepton_electron_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt())
                el_trg_sf, el_trg_unc = lepton_electron_trigger_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt())
                # on 0 position is the value, on 1 is uncertainty

                weight_lep_pu_elIDUp   = el_trg_sf * (el_sfs_reco[0] + el_sfs_reco[1]) * (el_sfs_id[0] + el_sfs_id[1]) * weight_pu
                weight_lep_pu_elIDDown = el_trg_sf * (el_sfs_reco[0] - el_sfs_reco[1]) * (el_sfs_id[0] - el_sfs_id[1]) * weight_pu

                weight_lep_pu_elTRGUp   = (el_trg_sf + el_trg_unc) * el_sfs_reco[0] * el_sfs_id[0] * weight_pu
                weight_lep_pu_elTRGDown = (el_trg_sf - el_trg_unc) * el_sfs_reco[0] * el_sfs_id[0] * weight_pu

                # test how much leptons SFs affect
                weight_lep_pu = el_trg_sf * el_sfs_reco[0] * el_sfs_id[0] * weight_pu

                weight_lep_pu_PUUp   = el_trg_sf * el_sfs_reco[0] * el_sfs_id[0] * weight_pu_up
                weight_lep_pu_PUDown = el_trg_sf * el_sfs_reco[0] * el_sfs_id[0] * weight_pu_dn

                weight_lep_EL_trg = el_trg_sf
                weight_lep_EL_id  = el_sfs_reco[0] * el_sfs_id[0]


            if isMC and pass_elel:
                el_sfs_reco1, el_sfs_id1 = lepton_electron_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt())
                el_sfs_reco2, el_sfs_id2 = lepton_electron_SF(ev.lep_p4[1].eta(), ev.lep_p4[1].pt())
                el_reco, el_reco_unc = dilepton_and_sfs(el_sfs_reco1[0], el_sfs_reco1[1], el_sfs_reco2[0], el_sfs_reco2[1])
                el_id,   el_id_unc   = dilepton_and_sfs(el_sfs_id1[0], el_sfs_id1[1], el_sfs_id2[0], el_sfs_id2[1])

                el_trg_sf1, el_trg_unc2 = lepton_electron_trigger_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt())
                el_trg_sf2, el_trg_unc2 = lepton_electron_trigger_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt())
                el_trg_sf, el_trg_unc = dilepton_or_sfs(el_trg_sf1, el_trg_unc2, el_trg_sf2, el_trg_unc2)

                # on 0 position is the value, on 1 is uncertainty
                weight_lep_pu_elIDUp   = el_trg_sf * (el_reco + el_reco_unc) * (el_id + el_id_unc) * weight_pu
                weight_lep_pu_elIDDown = el_trg_sf * (el_reco - el_reco_unc) * (el_id - el_id_unc) * weight_pu

                weight_lep_pu_elTRGUp   = (el_trg_sf + el_trg_unc) * el_reco * el_id * weight_pu
                weight_lep_pu_elTRGDown = (el_trg_sf - el_trg_unc) * el_reco * el_id * weight_pu

                # test how much leptons SFs affect
                weight_lep_pu = el_trg_sf * el_reco * el_id * weight_pu

                weight_lep_pu_PUUp   = el_trg_sf * el_reco * el_id * weight_pu_up
                weight_lep_pu_PUDown = el_trg_sf * el_reco * el_id * weight_pu_dn

                weight_lep_EL_trg = el_trg_sf
                weight_lep_EL_id  = el_reco * el_id

            weight_init *= weight_lep_pu if weight_lep_pu > 0.0001 else 0.

        control_counters.Fill(52)

        # -------------------- LEPTONS
        leps     = (ev.lep_p4,        ev.lep_relIso,        ev.lep_matching_gen,        ev.lep_matching_gen_dR,        ev.lep_id)
        leps_all = (ev.lep_alliso_p4, ev.lep_alliso_relIso, ev.lep_alliso_matching_gen, ev.lep_alliso_matching_gen_dR, ev.lep_alliso_id)
        #LeptonSelection(iso=, alliso=)

        # -------------------- TAUS

        # tau pt-s
        # ES correction
        # modes have different correction but the same uncertainty = +- 1.2% = 0.012
        # uncertainties are not correlated, but I'll do correlated UP/DOWN -- all modes UP or all modes DOWN
        #tau_pts_corrected = []
        #tau_pts_corrected_up = []
        #tau_pts_corrected_down = []
        #Mt_tau_met_nominal, Mt_tau_met_up, Mt_tau_met_down = None, None, None
        Mt_tau_met_nominal, Mt_tau_met_up, Mt_tau_met_down = 0, 0, 0

        # so, actually what I need from taus is
        # whether there is a medium tau with pt 30, eta 2.4
        # and if it is OS with the lepton
        # usually it is the tau on first position (0)
        # should I really loop?

        # each selection requires only 1 ID level of tau:
        # lowest and cuts -- tight
        # old -- medium
        # I also need preselection for WJets SS and preselection SS (for QCD)
        taus_main       = []
        taus_candidates = []
        taus_candidates_alliso = []

        # each tau is
        # p4, ES factor (ID lev is per collection), pdg ID (for OS, SS)
        # only top pt tau is treated but that's fine
        '''IDlev
        1 VLoose
        2 Loose
        3 Medium
        4 Tight
        5 VTight
        '''

        #tlep_p4 = TLorentzVector(ev.lep_p4[0].X(), ev.lep_p4[0].Y(), ev.lep_p4[0].Z(), ev.lep_p4[0].T())
        #if ev.tau_p4.size() > 0 and ev.tau_IDlev[0] > 1 and abs(ev.tau_p4[0].eta()) < 2.4:
        for i, (p4, tau_ID, tau_DM, tau_pdgID) in enumerate(zip(ev.tau_p4, ev.tau_IDlev, ev.tau_decayMode, ev.tau_id)):
            # discard taus if they match to lepton
            # TODO: turn it on in next NT
            #if ev.tau_matching_lep[i]:
            #    continue
            # the original index i is needed for the match to refitted values

            if abs(p4.eta()) > 2.4: continue
            # it should work like Python does and not copy these objects! (cast)
            #p4, DM, IDlev = ev.tau_p4[0], ev.tau_decayMode[0], ev.tau_IDlev[0]
            #if IDlev < 3 or abs(p4.eta()) < 2.4: continue # only Medium taus

            ## check dR to lepton
            # it's done in Ntupler -- not needed here
            #ttau_p4 = TLorentzVector(tau_p4.X(), tau_p4.Y(), tau_p4.Z(), tau_p4.T())
            #if not tlep_p4.DeltaR(ttau_p4) > 0.3:
            #    continue

            tau_pt      = p4.pt()
            tau_pt_up   = tau_pt
            tau_pt_down = tau_pt
            TES_factor    = 1. # this factor exists for nominal taus
            TES_factor_up = 1.
            TES_factor_dn = 1.

            jetmatched = ev.tau_dR_matched_jet[i] # > -1

            match_lep        = ev.tau_matching_lep[i]
            match_lep_alliso = ev.tau_matching_allIso_lep[i]

            # option of passing only 3PI taus with high significance
            if ONLY_3PI_TAUS:
                # (unnecessary if, but readable)
                #if not (ev.tau_hasSecondaryVertex[i] and ev.tau_refited_index[i] > -1 and ev.tau_SV_geom_flightLenSign[ev.tau_refited_index[i]] > SV_SIGN_CUT): continue
                if not (ev.tau_refited_index[i] > -1 and ev.tau_SV_geom_flightLenSign[ev.tau_refited_index[i]] > SV_SIGN_CUT): continue

            pass_pt = False
            # MC corrections
            if isMC: # and with_TauES_sys:
                if tau_DM == 0:
                    TES_factor = 0.995
                    TES_factor_up = 0.995 + 0.012
                    TES_factor_dn = 0.995 - 0.012
                elif tau_DM < 10:
                    TES_factor = 1.011
                    TES_factor_up = 1.011 + 0.012
                    TES_factor_dn = 1.011 - 0.012
                else:
                    TES_factor = 1.006
                    TES_factor_up = 1.006 + 0.012
                    TES_factor_dn = 1.006 - 0.012
                tau_pt_up   *= TES_factor_up
                tau_pt_down *= TES_factor_dn
                tau_pt      *= TES_factor
                pass_pt = any(pt_var > TAUS_PT_CUT for pt_var in (tau_pt, tau_pt_up, tau_pt_down))
            else:
                pass_pt = tau_pt > TAUS_PT_CUT

            # I store only 1 ID per selection
            # tight for lowest and cuts
            # medium for old
            # presel for presel

            # nominals
            if tau_pt > 20. and not match_lep:
                taus_candidates.append((p4, TES_factor, tau_pdgID, i, tau_ID, jetmatched))

            if tau_pt > 20. and not match_lep_alliso:
                taus_candidates_alliso.append((p4, TES_factor, tau_pdgID, i, tau_ID, jetmatched))

            if pass_pt and tau_ID > TAUS_ID_CUT and not match_lep:
                taus_main.append((p4, (TES_factor, TES_factor_up, TES_factor_dn), tau_pdgID, i, tau_ID, jetmatched))

            # all taus to met?
            if isMC and PROP_TAU and tau_pt > TAUS_PT_CUT:
                proc_met         -= p4 * (TES_factor - 1.)
                proc_met_JERUp   -= p4 * (TES_factor - 1.)
                proc_met_JERDown -= p4 * (TES_factor - 1.)
                proc_met_JESUp   -= p4 * (TES_factor - 1.)
                proc_met_JESDown -= p4 * (TES_factor - 1.)

            if isMC and tau_pt_up > TAUS_PT_CUT:
                proc_met_TESUp   -=  p4 * (TES_factor_up - 1.)
            if isMC and tau_pt_down > TAUS_PT_CUT:
                proc_met_TESDown -=  p4 * (TES_factor_dn - 1.)

        # sort by IDlev and pt
        taus_candidates_alliso.sort(key=lambda it: (it[4], it[1]   *it[0].pt()), reverse=True)
        taus_candidates       .sort(key=lambda it: (it[4], it[1]   *it[0].pt()), reverse=True)
        taus_main             .sort(key=lambda it: (it[4], it[1][0]*it[0].pt()), reverse=True)
        ## sort by pt only
        #taus_candidates_alliso.sort(key=lambda it: it[1]   *it[0].pt(), reverse=True)
        #taus_candidates       .sort(key=lambda it: it[1]   *it[0].pt(), reverse=True)
        #taus_main             .sort(key=lambda it: it[1][0]*it[0].pt(), reverse=True)


        # ---------- JETS

        # jets save (p4, factor) and (..., dR_far_of_tau) for adv selection -- counting tau-cleaned jets for Loose2
        # TODO: add jet b-discr or b-ID, also consider adding dR to tau (to which tau? loose? medium? save the ID of dR-ed tau?)
        ##SystematicJets = namedtuple('Jets', 'nom sys_JERUp sys_JERDown sys_JESUp sys_JESDown sys_bUp sys_bDown')
        #jets_all_optimized          = SystematicJets('nom'=[], sys_JERUp=[], sys_JERDown=[], sys_JESUp=[], sys_JESDown=[], sys_bUp=[], sys_bDown=[])
        #jets_all_optimized_cuts     = SystematicJets('nom'=[], sys_JERUp=[], sys_JERDown=[], sys_JESUp=[], sys_JESDown=[], sys_bUp=[], sys_bDown=[])
        #jets_all_optimized_old_cuts = SystematicJets('nom'=[], sys_JERUp=[], sys_JERDown=[], sys_JESUp=[], sys_JESDown=[], sys_bUp=[], sys_bDown=[])
        ## each systematic is a list of jets with
        ## p4, energy factor, b SF weight, b ID level

        # each selection needs b-jets (loose and medium) and other jets
        # -- for lj calculations and b-jet control
        # taumatched jets are needed jor lj too
        # they are split in b and non-b jets according to the selection
        '''
        jets_nom     = JetCutsPerSystematic(lowest=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])),
                                              cuts=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])),
                                               old=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])))
        jets_JERUp   = JetCutsPerSystematic(lowest=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])),
                                              cuts=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])),
                                               old=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])))
        jets_JERDown = JetCutsPerSystematic(lowest=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])),
                                              cuts=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])), 
                                               old=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])))
        jets_JESUp   = JetCutsPerSystematic(lowest=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])),
                                              cuts=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])), 
                                               old=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])))
        jets_JESDown = JetCutsPerSystematic(lowest=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])),
                                              cuts=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])), 
                                               old=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])))
        jets_bUp     = JetCutsPerSystematic(lowest=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])),
                                              cuts=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])), 
                                               old=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])))
        jets_bDown   = JetCutsPerSystematic(lowest=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])),
                                              cuts=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])), 
                                               old=JetBSplit(medium=[], loose=[], rest=[], taumatched=([], [])))
        # maybe it's worthwhile to make a custom object with propper defaults and fancy initialization without code
        '''

        jets_nom     = JetBSplit(medium=[], loose=[], rest=[], taumatched=([], []), lepmatched=[])
        jets_JERUp   = JetBSplit(medium=[], loose=[], rest=[], taumatched=([], []), lepmatched=[])
        jets_JERDown = JetBSplit(medium=[], loose=[], rest=[], taumatched=([], []), lepmatched=[])
        jets_JESUp   = JetBSplit(medium=[], loose=[], rest=[], taumatched=([], []), lepmatched=[])
        jets_JESDown = JetBSplit(medium=[], loose=[], rest=[], taumatched=([], []), lepmatched=[])
        jets_bUp     = JetBSplit(medium=[], loose=[], rest=[], taumatched=([], []), lepmatched=[])
        jets_bDown   = JetBSplit(medium=[], loose=[], rest=[], taumatched=([], []), lepmatched=[])
        # maybe it's worthwhile to make a custom object with propper defaults and fancy initialization without code


        #bjetsL_nom     = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsL_JERUp   = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsL_JERDown = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsL_JESUp   = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsL_JESDown = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsL_bUp     = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsL_bDown   = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsM_nom     = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsM_JERUp   = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsM_JERDown = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsM_JESUp   = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsM_JESDown = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsM_bUp     = JetCutsPerSystematic(lowest=[], cuts=[], old=[])
        #bjetsM_bDown   = JetCutsPerSystematic(lowest=[], cuts=[], old=[])

        # I forgot about the needed b-jet splitting (for lj and easy 1M, 2L1M selections)
        # and it got out of hand
        # I don't want to use dict-s of strings for performance sake
        # and trying to operate on direct pointers with named tuples (hopefully that's how it is implemented TODO check speed)
        # macro would be nice here
        # TODO is there another way?

        if len(ev.jet_b_discr) > 0:
            lead_jet_b_discr = ev.jet_b_discr[0]

        b_tag_wp_loose  = 0.5426 # 0.460 # 
        b_tag_wp_medium = 0.8484 # 0.8   # 
        b_tag_wp_tight  = 0.9535 # 0.935 # 
        b_tag_wp = b_tag_wp_medium

        tau_dR_jet_min = 0.4 # ak4 jets

        """
        # lists of "light jets", not passing b-tag
        jets_nominal = [] # nominal jet pts
        jets_jer_up, jets_jer_down = [], []
        jets_jes_up, jets_jes_down = [], []
        # lists of "heavy jets", passing b-tag
        jets_b_nominal = [] # nominal jet pts
        jets_b_jer_up, jets_b_jer_down = [], []
        jets_b_jes_up, jets_b_jes_down = [], []

        #has_2_loose_bjets
        jets_b_nominal_loose = [] # nominal jet pts
        jets_b_nominal_tight = [] # nominal jet pts
        jets_not_loose_b = [] # nominal jet pts
        jets_not_tight_b = [] # nominal jet pts

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

        ## for tt->elmu FAKERATES [TODO: add this as some ifs]
        #all_jets_for_fakes = []
        """

        # these collections should be cross-cleaned from loose+ tau
        #if taus_nominal_min:
        #    tau_p4 = taus_nominal_min[0][0]
        #    ttau_p4 = TLorentzVector(tau_p4.X(), tau_p4.Y(), tau_p4.Z(), tau_p4.T())
        # I'll do cross-clean with nominal loose taus
        #tausL_nom    = TauCutsPerSystematic(lowest=[], cuts=[], old=[])
        # of lowest cuts? -- no, of the same as in the selection
        #                 -- and this again messes up the code.......

        ##taus_nom_TLorentz    = TauCutsPerSystematic(lowest=[], loose=[], cuts=[], old=[], oldVloose=[], presel=[], presel_alliso=[]) # presel won't be needed
        #taus_nom_TLorentz    = TauSplit(candidates=[], medium=[])
        ## presel should also be cleaned of the tau -- of the one I do SS with
        ## in fact the cross-cleaning should be done only from the used tau
        ## I can sort the taus by pt and use the first one

        ## just lists of tlorentz vectors
        ## in lowest and cuts we use tight tau
        ## -- using only the first tau in the list -- the target one!
        #for tau_p4, _, _, _, _ in taus_candidates[:1]:
        #    taus_nom_TLorentz.candidates.append(TLorentzVector(tau_p4.X(), tau_p4.Y(), tau_p4.Z(), tau_p4.T()))

        ## in old we use medium tau
        #for tau_p4, _, _, _, _ in taus_main[:1]:
        #    taus_nom_TLorentz.medium.append(TLorentzVector(tau_p4.X(), tau_p4.Y(), tau_p4.Z(), tau_p4.T()))
        ##    taus_nom_TLorentz.oldVloose.append(TLorentzVector(tau_p4.X(), tau_p4.Y(), tau_p4.Z(), tau_p4.T()))
        ## in principle I should also change TES -- but this effect is very second order

        # I need to separate b-tagged jets and tau candidates
        # -- I need to do it in main selection with tau passing ID lev
        #    and in alliso selection, with the tau candidate
        # loose and tight b jet collections are cross cleaned from loose zero tau
        # only 1 jet is taken out
        #tau_match_lowest, tau_match_cuts = False, False
        taus_main_TLorentz = []
        taus_candidates_alliso_TLorentz = []

        for tau_p4, _, _, _, _, _ in taus_main[:1]:
            taus_main_TLorentz.append(TLorentzVector(tau_p4.X(), tau_p4.Y(), tau_p4.Z(), tau_p4.T()))
        for tau_p4, _, _, _, _, _ in taus_candidates_alliso[:1]:
            taus_candidates_alliso_TLorentz.append(TLorentzVector(tau_p4.X(), tau_p4.Y(), tau_p4.Z(), tau_p4.T()))

        tau_match_medium = False
        tau_match_candidate_alliso = False

        # transition to v20-v21
        #miniaod_jets = ev.jet_p4 if OLD_MINIAOD_JETS else ev.jet_initial_p4

        # MC jets are corrected for JER, which is propagated to met
        #jet_cor

        N_jets_nom_med     = 0 # b-tagged jets not matching medium tau
        N_jets_nom_med_all = 0 # all b-taged jets, regardless of taus
        N_jets_nom_all     = 0
        N_jets_JERUp_med   = 0
        N_jets_JERUp_all   = 0
        N_jets_JERDown_med = 0
        N_jets_JERDown_all = 0
        N_jets_JESUp_med   = 0
        N_jets_JESUp_all   = 0
        N_jets_JESDown_med = 0
        N_jets_JESDown_all = 0

        N_jets_nom_med_alliso = 0
        N_jets_nom_all_alliso = 0

        # -- I miss the bSF weight
        # TODO: in principle it should be added according to the b-s used in the selection and passed down to the channel
        #       notice the b SF are calculated in old scheme, for medium jets only
        #       so for now I'll use this weight everywhere
        weight_bSF     = 1.
        weight_bSFUp   = 1.
        weight_bSFDown = 1.
        weight_bSF_JERUp   = 1.
        weight_bSF_JERDown = 1.
        weight_bSF_JESUp   = 1.
        weight_bSF_JESDown = 1.

        sub_lep = True
        for i in xrange(ev.jet_p4.size()):
            # discard jets if they match to lepton
            # TODO: turn it on in next NT
            #if ev.jet_matching_lep[i]:
            #    continue

            #pfid, p4 = ev.jet_PFID[i], ev.jet_p4[i]
            # RECORRECTED jets
            #p4 = ev.jet_uncorrected_p4[i]
            # jet_p4 is fully corrected online
            # and whole correction is propagated to corrected met
            # as I did before I use the subset of passed jets, and propagate only their correction to met
            # TODO: add comparison to using all-corrected jets
            #p4 = ev.jet_initial_p4[i]
            #p4 = miniaod_jets[i]

            # nominally corrected jets in NTuple (v28 data and all previous for MC)
            #p4 = ev.jet_p4[i]

            # the default miniaod objects
            p4 = ev.jet_initial_p4[i]
            #p4 = ev.jet_uncorrected_p4[i]

            pfid = ev.jet_PFID[i]

            # actually the jet_p4 is the slimmed jet
            # -- works for slimmed met, which is met_init
            # met_x = ev.met_init.Px()

            # to apply JES I need uncorrected jet and to propagate it to MET I need the uncorrection factor included
            # = main thing is to start from the same point

            # but now I need full factor for MC
            if pfid < 1 or abs(p4.eta()) > 2.5: continue # Loose PFID and eta

            jet_b_discr = ev.jet_b_discr[i]
            #all_jets_b_discrs.append(jet_b_discr)

            #match_lep        = ev.jet_matching_lep[i]
            match_lep_alliso = ev.jet_matching_allIso_lep[i]
            match_lep        = ev.jet_matching_lep[i]

            # jet energy scale factor = JES factor * JER
            # But the JES is re-correction
            # Thus also * UncorFactor
            en_factor = 1.

            jet_factor_JERUp   = 1.
            jet_factor_JERDown = 1.

            jet_factor_JESUp   = 1.
            jet_factor_JESDown = 1.

            #jes_factor = ev.jet_jes_recorrection[i] # ALREADY APPLIED IN NTUPLER
            #en_factor *= jes_factor
            # the idea is to match these jets with the full_corr met

            HF = -1
            PF = -1
            jet_index = i
            #genmatch = 0
            jes_uncorFactor = ev.jet_uncorrected_jecFactor[i]
            if isMC:
                # JES is ALREADY APPLIED
                #en_factor *= jer_factor * jes_uncorFactor

                # jer factor is needed for the syst variation
                jer_factor = ev.jet_jer_factor[i]
                # and let's try the new prescription: miniaod objects + jer
                en_factor *= jer_factor
                # this jer is calculated on top of recorected jes, let's hope it's fine
                # (mostly depends on eta, dependence on pt is small)
                # and the jes variation (jes_shift used in the following)
                # is calculated in the recorrection procedure
                # -- let's hope it is also adequate
                HF = ev.jet_hadronFlavour[i]
                PF = ev.jet_partonFlavour[i]

                # JES must be reapplied for TT PS systematic datasets
                if isTT_systematic and ev.genjet_matched[i] and correction_jes_barrel_b:
                    pt_gen  = ev.genjet_pt[i]
                    #jes_uncorFactor = ev.jet_uncorrected_jecFactor[i]
                    jet_abs_eta = abs(p4.eta())
                    if HF == 5: # by hadron flavour it's a b jet
                        if jet_abs_eta < 1.5: # barrel
                            ps_jes_factor = correction_jes_barrel_b.GetBinContent(correction_jes_barrel_b.FindBin(pt_gen))
                        else:                 # endcap
                            ps_jes_factor = correction_jes_endcap_b.GetBinContent(correction_jes_endcap_b.FindBin(pt_gen))
                    elif HF == 4:
                        if jet_abs_eta < 1.5:
                            ps_jes_factor = correction_jes_barrel_c.GetBinContent(correction_jes_barrel_c.FindBin(pt_gen))
                        else:
                            ps_jes_factor = correction_jes_endcap_c.GetBinContent(correction_jes_endcap_c.FindBin(pt_gen))
                    else:
                        if jet_abs_eta < 1.5:
                            ps_jes_factor = correction_jes_barrel_udsg.GetBinContent(correction_jes_barrel_udsg.FindBin(pt_gen))
                        else:
                            ps_jes_factor = correction_jes_endcap_udsg.GetBinContent(correction_jes_endcap_udsg.FindBin(pt_gen))

                    # in the result we apply correction to the nominal JES in MC
                    # the correction is ratio of custom JES measured in nominal and PS TTbar MC
                    # so all possible bias due to the custom measurement method cancels out
                    en_factor *= ps_jes_factor

                # this check is done online with dR 0.4
                #if ev.jet_matching_gen_dR[i] < 0.3:
                    #genmatch = ev.jet_matching_gen[i]

            if sub_lep and match_lep:
                # sub jet by lep
                # works only for 1 lepton
                #proc_met_lepsub -= ev.lep_p4[0] - p4
                proc_met_lepsub += p4
                sub_lep = False
            else:
                proc_met_lepsub -= p4 * (en_factor - 1.)

            ## propagate the correction to met
            #proc_met -= p4 * (en_factor - 1.)

            ## and with systematic variations
            #if isMC and with_JER_sys:
            #    #jer_factor, up, down = ev.jet_jer_factor[i], ev.jet_jer_factor_up[i], ev.jet_jer_factor_down[i]
            #    jer_up, jer_down = ev.jet_jer_factor_up[i], ev.jet_jer_factor_down[i]
            #    proc_met_JERUp   -= p4 * (jer_up - 1.)
            #    proc_met_JERDown -= p4 * (jer_down - 1.)

            #if isMC and with_JES_sys:
            #    #jes_shift = ev.jet_jes_correction_relShift[i]
            #    jes_shift = ev.jet_jes_uncertainty[i]
            #    proc_met_JESUp   -= p4 * (en_factor*(1 + jes_shift) - 1.)
            #    proc_met_JESDown -= p4 * (en_factor*(1 - jes_shift) - 1.)


            # nominals
            jet_pt_init = p4.pt()

            # JER are multiplied by the separate JER factor
            jet_pt_JERUp   = jet_pt_init
            jet_pt_JERDown = jet_pt_init

            # nominal and JES are multiplied by the nominal JER
            jet_pt  = jet_pt_init * en_factor
            jet_eta = p4.eta()

            jet_pt_JESUp   = jet_pt
            jet_pt_JESDown = jet_pt

            possible_jet_pts = [jet_pt]

            #if p4.pt() > 30: # nominal jet
            # TODO: for optimization tests I reduced the cut, review it with the test results
            if jet_pt > 20: # nominal jet
                #for ttau_p4 in taus_nom_TLorentz.lowest:
                #    if tj_p4.DeltaR(ttau_p4) < 0.3:
                #        tau_match_lowest = True
                #        break

                #for ttau_p4 in taus_nom_TLorentz.cuts:
                #    if tj_p4.DeltaR(ttau_p4) < 0.3:
                #        tau_match_cuts = True
                #        break

                #for ttau_p4 in taus_nom_TLorentz.old:
                #    if tj_p4.DeltaR(ttau_p4) < 0.3:
                #        tau_match_medium = True
                #        break

                jet_tau_match_old    = False
                jet_tau_match_alliso = False
                if not tau_match_medium or not tau_match_candidate_alliso:
                    tj_p4 = TLorentzVector(p4.X(), p4.Y(), p4.Z(), p4.T())

                if not tau_match_medium:
                    tau_match_medium    = jet_tau_match_old    = any(tj_p4.DeltaR(ttau) < tau_dR_jet_min for ttau in taus_main_TLorentz)
                if not tau_match_candidate_alliso:
                    tau_match_candidate_alliso = jet_tau_match_alliso = any(tj_p4.DeltaR(ttau) < tau_dR_jet_min for ttau in taus_candidates_alliso_TLorentz)

                # jet passed ID, eta, lowest pt threshold and doesn't match to a Loose tau
                # now find its' systematic corrections and save
                flavId = ev.jet_hadronFlavour[i]
                b_tagged = b_tagged_medium = jet_b_discr > b_tag_wp_medium
                b_tagged_loose = jet_b_discr > b_tag_wp_loose
                b_tagged_tight = jet_b_discr > b_tag_wp_tight

                bID_lev = sum((b_tagged_loose, b_tagged_medium, b_tagged_tight))

                # jet at given systematic and selection is defined by p4, energy factor, bSF weight and b ID lev
                # JER and JES are already nominally corrected -- only Up/Down are needed
                jet_factor_JERUp = 1.
                jet_factor_JESUp = 1.
                jet_factor_JERDown = 1.
                jet_factor_JESDown = 1.
                # b SF needs to be calculated
                jet_weight_bSF_nom = 1.     # there might be several schemes for b SF weights
                jet_weight_bSFUp   = 1.
                jet_weight_bSFDown = 1.

                # so, for nominal jets I just need to check is bSF needs to be calculated and add it if so

                if isMC:
                    if with_bSF:
                        # for now I calculate only the old b SF scheme with only medium b WP
                        # TODO: add the 2L1M b SF scheme
                        jet_weight_bSF_nom, _, _ = calc_btag_sf_weight(b_tagged_medium, flavId, jet_pt, jet_eta)
                        if with_bSF_sys:
                            # again only old scheme
                            jet_weight_bSFUp  , _, _ = calc_btag_sf_weight(b_tagged, flavId, jet_pt, jet_eta, "up")
                            jet_weight_bSFDown, _, _ = calc_btag_sf_weight(b_tagged, flavId, jet_pt, jet_eta, "down")

                    if with_JER_sys:
                        #jer_factor, up, down = ev.jet_jer_factor[i], ev.jet_jer_factor_up[i], ev.jet_jer_factor_down[i]
                        up, down = ev.jet_jer_factor_up[i], ev.jet_jer_factor_down[i]
                        jet_factor_JERUp   = up    if jer_factor > 0 else 0
                        jet_factor_JERDown = down  if jer_factor > 0 else 0

                        #jet_pt_JERUp   = jet_pt * jet_factor_JERUp
                        #jet_pt_JERDown = jet_pt * jet_factor_JERDown
                        # --- this is correct, jet_pt is * jer_factor * jer_up / jer_factor
                        #     the jer_up factor must be * jer_nom
                        # for consistency let's make sys jet en corrections work like the nominal
                        jet_pt_JERUp    = jet_pt_init * jet_factor_JERUp
                        jet_pt_JERDown  = jet_pt_init * jet_factor_JERDown
                        possible_jet_pts.extend([jet_pt_JERUp, jet_pt_JERDown])

                    if with_JES_sys:
                        #jes_shift = ev.jet_jes_correction_relShift[i]
                        jes_shift = ev.jet_jes_uncertainty[i]
                        # for correct propagation to MET
                        # the factor of jet systematic must inlcude full correction from the initial jet
                        jet_factor_JESUp   = en_factor * (1 + jes_shift)
                        jet_factor_JESDown = en_factor * (1 - jes_shift)

                        #jet_pt_JESUp   = jet_pt * jet_factor_JESUp
                        #jet_pt_JESDown = jet_pt * jet_factor_JESDown
                        # --- these are correct, jet_pt = jet * jer_factor * jes_up
                        jet_pt_JESUp    = jet_pt_init * jet_factor_JESUp
                        jet_pt_JESDown  = jet_pt_init * jet_factor_JESDown
                        possible_jet_pts.extend([jet_pt_JESUp, jet_pt_JESDown])

                en_factors  = (en_factor, jet_factor_JERUp, jet_factor_JERDown, jet_factor_JESUp, jet_factor_JESDown)
                bSF_weights = (jet_weight_bSF_nom, jet_weight_bSFUp, jet_weight_bSFDown)

                pass_jet_pt_cut = any(a_pt > JETS_PT_CUT for a_pt in possible_jet_pts)
                if not pass_jet_pt_cut:
                    continue

                # match lep alliso for alliso study selection
                if not match_lep_alliso:
                    # only presel is relevant, therefore don't count tau dR
                    if jet_pt > JETS_PT_CUT:
                        N_jets_nom_all_alliso += 1
                        # emulate preselection, don't count the tau candidate
                        #if b_tagged_medium:
                        # no, there is a basic preselection without tau candidate and the special preselection with the candidate
                        if b_tagged_medium and not jet_tau_match_alliso:
                            N_jets_nom_med_alliso += 1

                # correct the met from lep-matched jet
                if match_lep:

                  if PROP_UNCORLEPJET:
                    # done for data and MC
                    #if jet_pt > JETS_PT_CUT:
                    # uncorrect the lep jet, return the correction to MET
                    proc_met -= p4 * (jes_uncorFactor - 1.)
                    # also propagate to TES variations
                    proc_met_TESUp   -= p4 * (jes_uncorFactor - 1.)
                    proc_met_TESDown -= p4 * (jes_uncorFactor - 1.)
                    #if jet_pt_JERUp > JETS_PT_CUT:
                    proc_met_JERUp   -= p4 * (jes_uncorFactor - 1.)
                    #if jet_pt_JERDown > JETS_PT_CUT:
                    proc_met_JERDown -= p4 * (jes_uncorFactor - 1.)
                    #if jet_pt_JESUp > JETS_PT_CUT:
                    proc_met_JESUp   -= p4 * (jes_uncorFactor - 1.)
                    #if jet_pt_JESDown > JETS_PT_CUT:
                    proc_met_JESDown -= p4 * (jes_uncorFactor - 1.)

                  if REMOVE_LEPJET:
                    # done for data and MC
                    # the idea: the lep jet is a phantom object, which does not represent any real part of the decay
                    #           therefore the MET must not be changed due to some corrections of this jet or by the jet itself
                    #if jet_pt > JETS_PT_CUT:
                    # uncorrect the lep jet, return the correction to MET
                    # lepjet is a phantom 
                    proc_met += p4
                    # also propagate to TES variations
                    proc_met_TESUp   += p4
                    proc_met_TESDown += p4
                    #if jet_pt_JERUp > JETS_PT_CUT:
                    proc_met_JERUp   += p4
                    #if jet_pt_JERDown > JETS_PT_CUT:
                    proc_met_JERDown += p4
                    #if jet_pt_JESUp > JETS_PT_CUT:
                    proc_met_JESUp   += p4
                    #if jet_pt_JESDown > JETS_PT_CUT:
                    proc_met_JESDown += p4

                  if isMC and PROP_LEPJET:
                    # this does not make sense:
                    # to correct a non existent jet and propagate the correction to MET in MC only
                    # -- there is no physical reason for this change of MET
                    # also it is done only if the jet passes the pT threshold?? what's the reason for this?
                    if jet_pt > JETS_PT_CUT:
                        proc_met -= p4 * (en_factor - 1.)
                        # also propagate to TES variations
                        proc_met_TESUp   -= p4 * (en_factor - 1.)
                        proc_met_TESDown -= p4 * (en_factor - 1.)
                    if jet_pt_JERUp > JETS_PT_CUT:
                        proc_met_JERUp   -= p4 * (jet_factor_JERUp - 1.)
                    if jet_pt_JERDown > JETS_PT_CUT:
                        proc_met_JERDown -= p4 * (jet_factor_JERDown - 1.)
                    if jet_pt_JESUp > JETS_PT_CUT:
                        proc_met_JESUp   -= p4 * (jet_factor_JESUp - 1.)
                    if jet_pt_JESDown > JETS_PT_CUT:
                        proc_met_JESDown -= p4 * (jet_factor_JESDown - 1.)
                    # these were broken! the nominal jer factor was not applied, but now it is.

                  if isMC and PROP_LEPJET_UNCLUSTER:
                    # return jes_cor to met
                    # in both data and MC
                    jes_uncorFactor = ev.jet_uncorrected_jecFactor[i]
                    if jet_pt > JETS_PT_CUT:
                        proc_met         -= p4 * (jes_uncorFactor - 1.)
                        proc_met_TESUp   -= p4 * (jes_uncorFactor - 1.)
                        proc_met_TESDown -= p4 * (jes_uncorFactor - 1.)
                    if jet_pt_JERUp > JETS_PT_CUT:
                        proc_met_JERUp   -= p4 * (jes_uncorFactor - 1.)
                    if jet_pt_JERDown > JETS_PT_CUT:
                        proc_met_JERDown -= p4 * (jes_uncorFactor - 1.)
                    if jet_pt_JESUp > JETS_PT_CUT:
                        proc_met_JESUp   -= p4 * (jes_uncorFactor - 1.)
                    if jet_pt_JESDown > JETS_PT_CUT:
                        proc_met_JESDown -= p4 * (jes_uncorFactor - 1.)

                # match lep
                if not match_lep:

                  # correct the met from all not lep-matched jets
                  if isMC and PROP_JETS:
                      if jet_pt > JETS_PT_CUT:
                          proc_met       -= p4 * (en_factor - 1.)
                          proc_met_TESUp   -= p4 * (en_factor - 1.)
                          proc_met_TESDown -= p4 * (en_factor - 1.)
                      if jet_pt_JERUp > JETS_PT_CUT:
                          proc_met_JERUp   -= p4 * (jet_factor_JERUp - 1.)
                      if jet_pt_JERDown > JETS_PT_CUT:
                          proc_met_JERDown -= p4 * (jet_factor_JERDown - 1.)
                      if jet_pt_JESUp > JETS_PT_CUT:
                          proc_met_JESUp   -= p4 * (jet_factor_JESUp - 1.)
                      if jet_pt_JESDown > JETS_PT_CUT:
                          proc_met_JESDown -= p4 * (jet_factor_JESDown - 1.)
                      # these were broken! the nominal jer factor was not applied, but now it is.

                  # multiply the b weights for all possible combinations
                  # and correct the met with all except lep-jet
                  if jet_pt > JETS_PT_CUT:
                      weight_bSF     *= jet_weight_bSF_nom
                      weight_bSFUp   *= jet_weight_bSFUp
                      weight_bSFDown *= jet_weight_bSFDown

                  if jet_pt_JERUp > JETS_PT_CUT:
                      weight_bSF_JERUp *= jet_weight_bSF_nom
                  if jet_pt_JERDown > JETS_PT_CUT:
                      weight_bSF_JERDown *= jet_weight_bSF_nom

                  if jet_pt_JESUp > JETS_PT_CUT:
                      weight_bSF_JESUp *= jet_weight_bSF_nom
                  if jet_pt_JESDown > JETS_PT_CUT:
                      weight_bSF_JESDown *= jet_weight_bSF_nom

                  # count the jets
                  if b_tagged_medium:
                      if jet_pt > JETS_PT_CUT:
                          N_jets_nom_med_all += 1
                  if not jet_tau_match_old and b_tagged_medium:
                      if jet_pt > JETS_PT_CUT:
                          N_jets_nom_med += 1
                      if jet_pt_JERUp > JETS_PT_CUT:
                          N_jets_JERUp_med += 1
                      if jet_pt_JERDown > JETS_PT_CUT:
                          N_jets_JERDown_med += 1
                      if jet_pt_JESUp > JETS_PT_CUT:
                          N_jets_JESUp_med += 1
                      if jet_pt_JESDown > JETS_PT_CUT:
                          N_jets_JESDown_med += 1

                  if jet_pt > JETS_PT_CUT:
                      N_jets_nom_all += 1
                  if jet_pt_JERUp > JETS_PT_CUT:
                      N_jets_JERUp_all += 1
                  if jet_pt_JERDown > JETS_PT_CUT:
                      N_jets_JERDown_all += 1
                  if jet_pt_JESUp > JETS_PT_CUT:
                      N_jets_JESUp_all += 1
                  if jet_pt_JESDown > JETS_PT_CUT:
                      N_jets_JESDown_all += 1

                  # save the nominal jets for lj cut etc
                  if not (jet_pt > JETS_PT_CUT):
                      continue
                  # match tau
                  if not jet_tau_match_old:
                    # the same, but the cuts and lowest -> cuts
                    # nominals
                    # match b-tag
                    if b_tagged_medium:
                        jets_nom.medium.append((p4, en_factors, bSF_weights, jet_b_discr, HF, PF, jet_index))

                    elif b_tagged_loose:
                        jets_nom.loose.append((p4, en_factors, bSF_weights, jet_b_discr, HF, PF, jet_index))

                    else:
                        jets_nom.rest.append((p4, en_factors, bSF_weights, jet_b_discr, HF, PF, jet_index))

                  else:
                    # this is the old selection, the medium b jets count
                    if b_tagged_medium:
                        jets_nom.taumatched[0].append((p4, en_factors, bSF_weights, jet_b_discr, HF, PF, jet_index))

                    else:
                        jets_nom.taumatched[1].append((p4, en_factors, bSF_weights, jet_b_discr, HF, PF, jet_index))

                else:
                    # matched lep
                    if jet_pt > JETS_PT_CUT:
                        jets_nom.lepmatched.append((p4, en_factors, bSF_weights, jet_b_discr, HF, PF, jet_index))
                    # TODO: notice, the b weights of these jets are not considered
                    # and these jets are not considered in number of jets
                    #if jet_pt > JETS_PT_CUT:
                    #    N_jets_nom_all += 1
                    #if jet_pt_JERUp > JETS_PT_CUT:
                    #    N_jets_JERUp_all += 1
                    #if jet_pt_JERDown > JETS_PT_CUT:
                    #    N_jets_JERDown_all += 1
                    #if jet_pt_JESUp > JETS_PT_CUT:
                    #    N_jets_JESUp_all += 1
                    #if jet_pt_JESDown > JETS_PT_CUT:
                    #    N_jets_JESDown_all += 1






        # sort jets by pt
        # mainly to record various "leading blah jet"
        # TODO: not sure if it is worthwhile -- shouldn't all this be already sorted in ntupler?

        """ trying without sort
        if isMC:
            if with_bSF_sys:
                jets_bUp   .lowest.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bUp     .cuts.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bUp      .old.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bDown .lowest.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bDown   .cuts.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bDown    .old.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)

                jets_bUp   .lowest.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bUp     .cuts.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bUp      .old.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bDown .lowest.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bDown   .cuts.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bDown    .old.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)

                jets_bUp   .lowest.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bUp     .cuts.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bUp      .old.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bDown .lowest.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bDown   .cuts.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_bDown    .old.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)

            if with_JER_sys:
                jets_JERUp   .lowest.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_JERUp     .cuts.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_JERUp      .old.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_JERDown .lowest.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_JERDown   .cuts.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_JERDown    .old.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)

            if with_JES_sys:
                jets_JESUp   .lowest.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_JESUp     .cuts.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_JESUp      .old.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_JESDown .lowest.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_JESDown   .cuts.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
                jets_JESDown    .old.sort(key=lambda it: it[1]*it[0].pt(), reverse=True)
        """ # trying without sort

        #has_medium_tau = any(IDlev > 2 and p4.pt() > 30 for IDlev, p4 in zip(ev.tau_IDlev, ev.tau_p4))
        #has_medium_tau = ev.tau_IDlev.size() > 0 and ev.tau_IDlev[0] > 2 and ev.tau_p4[0].pt() > 30
        #has_medium_tau = bool(tau_pts_corrected)
        #TODO: propagate TES to MET?

        '''
        # for tau FAKERATES
        # loop over all jets and taus
        # matching them in dR
        # and if matched adding into according collection
        tau_jets_candidates = []
        tau_jets_vloose = []
        tau_jets_loose  = []
        tau_jets_medium = []
        tau_jets_tight  = []
        tau_jets_vtight = []
        for jet_p4 in all_jets_for_fakes:
            for tau_cand_p4, cand_IDlev in zip(ev.tau_p4, ev.tau_IDlev):
                #if tau_cand_p4.eta() < 2.3 and tau_cand_p4.pt() > 20: was done in ntupler
                tj_p4   = TLorentzVector(jet_p4.X(), jet_p4.Y(), jet_p4.Z(), jet_p4.T())
                ttau_p4 = TLorentzVector(tau_cand_p4.X(), tau_cand_p4.Y(), tau_cand_p4.Z(), tau_cand_p4.T())
                if tj_p4.DeltaR(ttau_p4) < 0.3:
                    # add the jet p4 into collection according to IDlev
                    #Int_t IDlev = 0;
                    #if (tau.tauID(tau_VTight_ID)) IDlev = 5;
                    #else if (tau.tauID(tau_Tight_ID))  IDlev = 4;
                    #else if (tau.tauID(tau_Medium_ID)) IDlev = 3;
                    #else if (tau.tauID(tau_Loose_ID))  IDlev = 2;
                    #else if (tau.tauID(tau_VLoose_ID)) IDlev = 1;
                    tau_jets_candidates.append(jet_p4)
                    if cand_IDlev > 0:
                        if cand_IDlev == 1:
                            tau_jets_vloose.append(jet_p4)
                        elif cand_IDlev == 2:
                            tau_jets_vloose.append(jet_p4)
                            tau_jets_loose .append(jet_p4)
                        elif cand_IDlev == 3:
                            tau_jets_vloose.append(jet_p4)
                            tau_jets_loose .append(jet_p4)
                            tau_jets_medium.append(jet_p4)
                        elif cand_IDlev == 4:
                            tau_jets_vloose.append(jet_p4)
                            tau_jets_loose .append(jet_p4)
                            tau_jets_medium.append(jet_p4)
                            tau_jets_tight .append(jet_p4)
                        elif cand_IDlev == 5:
                            tau_jets_vloose.append(jet_p4)
                            tau_jets_loose .append(jet_p4)
                            tau_jets_medium.append(jet_p4)
                            tau_jets_tight .append(jet_p4)
                            tau_jets_vtight.append(jet_p4)
        '''



        # shape systematics are:
        # corrected jet pts and tau pts
        # weight variations
        #nominal systematics
        # jet pts, tau pts, b weight (=1 for data), pu weight (=1 for data)

        syst_weights_nominal = [weight_init, 1., 1.]
        syst_weights = {} # {'NOMINAL': syst_weights_nominal}
        syst_objects = {'NOMINAL': [jets_nom, taus_main, taus_candidates, proc_met]}
        #print syst_objects.keys()
        #print syst_weights.keys()

        # SYSTEMATIC LOOP, RECORD
        # for each systematic
        # pass one of the reco selections
        # check the subprocess
        # store distr

        control_counters.Fill(53)

        #for sys_i, (sys_name, (jets, taus, proc_met)) in enumerate(syst_objects.items()):
        sys_i, (sys_name, (jets, sel_taus, sel_tau_candidates, proc_met)) = 0, syst_objects.items()[0]
        #control_counters.Fill(4 + sys_i)
        '''
        1) for each variation of objects check if which channels pass
        2) then in the recording make a special section to fill histos with systematic weights
        3) for the rest use only the nominal weights

        from the object I only need the weight of bSF?
        routine to record a histo with all sys weights?
        and to create that histo with all sys weights?
        -- already the histos are created only for certain systs
        '''


        #sys_weight_min = weight * weight_bSF_min * weight_PU * weight_top_pt
        # pass reco selections

        # all the channel selections follow
        passed_channels = []

        if sys_i == 0:
            control_counters.Fill(100)

        if pass_mu:
            control_counters.Fill(101)

        if pass_el:
            control_counters.Fill(102)

        if pass_elmu:
            control_counters.Fill(103)

        if pass_elmu_el:
            control_counters.Fill(104)

        if pass_mu_all:
            control_counters.Fill(105)

        if pass_el_all:
            control_counters.Fill(106)

        #has_lowest_2L1M = len(jets.lowest.medium) > 0 and (len(jets.lowest.medium) + len(jets.lowest.loose)) > 1 

        weight_bSF_lowest = 1.
        # p4, energy factor, b SF weight, ID lev
        #for _, _, jet_weight, _, _, _, _ in jets.lowest.medium + jets.lowest.loose + jets.lowest.rest:
        #    weight_bSF_lowest *= jet_weight


        # separate b and tau tags
        old_jet_sel = len(jets.medium) > 0 and (len(jets.taumatched[0]) + len(jets.taumatched[1]) + len(jets.medium) + len(jets.loose) + len(jets.rest)) > 2

        pass_2b = len(jets.medium) > 1

        old_jet_sel_alliso = len(jets.medium) > 0 and (len(jets.taumatched[0]) + len(jets.taumatched[1]) + len(jets.medium) + len(jets.loose) + len(jets.rest)) > 2

        if old_jet_sel:
            control_counters.Fill(111)

        if len(sel_taus) > 0:
            control_counters.Fill(112)

        if pass_mu and old_jet_sel:
            control_counters.Fill(203)

        if pass_el and old_jet_sel:
            control_counters.Fill(204)

        if pass_mu and old_jet_sel and len(sel_taus) > 0:
            control_counters.Fill(304)

        if pass_el and old_jet_sel and len(sel_taus) > 0:
            control_counters.Fill(305)

        # 3 jets (2b and 1 tau) and 1 b-tagged
        # TODO: compare with previous result with only 2 jets and 1 b-tagged -- check njets distrs
        pass_old_mu_sel = pass_mu and old_jet_sel and len(sel_taus) > 0 # and no met cut
        pass_old_el_sel = pass_el and old_jet_sel and len(sel_taus) > 0 # and no met cut

        pass_old_mu_presel = pass_mu and old_jet_sel and len(sel_tau_candidates) > 0
        pass_old_el_presel = pass_el and old_jet_sel and len(sel_tau_candidates) > 0
        pass_old_lep_presel = pass_old_mu_presel or pass_old_el_presel

        pass_elmu_close = pass_elmu and (len(jets.medium) + len(jets.taumatched[0])) > 0 and (len(jets.taumatched[0]) + len(jets.taumatched[1]) + len(jets.medium) + len(jets.loose) + len(jets.rest)) > 1

        pass_elmu_el_close = pass_elmu_el and (len(jets.medium) + len(jets.taumatched[0])) > 0 and (len(jets.taumatched[0]) + len(jets.taumatched[1]) + len(jets.medium) + len(jets.loose) + len(jets.rest)) > 1

        sel_leps, sel_jets = leps, jets_nom

        # final CHANNEL SELECTION decision
        tt_channel_presel_stage = passes_tt_preselection_stages(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_med_all, N_jets_nom_all), taus_candidates, proc_met)

        # nominal
        tt_channel_sel_stage          = passes_tt_selection_stages(passed_triggers, leps, (N_jets_nom_med,     N_jets_nom_all),     [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)

        tt_channel_sel_stage_TESUp    = passes_tt_selection_stages(passed_triggers, leps, (N_jets_nom_med,     N_jets_nom_all),     [tau for tau in sel_taus if (tau[0]*tau[1][1]).pt() > TAUS_PT_CUT], proc_met)
        tt_channel_sel_stage_TESDown  = passes_tt_selection_stages(passed_triggers, leps, (N_jets_nom_med,     N_jets_nom_all),     [tau for tau in sel_taus if (tau[0]*tau[1][2]).pt() > TAUS_PT_CUT], proc_met)
        tt_channel_sel_stage_JERUp    = passes_tt_selection_stages(passed_triggers, leps, (N_jets_JERUp_med,   N_jets_JERUp_all),   [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        tt_channel_sel_stage_JERDown  = passes_tt_selection_stages(passed_triggers, leps, (N_jets_JERDown_med, N_jets_JERDown_all), [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        tt_channel_sel_stage_JESUp    = passes_tt_selection_stages(passed_triggers, leps, (N_jets_JESUp_med,   N_jets_JESUp_all),   [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        tt_channel_sel_stage_JESDown  = passes_tt_selection_stages(passed_triggers, leps, (N_jets_JESDown_med, N_jets_JESDown_all), [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)

        # alliso tt, the presel is most important here
        #tt_channel_presel_stage = passes_tt_selection_stages(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_all), taus_candidates, proc_met)
        tt_channel_stage_alliso = passes_tt_selection_stages_alliso(passed_triggers, leps_all, (N_jets_nom_med_alliso, N_jets_nom_all_alliso), taus_candidates_alliso, proc_met)

        # WJets control selection
        selection_wjets_control = passes_wjets_control_selection(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_all), taus_candidates, proc_met)

        # full syst for DY
        dy_channel_sel_stage         = passes_dy_tautau_selection_stages(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_all), [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        dy_channel_sel_stage_TESUp   = passes_dy_tautau_selection_stages(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_all), [tau for tau in sel_taus if (tau[0]*tau[1][1]).pt() > TAUS_PT_CUT], proc_met)
        dy_channel_sel_stage_TESDown = passes_dy_tautau_selection_stages(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_all), [tau for tau in sel_taus if (tau[0]*tau[1][2]).pt() > TAUS_PT_CUT], proc_met)
        dy_channel_sel_stage_JERUp   = passes_dy_tautau_selection_stages(passed_triggers, leps, (N_jets_JERUp_med,   N_jets_JERUp_all),   [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        dy_channel_sel_stage_JERDown = passes_dy_tautau_selection_stages(passed_triggers, leps, (N_jets_JERDown_med, N_jets_JERDown_med), [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        dy_channel_sel_stage_JESUp   = passes_dy_tautau_selection_stages(passed_triggers, leps, (N_jets_JESUp_med,   N_jets_JESUp_all),   [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        dy_channel_sel_stage_JESDown = passes_dy_tautau_selection_stages(passed_triggers, leps, (N_jets_JESDown_med, N_jets_JESDown_med), [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)

        # and DY mumu
        dy_mumu_channel_sel_stage         = passes_dy_mumu_selection_stages(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_all), [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        dy_mumu_channel_sel_stage_TESUp   = passes_dy_mumu_selection_stages(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_all), [tau for tau in sel_taus if (tau[0]*tau[1][1]).pt() > TAUS_PT_CUT], proc_met)
        dy_mumu_channel_sel_stage_TESDown = passes_dy_mumu_selection_stages(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_all), [tau for tau in sel_taus if (tau[0]*tau[1][2]).pt() > TAUS_PT_CUT], proc_met)
        dy_mumu_channel_sel_stage_JERUp   = passes_dy_mumu_selection_stages(passed_triggers, leps, (N_jets_JERUp_med,   N_jets_JERUp_all),   [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        dy_mumu_channel_sel_stage_JERDown = passes_dy_mumu_selection_stages(passed_triggers, leps, (N_jets_JERDown_med, N_jets_JERDown_med), [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        dy_mumu_channel_sel_stage_JESUp   = passes_dy_mumu_selection_stages(passed_triggers, leps, (N_jets_JESUp_med,   N_jets_JESUp_all),   [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        dy_mumu_channel_sel_stage_JESDown = passes_dy_mumu_selection_stages(passed_triggers, leps, (N_jets_JESDown_med, N_jets_JESDown_med), [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)

        # elmu dilepton selection
        em_channel_sel_stage = passes_elmu_selection_stages(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_all), [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        em_channel_sel_stage_TESUp   = passes_elmu_selection_stages(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_all), [tau for tau in sel_taus if (tau[0]*tau[1][1]).pt() > TAUS_PT_CUT], proc_met)
        em_channel_sel_stage_TESDown = passes_elmu_selection_stages(passed_triggers, leps, (N_jets_nom_med, N_jets_nom_all), [tau for tau in sel_taus if (tau[0]*tau[1][2]).pt() > TAUS_PT_CUT], proc_met)
        em_channel_sel_stage_JERUp   = passes_elmu_selection_stages(passed_triggers, leps, (N_jets_JERUp_med,   N_jets_JERUp_all),   [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        em_channel_sel_stage_JERDown = passes_elmu_selection_stages(passed_triggers, leps, (N_jets_JERDown_med, N_jets_JERDown_all), [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        em_channel_sel_stage_JESUp   = passes_elmu_selection_stages(passed_triggers, leps, (N_jets_JESUp_med,   N_jets_JESUp_all),   [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)
        em_channel_sel_stage_JESDown = passes_elmu_selection_stages(passed_triggers, leps, (N_jets_JESDown_med, N_jets_JESDown_all), [tau for tau in sel_taus if (tau[0]*tau[1][0]).pt() > TAUS_PT_CUT], proc_met)

        '''
        if tt_channel_sel_stage > 0:
            tt_channel_stage = 100 + tt_channel_sel_stage
        else:
            tt_channel_stage = tt_channel_presel_stage

        if tt_channel_sel_stage_TESUp > 0:
            tt_channel_sel_stage_TESUp += 100
        if tt_channel_sel_stage_TESDown > 0:
            tt_channel_sel_stage_TESDown += 100

        if tt_channel_sel_stage_JESUp > 0:
            tt_channel_sel_stage_JESUp += 100
        if tt_channel_sel_stage_JESDown > 0:
            tt_channel_sel_stage_JESDown += 100

        if tt_channel_sel_stage_JERUp > 0:
            tt_channel_sel_stage_JERUp += 100
        if tt_channel_sel_stage_JERDown > 0:
            tt_channel_sel_stage_JERDown += 100
        '''

        #passes = tt_channel_stage < 1 and tt_channel_sel_stage_TESUp < 1 and tt_channel_sel_stage_TESDown < 1 and tt_channel_sel_stage_JESUp < 1 and tt_channel_sel_stage_JESDown < 1 and tt_channel_sel_stage_JERUp < 1 and tt_channel_sel_stage_JERDown < 1 and dy_channel_sel_stage < 1
        notpasses_tt_leptau = tt_channel_sel_stage < 1 and tt_channel_sel_stage_TESUp < 1 and tt_channel_sel_stage_TESDown < 1 and tt_channel_sel_stage_JESUp < 1 and tt_channel_sel_stage_JESDown < 1 and tt_channel_sel_stage_JERUp < 1 and tt_channel_sel_stage_JERDown < 1
        # and em_channel_sel_stage < 1 and tt_channel_stage_alliso < 1
        #notpasses_tt_elmu   = em_channel_sel_stage < 1 and em_channel_sel_stage_TESUp < 1 and em_channel_sel_stage_TESDown < 1 and em_channel_sel_stage_JESUp < 1 and em_channel_sel_stage_JESDown < 1 and em_channel_sel_stage_JERUp < 1 and em_channel_sel_stage_JERDown < 1
        notpasses_em        = em_channel_sel_stage < 1 and em_channel_sel_stage_TESUp < 1 and em_channel_sel_stage_TESDown < 1 and em_channel_sel_stage_JERUp < 1 and em_channel_sel_stage_JERDown < 1 and em_channel_sel_stage_JESUp < 1 and em_channel_sel_stage_JESDown < 1
        notpasses_alliso = tt_channel_stage_alliso < 1

        notpasses_dy_tautau = dy_channel_sel_stage < 1 and dy_channel_sel_stage_TESUp < 1 and dy_channel_sel_stage_TESDown < 1 and dy_channel_sel_stage_JESUp < 1 and dy_channel_sel_stage_JESDown < 1 and dy_channel_sel_stage_JERUp < 1 and dy_channel_sel_stage_JERDown < 1
        notpasses_dy_mumu   = dy_mumu_channel_sel_stage < 1 and dy_mumu_channel_sel_stage_TESUp < 1 and dy_mumu_channel_sel_stage_TESDown < 1 and dy_mumu_channel_sel_stage_JESUp < 1 and dy_mumu_channel_sel_stage_JESDown < 1 and dy_mumu_channel_sel_stage_JERUp < 1 and dy_mumu_channel_sel_stage_JERDown < 1

        notpasses_presel = tt_channel_presel_stage < 1
        notpasses_wjets  = selection_wjets_control < 1

        if notpasses_tt_leptau and notpasses_em and notpasses_dy_tautau and notpasses_dy_mumu and notpasses_wjets and notpasses_alliso and notpasses_presel:
            continue

        # SAVE SELECTION, objects and weights

        selection_stage_presel[0] = tt_channel_presel_stage

        selection_stage[0] = tt_channel_sel_stage
        selection_stage_TESUp[0]   = tt_channel_sel_stage_TESUp
        selection_stage_TESDown[0] = tt_channel_sel_stage_TESDown

        selection_stage_JESUp[0]   = tt_channel_sel_stage_JESUp
        selection_stage_JESDown[0] = tt_channel_sel_stage_JESDown
        selection_stage_JERUp[0]   = tt_channel_sel_stage_JERUp
        selection_stage_JERDown[0] = tt_channel_sel_stage_JERDown

        selection_stage_tt_alliso[0] = tt_channel_stage_alliso

        selection_stage_wjets[0] = selection_wjets_control

        selection_stage_dy[0] = dy_channel_sel_stage
        selection_stage_dy_TESUp  [0] = dy_channel_sel_stage_TESUp   
        selection_stage_dy_TESDown[0] = dy_channel_sel_stage_TESDown 
        selection_stage_dy_JERUp  [0] = dy_channel_sel_stage_JERUp   
        selection_stage_dy_JERDown[0] = dy_channel_sel_stage_JERDown 
        selection_stage_dy_JESUp  [0] = dy_channel_sel_stage_JESUp   
        selection_stage_dy_JESDown[0] = dy_channel_sel_stage_JESDown 

        selection_stage_dy_mumu[0]         = dy_mumu_channel_sel_stage
        selection_stage_dy_mumu_TESUp  [0] = dy_mumu_channel_sel_stage_TESUp   
        selection_stage_dy_mumu_TESDown[0] = dy_mumu_channel_sel_stage_TESDown 
        selection_stage_dy_mumu_JERUp  [0] = dy_mumu_channel_sel_stage_JERUp   
        selection_stage_dy_mumu_JERDown[0] = dy_mumu_channel_sel_stage_JERDown 
        selection_stage_dy_mumu_JESUp  [0] = dy_mumu_channel_sel_stage_JESUp   
        selection_stage_dy_mumu_JESDown[0] = dy_mumu_channel_sel_stage_JESDown 

        selection_stage_em[0] = em_channel_sel_stage
        selection_stage_em_TESUp  [0] = em_channel_sel_stage_TESUp  
        selection_stage_em_TESDown[0] = em_channel_sel_stage_TESDown
        selection_stage_em_JERUp  [0] = em_channel_sel_stage_JERUp  
        selection_stage_em_JERDown[0] = em_channel_sel_stage_JERDown
        selection_stage_em_JESUp  [0] = em_channel_sel_stage_JESUp  
        selection_stage_em_JESDown[0] = em_channel_sel_stage_JESDown

        runNumber[0]   = ev.runNumber
        indexevents[0] = ev.indexevents
        lumi[0]        = ev.lumi

        # calc lj var
        requires_lj = True # pass_elmu_close or pass_old_lep_presel or pass_old_mu_sel or pass_old_el_sel
        jets_input_has = -11
        jets_found_has = -11
        lj_var = -11.
        if requires_lj:
            # all jets, without regard to tau in the event go into the calculation
            # (taumatched jets go too)
            with_all_permutation_masses = True
            # order: not-b-taged, b-taged
            # only medium b-tags p8 p9
            b_cand_jets     = jets.medium + jets.taumatched[0]
            light_cand_jets = jets.rest + jets.loose + jets.taumatched[1]
            lj_var, w_mass, t_mass, lj_gens, all_masses = calc_lj_var(ev, light_cand_jets, b_cand_jets, with_all_permutation_masses, isMC)
            # medium and loose b-tags
            # but tau-matches are done to Medium b!
            #lj_var, w_mass, t_mass, lj_gens, all_masses = calc_lj_var(ev, jets.rest + jets.taumatched[1], jets.medium + jets.loose + jets.taumatched[0], with_all_permutation_masses, isMC)
            n_bjets_used_in_lj = len(b_cand_jets)
            lj_cut = 60.
            large_lj = lj_var > lj_cut

            # for tt_lj MC truth find the jets in these groups
            if (abs(ev.gen_t_w_decay_id) == 13 or abs(ev.gen_t_w_decay_id) == 11) and abs(ev.gen_tb_w_decay_id) == 1:
                # in tb the jets' genmatches are -6 (b) and -5 (w)
                propper_b = -6
                propper_w = -5
            elif abs(ev.gen_t_w_decay_id) == 1 and (abs(ev.gen_tb_w_decay_id) == 13 or abs(ev.gen_tb_w_decay_id) == 11):
                propper_b = 6
                propper_w = 5
            else:
                propper_b = 0
                propper_w = 0

            jets_input_has = 0
            #3*(lj_b_gen == propper_b) + (lj_w1_gen == propper_w) + (lj_w2_gen == propper_w)
            #closest_pair_gens = (ev.jet_matching_gen[light_jets[i][6]], ev.jet_matching_gen[light_jets[u][6]])
            if propper_b:
                # what was in the input
                for b_cand in b_cand_jets:
                    if b_cand[6] == propper_b:
                        jets_input_has += 3
                        break

                found_light = False
                for cand in light_cand_jets:
                    if cand[6] == propper_w:
                        jets_input_has += 1
                        if found_light: break
                        found_light = True

                # what was found in lj minimization
                (found_w1, found_w2), found_b = lj_gens
                if found_b  == propper_b:
                    jets_found_has += 3
                if found_w1 == propper_w:
                    jets_found_has += 1
                if found_w2 == propper_w:
                    jets_found_has += 1

        event_jets_lj_var[0] = lj_var
        event_jets_input_has[0] = jets_input_has
        event_jets_found_has[0] = jets_found_has

        event_jets_n_jets[0]  = N_jets_nom_all
        event_jets_n_bjets[0] = N_jets_nom_med
        event_jets_n_jets_lepmatched[0] = len(jets_nom.lepmatched)

        event_jets_alliso_n_jets[0]  = N_jets_nom_med_alliso
        event_jets_alliso_n_bjets[0] = N_jets_nom_all_alliso

        # objects
        #selection_objects = leps, jets, taus.medium
        #for chan_i, (chan, apply_bSF, sel_leps, sel_jets, sel_taus) in enumerate((ch for ch in passed_channels if ch[0] in selected_channels)):

        '''
        leptons = ROOT.LorentzVectorS()
        ttree_out.Branch("leptons", leptons)

        taus = ROOT.LorentzVectorS()
        ttree_out.Branch("taus", taus)
        taus_l = ROOT.LorentzVectorS()
        ttree_out.Branch("taus_l", taus_l)

        jets_b = ROOT.LorentzVectorS()
        ttree_out.Branch("jets_b", jets_b)
        jets_r = ROOT.LorentzVectorS()
        ttree_out.Branch("jets_r", jets_r)
        jets_l = ROOT.LorentzVectorS()
        ttree_out.Branch("jets_l", jets_l)
        jets_t = ROOT.LorentzVectorS()
        ttree_out.Branch("jets_t", jets_t)
        '''

        lep_p4, lep_relIso, lep_matching_gen, lep_matching_gen_dR, lep_id = sel_leps
        for i, p4 in enumerate(lep_p4):
            event_leptons    .push_back(p4)
            event_leptons_ids.push_back(lep_id[i])
            event_leptons_dxy     .push_back(ev.lep_dxy[i])
            if isMC:
                event_leptons_genmatch.push_back(lep_matching_gen[i])

        # alliso leps
        #lep_alliso_p4, lep_relIso, lep_matching_gen, lep_matching_gen_dR, lep_id = leps_all
        #leps_all = (ev.lep_alliso_p4, ev.lep_alliso_relIso, ev.lep_alliso_matching_gen, ev.lep_alliso_matching_gen_dR, ev.lep_alliso_id)
        for i, relIso in enumerate(ev.lep_alliso_relIso):
            event_leptons_alliso_reliso .push_back(relIso)
            event_leptons_alliso_p4     .push_back(ev.lep_alliso_p4[i])
            event_leptons_alliso_pdgId  .push_back(ev.lep_alliso_id[i])

        for i, tau_cand in enumerate(taus_candidates_alliso):
            #taus_candidates_alliso.append((p4, TES_factor, tau_pdgID, i, tau_ID, jetmatched))
            event_taus_alliso_p4    .push_back(tau_cand[0]*tau_cand[1])
            event_taus_alliso_pdgId .push_back(tau_cand[2])
            event_taus_alliso_IDlev .push_back(tau_cand[4])

        # mass of the dilepton system, lep-tau or lep-lep
        if len(lep_p4)>1:
            event_dilep_mass[0] = (lep_p4[0] + lep_p4[1]).mass()
        elif len(lep_p4)>0 and len(sel_taus)>0:
            # nominal TES SF
            event_dilep_mass[0] = (lep_p4[0] + sel_taus[0][0] * sel_taus[0][1][0]).mass()
        else:
            event_dilep_mass[0] = -111.

        # calc the top mass from the met+l+b
        if REQUIRE_MLB:
            lep = lep_p4[0]
            met_e, met_z = solve_met_e_z(proc_met.pt(), lep.pt(), lep.z(), lep.e())
            full_met     = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(proc_met.x(), proc_met.y(), met_z, met_e)
            met_lep = full_met + lep
            # loop over b jets and save TODO: loose b-jets?
            #b_cand_jets     = jets.medium + jets.taumatched[0] + jets.loose
            closest_distance = 1000.
            for j, mult, _, _, _, _, jet_index in jets.medium + jets.taumatched[0]:
                m_t = (met_lep + j*mult[0]).mass()
                event_top_masses_medium .push_back(m_t)
                m_t_distance = abs(m_t - m_t_quark)
                if m_t_distance < closest_distance:
                    closest_distance = m_t_distance

            for j, mult, _, _, _, _, jet_index in jets.loose:
                m_t = (met_lep + j*mult[0]).mass()
                event_top_masses_loose  .push_back(m_t)
                m_t_distance = abs(m_t - m_t_quark)
                if m_t_distance < closest_distance:
                    closest_distance = m_t_distance

            event_top_masses_closest[0]  = closest_distance

        # if sel_jets.medium:
        #     bMjet0_pt = sel_jets.medium[0][0].pt() * sel_jets.medium[0][1]
        #  n_rest_jets   = len(sel_jets.rest) #  + len(sel_jets.taumatched[1])
        #  #n_medium_jets = len(sel_jets.medium)
        #  n_medium_jets = len(sel_jets.medium) # + len(sel_jets.taumatched[0])
        #  n_loose_jets  = len(sel_jets.loose)
        # jets_to_prop_met = sel_jets.medium + sel_jets.loose + sel_jets.rest + sel_jets.taumatched[0] + sel_jets.taumatched[1] + sel_jets.lepmatched if ALL_JETS else all_sel_jets
        # jets_nom.taumatched[1].append((p4, en_factor, jet_weight_bSF_nom, jet_b_discr, HF, PF, jet_index))
        # en_factors  = (en_factor, jet_factor_JERUp, jet_factor_JERDown, jet_factor_JESUp, jet_factor_JESDown)
        for jet in sel_jets.medium:
            corrected_jet_p4 = jet[0] * jet[1][0]
            event_jets_b.push_back(corrected_jet_p4)

            en_nom, JER_up, JER_down, JES_up, JES_down = jet[1]
            event_jets_JERUp  .push_back(JER_up)
            event_jets_JERDown.push_back(JER_down)
            event_jets_JESUp  .push_back(JES_up)
            event_jets_JESDown.push_back(JES_down)

            jet_bSFweight_nom, jet_bdiscr, jet_index = jet[2][0], jet[3], jet[6]
            event_jets_b_bdiscr   .push_back(jet_bdiscr)
            if isMC:
                event_jets_b_genmatch .push_back(ev.jet_matching_gen[jet_index])
                event_jets_b_bSFweight.push_back(jet_bSFweight_nom)

        for jet in sel_jets.rest + sel_jets.loose:
            corrected_jet_p4 = jet[0] * jet[1][0]
            event_jets_r.push_back(corrected_jet_p4)

            en_nom, JER_up, JER_down, JES_up, JES_down = jet[1]
            event_jets_JERUp  .push_back(JER_up)
            event_jets_JERDown.push_back(JER_down)
            event_jets_JESUp  .push_back(JES_up)
            event_jets_JESDown.push_back(JES_down)

            jet_bSFweight_nom, jet_bdiscr, jet_index = jet[2][0], jet[3], jet[6]
            event_jets_r_bdiscr   .push_back(jet_bdiscr)
            if isMC:
                event_jets_r_genmatch .push_back(ev.jet_matching_gen[jet_index])
                event_jets_r_bSFweight.push_back(jet_bSFweight_nom)

        for jet in sel_jets.taumatched[0] + sel_jets.taumatched[1]:
            corrected_jet_p4 = jet[0] * jet[1][0]
            event_jets_t.push_back(corrected_jet_p4)

            en_nom, JER_up, JER_down, JES_up, JES_down = jet[1]
            event_jets_JERUp  .push_back(JER_up)
            event_jets_JERDown.push_back(JER_down)
            event_jets_JESUp  .push_back(JES_up)
            event_jets_JESDown.push_back(JES_down)

            jet_bSFweight_nom, jet_bdiscr, jet_index = jet[2][0], jet[3], jet[6]
            event_jets_t_bdiscr   .push_back(jet_bdiscr)
            if isMC:
                event_jets_t_genmatch .push_back(ev.jet_matching_gen[jet_index])
                event_jets_t_bSFweight.push_back(jet_bSFweight_nom)

        for jet in sel_jets.lepmatched:
            corrected_jet_p4 = jet[0] * jet[1][0]
            event_jets_l.push_back(corrected_jet_p4)

            en_nom, JER_up, JER_down, JES_up, JES_down = jet[1]
            event_jets_JERUp  .push_back(JER_up)
            event_jets_JERDown.push_back(JER_down)
            event_jets_JESUp  .push_back(JES_up)
            event_jets_JESDown.push_back(JES_down)

            jet_bSFweight_nom, jet_bdiscr, jet_index = jet[2][0], jet[3], jet[6]
            event_jets_l_bdiscr   .push_back(jet_bdiscr)
            if isMC:
                event_jets_l_genmatch .push_back(ev.jet_matching_gen[jet_index])
                event_jets_l_bSFweight.push_back(jet_bSFweight_nom)


        # taus candidates
        for tau in sel_tau_candidates:
            #tau_index = tau[3]
            event_candidate_taus    .push_back(tau[0] * tau[1])
            event_candidate_taus_ids.push_back(tau[2])

        # taus
        # taus_nom.medium.append((p4, TES_factor, tau_pdgID, i, jetmatched))
        for tau in sel_taus:
            tau_index = tau[3]
            TES_nom, TES_up, TES_down = tau[1]
            event_taus          .push_back(tau[0] * TES_nom)
            event_taus_TES_up   .push_back(TES_up / TES_nom)
            event_taus_TES_down .push_back(TES_down / TES_nom)
            event_taus_ids.push_back(tau[2])
            event_taus_IDlev.push_back(tau[4])
            if isMC:
                event_taus_genmatch.push_back(ev.tau_matching_gen[tau_index])

            # save 3ch info if possible
            tau_pat_SV_sign = -11.
            tau_SV_sign = -11.
            tau_SV_leng = -11.
            dalitz_m1 = -11.
            dalitz_m2 = -11.
            track_energy = -11.
            tau_jet_bdiscr = -11.

            tau_refit_index = ev.tau_refited_index[tau_index] # tau id in the original vector
            # require refit and dR quality of refit
            refitted = tau_refit_index > -1 and ev.tau_SV_fit_track_OS_matched_track_dR[tau_refit_index] + ev.tau_SV_fit_track_SS1_matched_track_dR[tau_refit_index] + ev.tau_SV_fit_track_SS2_matched_track_dR[tau_refit_index] < 0.002
            if refitted:
              try:
                tau_pat_SV_sign = ev.tau_flightLengthSignificance[tau_index]
                tau_SV_sign     = ev.tau_SV_geom_flightLenSign [tau_refit_index]
                tau_SV_leng     = ev.tau_SV_geom_flightLen     [tau_refit_index]
                dalitz_m1 = (ev.tau_SV_fit_track_OS_p4[tau_refit_index] + ev.tau_SV_fit_track_SS1_p4[tau_refit_index]).mass()
                dalitz_m2 = (ev.tau_SV_fit_track_OS_p4[tau_refit_index] + ev.tau_SV_fit_track_SS2_p4[tau_refit_index]).mass()
                track_sum = ev.tau_SV_fit_track_OS_p4[tau_refit_index] + ev.tau_SV_fit_track_SS1_p4[tau_refit_index] + ev.tau_SV_fit_track_SS2_p4[tau_refit_index]
                track_energy = track_sum.energy()
              except IndexError:
                  logging.error("IndexError  : %d, (%d, %d, %d, %d)" % (tau_refit_index, ev.tau_SV_fit_track_OS_p4.size(), ev.tau_SV_fit_track_SS1_p4.size(), ev.tau_SV_fit_track_SS2_p4.size(), ev.tau_p4.size()))
                  logging.error("IndexError2 : %d, %d" % (ev.indexevents, iev))

            event_taus_pat_sv_sign.push_back(tau_pat_SV_sign)
            event_taus_sv_sign.push_back(tau_SV_sign)
            event_taus_sv_leng.push_back(tau_SV_leng)
            event_taus_sv_dalitz_m1.push_back(dalitz_m1)
            event_taus_sv_dalitz_m2.push_back(dalitz_m2)
            event_taus_track_energy.push_back(track_energy)

            # dR matched jet
            tau_jet_index   = ev.tau_dR_matched_jet[tau_index]
            if tau_jet_index > -1:
                tau_jet_bdiscr = ev.jet_b_discr[tau_jet_index]
                # can recalculate dR if needed

            event_taus_jet_bdiscr.push_back(tau_jet_bdiscr)

        # also if there are taus and it's tt_lj than try to specify the origin of the fake
        # if it is required by the procs
        if isTT and sel_taus and gen_proc_id[0] == genproc_tt_lj:

            '''
            # the 0 tau is used
            # sadly root...
            # I have to convert ALL lorentzvectors to TLorentzVectors to have convenient dR
            tau_p4 = TLorentzVector(sel_taus[0][0].X(), sel_taus[0][0].Y(), sel_taus[0][0].Z(), sel_taus[0][0].T())
            # all gen products are in
            # gen_t_w1_final_p4s gen_t_w2_final_p4s gen_tb_w1_final_p4s gen_tb_w2_final_p4s
            # gen_t_b_final_p4s gen_tb_b_final_p4s
            # but neutrinos are also there the ID to check:
            # gen_t_w1_final_pdgIds etc
            gen_match_cut = 0.4
            matched_w = False
            for w_prod, w_prod_ID in zip(ev.gen_t_w1_final_p4s, ev.gen_t_w1_final_pdgIds) + zip(ev.gen_t_w2_final_p4s, ev.gen_t_w2_final_pdgIds) + \
                                     zip(ev.gen_tb_w1_final_p4s, ev.gen_tb_w1_final_pdgIds) + zip(ev.gen_tb_w2_final_p4s, ev.gen_tb_w2_final_pdgIds):
                if abs(w_prod_ID) in (12, 14, 16):
                    continue
                w_prod_p4 = TLorentzVector(w_prod.X(), w_prod.Y(), w_prod.Z(), w_prod.T())
                if w_prod_p4.DeltaR(tau_p4) < gen_match_cut:
                    matched_w = True
                    break

            matched_b = False
            for b_prod, b_prod_ID in zip(ev.gen_t_b_final_p4s, ev.gen_t_b_final_pdgIds) + zip(ev.gen_tb_b_final_p4s, ev.gen_tb_b_final_pdgIds):
                if abs(b_prod_ID) in (12, 14, 16):
                    continue
                b_prod_p4 = TLorentzVector(b_prod.X(), b_prod.Y(), b_prod.Z(), b_prod.T())
                if b_prod_p4.DeltaR(tau_p4) < gen_match_cut:
                    matched_b = True
                    break
            '''
            # the same using genmatch info
            tau_index = sel_taus[0][3]
            gen_dR = ev.tau_matching_gen_dR[tau_index]
            gen_id = ev.tau_matching_gen[tau_index]

            matched_w = gen_dR < 0.3 and abs(gen_id) == 5
            matched_b = gen_dR < 0.3 and abs(gen_id) == 6

            if matched_w and matched_b:
                gen_proc_id[0] = genproc_tt_ljo
            elif matched_w:
                gen_proc_id[0] = genproc_tt_ljw
            elif matched_b:
                gen_proc_id[0] = genproc_tt_ljb

        #weight_bSF = weight_bSF_to_apply = 1.
        # consider tau-matched jets too as in:
        #old_jet_sel = (len(jets.old.medium) + len(jets.old.taumatched[0])) > 0 and (len(jets.old.taumatched[0]) + len(jets.old.taumatched[1]) + len(jets.old.medium) + len(jets.old.loose) + len(jets.old.rest)) > 2
        # -- TODO: by the way this is a weird part for b-tagging -- are there studies of b-taging for true taus?
        #for _, _, jet_weight, _, _, _, _ in sel_jets.medium + sel_jets.loose + sel_jets.rest + sel_jets.taumatched[0] + sel_jets.taumatched[1]:
        #    weight_bSF *= jet_weight
        selection_requires_b = True
        if selection_requires_b:
            weight_bSF_to_apply = weight_bSF

        #if roccor_factor != 1.:
        #    lep_p4[0] *= roccor_factor
        #    #lep_p4[0].SetPt(lep_p4.pt() * roccor_factor) # works only for 0
        #    proc_met -= lep_p4[0] * (roccor_factor - 1.)
        #    #met_x -=
        #    #met_x -=  * (en_factor - 1.)

        #control_counters.Fill(50 + 20*sys_i + 2*chan_i + 1)

        #record_weight = sys_weight if chan not in ('sel_mu_min', 'sel_mu_min_ss', 'sel_mu_min_medtau') else sys_weight_min

        # nominal systematics
        # weight init includes basic MC weights (amcatnlo -1, Z, recoil) and averaged per epoch PU and LEP
        weight_init, weight_top_pt_NOM, weight_th = syst_weights_nominal
        nom_sys_weight            = weight_init * weight_th * weight_top_pt_NOM
        nom_sys_weight_without_PU = weight_init * weight_th * weight_top_pt_NOM / weight_pu # for PU tests

        event_weight[0] = nom_sys_weight * weight_bSF_to_apply
        event_weight_init[0]  = weight_init

        event_weight_PU[0]    = weight_pu
        event_weight_PU_el[0] = weight_pu_el
        event_weight_PU_mu[0] = weight_pu_mu
        event_weight_PU_bcdef[0] = weight_pu_bcdef
        event_weight_PU_gh[0]    = weight_pu_gh
        event_weight_PU_per_epoch[0]    = ratio_bcdef*weight_pu_bcdef + ratio_gh*weight_pu_gh

        event_weight_th[0]    = weight_th
        event_weight_bSF[0]     = weight_bSF
        event_weight_bSFUp[0]   = weight_bSFUp
        event_weight_bSFDown[0] = weight_bSFDown
        event_weight_bSF_JERUp[0]   = weight_bSF_JERUp
        event_weight_bSF_JERDown[0] = weight_bSF_JERDown
        event_weight_bSF_JESUp[0]   = weight_bSF_JESUp
        event_weight_bSF_JESDown[0] = weight_bSF_JESDown

        # systematic variation of the event weight
        if isMC:
            event_weight_toppt[0]   = weight_top_pt if isTT else 1.
            event_weight_LEPall[0]  = weight_lepall_pu

            # leptons and pu with systematics
            event_weight_LEP_PU[0]     = weight_lep_pu if weight_lep_pu > 0.0001 else 0.

            max_weight_lep_pu = 2*weight_lep_pu
            event_weight_LEPmuIDUp[0]    = (weight_lep_pu_muIDUp    / weight_lep_pu if weight_lep_pu_muIDUp    < max_weight_lep_pu else max_weight_lep_pu) if weight_lep_pu > 0.0001 and weight_lep_pu_muIDUp    > 0.0001 else 0.
            event_weight_LEPmuIDDown[0]  = (weight_lep_pu_muIDDown  / weight_lep_pu if weight_lep_pu_muIDDown  < max_weight_lep_pu else max_weight_lep_pu) if weight_lep_pu > 0.0001 and weight_lep_pu_muIDDown  > 0.0001 else 0.
            event_weight_LEPmuTRGUp[0]   = (weight_lep_pu_muTRGUp   / weight_lep_pu if weight_lep_pu_muTRGUp   < max_weight_lep_pu else max_weight_lep_pu) if weight_lep_pu > 0.0001 and weight_lep_pu_muTRGUp   > 0.0001 else 0.
            event_weight_LEPmuTRGDown[0] = (weight_lep_pu_muTRGDown / weight_lep_pu if weight_lep_pu_muTRGDown < max_weight_lep_pu else max_weight_lep_pu) if weight_lep_pu > 0.0001 and weight_lep_pu_muTRGDown > 0.0001 else 0.
            event_weight_LEPelIDUp[0]    = (weight_lep_pu_elIDUp    / weight_lep_pu if weight_lep_pu_elIDUp    < max_weight_lep_pu else max_weight_lep_pu) if weight_lep_pu > 0.0001 and weight_lep_pu_elIDUp    > 0.0001 else 0.
            event_weight_LEPelIDDown[0]  = (weight_lep_pu_elIDDown  / weight_lep_pu if weight_lep_pu_elIDDown  < max_weight_lep_pu else max_weight_lep_pu) if weight_lep_pu > 0.0001 and weight_lep_pu_elIDDown  > 0.0001 else 0.
            event_weight_LEPelTRGUp[0]   = (weight_lep_pu_elTRGUp   / weight_lep_pu if weight_lep_pu_elTRGUp   < max_weight_lep_pu else max_weight_lep_pu) if weight_lep_pu > 0.0001 and weight_lep_pu_elTRGUp   > 0.0001 else 0.
            event_weight_LEPelTRGDown[0] = (weight_lep_pu_elTRGDown / weight_lep_pu if weight_lep_pu_elTRGDown < max_weight_lep_pu else max_weight_lep_pu) if weight_lep_pu > 0.0001 and weight_lep_pu_elTRGDown > 0.0001 else 0.

            event_weight_PUUp[0]    = (weight_lep_pu_PUUp / weight_lep_pu) if weight_lep_pu > 0.0001 else 0.
            event_weight_PUDown[0]  = (weight_lep_pu_PUDown / weight_lep_pu) if weight_lep_pu > 0.0001 else 0.

            # control and reconstruction of full weight

            event_weight_LEPmuBID[0]    = weight_lep_B_MU_id  if weight_lep_B_MU_id  > 0.0001 else 0.
            event_weight_LEPmuBTRG[0]   = weight_lep_B_MU_trg if weight_lep_B_MU_trg > 0.0001 else 0.
            event_weight_LEPmuHID[0]    = weight_lep_H_MU_id  if weight_lep_H_MU_id  > 0.0001 else 0.
            event_weight_LEPmuHTRG[0]   = weight_lep_H_MU_trg if weight_lep_H_MU_trg > 0.0001 else 0.
            event_weight_LEPelID[0]     = weight_lep_EL_id  if weight_lep_EL_id  > 0.0001 else 0.
            event_weight_LEPelTRG[0]    = weight_lep_EL_trg if weight_lep_EL_trg > 0.0001 else 0.

            #event_weight_LEPmu0ID[0]     = weight_lepMU0_id  if weight_lepMU0_id  > 0.0001 else 0.
            #event_weight_LEPmu0TRG[0]    = weight_lepMU0_trg if weight_lepMU0_trg > 0.0001 else 0.
            #event_weight_LEPmu0TRGUp[0]  = weight_lepMU0_trg_Up if weight_lepMU0_trg_Up > 0.0001 else 0.
            #event_weight_LEPmuTRGctrlUp[0]   = weight_lepMU_trg_Up   if weight_lepMU_trg_Up  > 0.0001 else 0.
            #event_weight_LEPmuTRGctrlDown[0] = weight_lepMU_trg_Down if weight_lepMU_trg_Down  > 0.0001 else 0.

            event_weight_PUUp_el[0]    = weight_pu_el_up
            event_weight_PUDown_el[0]  = weight_pu_el_dn
            event_weight_PUUp_mu[0]    = weight_pu_mu_up
            event_weight_PUDown_mu[0]  = weight_pu_mu_dn

            event_weight_PUUp_bcdef[0]    = weight_pu_bcdef_up
            event_weight_PUDown_bcdef[0]  = weight_pu_bcdef_dn

            event_weight_PUUp_gh[0]    = weight_pu_gh_up
            event_weight_PUDown_gh[0]  = weight_pu_gh_dn

        if with_Frag_sys and isTT:
            event_weight_PetersonUp[0] = weights_gen_weight_Peterson / weights_gen_weight_centralFrag
            event_weight_FragUp  [0]   = weights_gen_weight_Frag[0]  / weights_gen_weight_centralFrag
            event_weight_FragDown[0]   = weights_gen_weight_Frag[1]  / weights_gen_weight_centralFrag
            event_weight_SemilepBRUp  [0] = weights_gen_weight_semilepbr[0] / weights_gen_weight_centralFrag
            event_weight_SemilepBRDown[0] = weights_gen_weight_semilepbr[1] / weights_gen_weight_centralFrag
        else:
            event_weight_PetersonUp[0]    = 1.0
            event_weight_FragUp  [0]      = 1.0
            event_weight_FragDown[0]      = 1.0
            event_weight_SemilepBRUp  [0] = 1.0
            event_weight_SemilepBRDown[0] = 1.0

        if with_MEscale_sys and isTT:
            #event_weight_me_f_ = weights_gen_weight_nom   #= ev.gen_weights_renorm_fact[MUf_nom_MUr_nom] if ev.gen_weights_renorm_fact[MUf_nom_MUr_nom] > 0. else 0.00001
            event_weight_me_f_rUp[0] = weights_gen_weight_f_rUp #= ev.gen_weights_renorm_fact[MUf_nom_MUr_up]
            event_weight_me_f_rDn[0] = weights_gen_weight_f_rDn #= ev.gen_weights_renorm_fact[MUf_nom_MUr_down]
            event_weight_me_fUp_r[0] = weights_gen_weight_fUp_r #= ev.gen_weights_renorm_fact[MUf_up_MUr_nom]
            event_weight_me_fDn_r[0] = weights_gen_weight_fDn_r #= ev.gen_weights_renorm_fact[MUf_down_MUr_nom]
            event_weight_me_frUp [0] = weights_gen_weight_frUp  #= ev.gen_weights_renorm_fact[MUf_up_MUr_up]
            event_weight_me_frDn [0] = weights_gen_weight_frDn  #= ev.gen_weights_renorm_fact[MUf_down_MUr_down]

        if with_PDF_sys and isTT:
            for i in range(1,57):
                pdf_w = ev.gen_weights_pdf_hessians[i] / weights_gen_weight_norm
                # 1 event in muon selection has a number of PDF sets far above nominal value, weighted at 18 and more
                # row 92098 in
                # /gstore/t3cms/store/user/otoldaie/v23/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v23_MC2016_Summer16_TTJets_powheg/180317_001114/0000/MC2016_Summer16_TTJets_powheg_91.root
                # for PDF n2 (one of the problematic ones) there are 2 more events with pdf weight more than 2
                # the rest (far larger majority) are below
                if pdf_w > 2.:
                    pdf_w = 1.

                event_weight_pdf.push_back(pdf_w)  #= ev.gen_weights_renorm_fact[MUf_down_MUr_down]

                #pdf_sys_name = pdf_sys_name_up(i) #pdf_sys_names_Up[i-1]
                #control_hs[pdf_sys_name].Fill(pdf_w)
                #out_hs[(chan, record_proc, pdf_sys_name)]['Mt_lep_met_f'] .Fill(Mt_lep_met, record_weight * pdf_w)
                #pdf_sys_name = pdf_sys_name_down(i) #pdf_sys_names_Down[i-1] # down is nominal
                #out_hs[(chan, record_proc, pdf_sys_name)]['Mt_lep_met_f'] .Fill(Mt_lep_met, record_weight)

            event_weight_AlphaS_up [0] = weights_gen_weight_alphas_up
            event_weight_AlphaS_dn [0] = weights_gen_weight_alphas_dn

        # propagation of corrections to met

        met_x_prop_taus = proc_met.Px() #met_x
        met_y_prop_taus = proc_met.Py() #met_y

        met_x_prop_jets = proc_met.Px() # met_x
        met_y_prop_jets = proc_met.Py() # met_y

        met_x_prop = proc_met.Px() # met_x
        met_y_prop = proc_met.Py() # met_y

        # visible particles + met sum control
        sum_jets_corr_X = 0.
        sum_jets_corr_Y = 0.
        sum_jets_init_X = 0.
        sum_jets_init_Y = 0.
        #
        sum_tau_corr_X = 0.
        sum_tau_corr_Y = 0.
        sum_tau_init_X = 0.
        sum_tau_init_Y = 0.

        sum_tau_corr_pt = 0.
        substitution_pt = 0.

        # PROPAGATE tau and jet systematic variations to met
        # 
        # taus = [(p4, TES_factor, tau_pdgID)]
        # nominal taus do have a factor
        if isMC and sel_taus: # and 'TauES' in sys_name:
            #try:
            sum_tau_init = sel_taus[0][0]
            sum_tau_corr = sel_taus[0][0] * sel_taus[0][1][0]
            tau_cor = sel_taus[0][0] * (1. - sel_taus[0][1][0])
            for tau in sel_taus[1:]:
                tau_cor += tau[0] * (1. - tau[1][0])
                sum_tau_init += tau[0]
                sum_tau_corr += tau[0] * tau[1][0]

            sum_tau_corr_X = sum_tau_corr.Px()
            sum_tau_corr_Y = sum_tau_corr.Py()
            sum_tau_init_X = sum_tau_init.Px()
            sum_tau_init_Y = sum_tau_init.Py()

            sum_tau_corr_pt = tau_cor.pt()

            #except IndexError:
            #    print len(sel_taus), type(sel_taus)
            #    print sel_taus
            #    print len(sel_taus[0])
            #    raise IndexError
            met_x_prop_taus += tau_cor.X()
            met_y_prop_taus += tau_cor.Y()
            met_x_prop += tau_cor.X()
            met_y_prop += tau_cor.Y()

        #### and substitute the jet->tau in met p7 p9 -> 1) v25 p2_tt_jtau, 2) v25 p2_jes_recor
        ### this works very strangely: data is shifted to high Mt?
        ### but the study of jet/tau pt shows approximatly the same values in both MC and Data
        #if sel_taus and sel_taus[0][4] > -1:
        #    # only first tau is taken
        #    tau_index = sel_taus[0][3]
        #    the_tau_p4 = sel_taus[0][0] * sel_taus[0][1][0]
        #    tau_jet_index   = sel_taus[0][4]

        #    # substitute the nominal jet
        #    jer_factor = ev.jet_jer_factor[tau_jet_index] if isMC else 1.
        #    #jes_factor = ev.jet_jes_recorrection[tau_jet_index]
        #    #jes_uncorFactor = ev.jet_uncorrected_jecFactor[tau_jet_index]
        #    en_factor = jer_factor # * jes_factor * jes_uncorFactor
        #    the_jet_p4 = ev.jet_initial_p4[tau_jet_index] * en_factor #miniaod_jets[tau_jet_index]
        #    #substitution = the_tau_p4 - the_jet_p4
        #    # + what is to remove from met
        #    # - what is to include in met
        #    substitution = the_jet_p4 - the_tau_p4
        #    substitution_pt = substitution.pt()
        #    #met_x_prop += substitution.X()
        #    #met_y_prop += substitution.Y()

        # PROPAGATE jet correcions
        all_sel_jets = sel_jets.medium + sel_jets.loose + sel_jets.rest
        all_sel_jets_taumatched = sel_jets.medium + sel_jets.loose + sel_jets.rest + sel_jets.taumatched[0] + sel_jets.taumatched[1]
        jets_to_prop_met = sel_jets.medium + sel_jets.loose + sel_jets.rest + sel_jets.taumatched[0] + sel_jets.taumatched[1] + sel_jets.lepmatched if ALL_JETS else all_sel_jets
        # propagate above
        ## propagation of corrections of all jets
        ## nominal jets have JER factor
        #if isMC and jets_to_prop_met: # and ('JES' in sys_name or 'JER' in sys_name):
        #    #try:
        #    sum_jets_init = jets_to_prop_met[0][0]
        #    sum_jets_corr = jets_to_prop_met[0][0] * jets_to_prop_met[0][1]
        #    jet_cor       = jets_to_prop_met[0][0] * (1. - jets_to_prop_met[0][1])
        #    #for jet in all_sel_jets[1:] + sel_jets.taumatched[0] + sel_jets.taumatched[1]:
        #    for jet in jets_to_prop_met[1:]:
        #        jet_cor += jet[0] * (1. - jet[1])
        #        sum_jets_init += jet[0]
        #        sum_jets_corr += jet[0] * jet[1]

        #    sum_jets_corr_X = sum_jets_corr.Px()
        #    sum_jets_corr_Y = sum_jets_corr.Py()
        #    sum_jets_init_X = sum_jets_init.Px()
        #    sum_jets_init_Y = sum_jets_init.Py()

        #    met_x_prop_jets += jet_cor.X()
        #    met_y_prop_jets += jet_cor.Y()
        #    met_x_prop += jet_cor.X()
        #    met_y_prop += jet_cor.Y()
        #    # met prop = met nom + nom - factor
        #    # met prop + factor = met nom + nom

        # control over "objects' met"
        #all_sel_objects = LorentzVector(lep_p4[0].X(), lep_p4[0].Y(), lep_p4[0].Z(), lep_p4[0].T())
        #ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(100.,0.,1.,200.)
        all_sel_objects = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(0.,0.,0.,0.)
        for lep in lep_p4:
            all_sel_objects += lep
        for tau in sel_taus:
            all_sel_objects += tau[0] * tau[1][0]
        for jet in all_sel_jets:
            all_sel_objects += jet[0] * jet[1][0]

        '''
        Mt_lep_met_init = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), ev.met_init.Px(), ev.met_init.Py())
        # at NOMINAL these two should be the same
        # only systematic variations differ
        Mt_lep_met         = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), met_x_prop, met_y_prop)
        if isMC:
            Mt_lep_met_shifted = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), met_x_prop + 7.2, met_y_prop + 1.4)
        else:
             Mt_lep_met_shifted = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), met_x_prop, met_y_prop)
        Mt_lep_met_corr    = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), ev.met_corrected.Px(), ev.met_corrected.Py())
        # test on Mt fluctuation in mu_sel
        #if Mt_lep_met < 1.:
        #    continue
        Mt_lep_met_sublep = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), proc_met_lepsub.Px(), proc_met_lepsub.Py())

        met_pt = TMath.Sqrt(met_x_prop*met_x_prop + met_y_prop*met_y_prop)
        met_pt_taus = TMath.Sqrt(met_x_prop_taus*met_x_prop_taus + met_y_prop_taus*met_y_prop_taus)
        met_pt_jets = TMath.Sqrt(met_x_prop_jets*met_x_prop_jets + met_y_prop_jets*met_y_prop_jets)
        met_pt_init = proc_met.pt() # TMath.Sqrt(met_x*met_x + met_y*met_y)

        met_cancellation_x = met_x_prop + all_sel_objects.Px()
        met_cancellation_y = met_y_prop + all_sel_objects.Py()
        met_cancelation = TMath.Sqrt(met_cancellation_x*met_cancellation_x + met_cancellation_y*met_cancellation_y)
        '''

        #event_met.SetPx(met_x_prop)
        #event_met.SetPy(met_y_prop)

        event_met.SetPx(proc_met.Px())
        event_met.SetPy(proc_met.Py())

        event_met_JERUp.SetPx(proc_met_JERUp.Px())
        event_met_JERUp.SetPy(proc_met_JERUp.Py())
        event_met_JERDown.SetPx(proc_met_JERDown.Px())
        event_met_JERDown.SetPy(proc_met_JERDown.Py())

        event_met_JESUp.SetPx(proc_met_JESUp.Px())
        event_met_JESUp.SetPy(proc_met_JESUp.Py())
        event_met_JESDown.SetPx(proc_met_JESDown.Px())
        event_met_JESDown.SetPy(proc_met_JESDown.Py())

        event_met_TESUp.SetPx(proc_met_TESUp.Px())
        event_met_TESUp.SetPy(proc_met_TESUp.Py())
        event_met_TESDown.SetPx(proc_met_TESDown.Px())
        event_met_TESDown.SetPy(proc_met_TESDown.Py())

        event_met_lep_mt         [0] = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), proc_met.Px(), proc_met.Py())
        event_met_lep_mt_JERUp   [0] = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), proc_met_JERUp.Px(),   proc_met_JERUp.Py())
        event_met_lep_mt_JERDown [0] = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), proc_met_JERDown.Px(), proc_met_JERDown.Py())
        event_met_lep_mt_JESUp   [0] = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), proc_met_JESUp.Px(),   proc_met_JESUp.Py())
        event_met_lep_mt_JESDown [0] = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), proc_met_JESDown.Px(), proc_met_JESDown.Py())
        event_met_lep_mt_TESUp   [0] = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), proc_met_TESUp.Px(),   proc_met_TESUp.Py())
        event_met_lep_mt_TESDown [0] = transverse_mass_pts(lep_p4[0].Px(), lep_p4[0].Py(), proc_met_TESDown.Px(), proc_met_TESDown.Py())

        ttree_out.Fill()


        if save_weights:
          #weight_bSF = 1.
          #weight_bSF_up, weight_bSF_down = 1., 1.
          #weight_bSF_jer_up, weight_bSF_jer_down = 1., 1.
          #weight_bSF_jes_up, weight_bSF_jes_down = 1., 1.
          control_hs['weight_z_mass_pt'] .Fill(weight_z_mass_pt)
          control_hs['weight_bSF']    .Fill(weight_bSF)
          #control_hs['weight_top_pt'] .Fill(weight_top_pt) # done above

          if pass_mu: # CHECK: add these for dileptons if needed or pass_mumu: # or pass_elmu
            # bcdef_weight_trk, bcdef_weight_id, bcdef_weight_iso, gh_weight_trk, gh_weight_id, gh_weight_iso
            #mu_sfs = lepton_muon_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt())
            #mu_trg_sf = lepton_muon_trigger_SF(ev.lep_p4[0].eta(), ev.lep_p4[0].pt())

            control_hs['weight_mu_trk_bcdef']        .Fill(mu_sfs_b[0])
            control_hs['weight_mu_trk_bcdef_vtx_gen'].Fill(mu_sfs_b[1])
            control_hs['weight_mu_trk_bcdef_vtx']    .Fill(mu_sfs_b[2])
            control_hs['weight_mu_id_bcdef']         .Fill(mu_sfs_b[3][0])
            control_hs['weight_mu_iso_bcdef']        .Fill(mu_sfs_b[4][0])
            control_hs['weight_mu_trg_bcdef'].Fill(mu_trg_sf_b)
            control_hs['weight_mu_all_bcdef'].Fill(mu_trg_sf_b * mu_sfs_b[0] * mu_sfs_b[1] * mu_sfs_b[3][0] * mu_sfs_b[4][0])

            control_hs['weight_mu_trk_gh']        .Fill(mu_sfs_h[0])
            control_hs['weight_mu_trk_gh_vtx_gen'].Fill(mu_sfs_h[1])
            control_hs['weight_mu_trk_gh_vtx']    .Fill(mu_sfs_h[2])
            control_hs['weight_mu_id_gh']         .Fill(mu_sfs_h[3][0])
            control_hs['weight_mu_iso_gh']        .Fill(mu_sfs_h[4][0])
            control_hs['weight_mu_trg_gh'] .Fill(mu_trg_sf_h)
            control_hs['weight_mu_all_gh'] .Fill(mu_trg_sf_h * mu_sfs_h[0] * mu_sfs_h[1] * mu_sfs_h[3][0] * mu_sfs_h[4][0])

            control_hs['weight_mu_bSF']     .Fill(weight_bSF)
            #control_hs['weight_mu_bSF_up']  .Fill(weight_bSF_up)
            #control_hs['weight_mu_bSF_down'].Fill(weight_bSF_down)

          elif pass_el: # or pass_elel or pass_elmu_el:
            #el_sfs = lepton_electron_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())
            #el_trg_sf = lepton_electron_trigger_SF(abs(ev.lep_p4[0].eta()), ev.lep_p4[0].pt())

            control_hs['weight_el_trk'].Fill(el_sfs_reco[0])
            control_hs['weight_el_idd'].Fill(el_sfs_id[0])
            control_hs['weight_el_trg'].Fill(el_trg_sf)
            control_hs['weight_el_all'].Fill(el_trg_sf * el_sfs_reco[0] * el_sfs_id[0])

            control_hs['weight_el_bSF']     .Fill(weight_bSF)
            #control_hs['weight_el_bSF_up']  .Fill(weight_bSF_up)
            #control_hs['weight_el_bSF_down'].Fill(weight_bSF_down)

    profile.disable()
    #profile.print_stats()
    # there is no returning string
    #profile.dump_stats()

    return control_hs, control_counters, profile





#def main(input_dir, dtag, outdir, range_min, range_max):
def main(input_filename, fout_name, outdir, channels_to_select, lumi_bcdef=19714., lumi_gh=16146.):
    '''main(input_filename, outdir, range_min, range_max, lumi_bcdef=19252.03, lumi_gh=16290.02)

    lumi defaults are from _full_ golden json for muon
    -- the bad lumis don't reduce it that much, should not affect the ratio too much
    '''

    print " OLD_MINIAOD_JETS DO_W_STITCHING ALL_JETS with_bSF"
    print OLD_MINIAOD_JETS, DO_W_STITCHING, ALL_JETS, with_bSF

    start_time = datetime.now()

    input_tfile = TFile(input_filename)
    tree = input_tfile.Get('ntupler/reduced_ttree')

    """
    if not range_max: range_max = tree.GetEntries()

    fout_name = input_filename.split('/')[-1].split('.root')[0] + "_%d-%d.root" % (range_min, range_max)

    if isfile(outdir + '/' + fout_name):
        print "output file exists: %s" % (outdir + '/' + fout_name)
        return None
    """

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
    logger.write("%s\n" % str(start_time))
    logger.write("input file = %s\n" % input_filename)
    logger.write("output dir = %s\n" % outdir)
    #f = TFile('outdir/v12.3/merged-sets/MC2016_Summer16_TTJets_powheg.root')

    logger.write("N entries = %s\n" % tree.GetEntries())

    fout = TFile(outdir + '/' + fout_name, "RECREATE")
    ttree_out = TTree( 'ttree_out', 'tree with stage2 selection' ) # INPUT NOW
    ttree_out.SetDirectory(fout)

    #logger.write("range = %d, %d\n" % (range_min, range_max))

    logger.write("lumi BCDEF GH = %f %f\n" % (lumi_bcdef, lumi_gh))
    c_hs, control_counters, perf_profile = full_loop(tree, ttree_out, input_filename, lumi_bcdef, lumi_gh, logger, channels_to_select=channels_to_select)

    perf_profile.dump_stats(logger_file.split('.log')[0] + '.cprof')

    events_counter = input_tfile.Get('ntupler/events_counter')
    weight_counter = input_tfile.Get('ntupler/weight_counter')
    systematic_weights = input_tfile.Get('ntupler/systematic_weights')
    events_counter.SetDirectory(0)
    weight_counter.SetDirectory(0)
    systematic_weights.SetDirectory(0)

    input_tfile.Close()

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

    #fout = TFile("lets_test.root", "RECREATE")
    #fout = TFile(outdir + '/' + dtag + "_%d-%d.root" % (range_min, range_max), "RECREATE")
    #fout_name = outdir + '/' + 
    fout.Write()

    fout.cd()

    ttree_out.Write()

    #events_counter.SetDirectory()
    #weight_counter.SetDirectory()
    events_counter.Write() # hopefully these go to the root of the tfile
    weight_counter.Write()
    systematic_weights.Write()
    control_counters.Write()

    fout.Write()
    fout.Close()

    work_time = str(datetime.now() - start_time)
    print "output written, time elapsed %s" % work_time
    logger.write("time elapsed %s\n" % work_time)

    ##
    #print "trying to exit without segfaults"
    ##ROOT.gSystem.Unload("libUserCodettbar-leptons-80X.so")
    ##ROOT.gSystem.Unload("RoccoR_cc")
    #ROOT.gROOT.Reset()
    #print "all commands done"

