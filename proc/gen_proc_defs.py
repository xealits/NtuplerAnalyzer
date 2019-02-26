
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

# >>>>>>>>>>>>>> gen procs defs
genproc_dy_tautau = 1
genproc_dy_other  = 0

genproc_wjets_tauh = 2
genproc_wjets_taul = 1
genproc_wjets      = 0

genproc_qcd = genproc_data = 0

genproc_stop_el    = 20
genproc_stop_mu    = 10
genproc_stop_lj    = 2
genproc_stop_elmu  = 1
genproc_stop_other = 0

#procs_el_3ch_fbw = tt_procs_el_3ch_fbw =  (['tt_eltau3ch', 'tt_eltau', 'tt_ljb', 'tt_ljw', 'tt_ljo', 'tt_lj', 'tt_taultauh', 'tt_taulj', 'tt_other'], 'tt_other')
#procs_mu_3ch_fbw = tt_procs_mu_3ch_fbw =  (['tt_mutau3ch', 'tt_mutau', 'tt_ljb', 'tt_ljw', 'tt_ljo', 'tt_lj', 'tt_taultauh', 'tt_taulj', 'tt_other'], 'tt_other')
#procs_elmu   = tt_procs_elmu   =  (['tt_elmu', 'tt_ltaul', 'tt_taueltaumu', 'tt_other'], 'tt_other')

genproc_tt_eltau3ch  = 42
genproc_tt_eltau     = 41
genproc_tt_mutau3ch  = 32
genproc_tt_mutau     = 31

genproc_tt_ljb       = 24
genproc_tt_ljw       = 23
genproc_tt_ljo       = 22
genproc_tt_lj        = 21
genproc_tt_taultauh  = 12
genproc_tt_taulj     = 11

genproc_tt_elmu       = 3
genproc_tt_ltaul      = 2
genproc_tt_taueltaumu = 1

genproc_tt_other = 0

# <<<<<<<<<<<<<<


genprocs_mutau_dy       = ('dy',    [('tautau', [genproc_dy_tautau])])
genprocs_mutau_qcd      = ('qcd',   [])
genprocs_mutau_wjets    = ('wjets', [('taul', [genproc_wjets_taul]), ('tauh', [genproc_wjets_tauh])])
genprocs_mutau_dibosons = ('dibosons', [])

genprocs_mutau_s_top = ('s_top', [('mutau', [genproc_stop_mu]),  ('lj', [genproc_stop_lj])])
genprocs_mutau_tt    = ('tt',    [('mutau', [genproc_tt_mutau, genproc_tt_mutau3ch]),
                                  ('lj', [genproc_tt_lj, genproc_tt_ljb, genproc_tt_ljw, genproc_tt_ljo]),
                                  ('taultauh', [genproc_tt_taultauh]),
                                  ('taulj', [genproc_tt_taulj])])

genprocs_eltau_s_top = ('s_top', [('eltau', [genproc_stop_el]),  ('lj', [genproc_stop_lj])])
genprocs_eltau_tt    = ('tt',    [('eltau', [genproc_tt_eltau, genproc_tt_eltau3ch]),
                                  ('lj', [genproc_tt_lj, genproc_tt_ljb, genproc_tt_ljw, genproc_tt_ljo]),
                                  ('taultauh', [genproc_tt_taultauh]),
                                  ('taulj', [genproc_tt_taulj])])


genprocs_ltau_dy       = {'el': ('dy',    [('tautau', [genproc_dy_tautau])])}
genprocs_ltau_qcd      = {'el': ('qcd',   [])}
genprocs_ltau_wjets    = {'el': ('wjets', [('taul', [genproc_wjets_taul]), ('tauh', [genproc_wjets_tauh])])}
genprocs_ltau_dibosons = {'el': ('dibosons', [])}

genprocs_ltau_s_top = {'el': ('s_top', [('eltau', [genproc_stop_el]),  ('lj', [genproc_stop_lj])])}
genprocs_ltau_tt    = {'el': ('tt',    [('eltau', [genproc_tt_eltau, genproc_tt_eltau3ch]),
                                  ('lj', [genproc_tt_lj, genproc_tt_ljb, genproc_tt_ljw, genproc_tt_ljo]),
                                  ('taultauh', [genproc_tt_taultauh]),
                                  ('taulj', [genproc_tt_taulj])])}

genprocs_ltau_s_top['mu'] = ('s_top', [('mutau', [genproc_stop_mu]),  ('lj', [genproc_stop_lj])])
genprocs_ltau_tt   ['mu'] = ('tt',    [('mutau', [genproc_tt_mutau, genproc_tt_mutau3ch]),
                                  ('lj', [genproc_tt_lj, genproc_tt_ljb, genproc_tt_ljw, genproc_tt_ljo]),
                                  ('taultauh', [genproc_tt_taultauh]),
                                  ('taulj', [genproc_tt_taulj])])

dtags_procs = {
'MC2016_Summer16_DYJetsToLL_10to50_amcatnlo'   : genprocs_ltau_dy,
'MC2016_Summer16_DYJetsToLL_50toInf_madgraph'  : genprocs_ltau_dy,
'MC2016_Summer16_QCD_HT-100-200'    : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_HT-1000-1500'  : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_HT-1500-2000'  : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_HT-200-300'    : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_HT-2000-Inf'   : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_HT-300-500'    : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_HT-500-700'    : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_HT-700-1000'   : genprocs_ltau_qcd,

'MC2016_Summer16_QCD_MuEnriched_Pt15_20toInf'  : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_1000toInf' : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_120to170'  : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_15to20'    : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_170to300'  : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_20to30'    : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_300to470'  : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_30to50'    : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_470to600'  : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_50to80'    : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_600to800'  : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_800to1000' : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_80to120'   : genprocs_ltau_qcd,

'MC2016_Summer16_QCD_EMEnriched_Pt-120to170.root' : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-170to300.root' : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-20to30.root'   : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-300toInf.root' : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-30to40.root'   : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-30to50.root'   : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-30toInf.root'  : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-40toInf.root'  : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-50to80.root'   : genprocs_ltau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-80to120.root'  : genprocs_ltau_qcd,

'MC2016_Summer16_SingleT_tW_5FS_powheg'                     : genprocs_ltau_s_top,
'MC2016_Summer16_SingleTbar_tW_5FS_powheg'                  : genprocs_ltau_s_top,
'MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo'      : genprocs_ltau_s_top,
'MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg' : genprocs_ltau_s_top,
'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg'     : genprocs_ltau_s_top,

'MC2016_Summer16_TTJets_powheg'        : genprocs_ltau_tt,
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4down' : genprocs_ltau_tt,
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4up'   : genprocs_ltau_tt,
'MC2016_Summer16_TTJets_powheg_fsrdown'        : genprocs_ltau_tt,
'MC2016_Summer16_TTJets_powheg_fsrup'          : genprocs_ltau_tt,
'MC2016_Summer16_TTJets_powheg_hdampDOWN'      : genprocs_ltau_tt,
'MC2016_Summer16_TTJets_powheg_hdampUP'        : genprocs_ltau_tt,
'MC2016_Summer16_TTJets_powheg_isrdown'        : genprocs_ltau_tt,
'MC2016_Summer16_TTJets_powheg_isrup'          : genprocs_ltau_tt,

'MC2016_Summer16_WWTo2L2Nu_powheg'             : genprocs_ltau_dibosons,
'MC2016_Summer16_WWToLNuQQ_powheg'             : genprocs_ltau_dibosons,
'MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin' : genprocs_ltau_dibosons,
'MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin'   : genprocs_ltau_dibosons,
'MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin'    : genprocs_ltau_dibosons,
'MC2016_Summer16_WZTo3LNu_powheg'              : genprocs_ltau_dibosons,
'MC2016_Summer16_ZZTo2L2Nu_powheg'             : genprocs_ltau_dibosons,
'MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin'    : genprocs_ltau_dibosons,
'MC2016_Summer16_W1Jets_madgraph'  : genprocs_ltau_wjets,
'MC2016_Summer16_W2Jets_madgraph'  : genprocs_ltau_wjets,
'MC2016_Summer16_W3Jets_madgraph'  : genprocs_ltau_wjets,
'MC2016_Summer16_W4Jets_madgraph'  : genprocs_ltau_wjets,
'MC2016_Summer16_WJets_madgraph'   : genprocs_ltau_wjets,
'MC2016_Summer16_WJets_amcatnlo'   : genprocs_ltau_wjets,

'SingleMuon'     : {'el': ('data', [])},
'SingleElectron' : {'el': ('data', [])},

'Data13TeV_SingleElectron2016B_03Feb2017_ver2' : {'el': ('data', [])},
'Data13TeV_SingleElectron2016C_03Feb2017_v1'   : {'el': ('data', [])},
'Data13TeV_SingleElectron2016D_03Feb2017_v1'   : {'el': ('data', [])},
'Data13TeV_SingleElectron2016E_03Feb2017_v1'   : {'el': ('data', [])},
'Data13TeV_SingleElectron2016F_03Feb2017_v1'   : {'el': ('data', [])},
'Data13TeV_SingleElectron2016G_03Feb2017_v1'   : {'el': ('data', [])},
'Data13TeV_SingleElectron2016H_03Feb2017_ver2' : {'el': ('data', [])},
'Data13TeV_SingleElectron2016H_03Feb2017_ver3' : {'el': ('data', [])},
'Data13TeV_SingleMuon2016B_03Feb2017_ver2' : {'el': ('data', [])},
'Data13TeV_SingleMuon2016C_03Feb2017_v1'   : {'el': ('data', [])},
'Data13TeV_SingleMuon2016D_03Feb2017_v1'   : {'el': ('data', [])},
'Data13TeV_SingleMuon2016E_03Feb2017_v1'   : {'el': ('data', [])},
'Data13TeV_SingleMuon2016F_03Feb2017_v1'   : {'el': ('data', [])},
'Data13TeV_SingleMuon2016G_03Feb2017_v1'   : {'el': ('data', [])},
'Data13TeV_SingleMuon2016H_03Feb2017_ver2' : {'el': ('data', [])},
'Data13TeV_SingleMuon2016H_03Feb2017_ver3' : {'el': ('data', [])},
}

#dtags_procs_el = {
#'MC2016_Summer16_TTJets_powheg'               : genprocs_eltau_tt,
#'MC2016_Summer16_SingleT_tW_5FS_powheg'       : genprocs_eltau_s_top,
#'MC2016_Summer16_SingleTbar_tW_5FS_powheg'    : genprocs_eltau_s_top,
#'MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo' : genprocs_eltau_s_top,
#'MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg' : genprocs_eltau_s_top,
#'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg' : genprocs_eltau_s_top,
#'MC2016_Summer16_TTJets_powheg_CUETP8M2T4down':  genprocs_eltau_tt,
#'MC2016_Summer16_TTJets_powheg_CUETP8M2T4up'  :  genprocs_eltau_tt,
#'MC2016_Summer16_TTJets_powheg_fsrdown'       :  genprocs_eltau_tt,
#'MC2016_Summer16_TTJets_powheg_fsrup'         :  genprocs_eltau_tt,
#'MC2016_Summer16_TTJets_powheg_hdampDOWN'     :  genprocs_eltau_tt,
#'MC2016_Summer16_TTJets_powheg_hdampUP'       :  genprocs_eltau_tt,
#'MC2016_Summer16_TTJets_powheg_isrdown'       :  genprocs_eltau_tt,
#'MC2016_Summer16_TTJets_powheg_isrup'         :  genprocs_eltau_tt,
#}


