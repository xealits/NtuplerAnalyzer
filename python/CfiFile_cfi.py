import FWCore.ParameterSet.Config as cms

#demo = cms.EDAnalyzer('NtuplerAnalyzer',
#           minTracks = cms.untracked.uint32(0)
#)

ntupler = cms.EDAnalyzer('NtuplerAnalyzer' ,
    record_ElTau         = cms.bool(False)  ,
    record_MuTau         = cms.bool(False)  ,
    record_tauCands      = cms.bool(False) ,
    record_tauID         = cms.bool(True)  ,
    record_tauIDantiIso  = cms.bool(False) ,
    record_bPreselection = cms.bool(False) ,
    record_MonitorHLT    = cms.bool(False) ,
    record_ElMu          = cms.bool(False) ,
    record_Dilep         = cms.bool(False) ,
    record_jets          = cms.bool(False) ,
    record_signal        = cms.bool(False) ,
    record_all           = cms.bool(False) ,

    isMC = cms.bool(False),
    is2016legacy = cms.bool(False),
    isLocal = cms.bool(False),
    withHLT     = cms.bool(True  ),    # in 2015 some datasets (like QCD) don't have HLT
    HLT_source  = cms.string("HLT"  ), # in 2015 had HLT2 for reHLT datasets

    # same HLTs as in 2016, prescaled in 2017
    low_pt_muHLT_MC1   = cms.string("HLT_IsoMu24_eta2p1_v*"),
    low_pt_muHLT_MC2   = cms.string("HLT_IsoMu24_eta2p1_v*"),
    low_pt_muHLT_Data1 = cms.string("HLT_IsoMu24_eta2p1_v*"), # the HLT pattern match
    low_pt_muHLT_Data2 = cms.string("HLT_IsoMu24_eta2p1_v*"),

    low_pt_elHLT_Data  = cms.string("HLT_Ele27_WPTight_Gsf_v*"),
    low_pt_elHLT_MC    = cms.string("HLT_Ele27_WPTight_Gsf_v*"), # Ele27 is about 31 fb^-1 in 2017
    #elHLT_MC    = cms.string("HLT_Ele27_WPTight_Gsf_v*"), #MATTEO cambiato triggere elettroni MC, vediamo come va prima era v7

    # new HLTs

    low_pt32_elHLT  = cms.string("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v*"),

    # not prescaled
    # HT 150 is ok for ttbar, but not for DY tautau
    low_pt28_150HT_elHLT  = cms.string("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v*"),

    # not prescaled
    # PF jet 35 GeV is too high for DY tauh
    low_pt30_35PFJet_elHLT  = cms.string("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v*"),

    # main new HLTs: lowest unprescaled single-lepton
    muHLT_MC1   = cms.string("HLT_IsoMu27_v*"  ),
    muHLT_MC2   = cms.string("HLT_IsoMu27_v*"),
    muHLT_Data1 = cms.string("HLT_IsoMu27_v*"  ), # the HLT pattern match
    muHLT_Data2 = cms.string("HLT_IsoMu27_v*"),
    elHLT_Data  = cms.string("HLT_Ele35_WPTight_Gsf_v*"),
    elHLT_MC    = cms.string("HLT_Ele35_WPTight_Gsf_v*"), #TODO: figure out how to pass some parameters from command line?

    # dilepton HLTs
    # at least here we can get a lot of events
    elmuHLT_1 = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*"),
    elmuHLT_2 = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"),
    elmuHLT_3 = cms.string("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*"),
    elmuHLT_4 = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*"),

    # elel
    elelHLT_1 = cms.string("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*"),
    elelHLT_2 = cms.string("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"),

    # no mumu for now

    # leptau triggers
    eltauHLT  = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1*"),
    mutauHLT1 = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1*"),
    mutauHLT2 = cms.string("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1*"),


    lepMonitorHLT = cms.string("HLT_PFHT400_*"),

    #hlt_objects = cms.InputTag("selectedPatTrigger"), # 2016
    hlt_objects = cms.InputTag("slimmedPatTrigger"),  # legacy 2016 and 2017

    # files with JES corrections
    jecDir          = cms.string('${CMSSW_BASE}/src/UserCode/NtuplerAnalyzer/data/jec/25ns/'),
    # JER files
    # original 2016
    #resolutionFile  = cms.string('${CMSSW_BASE}/src/UserCode/NtuplerAnalyzer/data/jec/25ns/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt'),
    #scaleFactorFile = cms.string('${CMSSW_BASE}/src/UserCode/NtuplerAnalyzer/data/jec/25ns/Spring16_25nsV10_MC_SF_AK4PFchs.txt'),
    # legacy 2016
    resolutionFile  = cms.string('${CMSSW_BASE}/src/UserCode/NtuplerAnalyzer/data/jec/25ns/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt'),
    scaleFactorFile = cms.string('${CMSSW_BASE}/src/UserCode/NtuplerAnalyzer/data/jec/25ns/Summer16_25nsV1_MC_SF_AK4PFchs.txt'),

    #mc2hessianCSV = cms.FileInPath('PhysicsTools/HepMCCandAlgos/data/NNPDF30_lo_as_0130_hessian_60.csv'), #MC2Hessian transformation matrix
    mc2hessianCSV = cms.FileInPath('PhysicsTools/HepMCCandAlgos/data/NNPDF30_nlo_as_0118_hessian_60.csv'), #MC2Hessian transformation matrix

    btag_threshold = cms.double(0.2217),
    # X80 2016 ReReco, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
    # Loose  0.5426
    # Medium 0.8484
    # Tight  0.9535
    # 94X 2016 legacy ReReco, 17Jul2018 re-miniaod
    # Loose  0.2217
    # Medium 0.6321
    # Tight  0.8953
    el_kino_cuts_pt  = cms.double( 27.),
    el_kino_cuts_eta = cms.double( 2.5),
    el_veto_kino_cuts_pt  = cms.double( 15.),
    el_veto_kino_cuts_eta = cms.double( 2.5),
    mu_kino_cuts_pt  = cms.double( 24.),
    mu_kino_cuts_eta = cms.double( 2.5),
    mu_veto_kino_cuts_pt  = cms.double( 10.),
    mu_veto_kino_cuts_eta = cms.double( 2.5),

    jet_kino_cuts_pt  = cms.double( 20.),
    jet_kino_cuts_eta = cms.double( 2.5),
    tau_kino_cuts_pt  = cms.double( 20.),
    tau_kino_cuts_eta = cms.double( 2.4),

    # object names, variable among datasets and reco-s
    tau_objs_name = cms.string('slimmedTaus'),
    tau_VVLoose_ID = cms.string("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"),
    tau_VLoose_ID  = cms.string("byVLooseIsolationMVArun2017v2DBoldDMwLT2017" ),
    tau_Loose_ID   = cms.string("byLooseIsolationMVArun2017v2DBoldDMwLT2017"  ),
    tau_Medium_ID  = cms.string("byMediumIsolationMVArun2017v2DBoldDMwLT2017" ),
    tau_Tight_ID   = cms.string("byTightIsolationMVArun2017v2DBoldDMwLT2017"  ),
    tau_VTight_ID  = cms.string("byVTightIsolationMVArun2017v2DBoldDMwLT2017" ),
    tau_VVTight_ID = cms.string("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"),
)

# per-dataset-reco parameters defs
parameter_names_per_reco = {
'tau_objs_name':  {'2017legacy_2016legacy': cms.string('NewTauIDsEmbedded'), 'default': cms.string('slimmedTaus'),},
'tau_VVLoose_ID': {'2017legacy_2016legacy': cms.string("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"), 'default': ''},
'tau_VLoose_ID' : {'2017legacy_2016legacy': cms.string("byVLooseIsolationMVArun2017v2DBoldDMwLT2017" ), 'default': cms.string("byVLooseIsolationMVArun2v1DBoldDMwLT")},
'tau_Loose_ID'  : {'2017legacy_2016legacy': cms.string("byLooseIsolationMVArun2017v2DBoldDMwLT2017"  ), 'default': cms.string("byLooseIsolationMVArun2v1DBoldDMwLT") },
'tau_Medium_ID' : {'2017legacy_2016legacy': cms.string("byMediumIsolationMVArun2017v2DBoldDMwLT2017" ), 'default': cms.string("byMediumIsolationMVArun2v1DBoldDMwLT")},
'tau_Tight_ID'  : {'2017legacy_2016legacy': cms.string("byTightIsolationMVArun2017v2DBoldDMwLT2017"  ), 'default': cms.string("byTightIsolationMVArun2v1DBoldDMwLT") },
'tau_VTight_ID' : {'2017legacy_2016legacy': cms.string("byVTightIsolationMVArun2017v2DBoldDMwLT2017" ), 'default': cms.string("byVTightIsolationMVArun2v1DBoldDMwLT")},
'tau_VVTight_ID': {'2017legacy_2016legacy': cms.string("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"), 'default': ''},
}

