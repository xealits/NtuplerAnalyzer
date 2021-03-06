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
    is2017rereco = cms.bool(False),
    isLocal = cms.bool(False),
    withHLT     = cms.bool(True  ),    # in 2015 some datasets (like QCD) don't have HLT
    HLT_source  = cms.string("HLT"  ), # in 2015 had HLT2 for reHLT datasets
    muHLT_MC1   = cms.string("HLT_IsoMu24_v4"  ),
    muHLT_MC2   = cms.string("HLT_IsoTkMu24_v4"),
    muHLT_Data1 = cms.string("HLT_IsoMu24_v*"  ), # the HLT pattern match
    muHLT_Data2 = cms.string("HLT_IsoTkMu24_v*"),
    elHLT_Data  = cms.string("HLT_Ele27_WPTight_Gsf_v*"),
    elHLT_MC    = cms.string("HLT_Ele27_WPTight_Gsf_v7"),
    lepMonitorHLT = cms.string("HLT_PFHT400_*"),

    jecDir          = cms.string('${CMSSW_BASE}/src/UserCode/NtuplerAnalyzer/data/jec/25ns/'),
    resolutionFile  = cms.string('${CMSSW_BASE}/src/UserCode/NtuplerAnalyzer/data/jec/25ns/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt'),
    scaleFactorFile = cms.string('${CMSSW_BASE}/src/UserCode/NtuplerAnalyzer/data/jec/25ns/Spring16_25nsV10_MC_SF_AK4PFchs.txt'),

    #mc2hessianCSV = cms.FileInPath('PhysicsTools/HepMCCandAlgos/data/NNPDF30_lo_as_0130_hessian_60.csv'), #MC2Hessian transformation matrix
    mc2hessianCSV = cms.FileInPath('PhysicsTools/HepMCCandAlgos/data/NNPDF30_nlo_as_0118_hessian_60.csv'), #MC2Hessian transformation matrix

    btag_threshold = cms.double(0.5426),
    # X80 2016 ReReco, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
    # Loose  0.5426
    # Medium 0.8484
    # Tight  0.9535
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
    tau_kino_cuts_eta = cms.double( 2.4)
)

