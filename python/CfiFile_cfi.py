import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('NtuplerAnalyzer',
           minTracks = cms.untracked.uint32(0)
)

ntupler = cms.EDAnalyzer('NtuplerAnalyzer',
    isMC = cms.untracked.bool(True),
    muHLT_MC1   = cms.string("HLT_IsoMu24_v4"  ),
    muHLT_MC2   = cms.string("HLT_IsoTkMu24_v4"),
    muHLT_Data1 = cms.string("HLT_IsoMu24_v*"  ),
    muHLT_Data2 = cms.string("HLT_IsoTkMu24_v*"),
    elHLT_Data  = cms.string("HLT_PFJet140_v"  "HLT_Ele27_WPTight_Gsf_v*"),
    elHLT_MC    = cms.string("HLT_Ele27_WPTight_Gsf_v7")
)

