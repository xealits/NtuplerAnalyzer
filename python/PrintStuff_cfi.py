import FWCore.ParameterSet.Config as cms

printer = cms.EDAnalyzer('PrintStuff',
    isMC = cms.bool(False),
    HLT_source = cms.string("HLT"),
    )

