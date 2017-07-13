import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('NtuplerAnalyzer',
           minTracks = cms.untracked.uint32(0)
)

ntupler = cms.EDAnalyzer('NtuplerAnalyzer',
           minTracks = cms.untracked.uint32(0)
)

