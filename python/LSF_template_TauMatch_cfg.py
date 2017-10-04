import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# this ivars thing, whatever this is, will hold the names to the input/output files
ivars = VarParsing.VarParsing('analysis')

ivars.inputFiles=(
 {input}
)

ivars.outputFile = '{outfile}'
# get and parse the command line arguments
ivars.parseArguments()



# MC or Data is used everywhere
isMC = {isMC}
dtag = '{dtag}'

process = cms.Process("Demo")


# setting GlobalTag for MC Morion2017:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Run2_Moriond17_re_digi_reco_camp
# data:
# Produced with: 8_0_26_patch1; Global tag: 80X_dataRun2_2016SeptRepro_v7 (eras B-G) 80X_dataRun2_Prompt_v16 (era H); the global tags are an update the 23Sep20216 and PromptReco ones to includes the 23Sep20216 V3 JECs on top .
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('80X_mcRun2_asymptotic_2016_TrancheIV_v6')

# some feature for tracks?
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
#process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool( True ),
#    #SkipEvent = cms.untracked.vstring('ProductNotFound')
#)

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
##process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    #limit = cms.untracked.int32(-1)
#    limit = cms.untracked.int32(10000)
#)

#process.MessageLogger = cms.Service("MessageLogger",
#       destinations   = cms.untracked.vstring(
#                                             'detailedInfo'
#                                               ,'critical'
#                                               ,'cerr'
#                    ),
#       ,critical       = cms.untracked.PSet(
#                       , threshold = cms.untracked.string('ERROR') 
#        ),
#       detailedInfo   = cms.untracked.PSet(
#                      threshold = cms.untracked.string('INFO') 
#       ),
#       cerr           = cms.untracked.PSet(
#                       threshold  = cms.untracked.string('WARNING') 
#        )
#)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )





process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring(
    #    'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
    #)
    fileNames = cms.untracked.vstring(ivars.inputFiles)
)

process.load("UserCode.NtuplerAnalyzer.TauMatch_cfi")
#process.demo.minTracks=1000

#process.tauMatch.dtag = cms.string('MC2016_TT_powheg')
process.tauMatch.isMC = cms.bool(isMC)
process.tauMatch.isLocal = cms.bool(True)
process.tauMatch.dtag = cms.string(dtag)

# for LumiDump:
process.tauMatch.input = cms.untracked.vstring({input})
process.tauMatch.outfile = cms.string('{outfile}')

record_scheme = 'tauID'
if record_scheme:
    process.tauMatch.record_tauID         = cms.bool('tauID'         in record_scheme)
    process.tauMatch.record_bPreselection = cms.bool('bPreselection' in record_scheme)
    process.tauMatch.record_MonitorHLT    = cms.bool('MonitorHLT'    in record_scheme)
    process.tauMatch.record_ElMu          = cms.bool('ElMu'          in record_scheme)
    process.tauMatch.record_Dilep         = cms.bool('Dilep'         in record_scheme)

#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.Tracer = cms.Service("Tracer")

### Set output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(ivars.outputFile)
)


#Testing MET filters
#process.load("RecoMET.METFilters.metFilters_cff")
#process.load("PhysicsTools.PatAlgos.slimming.metFilterPaths_cff")
#from PhysicsTools.PatAlgos.slimming.metFilterPaths_cff import *
process.load("RecoMET.METFilters.metFilters_cff") # back to this

#process.p = cms.Path(Flag_BadChargedCandidateFilter * Flag_BadPFMuonFilter * process.tauMatch)
#process.p = cms.Path(BadChargedCandidateFilter * BadPFMuonFilter * process.tauMatch)
#process.p = cms.Path(process.BadChargedCandidateFilter * process.BadPFMuonFilter * process.tauMatch)
#process.p = cms.Sequence(process.BadChargedCandidateFilter * process.BadPFMuonFilter * process.tauMatch)
#process.p = cms.Sequence(process.metFilters * process.tauMatch)

# another way:
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_the_Bad_Charged_Hadro

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.taggingMode = cms.bool(True)
#process.BadPFMuonFilter.filter = cms.bool(True)

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)
#process.BadChargedCandidateFilter.filter = cms.bool(True)


process.p = cms.Path(
 #process.TransientTrackBuilderESProducer +
 process.BadPFMuonFilter *
 process.BadChargedCandidateFilter *
 process.tauMatch)

#process.p = cms.Path(
# process.filter * process.tauMatch)

#process.outpath = cms.EndPath(process.tauMatch)
#process.schedule = cms.Schedule(process.p)



