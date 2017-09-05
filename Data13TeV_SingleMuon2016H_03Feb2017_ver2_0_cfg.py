from os import path as path

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
#import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
#LumiList.LumiList().getVLuminosityBlockRange()


# this ivars thing, whatever this is, will hold the names to the input/output files
ivars = VarParsing.VarParsing('analysis')

ivars.inputFiles=(
 'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/00633FF0-85EA-E611-811C-001E674FB25C.root',
'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/006517CB-8AEA-E611-8CF6-0CC47AC08BD4.root',
'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/008BE852-89EA-E611-A265-001E67396AEA.root',
'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/009A586B-A6EA-E611-B612-001E674FD1DD.root',
'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/00B474D3-ADEA-E611-9E30-D067E5F910F5.root',
'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/00BC6AB7-94EA-E611-BC4E-0CC47AD9908C.root',
'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/0226548B-7AEA-E611-B312-001E674FBC11.root',
'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/06161A14-88EA-E611-8994-001E674FBC11.root',
'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/069ADB63-89EA-E611-8781-001E674FB153.root',
'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/081F8997-ADEA-E611-AE36-002481CFE1EC.root'
)

ivars.outputFile = '/afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_25/src/UserCode/NtuplerAnalyzer/Data13TeV_SingleMuon2016H_03Feb2017_ver2_0.root'
# get and parse the command line arguments
ivars.parseArguments()


# MC or Data is used everywhere
isMC = False
dtag = 'Data13TeV_SingleMuon2016H_03Feb2017_ver2'

process = cms.Process("Demo")

# setting GlobalTag for MC Morion2017:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Run2_Moriond17_re_digi_reco_camp
# data:
# Produced with: 8_0_26_patch1;
# Global tag: 80X_dataRun2_2016SeptRepro_v7 (eras B-G)
# 80X_dataRun2_Prompt_v16 (era H);
# the global tags are an update the 23Sep20216 and PromptReco ones to includes the 23Sep20216 V3 JECs on top .
# in principle, now MET filters work without the GlobalTag anyway...
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if isMC:
    process.GlobalTag.globaltag = cms.string('80X_mcRun2_asymptotic_2016_TrancheIV_v6')
else:
    # so, era H has different global tag
    process.GlobalTag.globaltag = cms.string('80X_dataRun2_Prompt_v16') if '2016H' in dtag else cms.string('80X_dataRun2_2016SeptRepro_v7')

# some feature for tracks?
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")


# initialize MessageLogger and output report
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    limit = cms.untracked.int32(-1)
#)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
# reports N events tried and passed for each module in the filter-analyzer path
# example:
# TrigReport ---------- Module Summary ------------
# TrigReport    Visited   Executed     Passed     Failed      Error Name
# TrigReport       1025       1025       1025          0          0 BadChargedCandidateFilter
# TrigReport       1025       1025       1025          0          0 BadPFMuonFilter
# TrigReport       1025       1025       1025          0          0 TriggerResults
# TrigReport       1025       1025       1025          0          0 ntupler



#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(ivars.inputFiles))

# nope, doing lumis another way
#theLumiMask = path.expandvars("") # for MC it defaults for "", and somehow the lumi-checker works well with that
#process.ntupler.lumisToProcess = LumiList.LumiList(filename = theLumiMask).getVLuminosityBlockRange()

# new lumi processing from
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePythonTips#Use_a_JSON_file_of_good_lumi_sec
# not sure how it works and how it counts lumisections processed
if not isMC:
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    JSONfile = '/afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_25/src/UserCode/ttbar-leptons-80X/analysis/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
    myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
    process.source.lumisToProcess.extend(myLumis)
# it's not clear how to get the output from this: which LumiSecs have actually been processed (and which were not due to job crashes etc)
# maybe should use the old utility from llvv_fwk

# NTUPLER
process.load("UserCode.NtuplerAnalyzer.CfiFile_cfi")
process.ntupler.isMC = cms.bool(isMC)
#process.ntupler.dtag = cms.string('MC2016_TT_powheg')
process.ntupler.dtag = cms.string(dtag)

theLumiMask = path.expandvars("/afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_25/src/UserCode/ttbar-leptons-80X/analysis/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt")
process.ntupler.lumisToProcess = LumiList.LumiList(filename = theLumiMask).getVLuminosityBlockRange()

record_scheme = 'Dilep'
if record_scheme:
    record_tauID         = cms.bool('tauID'         in record_scheme)
    record_bPreselection = cms.bool('bPreselection' in record_scheme)
    record_MonitorHLT    = cms.bool('MonitorHLT'    in record_scheme)
    record_ElMu          = cms.bool('ElMu'          in record_scheme)
    record_Dilep         = cms.bool('Dilep'         in record_scheme)



#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.Tracer = cms.Service("Tracer")

### Set output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(ivars.outputFile)
)


# and this is supposedly met filters
#process.load("RecoMET.METFilters.metFilters_cff") # this loads the recommended (hopefully) metFilters sequence, together with BadPFMuon and BadChHadron

#process.p = cms.Path(process.metFilters * process.ntupler)
#process.p = cms.Path(process.ntupler)
#process.p = cms.Sequence(process.metFilters * process.ntupler)

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)


process.p = cms.Path(
 process.BadPFMuonFilter *
 process.BadChargedCandidateFilter *
 process.ntupler)




