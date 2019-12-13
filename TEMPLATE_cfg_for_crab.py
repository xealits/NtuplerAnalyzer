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
 'file:165F54A0-A3BE-E611-B3F7-0025905A606A.root',
 #'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/00B474D3-ADEA-E611-9E30-D067E5F910F5.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0806AB92-99BE-E611-9ECD-0025905A6138.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/165F54A0-A3BE-E611-B3F7-0025905A606A.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/18E31463-B3BE-E611-B6A3-0CC47A4D7678.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/22B7C5F1-A4BE-E611-AC0A-0025905B85AE.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/26ABF488-A0BE-E611-BEEB-0CC47A4D7640.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/38AB0251-9FBE-E611-A8C2-0025905A60DE.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/3CCF34DF-9DBE-E611-9512-0025905B858E.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/3E18521B-A4BE-E611-8843-0025905A607E.root'
)


# MC or Data is used everywhere
isMC = {isMC}   # PARAMETER
dtag = '{dtag}' # PARAMETER
dataset_reco_name = '{dataset_reco_name}'

#ivars.outputFile = '/afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/MC2016_Summer16_TTJets_powheg_0.root'
ivars.outputFile = dtag + '.root'
# get and parse the command line arguments
ivars.parseArguments()

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

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
if dataset_reco_name == '2017legacy':
    if isMC:
        #process.GlobalTag.globaltag = cms.string('94X_mc2017_realistic_v17')
        process.GlobalTag.globaltag = cms.string('80X_mcRun2_asymptotic_2016_TrancheIV_v6')
    else:
        # for some reason this tag crashes TransientTrackBuilder
        #process.GlobalTag.globaltag = cms.string('94X_dataRun2_v11')
        process.GlobalTag.globaltag = cms.string('80X_dataRun2_2016SeptRepro_v7')

else:
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



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) # all events

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(ivars.inputFiles))

# nope, doing lumis another way
#theLumiMask = path.expandvars("") # for MC it defaults for "", and somehow the lumi-checker works well with that
#process.ntupler.lumisToProcess = LumiList.LumiList(filename = theLumiMask).getVLuminosityBlockRange()

# it is handled by CRAB
## new lumi processing from
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePythonTips#Use_a_JSON_file_of_good_lumi_sec
## not sure how it works and how it counts lumisections processed
#if not isMC:
#    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#    JSONfile = ''
#    myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
#    process.source.lumisToProcess.extend(myLumis)
## it's not clear how to get the output from this: which LumiSecs have actually been processed (and which were not due to job crashes etc)
## maybe should use the old utility from llvv_fwk


# the fragmentation systematics calculator
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

#from doc on particlelevel https://twiki.cern.ch/twiki/bin/viewauth/CMS/ParticleLevelProducer
process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi") 

# from b-frag example https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer
#process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
#                    inputPruned = cms.InputTag("prunedGenParticles"),
#                        inputPacked = cms.InputTag("packedGenParticles"),
#)
#from GeneratorInterface.RivetInterface.genParticles2HepMC_cfi import genParticles2HepMC
#process.genParticles2HepMC = genParticles2HepMC.clone( genParticles = cms.InputTag("mergedGenParticles") )
#process.load("GeneratorInterface.RivetInterface.particleLevel_cfi")
#process.particleLevel.excludeNeutrinosFromJetClustering = False
process.load('TopQuarkAnalysis.BFragmentationAnalyzer.bfragWgtProducer_cfi')


# NTUPLER
#process.load("UserCode.NtuplerAnalyzer.CfiFile_cfi")
import UserCode.NtuplerAnalyzer.CfiFile_cfi as std_config

process.ntupler = std_config.ntupler

# set name customizations per dataset per reco
parameter_names_per_reco = std_config.parameter_names_per_reco

for param_name, param_defs in parameter_names_per_reco.items():
    param_def = param_defs['default']
    for same_definition_group in param_defs:
        if dataset_reco_name in same_definition_group:
            param_def = param_defs[same_definition_group]

    setattr(process.ntupler, param_name, param_def)


process.ntupler.isMC = cms.bool(isMC)
process.ntupler.is2016legacy = cms.bool(dataset_reco_name == '2016legacy')
#process.ntupler.dtag = cms.string('MC2016_TT_powheg')
process.ntupler.dtag = cms.string(dtag)

process.ntupler.isLocal = cms.bool(False) # LSF submition is local

process.ntupler.withHLT = cms.bool({withHLT})

# the triggers in 2017:
#process.ntupler.elHLT_MC   = cms.string("HLT_Ele27_WPTight_Gsf_v14")
#process.ntupler.muHLT_MC1   = cms.string("HLT_IsoMu24_v*"  )
#process.ntupler.muHLT_MC2   = cms.string("HLT_IsoTkMu24_v*")

# for LumiDump (to be scraped):
process.ntupler.input = cms.untracked.vstring(
'file:165F54A0-A3BE-E611-B3F7-0025905A606A.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0693E0E7-97BE-E611-B32F-0CC47A78A3D8.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/0806AB92-99BE-E611-9ECD-0025905A6138.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/165F54A0-A3BE-E611-B3F7-0025905A606A.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/18E31463-B3BE-E611-B6A3-0CC47A4D7678.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/22B7C5F1-A4BE-E611-AC0A-0025905B85AE.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/26ABF488-A0BE-E611-BEEB-0CC47A4D7640.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/38AB0251-9FBE-E611-A8C2-0025905A60DE.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/3CCF34DF-9DBE-E611-9512-0025905B858E.root',
#'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/3E18521B-A4BE-E611-8843-0025905A607E.root'
)
process.ntupler.outfile = cms.string('/afs/cern.ch/user/o/otoldaie/work/private/16/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/MC2016_Summer16_TTJets_powheg_test.root')

## triggers
#if '2017legacy' in dataset_reco_name:
#    process.ntupler.muHLT_MC1   = cms.string("HLT_IsoMu24_v*"  )
#    process.ntupler.muHLT_MC2   = cms.string("HLT_IsoTkMu24_v*")
#    process.ntupler.muHLT_Data1 = cms.string("HLT_IsoMu24_v*"  )
#    process.ntupler.muHLT_Data2 = cms.string("HLT_IsoTkMu24_v*")
#    process.ntupler.elHLT_Data  = cms.string("HLT_Ele27_WPTight_Gsf_v*")
#    process.ntupler.elHLT_MC    = cms.string("HLT_Ele27_WPTight_Gsf_v*")


theLumiMask = path.expandvars("") # -- lumi should be handled via CRAB3, but for now I leave this config option available in Ntupler for local runs
process.ntupler.lumisToProcess = LumiList.LumiList(filename = theLumiMask).getVLuminosityBlockRange()

record_scheme = {record_scheme} # PARAMETER
if record_scheme:
    process.ntupler.record_ElTau         = cms.bool('ElTau'         in record_scheme)
    process.ntupler.record_MuTau         = cms.bool('MuTau'         in record_scheme)
    process.ntupler.record_tauCands      = cms.bool('tauCands'      in record_scheme)
    process.ntupler.record_tauID         = cms.bool('tauID'         in record_scheme)
    process.ntupler.record_tauIDantiIso  = cms.bool('tauIDantiIso'  in record_scheme)
    process.ntupler.record_bPreselection = cms.bool('bPreselection' in record_scheme)
    process.ntupler.record_MonitorHLT    = cms.bool('MonitorHLT'    in record_scheme)
    process.ntupler.record_ElMu          = cms.bool('ElMu'          in record_scheme)
    process.ntupler.record_Dilep         = cms.bool('Dilep'         in record_scheme)
    process.ntupler.record_jets          = cms.bool('jets'          in record_scheme)
    process.ntupler.record_signal        = cms.bool('signal'        in record_scheme)
    process.ntupler.record_lepTauVL_b    = cms.bool('lepTauVLid'    in record_scheme)
    process.ntupler.record_ElMu_b        = cms.bool('el_mu_b'       in record_scheme)



#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.Tracer = cms.Service("Tracer")

### Set output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(ivars.outputFile)
)


# and this is supposedly met filters
process.load("RecoMET.METFilters.metFilters_cff") # this loads the recommended (hopefully) metFilters sequence, together with BadPFMuon and BadChHadron

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


# for legacy 2016 the tau ID must be recomputed
# setting up the recomputation here

from UserCode.NtuplerAnalyzer.runTauIdMVA import *

na = TauIDEmbedder(process, cms, # pass tour process object
    #cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff90"), # <-- a tricky place, WP90 is ??
    #cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff70"),
    tauIdDiscrMVA_2017_version = 'v2',
    debug  = True,
    toKeep = ["2017v2"] # pick the one you need: ["2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1"]
)
na.runTauID()


"""
if isMC:
    print "MC"
    process.p = cms.Path(
     process.rerunMvaIsolationSequence *
     process.NewTauIDsEmbedded *
     process.BadPFMuonFilter *
     process.BadChargedCandidateFilter *
     process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel*process.bfragWgtProducer*
     process.ntupler)

else:
    print "data"
    process.p = cms.Path(
     process.rerunMvaIsolationSequence *
     process.NewTauIDsEmbedded *
     process.BadPFMuonFilter *
     process.BadChargedCandidateFilter *
     process.ntupler)
"""


if dataset_reco_name in '2017legacy_2016legacy':
    # on the fly tau ID in 2016legacy and 2017legacy
    run_path = process.rerunMvaIsolationSequence * process.NewTauIDsEmbedded * \
        process.BadPFMuonFilter * process.BadChargedCandidateFilter
else:
    run_path = process.BadPFMuonFilter * process.BadChargedCandidateFilter

if isMC:
    print "MC"
    # include bFrag systematics
    run_path = run_path * process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel*process.bfragWgtProducer

else:
    print "Data"

run_path = run_path * process.ntupler
process.p = cms.Path(run_path)

