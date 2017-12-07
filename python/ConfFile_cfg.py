import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# this ivars thing, whatever this is, will hold the names to the input/output files
ivars = VarParsing.VarParsing('analysis')


 # Data, file on netwokr -- let's see how cmsRun gets it
 #'root://cms-xrd-global.cern.ch///store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/001E3E7D-57EB-E611-8469-0CC47A7C35D2.root'
 #31 Aug present on CERN
 #'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/00B474D3-ADEA-E611-9E30-D067E5F910F5.root'
# single top
 #'file:ST_tW_top_0C2044DB-0EC2-E611-8567-0CC47A7FC378.root'
 # DY file
input_files, isMC, dtag = ('root://eoscms//eos/cms///store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/0EE8D393-D0DE-E611-9106-A4BF0101202F.root',), True, 'DYJets'
 # WJets
input_files, isMC, dtag = ('root://eoscms//eos/cms///store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/009CE684-45BB-E611-A261-001E67E6F8FA.root',), True, 'WJets'
input_files, isMC, dtag = ('root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/00B474D3-ADEA-E611-9E30-D067E5F910F5.root',), False, 'Data13TeV_SingleMuon2016B_03Feb2017_ver2'
 # TT for tau-rich events
input_files, isMC, dtag = ('file:165F54A0-A3BE-E611-B3F7-0025905A606A.root',), True, 'TTJets'
 #'root://eoscms//eos/cms///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/26ABF488-A0BE-E611-BEEB-0CC47A4D7640.root'
 #'root://eoscms//eos/cms///store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/165F54A0-A3BE-E611-B3F7-0025905A606A.root'
 #'root://eoscms//eos/cms///store/data/Run2016G/SingleMuon/AOD/23Sep2016-v1/120000/D400618A-839C-E611-8E83-008CFAF73190.root' # no CTPPS here

# Exception Message:
#Principal::getByToken: Found zero products matching all criteria
#Looking for type: std::vector<pat::Muon>
#Looking for module label: slimmedMuons
#Looking for productInstanceName: 
# aha
#'root://eoscms//eos/cms///store/data/Run2016E/SingleElectron/AOD/18Apr2017-v1/120000/D00E73C7-EA3A-E711-8E6C-02163E011598.root' # AOD file on eos, testing accessing the same data in AOD
# this one complains about absence of some wrapper of CTPPS Dimond Digi stuff...................
#'root://eoscms//eos/cms///store/data/Run2016G/SingleMuon/MINIAOD/03Feb2017-v1/100000/02382B19-D1EA-E611-B2F9-0CC47ABAC11C.root'
#'root://xrootd.unl.edu//store/data/Run2015D/Charmonium/AOD/PromptReco-v4/000/258/159/00000/02D2473D-E06B-E511-80DA-02163E01418B.root'
#), True

ivars.inputFiles = input_files

output_file = '/afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/NtuplerAnalyzer_test_Anti-iso_%s.root' % dtag
ivars.outputFile = output_file
# get and parse the command line arguments
ivars.parseArguments()


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
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )





process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring(
    #    'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
    #)
    fileNames = cms.untracked.vstring(ivars.inputFiles)
)


# the fragmentation systematics calculator
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                    inputPruned = cms.InputTag("prunedGenParticles"),
                        inputPacked = cms.InputTag("packedGenParticles"),
)
from GeneratorInterface.RivetInterface.genParticles2HepMC_cfi import genParticles2HepMC
process.genParticles2HepMC = genParticles2HepMC.clone( genParticles = cms.InputTag("mergedGenParticles") )
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi")
process.particleLevel.excludeNeutrinosFromJetClustering = False
process.load('TopQuarkAnalysis.BFragmentationAnalyzer.bfragWgtProducer_cfi')



process.load("UserCode.NtuplerAnalyzer.CfiFile_cfi")
#process.demo.minTracks=1000

#process.ntupler.dtag = cms.string('MC2016_TT_powheg')
process.ntupler.isMC = cms.bool(isMC)
process.ntupler.isLocal = cms.bool(True)
process.ntupler.dtag = cms.string(dtag)

# for LumiDump:
process.ntupler.input = cms.untracked.vstring(input_files)
process.ntupler.outfile = cms.string(output_file)

record_scheme = 'tauID Dilep MonitorHLT tauIDantiIso jets'
if record_scheme:
    process.ntupler.record_tauID         = cms.bool('tauID'         in record_scheme)
    process.ntupler.record_tauIDantiIso  = cms.bool('tauIDantiIso'  in record_scheme)
    process.ntupler.record_bPreselection = cms.bool('bPreselection' in record_scheme)
    process.ntupler.record_MonitorHLT    = cms.bool('MonitorHLT'    in record_scheme)
    process.ntupler.record_ElMu          = cms.bool('ElMu'          in record_scheme)
    process.ntupler.record_Dilep         = cms.bool('Dilep'         in record_scheme)
    process.ntupler.record_jets          = cms.bool('jets'          in record_scheme)

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

#process.p = cms.Path(Flag_BadChargedCandidateFilter * Flag_BadPFMuonFilter * process.ntupler)
#process.p = cms.Path(BadChargedCandidateFilter * BadPFMuonFilter * process.ntupler)
#process.p = cms.Path(process.BadChargedCandidateFilter * process.BadPFMuonFilter * process.ntupler)
#process.p = cms.Sequence(process.BadChargedCandidateFilter * process.BadPFMuonFilter * process.ntupler)
#process.p = cms.Sequence(process.metFilters * process.ntupler)

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
 # the calculators for b->hadron fragmentation systematic weights
 process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel*process.bfragWgtProducer*
 process.ntupler)

#process.p = cms.Path(
# process.filter * process.ntupler)

#process.outpath = cms.EndPath(process.ntupler)
#process.schedule = cms.Schedule(process.p)


