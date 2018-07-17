from os import path as path
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
#LumiList.LumiList().getVLuminosityBlockRange()

# this ivars thing, whatever this is, will hold the names to the input/output files
ivars = VarParsing.VarParsing('analysis')


HLT_source = 'HLT' # 'HLT2'
withHLT = True

 # Data, file on netwokr -- let's see how cmsRun gets it
 #'root://cms-xrd-global.cern.ch///store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/001E3E7D-57EB-E611-8469-0CC47A7C35D2.root'
 #31 Aug present on CERN
 #'root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/00B474D3-ADEA-E611-9E30-D067E5F910F5.root'
# single top
 #'file:ST_tW_top_0C2044DB-0EC2-E611-8567-0CC47A7FC378.root'

# testing MET filters
input_files, isMC, dtag = ('root://eoscms//eos/cms///store/data/Run2016H/SingleMuon/MINIAOD/03Feb2017_ver2-v1/110000/00B474D3-ADEA-E611-9E30-D067E5F910F5.root',), False, 'Data13TeV_SingleMuon2016H_03Feb2017_ver2'

# RunB: 'root://eoscms//eos/cms///store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver1-v1/50000/3AFB9551-E6EB-E611-8EDA-0025905C3D98.root'
input_files, isMC, dtag = ('root://eoscms//eos/cms///store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver1-v1/50000/3AFB9551-E6EB-E611-8EDA-0025905C3D98.root',), False, 'Data13TeV_SingleMuon2016B_03Feb2017_ver1'

# single top files, jobs segfault on them now
input_files, isMC, dtag = ('root://eoscms//eos/cms///store/mc/RunIISummer16MiniAODv2/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/AA6B7439-A1C1-E611-BB03-FA163E27A6DC.root',), True, 'ST_tW'

 # WJets
#input_files, isMC, dtag = ('root://eoscms//eos/cms///store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/009CE684-45BB-E611-A261-001E67E6F8FA.root',), True, 'WJets'
input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/MC2016_Summer16_W0Jets_amcatnlo_2475B5FF-D7BE-E611-8E3A-0CC47A7FC350.root',), True, 'WJets'

 # TT for tau-rich events
# /eos/user/o/otoldaie/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_D4F3B681-A3BE-E611-913C-0CC47AC08816.root
input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_06437FA0-B2D8-E611-923D-02163E019DBD.root',), True, 'WZTo2L2Q'


input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/Data_SingleElectron_G_FAAE8AFC-32EB-E611-88E3-0CC47A6C186C.root',), False, 'SingleElectron2016G'

input_files, isMC, is2017rereco, dtag = ('file:/eos/user/o/otoldaie/Data_SingleMuon_B_Aug_D83FEA89-CC81-E711-83A8-008CFA197E74.root',), False, True,  'SingleMuon2016BAug'
input_files, isMC, is2017rereco, dtag = ('file:/eos/user/o/otoldaie/Data_SingleMuon_B_3AFB9551-E6EB-E611-8EDA-0025905C3D98.root',),     False, False, 'SingleMuon2016B'
input_files, isMC, is2017rereco, dtag = ('file:/eos/user/o/otoldaie/Data_SingleMuon_B_Aug_5ADCFC9D-4B81-E711-A19C-0CC47A4C8E8A.root',), False, True,  'SingleMuon2016BAug2'

# DY file
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/MC2015_Spring16_reHLT_DY-50_amcatnlo_00200284-F15C-E611-AA9B-002590574776.root",), True, 'DYJets2015_amcatnlo'
input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/MC2016_Summer16_DY50_amcatnlo_0EE8D393-D0DE-E611-9106-A4BF0101202F.root',), True, 'DYJets'

# QCD
input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/QCD2015_HT100to200_48346CCA-B351-E611-8851-C4346BC8D390.root',), True, "QCD2015-100-200"
input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/QCD2016_HT100to200_E4183151-F8BD-E611-8228-002590207E3C.root',), True, "QCD2016-100-200"

input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/TT_TuneCUETP8M2T4_13TeV-powheg-fsrdown-pythia8_9EAC0046-06B7-E611-AF76-141877411FED.root',), True, 'TTJets_fsrdown'
input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/TT_165F54A0-A3BE-E611-B3F7-0025905A606A.root',), True, 'TTJets2'


if any('2015' in infile for infile in input_files) or '2015' in dtag:
    HLT_source = 'HLT2'

record_scheme = 'signal' # 'tauID' # 'tauCands' #  Dilep MonitorHLT tauIDantiIso jets'

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

to_tag = True

output_file = '/afs/cern.ch/work/o/otoldaie/private/16/CMSSW_8_0_29/src/UserCode/NtuplerAnalyzer/NtuplerAnalyzer_test_METfilters%s_%s_%s.root' % ('OFF' if to_tag else 'ON', dtag, record_scheme)
print output_file
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
process.MessageLogger.cerr.FwkReport.reportEvery = 100

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



process.load("UserCode.NtuplerAnalyzer.CfiFile_cfi")
#process.demo.minTracks=1000

#process.ntupler.dtag = cms.string('MC2016_TT_powheg')
process.ntupler.isMC = cms.bool(isMC)
process.ntupler.is2017rereco = cms.bool(is2017rereco)
process.ntupler.isLocal = cms.bool(False)
process.ntupler.withHLT = cms.bool(withHLT)
process.ntupler.dtag = cms.string(dtag)
process.ntupler.HLT_source = cms.string(HLT_source)
#process.ntupler.elHLT_MC   = cms.string("HLT_Ele32_eta2p1_WPTight_Gsf_v8")
#process.ntupler.elHLT_Data = cms.string("HLT_Ele32_eta2p1_WPTight_Gsf_v*")

## test of DY 2015
#process.ntupler.muHLT_MC1   = cms.string("HLT_IsoMu24_v*"  )
#process.ntupler.muHLT_MC2   = cms.string("HLT_IsoTkMu24_v*")
#process.ntupler.muHLT_Data1 = cms.string("HLT_IsoMu24_v*"  )
#process.ntupler.muHLT_Data2 = cms.string("HLT_IsoTkMu24_v*")
#process.ntupler.elHLT_Data  = cms.string("HLT_Ele27_WPTight_Gsf_v*")
#process.ntupler.elHLT_MC    = cms.string("HLT_Ele27_WPTight_Gsf_v*")

#theLumiMask = path.expandvars("") # -- lumi should be handled via CRAB3, but for now I leave this config option available in Ntupler for local runs
#process.ntupler.lumisToProcess = LumiList.LumiList(filename = theLumiMask).getVLuminosityBlockRange()

# for LumiDump:
process.ntupler.input = cms.untracked.vstring(input_files)
process.ntupler.outfile = cms.string(output_file)

print 'record_scheme', record_scheme
if record_scheme:
    process.ntupler.record_tauID         = cms.bool('tauID'         in record_scheme)
    process.ntupler.record_tauCands      = cms.bool('tauCands'      in record_scheme)
    process.ntupler.record_tauIDantiIso  = cms.bool('tauIDantiIso'  in record_scheme)
    process.ntupler.record_bPreselection = cms.bool('bPreselection' in record_scheme)
    process.ntupler.record_MonitorHLT    = cms.bool('MonitorHLT'    in record_scheme)
    process.ntupler.record_ElMu          = cms.bool('ElMu'          in record_scheme)
    process.ntupler.record_Dilep         = cms.bool('Dilep'         in record_scheme)
    process.ntupler.record_jets          = cms.bool('jets'          in record_scheme)
    process.ntupler.record_signal        = cms.bool('signal'        in record_scheme)

#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.Tracer = cms.Service("Tracer")

### Set output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(ivars.outputFile)
)


# this is loaded in ../../RecoMET/METFilters/test/test_badChargedCandidateFilter.py
# before met filters
process.load('Configuration.StandardSequences.Services_cff')
# .. and I still get MissingParameter: Parameter 'BadChargedCandidateFilter' not found.

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
process.BadPFMuonFilter.taggingMode = cms.bool(to_tag)
process.BadPFMuonFilter.filter = cms.bool(not to_tag)

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.taggingMode = cms.bool(to_tag)
process.BadChargedCandidateFilter.filter = cms.bool(not to_tag)


if isMC:
    print "MC"
    process.p = cms.Path(
     process.BadPFMuonFilter *
     process.BadChargedCandidateFilter *
     process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel*process.bfragWgtProducer*
     process.ntupler)

else:
    print "data"
    process.p = cms.Path(
     process.BadPFMuonFilter *
     process.BadChargedCandidateFilter *
     process.ntupler)



