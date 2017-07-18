import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
#import PhysicsTools.PythonAnalysis.LumiList as LumiList
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
#LumiList.LumiList().getVLuminosityBlockRange()

# this ivars thing, whatever this is, will hold the names to the input/output files
ivars = VarParsing.VarParsing('analysis')

ivars.inputFiles=(
 {input}
)

ivars.outputFile = '{outfile}'
# get and parse the command line arguments
ivars.parseArguments()


process = cms.Process("Demo")

# initialize MessageLogger and output report
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#    limit = cms.untracked.int32(-1)
#)
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(ivars.inputFiles))

# MC or Data is used everywhere
isMC = {isMC}

# nope, doing lumis another way
#theLumiMask = path.expandvars("") # for MC it defaults for "", and somehow the lumi-checker works well with that
#process.ntupler.lumisToProcess = LumiList.LumiList(filename = theLumiMask).getVLuminosityBlockRange()

# new lumi processing from
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePythonTips#Use_a_JSON_file_of_good_lumi_sec
# not sure how it works and how it counts lumisections processed
if not isMC:
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    JSONfile = '{lumiMask}'
    myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
    process.source.lumisToProcess.extend(myLumis)

# NTUPLER
process.load("UserCode.NtuplerAnalyzer.CfiFile_cfi")
process.ntupler.isMC = cms.bool({isMC})
#process.ntupler.dtag = cms.string('MC2016_TT_powheg')
process.ntupler.dtag = cms.string('{dtag}')


#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.Tracer = cms.Service("Tracer")

### Set output
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(ivars.outputFile)
)


# and this is supposedly met filters
#process.load("RecoMET.METFilters.metFilters_cff")
# ok, these are some old filters, they crash now

#process.p = cms.Path(process.metFilters * process.ntupler)
process.p = cms.Path(process.ntupler)



