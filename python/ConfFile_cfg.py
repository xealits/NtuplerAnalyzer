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

input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_366459FB-DFC5-E611-BE5E-D8D385AF8AE2.root",), True, 'TTJets_amcatnlo'
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TT_TuneCUETP8M2T4_13TeV-powheg-fsrdown-pythia8_3071D964-03B7-E611-98C4-B083FED18596.root",), True, 'TTJets_fsrdown'
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TT_TuneCUETP8M2T4_13TeV-powheg-fsrup-pythia8_BE2542A8-C2D6-E611-BCB2-0CC47A6C0682.root",),   True, 'TTJets_fsrup'
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TT_TuneCUETP8M2T4_13TeV-powheg-isrdown-pythia8_244A6BD6-CFB5-E611-98C3-008CFA16615C.root",), True, 'TTJets_isrdown'
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TT_TuneCUETP8M2T4_13TeV-powheg-isrup-pythia8_14D0A1F0-A0BE-E611-BF78-FA163E795A09.root",),   True, 'TTJets_isrup'
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TT_TuneCUETP8M2T4_GluonMoveCRTune_13TeV-powheg-pythia8_B6EEA189-C940-E711-B6FB-E0071B6CAD00.root",), True, 'TTJets_GluonMoveCRTune'
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TT_TuneCUETP8M2T4_QCDbasedCRTune_erdON_13TeV-powheg-pythia8_BE256AF2-B926-E711-A36C-A0369FD20608.root",), True, 'TTJets_QCDbasedCRTune'
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TT_TuneCUETP8M2T4down_13TeV-powheg-pythia8_B87402B7-52C0-E611-95D4-B8CA3A70A410.root",),   True, 'TTJets_Tunedown'
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TT_TuneCUETP8M2T4up_13TeV-powheg-pythia8_BAD9F05B-C6B6-E611-9C69-0023AEEEB226.root",),     True, 'TTJets_Tuneup'
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TT_TuneEE5C_13TeV-powheg-herwigpp_6820BBE0-01B5-E611-85C8-002590E7E01A.root",), True, 'TTJets_herwigpp'
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TT_hdampDOWN_TuneCUETP8M2T4_13TeV-powheg-pythia8_F8364FC9-A5F2-E611-90F2-FA163E73666D.root",), True, 'TTJets_HDAMPDown'
input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/TT_hdampUP_TuneCUETP8M2T4_13TeV-powheg-pythia8_1EA7E8C9-C0F2-E611-85C3-FA163E65B77A.root",),   True, 'TTJets_HDAMPUp'

input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/TT_TuneCUETP8M2T4_13TeV-powheg-fsrdown-pythia8_9EAC0046-06B7-E611-AF76-141877411FED.root',), True, 'TTJets_fsrdown'

input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/TT_165F54A0-A3BE-E611-B3F7-0025905A606A.root',), True, 'TTJets2_CMSSW94_test'
input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/test-data-files/RunIIFall17MiniAOD_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8_MINIAODSIM_94X_mc2017_realistic_v10-v1_00000_9E2812DC-B7F4-E711-A3AC-0025905C54BA.root',), True, 'TTJets_2017_test'
# same file in public directory:

input_files, isMC, dtag = ('file:/afs/cern.ch/user/m/mmagheri/workspace/public/FEC01104-3B42-E811-8349-001E6739BB01.root',), True, 'TTJets_2017_test_crashing'

input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/test-data-files/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_RunIISummer16MiniAODv3_PUMoriond17_94X_mcRun2_asymptotic_v3-v2_MINIAODSIM_100000_4AF402CC-38ED-E811-9602-AC1F6B23C806.root',), True, 'TTJets_ttHtranche3_2016legacy'

input_files, isMC, dtag = ("file:/eos/user/o/otoldaie/test-data-files/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_RunIISummer16MiniAODv3_PUMoriond17_94X_mcRun2_asymptotic_v3-v1_110000_B22320DB-EBBC-E811-A987-14187763B811.root",), True, 'MC_TTJets_2016legacy'
input_files, isMC, dtag = ('file:/afs/cern.ch/user/o/otoldaie/public/RunIIFall17MiniAOD_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8_MINIAODSIM_94X_mc2017_realistic_v10-v1_00000_9E2812DC-B7F4-E711-A3AC-0025905C54BA.root',), True, 'MC_TTJets_2017legacy_test2'
input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/test-data-files/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIFall17MiniAODv2_PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_00000_523E450B-CB41-E811-AACA-001E6739B849.root',), True, 'MC_TTJets_2017legacy_test2'

input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/test-data-files/Data_SingleElectron_Run2016G_17Jul2018-v1_MINIAOD_50000_02DC57F5-A48B-E811-A98E-0CC47AF9B496.root',), False, 'Data_SingleElectron_Run2016G_2016legacy'
input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/test-data-files/Data_SingleMuon_Run2017D_31Mar2018-v1_MINIAOD_80000_1E703527-F436-E811-80A7-E0DB55FC1055.root',),     False, 'Data_SingleMuon_Run2017D_2017legacy'

dataset_reco_name = '2017legacy'

'''
running on
input_files, isMC, dtag = ('file:/afs/cern.ch/user/m/mmagheri/workspace/public/FEC01104-3B42-E811-8349-001E6739BB01.root',), True, 'TTJets_2017_test_crashing'

Requested UserFloat energyScaleEtUp is not available! Possible UserFloats are: 
ElectronMVAEstimatorRun2Fall17IsoV1Values ElectronMVAEstimatorRun2Fall17NoIsoV1Values ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values ElectronMVAEstimatorRun2Spring15Trig25nsV1Values ElectronMVAEstimatorRun2Spring15Trig50nsV1Values ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values ElectronMVAEstimatorRun2Spring16HZZV1Values ecalEnergyErrPostCorr ecalEnergyErrPreCorr ecalEnergyPostCorr ecalEnergyPreCorr ecalTrkEnergyErrPostCorr ecalTrkEnergyErrPreCorr ecalTrkEnergyPostCorr ecalTrkEnergyPreCorr energyScaleDown energyScaleGainDown energyScaleGainUp energyScaleStatDown energyScaleStatUp energyScaleSystDown energyScaleSystUp energyScaleUp energyScaleValue energySigmaDown energySigmaPhiDown energySigmaPhiUp energySigmaRhoDown energySigmaRhoUp energySigmaUp energySigmaValue energySmearNrSigma heepTrkPtIso 

on
input_files, isMC, dtag = ('file:/eos/user/o/otoldaie/test-data-files/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_RunIISummer16MiniAODv3_PUMoriond17_94X_mcRun2_asymptotic_v3-v2_MINIAODSIM_100000_4AF402CC-38ED-E811-9602-AC1F6B23C806.root',), True, 'TTJets_2016legacy'

-- legacy 2016 was split into sub-processes???

pat::Tau: the ID byVTightIsolationMVArun2017v2DBoldDMwLT2017 can't be found in this pat::Tau.
The available IDs are: 'againstElectronLooseMVA6' 'againstElectronMVA6Raw' 'againstElectronMVA6category' 'againstElectronMediumMVA6' 'againstElectronTightMVA6' 'againstElectronVLooseMVA6' 'againstElectronVTightMVA6' 'againstMuonLoose3' 'againstMuonTight3' 'byCombinedIsolationDeltaBetaCorrRaw3Hits' 'byIsolationMVArun2v1DBdR03oldDMwLTraw' 'byIsolationMVArun2v1DBnewDMwLTraw' 'byIsolationMVArun2v1DBoldDMwLTraw' 'byIsolationMVArun2v1PWdR03oldDMwLTraw' 'byIsolationMVArun2v1PWnewDMwLTraw' 'byIsolationMVArun2v1PWoldDMwLTraw' 'byLooseCombinedIsolationDeltaBetaCorr3Hits' 'byLooseIsolationMVArun2v1DBdR03oldDMwLT' 'byLooseIsolationMVArun2v1DBnewDMwLT' 'byLooseIsolationMVArun2v1DBoldDMwLT' 'byLooseIsolationMVArun2v1PWdR03oldDMwLT' 'byLooseIsolationMVArun2v1PWnewDMwLT' 'byLooseIsolationMVArun2v1PWoldDMwLT' 'byMediumCombinedIsolationDeltaBetaCorr3Hits' 'byMediumIsolationMVArun2v1DBdR03oldDMwLT' 'byMediumIsolationMVArun2v1DBnewDMwLT' 'byMediumIsolationMVArun2v1DBoldDMwLT' 'byMediumIsolationMVArun2v1PWdR03oldDMwLT' 'byMediumIsolationMVArun2v1PWnewDMwLT' 'byMediumIsolationMVArun2v1PWoldDMwLT' 'byPhotonPtSumOutsideSignalCone' 'byTightCombinedIsolationDeltaBetaCorr3Hits' 'byTightIsolationMVArun2v1DBdR03oldDMwLT' 'byTightIsolationMVArun2v1DBnewDMwLT' 'byTightIsolationMVArun2v1DBoldDMwLT' 'byTightIsolationMVArun2v1PWdR03oldDMwLT' 'byTightIsolationMVArun2v1PWnewDMwLT' 'byTightIsolationMVArun2v1PWoldDMwLT' 'byVLooseIsolationMVArun2v1DBdR03oldDMwLT' 'byVLooseIsolationMVArun2v1DBnewDMwLT' 'byVLooseIsolationMVArun2v1DBoldDMwLT' 'byVLooseIsolationMVArun2v1PWdR03oldDMwLT' 'byVLooseIsolationMVArun2v1PWnewDMwLT' 'byVLooseIsolationMVArun2v1PWoldDMwLT' 'byVTightIsolationMVArun2v1DBdR03oldDMwLT' 'byVTightIsolationMVArun2v1DBnewDMwLT' 'byVTightIsolationMVArun2v1DBoldDMwLT' 'byVTightIsolationMVArun2v1PWdR03oldDMwLT' 'byVTightIsolationMVArun2v1PWnewDMwLT' 'byVTightIsolationMVArun2v1PWoldDMwLT' 'byVVLooseIsolationMVArun2v1DBoldDMwLT' 'byVVTightIsolationMVArun2v1DBdR03oldDMwLT' 'byVVTightIsolationMVArun2v1DBnewDMwLT' 'byVVTightIsolationMVArun2v1DBoldDMwLT' 'byVVTightIsolationMVArun2v1PWdR03oldDMwLT' 'byVVTightIsolationMVArun2v1PWnewDMwLT' 'byVVTightIsolationMVArun2v1PWoldDMwLT' 'chargedIsoPtSum' 'chargedIsoPtSumdR03' 'decayModeFinding' 'decayModeFindingNewDMs' 'footprintCorrection' 'footprintCorrectiondR03' 'neutralIsoPtSum' 'neutralIsoPtSumWeight' 'neutralIsoPtSumWeightdR03' 'neutralIsoPtSumdR03' 'photonPtSumOutsideSignalCone' 'photonPtSumOutsideSignalConedR03' 'puCorrPtSum' .

now for 2017legacy TTJets:

pat::Electron: the ID cutBasedElectronID-Fall17-94X-V1-tight can't be found in this pat::Electron.
The available IDs are: 'cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose' 'cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium' 'cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight' 'cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto' 'cutBasedElectronID-Spring15-25ns-V1-standalone-loose' 'cutBasedElectronID-Spring15-25ns-V1-standalone-medium' 'cutBasedElectronID-Spring15-25ns-V1-standalone-tight' 'cutBasedElectronID-Spring15-25ns-V1-standalone-veto' 'cutBasedElectronID-Spring15-50ns-V2-standalone-loose' 'cutBasedElectronID-Spring15-50ns-V2-standalone-medium' 'cutBasedElectronID-Spring15-50ns-V2-standalone-tight' 'cutBasedElectronID-Spring15-50ns-V2-standalone-veto' 'eidLoose' 'eidRobustHighEnergy' 'eidRobustLoose' 'eidRobustTight' 'eidTight' 'heepElectronID-HEEPV60' 'mvaEleID-Spring15-25ns-Trig-V1-wp80' 'mvaEleID-Spring15-25ns-Trig-V1-wp90' 'mvaEleID-Spring15-25ns-nonTrig-V1-wp80' 'mvaEleID-Spring15-25ns-nonTrig-V1-wp90' 'mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose' 'mvaEleID-Spring15-50ns-Trig-V1-wp80' 'mvaEleID-Spring15-50ns-Trig-V1-wp90' .

new 2017 legacy TT file:
Requested UserFloat energyScaleEtUp is not available! Possible UserFloats are: 
ElectronMVAEstimatorRun2Fall17IsoV1Values ElectronMVAEstimatorRun2Fall17NoIsoV1Values ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values ElectronMVAEstimatorRun2Spring15Trig25nsV1Values ElectronMVAEstimatorRun2Spring15Trig50nsV1Values ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values ElectronMVAEstimatorRun2Spring16HZZV1Values ecalEnergyErrPostCorr ecalEnergyErrPreCorr ecalEnergyPostCorr ecalEnergyPreCorr ecalTrkEnergyErrPostCorr ecalTrkEnergyErrPreCorr ecalTrkEnergyPostCorr ecalTrkEnergyPreCorr energyScaleDown energyScaleGainDown energyScaleGainUp energyScaleStatDown energyScaleStatUp energyScaleSystDown energyScaleSystUp energyScaleUp energyScaleValue energySigmaDown energySigmaPhiDown energySigmaPhiUp energySigmaRhoDown energySigmaRhoUp energySigmaUp energySigmaValue energySmearNrSigma heepTrkPtIso 

-- ScaleEt was 2016legacy only:
https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing

also:
pat::Tau: the ID byVTightIsolationMVArun2017v2DBoldDMwLT2017 can't be found in this pat::Tau.
The available IDs are: 'againstElectronLooseMVA6' 'againstElectronMVA6Raw' 'againstElectronMVA6category' 'againstElectronMediumMVA6' 'againstElectronTightMVA6' 'againstElectronVLooseMVA6' 'againstElectronVTightMVA6' 'againstMuonLoose3' 'againstMuonTight3' 'byCombinedIsolationDeltaBetaCorrRaw3Hits' 'byIsolationMVArun2v1DBdR03oldDMwLTraw' 'byIsolationMVArun2v1DBnewDMwLTraw' 'byIsolationMVArun2v1DBoldDMwLTraw' 'byIsolationMVArun2v1PWdR03oldDMwLTraw' 'byIsolationMVArun2v1PWnewDMwLTraw' 'byIsolationMVArun2v1PWoldDMwLTraw' 'byLooseCombinedIsolationDeltaBetaCorr3Hits' 'byLooseIsolationMVArun2v1DBdR03oldDMwLT' 'byLooseIsolationMVArun2v1DBnewDMwLT' 'byLooseIsolationMVArun2v1DBoldDMwLT' 'byLooseIsolationMVArun2v1PWdR03oldDMwLT' 'byLooseIsolationMVArun2v1PWnewDMwLT' 'byLooseIsolationMVArun2v1PWoldDMwLT' 'byMediumCombinedIsolationDeltaBetaCorr3Hits' 'byMediumIsolationMVArun2v1DBdR03oldDMwLT' 'byMediumIsolationMVArun2v1DBnewDMwLT' 'byMediumIsolationMVArun2v1DBoldDMwLT' 'byMediumIsolationMVArun2v1PWdR03oldDMwLT' 'byMediumIsolationMVArun2v1PWnewDMwLT' 'byMediumIsolationMVArun2v1PWoldDMwLT' 'byPhotonPtSumOutsideSignalCone' 'byTightCombinedIsolationDeltaBetaCorr3Hits' 'byTightIsolationMVArun2v1DBdR03oldDMwLT' 'byTightIsolationMVArun2v1DBnewDMwLT' 'byTightIsolationMVArun2v1DBoldDMwLT' 'byTightIsolationMVArun2v1PWdR03oldDMwLT' 'byTightIsolationMVArun2v1PWnewDMwLT' 'byTightIsolationMVArun2v1PWoldDMwLT' 'byVLooseIsolationMVArun2v1DBdR03oldDMwLT' 'byVLooseIsolationMVArun2v1DBnewDMwLT' 'byVLooseIsolationMVArun2v1DBoldDMwLT' 'byVLooseIsolationMVArun2v1PWdR03oldDMwLT' 'byVLooseIsolationMVArun2v1PWnewDMwLT' 'byVLooseIsolationMVArun2v1PWoldDMwLT' 'byVTightIsolationMVArun2v1DBdR03oldDMwLT' 'byVTightIsolationMVArun2v1DBnewDMwLT' 'byVTightIsolationMVArun2v1DBoldDMwLT' 'byVTightIsolationMVArun2v1PWdR03oldDMwLT' 'byVTightIsolationMVArun2v1PWnewDMwLT' 'byVTightIsolationMVArun2v1PWoldDMwLT' 'byVVLooseIsolationMVArun2v1DBoldDMwLT' 'byVVTightIsolationMVArun2v1DBdR03oldDMwLT' 'byVVTightIsolationMVArun2v1DBnewDMwLT' 'byVVTightIsolationMVArun2v1DBoldDMwLT' 'byVVTightIsolationMVArun2v1PWdR03oldDMwLT' 'byVVTightIsolationMVArun2v1PWnewDMwLT' 'byVVTightIsolationMVArun2v1PWoldDMwLT' 'chargedIsoPtSum' 'chargedIsoPtSumdR03' 'decayModeFinding' 'decayModeFindingNewDMs' 'footprintCorrection' 'footprintCorrectiondR03' 'neutralIsoPtSum' 'neutralIsoPtSumWeight' 'neutralIsoPtSumWeightdR03' 'neutralIsoPtSumdR03' 'photonPtSumOutsideSignalCone' 'photonPtSumOutsideSignalConedR03' 'puCorrPtSum' .

-- taus are calculated on the fly in 2016legacy and 2017legacy
https://twiki.cern.ch/twiki/bin/view/CMS/TauIDRecommendation13TeV#2017v2_discriminators ...
look TEMPLATE for how it was done for ntupler production

'''


from os import environ

if 'INFILE' in environ:
    input_files = ('file:%s' % environ['INFILE'],)
    dtag = environ['DTAG']
    isMC = True

if any('2015' in infile for infile in input_files) or '2015' in dtag:
    HLT_source = 'HLT2'

record_scheme = ['tauID','Dilep'] # ('signal','recordAll') # ('signal',) # ('tauID',) # ('tauCands',) #  Dilep MonitorHLT tauIDantiIso jets'

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

output_file = './NtuplerAnalyzer_test_%s_METfilters%s_%s_%s.root' % ('MC' if isMC else 'Data', 'OFF' if to_tag else 'ON', dtag, '_'.join(record_scheme))
ivars.outputFile = output_file
# get and parse the command line arguments
ivars.parseArguments()

print input_files
print output_file

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



#process.load("UserCode.NtuplerAnalyzer.CfiFile_cfi")
import UserCode.NtuplerAnalyzer.CfiFile_cfi as std_config

process.ntupler = std_config.ntupler
#process.demo.minTracks=1000

print process.ntupler.tau_objs_name
#print process.ntupler['tau_objs_name']
print getattr(process.ntupler, 'tau_objs_name')

parameter_names_per_reco = std_config.parameter_names_per_reco

for param_name, param_defs in parameter_names_per_reco.items():
    param_def = param_defs['default']
    for same_definition_group in param_defs:
        if dataset_reco_name in same_definition_group:
            param_def = param_defs[same_definition_group]

    setattr(process.ntupler, param_name, param_def)

print dataset_reco_name, process.ntupler.tau_objs_name, process.ntupler.tau_VLoose_ID


# for legacy 2016 and 2017 the tau ID must be recomputed
# setting up the recomputation here

if dataset_reco_name in ('2016legacy', '2017legacy'):
    from UserCode.NtuplerAnalyzer.runTauIdMVA import *
    na = TauIDEmbedder(process, cms, # pass tour process object
        #cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff90"), # <-- a tricky place, WP90 is ??
        #cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff70"),
        tauIdDiscrMVA_2017_version = 'v2',
        debug  = True,
        toKeep = ["2017v2"] # pick the one you need: ["2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1"]
    )
    na.runTauID()


#process.ntupler.dtag = cms.string('MC2016_TT_powheg')
process.ntupler.isMC = cms.bool(isMC)
process.ntupler.is2016legacy = cms.bool(dataset_reco_name == '2016legacy')
process.ntupler.isLocal = cms.bool(False)
process.ntupler.withHLT = cms.bool(withHLT)
process.ntupler.dtag = cms.string(dtag)
process.ntupler.HLT_source = cms.string(HLT_source)
#process.ntupler.elHLT_MC   = cms.string("HLT_Ele32_eta2p1_WPTight_Gsf_v8")
#process.ntupler.elHLT_Data = cms.string("HLT_Ele32_eta2p1_WPTight_Gsf_v*")

# triggers
if '2017legacy' in dataset_reco_name:
    process.ntupler.muHLT_MC1   = cms.string("HLT_IsoMu24_v*"  )
    process.ntupler.muHLT_MC2   = cms.string("HLT_IsoTkMu24_v*")
    process.ntupler.muHLT_Data1 = cms.string("HLT_IsoMu24_v*"  )
    process.ntupler.muHLT_Data2 = cms.string("HLT_IsoTkMu24_v*")
    process.ntupler.elHLT_Data  = cms.string("HLT_Ele27_WPTight_Gsf_v*")
    process.ntupler.elHLT_MC    = cms.string("HLT_Ele27_WPTight_Gsf_v*")

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
    process.ntupler.record_all           = cms.bool('recordAll'     in record_scheme)

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


if dataset_reco_name in '2017legacy_2016legacy':
    run_path = process.rerunMvaIsolationSequence * process.NewTauIDsEmbedded * \
        process.BadPFMuonFilter * process.BadChargedCandidateFilter
else:
    run_path = process.BadPFMuonFilter * process.BadChargedCandidateFilter

if isMC:
    print "MC"
    # include bFrag systematics
    run_path = run_path * process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel*process.bfragWgtProducer
    #run_path = run_path * process.NewTauIDsEmbedded
    #process.p = cms.Path(
    # process.BadPFMuonFilter *
    # process.BadChargedCandidateFilter *
    # process.NewTauIDsEmbedded *
    # #process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel*process.bfragWgtProducer*
    # process.ntupler)

    #process.p = cms.Path(
    # process.rerunMvaIsolationSequence *
    # process.NewTauIDsEmbedded *
    # process.BadPFMuonFilter *
    # process.BadChargedCandidateFilter *
    # #process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel*process.bfragWgtProducer*
    # process.ntupler)

else:
    print "data"
    #process.p = cms.Path(
    # process.BadPFMuonFilter *
    # process.BadChargedCandidateFilter *
    # process.ntupler)

run_path = run_path * process.ntupler
process.p = cms.Path(run_path)

