// -*- C++ -*-
//
// Package:    UserCode/NtuplerAnalyzer
// Class:      NtuplerAnalyzer
// 
/**\class NtuplerAnalyzer NtuplerAnalyzer.cc UserCode/NtuplerAnalyzer/plugins/NtuplerAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Oleksii Toldaiev
//         Created:  Thu, 13 Jul 2017 00:08:32 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// for TFileService -- to get to the output file
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


// lepton ID/Iso prescriptions
#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
// lumiUtils
#include "UserCode/llvv_fwk/interface/LumiUtils.h"
// couple functions processing leptons
#include "UserCode/ttbar-leptons-80X/interface/ProcessingMuons.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingElectrons.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingTaus.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingJets.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingBJets.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingHLT.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingGenParticles.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingDRCleaning.h"



// ptocedures for initializations
namespace utils
	{
	namespace cmssw
		{
		// TODO: it is the same jetCorrector as in MacroUtils, only Fall_ prefix is set
		// Fall15_25nsV2_
		FactorizedJetCorrector* getJetCorrector(TString baseDir, TString pf, bool isMC)
			{
			gSystem->ExpandPathName(baseDir);
			//TString pf(isMC ? "MC" : "DATA");
			// TString pf("Fall15_25nsV2_");
			//pf += (isMC ? "MC" : "DATA");

			//order matters: L1 -> L2 -> L3 (-> Residuals)
			std::vector<std::string> jetCorFiles;
			std::cout << baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt" << std::endl;
			jetCorFiles.push_back((baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt").Data());
			jetCorFiles.push_back((baseDir+"/"+pf+"_L2Relative_AK4PFchs.txt").Data());
			jetCorFiles.push_back((baseDir+"/"+pf+"_L3Absolute_AK4PFchs.txt").Data());
			if(!isMC) jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
			// now there is a practically empty file Fall15_25nsV2_MC_L2L3Residual_AK4PFchs.txt
			// adding the run on it anyway
			//jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
			// it is dummy/empty file for MC and apparently is is not used
			// but in v13.1 it seemed to influence selection a bit
			// adding it for v13.4 -- will test later without it
			// and removing in 13.7 test -- compare with 13.4 & 13.4_repeat

			//init the parameters for correction
			std::vector<JetCorrectorParameters> corSteps;
			for(size_t i=0; i<jetCorFiles.size(); i++) corSteps.push_back(JetCorrectorParameters(jetCorFiles[i]));
			//return the corrector
			return new FactorizedJetCorrector(corSteps);
			}
		}
	}


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class NtuplerAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit NtuplerAnalyzer(const edm::ParameterSet&);
      ~NtuplerAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      unsigned int minTracks_;

	/* the new, never documented in Workbook stuff with tikens
	 * delivered by this guy:
	 * https://twiki.cern.ch/twiki/bin/view/Main/HanDle
	 */
	/*declare in the class */
	//edm::EDGetTokenT<reco::TrackCollection> tracks_;
	edm::EDGetTokenT<pat::MuonCollection> muons_;
	edm::EDGetTokenT<pat::ElectronCollection> electrons_;
	edm::EDGetTokenT<pat::TauCollection> taus_;
	edm::EDGetTokenT<reco::VertexCollection> vtx_;
	edm::EDGetTokenT<double> rho_;
	edm::EDGetTokenT<edm::TriggerResults> trigResults_;
	//triggerObjectsHandle.getByLabel(ev, "selectedPatTrigger");
	//iEvent.getByToken(triggerObjects_, triggerObjectsHandle);
	edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> triggerObjects_;
	edm::EDGetTokenT<LHEEventProduct> lheEPToken_;
	edm::EDGetTokenT< std::vector < PileupSummaryInfo > > puInfo_;
	edm::EDGetTokenT< std::vector < PileupSummaryInfo > > puInfo2_;
	edm::EDGetTokenT<reco::GenParticleCollection> genParticle_;

	edm::EDGetTokenT<pat::METCollection> mets1_, mets2_, mets_uncorrected_;

	edm::EDGetTokenT<vector<reco::GenJet>> genJets_;
	//genJetsHandle.getByLabel(ev, "slimmedGenJets");
	edm::EDGetTokenT<pat::JetCollection> jets_;
	//jetsHandle.getByLabel(ev, "slimmedJets");

	TString dtag;
	bool isMC;
	string  muHLT_MC1  , muHLT_MC2  ,
		muHLT_Data1, muHLT_Data2,
		elHLT_Data , elHLT_MC   ;

	jet_id    jetID;
	pu_jet_id jetPUID;

	TString jecDir;
	TString TjetResolutionFileName;
	TString TjetResolutionSFFileName;

	FactorizedJetCorrector *jesCor;
	JetCorrectionUncertainty *totalJESUnc;
	JME::JetResolution jet_resolution_in_pt;
	JME::JetResolutionScaleFactor jet_resolution_sf_per_eta;

	double tau_kino_cuts_pt, tau_kino_cuts_eta;
	double jet_kino_cuts_pt, jet_kino_cuts_eta;
	double btag_threshold;

	edm::EDGetTokenT<bool> BadChCandFilterToken_;
	edm::EDGetTokenT<bool> BadPFMuonFilterToken_;

	//lumiUtils::GoodLumiFilter goodLumiFilter;

	// random numbers for corrections & uncertainties
	TRandom3 *r3;

	// gotta be a new option for definition of the interface
	//#include "ntupleOutput_leps.h"
	// so, I need to
	//  - declare the interface as class atributes
	//  - then make the TTree and connect branches in constructor
	//  - and reset/fill stuff in analyze
	// let's do it first manually for p4-s of leptons
	TTree* NT_output_ttree; 
	/*
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > NT_lep_p4;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* pt_lep_p4; // yep, vectors of complex objects require additional persistent pointers
	*/
	#define NTUPLE_INTERFACE_CLASS_DECLARE
	#include "UserCode/NtuplerAnalyzer/interface/ntupleOutput.h"
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
NtuplerAnalyzer::NtuplerAnalyzer(const edm::ParameterSet& iConfig) :
dtag       (iConfig.getParameter<std::string>("dtag")),
isMC       (iConfig.getParameter<bool>("isMC")),
muHLT_MC1  (iConfig.getParameter<std::string>("muHLT_MC1")),
muHLT_MC2  (iConfig.getParameter<std::string>("muHLT_MC2")),
muHLT_Data1(iConfig.getParameter<std::string>("muHLT_Data1")),
muHLT_Data2(iConfig.getParameter<std::string>("muHLT_Data2")),
elHLT_Data (iConfig.getParameter<std::string>("elHLT_Data")),
elHLT_MC   (iConfig.getParameter<std::string>("elHLT_MC")),
jecDir     (iConfig.getParameter<std::string>("jecDir")),
TjetResolutionFileName     (iConfig.getParameter<std::string>("resolutionFile")),
TjetResolutionSFFileName   (iConfig.getParameter<std::string>("scaleFactorFile")),
tau_kino_cuts_pt    (iConfig.getParameter<double>("tau_kino_cuts_pt")),
tau_kino_cuts_eta   (iConfig.getParameter<double>("tau_kino_cuts_eta")),
jet_kino_cuts_pt    (iConfig.getParameter<double>("jet_kino_cuts_pt")),
jet_kino_cuts_eta   (iConfig.getParameter<double>("jet_kino_cuts_eta")),
btag_threshold   (iConfig.getParameter<double>("btag_threshold"))
//goodLumiFilter(iConfig.getUntrackedParameter<std::vector<edm::LuminosityBlockRange>>("lumisToProcess", std::vector<edm::LuminosityBlockRange>()))

{
	r3 = new TRandom3();

	/* define in constructor via call to consumes (magic thingy) */
	//tracks_    = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
	muons_     = consumes<pat::MuonCollection>    (edm::InputTag("slimmedMuons"));
	electrons_ = consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"));
	taus_ = consumes<pat::TauCollection>(edm::InputTag("slimmedTaus"));
	vtx_ = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));
	rho_ = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
	trigResults_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
	triggerObjects_ = consumes<vector<pat::TriggerObjectStandAlone>>(edm::InputTag("selectedPatTrigger"));
	//fwlite::Handle<double> rhoHandle;
	//rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
	lheEPToken_ = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
	puInfo_  = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));
	puInfo2_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"));
	genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));

	mets1_ = consumes<pat::METCollection>(edm::InputTag("slimmedMETs"));
	mets2_ = consumes<pat::METCollection>(edm::InputTag("slimmedMETsMuEGClean"));
	mets_uncorrected_ = consumes<pat::METCollection>(edm::InputTag("slimmedMETsUncorrected"));

	genJets_ = consumes<vector<reco::GenJet>>(edm::InputTag("slimmedGenJets"));
	jets_    = consumes<pat::JetCollection>(edm::InputTag("slimmedJets"));

	//BadChCandFilterToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"));
	//BadPFMuonFilterToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuonFilter"));
	//BadChCandFilterToken_ = consumes<bool>(edm::InputTag("BadChargedCandidate"));
	/* try one of these strings:
	 * "BadParticleFilter",
	 *  PFCandidates  = cms.InputTag("particleFlow"),   # Collection to test
	 *  muons  = cms.InputTag("muons"),   # Collection to test
	 *  taggingMode   = cms.bool(False),
	 *  filterType  =cms.string("BadChargedCandidate"
	 *  
	 *  and for muons:
		BadPFMuonFilter = cms.EDFilter(
		    "BadParticleFilter",
		    PFCandidates  = cms.InputTag("particleFlow"),   # Collection to test
		    muons  = cms.InputTag("muons"),   # Collection to test 
		    taggingMode   = cms.bool(False),
		    filterType  =cms.string("BadPFMuon"),
		    maxDR         = cms.double(0.001),              # Maximum DR between reco::muon->innerTrack and pfCandidate 
		...
		BadChargedCandidateFilter = cms.EDFilter(
		    "BadParticleFilter",
		    PFCandidates  = cms.InputTag("particleFlow"),   # Collection to test
		    muons  = cms.InputTag("muons"),   # Collection to test
		    taggingMode   = cms.bool(False),
		    filterType  =cms.string("BadChargedCandidate"),
		    maxDR         = cms.double(0.00001),              # Maximum DR between reco::muon->innerTrack and pfCandidate 
		    minPtDiffRel = cms.double(0.00001),               # lower threshold on difference between pt of reco::muon->innerTrack and pfCandidate
		    minMuonTrackRelErr = cms.double(2.0),          # minimum ptError/pt on muon best track
		    innerTrackRelErr   = cms.double(1.0),          # minimum relPtErr on innerTrack
		    minMuonPt     = cms.double(100.0),               # minimum muon pt 
		    segmentCompatibility = cms.double(0.3),        # compatibility between the inner track and the segments in the muon spectrometer
		)
	 *
	 * -- there is no filter = True option!!
	 */
	//BadPFMuonFilterToken_ = consumes<bool>(edm::InputTag("BadPFMuon"));
	BadPFMuonFilterToken_ = consumes<bool>(edm::InputTag("BadPFMuonFilter"));
	BadChCandFilterToken_ = consumes<bool>(edm::InputTag("BadChargedCandidateFilter"));

	// dtag configs
	bool period_BCD = !isMC && (dtag.Contains("2016B") || dtag.Contains("2016C") || dtag.Contains("2016D"));
	bool period_EF  = !isMC && (dtag.Contains("2016E") || dtag.Contains("2016F"));
	bool period_G   = !isMC && (dtag.Contains("2016G"));
	bool period_H   = !isMC && (dtag.Contains("2016H"));


	// jet IDs, corrections, resolutions etc
	jetID = LooseJET; // TODO: move to Conf
	jetPUID = LoosePUJET;

	// JEC, JES, JER
	gSystem->ExpandPathName (jecDir);
	// v1
	// getJetCorrector(TString baseDir, TString pf, bool isMC)
	//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
	
	// in 2016 the corrections for data are per-period:
	// Summer16_23Sep2016BCDV4_DATA_        Summer16_23Sep2016EFV4_DATA_        Summer16_23Sep2016GV4_DATA_        Summer16_23Sep2016HV4_DATA_
	TString jet_corr_files;
	if (isMC)
		jet_corr_files = "/Summer16_23Sep2016V4_MC";
	else if (period_BCD)
		jet_corr_files = "/Summer16_23Sep2016BCDV4_DATA";
	else if (period_EF)
		jet_corr_files = "/Summer16_23Sep2016EFV4_DATA";
	else if (period_G)
		jet_corr_files = "/Summer16_23Sep2016GV4_DATA";
	else if (period_H)
		jet_corr_files = "/Summer16_23Sep2016HV4_DATA";
	jesCor = utils::cmssw::getJetCorrector (jecDir, jet_corr_files, isMC);
	totalJESUnc = new JetCorrectionUncertainty ((jecDir + jet_corr_files + "_Uncertainty_AK4PFchs.txt").Data());

	// resolution and scale-factors for the systematics
	gSystem->ExpandPathName(TjetResolutionFileName);
	gSystem->ExpandPathName(TjetResolutionSFFileName);

	string jetResolutionFileName   (TjetResolutionFileName);
	string jetResolutionSFFileName (TjetResolutionSFFileName);
	// <---- ROOT & CMSSW are best friends
	jet_resolution_in_pt = JME::JetResolution(jetResolutionFileName);
	jet_resolution_sf_per_eta = JME::JetResolutionScaleFactor(jetResolutionSFFileName);

	edm::Service<TFileService> fs;
	// init ttree
	NT_output_ttree = fs->make<TTree>("reduced_ttree", "TTree with reduced event data");
	// connect the branch
	/*
	// set the additional pointer:
	pt_lep_p4 = &NT_lep_p4;
	NT_output_ttree->Branch("lep_p4", "vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &pt_lep_p4);
	// should be ok
	*/

	// connect the branch with macro:
	#undef NTUPLE_INTERFACE_CLASS_DECLARE
	#define NTUPLE_INTERFACE_CLASS_INITIALIZE
	#include "UserCode/NtuplerAnalyzer/interface/ntupleOutput.h"

	//now do what ever initialization is needed
	usesResource("TFileService");

	/* output via EDM stuff
	 * does it work ok with skipping events?
	 * test the straight forward way, if doesn't work -- leave it, save in TTree manualy
	 *
	 * error: 'produces' was not declared in this scope
	 * -- so it's really just for producers
	 */
	//produces<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >>("leps_p4");

}


NtuplerAnalyzer::~NtuplerAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
NtuplerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	LogInfo ("Demo") << "entered event";

	// reset the output objects with macro
	#undef NTUPLE_INTERFACE_CLASS_DECLARE
	#undef NTUPLE_INTERFACE_CLASS_INITIALIZE
	#define NTUPLE_INTERFACE_CLASS_RESET
	#include "UserCode/NtuplerAnalyzer/interface/ntupleOutput.h"
	// if output contains stand-alone objects (not vector of LorentxVector-s, but just 1 LorentzVector, like MET or something)
	// you have to reset them yourself, since each object might have its' own method
	NT_met_init.SetXYZT(0,0,0,0);
	NT_met_uncorrected.SetXYZT(0,0,0,0);
	NT_met_corrected.SetXYZT(0,0,0,0);

	//Handle<reco::TrackCollection> tracks;
	//iEvent.getByToken( tracks_, tracks );
	//LogInfo("Demo") << "number of tracks "<<tracks->size();
	/*
	if( minTracks_ <= tracks->size() ) {
	   LogInfo("Demo") << "number of tracks "<<tracks->size();
	}
	*/

	// ------------------------------------------------- Apply MET FILTERS

	edm::Handle<bool> ifilterbadChCand;
	edm::Handle<bool> ifilterbadPFMuon;

	iEvent.getByToken(BadChCandFilterToken_, ifilterbadChCand);
	iEvent.getByToken(BadPFMuonFilterToken_, ifilterbadPFMuon);
	//bool  filterbadChCandidate = *ifilterbadChCand;
	//bool filterbadPFMuon = *ifilterbadPFMuon;
	NT_METfilterbadChCand = *ifilterbadChCand;
	NT_METfilterbadPFMuon = *ifilterbadPFMuon;

	// in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2017
	// they say the bool is false if rejected by event

	//if (!(filterbadChCandidate && filterbadPFMuon)) return;
	//if ((filterbadChCandidate && filterbadPFMuon)) return;

	/*
	 * MET filters are data-only thing -- remove events before passing and counting lumi, since MC is then normalized to data lumi
	 * thus after passing lumi data and MC should only have the same cuts
	 *
	 * info on MET filters and their presence in MINIAOD:
	 *   https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#ETmiss_filters
	 *   https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2017
	 *     filter                      location                     data	MC(fullSim)	MC(fastSim)	comment
	 *     primary vertex filter       available in miniAOD v2	DONE	suggested	suggested	 
	 *     beam halo filter            available in miniAOD v2	DONE	suggested	not suggested	Beam Halo Presentation
	 *     HBHE noise filter	   available in miniAOD v2	DONE	suggested	suggested	HCAL DPG Presentation
	 *     HBHEiso noise filter	   available in miniAOD v2	DONE	suggested	suggested	same as above
	 *     ECAL TP filter              available in miniAOD v2	DONE	suggested	suggested	ECAL DPG Presentation
	 *     Bad PF Muon Filter          to be run on the fly 	DONE	suggested	suggested	PPD presentation
	 *     Bad Charged Hadron Filter   to be run on the fly 	DONE	suggested	suggested	PPD presentation
	 *     ee badSC noise filter       available in miniAOD v2	DONE	not suggested	not suggested
	 *   https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
	 *   https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes#MET_Recipes
	 *
	 *   https://twiki.cern.ch/twiki/bin/view/CMS/MissingET
	 *   and their hypernews:
	 *   https://hypernews.cern.ch/HyperNews/CMS/get/met.html
	 */

	/*
	edm::TriggerResultsByName metFilters = iEvent.triggerResultsByName("PAT");   //is present only if PAT (and miniAOD) is not run simultaniously with RECO
	if(!metFilters.isValid()){metFilters = iEvent.triggerResultsByName("RECO");} //if not present, then it's part of RECO
	if(!metFilters.isValid()){       
		LogInfo("Demo") << "TriggerResultsByName for MET filters is not found in the process, as a consequence the MET filter is disabled for this event";
		return;
		}
	if (! isMC && metFilters.isValid())
		{
		// event is good if all filters ar true
		bool filters1 = utils::passTriggerPatterns(metFilters, "Flag_HBHENoiseFilter*", "Flag_HBHENoiseIsoFilter*", "Flag_EcalDeadCellTriggerPrimitiveFilter*");
		bool good_vertices = utils::passTriggerPatterns(metFilters, "Flag_goodVertices");
		bool eebad = utils::passTriggerPatterns(metFilters, "Flag_eeBadScFilter");
		bool halo  = utils::passTriggerPatterns(metFilters, "Flag_globalTightHalo2016Filter");
		// 2016 thing: bad muons
		bool flag_noBadMuons = utils::passTriggerPatterns(metFilters, "Flag_noBadMuons");
		//bool flag_duplicateMuons = utils::passTriggerPatterns(metFilters, "Flag_duplicateMuons");
		// from
		// https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes#Event_flags
		// Three flags are saved in the event:
		//    Flag_badMuons: the event contained at least one PF muon of pT > 20 GeV that is flagged as bad
		//    Flag_duplicateMuons: the event contained at least one PF muon of pT > 20 GeV that is flagged as duplicate
		//    Flag_noBadMuons: the event does not contain any PF muon of pT > 20 GeV flagged as bad or duplicate (i.e. the event is safe)

		// --- thus the Flag_noBadMuons should be enough
		///

		if (! (filters1 & good_vertices & eebad & halo & flag_noBadMuons)) return;
		// these Flag_noBadMuons/Flag_duplicateMuons are MET flags (the issue with bad muons in 2016),
		// they are true if the MET got corrected and event is fine

		// 
		// add: BadChHadron and BadPFMuon -- it seems their name should be Flag_BadChHadron etc
		//
		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2017
		// Bad PF Muon Filter	to be run on the fly	
		// -- "run on the fly", no this flag in the data itself
		//
		// but at the same time:
		//
		// Note that with the in the re-miniaod you will have (will rerun as pointed out below for) the following flags for the "bad muon" events:
		//    Bad PF Muon Filter
		//    Bad Charged Hadrons
		//    Flag_badMuons -> New Giovanni's Filter that the MET is corrected for (flag is set to true if the MET got corrected)
		//    Flag_duplicateMuons -> New Giovanni's Filter that the MET is corrected for (flag is set to true if the MET got corrected)

		// aha https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#ETmiss_filters
		// "Note that many of the current recommended filters can be accessed directly from Miniaod using the flag stored in the TriggerResults,
		//  with the exception of Bad Charged Hadron and Bad Muon Filters."
		// --- so, 2 vs 1 that there should be no Flags for these two in MINIAOD
		//  they should be run on the fly
		///


		//
		// MET POG gives some names to their filters instead of givin the name in code
		// apparently the actual name in the code can be found at:
		// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
		//
		// and there is no BadChHandron
		// the closes to their names are:
		// BadChargedCandidateFilter BadPFMuonFilter
		//
		// -- need to print out what actually is in 03Feb ReReco & ask on hypernews.
		//
		//  found these:
		//  root [7] metFilters.triggerNames()
		//  (const std::vector<std::string> &)
		//  { "Flag_duplicateMuons", "Flag_badMuons", "Flag_noBadMuons",
		//    "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter",
		//    "Flag_CSCTightHaloFilter", "Flag_CSCTightHaloTrkMuUnvetoFilter", "Flag_CSCTightHalo2015Filter",
		//    "Flag_globalTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter",
		//    "Flag_HcalStripHaloFilter", "Flag_hcalLaserEventFilter",
		//    "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellBoundaryEnergyFilter",
		//    "Flag_goodVertices",
		//    "Flag_eeBadScFilter",
		//    "Flag_ecalLaserCorrFilter",
		//    "Flag_trkPOGFilters",
		//    "Flag_chargedHadronTrackResolutionFilter",
		//    "Flag_muonBadTrackFilter",
		//    "Flag_trkPOG_manystripclus53X", "Flag_trkPOG_toomanystripclus53X", "Flag_trkPOG_logErrorTooManyClusters",
		//    "Flag_METFilters" }
		///
		//
		// it seems the track of these two filters goes to:
		// https://indico.cern.ch/event/591506/contributions/2387636/attachments/1381281/2099935/2016_12_01_MET_Scanning_Report_PPD.pdf
		// https://twiki.cern.ch/twiki/bin/view/CMS/MissingETScanners#More_info_on_filter_bad_ChargedC
		// and back to
		// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
		// -- Flag_BadChargedCandidateFilter
		// and Flag_BadPFMuonFilter
		}
		*/

	LogInfo ("Demo") << "passed MET filters";

	// PASS LUMI
	// done with some trick in crab/cmsRun python config
	//if (!isMC)
	//	if(!goodLumiFilter.isGoodLumi(iEvent.eventAuxiliary().run(), iEvent.eventAuxiliary().luminosityBlock())) return; 


	edm::Handle<double> rhoHandle;
	iEvent.getByToken(rho_, rhoHandle);
	if(rhoHandle.isValid() ) NT_fixedGridRhoFastjetAll = *rhoHandle;
	//NT_fixedGridRhoFastjetAll = NT_fixedGridRhoFastjetAll;

	Handle<pat::MuonCollection> muons_h;
	iEvent.getByToken( muons_, muons_h );
	pat::MuonCollection muons = *muons_h;
	Handle<pat::ElectronCollection> electrons_h;
	iEvent.getByToken( electrons_, electrons_h );
	pat::ElectronCollection electrons = *electrons_h;
	Handle<pat::TauCollection> taus_h;
	iEvent.getByToken( taus_, taus_h );
	pat::TauCollection taus = *taus_h;

	//LogInfo("Demo") << "number of muons "<< muons.size();


	// HLT TRIGGER
	bool eTrigger = false;
	bool muTrigger1, muTrigger2, muTrigger = false;

	// TriggerNames for TriggerObjects --------------------
	edm::Handle<edm::TriggerResults> trigResults; //our trigger result object
	//edm::InputTag * trigResultsTag; // the tag object, trigResults are extracted from the event via this tag

	string matched_elTriggerName("");
	string matched_muTriggerName1("");
	string matched_muTriggerName2("");
	//edm::TriggerResultsByName tr = ev.triggerResultsByName ("HLT");
	edm::TriggerResultsByName tr = iEvent.triggerResultsByName ("HLT");
	if (!tr.isValid ()){
		cout << "HLT is NOT valid!" << endl;
		return;
		// HLT2 was a quirk of Spring16 MC campaigns (noHLT/reHLT/withHLT thing)
		//tr = ev.triggerResultsByName ("HLT2");
		}
	else
		{
		LogInfo("Demo") << "Trigger HLT is valid";
		//trigResultsTag = new edm::InputTag("TriggerResults","","HLT"); //make sure have correct process on MC
		// pass trigger
		// using this:
		//   bool passTriggerPatternsAndGetName(edm::TriggerResultsByName& tr, std::string& pathName, std::string pattern)
		// -- pathName is the matched part of the trigger name (as I got it)
		//    it is passed to select trigger objects
		eTrigger   = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_elTriggerName, elHLT_MC)   : utils::passTriggerPatternsAndGetName(tr, matched_elTriggerName, elHLT_Data));
		muTrigger1 = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName1, muHLT_MC1) : utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName1, muHLT_Data1));
		muTrigger2 = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName2, muHLT_MC2) : utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName2, muHLT_Data2));
		muTrigger = muTrigger1 || muTrigger2;
		}

	if (!(eTrigger || muTrigger)) return; // orthogonalization is done afterwards

	LogInfo ("Demo") << "passed HLT " << eTrigger << ' ' << muTrigger << '(' << muTrigger1 << ',' << muTrigger2 << ')' << ';' << matched_elTriggerName << ' ' << matched_muTriggerName1 << ',' << matched_muTriggerName2;

	// names for trigger bits
	//edm::EDGetTokenT<edm::TriggerResults> trigResults_ = consumes<edm::TriggerResults>(trigResultsTag);
	//ev.getByLabel(*trigResultsTag, trigResults);
	iEvent.getByToken( trigResults_, trigResults );
	const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);

	//fwlite::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjectsHandle;
	edm::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjectsHandle;
	//triggerObjectsHandle.getByLabel(ev, "selectedPatTrigger");
	iEvent.getByToken(triggerObjects_, triggerObjectsHandle);
	if (!triggerObjectsHandle.isValid())
		{
		LogInfo("Demo") << "!triggerObjectsHandle.isValid()";
		return;
		}
	LogInfo ("Demo") << "got trigger objects";
	vector<pat::TriggerObjectStandAlone> trig_objs = *triggerObjectsHandle;

	// objects of our triggers
	vector<pat::TriggerObjectStandAlone> el_trig_objs;
	vector<pat::TriggerObjectStandAlone> mu_trig_objs, mu_trig_objs2;

	if (eTrigger)
		{
		NT_HLT_el = true;
		Processing_selectHLTobjects(trig_objs, trigNames, el_trig_objs, matched_elTriggerName);
		}
	if (muTrigger)
		{
		NT_HLT_mu = true;
		Processing_selectHLTobjects(trig_objs, trigNames, mu_trig_objs,  matched_muTriggerName1);
		Processing_selectHLTobjects(trig_objs, trigNames, mu_trig_objs2, matched_muTriggerName2);
		// vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
		mu_trig_objs.insert(mu_trig_objs.end(), mu_trig_objs2.begin(), mu_trig_objs2.end());
		}

	LogInfo ("Demo") << "our trigger objects: " << el_trig_objs.size() << ',' << mu_trig_objs.size();

	// PRIMARY VERTEX
	reco::VertexCollection vtx;
	edm::Handle<reco::VertexCollection> vtxHandle;
	iEvent.getByToken(vtx_, vtxHandle);
	if(vtxHandle.isValid() ) vtx = *vtxHandle;
	NT_nvtx = vtx.size();

	reco::Vertex goodPV;                                                                                                                                                                                               
	unsigned int nGoodPV(0);
	for(size_t ivtx=0; ivtx<vtx.size(); ++ivtx)
		{
		//if(utils::isGoodVertex(vtx[ivtx]))
		// directly from rumors
		// some use:
		// * at least 4 degrees of freedome (ndof) (>=4) (!)
		// * Rho < 2 (impact parameter to the beam spot
		// * z < 24
		bool its_good = (!vtx[ivtx].isFake()) && vtx[ivtx].ndof() > 4 && abs(vtx[ivtx].z()) < 24 && abs(vtx[ivtx].position().Rho()) < 2;
		// it should be equivalent to patUtils procedure
		// only they use reverse: ! > 24 etc -- but without the >=, thus there is 1 bit of discrepancy
		if (its_good)
			{
			if(nGoodPV==0) goodPV=vtx[ivtx];
			nGoodPV++;
			}
		}

	// couple things for MC:
	//  - gen nvtx
	//  - gen NUP (NUmber of Particles? needed for WNJets)
	//  - top pt-s, channel ID and other processing gen particles (save LorentzVector of generated taus, or non-neutrino part of generated tau)
	if(isMC)
		{
		LogInfo ("Demo") << "Processing MC";
		// ----------------------- gen nvtx
		int ngenITpu = 0;
		edm::Handle < std::vector<PileupSummaryInfo>> puInfoH;
		//puInfoH.getByLabel (ev, "slimmedAddPileupInfo");
		iEvent.getByToken(puInfo_, puInfoH);
		if (!puInfoH.isValid())
			{
			//puInfoH.getByLabel( ev, "addPileupInfo" );
			iEvent.getByToken(puInfo2_, puInfoH);
			if (!puInfoH.isValid()) {printf("collection PileupSummaryInfo with name slimmedAddPileupInfo or addPileupInfo does not exist\n"); exit(0);}
			}
		// so here we have valid puInfoH
		// otherwise exit was called
		for (std::vector < PileupSummaryInfo >::const_iterator it = puInfoH->begin (); it != puInfoH->end (); it++)
			{
			//if (it->getBunchCrossing () == 0) ngenITpu += it->getPU_NumInteractions ();
			// guys and Mara use getTrueNumInteractions :
			if (it->getBunchCrossing () == 0) ngenITpu += it->getTrueNumInteractions();
			}
		NT_nvtx_gen = ngenITpu;

		// ----------------------- gen NUP
		edm::Handle < LHEEventProduct > lheEPHandle;
		//lheEPHandle.getByLabel (ev, "externalLHEProducer");
		iEvent.getByToken(lheEPToken_, lheEPHandle);
		if (isMC && lheEPHandle.isValid()) NT_NUP_gen = lheEPHandle->hepeup().NUP;

		LogInfo ("Demo") << "Processing MC, gen particles";
		// ----------------------- gen particles
		// parse gen particles tree and get top pt-s and channel
		// channels are needed for:
		// TTbar, Single-top, DY
		// for now implement thed only TTbar and Single-Top
		reco::GenParticleCollection gen;
		edm::Handle<reco::GenParticleCollection> genHandle;
		//genHandle.getByLabel(ev, "prunedGenParticles");
		iEvent.getByToken(genParticle_, genHandle);
		if(genHandle.isValid())
			{
			gen = *genHandle;
			// For reference, some PDG IDs:
			// QUARKS
			// d  1
			// u  2
			// s  3
			// c  4
			// b  5
			// t  6
			// b' 7
			// t' 8
			// g 21
			// gamma 22
			// Z     23
			// W     24
			// h     25
			// e, ve     11, 12
			// mu, vmu   13, 14
			// tau, vtau 15, 16

			LogInfo ("Demo") << "Processing MC, gen particles, t decays and taus";
			NT_gen_pythia8_prompt_leptons_N = 0;
			NT_gen_N_wdecays = 0;
			NT_gen_N_zdecays = 0;
			for(size_t i = 0; i < gen.size(); ++ i)	
				{
				const reco::GenParticle & p = gen[i];
				int id = p.pdgId();
				int st = p.status(); // TODO: check what is status in decat simulation (pythia for our TTbar set)
				int n_daughters = p.numberOfDaughters();

				if (abs(id) == 6) // if it is a t quark
					{
					if (n_daughters == 2) // it is a decay vertes of t to something
						{
						// find the W decay channel in this top
						unsigned int d0_id = abs(p.daughter(0)->pdgId());
						unsigned int d1_id = abs(p.daughter(1)->pdgId());
						int W_num = d0_id == 24 ? 0 : (d1_id == 24 ? 1 : -1) ;
						if (W_num < 0) continue;
						const reco::Candidate * W = p.daughter( W_num );
						const reco::Candidate * W_final = find_W_decay(W);
						int decay_id = 1;
						// = id of lepton or 1 for quarks
						if (fabs(W_final->daughter(0)->pdgId()) == 11 || fabs(W_final->daughter(0)->pdgId()) == 13 || fabs(W_final->daughter(0)->pdgId()) == 15)
							{
							decay_id = W_final->daughter(0)->pdgId();
							if (fabs(W_final->daughter(0)->pdgId()) == 15)
								decay_id *= simple_tau_decay_id(W_final->daughter(0)); // = 11, 13 for leptons and 20 + 5*(Nch-1) + Npi0 for hadrons
							}
						else if (fabs(W_final->daughter(1)->pdgId()) == 11 || fabs(W_final->daughter(1)->pdgId()) == 13 || fabs(W_final->daughter(1)->pdgId()) == 15)
							{
							decay_id = W_final->daughter(1)->pdgId();
							if (fabs(W_final->daughter(1)->pdgId()) == 15)
								decay_id *= simple_tau_decay_id(W_final->daughter(1));
							}

						// save stuff, according to top Id, also save top p_T
						if (id>0)
							{
							NT_gen_t_pt  = p.pt();
							NT_gen_t_w_decay_id = decay_id;
							}
						else
							{
							NT_gen_tb_pt = p.pt();
							NT_gen_tb_w_decay_id = decay_id;
							}
						}
					}

				// if it is a tau -- save the non-neutrino part to the output
				//  the status is 1 or 2
				//  1. final state, not decays, so it should never happen for tau
				//  2. decayed or fragmented -- the case for tau
				if (abs(id) == 15 && st == 1)
					NT_gen_tau_p4.push_back(p.p4()); 
				else if (abs(id) == 15 && st == 2)
					{
					// it's a final state tau
					// select its' daughters, skipping neutrinos
					// add their momenta -- use the sum as a visible_gen_tau
					LorentzVector vis_ds(0,0,0,0);
					LogInfo ("Demo") << "N tau daughters: " << n_daughters;
					for (int j = 0; j < n_daughters; ++j)
						{
						const reco::Candidate * d = p.daughter(j);
						unsigned int d_id = abs(d->pdgId());
						LogInfo ("Demo") << j << " tau daughter ID = " << d->pdgId();
						if (d_id == 12 || d_id == 14 || d_id == 16) continue;
						vis_ds += d->p4();
						}
					NT_gen_tau_p4.push_back(vis_ds); 
					}

				// Prompt leptons ID for Z->LL
				// pythia 8 stores prompt particles with status 21-29 ("hardest subprocess", PYTHIA 8 Worksheet for tutorial at ASP 2012 Summer School)
				// -- checked DYJetsToLL -- there is no Z (pdgId 23) particles, but prompt leptons work fine
				// thus save N prompt leptons in the process
				// and their ID
				if ((abs(id) == 11 || abs(id) == 13 || abs(id) == 15) && st > 20 && st < 30)
					{
					NT_gen_pythia8_prompt_leptons_N += 1;
					int gen_prompt_lepton_ID = id;
					if (abs(id) == 15)
						gen_prompt_lepton_ID *= simple_tau_decay_id(&p);
					NT_gen_pythia8_prompt_leptons_IDs.push_back(gen_prompt_lepton_ID);
					LogInfo ("Demo") << "Found (pythia8) prompt lepton: " << id << ' ' << gen_prompt_lepton_ID << ' ' << NT_gen_pythia8_prompt_leptons_N;
					}

				// but madgraph DY (50-Inf, i.e. the main one) has the Z-s............
				if (abs(id) == 23 && n_daughters == 2)
					{
					int Zdecay_id = 1;
					int d0_id = abs(p.daughter(0)->pdgId());
					int d1_id = abs(p.daughter(1)->pdgId());
					int lep_daughter = (d0_id == 11 || d0_id == 13 || d0_id == 15 ? 0 : (d1_id == 11 || d1_id == 13 || d1_id == 15 ? 1 : -1));
					if (lep_daughter >= 0)
						{
						Zdecay_id = p.daughter(lep_daughter)->pdgId();
						if (abs(Zdecay_id) == 15)
							Zdecay_id *= simple_tau_decay_id(p.daughter(lep_daughter));
						}

					NT_gen_N_zdecays += 1;
					NT_gen_zdecays_IDs.push_back(Zdecay_id);
					}

				// and W->Lnu processes, apparently prompt leptons don't work there -- sometimes lepton goes directly to final state
				// search for first W decay -- it's supposedly prompt
				// could merge this with t decay procedure, or search from the final state leptons..
				if (abs(id) == 24 && n_daughters == 2)
					{
					int wdecay_id = 1;
					int d0_id = abs(p.daughter(0)->pdgId());
					int d1_id = abs(p.daughter(1)->pdgId());
					int lep_daughter = (d0_id == 11 || d0_id == 13 || d0_id == 15 ? 0 : (d1_id == 11 || d1_id == 13 || d1_id == 15 ? 1 : -1));
					if (lep_daughter >= 0)
						{
						wdecay_id = p.daughter(lep_daughter)->pdgId();
						if (abs(wdecay_id) == 15)
							wdecay_id *= simple_tau_decay_id(p.daughter(lep_daughter));
						}

					NT_gen_N_wdecays += 1;
					NT_gen_wdecays_IDs.push_back(wdecay_id);
					}
				}
			LogInfo ("Demo") << "Found: t decay = " << NT_gen_t_w_decay_id << " ; tb decay = " << NT_gen_tb_w_decay_id;
			}
		}

	double weight = 1;

	// MUONS
	LorentzVector muDiff(0., 0., 0., 0.);
	unsigned int nVetoMu(0);
	pat::MuonCollection selIDMuons, selMuons;
	processMuons_ID_ISO_Kinematics(muons, goodPV, weight, patUtils::llvvMuonId::StdTight, patUtils::llvvMuonId::StdLoose, patUtils::llvvMuonIso::Tight, patUtils::llvvMuonIso::Loose,               
		30., 2.4, 10., 2.5, selIDMuons, muDiff, nVetoMu, false, false);
	nVetoMu += processMuons_MatchHLT(selIDMuons, mu_trig_objs, 0.4, selMuons);

	// ELECTRONS
	pat::ElectronCollection selIDElectrons, selElectrons;
	unsigned int nVetoE(0);
	LorentzVector elDiff(0., 0., 0., 0.);
	processElectrons_ID_ISO_Kinematics(electrons, goodPV, NT_fixedGridRhoFastjetAll, weight, patUtils::llvvElecId::Tight, patUtils::llvvElecId::Loose, patUtils::llvvElecIso::Tight, patUtils::llvvElecIso::Loose,
		30., 2.4, 15., 2.5, selIDElectrons, elDiff, nVetoE, false, false);

	nVetoE += processElectrons_MatchHLT(selIDElectrons, el_trig_objs, 0.4, selElectrons);

	std::vector<patUtils::GenericLepton> selLeptons;
	for(size_t l=0; l<selElectrons.size(); ++l) selLeptons.push_back(patUtils::GenericLepton (selElectrons[l] ));
	for(size_t l=0; l<selMuons.size(); ++l)     selLeptons.push_back(patUtils::GenericLepton (selMuons[l]     ));
	std::sort(selLeptons.begin(), selLeptons.end(), utils::sort_CandidatesByPt);

	LogInfo ("Demo") << "selected leptons: " << '(' << selIDElectrons.size() << ',' << selIDMuons.size() << ')' <<  selLeptons.size() << ' ' << nVetoE << ',' << nVetoMu;

	bool clean_lep_conditions = nVetoE==0 && nVetoMu==0 && nGoodPV != 0;
	if (!(clean_lep_conditions && selLeptons.size() > 0 && selLeptons.size() < 3)) return; // exit now to reduce computation

	LogInfo ("Demo") << "passed lepton conditions ";

	NT_nleps = selLeptons.size();

	NT_leps_ID = 1;
	for (unsigned int i = 0; i<selLeptons.size(); i++)
		{
		NT_leps_ID *= selLeptons[i].pdgId();
		}
	//NT_leps_ID = NT_leps_ID;



	// MET
	pat::METCollection mets;
	edm::Handle<pat::METCollection> metsHandle;
	if (isMC)
		//metsHandle.getByLabel(ev, "slimmedMETs"); // 2016: slimmedMETs are METs corrected by muons
		iEvent.getByToken(mets1_, metsHandle);
	else // ReReco 03Feb data
		//metsHandle.getByLabel(ev, "slimmedMETsMuEGClean");
		iEvent.getByToken(mets2_, metsHandle);
	// 2016: slimmedMETsMuEGClean are corrected by muons and electrons, only in Data!
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes
	if(metsHandle.isValid() ) mets = *metsHandle;
	pat::MET MET = mets[0];
	// LorentzVector met = mets[0].p4 ();

	NT_met_init = MET.p4();

	// also for control let's get uncorrected met and compare the two:
	if (!isMC) // sadly this exists only in latest ReReco data made with 8.0.26 CMSSW, not in Summer16 MC
		{
		pat::METCollection mets_uncorrected;
		edm::Handle<pat::METCollection> mets_uncorrectedHandle;
		//mets_uncorrectedHandle.getByLabel(ev, "slimmedMETsUncorrected");
		iEvent.getByToken( mets_uncorrected_, mets_uncorrectedHandle);
		if(mets_uncorrectedHandle.isValid() ) mets_uncorrected = *mets_uncorrectedHandle;
		//pat::MET met_uncorrected = mets_uncorrected[0];
		//NT_met_uncorrected = met_uncorrected.p4();
		}


	// JETS
	/* jets additionally need initialization of:
	 * genJets
	 * jesCor, totalJESUnc,
	 * pass jet ID, PU jet ID (with/without PU),
	 * systematic variation (NOMINAL, the variation factors are saved per jet for offline)
	 * jet resolution in pt, eta (?)
	 * kinematic cuts
	 */

	// jets
	pat::JetCollection jets;
	edm::Handle<pat::JetCollection>jetsHandle;
	//jetsHandle.getByLabel(ev, "slimmedJets");
	iEvent.getByToken(jets_, jetsHandle);
	if(jetsHandle.isValid() ) jets = *jetsHandle;

	// get genJets from the event
	std::vector<reco::GenJet> genJets;
	edm::Handle<std::vector<reco::GenJet>> genJetsHandle;
	//genJetsHandle.getByLabel(ev, "slimmedGenJets");
	iEvent.getByToken( genJets_, genJetsHandle); // twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#GenJets
	if(genJetsHandle.isValid() ) genJets = *genJetsHandle;

	LorentzVector full_jet_corr(0., 0., 0., 0.);
	pat::JetCollection IDjets;
	//map<systematic_shift, pat::JetCollection> IDjets;
	// it's filled with jetSystematics by processJets_CorrectJES_SmearJERnJES_ID_ISO_with_systematics
	//string jetID("Loose");
	//string jetPUID("MediumPU");
	Variation jet_m_systematic_variation = Variation::NOMINAL;

	processJets_CorrectJES_SmearJERnJES_ID_ISO(jets, genJets, isMC, weight, NT_fixedGridRhoFastjetAll, nGoodPV, jesCor, totalJESUnc, 0.4/2,
		jet_resolution_in_pt, jet_resolution_sf_per_eta, jet_m_systematic_variation, jetID, jetPUID, /*with_PU*/ false, r3, full_jet_corr, IDjets, false, false);

	// ALSO MET
	//LorentzVector MET_corrected = MET.p4() - full_jet_corr;
	//float met_corrected = MET_corrected.pt();
	//NT_met_corrected = MET_corrected;

	pat::JetCollection selJets;
	processJets_Kinematics(IDjets, /*bool isMC,*/ weight, jet_kino_cuts_pt, jet_kino_cuts_eta, selJets, false, false);

	pat::JetCollection selJetsNoLep;
	crossClean_in_dR(selJets, selLeptons, 0.4, selJetsNoLep, weight, string("selJetsNoLep"), false, false);
	// and these are output jets for NTuple
	// they pass ID, corrected with JEC (smeared JES for MC)
	// pass kinematic cuts (pt, eta)
	// and dR-cleaned from selected leptons

	std::sort (selJetsNoLep.begin(),  selJetsNoLep.end(),  utils::sort_CandidatesByPt);

	/* What is saved for b-tagging?
	 * only N of our jets passing Medium WP
	 *
	 * in b-tagging MC is corrected with weight according to efficiency (measured in MC) of tagging and scale-factors (bCallibrator),
	 * the systematics = up/down shifts in calibrator
	 * the efficiencies and SF from the callibrator are done per pt-eta and hadronFlavour
	 * thus all of it can be done offline, on ntuples
	 *
	 * here lets just save nbjets (N of our jets passing medium CSV WP)
	 * for the record decision
	 */
	string btagger_label("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	NT_nbjets = 0;
	for (unsigned int i = 0; i<selJetsNoLep.size(); i++)
		{
		pat::Jet& jet = selJetsNoLep[i];
		float b_discriminator = jet.bDiscriminator(btagger_label);

		// dR match to ID jet -- it's the jet from MiniAOD, without redone corrections (needed for control of corrections)
		LorentzVector id_jet_p4(0,0,0,0);
		for (unsigned int j=0; j<IDjets.size(); j++)
			{
			if (reco::deltaR(jet, IDjets[j]) < 0.4)
				{
				id_jet_p4 = IDjets[j].p4();
				break;
				}
			}
		// find dR match to gen_jet
		LorentzVector gen_jet_p4(0,0,0,0);
		for (unsigned int j=0; j<genJets.size(); j++)
			{
			if (reco::deltaR(jet, genJets[j]) < 0.4)
				{
				gen_jet_p4 = genJets[j].p4();
				break;
				}
			}

		NT_jet_id.push_back(jet.pdgId());
		NT_jet_initial_p4.        push_back(id_jet_p4);
		NT_jet_p4.                push_back(jet.p4());
		NT_jet_uncorrected_p4.    push_back(jet.correctedP4("Uncorrected"));
		NT_jet_matched_genjet_p4. push_back(gen_jet_p4);
		NT_jet_jes_correction.         push_back(jet.userFloat("jes_correction"));
		NT_jet_jes_correction_relShift.push_back(jet.userFloat("jes_correction_relShift"));
		NT_jet_resolution.push_back(jet.userFloat("jet_resolution"));
		NT_jet_sf.       push_back(jet.userFloat("jer_sf"));
		NT_jet_sf_up.    push_back(jet.userFloat("jer_sf_up"));
		NT_jet_sf_down.  push_back(jet.userFloat("jer_sf_down"));
		NT_jet_jer_factor.      push_back(jet.userFloat("jer_factor"));
		NT_jet_jer_factor_up.   push_back(jet.userFloat("jer_factor_up"));
		NT_jet_jer_factor_down. push_back(jet.userFloat("jer_factor_down"));
		//NT_jet_rad.push_back(jet_radius(jet));
		NT_jet_etaetaMoment.push_back(jet.etaetaMoment());
		NT_jet_phiphiMoment.push_back(jet.phiphiMoment());
		NT_jet_pu_discr.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
		NT_jet_b_discr.push_back(b_discriminator);
		NT_jet_hadronFlavour.push_back(jet.hadronFlavour());
		NT_jet_partonFlavour.push_back(jet.partonFlavour());

		if (b_discriminator > btag_threshold) NT_nbjets += 1;
		}


	/*
	 * TAUS
	 */
	//LogInfo("Demo") << "taus.size() = "<< taus.size();
	//string tau_Loose_ID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
	string tau_Loose_ID  ("byLooseIsolationMVArun2v1DBoldDMwLT");
	string tau_Medium_ID ("byMediumIsolationMVArun2v1DBoldDMwLT");
	string tau_Tight_ID  ("byTightIsolationMVArun2v1DBoldDMwLT");
	string tau_decayMode       ("decayModeFinding");
	string tau_againstMuon     ("againstMuonTight3");
	string tau_againstElectron ("againstElectronTightMVA6");

	pat::TauCollection IDtaus, selTaus;
	processTaus_ID_ISO    (taus,   weight, tau_decayMode, tau_Loose_ID, tau_againstMuon, tau_againstElectron, IDtaus, false, false);
	processTaus_Kinematics(IDtaus, weight, tau_kino_cuts_pt, tau_kino_cuts_eta, selTaus,      false, false);

	pat::TauCollection selTausNoLep;
	crossClean_in_dR(selTaus,       selLeptons, 0.4, selTausNoLep,        weight, string("selTausNoLep"),        false, false);

	// and these are the NT output taus
	std::sort (selTausNoLep.begin(),  selTausNoLep.end(),  utils::sort_CandidatesByPt);

	for(size_t i=0; i<selTausNoLep.size(); ++i)
		{
		pat::Tau& tau = selTausNoLep[i];

		Int_t IDlev = 0;
		if (tau.tauID(tau_Tight_ID)) IDlev = 3;
		else if (tau.tauID(tau_Medium_ID)) IDlev = 2;
		else if (tau.tauID(tau_Loose_ID)) IDlev = 1;

		NT_tau_id.push_back(tau.pdgId());
		NT_tau_decayMode.push_back(tau.decayMode());
		NT_tau_p4.push_back(tau.p4());
		NT_tau_IDlev.push_back(IDlev);
		NT_tau_leading_track_pt.push_back(tau.userFloat("leading_track_pt"));
		NT_tau_leadChargedHadrCand_pt.push_back(tau.userFloat("leadChargedHadrCand_pt"));
		NT_tau_leadNeutralCand_pt.push_back(tau.userFloat("leadNeutralCand_pt"));
		NT_tau_leadCand_pt.push_back(tau.userFloat("leadCand_pt"));
		NT_tau_hasSecondaryVertex.push_back(tau.hasSecondaryVertex());
		//NT_tau_hcalEnergy = tau.hcalEnergy();
		//NT_tau_hcalEnergyLeadChargedHadrCand = tau.hcalEnergyLeadChargedHadrCand();

		if (tau.hasSecondaryVertex())
			{
			//const pat::tau::TauPFEssential::CovMatrix& flightCovMatr = taus[i].flightLengthCov();
			float x = tau.flightLength().x();
			float y = tau.flightLength().y();
			float z = tau.flightLength().z();
			NT_tau_flightLength.push_back(x*x + y*y + z*z);
			NT_tau_flightLengthSignificance.push_back(tau.flightLengthSig());
			/*
			LogInfo("Demo") << "flightLengthSig = "<< taus[i].flightLengthSig();
			LogInfo("Demo") << "flightLength    = "<< taus[i].flightLength().x() << ',' << taus[i].flightLength().y() << ',' << taus[i].flightLength().z();
			LogInfo("Demo") << "flightLength Covar = " << endl <<
				flightCovMatr(0,0) << ',' << flightCovMatr(0,1) << ',' << flightCovMatr(0,2) << endl << 
				flightCovMatr(1,0) << ',' << flightCovMatr(1,1) << ',' << flightCovMatr(2,2) << endl << 
				flightCovMatr(2,0) << ',' << flightCovMatr(2,1) << ',' << flightCovMatr(2,2) << endl;
			*/
			}
		else
			{
			NT_tau_flightLength.push_back(-111);
			NT_tau_flightLengthSignificance.push_back(-111);
			}

		// dR match to jet
		Int_t matched_jet_number = -1;
		for (unsigned int i=0; i<selJetsNoLep.size(); i++)
			{
			pat::Jet& jet = selJetsNoLep[i];
			if (reco::deltaR(jet, tau) < 0.4)
				{
				matched_jet_number = i;
				break;
				}
			}

		NT_tau_dR_matched_jet.push_back(matched_jet_number); // number of the jet in jet vectors, if no match = -1
		}

	NT_njets  = selJetsNoLep.size();

	LogInfo ("Demo") << "all particles/objects are selected, nbjets = " << NT_nbjets;

	//bool record_ntuple = (isSingleMu || isSingleE || pass_dileptons) && NT_nbjets > 0 && NT_tau_IDlev.size() > 0; // at least 1 b jet and 1 loose tau
	bool record_ntuple = clean_lep_conditions && selLeptons.size() > 0 && selLeptons.size() < 3 && NT_nbjets > 0; // leptons and at least 1 b jet

	if (record_ntuple)
		{
		LogInfo ("Demo") << "recording";

		for(size_t l=0; l<selLeptons.size(); ++l)
			{
			NT_lep_p4.push_back(selLeptons[l].p4());
			}

		NT_output_ttree->Fill();
		// EDM output
		//   but it's under if! not every event will store stuff -- see if it plays well with rest of EDM system
		//iEvent.put(NT_lep_p4, "leps_p4");
		// some use unknown C++ feature std::move :
		//iEvent.put(std::move(obj), "objname");
		}



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
NtuplerAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NtuplerAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtuplerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtuplerAnalyzer);
