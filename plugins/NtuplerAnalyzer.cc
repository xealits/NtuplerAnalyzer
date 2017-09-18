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

// the vertex fitter for tau secondary vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"//VertexCollection

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

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

double pu_vector_NOMINAL[] = {0, 0.360609416811339, 0.910848525427002, 1.20629960507795, 0.965997726573782, 1.10708082813183, 1.14843491548622, 0.786526251164482, 0.490577792661333, 0.740680941110478,
0.884048630953726, 0.964813189764159, 1.07045369167689, 1.12497267309738, 1.17367530613108, 1.20239808206413, 1.20815108390021, 1.20049333094509, 1.18284686347315, 1.14408796655615,
1.0962284704313, 1.06549162803223, 1.05151011089581, 1.05159666626121, 1.05064452078328, 1.0491726301522, 1.05772537082991, 1.07279673875566, 1.0837536468865, 1.09536667397119,
1.10934472980173, 1.09375894592864, 1.08263679568271, 1.04289345879947, 0.9851490341672, 0.909983816540809, 0.821346330143864, 0.71704523475871, 0.609800913869359, 0.502935245638477,
0.405579825620816, 0.309696044611377, 0.228191137503131, 0.163380359253309, 0.113368437957202, 0.0772279997453792, 0.0508111733313502, 0.0319007262683943, 0.0200879459309245, 0.0122753366005436,
0.00739933885813127, 0.00437426967257811, 0.00260473545284139, 0.00157047254226743, 0.000969500595715493, 0.000733193118123283, 0.000669817107713128, 0.000728548958604492, 0.000934559691182011, 0.00133719688378802,
0.00186652283903214, 0.00314422244976771, 0.00406954793369611, 0.00467888840511915, 0.00505224284441512, 0.00562827194936864, 0.0055889504870752, 0.00522867039470319, 0.00450752163476433, 0.00395300774604375,
0.00330577167682956, 0.00308353042577215, 0.00277846504893301, 0.00223943190687725, 0.00196650068765464, 0.00184742734258922, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0};


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
	edm::EDGetTokenT<edm::TriggerResults> trigResults_, trigResultsRECO_, trigResultsPAT_;
	//triggerObjectsHandle.getByLabel(ev, "selectedPatTrigger");
	//iEvent.getByToken(triggerObjects_, triggerObjectsHandle);
	edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> triggerObjects_;
	edm::EDGetTokenT<LHEEventProduct> lheEPToken_;
	edm::EDGetTokenT< std::vector < PileupSummaryInfo > > puInfo_;
	edm::EDGetTokenT< std::vector < PileupSummaryInfo > > puInfo2_;
	edm::EDGetTokenT<reco::GenParticleCollection> genParticle_;
	edm::EDGetTokenT<GenEventInfoProduct> evt_;

	edm::EDGetTokenT<edm::View<pat::PackedCandidate>> tracks_;
	edm::EDGetTokenT<reco::BeamSpot> beamSpot_;

	edm::EDGetTokenT<pat::METCollection> mets_slimmedMETs_, mets_slimmedMETsMuEGClean_, mets_uncorrected_;

	edm::EDGetTokenT<vector<reco::GenJet>> genJets_;
	//genJetsHandle.getByLabel(ev, "slimmedGenJets");
	edm::EDGetTokenT<pat::JetCollection> jets_;
	//jetsHandle.getByLabel(ev, "slimmedJets");

	bool record_tauID, record_bPreselection, record_MonitorHLT, record_ElMu, record_Dilep;

	TString dtag;
	bool isMC, aMCatNLO;
	bool isLocal;
	string  muHLT_MC1  , muHLT_MC2  ,
		muHLT_Data1, muHLT_Data2,
		elHLT_Data , elHLT_MC,
		lepMonitorHLT;

	jet_id    jetID;
	pu_jet_id jetPUID;

	TString jecDir;
	TString TjetResolutionFileName;
	TString TjetResolutionSFFileName;

	FactorizedJetCorrector *jesCor;
	JetCorrectionUncertainty *totalJESUnc;
	JME::JetResolution jet_resolution_in_pt;
	JME::JetResolutionScaleFactor jet_resolution_sf_per_eta;

	double  el_kino_cuts_pt, el_kino_cuts_eta, el_veto_kino_cuts_pt, el_veto_kino_cuts_eta,
		mu_kino_cuts_pt, mu_kino_cuts_eta, mu_veto_kino_cuts_pt, mu_veto_kino_cuts_eta,
		tau_kino_cuts_pt, tau_kino_cuts_eta;
	double jet_kino_cuts_pt, jet_kino_cuts_eta;
	double btag_threshold;

	edm::EDGetTokenT<bool> BadChCandFilterToken_;
	edm::EDGetTokenT<bool> BadPFMuonFilterToken_;

	lumiUtils::GoodLumiFilter goodLumiFilter;
	std::vector < std::string > urls; // = runProcess.getUntrackedParameter < std::vector < std::string > >("input");
	TString outUrl;

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
	TH1D *event_counter, *weight_counter; 
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
record_tauID         (iConfig.getParameter<bool>("record_tauID"))         ,
record_bPreselection (iConfig.getParameter<bool>("record_bPreselection")) ,
record_MonitorHLT    (iConfig.getParameter<bool>("record_MonitorHLT"))    ,
record_ElMu          (iConfig.getParameter<bool>("record_ElMu"))          ,
record_Dilep         (iConfig.getParameter<bool>("record_Dilep"))         ,
dtag       (iConfig.getParameter<std::string>("dtag")),
isMC       (iConfig.getParameter<bool>("isMC")),
isLocal    (iConfig.getParameter<bool>("isLocal")),
muHLT_MC1  (iConfig.getParameter<std::string>("muHLT_MC1")),
muHLT_MC2  (iConfig.getParameter<std::string>("muHLT_MC2")),
muHLT_Data1(iConfig.getParameter<std::string>("muHLT_Data1")),
muHLT_Data2(iConfig.getParameter<std::string>("muHLT_Data2")),
elHLT_Data (iConfig.getParameter<std::string>("elHLT_Data")),
elHLT_MC   (iConfig.getParameter<std::string>("elHLT_MC")),
lepMonitorHLT   (iConfig.getParameter<std::string>("lepMonitorHLT")),
jecDir     (iConfig.getParameter<std::string>("jecDir")),
TjetResolutionFileName     (iConfig.getParameter<std::string>("resolutionFile")),
TjetResolutionSFFileName   (iConfig.getParameter<std::string>("scaleFactorFile")),
el_kino_cuts_pt    (iConfig.getParameter<double>("el_kino_cuts_pt")),
el_kino_cuts_eta   (iConfig.getParameter<double>("el_kino_cuts_eta")),
el_veto_kino_cuts_pt    (iConfig.getParameter<double>("el_veto_kino_cuts_pt")),
el_veto_kino_cuts_eta   (iConfig.getParameter<double>("el_veto_kino_cuts_eta")),
mu_kino_cuts_pt    (iConfig.getParameter<double>("mu_kino_cuts_pt")),
mu_kino_cuts_eta   (iConfig.getParameter<double>("mu_kino_cuts_eta")),
mu_veto_kino_cuts_pt    (iConfig.getParameter<double>("mu_veto_kino_cuts_pt")),
mu_veto_kino_cuts_eta   (iConfig.getParameter<double>("mu_veto_kino_cuts_eta")),
tau_kino_cuts_pt    (iConfig.getParameter<double>("tau_kino_cuts_pt")),
tau_kino_cuts_eta   (iConfig.getParameter<double>("tau_kino_cuts_eta")),
jet_kino_cuts_pt    (iConfig.getParameter<double>("jet_kino_cuts_pt")),
jet_kino_cuts_eta   (iConfig.getParameter<double>("jet_kino_cuts_eta")),
btag_threshold   (iConfig.getParameter<double>("btag_threshold")),
goodLumiFilter   (iConfig.getUntrackedParameter<std::vector<edm::LuminosityBlockRange>>("lumisToProcess", std::vector<edm::LuminosityBlockRange>())),
urls   (iConfig.getUntrackedParameter <std::vector <std::string> >("input")),
outUrl (iConfig.getParameter<std::string>("outfile"))

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
	trigResultsRECO_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","RECO"));
	trigResultsPAT_     = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","PAT"));
	triggerObjects_ = consumes<vector<pat::TriggerObjectStandAlone>>(edm::InputTag("selectedPatTrigger"));
	//fwlite::Handle<double> rhoHandle;
	//rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
	lheEPToken_ = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
	puInfo_  = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));
	puInfo2_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"));
	genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));

	evt_  = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

	tracks_ = consumes<edm::View<pat::PackedCandidate>> (edm::InputTag("packedPFCandidates"));
	beamSpot_ = consumes<reco::BeamSpot> (edm::InputTag("offlineBeamSpot"));

	mets_slimmedMETs_ = consumes<pat::METCollection>(edm::InputTag("slimmedMETs"));
	mets_slimmedMETsMuEGClean_ = consumes<pat::METCollection>(edm::InputTag("slimmedMETsMuEGClean"));
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

	aMCatNLO = dtag.Contains("amcatnlo");


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
	event_counter  = fs->make<TH1D>( "events_counter"  , "pass category", 100,  0, 100);
	weight_counter = fs->make<TH1D>( "weight_counter"  , "pass category", 100,  0, 100); // for control of effect from PU reweighting, aMCatNLO gen -1 weights, top pt
	// in principle it should be orthogonal to the ntuple preselection
	// and be possible to do it after the preselection
	// N_presel / weighted N_presel = N_all / weighted N_all

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
	LogInfo ("Demo") << "entered event " << iEvent.eventAuxiliary().run() << ' ' << iEvent.eventAuxiliary().luminosityBlock();

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
	NT_met_slimmedMets.SetXYZT(0,0,0,0);
	NT_met_slimmedMETsMuEGClean.SetXYZT(0,0,0,0);

	NT_jets_full_correction.SetXYZT(0,0,0,0);

	NT_gen_t_w_final_p4.SetXYZT(0,0,0,0);
	NT_gen_t_b_final_p4.SetXYZT(0,0,0,0);
	NT_gen_tb_w_final_p4.SetXYZT(0,0,0,0);
	NT_gen_tb_b_final_p4.SetXYZT(0,0,0,0);

	math::Error<3>::type pvCov;
	pvCov(0,0) = 999;
	pvCov(1,1) = 999;
	pvCov(2,2) = 999;
	NT_PV_cov = pvCov;


	unsigned int event_checkpoint = 0;
	event_counter->Fill(event_checkpoint++);
	double weight = 1;

	NT_indexevents = iEvent.id().event();
	NT_runNumber   = iEvent.id().run();
	NT_lumi        = iEvent.luminosityBlock();

	// couple things for MC:
	//  - gen aMCatNLO
	//  - gen nvtx
	//  - gen NUP (NUmber of Particles? needed for WNJets)
	//  - top pt-s, channel ID and other processing gen particles (save LorentzVector of generated taus, or non-neutrino part of generated tau)
	if(isMC)
		{
		LogInfo ("Demo") << "Processing MC";
		double weight_TopPT = 1, weight_Gen = 1;
		// ----------------------- aMCatNLO -1 weights
		//fwlite::Handle<GenEventInfoProduct> evt;
		//evt.getByLabel(ev, "generator");
		edm::Handle<GenEventInfoProduct> evt;
		iEvent.getByToken(evt_, evt);
		if(evt.isValid())
			{
			//if (debug) cout << "evt is valid, evt->weight() = " << evt->weight() << "\n";
			LogInfo("Demo") << "evt is valid, evt->weight() = " << evt->weight();
			weight_Gen = (evt->weight() > 0 ) ? 1. : -1. ;
			NT_aMCatNLO_weight = evt->weight();
			}

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

        		vector<const reco::Candidate*> t_b_parts, tb_b_parts, t_W_parts, tb_W_parts;
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

				if (abs(id) == 6 && n_daughters == 2) // if it is a t quark
					{ // it is a decay vertes of t to something
					// calculate top_pt weights:
					weight_TopPT *= TMath::Sqrt(top_pT_SF(p.pt()));
					// find the W decay channel in this top
					unsigned int d0_id = abs(p.daughter(0)->pdgId());
					unsigned int d1_id = abs(p.daughter(1)->pdgId());
					int W_num = d0_id == 24 ? 0 : (d1_id == 24 ? 1 : -1) ;
					if (W_num < 0) continue;
					const reco::Candidate * W = p.daughter( W_num );
					const reco::Candidate * b = p.daughter( 1 - W_num );
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
						save_final_states(W, NT_gen_t_w_final_p4s, NT_gen_t_w_final_pdgIds, NT_gen_t_w_final_statuses, t_W_parts);
						save_final_states(b, NT_gen_t_b_final_p4s, NT_gen_t_b_final_pdgIds, NT_gen_t_b_final_statuses, t_b_parts);
						for (unsigned int i=0; i<NT_gen_t_w_final_p4s.size(); i++)
							{
							// skip neutrinos
							if (abs(NT_gen_t_w_final_pdgIds[i]) == 12 || abs(NT_gen_t_w_final_pdgIds[i]) == 14 || abs(NT_gen_t_w_final_pdgIds[i]) == 16) continue;
							NT_gen_t_w_final_p4 += NT_gen_t_w_final_p4s[i];
							}
						for (unsigned int i=0; i<NT_gen_t_b_final_p4s.size(); i++)
							{
							// skip neutrinos
							if (abs(NT_gen_t_b_final_pdgIds[i]) == 12 || abs(NT_gen_t_b_final_pdgIds[i]) == 14 || abs(NT_gen_t_b_final_pdgIds[i]) == 16) continue;
							NT_gen_t_b_final_p4 += NT_gen_t_b_final_p4s[i];
							}
						}
					else
						{
						NT_gen_tb_pt = p.pt();
						NT_gen_tb_w_decay_id = decay_id;
						save_final_states(W, NT_gen_tb_w_final_p4s, NT_gen_tb_w_final_pdgIds, NT_gen_tb_w_final_statuses, tb_W_parts);
						save_final_states(b, NT_gen_tb_b_final_p4s, NT_gen_tb_b_final_pdgIds, NT_gen_tb_b_final_statuses, tb_b_parts);
						for (unsigned int i=0; i<NT_gen_tb_w_final_p4s.size(); i++)
							{
							// skip neutrinos
							if (abs(NT_gen_tb_w_final_pdgIds[i]) == 12 || abs(NT_gen_tb_w_final_pdgIds[i]) == 14 || abs(NT_gen_tb_w_final_pdgIds[i]) == 16) continue;
							NT_gen_tb_w_final_p4 += NT_gen_tb_w_final_p4s[i];
							}
						for (unsigned int i=0; i<NT_gen_tb_b_final_p4s.size(); i++)
							{
							// skip neutrinos
							if (abs(NT_gen_tb_b_final_pdgIds[i]) == 12 || abs(NT_gen_tb_b_final_pdgIds[i]) == 14 || abs(NT_gen_tb_b_final_pdgIds[i]) == 16) continue;
							NT_gen_tb_b_final_p4 += NT_gen_tb_b_final_p4s[i];
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
					NT_gen_N_zdecays += 1;
					int d0_id = p.daughter(0)->pdgId();
					int d1_id = p.daughter(1)->pdgId();
					int a_d0_id = abs(d0_id);
					int a_d1_id = abs(d1_id);
					int lep_daughter = (a_d0_id == 11 || a_d0_id == 13 || a_d0_id == 15 ? 0 : (a_d1_id == 11 || a_d1_id == 13 || a_d1_id == 15 ? 1 : -1));
					if (lep_daughter >= 0)
						{
						if (a_d0_id == 15) d0_id *= simple_tau_decay_id(p.daughter(0));
						if (a_d1_id == 15) d1_id *= simple_tau_decay_id(p.daughter(1));
						}
					NT_gen_zdecays_IDs.push_back(d0_id);
					NT_gen_zdecays_IDs.push_back(d1_id);
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

		weight = pu_vector_NOMINAL[NT_nvtx_gen] * (aMCatNLO? weight_Gen : 1) * weight_TopPT;
		}
	weight_counter->Fill(event_checkpoint, weight);

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

	edm::TriggerResultsByName metFilters = iEvent.triggerResultsByName("RECO"); //is present only if PAT (and miniAOD) is not run simultaniously with RECO
	if(!metFilters.isValid()){metFilters = iEvent.triggerResultsByName("PAT");} //if not present, then it's part of RECO
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
		//bool flag_noBadMuons = utils::passTriggerPatterns(metFilters, "Flag_noBadMuons"); // <---- the bad muons are done on the fly with cfg.py thingy
		//bool flag_duplicateMuons = utils::passTriggerPatterns(metFilters, "Flag_duplicateMuons");
		// from
		// https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes#Event_flags
		// Three flags are saved in the event:
		//    Flag_badMuons: the event contained at least one PF muon of pT > 20 GeV that is flagged as bad
		//    Flag_duplicateMuons: the event contained at least one PF muon of pT > 20 GeV that is flagged as duplicate
		//    Flag_noBadMuons: the event does not contain any PF muon of pT > 20 GeV flagged as bad or duplicate (i.e. the event is safe)

		// --- thus the Flag_noBadMuons should be enough
		///

		//if (! (filters1 & good_vertices & eebad & halo & flag_noBadMuons)) return;
		// save for study
		NT_pass_basic_METfilters = filters1 & good_vertices & eebad & halo;
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

	LogInfo ("Demo") << "passed MET filters";

	// PASS LUMI
	// done with some trick in crab/cmsRun python config
	if (!isMC && isLocal)
		{
		if (!goodLumiFilter.isGoodLumi(iEvent.eventAuxiliary().run(), iEvent.eventAuxiliary().luminosityBlock())) return; 
		}
	event_counter->Fill(event_checkpoint++);
	weight_counter->Fill(event_checkpoint, weight);


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
	bool lepMonitorTrigger = false;
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
		lepMonitorTrigger   = utils::passTriggerPatterns(tr, lepMonitorHLT);
		eTrigger   = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_elTriggerName, elHLT_MC)   : utils::passTriggerPatternsAndGetName(tr, matched_elTriggerName, elHLT_Data));
		muTrigger1 = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName1, muHLT_MC1) : utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName1, muHLT_Data1));
		muTrigger2 = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName2, muHLT_MC2) : utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName2, muHLT_Data2));
		muTrigger = muTrigger1 || muTrigger2;
		}

	if (!(eTrigger || muTrigger || lepMonitorTrigger)) return; // orthogonalization is done afterwards
	event_counter ->Fill(event_checkpoint++);
	weight_counter->Fill(event_checkpoint, weight);

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

	NT_HLT_lepMonitor = lepMonitorTrigger;
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
			// save info on the vertex
			NT_PV_x.push_back(vtx[ivtx].x());
			NT_PV_y.push_back(vtx[ivtx].y());
			NT_PV_z.push_back(vtx[ivtx].z());
			NT_PV_x_err.push_back(vtx[ivtx].xError());
			NT_PV_y_err.push_back(vtx[ivtx].yError());
			NT_PV_z_err.push_back(vtx[ivtx].zError());
			}
		}

	//weight = 1; // reset weights?

	// MUONS
	LorentzVector muDiff(0., 0., 0., 0.);
	unsigned int nVetoMu(0);
	//pat::MuonCollection selIDMuons, selMuons;
	pat::MuonCollection selMuons;
	processMuons_ID_ISO_Kinematics(muons, goodPV, weight, patUtils::llvvMuonId::StdTight, patUtils::llvvMuonId::StdLoose, patUtils::llvvMuonIso::Tight, patUtils::llvvMuonIso::Loose,               
		mu_kino_cuts_pt, mu_kino_cuts_eta, mu_veto_kino_cuts_pt, mu_veto_kino_cuts_eta, selMuons, muDiff, nVetoMu, false, false);

	//nVetoMu += processMuons_MatchHLT(selIDMuons, mu_trig_objs, 0.4, selMuons);

	// ELECTRONS
	//pat::ElectronCollection selIDElectrons, selElectrons;
	pat::ElectronCollection selElectrons;
	unsigned int nVetoE(0);
	LorentzVector elDiff(0., 0., 0., 0.);
	processElectrons_ID_ISO_Kinematics(electrons, goodPV, NT_fixedGridRhoFastjetAll, weight, patUtils::llvvElecId::Tight, patUtils::llvvElecId::Loose, patUtils::llvvElecIso::Tight, patUtils::llvvElecIso::Loose,
		el_kino_cuts_pt, el_kino_cuts_eta, el_veto_kino_cuts_pt, el_veto_kino_cuts_eta, selElectrons, elDiff, nVetoE, false, false);

	//nVetoE += processElectrons_MatchHLT(selIDElectrons, el_trig_objs, 0.4, selElectrons);

	for(size_t l=0; l<selMuons.size(); ++l)
		{
		NT_lep_p4.push_back(selMuons[l].p4());
		NT_lep_id.push_back(selMuons[l].pdgId());
		// mu_trig_objs or el_trig_objs
		NT_lep_matched_HLT.push_back(processMuon_matchesHLTs(selMuons[l], mu_trig_objs, 0.4));
		NT_lep_dz  .push_back(selMuons[l].muonBestTrack()->dz (goodPV.position()));
		NT_lep_dxy .push_back(selMuons[l].muonBestTrack()->dxy (goodPV.position()));
		NT_lep_dB.push_back(selMuons[l].dB());
		NT_lep_relIso.push_back(relIso(selMuons[l], NT_fixedGridRhoFastjetAll));
		}

	for(size_t l=0; l<selElectrons.size(); ++l)
		{
		NT_lep_p4.push_back(selElectrons[l].p4());
		NT_lep_id.push_back(selElectrons[l].pdgId());
		// mu_trig_objs or el_trig_objs
		NT_lep_matched_HLT.push_back(processElectron_matchesHLTs(selElectrons[l], el_trig_objs, 0.4));
		NT_lep_dz  .push_back(selElectrons[l].gsfTrack()->dz (goodPV.position()));
		NT_lep_dxy .push_back(selElectrons[l].gsfTrack()->dxy (goodPV.position()));
		NT_lep_dB.push_back(selElectrons[l].dB());
		NT_lep_relIso.push_back(relIso(selElectrons[l], NT_fixedGridRhoFastjetAll));
		}

	std::vector<patUtils::GenericLepton> selLeptons;
	for(size_t l=0; l<selElectrons.size(); ++l) selLeptons.push_back(patUtils::GenericLepton (selElectrons[l] ));
	for(size_t l=0; l<selMuons.size(); ++l)     selLeptons.push_back(patUtils::GenericLepton (selMuons[l]     ));
	std::sort(selLeptons.begin(), selLeptons.end(), utils::sort_CandidatesByPt);

	//LogInfo ("Demo") << "selected leptons: " << '(' << selIDElectrons.size() << ',' << selIDMuons.size() << ')' <<  selLeptons.size() << ' ' << nVetoE << ',' << nVetoMu;
	LogInfo ("Demo") << "selected leptons: " << '(' << selElectrons.size() << ',' << selMuons.size() << ')' <<  selLeptons.size() << ' ' << nVetoE << ',' << nVetoMu;

	bool clean_lep_conditions = nVetoE==0 && nVetoMu==0 && nGoodPV != 0;
	//if (!(clean_lep_conditions && selLeptons.size() > 0 && selLeptons.size() < 3)) return; // exit now to reduce computation -- nope, there are other records too
	event_counter ->Fill(event_checkpoint++);
	weight_counter->Fill(event_checkpoint, weight);

	LogInfo ("Demo") << "passed lepton conditions ";

	NT_nleps = selLeptons.size();

	NT_leps_ID = 1;
	for (unsigned int i = 0; i<selLeptons.size(); i++)
		{
		NT_leps_ID *= selLeptons[i].pdgId();
		}
	//NT_leps_ID = NT_leps_ID;



	// MET
	pat::METCollection mets_slimmedMETs, mets_slimmedMETsMuEGClean;
	edm::Handle<pat::METCollection> metsHandle_slimmedMETs;
	edm::Handle<pat::METCollection> metsHandle_slimmedMETsMuEGClean;
	iEvent.getByToken(mets_slimmedMETs_, metsHandle_slimmedMETs);
	iEvent.getByToken(mets_slimmedMETsMuEGClean_, metsHandle_slimmedMETsMuEGClean);

	// 2016: slimmedMETsMuEGClean are corrected by muons and electrons, only in Data!
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes
	// slimmedMets should be there for both Data and MC
	if(metsHandle_slimmedMETs.isValid() )
		{
		const pat::MET& MET = metsHandle_slimmedMETs->front();
		NT_met_slimmedMets = MET.p4();
		NT_met_init = MET.p4();
		}
	// LorentzVector met = mets_slimmedMETs[0].p4 ();

	// these are valid only for data
	if(metsHandle_slimmedMETsMuEGClean.isValid() )
		{
		const pat::MET& MET2 = metsHandle_slimmedMETsMuEGClean->front();
		NT_met_slimmedMETsMuEGClean = MET2.p4();
		}

	// also for control let's get uncorrected met and compare the two:
	if (!isMC) // sadly this exists only in latest ReReco data made with 8.0.26 CMSSW, not in Summer16 MC
		{
		pat::METCollection mets_uncorrected;
		edm::Handle<pat::METCollection> mets_uncorrectedHandle;
		//mets_uncorrectedHandle.getByLabel(ev, "slimmedMETsUncorrected");
		iEvent.getByToken( mets_uncorrected_, mets_uncorrectedHandle);
		if(mets_uncorrectedHandle.isValid() ) mets_uncorrected = *mets_uncorrectedHandle;
		pat::MET met_uncorrected = mets_uncorrected[0];
		NT_met_uncorrected = met_uncorrected.p4();
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
	//pat::JetCollection IDjets;
	//pat::JetCollection selJets;
	pat::JetCollection selJetsNoLep;
	//map<systematic_shift, pat::JetCollection> IDjets;
	// it's filled with jetSystematics by processJets_CorrectJES_SmearJERnJES_ID_ISO_with_systematics
	//string jetID("Loose");
	//string jetPUID("MediumPU");
	//Variation jet_m_systematic_variation = Variation::NOMINAL;

	//processJets_CorrectJES_SmearJERnJES_ID_ISO(jets, genJets, isMC, weight, NT_fixedGridRhoFastjetAll, nGoodPV, jesCor, totalJESUnc, 0.4/2,
	//	jet_resolution_in_pt, jet_resolution_sf_per_eta, jet_m_systematic_variation, jetID, jetPUID, /*with_PU*/ false, r3, full_jet_corr, IDjets, false, false);
	// selecting jets in the following loop
	// similar to https://github.com/LLRCMS/LLRHiggsTauTau/blob/b8bc9146cab462fdaf8e4161761b5e70a08f4a65/NtupleProducer/plugins/HTauTauNtuplizer.cc

	NT_nbjets = 0;
	NT_nallbjets = 0;
	NT_njets = 0;
	NT_nalljets = 0;
	string btagger_label("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	for (unsigned int ijet=0; ijet<jets.size(); ijet++)
		{
		pat::Jet& jet = jets[ijet];
		// the only requirements are pt and no-leptons in dR:
		if (jet.pt() < jet_kino_cuts_pt) continue;
		// all eta pass -- forward jets too, usefull for WJets control region

		//crossClean_in_dR(selJets, selLeptons, 0.4, selJetsNoLep, weight, string("selJetsNoLep"), false, false);
	        bool overlapWithLepton(false);
		float min_dR = 0.4;
	        for(unsigned int l=0; l<(unsigned int)selLeptons.size();++l)
	                {
	                if (reco::deltaR(jet, selLeptons[l])<min_dR)
	                        { overlapWithLepton=true; break; }
	                }
	        if (overlapWithLepton) continue;

		NT_jet_id.push_back(jet.pdgId());
		NT_jet_p4.                push_back(jet.p4());

		LorentzVector jet_initial_p4 = jet.p4();

		//NT_jet_rad.push_back(jet_radius(jet));
		NT_jet_etaetaMoment.push_back(jet.etaetaMoment());
		NT_jet_phiphiMoment.push_back(jet.phiphiMoment());
		NT_jet_pu_discr.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
		//NT_jet_pu_discr_updated.push_back(ijet->hasUserFloat("pileupJetIdUpdated:fullDiscriminant") ? ijet->userFloat("pileupJetIdUpdated:fullDiscriminant") : -999);
		float b_discriminator = jet.bDiscriminator(btagger_label);
		NT_jet_b_discr.push_back(b_discriminator);
		if (b_discriminator > btag_threshold) NT_nallbjets += 1;

		NT_jet_hadronFlavour.push_back(jet.hadronFlavour());
		NT_jet_partonFlavour.push_back(jet.partonFlavour());

		// PF jet ID
	        // from the twiki:
	        // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
	        float NHF = jet.neutralHadronEnergyFraction();
	        float NEMF = jet.neutralEmEnergyFraction();
	        float CHF = jet.chargedHadronEnergyFraction();
	        float MUF = jet.muonEnergyFraction();
	        float CEMF = jet.chargedEmEnergyFraction(); // chargedEmEnergyFraction (relative to uncorrected jet energy)
	        float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	        float NumNeutralParticles =jet.neutralMultiplicity();
	        float CHM = jet.chargedMultiplicity();

		double abseta = fabs(jet.eta());

		bool looseJetID = false;
		bool tightJetID = false;
		bool tightLepVetoJetID = false;

		// and latest (at Moriond17) stuff:
		// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
		// quote:
		// > Note: All fractions are calculated with the raw/uncorrected energy of the jet (only then they add up to unity). So the PF JetID has to be applied before the jet energy corrections.
		// --- khmm....
		// in MINIAOD the jets are already corrected  --> one gets the uncorrected jet and reapplies the corrections
		// > ... collection of AK4 jets slimmedJets, made from ak4PFJetsCHS ... "These have standard jet energy corrections applied (L1FastJet, L2, L3), and a pT cut at 10 GeV"
		// so now one need to get the uncorrected jet --> 
		// TODO <-- it doesn't make sense: corrected jet differs only in pt, which doesn't matter in these cuts

		if (abseta <= 2.7)
			{
			looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abseta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abseta>2.4);
			tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abseta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abseta>2.4);
			tightLepVetoJetID = ((NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abseta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abseta>2.4) );
			}
		else if (abseta <= 3.0)
			{
			looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
			tightJetID = looseJetID;
			}
		else
			{
			looseJetID = (NEMF<0.90 && NumNeutralParticles>10);
			tightJetID = looseJetID;
			}

		int jetid = 0;
		if (looseJetID)
			{
			++jetid;
			if (abseta < jet_kino_cuts_eta)
				{
				NT_njets += 1;
				if (b_discriminator > btag_threshold) NT_nbjets += 1;
				}
			// counter of our old jets: cut on pt, on eta (small eta, central jet) and PFID Loose
			}
		if (tightJetID) ++jetid;
		if (tightLepVetoJetID) ++jetid;

		NT_jet_PFID.push_back(jetid);

		// corrections
		NT_jet_area.push_back (jet.jetArea());

		float jecFactor = jet.jecFactor("Uncorrected") ;
	 	//float jetRawPt = jecFactor * jet.pt();
		//float jetRawPt2 = ijet->pt() / jecFactor; // this is wrong
		NT_jet_uncorrected_jecFactor.push_back (jecFactor);
		LorentzVector rawJet = jet.correctedP4("Uncorrected");
		NT_jet_uncorrected_p4.    push_back(rawJet);
		// TODO: compare the two offline

		// save FactorizedCorrector factor
		//double toRawSF=jet.correctedJet("Uncorrected").pt()/jet.pt();
		//LorentzVector rawJet(jet*toRawSF);
		// FactorizedCorrector with all jet correction files
		jesCor->setJetEta(rawJet.eta());
		jesCor->setJetPt(rawJet.pt());
		jesCor->setJetA(jet.jetArea());
		jesCor->setRho(NT_fixedGridRhoFastjetAll);
		//jesCor->setNPV(nGoodPV); // not used in PU Jet ID example, shouldn't matter
		float jes_correction = jesCor->getCorrection();
		//jet.setP4(rawJet*jes_correction);
		//LorentzVector jesCorJet = (rawJet*jes_correction);
		//jet.addUserFloat("jes_correction", jes_correction);
		// TODO: compare jet_p4 (default MiniAOD jet) and uncorrected_p4 * jes_correction <- re-corrected jet
		NT_jet_jes_recorrection.push_back(jes_correction);

		float dR_max = 0.4/2;
		double jet_resolution = -1;
		double jer_sf = -1;
		double jer_sf_up = -1;
		double jer_sf_down = -1;
		// and the final factors from these SFs
		double jer_factor = -1, jer_factor_up = -1, jer_factor_down = -1;
		LorentzVector gen_jet_p4(0,0,0,0);
		// here is the matching of the jet:
		if(isMC)
			{
			// the JER SF and resolution for the jet:
			//std::vector<double> jer_sf_pair = JER_SF(jet.eta());
			//double jer_sf = jer_sf_pair[0];
			//double jer_resolution = jer_sf_pair[1]; // TODO: not sure about this -- is the table the same as what their tool returns?
			// getting it with the tool from files:
			jet_resolution = jet_resolution_in_pt.getResolution({{JME::Binning::JetPt, jet.pt()},
				{JME::Binning::JetEta, jet.eta()},
				{JME::Binning::Rho, NT_fixedGridRhoFastjetAll}});
			//jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, m_systematic_variation);
			jer_sf      = jet_resolution_sf_per_eta.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::NOMINAL);
			jer_sf_up   = jet_resolution_sf_per_eta.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::UP);
			jer_sf_down = jet_resolution_sf_per_eta.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::DOWN);

			// matching to generation jet:
			//const reco::GenJet* genJet=jet.genJet();
			// the PAT doc describes it as "return the matched generated jet"
			// what's the matching procedure?
			// for now I'll do it manually in dR, as shown in Jet POG example
			// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
			const reco::GenJet* matched_genJet = nullptr;
			//double dR_max = 0.4/2; // 0.4 is the jet cone parameter of AK4 jets, which I use
			// moved it to parameters of the procedure
			for (unsigned int i=0; i<genJets.size(); i++)
				{
				reco::GenJet& genJet = genJets[i];
				double dR = reco::deltaR(jet, genJet);

				if (dR > dR_max) continue;

				double dPt = std::abs(genJet.pt() - jet.pt());
				double dPt_max_factor = 3*jet.pt(); // from twiki
				if (dPt > dPt_max_factor * jet_resolution) continue;

				matched_genJet = &genJet;
				}

			// calculate and apply the smearing factor to jet energy
			// with one of the two algorithms, if jet matched to genJet or not
			if (matched_genJet)
				{ // the scaling is ok
				gen_jet_p4 = matched_genJet->p4();

				double dPt = jet.pt() - matched_genJet->pt();
				//double genjetpt( genJet ? genJet->pt(): 0.);                    
				//std::vector<double> smearJER=utils::cmssw::smearJER(jet.pt(),jet.eta(),genjetpt);
				// using the local smear:
				//std::vector<double> JER_smearing_factor = smearJER(jet.pt(),jet.eta(),genjetpt);
				//double jer_smearing = JER_smearing_factor[0];
				jer_factor = TMath::Max(0., 1. + (jer_sf - 1) * dPt / jet.pt());
				jer_factor_up   = TMath::Max(0., 1. + (jer_sf_up   - 1) * dPt / jet.pt());
				jer_factor_down = TMath::Max(0., 1. + (jer_sf_down - 1) * dPt / jet.pt());
				jet.setP4(jet.p4()*jer_factor); // same as scaleEnergy in the Jet POG example
				// but they also do MIN_ENERGY thing
				// which is static constexpr const double MIN_JET_ENERGY = 1e-2;
				}
			else
				{ // the smearing with gaussian randomizer

				// this is the example:
				//double sigma = jet_resolution * std::sqrt(jer_sf * jer_sf - 1);
				//double smearFactor = 1 + r3->Gaus(0, sigma);
				// this is the twiki:
				double smearFactor      = 1. + r3->Gaus(0, jet_resolution) * std::sqrt(TMath::Max(0., jer_sf*jer_sf - 1.));
				double smearFactor_up   = 1. + r3->Gaus(0, jet_resolution) * std::sqrt(TMath::Max(0., jer_sf_up*jer_sf_up - 1.));
				double smearFactor_down = 1. + r3->Gaus(0, jet_resolution) * std::sqrt(TMath::Max(0., jer_sf_down*jer_sf_down - 1.));
				jer_factor = TMath::Max(0., smearFactor);
				jer_factor_up   = TMath::Max(0., smearFactor_up);
				jer_factor_down = TMath::Max(0., smearFactor_down);

				// multiplying a Gaussian should = to multiplying the sigma
				jet.setP4(jet.p4()*jer_factor);
				}

			// THUS MC was shifted by SF jer_factor -- it got corrected and the correction has to be propagated to MET
			full_jet_corr += jet.p4() - jet_initial_p4; // initial jet + this difference = corrected jet
			}

		NT_jet_resolution.push_back(jet_resolution);
		NT_jet_jer_sf.        push_back(jer_sf);
		NT_jet_jer_sf_up.     push_back(jer_sf_up);
		NT_jet_jer_sf_down.   push_back(jer_sf_down);

		NT_jet_matched_genjet_p4.push_back(gen_jet_p4);
		NT_jet_jer_factor.       push_back(jer_factor);
		NT_jet_jer_factor_up.    push_back(jer_factor_up);
		NT_jet_jer_factor_down.  push_back(jer_factor_down);

		// jet energy scale has uncertainty
		totalJESUnc->setJetEta(jet.eta());
		totalJESUnc->setJetPt(jet.pt()); // should be the corrected jet pt <- de hell this means? do it after MC SF?
		float relShift = fabs(totalJESUnc->getUncertainty(true));
		//jet.addUserFloat("jes_correction_relShift", relShift);
		// use it with rawJet*jes_correction*(1 +- relShift)
		// since all these corrections are multiplication of p4
		// I can do this shift whenever I want
		// uncertainty shift is saved only for the NOMINAL jet, which is default MiniAOD one now
		NT_jet_jes_uncertainty.push_back(relShift);

		/* just a note:
		 * so, my jets are default MiniAOD jets + MC SF for energy resolution
		 * the energy resolution UP/DOWN is saved and the energy scale uncertainty
		 *
		 * also I save JES factors: for reapplied Factorized Corrector and for uncorrected jet
		 * in principle I can jet nominal jets with reapplied Factorized Corrector, just without systematic uncertainties
		 * I save slimmedMETs[0] as main MET and propagate jet MC SF to it as met_corrected
		 * it corresponds to nominal jets
		 * there are other control mets
		 *
		 * it corresponds to LLRHiggs mets and jets
		 * and it is simple separation of different corrections and mets
		 */

		// save this jet -- it's used in tau dR match to jets
		selJetsNoLep.push_back(jet);
		/* requirements for there jets: pt and no leptons in dR
		 * no PFID requirement and no eta
		 */
		}
	NT_nalljets  = selJetsNoLep.size();

	NT_jets_full_correction = full_jet_corr; // initial jets + this correction = corrected jets
	// ALSO MET
	LorentzVector MET_corrected = NT_met_init - full_jet_corr;
	//float met_corrected = MET_corrected.pt();
	NT_met_corrected = MET_corrected;
	// TODO: I don't correct NT_met_slimmedMets and the other -- the correction should be applied offline if needed
	// in principle NT_met_init = NT_met_slimmedMets -- so NT_met_corrected saves the applied correction

	//pat::JetCollection selJets;
	//processJets_Kinematics(IDjets, /*bool isMC,*/ weight, jet_kino_cuts_pt, jet_kino_cuts_eta, selJets, false, false);

	//pat::JetCollection selJetsNoLep;
	//crossClean_in_dR(selJets, selLeptons, 0.4, selJetsNoLep, weight, string("selJetsNoLep"), false, false);
	// and these are output jets for NTuple
	// they pass ID, corrected with JEC (smeared JES for MC)
	// pass kinematic cuts (pt, eta)
	// and dR-cleaned from selected leptons

	//std::sort (selJetsNoLep.begin(),  selJetsNoLep.end(),  utils::sort_CandidatesByPt);
	// no need for sort


	/*
	 * TAUS
	 */
	//LogInfo("Demo") << "taus.size() = "<< taus.size();
	//string tau_Loose_ID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
	string tau_VLoose_ID ("byVLooseIsolationMVArun2v1DBoldDMwLT");
	string tau_Loose_ID  ("byLooseIsolationMVArun2v1DBoldDMwLT");
	string tau_Medium_ID ("byMediumIsolationMVArun2v1DBoldDMwLT");
	string tau_Tight_ID  ("byTightIsolationMVArun2v1DBoldDMwLT");
	string tau_VTight_ID ("byVTightIsolationMVArun2v1DBoldDMwLT");
	string tau_decayMode       ("decayModeFinding");
	string tau_againstMuon     ("againstMuonTight3");
	string tau_againstElectron ("againstElectronTightMVA6");

	pat::TauCollection IDtaus, selTaus;
	processTaus_ID    (taus,   weight, tau_decayMode, tau_againstMuon, tau_againstElectron, IDtaus, false, false);
	//processTaus_ID_ISO    (taus,   weight, tau_decayMode, tau_VLoose_ID, tau_againstMuon, tau_againstElectron, IDtaus, false, false);
	processTaus_Kinematics(IDtaus, weight, tau_kino_cuts_pt, tau_kino_cuts_eta, selTaus,      false, false);

	pat::TauCollection selTausNoLep;
	crossClean_in_dR(selTaus,       selLeptons, 0.4, selTausNoLep,        weight, string("selTausNoLep"),        false, false);

	// and these are the NT output taus
	std::sort (selTausNoLep.begin(),  selTausNoLep.end(),  utils::sort_CandidatesByPt);

	// PREPARE TRACKS AND STUFF for REFITTING TAU SV AND PV

	// BEAMSPOT
	//
	reco::BeamSpot beamSpot;
	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByToken(beamSpot_, beamSpotHandle);
	if (beamSpotHandle.isValid()) beamSpot = *beamSpotHandle;

	// TRACKS
	//
	edm::Handle<edm::View<pat::PackedCandidate>> tracksHandle;
	iEvent.getByToken(tracks_, tracksHandle);
	//if (tracksHandle.isValid()) tracks = *tracksHandle;
	const edm::View<pat::PackedCandidate>* track_cands = tracksHandle.product();

	// some more feature for tracks
	// the python conf load some module TransientTrackBuilder from cmssw package TrackingTools.TransientTrack
	// probably the module stores this stuff to the event..
	edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);

	// select tracks associated with PV
	reco::TrackCollection pvTracks;
	reco::TrackCollection allTracks; // for taus (with possible SV) (testing now)

	//TLorentzVector aTrack;
	for(size_t i=0; i<track_cands->size(); ++i)
		{
		if((*track_cands)[i].charge()==0 || (*track_cands)[i].vertexRef().isNull()) continue;
		if(!(*track_cands)[i].bestTrack()) continue;

		unsigned int key = (*track_cands)[i].vertexRef().key();
		int quality = (*track_cands)[i].pvAssociationQuality();

		// here I need to select "good" tracks
		// save them to all tracks
		// and if they belong to PV save them to pv tracks
		if (!(key!=0 ||
			(quality!=pat::PackedCandidate::UsedInFitTight
			 && quality!=pat::PackedCandidate::UsedInFitLoose)))// continue;
			{
			pvTracks.push_back(*((*track_cands)[i].bestTrack()));
			}

		// TODO: add requirement of "goodness"?
		allTracks.push_back(*((*track_cands)[i].bestTrack()));
		}

	//NT_ntaus = 0;
	NT_ntaus = selTausNoLep.size(); // all tau before MVA anti-jet iso
	std::vector<double > tracksToBeRemoved_PV; // compared by Pt due to the conflict of comparing const and not const iterators
	// these tracks correspond to all 3pi tau tracks in the event
	// without considering the ID of these taus
	// whenever refit works for a tau the corresponding tracks are removed
	for(size_t i=0; i<selTausNoLep.size(); ++i)
		{
		pat::Tau& tau = selTausNoLep[i];

		Int_t IDlev = 0;
		if (tau.tauID(tau_VTight_ID)) IDlev = 5;
		else if (tau.tauID(tau_Tight_ID))  IDlev = 4;
		else if (tau.tauID(tau_Medium_ID)) IDlev = 3;
		else if (tau.tauID(tau_Loose_ID))  IDlev = 2;
		else if (tau.tauID(tau_VLoose_ID)) IDlev = 1;

		//if (IDlev > 0) NT_ntaus += 1;

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

		// Re-reconstructed Secondary Vertex (SV) for taus with 3 pions
		// for each tau save the bool if fit is ok
		// and parameters,
		// default value for the parameters of not good fit are 999.
		// Also save matching quality (= sum of dR of the matching to general tracks).
		bool fitOk_SV = false;  
		std::vector<double > tracksToBeRemoved; // compared by Pt due to the conflict of comparing const and not const iterators
		// if fit ok some tracks have to be removed from PV refit
		TransientVertex transVtx_SV;
		double matchingQuality(0);
		if (tau.decayMode() > 9)
			{
			/*
			 * 1) unpack things necessary for the fitter, tau candidate tracks, beamSpot & general tracks
			 * 2) match tracks it tau candidates
			 * 3) fit adaptivefitter to the matched tracks
			 * 4) save with "quality score" = sum of dR of the matching
			 */

			// PAT tau tracks
			reco::CandidatePtrVector sigCands = tau.signalChargedHadrCands();//signalCands();

			// rebuild tau tracks (supposed o be more precise)
			// match allTracks (not pvTracks) to tau tracks
			std::vector<reco::TransientTrack> transTracks_tau;  

			for (reco::CandidatePtrVector::const_iterator itr = sigCands.begin(); itr != sigCands.end(); ++itr)
				{
				double deR(999.); 
				double checkqual(0);
				reco::Track closestTrack;

				//for(auto iter: pvTracks)
				for(auto iter: allTracks)
					{
					if(std::find(tracksToBeRemoved.begin(), tracksToBeRemoved.end(), iter.pt())!=tracksToBeRemoved.end())
						continue;
					if( sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2))  < deR)
						{
						deR = sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2));
						checkqual=deR;
						closestTrack = iter;
						}
					}

				matchingQuality+=checkqual;
				tracksToBeRemoved.push_back(closestTrack.pt());
				transTracks_tau.push_back(transTrackBuilder->build(closestTrack));
				//cout<<"  closestTrackiter eta  :  "<<   closestTrack.eta() << "   phi   " << closestTrack.phi() << "    pt  "<< closestTrack.pt() <<endl;
				}

			// do the tau SV refit
			if (transTracks_tau.size() >= 2 )
				{
				AdaptiveVertexFitter avf;
				avf.setWeightThreshold(0.001); 
				try
					{
					transVtx_SV = avf.vertex(transTracks_tau, beamSpot);
					fitOk_SV = true; 
					}
				catch (...)
					{
					fitOk_SV = false; 
					std::cout<<"Vtx fit failed!"<<std::endl;
					}
				}

			fitOk_SV = fitOk_SV && transVtx_SV.isValid() && fabs(transVtx_SV.position().x())<1 && fabs(transVtx_SV.position().y())<1;
			}

		// save the vertex
		//NT_tau_SV_fit_isOk.push_back(fitOk_SV);
		if(fitOk_SV)
			{
			// store the index of refited taus to the all taus vectors
			NT_tau_refited_index.push_back(NT_tau_SV_fit_matchingQuality.size()); // the last number in current vectors of refited taus

			// now add the new tau to refitted taus
			// set the tracks corresponding to tau to be removed for PV refit
			// tracksToBeRemoved_PV
			for (unsigned int i=0; i<tracksToBeRemoved.size(); i++)
				tracksToBeRemoved_PV.push_back(tracksToBeRemoved[i]);

			///NOTE: we take original vertex z position, as this gives the best reults on CP
			///variables. To be understood; probable reason are missing tracks with Pt<0.95GeV
			NT_tau_SV_fit_matchingQuality.push_back(matchingQuality);
			NT_tau_SV_fit_x.push_back(transVtx_SV.position().x());
			NT_tau_SV_fit_y.push_back(transVtx_SV.position().y());
			NT_tau_SV_fit_z.push_back(transVtx_SV.position().z());
			reco::Vertex secondaryVertex = transVtx_SV;
			// covariance matrix
			math::Error<3>::type svCov;
			secondaryVertex.fill(svCov);
			NT_tau_SV_cov.push_back(svCov);

			// and save the tracks (for kinematic fits on them)
			//transTracks_tau // <- not sure about using these for momentum TODO: check
			// using Signal Candidates of tau
			reco::CandidatePtrVector sigCands = tau.signalChargedHadrCands();
			// for control number of signal candidates (should only be = 3, the 3pi decay):
			NT_tau_SV_fit_ntracks.push_back(sigCands.size());

			int i = 0;
			for (reco::CandidatePtrVector::const_iterator itr = sigCands.begin(); itr != sigCands.end(); ++itr, i++)
				{
				// save as Same Sign track
				if ((*itr)->charge() * tau.pdgId() > 0)
					{
					NT_tau_SV_fit_track_SS_p4.push_back((*itr)->p4());
					}
				else if (i == 0) // first OS track
					{
					NT_tau_SV_fit_track_OS1_p4.push_back((*itr)->p4());
					}
				else // second OS track
					NT_tau_SV_fit_track_OS2_p4.push_back((*itr)->p4());
				}
			}
		else
			{
			// store default index of refited taus to the all taus vectors
			NT_tau_refited_index.push_back(-1);
			}
		}

	// if some tracks were removed (SV was refitted for a tau)
	// refit PV
	// build transient tracks of the rest of the tracks
	// for PV refitting
	if (tracksToBeRemoved_PV.size()>0)
		{
		TransientVertex transVtx_PV;
		std::vector<reco::TransientTrack> transTracks_nottau;

		for(auto iter: pvTracks)
			{
			if(std::find(tracksToBeRemoved_PV.begin(), tracksToBeRemoved_PV.end(), iter.pt())!=tracksToBeRemoved_PV.end())
				continue;
			transTracks_nottau.push_back(transTrackBuilder->build(iter));
			}

		bool fitOk_PV = false;  
		if (transTracks_nottau.size() >= 2 ) {
			AdaptiveVertexFitter avf;
			avf.setWeightThreshold(0.001); 
			try {
				transVtx_PV = avf.vertex(transTracks_nottau, beamSpot);
				fitOk_PV = true; 
				}
			catch (...) {
				fitOk_PV = false; 
				std::cout<<"Vtx fit failed!"<<std::endl;
				}
			}

		fitOk_PV = fitOk_PV && transVtx_PV.isValid() && fabs(transVtx_PV.position().x())<1 && fabs(transVtx_PV.position().y())<1;

		// save the vertex
		NT_PV_fit_isOk = fitOk_PV;
		if(fitOk_PV)
			{
			// save position and cov of vertex
			NT_PV_fit_x = transVtx_PV.position().x();
			NT_PV_fit_y = transVtx_PV.position().y();
			NT_PV_fit_z = transVtx_PV.position().z();
			reco::Vertex pvVertex = transVtx_PV;
			// covariance matrix
			//math::Error<3>::type pvCov;
			//pvVertex.fill(pvCov);
			//NT_PV_cov = pvCov;
			pvVertex.fill(NT_PV_cov);
			}
		}



	LogInfo ("Demo") << "all particles/objects are selected, nbjets = " << NT_nbjets;

	//bool record_ntuple = (isSingleMu || isSingleE || pass_dileptons) && NT_nbjets > 0 && NT_tau_IDlev.size() > 0; // at least 1 b jet and 1 loose tau
	bool pass_leptons = clean_lep_conditions && selLeptons.size() > 0 && selLeptons.size() < 3;
	bool record_ntuple = false;

	if (record_tauID)
		{
		/*
		 * tau ID preselection
		 * about twice less events then in our preselection with b jets
		 *
		 * should contain good WJets control sample
		 */
		record_ntuple |= pass_leptons && NT_ntaus > 0;
		}
	if (record_bPreselection)
		{
		/* leptons and at least 1 b jet
		 * our old preselection
		 */
		record_ntuple |= pass_leptons && NT_nbjets > 0;
		}
	if (record_MonitorHLT)
		{
		/* the HLT efficiency study
		 * only few events in lepMonitorTrigger
		 */
		record_ntuple |= pass_leptons && lepMonitorTrigger;
		}
	if (record_ElMu)
		{
		// all el-mu events (it's mainly TTbar, so should be few)
		record_ntuple |= clean_lep_conditions && abs(NT_leps_ID) == 143;
		}
	if (record_Dilep)
		{
		/* just all dileptons -- lots of DY here
		 * control for lepton IDs
		 * about = to tau ID preselection
		 */
		record_ntuple |= clean_lep_conditions && selLeptons.size() == 2;
		}

	if (record_ntuple)
		{
		LogInfo ("Demo") << "recording " << record_tauID << record_bPreselection << record_MonitorHLT << record_ElMu << record_Dilep;
		event_counter ->Fill(event_checkpoint++);
		weight_counter->Fill(event_checkpoint, weight);

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
if(!isMC){
	goodLumiFilter.FindLumiInFiles(urls); // urls! why are they here at all? no even 1 comment in that "Utilities!"
	goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
	}
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
