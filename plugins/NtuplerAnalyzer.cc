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
// couple functions processing leptons
#include "UserCode/ttbar-leptons-80X/interface/ProcessingMuons.h"                                                                                                                                                                  
#include "UserCode/ttbar-leptons-80X/interface/ProcessingElectrons.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingTaus.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingJets.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingBJets.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingHLT.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingGenParticles.h"

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
	edm::EDGetTokenT<reco::TrackCollection> tracks_;
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

	bool isMC;
	string  muHLT_MC1  , muHLT_MC2  ,
		muHLT_Data1, muHLT_Data2,
		elHLT_Data , elHLT_MC   ;

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
isMC(iConfig.getUntrackedParameter<bool>("isMC", true)),
muHLT_MC1  (iConfig.getParameter<std::string>("muHLT_MC1")),
muHLT_MC2  (iConfig.getParameter<std::string>("muHLT_MC2")),
muHLT_Data1(iConfig.getParameter<std::string>("muHLT_Data1")),
muHLT_Data2(iConfig.getParameter<std::string>("muHLT_Data2")),
elHLT_Data (iConfig.getParameter<std::string>("elHLT_Data")),
elHLT_MC   (iConfig.getParameter<std::string>("elHLT_MC"))

{
	/* define in constructor via call to consumes (magic thingy) */
	tracks_    = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
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

	Float_t NT_fixedGridRhoFastjetAll = -1;

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
	bool muTrigger = false;

	// TriggerNames for TriggerObjects --------------------
	edm::Handle<edm::TriggerResults> trigResults; //our trigger result object
	//edm::InputTag * trigResultsTag; // the tag object, trigResults are extracted from the event via this tag

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
		eTrigger =    ( isMC ?  utils::passTriggerPatterns(tr, elHLT_MC) : utils::passTriggerPatterns(tr, elHLT_Data));
		muTrigger =   ( isMC ?  utils::passTriggerPatterns(tr, muHLT_MC1, muHLT_MC2) : utils::passTriggerPatterns (tr, muHLT_Data1, muHLT_Data2));
		}

	if (!(eTrigger || muTrigger)) return; // orthogonalization is done afterwards

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
	vector<pat::TriggerObjectStandAlone> trig_objs = *triggerObjectsHandle;

	// objects of our triggers
	vector<pat::TriggerObjectStandAlone> el_trig_objs;
	vector<pat::TriggerObjectStandAlone> mu_trig_objs, mu_trig_objs2;

	if (eTrigger)
		{
		NT_HLT_el = true;
		Processing_selectHLTobjects(trig_objs, trigNames, el_trig_objs, (isMC? elHLT_MC : elHLT_Data));
		}
	if (muTrigger)
		{
		NT_HLT_mu = true;
		Processing_selectHLTobjects(trig_objs, trigNames, mu_trig_objs,  (isMC? muHLT_MC1 : muHLT_Data1));
		Processing_selectHLTobjects(trig_objs, trigNames, mu_trig_objs2, (isMC? muHLT_MC2 : muHLT_Data2));
		// vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
		mu_trig_objs.insert(mu_trig_objs.end(), mu_trig_objs2.begin(), mu_trig_objs2.end());
		}

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
						int d0_id = p.daughter(0)->pdgId();
						int d1_id = p.daughter(1)->pdgId();
						int W_num = d0_id == 24 ? 0 : (d1_id == 24 ? 1 : -1) ;
						if (W_num < 0) continue;
						const reco::Candidate * W = p.daughter( W_num );
						const reco::Candidate * W_final = find_W_decay(W);
						int decay_id = 1;
						// = id of lepton or 1 for quarks
						if (fabs(W_final->daughter(0)->pdgId()) == 11 || fabs(W_final->daughter(0)->pdgId()) == 13 || fabs(W_final->daughter(0)->pdgId()) == 15)
							decay_id = W_final->daughter(0)->pdgId();
						else if (fabs(W_final->daughter(1)->pdgId()) == 11 || fabs(W_final->daughter(1)->pdgId()) == 13 || fabs(W_final->daughter(1)->pdgId()) == 15)
							decay_id = W_final->daughter(1)->pdgId();

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
				if (id == 15 && st == 1)
					NT_gen_tau_p4.push_back(p.p4()); 
				else if (id == 15 && st == 2)
					{
					// it's a final state tau
					// select its' daughters, skipping neutrinos
					// add their momenta -- use the sum as a visible_gen_tau
					LorentzVector vis_ds(0,0,0,0);
					for (int j = 0; j < n_daughters; ++j)
						{
						const reco::Candidate * d = p.daughter(j);
						unsigned int d_id = abs(d->pdgId());
						if (d_id == 12 || d_id == 14 || d_id == 16) continue;
						vis_ds += d->p4();
						}
					NT_gen_tau_p4.push_back(vis_ds); 
					}
				}
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

	bool clean_lep_conditions = nVetoE==0 && nVetoMu==0 && nGoodPV != 0;
	if (!(clean_lep_conditions && selLeptons.size() > 0 && selLeptons.size() < 3)) return; // exit now to reduce computation

	/*
	 * TAUS
	 */
	//LogInfo("Demo") << "taus.size() = "<< taus.size();
	//string tau_Loose_ID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
	string tau_Loose_ID("byLooseIsolationMVArun2v1DBoldDMwLT");
	string tau_decayMode       ("decayModeFinding");
	string tau_againstMuon     ("againstMuonTight3");
	string tau_againstElectron ("againstElectronTightMVA6");
	double tau_kino_cuts_pt  = 30.;
	double tau_kino_cuts_eta = 2.3;

	pat::TauCollection IDtaus, selTaus;
	processTaus_ID_ISO    (taus,   weight, tau_decayMode, tau_Loose_ID, tau_againstMuon, tau_againstElectron, IDtaus, false, false);
	processTaus_Kinematics(IDtaus, weight, tau_kino_cuts_pt, tau_kino_cuts_eta, selTaus,      false, false);

	pat::TauCollection selTausNoLep;
	crossClean_in_dR(selTaus,       selLeptons, 0.4, selTausNoLep,        weights_FULL[SYS_NOMINAL], string("selTausNoLep"),        false, debug);

	// and these are the NT output taus
	std::sort (selTausNoLep.begin(),  selTausNoLep.end(),  utils::sort_CandidatesByPt);

	for(size_t i=0; i<selTausNoLep.size(); ++i)
		{
		pat::Tau& tau = selTausNoLep[i];

		Int_t IDlev = 0;
		if (tau.tauID(tau_Tight_ID)) IDlev = 3;
		else if (tau.tauID(tau_ID)) IDlev = 2;
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
			NT_tau_flightLength.push_back(-1);
			NT_tau_flightLengthSignificance.push_back(-1);
			}
		}

	// JETS

	//bool record_ntuple = (isSingleMu || isSingleE || pass_dileptons) && NT_nbjets > 0 && NT_tau_IDlev.size() > 0; // at least 1 b jet and 1 loose tau
	bool record_ntuple = clean_lep_conditions && selLeptons.size() > 0 && selLeptons.size() < 3 && NT_nbjets > 0; // leptons and at least 1 b jet

	if (record_ntuple)
		{
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
