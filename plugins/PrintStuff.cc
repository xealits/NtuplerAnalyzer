// -*- C++ -*-
//
// Package:    UserCode/NtuplerAnalyzer
// Class:      PrintStuff
// 
/**\class PrintStuff PrintStuff.cc UserCode/NtuplerAnalyzer/plugins/PrintStuff.cc

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
//#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

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

#include "PhysicsTools/HepMCCandAlgos/interface/PDFWeightsHelper.h"

#include "TRandom3.h"

// MET recoil corrections for DY and WJets, from higgs->tautau group
// usage: https://github.com/CMS-HTT/RecoilCorrections/blob/master/instructions.txt
//#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
// do correction off-line in processing

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

static double pu_vector_NOMINAL[] = {0, 0.360609416811339, 0.910848525427002, 1.20629960507795, 0.965997726573782, 1.10708082813183, 1.14843491548622, 0.786526251164482, 0.490577792661333, 0.740680941110478,
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


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

//class PrintStuff : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
//class PrintStuff : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
//class PrintStuff : public edm::one::EDAnalyzer  { // didn't try this one!
class PrintStuff : public edm::EDAnalyzer  {
   public:
      explicit PrintStuff(const edm::ParameterSet&);
      ~PrintStuff();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void beginRun(edm::Run const &, edm::EventSetup const &); // override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endRun(edm::Run const& iRun, edm::EventSetup const&); // override;
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
	edm::EDGetTokenT<LHEEventProduct>   lheEPToken_;
	edm::EDGetTokenT<LHERunInfoProduct> lheRPToken_;
	edm::EDGetTokenT< std::vector < PileupSummaryInfo > > puInfo_;
	edm::EDGetTokenT< std::vector < PileupSummaryInfo > > puInfo2_;
	edm::EDGetTokenT<reco::GenParticleCollection> genParticle_;
	edm::EDGetTokenT<GenEventInfoProduct> evt_;
	PDFWeightsHelper pdfweightshelper_;

	edm::EDGetTokenT<edm::View<pat::PackedCandidate>> tracks_;
	edm::EDGetTokenT<reco::BeamSpot> beamSpot_;

	edm::EDGetTokenT<pat::METCollection> mets_slimmedMETs_, mets_slimmedMETsMuEGClean_, mets_uncorrected_;

	edm::EDGetTokenT<vector<reco::GenJet>> genJets_;
	//genJetsHandle.getByLabel(ev, "slimmedGenJets");
	edm::EDGetTokenT<pat::JetCollection> jets_;
	//jetsHandle.getByLabel(ev, "slimmedJets");
	edm::EDGetTokenT<std::vector<reco::GenJet>  > genJetsToken_; // there is a twist why

	//edm::EDGetTokenT<edm::ValueMap<float> > petersonFragToken_;
	edm::EDGetTokenT<edm::ValueMap<float> > upFragToken_, centralFragToken_, downFragToken_, PetersonFragToken_, semilepbrUpToken_, semilepbrDownToken_;

	//RecoilCorrector* recoilPFMetCorrector;
	//TH2D* zPtMass_histo;

	bool record_tauID, record_tauIDantiIso, record_bPreselection, record_MonitorHLT, record_ElMu, record_Dilep, record_jets;

	TString dtag;
	bool isMC, aMCatNLO, isWJets, isDY;
	bool isLocal;
	string  HLT_source,
		muHLT_MC1  , muHLT_MC2  ,
		muHLT_Data1, muHLT_Data2,
		elHLT_Data , elHLT_MC,
		lepMonitorHLT;

	jet_id    jetID;
	pu_jet_id jetPUID;

	TString jecDir;
	TString TjetResolutionFileName;
	TString TjetResolutionSFFileName;

	double  el_kino_cuts_pt, el_kino_cuts_eta, el_veto_kino_cuts_pt, el_veto_kino_cuts_eta,
		mu_kino_cuts_pt, mu_kino_cuts_eta, mu_veto_kino_cuts_pt, mu_veto_kino_cuts_eta,
		tau_kino_cuts_pt, tau_kino_cuts_eta;
	double jet_kino_cuts_pt, jet_kino_cuts_eta;
	double btag_threshold;

	edm::EDGetTokenT<bool> BadChCandFilterToken_;
	edm::EDGetTokenT<bool> BadPFMuonFilterToken_;


	// gotta be a new option for definition of the interface
	//#include "ntupleOutput_leps.h"
	// so, I need to
	//  - declare the interface as class atributes
	//  - then make the TTree and connect branches in constructor
	//  - and reset/fill stuff in analyze
	// let's do it first manually for p4-s of leptons
	/*
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > NT_lep_p4;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* pt_lep_p4; // yep, vectors of complex objects require additional persistent pointers
	*/
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
PrintStuff::PrintStuff(const edm::ParameterSet& iConfig) :
isMC       (iConfig.getParameter<bool>("isMC")),
HLT_source (iConfig.getParameter<string>("HLT_source"))

{

	// lheRPToken_ = consumes<LHERunInfoProduct, edm::InRun>(edm::InputTag("externalLHEProducer"));
	// this results in:
	// ----- Begin Fatal Exception 14-Dec-2017 18:36:01 CET-----------------------
	// An exception of category 'BranchTypeMismatch' occurred while
	//    [0] Processing run: 1
	//    [1] Calling global endRun for module PrintStuff/'printer'
	// Exception Message:
	// A get using a EDGetToken was done in Run but the consumes call was for Event.
	//  Please modify the consumes call to use the correct branch type.
	// ----- End Fatal Exception -------------------------------------------------
	// in
	// https://github.com/KappaAnalysis/Kappa/blob/master/Producers/interface/KGenInfoProducer.h#L51-L52
	// found this:

	lheRPToken_ = consumes<LHERunInfoProduct, edm::InRun>(edm::InputTag("externalLHEProducer"));

	trigResults_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", HLT_source));
	trigResultsRECO_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","RECO"));
	trigResultsPAT_     = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","PAT"));

	//lheEPToken_ = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));

	//now do what ever initialization is needed
	//usesResource("TFileService");

	/* output via EDM stuff
	 * does it work ok with skipping events?
	 * test the straight forward way, if doesn't work -- leave it, save in TTree manualy
	 *
	 * error: 'produces' was not declared in this scope
	 * -- so it's really just for producers
	 */
	//produces<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >>("leps_p4");

}


PrintStuff::~PrintStuff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PrintStuff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	LogInfo ("Demo") << "entered event " << iEvent.eventAuxiliary().run() << ' ' << iEvent.eventAuxiliary().luminosityBlock();

	edm::TriggerResultsByName tr = iEvent.triggerResultsByName (HLT_source);
	if (tr.isValid())
		{ // TODO: make a separate executable
		cout << "Printing " << HLT_source << " trigger list" << endl;
		int i = 0;
		for(edm::TriggerNames::Strings::const_iterator trnames = tr.triggerNames().begin(); trnames!=tr.triggerNames().end(); ++trnames, ++i)
			cout << i << "\t" << *trnames << endl;
		cout << "----------- End of trigger list ----------" << endl;
		}
	else
		{
		cout << "HLT trigger !isValid()" << endl;
		}

	// print MET Filters
	//edm::TriggerResultsByName metFilters = iEvent.triggerResultsByName("RECO"); //is present only if PAT (and miniAOD) is not run simultaniously with RECO
	//if(!metFilters.isValid()){metFilters = iEvent.triggerResultsByName("PAT");} //if not present, then it's part of RECO
	edm::TriggerResultsByName  patFilters = iEvent.triggerResultsByName("PAT");
	edm::TriggerResultsByName recoFilters = iEvent.triggerResultsByName("RECO");

	//LogInfo("Demo") << "Printing HLT trigger list";
	//LogInfo("Demo") << "-- Commented out --";
	std::cout << "Printing PAT and RECO filters" << std::endl;
	std::cout << "-- Cout --" << std::endl;

	if (patFilters.isValid())
		{
		//LogInfo("Demo") << "patFilters :";
		std::cout << "patFilters :" << std::endl;
		int i = 0;
		for (edm::TriggerNames::Strings::const_iterator trnames = patFilters.triggerNames().begin(); trnames!=patFilters.triggerNames().end(); ++trnames, ++i)
			std::cout << i << "\t" << *trnames << std::endl;
			//LogInfo("Demo") << i << "\t" << *trnames;
		}
	else
		{
		LogInfo("Demo") << "patFilters is not isValid()";
		}

	if (recoFilters.isValid())
		{
		std::cout << "recoFilters :" << std::endl;
		int i = 0;
		for (edm::TriggerNames::Strings::const_iterator trnames = recoFilters.triggerNames().begin(); trnames!=recoFilters.triggerNames().end(); ++trnames, ++i)
			std::cout << i << "\t" << *trnames << std::endl;
			//LogInfo("Demo") << i << "\t" << *trnames;
		}
	else
		{
		LogInfo("Demo") << "recoFilters is not isValid()";
		}

	//LogInfo("Demo") << "----------- End of trigger list ----------";
	std::cout << "----------- End of trigger list ----------" << std::endl;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ info on PDF weights content for MC  ------------
// the only example of this is:
// http://www-ekp.physik.uni-karlsruhe.de/~berger/kappa/doc/class_k_tuple.html#ab825d473aac95b1bb2872bf5624a01c9
// http://www-ekp.physik.uni-karlsruhe.de/~berger/kappa/doc/_k_tuple_8cc.html
// source
// http://www-ekp.physik.uni-karlsruhe.de/~berger/kappa/doc/_k_tuple_8cc_source.html
void PrintStuff::beginRun(edm::Run const & iRun, edm::EventSetup const & setup) 		
{
std::cout << "beginRun" << std::endl;
}

void PrintStuff::endRun(edm::Run const & iRun, edm::EventSetup const & setup) 		
{
std::cout << "endRun" << std::endl;

if (isMC)
	{
	edm::Handle<LHERunInfoProduct> run;
	iRun.getByToken(lheRPToken_, run);
	//iEvent.getByToken(lheRPToken_, run);

	//edm::Handle<LHERunInfoProduct> run; 
	//iRun.getByLabel("externalLHEProducer", run);

	LHERunInfoProduct myLHERunInfoProduct = *(run.product());

	std::cout << "labels of PDF weights" << std::endl;
	typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
	for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++)
		{
		std::cout << iter->tag() << std::endl;
		std::vector<std::string> lines = iter->lines();
		for (unsigned int iLine = 0; iLine<lines.size(); iLine++)
			{
			std::cout << lines.at(iLine);
			}
		}
	}
else
	{
	std::cout << "data input" << std::endl;
	}

}

// ------------ method called once each job just before starting event loop  ------------
void 
PrintStuff::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PrintStuff::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PrintStuff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrintStuff);
