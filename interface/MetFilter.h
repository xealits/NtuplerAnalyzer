#ifndef METFILTER_h
#define METFILTER_h


#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//Load here all the dataformat that we will need
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/PtrVector.h"

//need for the good lumi filter
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"
#include "FWCore/Utilities/interface/Algorithms.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "UserCode/NtuplerAnalyzer/interface/MacroUtils.h"
#include "UserCode/NtuplerAnalyzer/interface/LumiUtils.h"

// Electron ID
#include "RecoEgamma/ElectronIdentification/interface/VersionedPatElectronSelector.h"

#include <vector>
#include "TVector3.h"
#include "TMath.h"
#include "TGraph.h"
#include <Math/VectorUtil.h>


namespace patUtils
{
   class MetFilter{
    private :
     struct RuLuEv {
        unsigned int Run;  unsigned int Lumi;  unsigned int Event;
        RuLuEv(unsigned int Run_, unsigned int Lumi_, unsigned int Event_){ Run = Run_; Lumi = Lumi_; Event = Event_;}
        bool operator==(const RuLuEv &other) const { return (Run == other.Run && Lumi == other.Lumi && Event == other.Event); }
     };
     struct RuLuEvHasher{
         std::size_t operator()(const RuLuEv& k) const{ using std::size_t; using std::hash;  using std::string;
            return ((hash<unsigned int>()(k.Run) ^ (hash<unsigned int>()(k.Lumi) << 1)) >> 1) ^ (hash<unsigned int>()(k.Event) << 1);
         }
     };

     typedef std::unordered_map<RuLuEv, int, RuLuEvHasher> MetFilterMap;
     MetFilterMap map;
    public :
     MetFilter(){}
     ~MetFilter(){}
     void Clear(){map.clear();}
     void FillBadEvents(std::string path);

     //     New Met Filters for 2016 Run II:    
     bool passBadPFMuonFilter(const fwlite::Event& ev);
     bool passBadChargedCandidateFilter(const fwlite::Event& ev);   
     
     int  passMetFilterInt(const fwlite::Event& ev);
     int  passMetFilterInt(const fwlite::Event& ev, bool is2016); 
     bool passMetFilter(const fwlite::Event& ev);
     bool BadGlobalMuonTaggerFilter(const fwlite::Event& ev,std::unique_ptr<edm::PtrVector<reco::Muon>> &out, bool selectClones=false);
     bool BadGlobalMuonTaggerFilter(const fwlite::Event& ev,std::unique_ptr<std::vector<reco::Muon*>> &out, bool selectClones=false);
   };
}


#endif
