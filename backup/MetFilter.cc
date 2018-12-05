//#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/NtuplerAnalyzer/interface/MetFilter.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

namespace patUtils
{

void MetFilter::FillBadEvents(std::string path){
     unsigned int Run=0; unsigned int Lumi=1; unsigned int Event=2;
     //LOOP on the files and fill the map
     std::ifstream File;
     File.open(path.c_str());
     if(!File.good() ){
       std::cout<<"ERROR:: File "<<path<<" NOT FOUND!!"<<std::endl; 
       return;
     }

     std::string line;
     while ( File.good() ){
         getline(File, line);
         std::istringstream stringfile(line);
         stringfile >> Run >> Lumi >> Event >> std::ws;
         map[RuLuEv(Run, Lumi, Event)] += 1;
    }
    File.close();
}
  

  bool MetFilter::passBadPFMuonFilter(const fwlite::Event& ev) {

    const bool debug_ = false; 

    const int algo_ = 14;
    //    const double minDZ_ = 0.1; 
    const double minTrkPtError_ = 0.5; 
    const double minMuPt_ = 100.; 

    pat::MuonCollection muons; 
    fwlite::Handle<pat::MuonCollection> muonsHandle; 
    muonsHandle.getByLabel(ev, "slimmedMuons"); 
    if(muonsHandle.isValid()){ muons = *muonsHandle;}

    pat::PackedCandidateCollection pfCandidates; 
    fwlite::Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
    pfCandidatesHandle.getByLabel(ev,"packedPFCandidates"); 
    if (pfCandidatesHandle.isValid()) { pfCandidates = *pfCandidatesHandle; }

    bool foundBadPFMuon = false;

    for ( unsigned i=0; i < muons.size(); ++i ) { // loop over all muons
    
      const reco::Muon & muon = muons[i];

      reco::TrackRef innerMuonTrack = muon.innerTrack();

      if ( innerMuonTrack.isNull() ) { 
	if (debug_) cout<<"Skipping this muon because it has no inner track"<<endl; 
	continue; 
      };

      if ( innerMuonTrack->pt() < minMuPt_) {
	if (debug_) cout<<"Skipping this muon because inner track pt."<<endl; 
	continue;
      }

      if ( innerMuonTrack->quality(reco::TrackBase::highPurity) ) { 
	if (debug_) cout<<"Skipping this muon because inner track is high purity."<<endl; 
	continue;
      }

      // Consider only muons with large relative pt error
      if (debug_) cout<<"Muon inner track pt rel err: "<<innerMuonTrack->ptError()/innerMuonTrack->pt()<<endl;
      if (not ( innerMuonTrack->ptError()/innerMuonTrack->pt() > minTrkPtError_ ) ) {
	if (debug_) cout<<"Skipping this muon because seems well measured."<<endl; 
	continue;
      }

      // Consider only muons from muonSeededStepOutIn algo
      if (debug_) cout<<"Muon inner track original algo: "<<innerMuonTrack->originalAlgo() << endl;
      if (not ( innerMuonTrack->originalAlgo() == algo_  && innerMuonTrack->algo() == algo_ ) ) {
	if (debug_) cout<<"Skipping this muon because is not coming from the muonSeededStepOutIn"<<endl; 
	continue;
      }
    
      for ( unsigned j=0; j < pfCandidates.size(); ++j ) {
	const reco::Candidate & pfCandidate = pfCandidates[j];
	// look for pf muon
	if ( not ( ( abs(pfCandidate.pdgId()) == 13) and (pfCandidate.pt() > minMuPt_) ) ) continue;
	// require small dR
	float dr = deltaR( muon.eta(), muon.phi(), pfCandidate.eta(), pfCandidate.phi() );
	if( dr < 0.001 ) {
	  foundBadPFMuon=true;
	  if (debug_) cout <<"found bad muon!"<<endl;
	  break;
	}
      }

      if (foundBadPFMuon) { break; };

    }

    bool pass = !foundBadPFMuon;
    if (debug_) cout<<"pass: "<<pass<<endl;

    return pass;

  }

  bool MetFilter::passBadChargedCandidateFilter(const fwlite::Event& ev) {

    const bool debug_ = false;

    const double maxDR_ = 0.001;
    const double minPtDiffRel_ = -0.5;
    const double minMuonTrackRelErr_ = 0.5;
    const double minMuonPt_ = 100.;

    //    reco::Muon MuonView;
    //    typedef View<reco::Muon> MuonView;
    pat::MuonCollection muons;
    fwlite::Handle<pat::MuonCollection> muonsHandle;
    muonsHandle.getByLabel(ev, "slimmedMuons");                                                                                          
    if(muonsHandle.isValid()){ muons = *muonsHandle;}  

    //reco::Candidate CandidateView;
    //    typedef View<reco::Candidate> CandidateView;
    pat::PackedCandidateCollection pfCandidates; 
    fwlite::Handle<pat::PackedCandidateCollection> pfCandidatesHandle;
    pfCandidatesHandle.getByLabel(ev,"packedPFCandidates");
    if (pfCandidatesHandle.isValid()) { pfCandidates = *pfCandidatesHandle; }

    bool foundBadChargedCandidate = false;

    for ( unsigned i=0; i < muons.size(); ++i ) { // loop over all muons

      const reco::Muon & muon = muons[i];

      if ( muon.pt() > minMuonPt_) {
	reco::TrackRef innerMuonTrack = muon.innerTrack();
        if (debug_) cout<<"muon "<<muon.pt()<<endl;

        if ( innerMuonTrack.isNull() ) { 
	  if (debug_) cout<<"Skipping this muon because it has no inner track"<<endl; 
	  continue; 
	};
        if ( innerMuonTrack->quality(reco::TrackBase::highPurity) ) { 
	  if (debug_) cout<<"Skipping this muon because inner track is high purity."<<endl; 
	  continue;
        }
        // Consider only muons with large relative pt error
        if (debug_) cout<<"Muon inner track pt rel err: "<<innerMuonTrack->ptError()/innerMuonTrack->pt()<<endl;
        if (not ( innerMuonTrack->ptError()/innerMuonTrack->pt() > minMuonTrackRelErr_ ) ) {
	  if (debug_) cout<<"Skipping this muon because seems well measured."<<endl; 
	  continue;
        }
        for ( unsigned j=0; j < pfCandidates.size(); ++j ) {
	  const reco::Candidate & pfCandidate = pfCandidates[j];
	  // look for charged hadrons
	  if (not ( abs(pfCandidate.pdgId()) == 211) ) continue;
	  float dr = deltaR( innerMuonTrack->eta(), innerMuonTrack->phi(), pfCandidate.eta(), pfCandidate.phi() );
	  float dpt = ( pfCandidate.pt() - innerMuonTrack->pt())/(0.5*(innerMuonTrack->pt() + pfCandidate.pt()));
	  if ( (debug_)  and (dr<0.5) ) cout<<" pt(it) "<<innerMuonTrack->pt()<<" candidate "<<pfCandidate.pt()<<" dr "<< dr
					    <<" dpt "<<dpt<<endl;
	  // require similar pt ( one sided ) and small dR
	  if ( ( deltaR( innerMuonTrack->eta(), innerMuonTrack->phi(), pfCandidate.eta(), pfCandidate.phi() ) < maxDR_ ) 
	       and ( ( pfCandidate.pt() - innerMuonTrack->pt())/(0.5*(innerMuonTrack->pt() + pfCandidate.pt())) > minPtDiffRel_ ) ) {
	    foundBadChargedCandidate = true;
	    cout <<"found bad track!"<<endl; 
	    break;
	  }
        }
        if (foundBadChargedCandidate) { break; };
      }
    } // end loop over muonss

    bool pass = !foundBadChargedCandidate;
    if (debug_) cout<<"pass: "<<pass<<endl;

    return pass;

  }
  

  int MetFilter::passMetFilterInt(const fwlite::Event& ev, bool is2016){
  
    if(map.find( RuLuEv(ev.eventAuxiliary().run(), ev.eventAuxiliary().luminosityBlock(), ev.eventAuxiliary().event()))!=map.end())return 1;                   

    edm::TriggerResultsByName metFilters = ev.triggerResultsByName("PAT");   //is present only if PAT (and miniAOD) is not run simultaniously with RECO       
    if(!metFilters.isValid()){metFilters = ev.triggerResultsByName("RECO");} //if not present, then it's part of RECO                                         
    if(!metFilters.isValid()){                                                                                                                                
      printf("TriggerResultsByName for MET filters is not found in the process, as a consequence the MET filter is disabled for this event\n");               
    }

    if(!utils::passTriggerPatterns(metFilters, "Flag_globalTightHalo2016Filter")) return 2;                                                                   
    if(!utils::passTriggerPatterns(metFilters, "Flag_goodVertices"              )) return 3;                                                                  
    if(!utils::passTriggerPatterns(metFilters, "Flag_eeBadScFilter"         )) return 4;                                                                      
    if(!utils::passTriggerPatterns(metFilters, "Flag_EcalDeadCellTriggerPrimitiveFilter")) return 5;                                                          
    if(!utils::passTriggerPatterns(metFilters, "Flag_HBHENoiseFilter"  )) return 6;                                                                           
    if(!utils::passTriggerPatterns(metFilters, "Flag_HBHENoiseIsoFilter"  )) return 7;          

    return 0;
  }    

  int MetFilter::passMetFilterInt(const fwlite::Event& ev){  
     if(map.find( RuLuEv(ev.eventAuxiliary().run(), ev.eventAuxiliary().luminosityBlock(), ev.eventAuxiliary().event()))!=map.end())return 1;

    // Legacy: the different collection name was necessary with the early 2015B prompt reco        
//    edm::TriggerResultsByName metFilters = isPromptReco ? ev.triggerResultsByName("RECO") : ev.triggerResultsByName("PAT");
    edm::TriggerResultsByName metFilters = ev.triggerResultsByName("PAT");   //is present only if PAT (and miniAOD) is not run simultaniously with RECO
    if(!metFilters.isValid()){metFilters = ev.triggerResultsByName("RECO");} //if not present, then it's part of RECO
    if(!metFilters.isValid()){       
      printf("TriggerResultsByName for MET filters is not found in the process, as a consequence the MET filter is disabled for this event\n");    
    }


    // Documentation:    
    // -------- Full MET filters list (see bin/chhiggs/runAnalysis.cc for details on how to print it out ------------------
    // Flag_trackingFailureFilter
    // Flag_goodVertices                        -------> Recommended by PAG
    // Flag_CSCTightHaloFilter                  -------> Recommended by PAG
    // Flag_trkPOGFilters
    // Flag_trkPOG_logErrorTooManyClusters
    // Flag_EcalDeadCellTriggerPrimitiveFilter  -------> Recommended by PAG
    // Flag_ecalLaserCorrFilter
    // Flag_trkPOG_manystripclus53X
    // Flag_eeBadScFilter                       -------> Recommended by PAG
    // Flag_METFilters
    // Flag_HBHENoiseFilter                     -------> Recommended by PAG
    // Flag_trkPOG_toomanystripclus53X
    // Flag_hcalLaserEventFilter
 

    // Notes (from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss_filters ):
    // - For the RunIISpring15DR74 MC campaing, the process name in PAT.
    // - For Run2015B PromptReco Data, the process name is RECO.
    // - For Run2015B re-MiniAOD Data 17Jul2015, the process name is PAT.
    // - MET filters are available in PromptReco since run 251585; for the earlier run range (251162-251562) use the re-MiniAOD 17Jul2015
    // - Recommendations on how to use MET filers are given in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2 . Note in particular that the HBHO noise filter must be re-run from MiniAOD instead of using the flag stored in the TriggerResults; this applies to all datasets (MC, PromptReco, 17Jul2015 re-MiniAOD)
    // -------------------------------------------------

    if(!utils::passTriggerPatterns(metFilters, "Flag_CSCTightHaloFilter"                     )) return 2; 
    if(!utils::passTriggerPatterns(metFilters, "Flag_goodVertices"                           )) return 3;
    if(!utils::passTriggerPatterns(metFilters, "Flag_eeBadScFilter"                          )) return 4;
    if(!utils::passTriggerPatterns(metFilters, "Flag_EcalDeadCellTriggerPrimitiveFilter"     )) return 5;
    if(!utils::passTriggerPatterns(metFilters, "Flag_HBHENoiseFilter"                        )) return 6;

/*    
    // HBHE filter needs to be complemented with , because it is messed up in data (see documentation below)
    //if(!utils::passTriggerPatterns(metFilters, "Flag_HBHENoiseFilter"   )) return false; // Needs to be rerun for both data (prompt+reReco) and MC, for now.
    // C++ conversion of the python FWLITE example: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2

    HcalNoiseSummary summary;
    fwlite::Handle <HcalNoiseSummary> summaryHandle;
    summaryHandle.getByLabel(ev, "hcalnoise");
    if(summaryHandle.isValid()) summary=*summaryHandle;

    //HBHE NOISE
    if(summary.maxHPDHits() >= 17)return 6;
    if(summary.maxHPDNoOtherHits() >= 10)return 7;
    if( summary.maxZeros() >= 9e9) return 8;
    if(summary.HasBadRBXRechitR45Loose())return 9;  //for 25ns only!  CHeck the recipe for 50ns.
//    bool failCommon(summary.maxHPDHits() >= 17  || summary.maxHPDNoOtherHits() >= 10 || summary.maxZeros() >= 9e9 );
    // IgnoreTS4TS5ifJetInLowBVRegion is always false, skipping.
//    if((failCommon || summary.HasBadRBXRechitR45Loose() )) return 5;  //for 25ns only
   
     //HBHE ISO NOISE
     if(summary.numIsolatedNoiseChannels() >=10)return 10;
     if(summary.isolatedNoiseSumE() >=50) return 11;
     if(summary.isolatedNoiseSumEt() >=25) return 12;
*/
     return 0;
}

  bool MetFilter::passMetFilter(const fwlite::Event& ev){
    return passMetFilterInt(ev)==0;
  }

bool MetFilter::BadGlobalMuonTaggerFilter(const fwlite::Event& ev,std::unique_ptr<std::vector<reco::Muon*>> &out1, bool selectClones_){

        reco::VertexCollection vtx;
        fwlite::Handle< reco::VertexCollection > vtxHandle; 
        vtxHandle.getByLabel(ev, "offlineSlimmedPrimaryVertices");
        if(vtxHandle.isValid()){ vtx = *vtxHandle;}	

	pat::MuonCollection muons; 
	fwlite::Handle< pat::MuonCollection > muonsHandle;
        muonsHandle.getByLabel(ev, "slimmedMuons");
	if(muonsHandle.isValid()){ muons = *muonsHandle;}

	std::vector<int> goodMuon;
 	const auto &PV = vtx.front().position();
	std::unique_ptr<std::vector<reco::Muon*>> out(new std::vector<reco::Muon*>());
        out = (std::move(out1));

       double ptCut_=20;
       for ( unsigned i=0; i < muons.size(); ++i ) { // loop over all muons
       const reco::Muon &mu = muons[i];

        if (!mu.isPFMuon() || mu.innerTrack().isNull()) {
            goodMuon.push_back(-1); // bad but we don't care
            continue;
        } 
        if (preselection(mu,selectClones_)) {
            float dxypv = std::abs(mu.innerTrack()->dxy(PV));
            float dzpv  = std::abs(mu.innerTrack()->dz(PV));
            if (tighterId(mu)) {
                bool ipLoose = ((dxypv < 0.5 && dzpv < 2.0) || mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() >= 2);
                goodMuon.push_back(ipLoose || (!selectClones_ && tightGlobal(mu)));
            } else if (safeId(mu)) {
                bool ipTight = (dxypv < 0.2 && dzpv < 0.5);
                goodMuon.push_back(ipTight);
           } else {
                goodMuon.push_back(0);
            }
        } else {
            goodMuon.push_back(3); // maybe good, maybe bad, but we don't care
        }
    }

    bool found = false;
    bool verbose_ = false;
    for (unsigned int i = 0, n = muons.size(); i < n; ++i) {
        if (muons[i].pt() < ptCut_ || goodMuon[i] != 0) continue;
        if (verbose_) printf("potentially bad muon %d of pt %.1f eta %+.3f phi %+.3f\n", int(i+1), muons[i].pt(), muons[i].eta(), muons[i].phi());
        bool bad = true;
        if (selectClones_) {
            bad = false; // unless proven otherwise
            unsigned int n1 = muons[i].numberOfMatches(reco::Muon::SegmentArbitration);
            for (unsigned int j = 0; j < n; ++j) {
                if (j == i || goodMuon[j] <= 0 || !partnerId(muons[j])) continue;
                unsigned int n2 = muons[j].numberOfMatches(reco::Muon::SegmentArbitration);
                if (deltaR2(muons[i],muons[j]) < 0.16 || (n1 > 0 && n2 > 0 && muon::sharedSegments(muons[i],muons[j]) >= 0.5*std::min(n1,n2))) {
                    if (verbose_) printf("     tagged as clone of muon %d of pt %.1f eta %+.3f phi %+.3f\n", int(j+1), muons[j].pt(), muons[j].eta(), muons[j].phi());
                    bad = true;
                    break;
                } 
            }
        }
        if (bad) {
            found = true;
            out->push_back(std::addressof(muons[i]));
        }
}
	return found;
}

}




