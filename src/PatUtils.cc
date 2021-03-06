//#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/NtuplerAnalyzer/interface/PatUtils.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

namespace patUtils
{

  bool passId (VersionedPatElectronSelector id, edm::EventBase const & event, pat::Electron el){
    // This assumes an object to be created ( *before the event loop* ):
    // VersionedPatElectronSelector loose_id("some_VID_tag_including_the_WP");
    return id(el, event);
  }
  
  bool passId(pat::Electron& el,  reco::Vertex& vtx, int IdLevel, int cutVersion){
    
    //for electron Id look here: //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns 
    //for the meaning of the different cuts here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
    float dEtaln         = fabs(el.deltaEtaSuperClusterTrackAtVtx());
    float dPhiln         = fabs(el.deltaPhiSuperClusterTrackAtVtx());
    float sigmaletaleta  = el.sigmaIetaIeta();
    float hem            = el.hadronicOverEm();
    double resol         = fabs((1/el.ecalEnergy())-(el.eSuperClusterOverP()/el.ecalEnergy()));
    double dxy           = fabs(el.gsfTrack()->dxy(vtx.position()));
    double dz            = fabs(el.gsfTrack()->dz(vtx.position())); 
    //double mHits         = el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    double mHits         = el.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
    bool conversionVeto = el.passConversionVeto(); 
    bool barrel = (fabs(el.superCluster()->eta()) <= 1.479);
    bool endcap = (!barrel && fabs(el.superCluster()->eta()) < 2.5);
    
    // Spring15 selection 
    switch(cutVersion){
    case CutVersion::Spring15Cut25ns :

            switch(IdLevel){
            case llvvElecId::Veto :
                 if(barrel                   &&
                 dEtaln        < 0.0152   &&
                 dPhiln        < 0.216    &&
                 sigmaletaleta < 0.0114   &&
                 hem           < 0.181    &&
                 dxy           < 0.0564   &&
                 dz            < 0.472    &&
                 resol         < 0.207    &&
                 mHits         <=2        &&
                 conversionVeto            )
        	return true;
                 if(endcap                   &&
                 dEtaln        < 0.0113   &&
                 dPhiln        < 0.237    &&
                 sigmaletaleta < 0.0352   &&
                 hem           < 0.116    &&
                 dxy           < 0.222    &&
                 dz            < 0.921    &&
                 resol         < 0.174    &&
                 mHits <= 3               &&
                 conversionVeto            )
                return true;
        	break;
              
              case llvvElecId::Loose :
              if(barrel                   &&
                 dEtaln        < 0.0105   &&
                 dPhiln        < 0.115    &&
                 sigmaletaleta < 0.0103   &&
                 hem           < 0.104    &&
                 dxy           < 0.0261   &&
                 dz            < 0.41     &&
                 resol         < 0.102    &&
                 mHits         <= 2       &&
                 conversionVeto            )
                return true;
              if(endcap                   &&
                 dEtaln        < 0.00814  &&
                 dPhiln        < 0.182    &&
                 sigmaletaleta < 0.0301   &&
                 hem           < 0.0897   &&
                 dxy           < 0.118    &&
                 dz            < 0.822    &&
                 resol         < 0.126    &&
                 mHits         <= 1       &&
                 conversionVeto            )
                return true;
              break;
        
            case llvvElecId::Medium :
              if(barrel                     &&
                 dEtaln          < 0.0103   &&
                 dPhiln          < 0.0336   &&
                 sigmaletaleta   < 0.0101   &&
                 hem             < 0.0876   &&
                 dxy             < 0.0118   &&
                 dz              < 0.373    &&
                 resol           < 0.0174   &&
                 mHits           <= 2       &&
                 conversionVeto             )
                return true;
              if(endcap                     &&
                 dEtaln          < 0.00733  &&
                 dPhiln          < 0.114    &&
                 sigmaletaleta   < 0.0283   &&
                 hem             < 0.0678   &&
                 dxy             < 0.0739   &&
                 dz              < 0.602    &&
                 resol           < 0.0898   &&
                 mHits            <= 1      &&
                 conversionVeto             )
                return true;
              break;
        
            case llvvElecId::Tight :
              if(barrel                     &&
                 dEtaln          < 0.00926  &&
                 dPhiln          < 0.0336   &&
                 sigmaletaleta   < 0.0101   &&
                 hem             < 0.0597   &&
                 dxy             < 0.0111   &&
                 dz              < 0.0466   &&
                 resol           < 0.012    &&
                 mHits           <= 2       &&
                 conversionVeto             )
                return true;
              if(endcap                     &&
                 dEtaln          < 0.00724  &&
                 dPhiln          < 0.0918   &&
                 sigmaletaleta   < 0.0279   &&
                 hem             < 0.0615   &&
                 dxy             < 0.0351   &&
                 dz              < 0.417    &&
                 resol           < 0.00999  &&
                 mHits           <= 1       &&
                 conversionVeto             )
                return true;
              break;
        
            case llvvElecId::LooseMVA :
            case llvvElecId::MediumMVA :
            case llvvElecId::TightMVA :
              printf("FIXME: MVA ID not yet implemented for the electron\n");
              return false;
                          break;
                          
            default:
              printf("FIXME ElectronId llvvElecId::%i is unkown\n", IdLevel);
              return false;
              break;
            }
	break;

    case CutVersion::Moriond17Cut :
    case CutVersion::ICHEP16Cut :

            switch(IdLevel){
            case llvvElecId::Veto :
                 if(barrel                   &&
                 dEtaln        < 0.00749  &&
                 dPhiln        < 0.228    &&
                 sigmaletaleta < 0.0115   &&
                 hem           < 0.356    &&
                 //dxy           < 0.050    &&
                 //dz            < 0.10     &&
                 resol         < 0.299    &&
                 mHits         <=2        &&
                 conversionVeto            )
        	return true;
                 if(endcap                   &&
                 dEtaln        < 0.00895  &&
                 dPhiln        < 0.213    &&
                 sigmaletaleta < 0.037    &&
                 hem           < 0.211    &&
                 //dxy           < 0.10     &&
                 //dz            < 0.20     &&
                 resol         < 0.15     &&
                 mHits <= 3               &&
                 conversionVeto            )
                return true;
        	break;
              
              case llvvElecId::Loose :
              if(barrel                   &&
                 dEtaln        < 0.00477  &&
                 dPhiln        < 0.222    &&
                 sigmaletaleta < 0.011    &&
                 hem           < 0.298    &&
                 //dxy           < 0.05     &&
                 //dz            < 0.10     &&
                 resol         < 0.241    &&
                 mHits         <= 1       &&
                 conversionVeto            )
                return true;
              if(endcap                   &&
                 dEtaln        < 0.00868  &&
                 dPhiln        < 0.213    &&
                 sigmaletaleta < 0.0314   &&
                 hem           < 0.101    &&
                 //dxy           < 0.10     &&
                 //dz            < 0.20     &&
                 resol         < 0.14     &&
                 mHits         <= 1       &&
                 conversionVeto            )
                return true;
              break;
        
            case llvvElecId::Medium :
              if(barrel                     &&
                 dEtaln          < 0.00311  &&
                 dPhiln          < 0.103    &&
                 sigmaletaleta   < 0.00998  &&
                 hem             < 0.253    &&
                 //dxy             < 0.05     &&
                 //dz              < 0.10     &&
                 resol           < 0.134    &&
                 mHits           <= 1       &&
                 conversionVeto             )
                return true;
              if(endcap                     &&
                 dEtaln          < 0.00609  &&
                 dPhiln          < 0.045    &&
                 sigmaletaleta   < 0.0298   &&
                 hem             < 0.0878   &&
                 //dxy             < 0.10     &&
                 //dz              < 0.20     &&
                 resol           < 0.13     &&
                 mHits            <= 1      &&
                 conversionVeto             )
                return true;
              break;
        
            case llvvElecId::Tight :
              if(barrel                     &&
                 dEtaln          < 0.00308  &&
                 dPhiln          < 0.0816   &&
                 sigmaletaleta   < 0.00998  &&
                 hem             < 0.0414   &&
                 //dxy             < 0.05     &&
                 //dz              < 0.10     &&
                 resol           < 0.0129   &&
                 mHits           <= 1       &&
                 conversionVeto             )
                return true;
              if(endcap                     &&
                 dEtaln          < 0.00605  &&
                 dPhiln          < 0.0394   &&
                 sigmaletaleta   < 0.0292   &&
                 hem             < 0.0641   &&
                 //dxy             < 0.10     &&
                 //dz              < 0.20     &&
                 resol           < 0.0129   &&
                 mHits           <= 1       &&
                 conversionVeto             )
                return true;
              break;
        
            case llvvElecId::LooseMVA :
            case llvvElecId::MediumMVA :
            case llvvElecId::TightMVA :
              printf("FIXME: MVA ID not yet implemented for the electron\n");
              return false;
                          break;
                          
            default:
              printf("FIXME ElectronId llvvElecId::%i is unkown\n", IdLevel);
              return false;
              break;
            }
      break;

    default:
      printf("FIXME CutVersion::%i is unkown\n", cutVersion);
      return false;
      break;
    }
    return false;
  }
  
  bool passId(pat::Muon& mu,  reco::Vertex& vtx, int IdLevel, int cutVersion){
    //for muon Id look here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#LooseMuon
   
    // Spring15 selection 
    switch(cutVersion){
    case CutVersion::Spring15Cut25ns :

            switch(IdLevel){
              
            case llvvMuonId::Loose :
              if(mu.isPFMuon() && (mu.isGlobalMuon() || mu.isTrackerMuon()))return true;
              break;
              
            case llvvMuonId::Soft :
              if(mu.isPFMuon() && mu.isTrackerMuon() && mu.muonID("TMOneStationTight") && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 1 &&
                 fabs(mu.innerTrack()->dxy(vtx.position())) < 0.3 && fabs(mu.innerTrack()->dz(vtx.position())) < 20. && mu.innerTrack()->normalizedChi2() < 1.8) return true;
              break;
              
            case llvvMuonId::Tight :
              if( mu.isPFMuon() && mu.isGlobalMuon() && mu.globalTrack()->normalizedChi2() < 10. && mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0. && mu.numberOfMatchedStations() > 1 &&
                  fabs(mu.muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(mu.muonBestTrack()->dz(vtx.position())) < 0.5 && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
                  mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5)return true;
              break;
              
            case llvvMuonId::StdLoose :
              if(mu.isLooseMuon()) return true;
              break;
              
            case llvvMuonId::StdSoft :
              if(mu.isSoftMuon(vtx)) return true;
              break;
              
            case llvvMuonId::StdTight :
              if(mu.isTightMuon(vtx)) return true;
              break;
              
            default:
              printf("FIXME MuonId llvvMuonId::%i is unkown\n", IdLevel);
              return false;
              break;
            }
        break;
	// ICHEP16 or Moriond17 selection
    case CutVersion::Moriond17Cut :
    case CutVersion::ICHEP16Cut :
        
            switch(IdLevel){
              
            case llvvMuonId::Loose :
              if(mu.isPFMuon() && (mu.isGlobalMuon() || mu.isTrackerMuon()))return true;
              break;
              
            case llvvMuonId::Soft :
              if(muon::isGoodMuon(mu, muon::TMOneStationTight) && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 &&
                 fabs(mu.innerTrack()->dxy(vtx.position())) < 0.3 && fabs(mu.innerTrack()->dz(vtx.position())) < 20.) return true;
              break;
              
            case llvvMuonId::Tight :
              if( mu.isPFMuon() && mu.isGlobalMuon() && mu.globalTrack()->normalizedChi2() < 10. && mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0. && mu.numberOfMatchedStations() > 1 &&
                  fabs(mu.muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(mu.muonBestTrack()->dz(vtx.position())) < 0.5 && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 &&
                  mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5)return true;
              break;

            case llvvMuonId::tkHighPT :
              if(mu.isTrackerMuon() && mu.track().isNonnull() && mu.numberOfMatchedStations() > 1 && (mu.muonBestTrack()->ptError()/mu.muonBestTrack()->pt()) < 0.3 && fabs(mu.muonBestTrack()->dxy(vtx.position()))<0.2 && fabs(mu.muonBestTrack()->dz(vtx.position())) < 0.5 && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement()>5) return true;
              break;
            case llvvMuonId::TightAndTlkHighPt :
	      {
	      bool tightID =  mu.isPFMuon() && mu.isGlobalMuon() && mu.globalTrack()->normalizedChi2() < 10. && mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0. && mu.numberOfMatchedStations() > 1 && fabs(mu.muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(mu.muonBestTrack()->dz(vtx.position())) < 0.5 && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5; 
              bool tkHighPt = mu.isTrackerMuon() && mu.track().isNonnull() && mu.numberOfMatchedStations() > 1 && (mu.muonBestTrack()->ptError()/mu.muonBestTrack()->pt()) < 0.3 && fabs(mu.muonBestTrack()->dxy(vtx.position()))<0.2 && fabs(mu.muonBestTrack()->dz(vtx.position())) < 0.5 && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement()>5;
              if (tightID||(mu.pt()>200&&tkHighPt)) return true;
              }
            case llvvMuonId::StdLoose :
              if(mu.isLooseMuon()) return true;
              break;
              
            case llvvMuonId::StdSoft :
              if(mu.isSoftMuon(vtx)) return true;
              break;
              
            case llvvMuonId::StdTight :
              if(mu.isTightMuon(vtx)) return true;
              break;
              
            default:
              printf("FIXME MuonId llvvMuonId::%i is unkown\n", IdLevel);
              return false;
              break;
            }
	break;
    default:
        printf("FIXME MuonID CutVersion::%i is unkown\n", cutVersion);
        return false;
        break;
    }

    return false;
  }  
  
  bool passId(pat::Photon& photon, double rho, int IdLevel){
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2
    // CSA14 selection, conditions: 25ns, better detector alignment. 
    // Used Savvas Kyriacou's slides, mailed from Ilya. 
    
    bool elevto = photon.hasPixelSeed();  //LQ  REACTIVATED FOR TIGHT ID, OTHERWISE MANY ELECtRONS pass the photon Id
    
    // sigma ieta ieta in 5x5 shower block
     float sigmaIetaIeta = photon.full5x5_sigmaIetaIeta();

    // H/E 
    float hoe = photon.hadTowOverEm();

    // isolation
    double pt=photon.pt();
    double eta=photon.superCluster()->eta();

    float chIso = photon.chargedHadronIso(); 
    float chArea = utils::cmssw::getEffectiveArea(22,eta,"chIso"); 

    float nhIso = photon.neutralHadronIso();
    float nhArea = utils::cmssw::getEffectiveArea(22,eta,"nhIso");

    float gIso = photon.photonIso();
    float gArea = utils::cmssw::getEffectiveArea(22,eta,"gIso");

    bool barrel = (fabs(eta) <= 1.479);
    bool endcap = (!barrel && fabs(eta) < 2.5);

    //SPRING16 selection :  https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Recommended_Working_points_for_2
          switch(IdLevel){
          case llvvPhotonId::Loose :
              
            if ( barrel
      	   //&& !elevto
      	   && hoe < 0.0597
      	   && sigmaIetaIeta < 0.01031
      	   && TMath::Max(chIso-chArea*rho,0.0) < 1.295
      	   && TMath::Max(nhIso-nhArea*rho,0.0) < 10.910 + 0.0148*pt + 0.000017*pt*pt
      	   && TMath::Max(gIso-gArea*rho,  0.0) < 3.630 + 0.0047*pt  )
      	return true; 
            if ( endcap
      	   //&& !elevto
      	   && hoe < 0.0481
      	   && sigmaIetaIeta < 0.03013
      	   && TMath::Max(chIso-chArea*rho,0.0) < 1.011
      	   && TMath::Max(nhIso-nhArea*rho,0.0) < 5.931 + 0.0163*pt+0.000014*pt*pt
      	   && TMath::Max(gIso-gArea*rho,  0.0) < 6.641 + 0.0034*pt  )
      	return true; 
                  
            break;
            
          case llvvPhotonId::Medium :
      
            if ( barrel
      	   //&& !elevto
      	   && hoe < 0.0396
      	   && sigmaIetaIeta < 0.01022
      	   && TMath::Max(chIso-chArea*rho,0.0) < 0.441 
      	   && TMath::Max(nhIso-nhArea*rho,0.0) < 2.725 + 0.0148*pt + 0.000017*pt*pt 
      	   && TMath::Max(gIso-gArea*rho,  0.0) < 2.571 + 0.0047*pt  )
      	return true; 
            if ( endcap
      	   //&& !elevto
      	   && hoe < 0.0219
      	   && sigmaIetaIeta < 0.03001
      	   && TMath::Max(chIso-chArea*rho,0.0) < 0.442
      	   && TMath::Max(nhIso-nhArea*rho,0.0) < 1.715 + 0.0163*pt+0.000014*pt*pt
      	   && TMath::Max(gIso-gArea*rho,  0.0) < 3.863 + 0.0034*pt  )
      	return true; 
                  
            break;
          case llvvPhotonId::Tight :
            
            if ( barrel   
           && !elevto 
      	   && hoe < 0.0269      
           && sigmaIetaIeta < 0.00994
      	   && TMath::Max(chIso-chArea*rho,0.0) < 0.202
      	   && TMath::Max(nhIso-nhArea*rho,0.0) < 0.264 + 0.0148*pt+0.000017*pt*pt
      	   && TMath::Max(gIso-gArea*rho,  0.0) < 2.362 + 0.0047*pt )
      	return true;
      
            if ( endcap
           && !elevto 
      	   && hoe < 0.0213 
      	   && sigmaIetaIeta < 0.03000
      	   && TMath::Max(chIso-chArea*rho,0.0) < 0.034 
      	   && TMath::Max(nhIso-nhArea*rho,0.0) < 0.586 + 0.0163*pt + 0.000014*pt*pt
      	   && TMath::Max(gIso-gArea*rho,  0.0) < 2.617 + 0.0034*pt ) 
      	return true; 
      break;          

    default:
      printf("FIXME PhotonId llvvPhotonId::%i is unkown\n", IdLevel);
      return false;
      break;
      
    }    
    
    return false; 
  }

  float relIso(patUtils::GenericLepton& lep, double rho){
    
    int lid=lep.pdgId();
    float relIso = 0.0; 
      
    if(lid==13){

      float  chIso   = lep.mu.pfIsolationR04().sumChargedHadronPt;
      float  nhIso   = lep.mu.pfIsolationR04().sumNeutralHadronEt;
      float  gIso    = lep.mu.pfIsolationR04().sumPhotonEt;
      float  puchIso = lep.mu.pfIsolationR04().sumPUPt;
      
      relIso  = (chIso + TMath::Max(0.,nhIso+gIso-0.5*puchIso)) / lep.mu.pt();
      
    } else if (lid==11){ 
      
      //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
      float  chIso   = lep.el.pfIsolationVariables().sumChargedHadronPt;
      float  nhIso   = lep.el.pfIsolationVariables().sumNeutralHadronEt;
      float  gIso    = lep.el.pfIsolationVariables().sumPhotonEt;
        
      if (rho == 0) {
	float  puchIso = lep.el.pfIsolationVariables().sumPUPt; 
	relIso  = (chIso + TMath::Max(0.,nhIso+gIso-0.5*puchIso)) / lep.el.pt();
      }
      else {
	float effArea = utils::cmssw::getEffectiveArea(11,lep.el.superCluster()->eta());
	relIso  = (chIso + TMath::Max(0.,nhIso+gIso-rho*effArea)) / lep.el.pt();
      }
      
    }
    
    return relIso;

  }


  bool passIso (VersionedPatElectronSelector id, pat::Electron& el){
    // This assumes an object to be created ( *before the event loop* ):
    // VersionedPatElectronSelector loose_id("some_VID_tag_including_the_WP");
    return true; // Isolation is now embedded into the ID.
  }
  
  bool passIso(pat::Electron& el, int IsoLevel, int cutVersion, double rho  ){
         //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
	  float  chIso   = el.pfIsolationVariables().sumChargedHadronPt;
          float  nhIso   = el.pfIsolationVariables().sumNeutralHadronEt;
          float  gIso    = el.pfIsolationVariables().sumPhotonEt;

	  float  relIso = 0.0; 
	  
	  if (rho == 0) {
	    float  puchIso = el.pfIsolationVariables().sumPUPt; 
	    relIso  = (chIso + TMath::Max(0.,nhIso+gIso-0.5*puchIso)) / el.pt();
	  }
	  else {
	    float effArea = utils::cmssw::getEffectiveArea(11,el.superCluster()->eta());
	    relIso  = (chIso + TMath::Max(0.,nhIso+gIso-rho*effArea)) / el.pt();
	  }
	  
          bool barrel = (fabs(el.superCluster()->eta()) <= 1.479);
          bool endcap = (!barrel && fabs(el.superCluster()->eta()) < 2.5);
          // Spring15 selection 
          switch(cutVersion){
          case CutVersion::Spring15Cut25ns :
      	  // Spring15 selection, conditions: PU20 bx25
                switch(IsoLevel){
                     case llvvElecIso::Veto :
                        if( barrel && relIso < 0.126    ) return true;
                        if( endcap && relIso < 0.144    ) return true;
                        break;
      
                     case llvvElecIso::Loose :
                        if( barrel && relIso < 0.0893   ) return true;
                        if( endcap && relIso < 0.121    ) return true;
                        break;
      
                     case llvvElecIso::Medium :
                        if( barrel && relIso < 0.0766   ) return true;
                        if( endcap && relIso < 0.0678   ) return true;
                        break;
      
                     case llvvElecIso::Tight :
                        if( barrel && relIso < 0.0354   ) return true;
                        if( endcap && relIso < 0.0646   ) return true;
                        break;
      
                     default:
                        printf("FIXME ElectronIso llvvElectronIso::%i is unkown\n", IsoLevel);
                        return false;
                        break;
                }
             break;
      
          case CutVersion::Moriond17Cut :
          case CutVersion::ICHEP16Cut :
      	  // ICHEP16 or Moriond17 selection, conditions: PU20 bx25
               switch(IsoLevel){
                     case llvvElecIso::Veto :
                        if( barrel && relIso < 0.175    ) return true;
                        if( endcap && relIso < 0.159    ) return true;
                        break;
      
                     case llvvElecIso::Loose :
                        if( barrel && relIso < 0.0994   ) return true;
                        if( endcap && relIso < 0.107    ) return true;
                        break;
      
                     case llvvElecIso::Medium :
                        if( barrel && relIso < 0.0695   ) return true;
                        if( endcap && relIso < 0.0821   ) return true;
                        break;
      
                     case llvvElecIso::Tight :
                        if( barrel && relIso < 0.0588   ) return true;
                        if( endcap && relIso < 0.0571   ) return true;
                        break;
      
                     default:
                        printf("FIXME ElectronIso llvvElectronIso::%i is unkown\n", IsoLevel);
                        return false;
                        break;
                }
             break;

          default:
             printf("FIXME ElectronIsolation  CutVersion::%i is unkown\n", cutVersion);
             return false;
             break;

	  }

          return false;  
   }
  
  bool passIso(pat::Muon& mu, int IsoLevel, int cutVersion){
    //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
    float  chIso   = mu.pfIsolationR04().sumChargedHadronPt;
    float  nhIso   = mu.pfIsolationR04().sumNeutralHadronEt;
    float  gIso    = mu.pfIsolationR04().sumPhotonEt;
    float  puchIso = mu.pfIsolationR04().sumPUPt;
    float  relIso  = (chIso + TMath::Max(0.,nhIso+gIso-0.5*puchIso)) / mu.pt();
    float  chIso03   = mu.pfIsolationR03().sumChargedHadronPt;
    float  nhIso03   = mu.pfIsolationR03().sumNeutralHadronEt;
    float  gIso03    = mu.pfIsolationR03().sumPhotonEt;
    float  puchIso03 = mu.pfIsolationR03().sumPUPt;
    float  relIso03  = (chIso03 + TMath::Max(0.,nhIso03+gIso03-0.5*puchIso03)) / mu.pt();
    float  trkrelIso = mu.isolationR03().sumPt/mu.pt(); // no PU correction
    
    switch(cutVersion){
       case CutVersion::Spring15Cut25ns :
           switch(IsoLevel){
              case llvvMuonIso::Loose : 
                 if( relIso < 0.20 ) return true;
                 break;
               	 
              case llvvMuonIso::Tight :
                if( relIso < 0.15 ) return true;
                break;
             
              default:
                printf("FIXME MuonIso llvvMuonIso::%i is unkown\n", IsoLevel);
                return false;
                break;
           }
           break;
       case CutVersion::ICHEP16Cut :
           switch(IsoLevel){
              case llvvMuonIso::Loose : 
                 if( relIso < 0.20 && trkrelIso < 0.1) return true;
                 break;
               	 
              case llvvMuonIso::Tight :
                if( relIso < 0.15 && trkrelIso < 0.1) return true;
                break;
             
              default:
                printf("FIXME MuonIso llvvMuonIso::%i is unkown\n", IsoLevel);
                return false;
                break;
           }
           break;

       case CutVersion::Moriond17Cut :
           switch(IsoLevel){
              case llvvMuonIso::Loose : 
                 if( relIso < 0.25 ) return true;
                 break;
               	 
              case llvvMuonIso::Tight :
                if( relIso < 0.15 ) return true;
                break;
              case llvvMuonIso::H4lWP :
                if( relIso03 < 0.35 ) return true;
                break;
              default:
                printf("FIXME MuonIso llvvMuonIso::%i is unkown\n", IsoLevel);
                return false;
                break;
           }
           break;

       default:
           printf("FIXME MuonIsolation  CutVersion::%i is unkown\n", cutVersion);
           return false;
           break;
    }

    return false;          
  }

/*
  bool passPhotonTrigger(fwlite::Event &ev, float &triggerThreshold, float &triggerPrescale, float& triggerThresholdHigh ){
    edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
    if( !tr.isValid() ) return false;

    bool hasPhotonTrigger(false);

    triggerPrescale = 1.0; 
    triggerThreshold = 0.0;
    triggerThresholdHigh = 999999;
    

    std::string successfulPath="";
    if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon300_*")){
      hasPhotonTrigger=true;
      triggerThreshold=300;
      triggerThresholdHigh=999999;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon250_*")){
      hasPhotonTrigger=true;
      triggerThreshold=250;
      triggerThresholdHigh=300;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon165_R9Id90_HE10_IsoM_*")){
      hasPhotonTrigger=true;
     triggerThreshold=165;
     triggerThresholdHigh=250;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon120_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=120;
      triggerThresholdHigh=165;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon90_R9Id90_HE10_IsoM_*")){ // HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=90;
      triggerThresholdHigh=120;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon75_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=75;
      triggerThresholdHigh=90;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon50_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=50;
      triggerThresholdHigh=75;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon36_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=36;
      triggerThresholdHigh=50;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon30_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){                       
      hasPhotonTrigger=true;                                                                                                                     
      triggerThreshold=30;    
      triggerThresholdHigh=36;
    } 
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon22_R9Id90_HE10_IsoM_*")){ //HE10_Iso40_EBOnly_*")){
      hasPhotonTrigger=true;
      triggerThreshold=22;
      triggerThresholdHigh=30;
    }
      
    if(successfulPath!=""){ //get the prescale associated to it
      fwlite::Handle< pat::PackedTriggerPrescales > prescalesHandle;
      prescalesHandle.getByLabel(ev, "patTrigger");
      pat::PackedTriggerPrescales prescales = *prescalesHandle;
      const edm::TriggerResults& trResults =  prescales.triggerResults();
      prescales.setTriggerNames( ev.triggerNames(trResults) );
      triggerPrescale = prescales.getPrescaleForName(successfulPath);
    }

    return hasPhotonTrigger; 
  }

  
  bool passVBFPhotonTrigger(fwlite::Event &ev, float &triggerThreshold,
			    float &triggerPrescale, float &triggerThresholdHigh ){
    edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
    if( !tr.isValid() ) return false;

    bool hasPhotonTrigger(false);
    // float triggerPrescale(1.0); 
    // float triggerThreshold(0);
    triggerPrescale = 1.0; 
    triggerThreshold = 0.0;
    triggerThresholdHigh = 999999; 

    std::string successfulPath="";
    if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon300_*")){
      hasPhotonTrigger=true;
      triggerThreshold=300;
      triggerThresholdHigh=999999;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon250_*")){
      hasPhotonTrigger=true;
      triggerThreshold=250;
      triggerThresholdHigh=300; 
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF_*")){
      hasPhotonTrigger=true;
      triggerThreshold=120;
      triggerThresholdHigh=165;
    }
    else if( utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF_*")){
      hasPhotonTrigger=true;
      triggerThreshold=90;
      triggerThresholdHigh=120;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_*")){
      hasPhotonTrigger=true;
      triggerThreshold=75;
      triggerThresholdHigh=90;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_*")){
      hasPhotonTrigger=true;
      triggerThreshold=50;
      triggerThresholdHigh=75;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_*")){
      hasPhotonTrigger=true;
      triggerThreshold=36;
      triggerThresholdHigh=50;
    }
    else if(utils::passTriggerPatternsAndGetName(tr, successfulPath, "HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_*")){
      hasPhotonTrigger=true;
      triggerThreshold=22;
      triggerThresholdHigh=36;
    }
      
    if(successfulPath!=""){ //get the prescale associated to it
      fwlite::Handle< pat::PackedTriggerPrescales > prescalesHandle;
      prescalesHandle.getByLabel(ev, "patTrigger");
      pat::PackedTriggerPrescales prescales = *prescalesHandle;
      const edm::TriggerResults& trResults =  prescales.triggerResults();
      prescales.setTriggerNames( ev.triggerNames(trResults) );
      triggerPrescale = prescales.getPrescaleForName(successfulPath);
    }

    return hasPhotonTrigger; 
  }
*/

  bool passPFJetID(std::string label,
                  pat::Jet jet){

    bool passID = false;

    float rawJetEn(jet.correctedJet("Uncorrected").energy() );
    // Note: All fractions are calculated with the raw/uncorrected energy of the jet (only then they add up to unity). So the PF JetID has to be applied before the jet energy corrections. 

    float nhf( (jet.neutralHadronEnergy() + jet.HFHadronEnergy())/rawJetEn );
    float nef( jet.neutralEmEnergy()/rawJetEn );
    float cef( jet.chargedEmEnergy()/rawJetEn );
    float chf( jet.chargedHadronEnergy()/rawJetEn );
    float nch    = jet.chargedMultiplicity();
    float nconst = jet.numberOfDaughters();
    //float muf( jet.muonEnergy()/rawJetEn);
    float NumNeutralParticles = jet.neutralMultiplicity();

    //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016 (27-dec-2016)
    if (label == "Loose"){
      if( fabs(jet.eta()) <= 2.7) passID = ( (nhf<0.99  && nef<0.99 && nconst>1) && ( fabs(jet.eta())>2.4 || (fabs(jet.eta()) <= 2.4 && chf>0 && nch>0 && cef<0.99) ) );
      if( fabs(jet.eta()) > 2.7 && fabs(jet.eta()) <= 3.0) passID = ( nef>0.01 && nhf<0.98 && NumNeutralParticles > 2);
      if( fabs(jet.eta()) > 3.0) passID = (nef<0.90 && NumNeutralParticles > 10);
    }
    if (label == "Tight"){ 
      if( fabs(jet.eta()) <= 2.7) passID = ( (nhf<0.90  && nef<0.90 && nconst>1) && ( fabs(jet.eta())>2.4 || (fabs(jet.eta()) <= 2.4 && chf>0 && nch>0 && cef<0.99) ) );
      if( fabs(jet.eta()) > 2.7 && fabs(jet.eta()) <= 3.0) passID = ( nef>0.01 && nhf<0.98 && NumNeutralParticles > 2);
      if( fabs(jet.eta()) > 3.0) passID = (nef<0.90 && NumNeutralParticles > 10);
    }
    return passID;
  }


  bool passPUJetID(pat::Jet j){

    double jpt = j.pt();
    double jeta = j.eta();

    //Recommendation of HZZ :https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2015#Jets
    //FIXME recommendation of HZZ are oudated (run 1)--->update them
    //FIXME Progress: Partially done (removed the cuts at eta>2.75) but have to wait for improvements from JetMET group
    //FIXME In miniAOD V2, the variable used is bugged ! Please avoir using it!
    float jpumva=0.;
    jpumva=j.userFloat("pileupJetId:fullDiscriminant");

    bool passPU = true;
    if(jpt>20){
      if(jeta>3.){
        //if(jpumva<=-0.45)passPU=false;
      }else if(jeta>2.75){
        //if(jpumva<=-0.55)passPU=false;
      }else if(jeta>2.5){
        if(jpumva<=-0.6)passPU=false;
      }else if(jpumva<=-0.63)passPU=false;
    }else{ //pt<20 : in the 2l2nu analysis, this means 15<pt<20
      if(jeta>3.){
        //if(jpumva<=-0.95)passPU=false;
      }else if(jeta>2.75){
        //if(jpumva<=-0.94)passPU=false;
      }else if(jeta>2.5){
        if(jpumva<=-0.96)passPU=false;
      }else if(jpumva<=-0.95)passPU=false;
    }
    return passPU;
  }

  bool exclusiveDataEventFilter(const double& run, const bool& isMC, const bool& isPromptReco)
  {
    bool passExclusiveDataEventFilter(false);
    
    if(isMC)
      passExclusiveDataEventFilter=true;
    else
      {
        bool isForPromptReco(run<251162 || run>251562);
        if(isPromptReco)  // Prompt reco keeps events outside of that range
          {
            if(isForPromptReco) passExclusiveDataEventFilter=true;
            else                passExclusiveDataEventFilter=false;
          }
        else // 17Jul15 ReReco keeps event inside that range
          {
            if(isForPromptReco) passExclusiveDataEventFilter=false;
            else                passExclusiveDataEventFilter=true;
          }
      }
    return passExclusiveDataEventFilter;
  }


  bool outInOnly(const reco::Muon &mu)  {
      const reco::Track &tk = *mu.innerTrack();
      return tk.algoMask().count() == 1 && tk.isAlgoInMask(reco::Track::muonSeededStepOutIn);
  }
  bool preselection(const reco::Muon &mu, bool selectClones_)  { 
      return mu.isGlobalMuon() && (!selectClones_ || outInOnly(mu));
  }
  bool tighterId(const reco::Muon &mu)  { 
      return muon::isMediumMuon(mu) && mu.numberOfMatchedStations() >= 2; 
  }
  bool tightGlobal(const reco::Muon &mu)  {
      return (mu.globalTrack()->hitPattern().muonStationsWithValidHits() >= 3 && mu.globalTrack()->normalizedChi2() <= 20);
  }
  bool safeId(const reco::Muon &mu)  { 
      if (mu.muonBestTrack()->ptError() > 0.2 * mu.muonBestTrack()->pt()) { return false; }
      return mu.numberOfMatchedStations() >= 1 || tightGlobal(mu);
  }
  bool partnerId(const reco::Muon &mu)  {
            return mu.pt() >= 10 && mu.numberOfMatchedStations() >= 1;
  }


std::pair<double, double> scaleVariation(const fwlite::Event& ev){
	//std::cout << " " << std::endl;
	//std::cout << "STARTING SCALE-VARIATION Estimation" << std::endl;
        fwlite::Handle<LHEEventProduct> lheEPHandle;
        lheEPHandle.getByLabel( ev, "externalLHEProducer");
        double scaleUp = 1.;
        double scaleDw = 1.;
        bool check_in = false;
        std::vector<int> idVect;
	std::vector<int>::iterator it;
        if( lheEPHandle.isValid() ){
                for (unsigned int i=0; i<lheEPHandle->weights().size(); i++) {
                        std::string::size_type sz;
                        double id = std::stod( lheEPHandle->weights()[i].id, &sz);
                        idVect.push_back( id );
                }
		for( unsigned int k=1001; k<1010; k++){
			it = find( idVect.begin(), idVect.end(), k);
			if( it != idVect.end() ){ check_in = true; }
			else{ check_in = false; }
		}
		if( check_in ){
                	for (unsigned int i=0; i<lheEPHandle->weights().size(); i++) {
                        	if( lheEPHandle->weights()[i].id != "1001" || lheEPHandle->weights()[i].id != "1006" || lheEPHandle->weights()[i].id != "1008" ){
                                	double local_weight = 0;
                                	local_weight = ( lheEPHandle->weights()[i].wgt / lheEPHandle->originalXWGTUP() );
					//std::cout << "Local weight: " << local_weight << std::endl; 
                                	scaleUp = std::max(scaleUp, local_weight);
                                	scaleDw = std::min(scaleDw, local_weight);
                        	}
                 	}
			
		 } else { scaleUp = 1.;  scaleDw = 1.;}
         }
	 //std::cout << "ScaleUp Value: " << scaleUp << "; ScaleDwn Value: " << scaleDw << std::endl;
	 return std::make_pair(scaleUp, scaleDw);
}

double pdfVariation(const fwlite::Event& ev){
	//std::cout << "  " << std::endl;
	//std::cout << "STARTING PDF Estimation" << std::endl;
        fwlite::Handle<LHEEventProduct> lheEPHandle;
        lheEPHandle.getByLabel( ev, "externalLHEProducer");
        int N = 0;
        double pdfVar = 0;
        double sum = 0;
        bool check_in = false;
        std::vector<int> idVect;
        std::vector<int>::iterator it;
        if( lheEPHandle.isValid() ){
                for (unsigned int i=0; i<lheEPHandle->weights().size(); i++) {
                        std::string::size_type sz;
                        double id = std::stod( lheEPHandle->weights()[i].id, &sz);
                        idVect.push_back( id );
                }
                for( unsigned int k=2001; k<2101; k++){
                        it = find( idVect.begin(), idVect.end(), k);
                        if( it != idVect.end() ){ check_in = true; }
                        else{ check_in = false; }
                }
                if( check_in ){
                	for (unsigned int i=0; i<lheEPHandle->weights().size(); i++) {
                        	std::string::size_type sz;
                        	double id = std::stod( lheEPHandle->weights()[i].id, &sz); 
                                if( id<2001 || id>2100 ) continue;	
				//std::cout << "Weight: " << lheEPHandle->weights()[i].wgt << "; Nominal Weight: " << lheEPHandle->originalXWGTUP() << std::endl;
                                sum += std::pow( (lheEPHandle->weights()[i].wgt / lheEPHandle->originalXWGTUP() - 1 ), 2);
                                N++; 
                                	
                        }	
			pdfVar = 1+ std::sqrt( sum/ ( N -1 ) ); //+1 variation
        	} else { pdfVar = 1.; }
	}
        return pdfVar;
}

double alphaVariation(const fwlite::Event& ev){
        //std::cout << "  " << std::endl;
        //std::cout << "STARTING ALPHA Estimation" << std::endl;
        fwlite::Handle<LHEEventProduct> lheEPHandle;
        lheEPHandle.getByLabel( ev, "externalLHEProducer");
        double alphaVar = 0;
        double local_alpha_one = 0;
        double local_alpha_two = 0;
	std::vector<int> idVect;
        std::vector<int>::iterator itone, ittwo;
        if( lheEPHandle.isValid() ){
                for (unsigned int i=0; i<lheEPHandle->weights().size(); i++) {
                        std::string::size_type sz;
                        double id = std::stod( lheEPHandle->weights()[i].id, &sz);		
			idVect.push_back( id );
		}
		itone = find( idVect.begin(), idVect.end(), 2101);
		ittwo = find( idVect.begin(), idVect.end(), 2102);
		if( itone != idVect.end() && ittwo != idVect.end() ){
			for (unsigned int i=0; i<lheEPHandle->weights().size(); i++) {
                        	std::string::size_type sz;
                        	double id = std::stod( lheEPHandle->weights()[i].id, &sz);	
                        	if( ( id == 2101 ) ){
                                	local_alpha_one = ( lheEPHandle->weights()[i].wgt / lheEPHandle->originalXWGTUP() );
                        	} else if( ( id == 2102 ) ){
                                	local_alpha_two = ( lheEPHandle->weights()[i].wgt / lheEPHandle->originalXWGTUP() );
                        	}
				//std::cout << "alpha one: " << local_alpha_one << "; alpha two: " << local_alpha_two << "; Nominal: " << lheEPHandle->originalXWGTUP() << std::endl;
                	}
			alphaVar = 1+std::sqrt(0.75)*std::abs( local_alpha_one - local_alpha_two )*0.5; //+1 variation
		} else { alphaVar = 1; }
        }
	//std::cout << "Alpha Uncertainties: " << alphaVar << std::endl;
        return alphaVar;
}
  
                                                                                                                                              
  double getHTScaleFactor(TString dtag, double lheHt)
  {
    // NNLO per-event weights as a function of generator level HT
    // (from https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2015#MC_and_data_samples )
    // For DY-5to50, only LO is available.   

    // Please look below for the code snippet necessary to run this.

    double htScaleFactor(1.0);              
    if(dtag.Contains("WJetsToLNu"))
      {
        // NNLO
        // Valid for:                                                       xsection [pb]
        // /WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8                  50690
        // /WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8       1345
        // /WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8        359.7
        // /WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8         48.91
        // /WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8         18.77
        if     (lheHt<100              ) htScaleFactor=0.8520862372;
        else if(lheHt>=100 && lheHt<200) htScaleFactor=0.1352710705;
        else if(lheHt>=200 && lheHt<400) htScaleFactor=0.076142149 ;                                                
        else if(lheHt>=400 && lheHt<600) htScaleFactor=0.0326980819;
        else if(lheHt>=600             ) htScaleFactor=0.0213743732;
      }                                     
    else if(dtag.Contains("DYJetsToLL_M-50"))
      {                                     
        // NNLO                             
        // Valid for:                                                               xsection [pb]
        // /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8                    4895
        // /DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8         139.4
        // /DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8          42.75
        // /DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8           5.497
        // /DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8           2.21
        if     (lheHt<100              ) htScaleFactor=0.6655715203;
        else if(lheHt>=100 && lheHt<200) htScaleFactor=0.0575124298;
        else if(lheHt>=200 && lheHt<400) htScaleFactor=0.049972089 ;
        else if(lheHt>=400 && lheHt<600) htScaleFactor=0.0062770613;
        else if(lheHt>=600             ) htScaleFactor=0.00271213  ;
      }                                     
    else if(dtag.Contains("DYJetsToLL_M-5to50"))
      {
        // LO                               
        // Valid for:                                                               xsection [pb]
        // /DYJetsToLL_M-5to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8                  71310
        // /DYJetsToLL_M-5to50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8        224.2
        // /DYJetsToLL_M-5to50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8         37.2
        // /DYJetsToLL_M-5to50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8          3.581
        // /DYJetsToLL_M-5to50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8          1.124
        if     (lheHt<100              ) htScaleFactor=7.5826225134;
        else if(lheHt>=100 && lheHt<200) htScaleFactor=0.2149472503;
        else if(lheHt>=200 && lheHt<400) htScaleFactor=0.0365903335;
        else if(lheHt>=400 && lheHt<600) htScaleFactor=0.0035837837;
        else if(lheHt>=600             ) htScaleFactor=0.0011156801;
      }                                     

    htScaleFactor /= 1000; // The scale factors are derived for 1/fb, whereas we normalize in picobarns
    return htScaleFactor;                   

    // In order to use this stitching, you need to run a full set of HT-binned samples plus the inclusive one.
    // You should *NOT* have the nhepup cut that was used for the stitching of jet-binned samples. Either one or the other.
    
    // HT-binned samples stitching: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2015#MC_and_data_samples
    /*
    double weightGen(1.0);
    // do not correct for negative weights: these are LO samples
      
    bool isV0JetsMC   (isMC && (dtag.Contains ("DYJetsToLL") || dtag.Contains ("WJets")));
    if(isV0JetsMC)
      {
        // access generator level HT               
        fwlite::Handle<LHEEventProduct> lheEventProduct;
        lheEventProduct.getByLabel(ev, "externalLHEProducer");
        //edm::Handle<LHEEventProduct> lheEventProduct;
        //ev.getByLabel( 'externalLHEProducer', lheEventProduct);
        const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup(); 
        std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
        double lheHt = 0.;
        size_t numParticles = lheParticles.size();
        for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
          int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
          int status = lheEvent.ISTUP[idxParticle];
          if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) { // quarks and gluons
            lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
          }                                        
        }
        if(debug) cout << "Sample: " << dtag << ", lheHt: " << lheHt << ", scale factor from spreadsheet: " << patUtils::getHTScaleFactor(dtag, lheHt) << endl;
        weightGen *=   patUtils::getHTScaleFactor(dtag, lheHt);
      }         

      weight *= weightGen;
    */
  }                                               

}
