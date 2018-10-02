#include "UserCode/NtuplerAnalyzer/interface/ProcessingTaus.h"
#include "UserCode/NtuplerAnalyzer/interface/recordFuncs.h"



int processTaus_ID(pat::TauCollection& taus, double weight, // input
	string& tauID_decayMode, //string& tauID,               // config/cuts
	string& tauID_IsoMuons,  string& tauID_IsoElectrons,
	pat::TauCollection& selTaus,                          // output
	bool record, bool debug) // more output

{
for (unsigned int count_ided_taus = 0, n = 0; n < taus.size(); ++n)
	{
	pat::Tau& tau = taus[n];

	if (record)
		{
		fill_2d(string("control_tau_slimmedtaus_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
		fill_1d(string("control_tau_slimmedtaus_phi"), 128, -3.2, 3.2, tau.phi(), weight);
		}

	// ---------- IDs
			
	// if(!tau.isPFTau()) continue; // Only PFTaus // It should be false for slimmedTaus
	// if(tau.emFraction() >=2.) continue;
			
	// Discriminators from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
	// "The tau passes the discriminator if pat::Tau::tauID("name") returns a value of 0.5 or greater"
	//if(tau.tauID("decayModeFindingNewDMs")<0.5) continue; // High pt tau. Otherwise, OldDMs
	// decayModeFindingOldDMs is synonim for decayModeFinding
	// Anyways, the collection of taus from miniAOD should be already afer decayModeFinding cut (the tag - Old or New - is unspecified in the twiki, though).
	// Consequently, there might be a small bias due to events that are cut by the OldDM and would not be cut by the NewDM

	if (tau.tauID(tauID_decayMode) < 0.5) continue;
	//if (tau.tauID(tauID) < 0.5)           continue; // See whether to us the new byMediumPileupWeightedIsolation3Hits that is available only for dynamic strip reconstruction (default in CMSSW_7_4_14)
	if (tau.tauID(tauID_IsoMuons) < 0.5)  continue; // Medium working point not usable. Available values: Loose, Tight
	if (tau.tauID(tauID_IsoElectrons) < 0.5) continue;

	selTaus.push_back(tau);

	//if (tau.tauID ("againstElectronMediumMVA5")                  <0.5) continue; // Tight working point not usable. Avaiable values: VLoose, Loose, Medium
	// if (tau.tauID ("againstElectronMediumMVA6")                  <0.5) continue;
	// if (tau.tauID ("byTightCombinedIsolationDeltaBetaCorr3Hits")<0.5) continue;
	//if (tau.tauID ("byMediumIsolationMVArun2v1DBoldDMwLT")<0.5) continue;
	// maybe use MVA-based discriminators:
	// byVLooseIsolationMVArun2v1DBoldDMwLT
	// byLooseIsolationMVArun2v1DBoldDMwLT
	// byMediumIsolationMVArun2v1DBoldDMwLT
	// byTightIsolationMVArun2v1DBoldDMwLT
	// byVTightIsolationMVArun2v1DBoldDMwLT

	// or the new discriminators:
	// byLooseCombinedIsolationDeltaBetaCorr3HitsdR03
	// byMediumCombinedIsolationDeltaBetaCorr3HitsdR03
	// byTightCombinedIsolationDeltaBetaCorr3HitsdR03
	// -- recommended for multi-object final states (ttH, H->tau-tau)
	// -- not found in noHLT TTbar

	if (record)
		{
		fill_2d(string("control_tau_idedtaus_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
		fill_1d(string("control_tau_idedtaus_phi"), 128, -3.2, 3.2, tau.phi(), weight);
		}

	// Pixel hits cut (will be available out of the box in new MINIAOD production)
	//{
	//int nChHadPixelHits = 0;
	//reco::CandidatePtrVector chCands = tau.signalChargedHadrCands();
	//for(reco::CandidatePtrVector::const_iterator iter = chCands.begin(); iter != chCands.end(); iter++)
	//	{
	//	pat::PackedCandidate const* packedCand = dynamic_cast<pat::PackedCandidate const*>(iter->get());
	//	int pixelHits = packedCand->numberOfPixelHits();
	//	if(pixelHits > nChHadPixelHits) nChHadPixelHits = pixelHits;
	//	}
	//if(nChHadPixelHits==0) continue;
	//}

	/*
	fill_pt_e( string("all_taus_ided_pt"), tau.pt(), weight);
	if (count_ided_taus<1)
		{
		fill_pt_e( string("top1pt_taus_ided_pt"), tau.pt(), weight);
		count_ided_taus += 1;
		}
	*/
	}

std::sort (selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);

return 0;
}
int processTaus_ID_ISO(pat::TauCollection& taus, double weight, // input
	string& tauID_decayMode, string& tauID,               // config/cuts
	string& tauID_IsoMuons,  string& tauID_IsoElectrons,
	pat::TauCollection& selTaus,                          // output
	bool record, bool debug) // more output

{
for (unsigned int count_ided_taus = 0, n = 0; n < taus.size(); ++n)
	{
	pat::Tau& tau = taus[n];

	if (record)
		{
		fill_2d(string("control_tau_slimmedtaus_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
		fill_1d(string("control_tau_slimmedtaus_phi"), 128, -3.2, 3.2, tau.phi(), weight);
		}

	// ---------- IDs
			
	// if(!tau.isPFTau()) continue; // Only PFTaus // It should be false for slimmedTaus
	// if(tau.emFraction() >=2.) continue;
			
	// Discriminators from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
	// "The tau passes the discriminator if pat::Tau::tauID("name") returns a value of 0.5 or greater"
	//if(tau.tauID("decayModeFindingNewDMs")<0.5) continue; // High pt tau. Otherwise, OldDMs
	// decayModeFindingOldDMs is synonim for decayModeFinding
	// Anyways, the collection of taus from miniAOD should be already afer decayModeFinding cut (the tag - Old or New - is unspecified in the twiki, though).
	// Consequently, there might be a small bias due to events that are cut by the OldDM and would not be cut by the NewDM

	if (tau.tauID(tauID_decayMode) < 0.5) continue;
	if (tau.tauID(tauID) < 0.5)           continue; // See whether to us the new byMediumPileupWeightedIsolation3Hits that is available only for dynamic strip reconstruction (default in CMSSW_7_4_14)
	if (tau.tauID(tauID_IsoMuons) < 0.5)  continue; // Medium working point not usable. Available values: Loose, Tight
	if (tau.tauID(tauID_IsoElectrons) < 0.5) continue;

	selTaus.push_back(tau);

	//if (tau.tauID ("againstElectronMediumMVA5")                  <0.5) continue; // Tight working point not usable. Avaiable values: VLoose, Loose, Medium
	// if (tau.tauID ("againstElectronMediumMVA6")                  <0.5) continue;
	// if (tau.tauID ("byTightCombinedIsolationDeltaBetaCorr3Hits")<0.5) continue;
	//if (tau.tauID ("byMediumIsolationMVArun2v1DBoldDMwLT")<0.5) continue;
	// maybe use MVA-based discriminators:
	// byVLooseIsolationMVArun2v1DBoldDMwLT
	// byLooseIsolationMVArun2v1DBoldDMwLT
	// byMediumIsolationMVArun2v1DBoldDMwLT
	// byTightIsolationMVArun2v1DBoldDMwLT
	// byVTightIsolationMVArun2v1DBoldDMwLT

	// or the new discriminators:
	// byLooseCombinedIsolationDeltaBetaCorr3HitsdR03
	// byMediumCombinedIsolationDeltaBetaCorr3HitsdR03
	// byTightCombinedIsolationDeltaBetaCorr3HitsdR03
	// -- recommended for multi-object final states (ttH, H->tau-tau)
	// -- not found in noHLT TTbar

	if (record)
		{
		fill_2d(string("control_tau_idedtaus_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
		fill_1d(string("control_tau_idedtaus_phi"), 128, -3.2, 3.2, tau.phi(), weight);
		}

	// Pixel hits cut (will be available out of the box in new MINIAOD production)
	//{
	//int nChHadPixelHits = 0;
	//reco::CandidatePtrVector chCands = tau.signalChargedHadrCands();
	//for(reco::CandidatePtrVector::const_iterator iter = chCands.begin(); iter != chCands.end(); iter++)
	//	{
	//	pat::PackedCandidate const* packedCand = dynamic_cast<pat::PackedCandidate const*>(iter->get());
	//	int pixelHits = packedCand->numberOfPixelHits();
	//	if(pixelHits > nChHadPixelHits) nChHadPixelHits = pixelHits;
	//	}
	//if(nChHadPixelHits==0) continue;
	//}

	/*
	fill_pt_e( string("all_taus_ided_pt"), tau.pt(), weight);
	if (count_ided_taus<1)
		{
		fill_pt_e( string("top1pt_taus_ided_pt"), tau.pt(), weight);
		count_ided_taus += 1;
		}
	*/
	}

std::sort (selTaus.begin(), selTaus.end(), utils::sort_CandidatesByPt);

return 0;
}



int processTaus_Kinematics(pat::TauCollection& taus,          // input
	double weight,
	double pt_cut, double eta_cut,
	pat::TauCollection& selTaus,                          // output
	bool record, bool debug) // more output

{

for (unsigned int count_ided_taus = 0, n = 0; n < taus.size(); ++n)
	{
	pat::Tau& tau = taus[n];

	// --------- Tau Kinematics
	if (tau.pt() < pt_cut || fabs (tau.eta()) > eta_cut) continue;

	// add leading track Pt as userFloat
	//reco::TrackRef pat::Tau::leadTrack	(		)	const
	reco::TrackRef lead_track = tau.leadTrack();
	if (lead_track.isNonnull())
		{
		// have to do separate if here!
		if (lead_track->outerOk())
			tau.addUserFloat("leading_track_pt", lead_track->outerPt());
		else
			tau.addUserFloat("leading_track_pt", -1);
		}
	else
		tau.addUserFloat("leading_track_pt", -1);
	// also add these
	// leadChargedHadrCand(), leadNeutralCand(), leadCand(), 
	// they should always be present in MINIAOD taus

	tau.addUserFloat("leadChargedHadrCand_pt", tau.leadChargedHadrCand().isNonnull() ? tau.leadChargedHadrCand()->pt() : -1);
	tau.addUserFloat("leadNeutralCand_pt",     tau.leadNeutralCand().isNonnull() ? tau.leadNeutralCand()->pt() : -1);
	tau.addUserFloat("leadCand_pt",            tau.leadCand().isNonnull() ? tau.leadCand()->pt() : -1);
	/*
	tau.addUserFloat("leadChargedHadrCand_pt", -1);
	tau.addUserFloat("leadNeutralCand_pt", -1);
	tau.addUserFloat("leadCand_pt", -1);
	*/

	selTaus.push_back(tau);
	if (record)
		{
		fill_2d(string("control_tau_selTaus_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
		fill_1d(string("control_tau_selTaus_phi"), 128, -3.2, 3.2, tau.phi(), weight);
		}

	}
return 0;
}
