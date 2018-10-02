#ifndef PROCESSINGMUONS_H
#define PROCESSINGMUONS_H

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "UserCode/NtuplerAnalyzer/interface/PatUtils.h"
//#include "UserCode/NtuplerAnalyzer/interface/recordFuncs.h"


int processMuons_ID_ISO_Kinematics(pat::MuonCollection& muons, reco::Vertex goodPV, double weight, // input
	patUtils::llvvMuonId::MuonId mu_ID, patUtils::llvvMuonId::MuonId veto_mu_ID,               // config/cuts
	patUtils::llvvMuonIso::MuonIso mu_ISO, patUtils::llvvMuonIso::MuonIso veto_mu_ISO,
	double pt_cut, double eta_cut, double veto_pt_cut, double veto_eta_cut,
	pat::MuonCollection& selMuons, LorentzVector& muDiff, unsigned int& nVetoMu, unsigned int& nVetoMu_all,               // output
	bool record, bool debug); // more output

int processMuons_MatchHLT(
	pat::MuonCollection& muons,
	vector<pat::TriggerObjectStandAlone>& trig_objs,    // input: trigger objects to match against (so, these should match HLT of interest)
	float min_dR,
	pat::MuonCollection& muons_matched
	);

bool processMuon_matchesHLTs(
	pat::Muon& muon,
	vector<pat::TriggerObjectStandAlone>& trig_objs,    // input: trigger objects to match against (so, these should match HLT of interest)
	float min_dR
	);

float relIso(pat::Muon& lep, double rho);

#endif /* PROCESSINGMUONS_H */


