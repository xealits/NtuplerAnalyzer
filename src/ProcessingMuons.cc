#include "UserCode/NtuplerAnalyzer/interface/ProcessingMuons.h"
//#include "UserCode/NtuplerAnalyzer/interface/recordFuncs.h"



int processMuons_ID_ISO_Kinematics(
	pat::MuonCollection& muons, reco::Vertex goodPV, double weight, // input
	patUtils::llvvMuonId::MuonId   mu_ID, patUtils::llvvMuonId::MuonId veto_mu_ID,                          // config/cuts
	patUtils::llvvMuonIso::MuonIso mu_ISO, patUtils::llvvMuonIso::MuonIso veto_mu_ISO,
	double pt_cut, double eta_cut, double veto_pt_cut, double veto_eta_cut,
	// output
	pat::MuonCollection& selMuons, pat::MuonCollection& selMuons_allIso,
	LorentzVector& muDiff,
	//unsigned int& nVetoMu_IsoImp, // muons have impact embedded into ID, no need to account for it
	unsigned int& nVetoMu_Iso, unsigned int& nVetoMu_all, // all the veto counters
	bool record, bool debug) // more output

{
//LorentzVector muDiff(0., 0., 0., 0.);
// std::vector<patUtils::GenericLepton> selLeptons;
//pat::MuonCollection selMuons;
// unsigned int count_idiso_muons = 0;

for(unsigned int count_idiso_muons = 0, n=0; n<muons.size (); ++n)
	{
	// patUtils::GenericLepton& lepton = leptons[n];
	pat::Muon& muon = muons[n];

	bool
		passKin(true),     passId(true),     passIso(true),
		passVetoKin(true), passVetoId(true), passVetoIso(true);

	int lid(muon.pdgId()); // should always be 13

	//apply muon corrections
	/* no lepcorrs in 13.6
	if(abs(lid) == 13 && muCor)
		{
		float qter;
		TLorentzVector p4(muon.px(), muon.py(), muon.pz(), muon.energy());
		// old corrections:
		// muCor->applyPtCorrection(p4, lid < 0 ? -1 : 1);
		// if(isMC) muCor->applyPtSmearing(p4, lid < 0 ? -1 : 1, false);
		// roch-cor (rochcor) corrections:
		if (isMC) muCor->momcor_mc  (p4, lid<0 ? -1 :1, 0, qter);
		else muCor->momcor_data(p4, lid<0 ? -1 :1, 0, qter);
		muDiff -= muon.p4();
		muon.setP4(LorentzVector(p4.Px(), p4.Py(), p4.Pz(), p4.E()));
		muDiff += muon.p4();
		}
	*/


	//no need for charge info any longer
	//lid = abs(lid);
	//TString lepStr(lid == 13 ? "mu" : "e");
			
	// no need to mess with photon ID // //veto nearby photon (loose electrons are many times photons...)
	// no need to mess with photon ID // double minDRlg(9999.);
	// no need to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	// no need to mess with photon ID //   minDRlg=TMath::Min(minDRlg,deltaR(leptons[ilep].p4(),selPhotons[ipho].p4()));
	// no need to mess with photon ID // if(minDRlg<0.1) continue;

	//Cut based identification

	// ------------------------- muon IDs
	passId     = patUtils::passId(muon, goodPV, mu_ID,      patUtils::CutVersion::Legacy2016_07Aug17Jul);
	passVetoId = patUtils::passId(muon, goodPV, veto_mu_ID, patUtils::CutVersion::Legacy2016_07Aug17Jul);

	// ------------------------- muon isolation
	passIso     = patUtils::passIso(muon, mu_ISO,      patUtils::CutVersion::Legacy2016_07Aug17Jul);
	passVetoIso = patUtils::passIso(muon, veto_mu_ISO, patUtils::CutVersion::Legacy2016_07Aug17Jul);


	// ---------------------------- Muon Kinematics
	double leta( muon.eta());

	// ---------------------- Main muon kin
	if(muon.pt() < pt_cut)   passKin = false;
	// if(leta > 2.1)        passKin = false;
	//if(muon.pt() < 26.)   passKin = false;
	if(leta > eta_cut)        passKin = false;

	// ---------------------- Veto muon kin
	//if (muon.pt () < 20)  passVetoKin = false;
	// if (leta > 2.1)       passVetoKin = false;
	if (muon.pt() < veto_pt_cut)   passVetoKin = false;
	if (leta > veto_eta_cut)       passVetoKin = false;

	//if     (passKin     && passId     && passIso)
	if (passKin && passId) // all iso muons for anti-iso QCD region
		{
		selMuons_allIso.push_back(muon);
		if (passIso)
			selMuons.push_back(muon);
		}
	else
		{
		//if(passVetoKin && passVetoId && passVetoIso && passImp) nVetoMu_IsoImp++;
		if(passVetoKin && passVetoId && passVetoIso) nVetoMu_Iso++;
		if(passVetoKin && passVetoId) nVetoMu_all++;
		}
	}


// TODO: there should be no need to sort selected muons here again -- they are in order of Pt
std::sort(selMuons.begin(),   selMuons.end(),   utils::sort_CandidatesByPt);

return 0;
}





// old compatibility call

int processMuons_ID_ISO_Kinematics(
	pat::MuonCollection& muons, reco::Vertex goodPV, double weight, // input
	patUtils::llvvMuonId::MuonId   mu_ID, patUtils::llvvMuonId::MuonId veto_mu_ID,                          // config/cuts
	patUtils::llvvMuonIso::MuonIso mu_ISO, patUtils::llvvMuonIso::MuonIso veto_mu_ISO,
	double pt_cut, double eta_cut, double veto_pt_cut, double veto_eta_cut,
	// output
	pat::MuonCollection& selMuons,
	LorentzVector& muDiff,
	//unsigned int& nVetoMu_IsoImp, // muons have impact embedded into ID, no need to account for it
	unsigned int& nVetoMu_Iso, unsigned int& nVetoMu_all, // all the veto counters
	bool record, bool debug) // more output

{
pat::MuonCollection selMuons_allIso;
return processMuons_ID_ISO_Kinematics(
	muons, goodPV, weight, // input
	mu_ID,  veto_mu_ID,                          // config/cuts
	mu_ISO, veto_mu_ISO,
	pt_cut, eta_cut, veto_pt_cut, veto_eta_cut,
	// output
	selMuons, selMuons_allIso,
	muDiff,
	nVetoMu_Iso, nVetoMu_all, // all the veto counters
	record, debug);
}


/*
 * loop over all HLT objects, find if anything matches to given muon
 */
bool processMuon_matchesHLTs(
	pat::Muon& muon,
	vector<pat::TriggerObjectStandAlone>& trig_objs,    // input: trigger objects to match against (so, these should match HLT of interest)
	float min_dR
	)
{
bool matched = false;
for (size_t i = 0; i < trig_objs.size(); i++)
	{
	pat::TriggerObjectStandAlone& obj = trig_objs[i];
	if (reco::deltaR(muon, obj) < min_dR)
		{
		matched = true;
		break;
		}
	}
return matched;
}

/*
 */
int processMuons_MatchHLT(
	pat::MuonCollection& muons,
	vector<pat::TriggerObjectStandAlone>& trig_objs,    // input: trigger objects to match against (so, these should match HLT of interest)
	float min_dR,
	pat::MuonCollection& muons_matched
	)
{
// loop over all muons
// for each muon check if it matches in dR (< min_dR) to any of given trigger objects
// save it if yes
// return amount of discarded leptons

unsigned int nDiscarded = 0;
for(unsigned int n=0; n<muons.size(); ++n)
	{
	pat::Muon& muon = muons[n];
	bool matched = false;
	for (size_t i = 0; i < trig_objs.size(); i++)
		{
		pat::TriggerObjectStandAlone& obj = trig_objs[i];
		if (reco::deltaR(muon, obj) < min_dR)
			{
			matched = true;
			muons_matched.push_back(muon);
			break;
			}
		}
	if (!matched) nDiscarded += 1;
	}

return nDiscarded;
}

float relIso(pat::Muon& muon, double rho)
{
float relIso = 0.0; 

float  chIso   = muon.pfIsolationR04().sumChargedHadronPt;
float  nhIso   = muon.pfIsolationR04().sumNeutralHadronEt;
float  gIso    = muon.pfIsolationR04().sumPhotonEt;
float  puchIso = muon.pfIsolationR04().sumPUPt;
      
relIso  = (chIso + TMath::Max(0.,nhIso+gIso-0.5*puchIso)) / muon.pt();

return relIso;
}


