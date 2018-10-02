#include "UserCode/NtuplerAnalyzer/interface/ProcessingDRCleaning.h"
#include "UserCode/NtuplerAnalyzer/interface/recordFuncs.h"


/* clean taus of leptons
 */
int crossClean_in_dR(pat::TauCollection& selTaus, std::vector<patUtils::GenericLepton>& leptons,
	float min_dR,
	pat::TauCollection& selTausNoLep, // output
	double weight,
	string control_name,
	bool record, bool debug) // more output
{

//pat::TauCollection selTausNoLep;
int closest_totaunolep_particle_id = 0; // wonder what is 0 particle
for (size_t itau = 0; itau < selTaus.size(); ++itau)
	{
	pat::Tau& tau = selTaus[itau];
	if(debug) cout << "selecting NoLep taus " << itau << endl;

	// cross-cleaning taus with leptons
	bool overlapWithLepton(false);
	for(int l=0; l<(int)leptons.size();++l)
		{
		if (reco::deltaR(tau, leptons[l])<min_dR)
			{ overlapWithLepton=true; break; }
		}
	if (overlapWithLepton) continue;

	selTausNoLep.push_back(tau);
	// so these are the final taus we use in the selection
	if (record)
		{
		fill_2d(control_name + string("_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
		fill_1d(control_name + string("_phi"), 128, -3.2, 3.2, tau.phi(), weight);
		}

	// for the fake-rate counts (in MC)
	// let's save how many taus we find:
	//increment(string("number_of_tausnolep_found"), 1);
	}

return 0;
}

/* clean taus of electrons
 */
int crossClean_in_dR(pat::TauCollection& selTaus, pat::ElectronCollection leptons,
	float min_dR,
	pat::TauCollection& selTausNoLep, // output
	double weight,
	string control_name,
	bool record, bool debug) // more output
{

//pat::TauCollection selTausNoLep;
int closest_totaunolep_particle_id = 0; // wonder what is 0 particle
for (size_t itau = 0; itau < selTaus.size(); ++itau)
	{
	pat::Tau& tau = selTaus[itau];
	if(debug) cout << "selecting NoLep taus " << itau << endl;

	// cross-cleaning taus with leptons
	bool overlapWithLepton(false);
	for(int l=0; l<(int)leptons.size();++l)
		{
		if (reco::deltaR(tau, leptons[l])<min_dR)
			{ overlapWithLepton=true; break; }
		}
	if (overlapWithLepton) continue;

	selTausNoLep.push_back(tau);
	// so these are the final taus we use in the selection
	if (record)
		{
		fill_2d(control_name + string("_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
		fill_1d(control_name + string("_phi"), 128, -3.2, 3.2, tau.phi(), weight);
		}

	// for the fake-rate counts (in MC)
	// let's save how many taus we find:
	//increment(string("number_of_tausnolep_found"), 1);
	}

return 0;
}


/* clean taus of muons
 */
int crossClean_in_dR(pat::TauCollection& selTaus, pat::MuonCollection leptons,
	float min_dR,
	pat::TauCollection& selTausNoLep, // output
	double weight,
	string control_name,
	bool record, bool debug) // more output
{

//pat::TauCollection selTausNoLep;
int closest_totaunolep_particle_id = 0; // wonder what is 0 particle
for (size_t itau = 0; itau < selTaus.size(); ++itau)
	{
	pat::Tau& tau = selTaus[itau];
	if(debug) cout << "selecting NoLep taus " << itau << endl;

	// cross-cleaning taus with leptons
	bool overlapWithLepton(false);
	for(int l=0; l<(int)leptons.size();++l)
		{
		if (reco::deltaR(tau, leptons[l])<min_dR)
			{ overlapWithLepton=true; break; }
		}
	if (overlapWithLepton) continue;

	selTausNoLep.push_back(tau);
	// so these are the final taus we use in the selection
	if (record)
		{
		fill_2d(control_name + string("_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
		fill_1d(control_name + string("_phi"), 128, -3.2, 3.2, tau.phi(), weight);
		}

	// for the fake-rate counts (in MC)
	// let's save how many taus we find:
	//increment(string("number_of_tausnolep_found"), 1);
	}

return 0;
}




/* clean jets of leptons
 */
int crossClean_in_dR(pat::JetCollection& selJets, std::vector<patUtils::GenericLepton>& selLeptons,
	float min_dR,
	pat::JetCollection& selJetsOut, // output
	double weight,
	string control_name,
	bool record, bool debug) // more output

{
//pat::JetCollection selJetsNoLep;
for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
	{
	pat::Jet& jet = selJets[ijet];

	double minDRlj (9999.);

	for (size_t ilep = 0; ilep < selLeptons.size(); ilep++)
		minDRlj = TMath::Min(minDRlj, reco::deltaR (jet, selLeptons[ilep]));

	if (minDRlj < min_dR) continue;

	selJetsOut.push_back(jet);

	if (record)
		{
		fill_2d(control_name + string("_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
		fill_1d(control_name + string("_phi"), 128, -3.2, 3.2, jet.phi(), weight);
		}
	}

return 0;
}


/* clean jets of taus
 */
int crossClean_in_dR(pat::JetCollection& selJets, pat::TauCollection& selTaus,
	float min_dR,
	pat::JetCollection& selJetsOut, // output
	double weight,
	string control_name,
	bool record, bool debug) // more output

{
for (size_t ijet = 0; ijet < selJets.size(); ++ijet)
	{
	pat::Jet& jet = selJets[ijet];

	double minDRlj (9999.);

	for (size_t ilep = 0; ilep < selTaus.size(); ilep++)
		minDRlj = TMath::Min(minDRlj, reco::deltaR (jet, selTaus[ilep]));

	if (minDRlj < min_dR) continue;

	selJetsOut.push_back(jet);

	if (record)
		{
		fill_2d(control_name + string("_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
		fill_1d(control_name + string("_phi"), 128, -3.2, 3.2, jet.phi(), weight);
		}
	}

return 0;
}



/* clean taus of jets
 */
int crossClean_in_dR(pat::TauCollection& selTaus, pat::JetCollection& selJets, 
	float min_dR,
	pat::TauCollection& selTausOut, // output
	double weight,
	string control_name,
	bool record, bool debug) // more output

{
//pat::TauCollection selTausNoLepNoJet;
for (size_t itau = 0; itau < selTaus.size(); ++itau)
	{
	pat::Tau& tau = selTaus[itau];

	// cross-cleaning taus with jets
	bool overlapWithJet(false);
	for(int n=0; n<(int)selJets.size();++n)
		{
		if (reco::deltaR(tau, selJets[n])<0.4)
			{ overlapWithJet=true; break; }
		}
	if (overlapWithJet) continue;

	// the taus, far from leptons and jets
	selTausOut.push_back(tau);

	if (record)
		{
		fill_2d(control_name + string("_pt_eta"), 250, 0., 500., 200, -4., 4., tau.pt(), tau.eta(), weight);
		fill_1d(control_name + string("_phi"), 128, -3.2, 3.2, tau.phi(), weight);
		}

	// for the fake-rate counts (in MC)
	// let's save how many taus we find:
	// increment(string("number_of_tausNoLepNoJet_found"), weight);
	}

return 0;
}


