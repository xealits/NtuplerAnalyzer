#include "UserCode/NtuplerAnalyzer/interface/ProcessingGenParticles.h"
#include "TMath.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;

/*
 * const reco::Candidate* find_W_decay(const reco::Candidate * W) {
 *
 * returns the final W, daughters of which are decay products
 * or returns NULL if something weird happend
 *
 * Usage:
 * const reco::Candidate * W = p.daughter( W_num );
 * const reco::Candidate * W_final = find_W_decay(W);
 * W_final->daughter(0)->pdgId();
 * W_final->daughter(1)->pdgId();
*/
const reco::Candidate* find_W_decay(const reco::Candidate * W) {
	const reco::Candidate * p = W; // eeh.. how C passes const arguments? are they local to function
	int d0_id, d1_id; // ids of decay daughters
	//bool found_decay=false;
	// follow W, whatever happens to it (who knows! defensive programming here)
	// until leptonic/quarkonic decay is found
	// 
	// assume there can be only W->W transitions inbetween (TODO: to check actually)
	//while (!found_decay) {
	int hops_count = 0;
	//while (p->numberOfDaughters() == 1 && hops_count < 100)
	while (hops_count < 100)
		{
		// if there are many daughters
		// then it is a W->W transition, probably W->gamme+W
		// follow Ws
		// if there are no W among daughters -- break
		unsigned int n_daughters = p->numberOfDaughters();
		if (n_daughters > 1)
			{
			// loop over daughters looking for W
			bool found_W_transition = false;
			for (unsigned int i = 0; i < n_daughters; ++i)
				{
				if (abs(p->daughter(i)->pdgId()) == 24)
					{
					p = p->daughter(i);
					found_W_transition = true;
					break;
					}
				}
			// if there were no W among daughters, then it is the decay vertex, get out of the loop and return
			if (!found_W_transition) break;
			}
		else
			{
			p = p->daughter(0);
			}
		hops_count++;
		}

	return p;

	/*
	while (true) {
		int n = p->numberOfDaughters();
		switch(n) {
		case 0: return NULL;
		case 1: // it should be another W->W transition
			p = p->daughter(0);
			break; // there should not be no infinite loop here!
		case 2: // so, it should be the decay
			//found_decay = true;
			return p;
		default: // probably, it can happen for hadronic W?
			return W;
		}
	}
	*/
}

/*
 * the same, but different name
 */
const reco::Candidate* find_part_decay(const reco::Candidate * W) {
	const reco::Candidate * p = W; // eeh.. how C passes const arguments? are they local to function
	int d0_id, d1_id; // ids of decay daughters
	//bool found_decay=false;
	// follow W, whatever happens to it (who knows! defensive programming here)
	// until leptonic/quarkonic decay is found
	// 
	// assume there can be only W->W transitions inbetween (TODO: to check actually)
	//while (!found_decay) {
	while (true) {
		int n = p->numberOfDaughters();
		switch(n) {
		case 0: return p;
		case 1: // it should be another W->W transition
			p = p->daughter(0);
			break; // there should not be no infinite loop here!
		case 2: // so, it should be the decay
			//found_decay = true;
			return p;
		default: // some decays are simulated in multipoint way (b decay, tau)
			return p;
		}
	}
}


string parse_W_decay(const reco::Candidate * W) {
	int d0_id = fabs(W->daughter(0)->pdgId());
	int d1_id = fabs(W->daughter(1)->pdgId());
	if (d0_id == 15 || d1_id == 15 ) return string("tau");
	if (d0_id == 11 || d1_id == 11 ) return string("el");
	if (d0_id == 13 || d1_id == 13 ) return string("mu");
	return string("q"); // FiXME: quite dangerous control-flow!
}

/*
 * works as:
 * decay_id *= simple_tau_decay_id(W_final->daughter(0));
 * = 11, 13 for leptons and 20 + 5*Nch + Npi0 for hadrons
 */
unsigned int simple_tau_decay_id(const reco::Candidate * tau)
{
int id = tau->pdgId();
int st = tau->status(); // TODO: check what is status in decat simulation (pythia for our TTbar set)
int n_daughters = tau->numberOfDaughters();

unsigned int tau_decay_id = 0;
// = 11, 13 for leptons
//   for hadronic decays 20 + 5*(Nch - 1) + Npi0  = 15 + 5*Nch + Npi0
//   pi0 pdgID = 111

// 2 is "fractured state", which the final for tau
while (st != 2)
	{
	id = tau->pdgId();
	st = tau->status(); // TODO: check what is status in decat simulation (pythia for our TTbar set)
	n_daughters = tau->numberOfDaughters();
	if (st == 1) // the "final state" tau, shouldn't happen
		return 1;
	else // loop to the tau daughter of the tau, supposedly there is only 1 daughter, but let's loop just in case (gen photons?)
		{
		for (unsigned int j = 0; j < n_daughters; ++j)
			{
			const reco::Candidate * d = tau->daughter(j);
			unsigned int d_id = abs(d->pdgId());
			if (d_id == 15)
				{
				tau = d;
				break;
				}
			}
		}
	}

//if (st == 2)
for (unsigned int j = 0; j < n_daughters; ++j)
	{
	const reco::Candidate * d = tau->daughter(j);
	unsigned int d_id = abs(d->pdgId());
	//LogInfo ("Demo") << j << " tau daughter ID = " << d->pdgId();
	if (d_id == 12 || d_id == 14 || d_id == 16) continue;    // neutrinos
	if (d_id == 11 || d_id == 13 || d_id == 15) return d_id; // leptons (redundant 15)
	// process hadrons, pi0 and the rest
	if (d_id == 111) // pi0
		tau_decay_id += 1;
	else // charged hadron
		tau_decay_id += 5;
	}

tau_decay_id += 15; // the shift to avoid collision with leptons

return tau_decay_id;
}


/* 
 * status codes:
 *   all generators
 *   1 -- final state
 *   2 -- unstable particle
 *   11-200 -- used by generators
 *
 * https://twiki.cern.ch/twiki/bin/view/Sandbox/KevinSappSandbox
 * PYTHIA 8 Worksheet, Torbj¨orn Sj¨ostrand, Richard Corke
 * pythia 8
 *    status 21-29: Particles from the hardest subprocess
 * -- so, these should be prompt
 */

/*
 * save all no-daughters particles, uniq by their pointer
 * save vectors with p4, pdgId and status
 */
void save_final_states(
	const reco::Candidate * part,
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>& p4_s,
	vector<Int_t>& pdgId_s,
	vector<Int_t>& status_s,
	vector<const reco::Candidate*>& saved_particles)

{
if (part->numberOfDaughters() == 0) // it's a final state
	{
	// check if it is already not saved
	if (std::find(saved_particles.begin(), saved_particles.end(), part) == saved_particles.end())
		{
		// then save it:
		saved_particles.push_back(part);
		p4_s.push_back(part->p4());
		pdgId_s.push_back(part->pdgId());
		status_s.push_back(part->status());
		}
	}
else
	{
	// loop through daughters
	for (int d_i=0; d_i < part->numberOfDaughters(); d_i++)
		{
		const reco::Candidate * daughter = part->daughter(d_i);
		save_final_states(daughter, p4_s, pdgId_s, status_s, saved_particles);
		}
	}

// strictly sequencial algorithm
}

/*
 * save all visible no-daughters candidates (leaf candidates)
 * (no gen neutrinos)
 */
void save_final_cands(
	const reco::Candidate * part,
	vector<const reco::Candidate*>& saved_particles)

{
unsigned int pdgId = abs(part->pdgId());
if (!(pdgId == 12 || pdgId == 14 || pdgId == 16))
	{
	if (part->numberOfDaughters() == 0) // it's a final state
		{
		// check if it is not already saved
		if (std::find(saved_particles.begin(), saved_particles.end(), part) == saved_particles.end())
			{
			// then save it:
			saved_particles.push_back(part);
			}
		}
	else
		{
		// loop through daughters
		for (int d_i=0; d_i < part->numberOfDaughters(); d_i++)
			{
			const reco::Candidate * daughter = part->daughter(d_i);
			save_final_cands(daughter, saved_particles);
			}
		}
	}
// strictly sequencial algorithm
}


void save_final_cands(
	const reco::Candidate * part,
	vector<LorentzVector>& saved_particles)

{
unsigned int pdgId = abs(part->pdgId());
if (!(pdgId == 12 || pdgId == 14 || pdgId == 16))
	{
	if (part->numberOfDaughters() == 0) // it's a final state
		{
		// check if it is not already saved
		LorentzVector part_p4 = part->p4();
		//if (std::find(saved_particles.begin(), saved_particles.end(), part) == saved_particles.end())
		bool dR_match_to_saved = false;
		for (int saved_p=0; saved_p<saved_particles.size(); saved_p++)
			{
			dR_match_to_saved = reco::deltaR(part_p4, saved_particles[saved_p]) < 0.0001;
			if (dR_match_to_saved) break;
			}
		// then save it:
		if (! dR_match_to_saved)
			saved_particles.push_back(part_p4);
		}
	else
		{
		// loop through daughters
		for (int d_i=0; d_i < part->numberOfDaughters(); d_i++)
			{
			const reco::Candidate * daughter = part->daughter(d_i);
			save_final_cands(daughter, saved_particles);
			}
		}
	}
// strictly sequencial algorithm
}

// sumup visible momenta of gen tau products
void sum_final_cands(
	const reco::Candidate * part,
	vector<LorentzVector>& saved_particles,
	LorentzVector& sum_parts,
	bool save_only_visible)

{
unsigned int pdgId = abs(part->pdgId());
if (save_only_visible && !(pdgId == 12 || pdgId == 14 || pdgId == 16))
	{
	if (part->numberOfDaughters() == 0) // it's a final state
		{
		// check if it is not already saved
		LorentzVector part_p4 = part->p4();
		//if (std::find(saved_particles.begin(), saved_particles.end(), part) == saved_particles.end())
		bool dR_match_to_saved = false;
		for (int saved_p=0; saved_p<saved_particles.size(); saved_p++)
			{
			dR_match_to_saved = reco::deltaR(part_p4, saved_particles[saved_p]) < 0.0001;
			if (dR_match_to_saved) break;
			}
		// then save it:
		if (! dR_match_to_saved)
			{
			sum_parts += part_p4;
			saved_particles.push_back(part_p4);
			}

		}
	else
		{
		// loop through daughters
		for (int d_i=0; d_i < part->numberOfDaughters(); d_i++)
			{
			const reco::Candidate * daughter = part->daughter(d_i);
			sum_final_cands(daughter, saved_particles, sum_parts, save_only_visible);
			}
		}
	}

if (!save_only_visible)
	{
	if (pdgId == 12 || pdgId == 14 || pdgId == 16)
		{
		// it's a neutrino, final state, check if it is not saved and sum it
		LorentzVector part_p4 = part->p4();
		//if (std::find(saved_particles.begin(), saved_particles.end(), part) == saved_particles.end())
		bool dR_match_to_saved = false;
		for (int saved_p=0; saved_p<saved_particles.size(); saved_p++)
			{
			dR_match_to_saved = reco::deltaR(part_p4, saved_particles[saved_p]) < 0.0001;
			if (dR_match_to_saved) break;
			}
		// then save it:
		if (! dR_match_to_saved)
			{
			sum_parts += part_p4;
			saved_particles.push_back(part_p4);
			}
		}
	else
		{
		// loop through daughters
		for (int d_i=0; d_i < part->numberOfDaughters(); d_i++)
			{
			const reco::Candidate * daughter = part->daughter(d_i);
			sum_final_cands(daughter, saved_particles, sum_parts, save_only_visible);
			}
		}
	}

// strictly sequencial algorithm
}

/** sumup visible momenta of gen tau products
*/
void sum_separate_final_cands(
	const reco::Candidate * part,
	LorentzVector& sum_parts,
	bool save_only_visible)

{
unsigned int pdgId = abs(part->pdgId());
unsigned int ndaughters = part->numberOfDaughters();
if (ndaughters == 0) // it's a final state
	{
	// handle neutrinos
	if (pdgId == 12 || pdgId == 14 || pdgId == 16)
		{
		if (!save_only_visible) sum_parts += part->p4();
		}
	else // add up everything else
		{
		sum_parts += part->p4();
		}
	}

	else
		{
		// loop through daughters
		for (int d_i=0; d_i < ndaughters; d_i++)
			{
			const reco::Candidate * daughter = part->daughter(d_i);
			sum_separate_final_cands(daughter, sum_parts, save_only_visible);
			}
		}
}

/*
 * save all visible no-daughters candidates (leaf candidates)
 * (no gen neutrinos)
 * with parallel ids
 */
void save_final_cands(
	const reco::Candidate * part,
	vector<LorentzVector>& saved_particles,
	vector<int>& provenance_ids,
	int provenance_id
	)

{
unsigned int pdgId = abs(part->pdgId());
if (!(pdgId == 12 || pdgId == 14 || pdgId == 16))
	{
	if (part->numberOfDaughters() == 0) // it's a final state
		{
		// check if it is not already saved
		LorentzVector part_p4 = part->p4();
		//if (std::find(saved_particles.begin(), saved_particles.end(), part) == saved_particles.end())
		bool dR_match_to_saved = false;
		for (int saved_p=0; saved_p<saved_particles.size(); saved_p++)
			{
			dR_match_to_saved = reco::deltaR(part_p4, saved_particles[saved_p]) < 0.0001;
			if (dR_match_to_saved) break;
			}
		// then save it:
		if (! dR_match_to_saved)
			{
			saved_particles.push_back(part_p4);
			provenance_ids.push_back(provenance_id);
			}
		}
	else
		{
		// loop through daughters
		for (int d_i=0; d_i < part->numberOfDaughters(); d_i++)
			{
			const reco::Candidate * daughter = part->daughter(d_i);
			save_final_cands(daughter, saved_particles, provenance_ids, provenance_id);
			}
		}
	}
// strictly sequencial algorithm
}



// Top pT reweighting
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
// -- the study is done for 7-8 Tev
//    but still recommended for 13TeV

double top_pT_SF(double x)
	{
	// the SF function is SF(x)=exp(a+bx)
	// where x is pT of the top quark (at generation?)
	// sqrt(s) 	channel     	a     	b
	// 7 TeV 	all combined 	0.199 	-0.00166
	// 7 TeV 	l+jets      	0.174 	-0.00137
	// 7 TeV 	dilepton    	0.222 	-0.00197
	// 8 TeV 	all combined 	0.156 	-0.00137
	// 8 TeV 	l+jets       	0.159 	-0.00141
	// 8 TeV 	dilepton     	0.148 	-0.00129
	// 13 TeV	all combined	0.0615	-0.0005
	// -- taking all combined 13 TeV
	double a = 0.0615;
	double b = -0.0005;
	return TMath::Exp(a + b*x);
	}


int parse_chain_id(const reco::Candidate& part)
{
unsigned int pdgId = abs(part.pdgId());
unsigned int n_mothers = part.numberOfMothers();
edm::LogInfo ("Demo") << "parse_chain_id: " << pdgId << ' ' << n_mothers;

// check pdgIDs of all mothers
for (unsigned int m_i=0; m_i<n_mothers; m_i++)
	{
	unsigned int m_pdg = abs(part.mother(m_i)->pdgId());
	edm::LogInfo ("Demo") << "has_mother: m " << m_pdg;
	if      (m_pdg == 24) return m_pdg;
	else if (m_pdg == 23) return m_pdg;
	else if (m_pdg == 6)  return m_pdg;
	}
// if still no match -- follow the first mother
if (n_mothers > 0)
	return parse_chain_id(*(part.mother(0)));
else
	return 0;
}

