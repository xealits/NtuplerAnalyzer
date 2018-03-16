
#include "UserCode/NtuplerAnalyzer/interface/handy.h"

Int_t parse_chain_id_old(const reco::Candidate* part)

{
Int_t chain_id = 0;
Int_t part_charge = (part->pdgId() > 0 ? 1 : -1);

unsigned int n_mothers = part->numberOfMothers();

unsigned int order_multiplier = 1;

bool transition = true;

// to work around collor reconnection I'll have to save all mothers of a particle and
// scan through them if there are target resonanses amoung their mothers
// etc
//std::vector<const reco::Candidate*> particle_pool;
//particle_pool.push_back(part);
// and I just want to know what kind of decay happened here, why all this info is mixed together in one mass?
// a dump from the generator
// really idiotic system
// // it becomes messy
// let's just mix all old crap together somehow

while (n_mothers > 0 && transition)
	{
	//for (unsigned int p_i=0; p_i<particle_pool.size(); p_i++)
	n_mothers = part->numberOfMothers();

	unsigned int part_pdg = abs(part->pdgId());

	transition = false;

	for (unsigned int d_i=0; d_i < n_mothers; d_i++)
		{
		const reco::Candidate * mother = part->mother(d_i);
		// mother must be different and one of the targets
		unsigned int m_pdg = abs(mother->pdgId());
		if (part_pdg == m_pdg)
			{
			part = mother;
			transition = true;
			break;
			}
		// 15 tau 24 W 5 b 6 top
		// so, on of these has to be in the imediate top level
		else if (m_pdg == 15 || m_pdg == 24 || m_pdg == 5 || m_pdg == 6)
			{
			// add to the chain id
			int id = 0;
			switch(m_pdg) {
			case 6:
				id = 4;
				part_charge = (part->pdgId() > 0 ? 1 : -1);
				break;
			case 5:
				id = 3;
				part_charge = (part->pdgId() > 0 ? 1 : -1);
				break;
			case 24:
				id = 2;
				part_charge = (part->pdgId() > 0 ? 1 : -1);
				break;
			case 15:
				id = 1;
				part_charge = (part->pdgId() > 0 ? 1 : -1);
				break;
			default:
				id = 0;
				break;
				}
			chain_id += order_multiplier * id;
			order_multiplier *= 10;
			part = mother;
			transition = true;
			break;
			}
		//else
		//	{
		//	//return chain_id;
		//	part = mother; // let's have the rule of following the first mother -- it should not matter for main chain
		//	// no, apparently here some kind of collor-reconnect may happen
		//	break;
		//	}
		}

	if (!transition)
		//edge case;
		{
		part = part->mother(0); // follow first mother
		// in principle I need to go recursive here to avoid more reconnection
		}
	}

chain_id *= part_charge;

return chain_id;
}

/*
 * save final state info for a decaying particle
 * -- I had problems with saving pointers to candidates
 *    it works for tt, but it seems to crash in other datasets
 */
void save_final_states(
	const reco::Candidate * part,
	std::vector<const reco::Candidate*>& saved_particles)

{
unsigned int pdgId = abs(part->pdgId());
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
		save_final_states(daughter, saved_particles);
		// the recursion is not very deep here, 1-2 leaps
		}
	}
// strictly sequencial algorithm
}



