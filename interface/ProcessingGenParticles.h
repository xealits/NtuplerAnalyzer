#ifndef PROCESSINGGENPARTICLES_H
#define PROCESSINGGENPARTICLES_H

#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include <string>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;


const reco::Candidate* find_W_decay(const reco::Candidate * W);
const reco::Candidate* find_part_decay(const reco::Candidate * part);
std::string parse_W_decay(const reco::Candidate * W);
unsigned int simple_tau_decay_id(const reco::Candidate * tau);
double top_pT_SF(double x);

void save_final_states(const reco::Candidate * part,
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>& p4_s,
	std::vector<Int_t>& pdgId_s,
	std::vector<Int_t>& status_s,
	std::vector<const reco::Candidate*>& saved_particles);

void save_final_cands(const reco::Candidate * part,
	std::vector<const reco::Candidate*>& saved_particles);

void save_final_cands(const reco::Candidate * part,
	std::vector<LorentzVector>& saved_particles);


#endif /* PROCESSINGGENPARTICLES_H */

