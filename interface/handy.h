//// somehow all these files are not found in compilation????
// but Int_t ctypes are found in compiling the plugin -- ?? how?
//#include "TClassTable.h"
//#include "TSystem.h"
//#include "TROOT.h"

// copied the buildfile from ttbar to root of this project -- it compiled

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

Int_t parse_chain_id(const reco::Candidate* part);
void save_final_states(
	const reco::Candidate * part,
	std::vector<const reco::Candidate*>& saved_particles);

