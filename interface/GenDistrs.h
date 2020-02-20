#ifndef GENDISTRS_H
#define GENDISTRS_H

#include "Math/Vector4D.h"
#include "TH1D.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

/** gen level kinematic histograms

The final state process is defined by the final state leptons.
The ordering of leptons is as follows:
1. order by PDG ID
2. if the ID is the same -- order by pT
3. in case of taus the leptonic and hadronic taus are separated,
   but different hadronic modes (1ch or 3ch etc) are treated as the same

The kinematic distributions saved for each object:
pt and eta of the leading and subleading objects
*/
#define GENDISTR_kinem_histos(name) TH1D \
	*name##_l1_pt,  *name##_l2_pt,  *name##_b1_pt,  *name##_b2_pt, \
	*name##_l1_eta, *name##_l2_eta, *name##_b1_eta, *name##_b2_eta

GENDISTR_kinem_histos(elel);
GENDISTR_kinem_histos(eltaul);
GENDISTR_kinem_histos(eltau1h);
GENDISTR_kinem_histos(eltau3h);

int GenDistrs_make_histos_in_FileService(edm::Service<TFileService>& fs);

struct GenDistrs_general_gen_params {
	Int_t gen_decay_lep1_id;
	Int_t gen_decay_lep2_id;
	const LorentzVector* NT_gen_decay_lep1_p4;
	const LorentzVector* NT_gen_decay_lep2_p4;
	const LorentzVector* NT_gen_decay_jet1_p4;
	const LorentzVector* NT_gen_decay_jet2_p4;
	const LorentzVector* NT_gen_decay_bjet1_p4;
	const LorentzVector* NT_gen_decay_bjet2_p4;
};
int GenDistrs_record(struct GenDistrs_general_gen_params);

#endif /* GENDISTRS_H */
