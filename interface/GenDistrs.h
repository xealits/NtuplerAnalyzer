#ifndef GENDISTRS_H
#define GENDISTRS_H

#include "Math/Vector4D.h"
#include "TH1D.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

struct LorentzVector_pointer_pair {
	const LorentzVector* first;
	const LorentzVector* second;
};

struct LorentzVector_pointer_pair sorted_byPt_LorentzVectors(const LorentzVector& v1, const LorentzVector& v2);

void genDistrs_fill_sorted_pair(const LorentzVector& v1, const LorentzVector& v2);

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

int genDistrs_make_histos_in_FileService(edm::Service<TFileService>& fs);

#endif /* GENDISTRS_H */
