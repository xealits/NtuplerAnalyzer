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

// gen level kinematic histograms
TH1D *elel_el1_pt, *elel_el2_pt, *elel_b1_pt, *elel_b2_pt,
	*elel_el1_eta, *elel_el2_eta, *elel_b1_eta, *elel_b2_eta,
	*eltaul_el_pt,  *eltaul_tau_pt,  *eltaul_b1_pt,  *eltaul_b2_pt,
	*eltaul_el_eta, *eltaul_tau_eta, *eltaul_b1_eta, *eltaul_b2_eta,
	*eltau1h_el_pt, *eltau1h_tau_pt, *eltau1h_b1_pt, *eltau1h_b2_pt,
	*eltau1h_el_eta, *eltau1h_tau_eta, *eltau1h_b1_eta, *eltau1h_b2_eta,
	*eltau3h_el_pt,  *eltau3h_tau_pt,  *eltau3h_b1_pt,  *eltau3h_b2_pt,
	*eltau3h_el_eta, *eltau3h_tau_eta, *eltau3h_b1_eta, *eltau3h_b2_eta;

int genDistrs_make_histos_in_FileService(edm::Service<TFileService>& fs);

#endif /* GENDISTRS_H */
