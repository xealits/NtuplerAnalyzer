#ifndef GENDISTRS_H
#define GENDISTRS_H

#define define_gen_level_pt(name)  name ##_pt  = fs->make<TH1D>( #name "_pt",  #name "_pt",  200,  0, 200)
#define define_gen_level_eta(name) name ##_eta = fs->make<TH1D>( #name "_eta", #name "_eta", 200,  -3.0, 3.0)

#include "Math/Vector4D.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

struct LorentzVector_pointer_pair {
	const LorentzVector* first;
	const LorentzVector* second;
};

struct LorentzVector_pointer_pair sorted_byPt_LorentzVectors(const LorentzVector& v1, const LorentzVector& v2);

#endif /* GENDISTRS_H */
