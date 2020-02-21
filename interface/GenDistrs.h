#ifndef GENDISTRS_H
#define GENDISTRS_H

#include "Math/Vector4D.h"
#include "TH1D.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

int GenDistrs_make_histos_in_FileService(edm::Service<TFileService>& fs);
int GenDistrs_cleanup(void);

/** the gen object in the event
*/
struct GenDistrs_recorded_gen_objects {
	Int_t lep_ids[2];
	LorentzVector* leps[2];
	LorentzVector* jets[2];
	LorentzVector* bjets[2];
};
int GenDistrs_record(struct GenDistrs_recorded_gen_objects);

#endif /* GENDISTRS_H */
