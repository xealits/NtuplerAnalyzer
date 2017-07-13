#ifndef NTUPLEOUTPUT_LEPS_H
#define NTUPLEOUTPUT_LEPS_H
#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput_interface.h"


#define LEPS_OUTPUT \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, lep_id); \
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, std::vector<LorentzVector>, lep_p4); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_dxy); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_dz); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_relIso); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_dB);

LEPS_OUTPUT

#define RESET_LEPS \
NT_lep_id.clear(); \
NT_lep_p4.clear(); \
NT_lep_dxy.clear(); \
NT_lep_dz.clear(); \
NT_lep_relIso.clear(); \
NT_lep_dB.clear();


#endif /* NTUPLEOUTPUT_LEPS_H */

