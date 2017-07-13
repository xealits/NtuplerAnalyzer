#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput_interface.h"

#ifndef NTUPLEOUTPUT_TAUS_H
#define NTUPLEOUTPUT_TAUS_H

#define NT_TAUS_N 2

#define TAU_OUTPUT(num) \
Int_t_in_NTuple(OUTNTUPLE, tau##num##_id); \
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, tau##num##_p4); \
Int_t_in_NTuple(OUTNTUPLE, tau##num##_IDlev); \
Float_t_in_NTuple(OUTNTUPLE, tau##num##_leading_track_pt); \
Float_t_in_NTuple(OUTNTUPLE, tau##num##_leadChargedHadrCand_pt); \
Float_t_in_NTuple(OUTNTUPLE, tau##num##_leadNeutralCand_pt); \
Float_t_in_NTuple(OUTNTUPLE, tau##num##_leadCand_pt);

//TAU_OUTPUT(0)
//TAU_OUTPUT(1)

#define TAUS_OUTPUT \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_id); \
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, std::vector<LorentzVector>, tau_p4); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_IDlev); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_leading_track_pt); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_leadChargedHadrCand_pt); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_leadNeutralCand_pt); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_leadCand_pt);

TAUS_OUTPUT

#define NT_tau(i, tau, IDlev) case i: { \
NT_tau ##i ##_id    = tau.pdgId(); \
NT_tau ##i ##_p4    = tau.p4(); \
NT_tau ##i ##_IDlev = IDlev; \
NT_tau ##i ##_leading_track_pt = tau.userFloat("leading_track_pt"); \
NT_tau ##i ##_leadChargedHadrCand_pt = tau.userFloat("leadChargedHadrCand_pt"); \
NT_tau ##i ##_leadNeutralCand_pt     = tau.userFloat("leadNeutralCand_pt"); \
NT_tau ##i ##_leadCand_pt            = tau.userFloat("leadCand_pt"); \
break; }

#define RESET_TAU(num) \
NT_tau##num##_id = -1; \
NT_tau##num##_p4.SetXYZT(0,0,0,0); \
NT_tau##num##_IDlev = -1; \
NT_tau##num##_leading_track_pt = -1; \
NT_tau##num##_leadChargedHadrCand_pt = -1; \
NT_tau##num##_leadNeutralCand_pt = -1; \
NT_tau##num##_leadCand_pt = -1;

#define RESET_TAUS \
NT_tau_id.clear(); \
NT_tau_p4.clear(); \
NT_tau_IDlev.clear(); \
NT_tau_leading_track_pt.clear(); \
NT_tau_leadChargedHadrCand_pt.clear(); \
NT_tau_leadNeutralCand_pt.clear(); \
NT_tau_leadCand_pt.clear();

#define RESET_TAUS_old \
RESET_TAU(0) \
RESET_TAU(1)

#endif /* NTUPLEOUTPUT_TAUS_H */

