#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput_interface.h"

#ifndef NTUPLEOUTPUT_JETS_H
#define NTUPLEOUTPUT_JETS_H

#define NT_JETS_N 5

#define JET_OUTPUT(num) \
Int_t_in_NTuple(OUTNTUPLE, jet##num##_id); \
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, jet##num##_initial_p4); \
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, jet##num##_p4);  \
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, jet##num##_uncorrected_p4); \
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, jet##num##_matched_genjet_p4); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_jes_correction); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_jes_correction_relShift); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_resolution); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_sf); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_sf_up); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_sf_down); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_jer_factor); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_jer_factor_up); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_jer_factor_down); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_rad); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_pu_discr); \
Float_t_in_NTuple(OUTNTUPLE, jet##num##_b_discr); \
Int_t_in_NTuple(OUTNTUPLE, jet##num##_hadronFlavour); \
Int_t_in_NTuple(OUTNTUPLE, jet##num##_partonFlavour);

//JET_OUTPUT(0)
//JET_OUTPUT(1)
//JET_OUTPUT(2)
//JET_OUTPUT(3)
//JET_OUTPUT(4)

// #define VECTOR_PARAMs_in_NTuple (NTuple, TYPE, Name)   std::vector<TYPE> NT_##Name; NTuple.Branch(#Name, &NT_##Name)
// #define VECTOR_OBJECTs_in_NTuple(NTuple, VECTOR_CLASS, Name)   VECTOR_CLASS NT_##Name; VECTOR_CLASS* pt_NT_##Name ; NTuple.Branch(#Name, #VECTOR_CLASS, &pt_NT_##Name)

#define JETS_OUTPUT \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, jet_id); \
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, std::vector<LorentzVector>, jet_initial_p4); \
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, std::vector<LorentzVector>, jet_p4);  \
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, std::vector<LorentzVector>, jet_uncorrected_p4); \
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, std::vector<LorentzVector>, jet_matched_genjet_p4); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jes_correction); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jes_correction_relShift); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_resolution); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_sf); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_sf_up); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_sf_down); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jer_factor); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jer_factor_up); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jer_factor_down); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_rad); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_pu_discr); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_b_discr); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   jet_hadronFlavour); \
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   jet_partonFlavour);

JETS_OUTPUT

#define NT_jet(i, jet, id_jet_p4, matched_genjet_p4, jet_radius_func, btagger_label) case i: { \
NT_jet ##i ##_id                = jet.pdgId(); \
NT_jet ##i ##_initial_p4        = id_jet_p4; \
NT_jet ##i ##_p4                = jet.p4(); \
NT_jet ##i ##_uncorrected_p4    = jet.correctedP4("Uncorrected"); \
NT_jet ##i ##_matched_genjet_p4 = matched_genjet_p4; \
NT_jet ##i ##_jes_correction = jet.userFloat("jes_correction"); \
NT_jet ##i ##_jes_correction_relShift = jet.userFloat("jes_correction_relShift"); \
NT_jet ##i ##_resolution  = jet.userFloat("jet_resolution"); \
NT_jet ##i ##_sf          = jet.userFloat("jer_sf"); \
NT_jet ##i ##_sf_up       = jet.userFloat("jer_sf_up"); \
NT_jet ##i ##_sf_down     = jet.userFloat("jer_sf_down"); \
NT_jet ##i ##_jer_factor          = jet.userFloat("jer_factor"); \
NT_jet ##i ##_jer_factor_up       = jet.userFloat("jer_factor_up"); \
NT_jet ##i ##_jer_factor_down     = jet.userFloat("jer_factor_down"); \
NT_jet ##i ##_rad               = jet_radius_func(jet); \
NT_jet ##i ##_pu_discr          = jet.userFloat("pileupJetId:fullDiscriminant"); \
NT_jet ##i ##_b_discr           = jet.bDiscriminator(btagger_label); \
NT_jet ##i ##_hadronFlavour     = jet.hadronFlavour(); \
NT_jet ##i ##_partonFlavour     = jet.partonFlavour(); \
break; }

// the automatic reset of all parameters for now
#define RESET_JET(num) \
NT_jet##num##_id = -1; \
NT_jet##num##_initial_p4.SetXYZT(0,0,0,0); \
NT_jet##num##_p4.SetXYZT(0,0,0,0);  \
NT_jet##num##_uncorrected_p4.SetXYZT(0,0,0,0); \
NT_jet##num##_matched_genjet_p4.SetXYZT(0,0,0,0); \
NT_jet##num##_jes_correction = -1; \
NT_jet##num##_jes_correction_relShift = -1; \
NT_jet##num##_resolution = -1; \
NT_jet##num##_sf = -1; \
NT_jet##num##_sf_up = -1; \
NT_jet##num##_sf_down = -1; \
NT_jet##num##_jer_factor = -1; \
NT_jet##num##_jer_factor_up = -1; \
NT_jet##num##_jer_factor_down = -1; \
NT_jet##num##_rad = -1; \
NT_jet##num##_pu_discr = -1; \
NT_jet##num##_b_discr = -1; \
NT_jet##num##_hadronFlavour = -1; \
NT_jet##num##_partonFlavour = -1;

#define RESET_JETS_old \
RESET_JET(0) \
RESET_JET(1) \
RESET_JET(2) \
RESET_JET(3) \
RESET_JET(4)

// the automatic reset of all parameters for now
#define RESET_JETS \
NT_jet_id.clear(); \
NT_jet_initial_p4.clear(); \
NT_jet_p4.clear();  \
NT_jet_uncorrected_p4.clear(); \
NT_jet_matched_genjet_p4.clear(); \
NT_jet_jes_correction.clear(); \
NT_jet_jes_correction_relShift.clear(); \
NT_jet_resolution.clear(); \
NT_jet_sf.clear(); \
NT_jet_sf_up.clear(); \
NT_jet_sf_down.clear(); \
NT_jet_jer_factor.clear(); \
NT_jet_jer_factor_up.clear(); \
NT_jet_jer_factor_down.clear(); \
NT_jet_rad.clear(); \
NT_jet_pu_discr.clear(); \
NT_jet_b_discr.clear(); \
NT_jet_hadronFlavour.clear(); \
NT_jet_partonFlavour.clear();

#endif /* NTUPLEOUTPUT_JETS_H */

