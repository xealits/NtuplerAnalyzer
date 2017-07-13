#include "UserCode/ttbar-leptons-80X/interface/ntupleOutput_interface.h"

#ifndef NTUPLEOUTPUT_COMMON_H
#define NTUPLEOUTPUT_COMMON_H

Float_t_in_NTuple(OUTNTUPLE, aMCatNLO_weight);
Float_t_in_NTuple(OUTNTUPLE, gen_t_pt);
Float_t_in_NTuple(OUTNTUPLE, gen_tb_pt);
Int_t_in_NTuple(OUTNTUPLE, gen_t_w_decay_id); // = id of lepton (11/13/15, but the sign means which product is lepton: minus=1, plus=2) or 1 for quarks
Int_t_in_NTuple(OUTNTUPLE, gen_tb_w_decay_id);
Int_t_in_NTuple(OUTNTUPLE, NUP_gen); // TODO: add gen info from TTbar
Int_t_in_NTuple(OUTNTUPLE, nvtx);

Bool_t_in_NTuple(OUTNTUPLE, HLT_el);
Bool_t_in_NTuple(OUTNTUPLE, HLT_mu);

Int_t_in_NTuple(OUTNTUPLE, leps_ID);
Int_t_in_NTuple(OUTNTUPLE, nleps);
Int_t_in_NTuple(OUTNTUPLE, njets);
Int_t_in_NTuple(OUTNTUPLE, nbjets);
Int_t_in_NTuple(OUTNTUPLE, ntaus);

OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, met_init);
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, met_uncorrected);
OBJECT_in_NTuple(OUTNTUPLE, LorentzVector, met_corrected);

#define RESET_COMMON \
NT_aMCatNLO_weight = -1; \
NT_gen_t_pt = -1; \
NT_gen_tb_pt = -1; \
NT_gen_t_w_decay_id  = -1; \
NT_gen_tb_w_decay_id = -1; \
NT_NUP_gen = -1; \
NT_nvtx = -1; \
NT_HLT_el = false; \
NT_HLT_mu = false; \
NT_leps_ID = -1; \
NT_nleps = -1; \
NT_njets = -1; \
NT_nbjets = -1; \
NT_ntaus = -1; \
NT_met_init.SetXYZT(0,0,0,0); \
NT_met_uncorrected.SetXYZT(0,0,0,0); \
NT_met_corrected.SetXYZT(0,0,0,0); \

#endif /* NTUPLEOUTPUT_COMMON_H */

