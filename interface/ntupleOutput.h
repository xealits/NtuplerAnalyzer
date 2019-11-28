//#ifndef NTUPLEOUTPUT_INTERFACE_H
//#define NTUPLEOUTPUT_INTERFACE_H
// hemmm... it's for the sake of redefining the macro.. kind of messy and will through warnings
// if necessary to not redefine them -- need to split class declaration and initialization in different files, which I don't do at the moment

/*
 * Useful info
 *
 * TNtuple is:
 *   branches of Float_t parameters, named somehow
 *
 * TTree is:
 *   branches of simple Float_t, Int_t parameters or complex classes, known to ROOT
 *   (new classes are given to ROOT by call gROOT->ProcessLine("#include <vector>"))
 *   the creation method for these 2 types is different,
 *   also these methods are different from from TNtuple ones (in better way)
 *
 */

/*
 * What is wanted:
 *   keep the output TTree interface (i.e. the definitions of Branches, their classes and names) in 1 file
 *   easily create or open TTree of this interface in main process
 *   loop over Entires
 *   and have full access to all the branches
 *
 * To do it:
 *   there are macro creating this interface
 *   and the list of them with current definition of the interface is in this file
 *   there are 2 bunches of macro -- for creating Branches of new TTree or for opening existing one
 *   which unfold into commands like
 *      Class_X NT_branchFoo; outputTTreeObject.Branch("branchFoo", "Class_X", &NT_branchFoo);
 *      or
 *      Class_X NT_branchFoo; outputTTreeObject.SetBranchAddress("branchFoo", &NT_branchFoo);
 *
 * the outputTTreeObject is defined in OUTNTUPLE
 * the mode of the interface (create or open ttree) is defined with NTUPLE_INTERFACE_CREATE or NTUPLE_INTERFACE_OPEN
 *
 * there are also a bunch of convenience macro for handling the TNtuple legacy bunches of Float_t parameters
 * -- they mostly should go away when propper objects are used
 *  and there is a macro reseting all the branch parameters -- it's ad-hoc, TODO: do it somehow in more automated, convenient way
 *
 * branch object names are prepended with NT_
 * so a branch named "foo" in the program namespace will have the object named NT_foo
 *
 * also default name of the TTree is NT_output_ttree
 *
 * there is a usage example in a comment further
 */

// default name of the output
#ifndef OUTNTUPLE
	#define OUTNTUPLE NT_output_ttree
#endif

#ifndef COMMA
#define COMMA ,
#endif

// another way of passing a comma, since in CMSSW94 the COMMA does not work anymore
#define UNPACK( ... ) __VA_ARGS__

#define STRINGIFY(x) #x
#define UNPACKS( ... ) #__VA_ARGS__

#ifndef DEFAULT_PARAMETER
#define DEFAULT_PARAMETER -111
#endif

/* macro declaring the object and setting a branch with its' pointer --- all in current, __not_global__ space (in main space)
 *
 * Notice the protocol:
 *    1) the object name in current namespace is `NT_Name`
 *    2) the branch name in the ntuple is `Name`
 */

#if defined(NTUPLE_INTERFACE_CLASS_DECLARE) // ALSO: in EDM Classes NTuple is pointer to TFile Service!
	// to declare vector of types (int/float etc): declate just vector
	#define VECTOR_PARAMs_in_NTuple(NTuple, TYPE, Name)   std::vector<TYPE> NT_##Name;
	// to declare vector of objects: declate vector and a pointer to it
	#define VECTOR_OBJECTs_in_NTuple(NTuple, Name, ...)   __VA_ARGS__ NT_##Name; __VA_ARGS__* pt_NT_##Name;
	// objects and types (simple parameters)
	#define OBJECT_in_NTuple(NTuple, Name, ...)     __VA_ARGS__   NT_##Name;
	#define Float_t_in_NTuple(NTuple, Name)         Float_t NT_##Name;
	#define Int_t_in_NTuple(NTuple, Name)           Int_t   NT_##Name;
	#define ULong64_t_in_NTuple(NTuple, Name)       ULong64_t   NT_##Name;
	#define Bool_t_in_NTuple(NTuple, Name)          Bool_t  NT_##Name;

#elif defined(NTUPLE_INTERFACE_CLASS_INITIALIZE)
	// hook up branch
	#define VECTOR_PARAMs_in_NTuple(NTuple, TYPE, Name)   NTuple->Branch(#Name, &NT_##Name);
	// hook up pointer for vectors of objects
	#define VECTOR_OBJECTs_in_NTuple(NTuple, Name, ...)   pt_NT_##Name = & NT_##Name; NTuple->Branch(#Name, #__VA_ARGS__, &pt_NT_##Name);
	// handle possible commas
	//#define VECTOR_OBJECTs_in_NTuple(NTuple, VECTOR_CLASS_MACRO, Name)   pt_NT_##Name = & NT_##Name; NTuple->Branch(#Name, VECTOR_CLASS_MACRO, &pt_NT_##Name);
	// objects and types (simple parameters)
	#define OBJECT_in_NTuple(NTuple, Name, ...)     NTuple->Branch(#Name, #__VA_ARGS__, &NT_##Name);
	#define Float_t_in_NTuple(NTuple, Name)         NTuple->Branch(#Name, &NT_##Name, #Name "/F");
	#define Int_t_in_NTuple(NTuple, Name)           NTuple->Branch(#Name, &NT_##Name, #Name "/I");
	#define ULong64_t_in_NTuple(NTuple, Name)       NTuple->Branch(#Name, &NT_##Name, #Name "/l");
	#define Bool_t_in_NTuple(NTuple, Name)          NTuple->Branch(#Name, &NT_##Name, #Name "/O");
#elif defined(NTUPLE_INTERFACE_CLASS_RESET)
	#define VECTOR_PARAMs_in_NTuple(NTuple, TYPE, Name)            NT_##Name.clear();
	#define VECTOR_OBJECTs_in_NTuple(NTuple, Name, ...)   NT_##Name.clear();
	// objects and types (simple parameters)
	//#define OBJECT_in_NTuple(NTuple, CLASS, Name)   CLASS   NT_##Name; NTuple.Branch(#Name, #CLASS, &NT_##Name);
	#define OBJECT_in_NTuple(NTuple, Name, ...)
	// WARNING: you'll have to reset your object yourself!
	// and these are defaults:
	#define Float_t_in_NTuple(NTuple, Name)         NT_##Name = DEFAULT_PARAMETER ;
	#define Int_t_in_NTuple(NTuple, Name)           NT_##Name = DEFAULT_PARAMETER ;
	#define ULong64_t_in_NTuple(NTuple, Name)       NT_##Name = DEFAULT_PARAMETER ;
	#define Bool_t_in_NTuple(NTuple, Name)          NT_##Name = false;
// ok, the classes should work
// the following are outdated at the moment:
#elif defined(NTUPLE_INTERFACE_CREATE)
	#define VECTOR_PARAMs_in_NTuple(NTuple, TYPE, Name)   std::vector<TYPE> NT_##Name; NTuple.Branch(#Name, &NT_##Name);
	#define VECTOR_OBJECTs_in_NTuple(NTuple, Name, ...)   __VA_ARGS__ NT_##Name; __VA_ARGS__* pt_NT_##Name ; NTuple.Branch(#Name, #__VA_ARGS__, &pt_NT_##Name);
	#define OBJECT_in_NTuple(NTuple, Name, ...)     __VA_ARGS__   NT_##Name; NTuple.Branch(#Name, #__VA_ARGS__, &NT_##Name);
	#define Float_t_in_NTuple(NTuple, Name)         Float_t NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/F");
	#define Int_t_in_NTuple(NTuple, Name)           Int_t   NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/I");
	#define ULong64_t_in_NTuple(NTuple, Name)       ULong64_t   NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/l");
	#define Bool_t_in_NTuple(NTuple, Name)          Bool_t  NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/O");
#elif defined(NTUPLE_INTERFACE_OPEN)
	#define VECTOR_PARAMs_in_NTuple(NTuple, TYPE, Name)   std::vector<TYPE> NT_##Name; std::vector<TYPE>* pt_NT_##Name = &NT_##Name; NTuple->SetBranchAddress(#Name, &pt_NT_##Name);
	#define VECTOR_OBJECTs_in_NTuple(NTuple, Name, ...)   __VA_ARGS__ NT_##Name; __VA_ARGS__* pt_NT_##Name = &NT_##Name ; NTuple->SetBranchAddress(#Name, &pt_NT_##Name);
	#define OBJECT_in_NTuple(NTuple, Name, ...)     __VA_ARGS__*  NT_##Name = 0; NTuple->SetBranchAddress(#Name, &NT_##Name);
	#define PARAMETER_in_NTuple(NTuple, TYPE, Name)  TYPE   NT_##Name; NTuple->SetBranchAddress(#Name, &NT_##Name);
	#define Float_t_in_NTuple(NTuple, Name)         PARAMETER_in_NTuple(NTuple, Float_t, Name);
	#define Int_t_in_NTuple(NTuple, Name)           PARAMETER_in_NTuple(NTuple, Int_t, Name);
	#define ULong64_t_in_NTuple(NTuple, Name)       PARAMETER_in_NTuple(NTuple, ULong64_t, Name);
	#define Bool_t_in_NTuple(NTuple, Name)          PARAMETER_in_NTuple(NTuple, Bool_t, Name);
#else
	error: set ntuple interface mode
#endif

/*
 * complex objects in vectors
 * from https://root.cern.ch/root/html/tutorials/math/mathcoreVectorCollection.C.html
 *   std::vector<ROOT::Math::XYZTVector>  tracks;
 *   std::vector<ROOT::Math::XYZTVector> * pTracks = &tracks;
 *   t1.Branch("tracks","std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&pTracks);
 *
 * simple objects (float, int) in vectors)
 * from
 *
 *   std::vector<float> vpx;
 *
 *   // Create a TTree
 *   TTree *t = new TTree("tvec","Tree with vectors");
 *   t->Branch("vpx",&vpx);
 *
 */

/*
 * this file ties the interface to our ntuple (TTree) in the current namespace
 *
 * Usage
 * 
 * actual files using it are:
 *   bin/ntupleEvents.cc           (creates the TTree, at line 1542 in main)
 *   test/likelihood_regions.cc    (uses existing TTree, at line 51 in pull_likelihood_regions)
 *
 * Roughly the idea is as follows
 *
 * declare your ntuple:
 *
 *     TTree output("reduced_ttree", "TTree with reduced event data");
 *     // if the name is not `NT_output_ttree` (which is assumed here in the interface)
 *     // define your name for preprocessor:
 *     #define OUTNTUPLE output
 *
 *     // set the mode of the interface to branches (create branches in a new TTree or open branches of existing TTree):
 *     #define NTUPLE_INTERFACE_CREATE
 * 
 *     // load this interface:
 *     #include "ntupleOutput.h"
 *
 * now you have NT_Name objects in the name space and the ntuple has branches "Name" with pointers to these objects
 * you can loop over TTree:
 *
 *     for (Long64_t i=0; i<NT_output_ttree.GetEntries(); i++)
 *         {
 *         NT_output_ttree.GetEntry(i);
 *         ...
 *         }
 *
 * copy objects from event (pseudocode):
 *
 *     NT_foo = events[i]["foo"])
 *
 * or actual example from ntupleEvents:
 *     NT_aMCatNLO_weight = evt->weight();
 * or from likelihood_regions:
 *     if (NT_tau_IDlev_0 != 3. && NT_tau_IDlev_1 != 3.) continue;
 *     
 * when done fill the ntuple:
 *
 *     output.Fill();
 *
 * clearing/reseting of the objects for each event -- currently it is responsibility of the programmer
 * but there is a sketchy macro for this now, it is in development
 * used as in ntupleEvents:
 *     RESET_NTUPLE_PARAMETERS // defaults all parameters
 *
 * -- with no ;
 * that's how sketchy it is
 *
 */


/* TODO: add what's missing from eventSelection
 *       and use proper TLorentzVector etc classes
 *
 * MET filters and lumisection certificate are done on the fly at ntuple production
 * lumi passes after MET filters -- to properly account for it in luminosity
 */

#ifndef NTUPLEOUTPUT_LORENTZVECTOR_H
#define NTUPLEOUTPUT_LORENTZVECTOR_H
// the exact LorentzVector declaration
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Vector_3D;
typedef ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3>>  Error_3D;
#endif /* NTUPLEOUTPUT_LORENTZVECTOR_H */

/* Then:
#define UNPACK(...) __VA_ARGS__
#define FOO(type,name) type name

FOO(UNPACK(std::map<int, int>), map_var);
*/

//#endif /* NTUPLEOUTPUT_INTERFACE_H */

/*
 * And declarations of the whole output
 * (tried splitting it into several files, for lighter/faster processing -- no significant speed-up, and not worth trouble for now)
 */

// COMMON OUTPUT
ULong64_t_in_NTuple(OUTNTUPLE, indexevents)
Int_t_in_NTuple(OUTNTUPLE, runNumber)
Int_t_in_NTuple(OUTNTUPLE, lumi)

Float_t_in_NTuple(OUTNTUPLE, aMCatNLO_weight)
Float_t_in_NTuple(OUTNTUPLE, gen_t_pt)
Float_t_in_NTuple(OUTNTUPLE, gen_tb_pt)
Int_t_in_NTuple(OUTNTUPLE, gen_t_w_decay_id) // = id of lepton (+-11/13/15, sign = sign top?) or 1 for quarks
// if it is tau the ID is multiplied by ID of the decay:
//     11, 13 for leptons (positive, no sign -- the sign is tau's)
//     the 20 + 5*(Nch-1) + Npi0 (just 5*Nch + Npi0 should be 1, 2, 10 etc -- the possible overlap is 3 charged + 1 neutral pi = 11, it's a rare decay, maybe negligible, but still let's add shift by 20)
//     (no overlaps with lepton id-s)
Int_t_in_NTuple(OUTNTUPLE, gen_tb_w_decay_id)

// final states of t/tb b and W:
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_t_w1_final_pdgIds)
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_t_w1_final_statuses)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_t_w1_final_p4s, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_t_w2_final_pdgIds)
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_t_w2_final_statuses)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_t_w2_final_p4s, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_t_b_final_pdgIds)
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_t_b_final_statuses)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_t_b_final_p4s, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_tb_w1_final_pdgIds)
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_tb_w1_final_statuses)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_tb_w1_final_p4s, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_tb_w2_final_pdgIds)
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_tb_w2_final_statuses)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_tb_w2_final_p4s, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_tb_b_final_pdgIds)
VECTOR_PARAMs_in_NTuple (OUTNTUPLE, Int_t, gen_tb_b_final_statuses)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_tb_b_final_p4s, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
// sums of p4-s
OBJECT_in_NTuple(OUTNTUPLE, gen_t_w1_final_p4 , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>)
OBJECT_in_NTuple(OUTNTUPLE, gen_t_w2_final_p4 , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>)
OBJECT_in_NTuple(OUTNTUPLE, gen_t_b_final_p4  , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>)
OBJECT_in_NTuple(OUTNTUPLE, gen_tb_w1_final_p4, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>)
OBJECT_in_NTuple(OUTNTUPLE, gen_tb_w2_final_p4, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>)
OBJECT_in_NTuple(OUTNTUPLE, gen_tb_b_final_p4 , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>)

Int_t_in_NTuple(OUTNTUPLE, gen_pythia8_prompt_leptons_N)  // N leptons with status = 21-29 (pythia 8, "particles from hardest subprocess", "Pythia 8 worksheet" for tutorial at the ASP 2012 Summer School)
//Int_t_in_NTuple(OUTNTUPLE, gen_prompt_leptons_ID) // product of their ID-s, tau ID = pdgID * by (20 + 5*...)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, gen_pythia8_prompt_leptons_IDs) // ID-s of the prompt leptons, tau ID = pdgID * by (20 + 5*...)
// for WJets MC
Int_t_in_NTuple(OUTNTUPLE, gen_N_wdecays)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, gen_wdecays_IDs) // ID-s of the prompt leptons, tau ID = pdgID * by (20 + 5*...)
Int_t_in_NTuple(OUTNTUPLE, gen_N_zdecays)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, gen_zdecays_IDs)

// gen2 (from rivet this with particlelevel thing) leptons and jets
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   gen2_jets_pdgId)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t,  gen2_jets_lep_dR_matched)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, gen2_jets_lep_dR)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   gen2_leptons_pdgId)

VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen2_jets_p4   , std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen2_leptons_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)

//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, gen2_fatjets_pdgId)
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, gen2_fatjets_p4)

/*
 * new homogeneous gen final states
 * should substitute all the above
 */

VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_final_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t, gen_final_PromptFinal)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t, gen_final_fromHardFinal)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t, gen_final_PromptTauDecay)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,  gen_final_pdgId)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,  gen_final_status)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,  gen_final_ndaughters)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,  gen_final_chainId)

/*
// custom ID: +-4 t,  +-3 b,  +- 2W,  +-1 tau
// multiplied by power of 10 in order of the decay chain

VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_PromptDecayed_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, gen_PromptDecayed_pdgId)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, gen_PromptDecayed_status)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, gen_PromptDecayed_chainId)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, gen_PromptDecayed_ndaughters)
*/

// for recoil corrections:
Float_t_in_NTuple(OUTNTUPLE, gen_genPx)
Float_t_in_NTuple(OUTNTUPLE, gen_genPy)
Float_t_in_NTuple(OUTNTUPLE, gen_visPx)
Float_t_in_NTuple(OUTNTUPLE, gen_visPy)

// for Z pt-mass correction:
Float_t_in_NTuple(OUTNTUPLE, genPt)
Float_t_in_NTuple(OUTNTUPLE, genMass)

VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   gen_match_lep_id) // IDs of matched gen particles final states for control
// it just shows number of leptons (el, mu, tau) found in event
// or for tau decay products passing p.isDirectHardProcessTauDecayProductFinalState
// it saves IDs of these products separately for 3ch taus and 1ch (including leptonic taus)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   gen_match_tau_id)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   gen_match_tau3ch_id)
// -- should be enough to distinguish number of taus etc
// it just confirms that the DY true gen taus are correctly found
//Int_t_in_NTuple(OUTNTUPLE, gen_match_w_id)
//Int_t_in_NTuple(OUTNTUPLE, gen_match_b_id)

Int_t_in_NTuple(OUTNTUPLE, gen_NUP) // TODO: add gen info from TTbar
Int_t_in_NTuple(OUTNTUPLE, gen_n_PUP)
Int_t_in_NTuple(OUTNTUPLE, nvtx)
Int_t_in_NTuple(OUTNTUPLE, nvtx_gen)
Float_t_in_NTuple(OUTNTUPLE, fixedGridRhoFastjetAll)

// renormalization/factorization scale weights (normalized already)
Float_t_in_NTuple(OUTNTUPLE, gen_weight_norm)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, gen_weights_renorm_fact)   // the envelope weights
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, gen_weights_pdf_hessians)  // the PDF weights
Float_t_in_NTuple(OUTNTUPLE, gen_weight_pdf_nn30_nominal)
Float_t_in_NTuple(OUTNTUPLE, gen_weight_alphas_1)
Float_t_in_NTuple(OUTNTUPLE, gen_weight_alphas_2)
Float_t_in_NTuple(OUTNTUPLE, gen_weight_too)

// fragmentation and lepton decay tables (for b->hadron processes)
Float_t_in_NTuple(OUTNTUPLE, gen_weight_centralFrag   )
Float_t_in_NTuple(OUTNTUPLE, gen_weight_FragUp        )
Float_t_in_NTuple(OUTNTUPLE, gen_weight_FragDown      )
Float_t_in_NTuple(OUTNTUPLE, gen_weight_PetersonFrag  )
Float_t_in_NTuple(OUTNTUPLE, gen_weight_semilepbrUp   )
Float_t_in_NTuple(OUTNTUPLE, gen_weight_semilepbrDown )

Bool_t_in_NTuple(OUTNTUPLE, METfilterbadChCand) // these two seem to not work
Bool_t_in_NTuple(OUTNTUPLE, METfilterbadPFMuon)
Bool_t_in_NTuple(OUTNTUPLE, pass_basic_METfilters)
Bool_t_in_NTuple(OUTNTUPLE, filters_hbhe             )
Bool_t_in_NTuple(OUTNTUPLE, filters_ecalDeadCellTrig )
Bool_t_in_NTuple(OUTNTUPLE, filters_good_vertices    )
Bool_t_in_NTuple(OUTNTUPLE, filters_eebad            )
Bool_t_in_NTuple(OUTNTUPLE, filters_halo             )
Bool_t_in_NTuple(OUTNTUPLE, filters_halo_super       )

Bool_t_in_NTuple(OUTNTUPLE, filters_noBadMuons       )
Bool_t_in_NTuple(OUTNTUPLE, filters_duplicateMuons   )

Bool_t_in_NTuple(OUTNTUPLE, BadChargedCandidateFilter )
Bool_t_in_NTuple(OUTNTUPLE, BadPFMuonFilter           )

Bool_t_in_NTuple(OUTNTUPLE, HLT_lepMonitor)
Bool_t_in_NTuple(OUTNTUPLE, HLT_el)
Bool_t_in_NTuple(OUTNTUPLE, HLT_mu)
Bool_t_in_NTuple(OUTNTUPLE, HLT_el_low_pt)
Bool_t_in_NTuple(OUTNTUPLE, HLT_mu_low_pt)

Bool_t_in_NTuple(OUTNTUPLE, HLT_jets140)
Bool_t_in_NTuple(OUTNTUPLE, HLT_jets400)

// additional triggers for testing in 2017

Bool_t_in_NTuple(OUTNTUPLE, HLT_el_low_pt32)
Bool_t_in_NTuple(OUTNTUPLE, HLT_el_low_pt28_150HT)
Bool_t_in_NTuple(OUTNTUPLE, HLT_el_low_pt30_35PFJet)

// and dilepton HLTs
Bool_t_in_NTuple(OUTNTUPLE, HLT_elmu_1)
Bool_t_in_NTuple(OUTNTUPLE, HLT_elmu_2)
Bool_t_in_NTuple(OUTNTUPLE, HLT_elmu_3)
Bool_t_in_NTuple(OUTNTUPLE, HLT_elmu_4)

Bool_t_in_NTuple(OUTNTUPLE, HLT_elel_1)
Bool_t_in_NTuple(OUTNTUPLE, HLT_elel_2)

// lep-tau HLTs in 2017

Bool_t_in_NTuple(OUTNTUPLE, HLT_eltau)
Bool_t_in_NTuple(OUTNTUPLE, HLT_mutau1)
Bool_t_in_NTuple(OUTNTUPLE, HLT_mutau2)

Int_t_in_NTuple(OUTNTUPLE, leps_ID)
Int_t_in_NTuple(OUTNTUPLE, nleps)
Int_t_in_NTuple(OUTNTUPLE, njets)  // small eta and Loose PFID
Int_t_in_NTuple(OUTNTUPLE, nbjets) // N b-jets among these
Int_t_in_NTuple(OUTNTUPLE, nbjets_noVLooseTau)
Int_t_in_NTuple(OUTNTUPLE, nbjets_noMediumTau)
Int_t_in_NTuple(OUTNTUPLE, nbjets_noTightTau )

Int_t_in_NTuple(OUTNTUPLE, nalljets)  // all jets, only pt and no leptons in dR
Int_t_in_NTuple(OUTNTUPLE, nallbjets) // N b-jets among them
Int_t_in_NTuple(OUTNTUPLE, ntaus)

// info on veto leps
// ele impact is ON by default since v22
//Bool_t_in_NTuple(OUTNTUPLE, no_std_veto_leps) // std veto lep bool for proc, includes iso of leptons and the suggested impact of electrons (muon impact is embedded into ID)
Bool_t_in_NTuple(OUTNTUPLE, no_iso_veto_leps) // only iso, no suggested impacts of electrons
//Int_t_in_NTuple(OUTNTUPLE, nleps_veto_el_isoimp)
Int_t_in_NTuple(OUTNTUPLE, nleps_veto_el_iso)
Int_t_in_NTuple(OUTNTUPLE, nleps_veto_el_all)
Int_t_in_NTuple(OUTNTUPLE, nleps_veto_mu_iso)
Int_t_in_NTuple(OUTNTUPLE, nleps_veto_mu_all)

// Primary Vertices
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, PV_x)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, PV_y)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, PV_z)
// just errors instead of full covar matrices for now:
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, PV_x_err)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, PV_y_err)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, PV_z_err)

// MET OUTPUT
OBJECT_in_NTuple(OUTNTUPLE, met_init                , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)

OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_UnclusteredEnUp   , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_UnclusteredEnDown , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_JetEnUp           , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_JetEnDown         , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_JetResUp          , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_JetResDown        , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_MuonEnUp          , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_MuonEnDown        , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_ElectronEnUp      , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_ElectronEnDown    , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_TauEnUp           , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_TauEnDown         , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)

//OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_JetResUpSmear     , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
//OBJECT_in_NTuple(OUTNTUPLE, met_init_shift_JetResDownSmear   , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)

OBJECT_in_NTuple(OUTNTUPLE, met_uncorrected         , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_corrected           , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_slimmedMets         , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, met_slimmedMETsMuEGClean, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
OBJECT_in_NTuple(OUTNTUPLE, jets_full_correction    , ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >)
//Float_t_in_NTuple(OUTNTUPLE, pfmetcorr_ex) // corrected with recoil corrections
//Float_t_in_NTuple(OUTNTUPLE, pfmetcorr_ey)

// LEPTONS OUTPUT
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, lep_id)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, lep_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_dxy)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_dz)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_relIso)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t,  lep_passIso)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_dB)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t,  lep_matched_HLT)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_matched_HLT_dR)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   lep_matching_gen)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_matching_gen_dR)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   lep_N_trackerLayersWithMeasurement)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_correction)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_energy_parameter1_r9) // a parameter for energy corrections, r9 in electrons, 0 in muons
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_energy_ScaleUp)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_energy_ScaleDown)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_energy_SmearUp)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_energy_SmearDown)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_energy_2016legacy_ScaleEtUp)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_energy_2016legacy_ScaleEtDown)

// tau gen matches for leptons from taus (in DY and TT)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   lep_matching_gen_collection)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_matching_gen_collection_dR)

// allIso LEPTONS
Int_t_in_NTuple(OUTNTUPLE, leps_ID_allIso)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, lep_alliso_id)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, lep_alliso_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_alliso_relIso)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t,  lep_alliso_matched_HLT)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_alliso_matched_HLT_dR)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   lep_alliso_matching_gen)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, lep_alliso_matching_gen_dR)

// JETS OUTPUT
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, jet_id)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, jet_initial_p4    , std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, jet_p4            , std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, jet_uncorrected_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)

VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   jet_matching_gen)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_matching_gen_dR)

// info on lep-dR-match
// skip or study these jets in proc
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t,  jet_matching_lep)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_matching_lep_dR)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t,  jet_matching_allIso_lep)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_matching_allIso_lep_dR)

//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, jet_matched_genjet_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t,  genjet_matched)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, genjet_pt)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, genjet_dR)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   genjet_i)

VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   jet_PFID)

VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jetCharge) // whatever, to check the idea of b-antib assymetry on charge, also the mis-modelling

VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_area)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_uncorrected_jecFactor)

VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jes_recorrection) // jet_uncorrected_p4 * by this factor should give recorrected jet
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jes_uncertainty)

VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_resolution)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jer_sf)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jer_sf_up)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jer_sf_down)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jer_factor)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jer_factor_up)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_jer_factor_down)

//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_rad)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_etaetaMoment)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_phiphiMoment)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_pu_discr)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_b_discr_original2016)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_b_discr)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_b_DeepCSV_legacy2016_1)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, jet_b_DeepCSV_legacy2016_2)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   jet_hadronFlavour)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   jet_partonFlavour)

// GEN TAUS
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_tau_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)

// for tau momentum transitions study
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_tt_tau_vis_p4  , std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_tt_tau_invis_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, gen_tt_tau_orig_p4 , std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, gen_tt_tau_simpleID)

// TAUS OUTPUT
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_id)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_IDlev)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_decayMode)

VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   tau_matching_gen)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_matching_gen_dR)

// tau gen matches in DY and TT for studying the momentum transition at different stages
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   tau_matching_gen_collection)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_matching_gen_collection_dR)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_matching_gen_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)

// info on lep matching
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t,  tau_matching_lep)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_matching_lep_dR)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t,  tau_matching_allIso_lep)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_matching_allIso_lep_dR)

VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_leading_track_pt)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_leadChargedHadrCand_pt)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_leadNeutralCand_pt)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_leadCand_pt)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t,  tau_hasSecondaryVertex)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_flightLength)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_flightLengthSignificance)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t,   tau_dR_matched_jet) // number of the jet in jet vectors, if no match = -1
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Bool_t, tau_SV_fit_isOk)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_refited_index) // number in the vectors of refited tau

// isolation parameters
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_isolationCands)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_isolationChargedHadrCands)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_isolationGammaCands)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_isolationNeutrHadrCands)

VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_signalCands)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_signalChargedHadrCands)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_signalGammaCands)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_signalNeutrHadrCands)

//// these have 0 in the output
//// isolation parameters
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_isolationPFCands)
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_isolationPFChargedHadrCands)
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_isolationPFGammaCands)
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_n_isolationPFNeutrHadrCands)

// TAUS refit OUTPUT
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_matchingQuality)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_x)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_y)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_z)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_cov, std::vector<ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>>>)
// info on tracks of the tau
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_ntracks) // for control
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS_p4 , std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_OS1_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_OS2_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)

// these are from sigCands
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_OS_p4 , std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS1_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS2_p4, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>)
// closest tracks, indexes in the tracks vectors
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_OS_matched_track)
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_SS1_matched_track)
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_SS2_matched_track)
// assume momentum from sigCands is correct
// get impact parameters (b) from general tracks
// save as ROOT::Math::XYZPoint?
// first -- match quality, dR
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_track_OS_matched_track_dR)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_track_SS1_matched_track_dR)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_track_SS2_matched_track_dR)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_track_matchQ) // just sum of all dR
// track parameters
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_OS_matched_track_vtxkey)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_SS1_matched_track_vtxkey)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_SS2_matched_track_vtxkey)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_OS_matched_track_vtxQ)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_SS1_matched_track_vtxQ)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_SS2_matched_track_vtxQ)

VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_OS_matched_track_pdgId)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_SS1_matched_track_pdgId)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_SS2_matched_track_pdgId)

// the impact parameter is ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>
// I'd like to save TVector3.. TVector3 is not writable to TTree branches without magic
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_OS_matched_track_b  , std::vector<TVector3>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS1_matched_track_b , std::vector<TVector3>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS2_matched_track_b , std::vector<TVector3>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_OS_matched_track_p3 , std::vector<TVector3>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS1_matched_track_p3, std::vector<TVector3>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS2_matched_track_p3, std::vector<TVector3>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_OS_matched_track_b  , std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS1_matched_track_b , std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS2_matched_track_b , std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_OS_matched_track_p3 , std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS1_matched_track_p3, std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>>)
VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS2_matched_track_p3, std::vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>>)

//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_OS_matched_track_b  , std::vector<Vector_3D>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS1_matched_track_b , std::vector<Vector_3D>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS2_matched_track_b , std::vector<Vector_3D>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_OS_matched_track_p3 , std::vector<Vector_3D>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS1_matched_track_p3, std::vector<Vector_3D>)
//VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, tau_SV_fit_track_SS2_matched_track_p3, std::vector<Vector_3D>)

// closest gen tau products, indexes in the products vectors
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_OS_matched_gen)
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_SS1_matched_gen)
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Int_t, tau_SV_fit_track_SS2_matched_gen)
//
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_track_OS_matched_gen_dR)
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_track_SS1_matched_gen_dR)
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_track_SS2_matched_gen_dR)
//VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_fit_track_matchQ_gen) // just sum of all dR

// finally, the geometrical SV
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_geom_flightLen)
VECTOR_PARAMs_in_NTuple(OUTNTUPLE, Float_t, tau_SV_geom_flightLenSign)

// PV REFITTED
Bool_t_in_NTuple(OUTNTUPLE, PV_fit_isOk)
Float_t_in_NTuple(OUTNTUPLE, PV_fit_x)
Float_t_in_NTuple(OUTNTUPLE, PV_fit_y)
Float_t_in_NTuple(OUTNTUPLE, PV_fit_z)
OBJECT_in_NTuple(OUTNTUPLE, PV_cov, ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>>)
//OBJECT_in_NTuple(OUTNTUPLE, UNPACK(ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>>), PV_cov)

// The class requested (vector<math::Error<3>::type>) for "tau_SV_cov" is different from the type of the pointer passed (vector<ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > >)


// undefine all macro
#undef VECTOR_PARAMs_in_NTuple
#undef VECTOR_OBJECTs_in_NTuple
#undef OBJECT_in_NTuple
#undef Float_t_in_NTuple
#undef Int_t_in_NTuple
#undef ULong64_t_in_NTuple
#undef Bool_t_in_NTuple

