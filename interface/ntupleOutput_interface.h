#ifndef NTUPLEOUTPUT_INTERFACE_H
#define NTUPLEOUTPUT_INTERFACE_H

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

/* macro declaring the object and setting a branch with its' pointer --- all in current, __not_global__ space (in main space)
 *
 * Notice the protocol:
 *    1) the object name in current namespace is `NT_Name`
 *    2) the branch name in the ntuple is `Name`
 */
#if defined(NTUPLE_INTERFACE_CREATE)
	#define VECTOR_PARAMs_in_NTuple(NTuple, TYPE, Name)   std::vector<TYPE> NT_##Name; NTuple.Branch(#Name, &NT_##Name)
	#define VECTOR_OBJECTs_in_NTuple(NTuple, VECTOR_CLASS, Name)   VECTOR_CLASS NT_##Name; VECTOR_CLASS* pt_NT_##Name ; NTuple.Branch(#Name, #VECTOR_CLASS, &pt_NT_##Name)
	#define OBJECT_in_NTuple(NTuple, CLASS, Name)   CLASS   NT_##Name; NTuple.Branch(#Name, #CLASS, &NT_##Name)
	#define Float_t_in_NTuple(NTuple, Name)         Float_t NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/F")
	#define Int_t_in_NTuple(NTuple, Name)           Int_t   NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/I")
	#define Bool_t_in_NTuple(NTuple, Name)          Bool_t  NT_##Name; NTuple.Branch(#Name, &NT_##Name, #Name "/O")
#elif defined(NTUPLE_INTERFACE_OPEN)
	#define OBJECT_in_NTuple(NTuple, CLASS, Name)   CLASS*  NT_##Name = 0; NTuple->SetBranchAddress(#Name, &NT_##Name)
	#define PARAMETER_in_NTuple(NTuple, TYPE, Name)  TYPE   NT_##Name; NTuple->SetBranchAddress(#Name, &NT_##Name)
	#define Float_t_in_NTuple(NTuple, Name)         PARAMETER_in_NTuple(NTuple, Float_t, Name)
	#define Int_t_in_NTuple(NTuple, Name)           PARAMETER_in_NTuple(NTuple, Int_t, Name)
	#define Bool_t_in_NTuple(NTuple, Name)          PARAMETER_in_NTuple(NTuple, Bool_t, Name)
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

// the exact LorentzVector declaration
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

#endif /* NTUPLEOUTPUT_INTERFACE_H */

