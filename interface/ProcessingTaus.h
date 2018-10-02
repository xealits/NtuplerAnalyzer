#ifndef PROCESSINGTAUS_H
#define PROCESSINGTAUS_H

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "UserCode/NtuplerAnalyzer/interface/PatUtils.h"
//#include "UserCode/NtuplerAnalyzer/interface/recordFuncs.h"



int processTaus_ID(pat::TauCollection& taus, double weight, // input
	string& tauID_decayMode, //string& tauID,               // config/cuts
	string& tauID_IsoMuons,  string& tauID_IsoElectrons,
	pat::TauCollection& selTaus,                          // output
	bool record, bool debug); // more output

int processTaus_ID_ISO(pat::TauCollection& taus, double weight, // input
	string& tauID_decayMode, string& tauID,               // config/cuts
	string& tauID_IsoMuons,  string& tauID_IsoElectrons,
	pat::TauCollection& selTaus,                          // output
	bool record, bool debug); // more output


int processTaus_Kinematics(pat::TauCollection& taus,          // input
	double weight,
	double pt_cut, double eta_cut,
	pat::TauCollection& selTaus,                          // output
	bool record, bool debug); // more output

#endif /* PROCESSINGTAUS_H */

