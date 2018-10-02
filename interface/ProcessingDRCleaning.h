#ifndef PROCESSINGDRCLEANING_H
#define PROCESSINGDRCLEANING_H

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "UserCode/NtuplerAnalyzer/interface/PatUtils.h"

//#include "UserCode/ttbar-leptons-80X/interface/recordFuncs.h"


/* clean taus of leptons
 */
int crossClean_in_dR(pat::TauCollection& selTaus, std::vector<patUtils::GenericLepton>& leptons,
	float min_dR,
	pat::TauCollection& selTausNoLep, // output
	double weight,
	string control_name,
	bool record, bool debug); // more output

/* clean taus of electrons
 */
int crossClean_in_dR(pat::TauCollection& selTaus, pat::ElectronCollection leptons,
	float min_dR,
	pat::TauCollection& selTausNoLep, // output
	double weight,
	string control_name,
	bool record, bool debug); // more output


/* clean taus of muons
 */
int crossClean_in_dR(pat::TauCollection& selTaus, pat::MuonCollection leptons,
	float min_dR,
	pat::TauCollection& selTausNoLep, // output
	double weight,
	string control_name,
	bool record, bool debug); // more output




/* clean jets of leptons
 */
int crossClean_in_dR(pat::JetCollection& selJets, std::vector<patUtils::GenericLepton>& selLeptons,
	float min_dR,
	pat::JetCollection& selJetsOut, // output
	double weight,
	string control_name,
	bool record, bool debug); // more output



/* clean jets of taus
 */
int crossClean_in_dR(pat::JetCollection& selJets, pat::TauCollection& selTaus,
	float min_dR,
	pat::JetCollection& selJetsOut, // output
	double weight,
	string control_name,
	bool record, bool debug); // more output




/* clean taus of jets
 */
int crossClean_in_dR(pat::TauCollection& selTaus, pat::JetCollection& selJets, 
	float min_dR,
	pat::TauCollection& selTausOut, // output
	double weight,
	string control_name,
	bool record, bool debug); // more output


#endif /* PROCESSINGDRCLEANING_H */

