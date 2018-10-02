#ifndef PROCESSINGJETS_H
#define PROCESSINGJETS_H

#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>
#include <JetMETCorrections/Modules/interface/JetResolution.h>

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TRandom3.h"

#include "UserCode/NtuplerAnalyzer/interface/SystematicShifts.h"
//#include "UserCode/NtuplerAnalyzer/interface/recordFuncs.h"

/*
 * Jet IDs
 * and PU Jet IDs
 */

typedef enum {LooseJET, TightJET} jet_id;
typedef enum {LoosePUJET, MediumPUJET, TightPUJET} pu_jet_id;

std::vector<double> smearJER(double pt, double eta, double genPt);

std::vector<double> JER_SF(double eta);




/* TODO: taken from the latest MacroUtils of llvv
 * need to check it
 * TODO: it seems this isn't used, done inline later
 */
std::vector<float> smearJES(float pt, float eta, JetCorrectionUncertainty *jecUnc);



bool passPFJetID(jet_id, pat::Jet);




int processJets_CorrectJES_SmearJERnJES_ID_ISO(pat::JetCollection& jets, std::vector<reco::GenJet>& genJets, // input
	bool isMC, double weight,
	double rho, unsigned int nGoodPV,
	FactorizedJetCorrector *jesCor,
	JetCorrectionUncertainty *totalJESUnc,
	double dR_max, // for jet matching in jet corrections smearing for MC
	JME::JetResolution& resolution, JME::JetResolutionScaleFactor& resolution_sf, Variation& m_systematic_variation,
	jet_id   & jetID,
	pu_jet_id& jetPUID,
	bool with_PUID,
	//double pt_cut, double eta_cut,
	TRandom3 *r3,   // the randomizer for the smearing
	LorentzVector& full_jet_corr, pat::JetCollection& selJets,                          // output
	bool record, bool debug); // more output



int processJets_CorrectJES_SmearJERnJES_ID_ISO_with_systematics(pat::JetCollection& jets, std::vector<reco::GenJet>& genJets, // input
	bool isMC, double weight,
	double rho, unsigned int nGoodPV,
	FactorizedJetCorrector *jesCor,
	JetCorrectionUncertainty *totalJESUnc,
	double dR_max, // for jet matching in jet corrections smearing for MC
	JME::JetResolution& resolution, JME::JetResolutionScaleFactor& resolution_sf, //Variation& m_systematic_variation,
	jet_id   & jetID,
	pu_jet_id& jetPUID,
	bool with_PUID,
	//double pt_cut, double eta_cut,
	TRandom3 *r3,   // the randomizer for the smearing
	//LorentzVector& full_jet_corr,
	map<systematic_shift, LorentzVector>& full_jet_corr,
	map<systematic_shift, pat::JetCollection>& selJets,                          // output
	bool record, bool debug); // more output


int processJets_Kinematics(pat::JetCollection& jets, // input
	//bool isMC,
	double weight,
	double pt_cut, double eta_cut,
	pat::JetCollection& selJets,                 // output
	bool record, bool debug); // more output

#endif /* PROCESSINGJETS_H */


