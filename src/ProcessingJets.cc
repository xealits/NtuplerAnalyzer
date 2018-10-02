#include "UserCode/NtuplerAnalyzer/interface/ProcessingJets.h"
//#include "UserCode/NtuplerAnalyzer/interface/recordFuncs.h"



std::vector<double> smearJER(double pt, double eta, double genPt)
	{
	std::vector<double> toReturn(3,pt);
	if(genPt<=0) return toReturn;

	// 2016 80X
	// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
	//abs(eta) region 	0.0–0.5 	0.5-0.8 	0.8–1.1 	1.1-1.3 	1.3–1.7 	1.7 - 1.9 	1.9–2.1 	2.1 - 2.3 	2.3 - 2.5 	2.5–2.8 	2.8-3.0 	3.0-3.2 	3.2-5.0
	//Data/MC SF
	//1.122 +-0.026
	//1.167 +-0.048
	//1.168 +-0.046
	//1.029 +-0.066
	//1.115 +-0.03
	//1.041 +-0.062
	//1.167 +-0.086
	//1.094 +-0.093
	//1.168 +-0.120
	//1.266 +-0.132
	//1.595 +-0.175
	//0.998 +-0.066
	//1.226 +-0.145
	//
	//Moriond:
	//abs(eta) region	0.0–0.5	0.5-0.8	0.8–1.1	1.1-1.3	1.3–1.7	1.7 - 1.9	1.9–2.1	2.1 - 2.3	2.3 - 2.5	2.5–2.8	2.8-3.0	3.0-3.2	3.2-5.0
	//Data/MC SF	1.109 +-0.008	1.138 +-0.013	1.114 +-0.013	1.123 +-0.024	1.084 +-0.011	1.082 +-0.035	1.140 +-0.047	1.067 +-0.053	1.177 +-0.041	1.364 +-0.039	1.857 +-0.071	1.328 +-0.022	1.16 +-0.029
	//
	//abs(eta) region	0.0–0.5	0.5-0.8	0.8–1.1	1.1-1.3	1.3–1.7	1.7 - 1.9	1.9–2.1	2.1 - 2.3	2.3 - 2.5	2.5–2.8	2.8-3.0	3.0-3.2	3.2-5.0
	//Data/MC SF
	// 1.109 +-0.008
	// 1.138 +-0.013
	// 1.114 +-0.013
	// 1.123 +-0.024
	// 1.084 +-0.011
	// 1.082 +-0.035
	// 1.140 +-0.047
	// 1.067 +-0.053
	// 1.177 +-0.041
	// 1.364 +-0.039
	// 1.857 +-0.071
	// 1.328 +-0.022
	// 1.16  +-0.029
	eta=fabs(eta);
	double ptSF(1.0), ptSF_err(0.06);
	if      (eta<0.5) { ptSF=1.109; ptSF_err=0.008; }
	else if (eta<0.8) { ptSF=1.138; ptSF_err=0.013; }
	else if (eta<1.1) { ptSF=1.114; ptSF_err=0.013; }
	else if (eta<1.3) { ptSF=1.123; ptSF_err=0.024; }
	else if (eta<1.7) { ptSF=1.084; ptSF_err=0.011; }
	else if (eta<1.9) { ptSF=1.082; ptSF_err=0.035; }
	else if (eta<2.1) { ptSF=1.140; ptSF_err=0.047; }
	else if (eta<2.3) { ptSF=1.067; ptSF_err=0.053; }
	else if (eta<2.5) { ptSF=1.177; ptSF_err=0.041; }
	else if (eta<2.8) { ptSF=1.364; ptSF_err=0.039; }
	else if (eta<3.0) { ptSF=1.857; ptSF_err=0.071; }
	else if (eta<3.2) { ptSF=1.328; ptSF_err=0.022; }
	else if (eta<5.0) { ptSF=1.16 ; ptSF_err=0.029; }

	/*
	toReturn[0]=TMath::Max(0., (genPt+ptSF*(pt-genPt))/pt );
	toReturn[1]=TMath::Max(0., (genPt+(ptSF+ptSF_err)*(pt-genPt))/pt );
	toReturn[2]=TMath::Max(0., (genPt+(ptSF-ptSF_err)*(pt-genPt))/pt );
	*/

	// TODO: check new SF-scaling application:
	// from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
	// (called c_JER there)
	toReturn[0]=TMath::Max(0., (1 + (ptSF          - 1)*(pt - genPt)/pt) );
	toReturn[1]=TMath::Max(0., (1 + (ptSF+ptSF_err - 1)*(pt - genPt)/pt) );
	toReturn[2]=TMath::Max(0., (1 + (ptSF-ptSF_err - 1)*(pt - genPt)/pt) );

	return toReturn;
	}

std::vector<double> JER_SF(double eta)
	{
	eta=fabs(eta);
	double ptSF(1.0), ptSF_err(0.06);
	std::vector<double> toReturn(2, ptSF);

	//Moriond:
	//abs(eta) region	0.0–0.5	0.5-0.8	0.8–1.1	1.1-1.3	1.3–1.7	1.7 - 1.9	1.9–2.1	2.1 - 2.3	2.3 - 2.5	2.5–2.8	2.8-3.0	3.0-3.2	3.2-5.0
	//Data/MC SF	1.109 +-0.008	1.138 +-0.013	1.114 +-0.013	1.123 +-0.024	1.084 +-0.011	1.082 +-0.035	1.140 +-0.047	1.067 +-0.053	1.177 +-0.041	1.364 +-0.039	1.857 +-0.071	1.328 +-0.022	1.16 +-0.029
	//
	//abs(eta) region	0.0–0.5	0.5-0.8	0.8–1.1	1.1-1.3	1.3–1.7	1.7 - 1.9	1.9–2.1	2.1 - 2.3	2.3 - 2.5	2.5–2.8	2.8-3.0	3.0-3.2	3.2-5.0
	//Data/MC SF
	// 1.109 +-0.008
	// 1.138 +-0.013
	// 1.114 +-0.013
	// 1.123 +-0.024
	// 1.084 +-0.011
	// 1.082 +-0.035
	// 1.140 +-0.047
	// 1.067 +-0.053
	// 1.177 +-0.041
	// 1.364 +-0.039
	// 1.857 +-0.071
	// 1.328 +-0.022
	// 1.16  +-0.029
	if      (eta<0.5) { ptSF=1.109; ptSF_err=0.008; }
	else if (eta<0.8) { ptSF=1.138; ptSF_err=0.013; }
	else if (eta<1.1) { ptSF=1.114; ptSF_err=0.013; }
	else if (eta<1.3) { ptSF=1.123; ptSF_err=0.024; }
	else if (eta<1.7) { ptSF=1.084; ptSF_err=0.011; }
	else if (eta<1.9) { ptSF=1.082; ptSF_err=0.035; }
	else if (eta<2.1) { ptSF=1.140; ptSF_err=0.047; }
	else if (eta<2.3) { ptSF=1.067; ptSF_err=0.053; }
	else if (eta<2.5) { ptSF=1.177; ptSF_err=0.041; }
	else if (eta<2.8) { ptSF=1.364; ptSF_err=0.039; }
	else if (eta<3.0) { ptSF=1.857; ptSF_err=0.071; }
	else if (eta<3.2) { ptSF=1.328; ptSF_err=0.022; }
	else if (eta<5.0) { ptSF=1.16 ; ptSF_err=0.029; }

	toReturn[0]=TMath::Max(0., ptSF);
	toReturn[1]=TMath::Max(0., ptSF_err);

	return toReturn;
	}





/* 
 * Jet Energy Scale (JES)
 *
 * taken from the latest MacroUtils of llvv
 * corresponds to
 * https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties
 * and here
 * https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Example_implementation
 * they say
 * > It is left up to the user to create separate instances of JetCorrectionUncertainty for any and each of the uncertainty sources.
 *
 * also in https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties
 * > Note: JEC uncertainty sources (i.e. correlations, for advanced users)
 * and they show example with
 * JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("Jec11_V12_Uncertainty_AK5PF.txt")
 */
std::vector<float> smearJES(float pt, float eta, JetCorrectionUncertainty *jecUnc)
	{
	jecUnc->setJetEta(eta);
	jecUnc->setJetPt(pt); // should be the corrected jet pt
	float relShift=fabs(jecUnc->getUncertainty(true));
	std::vector<float> toRet;
	toRet.push_back((1.0+relShift)*pt);
	toRet.push_back((1.0-relShift)*pt);
	return toRet;
	}



bool passPFJetID(jet_id label, pat::Jet jet)
	{
	// Set of cuts from the POG group: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data

	bool looseJetID = false;
	bool tightJetID = false;
	bool passID(false); 

	// float rawJetEn(jet.correctedJet("Uncorrected").energy() );

	double eta = fabs(jet.eta());
 
 	// from the twiki:
 	// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
	float NHF = jet.neutralHadronEnergyFraction();
	float NEMF = jet.neutralEmEnergyFraction();
	float CHF = jet.chargedHadronEnergyFraction();
	float MUF = jet.muonEnergyFraction();
	float CEMF = jet.chargedEmEnergyFraction(); // chargedEmEnergyFraction (relative to uncorrected jet energy)
	float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	float NumNeutralParticles =jet.neutralMultiplicity();
	float CHM = jet.chargedMultiplicity(); 
	// TODO: check if these change corresponding to jet corrections? Or apply corrections after passing the selection?
	//
	// the doc https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_14/doc/html/d6/d00/classpat_1_1Jet.html
	// says these `...Fraction` are respective uncorrected jets
	// which is exactly what is needed by their note:
	// > Note: All fractions are calculated with the raw/uncorrected energy of the jet (only then they add up to unity). So the PF JetID has to be applied before the jet energy corrections.

	// float nhf( (jet.neutralHadronEnergy() + jet.HFHadronEnergy())/rawJetEn );
	// float nef( jet.neutralEmEnergy()/rawJetEn );
	// float cef( jet.chargedEmEnergy()/rawJetEn );
	// float chf( jet.chargedHadronEnergy()/rawJetEn );
	// float nch    = jet.chargedMultiplicity();
	// float nconst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
	// float muf(jet.muonEnergy()/rawJetEn); 

	// at time of 80X all 13TeV (74X, 76X, 80X) recommendations got new specification for |eta| < 2.7
	/*
	if (label=="Loose")
		{
		// passID = ( ((nhf<0.99 && nef<0.99 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.99) || abs(eta)>2.4)) && abs(eta)<=3.0 );
		// the same, just with new names:
		// passID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 
		if (fabs(eta) <= 3.0)
			{
			if (fabs(eta) <= 2.7) passID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4);
			else passID = NEMF<0.90 && NumNeutralParticles > 2;
			}
		else passID = NEMF<0.90 && NumNeutralParticles > 10;
		}
	if (label=="Tight")
		{
		// passID = ( ((nhf<0.90 && nef<0.90 && nconst>1) && ((abs(eta)<=2.4 && chf>0 && nch>0 && cef<0.90) || abs(eta)>2.4)) && abs(eta) <=3.0);
		// passID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=3.0 ;
		// the same, only NHF and NEMF are 0.90 at |eta| < 2.7
		if (fabs(eta) <= 3.0)
			{
			if (fabs(eta) <= 2.7) passID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4);
			else passID = NEMF<0.90 && NumNeutralParticles > 2;
			}
		else passID = NEMF<0.90 && NumNeutralParticles > 10;
		}
	*/

	// and latest (at Moriond17) stuff:
	// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
	// quote:
	// > Note: All fractions are calculated with the raw/uncorrected energy of the jet (only then they add up to unity). So the PF JetID has to be applied before the jet energy corrections.
	// --- khmm....
	// in MINIAOD the jets are already corrected  --> one gets the uncorrected jet and reapplies the corrections
	// > ... collection of AK4 jets slimmedJets, made from ak4PFJetsCHS ... "These have standard jet energy corrections applied (L1FastJet, L2, L3), and a pT cut at 10 GeV"
	// so now one need to get the uncorrected jet --> 
	if (eta <= 2.7)
		{
		looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7 ;
		tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) && abs(eta)<=2.7 ;
		}
	else if (eta <= 3.0)
		{
		looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
		tightJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && abs(eta)>2.7 && abs(eta)<=3.0 );
		}
	else
		{
		looseJetID = (NEMF<0.90 && NumNeutralParticles>10 && abs(eta)>3.0 );
		tightJetID = (NEMF<0.90 && NumNeutralParticles>10 && abs(eta)>3.0 );
		}

	// there is also:
	//tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4) && abs(eta)<=3.0

	// And the abs(eta)>3.0 part, but we never consider such jets, so... Meh!

	if (label == LooseJET)
		return looseJetID;
	if (label == TightJET)
		return tightJetID;

	return passID; 
	}




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
	bool record, bool debug) // more output

{
// the PF ID (in Moriond17 recommended to be applied before corrections)
/*
JME::JetResolution resolution = JME::JetResolution(jetResolutionFileName);
JME::JetResolutionScaleFactor resolution_sf = JME::JetResolutionScaleFactor(jetResolutionSFFileName);
Variation m_systematic_variation = Variation::NOMINAL; // FIXME: it should be in some headers, included before... but remake it somehow
*/

for(size_t ijet=0; ijet<jets.size(); ijet++)
	{
	// the input jets here are, supposedly, slimmedJets from MINIAODs -- which are already corrected
	// -- apparently, the parameters that are used in PF ID calculation use the uncorrected values stored in the jets
	// so, everything is fine anyway
	pat::Jet& jet = jets[ijet];
	/*
	if (record)
		{
		fill_1d(string("control_jet_jets_forID_per_eta"), 100, -3., 3., jet.eta(), weight);
		fill_1d(string("control_jet_jets_forID_per_pt"),  100,  0., 400., jet.pt(), weight);
		}
	*/

	bool passID = passPFJetID(jetID, jet);

	// PUJet ID
	float PUJetID_descriminant = jet.userFloat("pileupJetId:fullDiscriminant");

	/*
	if (record)
		{
		fill_2d(string("control_jet_PUJetID_discriminant_vs_eta"), 100, -1.5, 1.5, 100, -3., 3., PUJetID_descriminant, jet.eta(), weight);
		}
	*/

	bool passPUJetID_Loose  = false;
	bool passPUJetID_Medium = false;
	bool passPUJetID_Tight  = false;
	// assume eta < 2.5 always
	if (jet.pt() > 30.)
		{
		passPUJetID_Tight  = PUJetID_descriminant >  0.86;
		passPUJetID_Medium = PUJetID_descriminant >  0.61;
		passPUJetID_Loose  = PUJetID_descriminant > -0.89;
		}
	else if (jet.pt() > 20.)
		{
		passPUJetID_Tight  = PUJetID_descriminant >  0.69; 
		passPUJetID_Medium = PUJetID_descriminant >  0.18; 
		passPUJetID_Loose  = PUJetID_descriminant > -0.97; 
		}
	else if (jet.pt() > 10.)
		{
		passPUJetID_Tight  = PUJetID_descriminant >  0.69; 
		passPUJetID_Medium = PUJetID_descriminant >  0.18; 
		passPUJetID_Loose  = PUJetID_descriminant > -0.97; 
		}
	else
		{
		passPUJetID_Tight  = PUJetID_descriminant >  0.69; 
		passPUJetID_Medium = PUJetID_descriminant >  0.18; 
		passPUJetID_Loose  = PUJetID_descriminant > -0.97; 
		}

	/*
	if (record)
		{
		if (passPUJetID_Loose)
			{
			fill_1d(string("control_jet_jets_passPUJetID_Loose_per_eta"), 100, -3., 3., jet.eta(), weight);
			fill_1d(string("control_jet_jets_passPUJetID_Loose_per_pt"),  100,  0., 400., jet.pt(), weight);
			}
		if (passPUJetID_Medium)
			{
			fill_1d(string("control_jet_jets_passPUJetID_Medium_per_eta"), 100, -3., 3., jet.eta(), weight);
			fill_1d(string("control_jet_jets_passPUJetID_Medium_per_pt"),  100,  0., 400., jet.pt(), weight);
			}
		if (passPUJetID_Tight)
			{
			fill_1d(string("control_jet_jets_passPUJetID_Tight_per_eta"), 100, -3., 3., jet.eta(), weight);
			fill_1d(string("control_jet_jets_passPUJetID_Tight_per_pt"),  100,  0., 400., jet.pt(), weight);
			}
		}
	*/

	bool passPUJetID        = false;
	if      (jetPUID == LoosePUJET)
		passPUJetID = passPUJetID_Loose;
	else if (jetPUID == MediumPUJET)
		passPUJetID = passPUJetID_Medium;
	else if (jetPUID == TightPUJET)
		passPUJetID = passPUJetID_Tight;

	if (with_PUID)
		passID &= (passPUJetID);

	if (passID)
		{
		selJets.push_back(jet);
		/*
		if (record)
			{
			fill_2d(string("control_jet_jetsIDed_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("control_jet_jetsIDed_phi"), 128, -3.2, 3.2, jet.phi(), weight);
			}
		*/
		}
	}

// The jet corrections

//LorentzVector full_jet_corr(0., 0., 0., 0.);
for(size_t ijet=0; ijet<selJets.size(); ijet++)
	{
	// TODO: so does this mean "in place"?
	pat::Jet& jet = selJets[ijet];

	// for MC smearing
	//TRandom *r3 = new TRandom3();

	/*
	if (record)
		{
		//fill_2d(string("slimmedjet_pt_eta"), 400, 0., 400., 200, -4., 4., jet.pt(), jet.eta(), weight);
		fill_2d(string("control_jet_slimmedjet_pt_eta"), 200, 0., 400., 100, -4., 4., jet.pt(), jet.eta(), weight);
		fill_1d(string("control_jet_slimmedjet_phi"), 128, -3.2, 3.2, jet.phi(), weight);
		}
	*/

	LorentzVector jet_initial_momentum = jet.p4();

	if(debug) cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;

	// no-corrections jet
	LorentzVector rawJet = jet.correctedP4("Uncorrected");

	/*
	if (record)
		fill_2d(string("control_jet_slimmedjet_uncorrected_pt_eta"), 200, 0., 400., 100, -4., 4., rawJet.pt(), rawJet.eta(), weight);
	*/

	if(debug) cout << rawJet.eta() << " " << rawJet.pt() << " " << rawJet.energy() << endl;

	//double toRawSF=jet.correctedJet("Uncorrected").pt()/jet.pt();
	//LorentzVector rawJet(jet*toRawSF);
	// FactorizedCorrector with all jet correction files
	jesCor->setJetEta(rawJet.eta());
	jesCor->setJetPt(rawJet.pt());
	jesCor->setJetA(jet.jetArea());
	jesCor->setRho(rho);
	//jesCor->setNPV(nGoodPV); // not used in PU Jet ID example, shouldn't matter
	float jes_correction = jesCor->getCorrection();
	jet.setP4(rawJet*jes_correction);
	//LorentzVector jesCorJet = (rawJet*jes_correction);
	jet.addUserFloat("jes_correction", jes_correction);

	// jet energy scale has uncertainty
	totalJESUnc->setJetEta(jet.eta());
	totalJESUnc->setJetPt(jet.pt()); // should be the corrected jet pt
	float relShift = fabs(totalJESUnc->getUncertainty(true));
	jet.addUserFloat("jes_correction_relShift", relShift);
	// use it with rawJet*jes_correction*(1 +- relShift)
	// since all these corrections are multiplication of p4
	// I can do this shift whenever I want

	/*
	fill_1d(string("control_jet_slimmedjet_jescorrection"), 400, 0., 2., jes_correction, weight);
	if (record)
		fill_2d(string("control_jet_slimmedjet_jescor_pt_eta"), 200, 0., 400., 100, -4., 4., jet.pt(), jet.eta(), weight);
	*/

	if(debug) cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;

	//smear JER
	//https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
	// if a jet matches well to a genJet -- it can be just scaled and that's fine
	// if id doesn't -- one needs to smear it with a randomizer according to the resolution and other stuff

	double jet_resolution = -1;
	double jer_sf = -1;
	double jer_sf_up = -1;
	double jer_sf_down = -1;
	// and the final factors from these SFs
	double jer_factor = -1, jer_factor_up = -1, jer_factor_down = -1;
	// here is the matching of the jet:
	if(isMC)
		{

		/*
		if (record)
		       {
		       fill_2d(string("control_jet_slimmedjet_mc_jerSmearing_before"), 200, 0., 400., 100, -4., 4., jet.pt(), jet.eta(), weight);
		       }
		*/

		// the JER SF and resolution for the jet:
		//std::vector<double> jer_sf_pair = JER_SF(jet.eta());
		//double jer_sf = jer_sf_pair[0];
		//double jer_resolution = jer_sf_pair[1]; // TODO: not sure about this -- is the table the same as what their tool returns?
		// getting it with the tool from files:
		jet_resolution = resolution.getResolution({{JME::Binning::JetPt, jet.pt()}, {JME::Binning::JetEta, jet.eta()}, {JME::Binning::Rho, rho}});
		//jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, m_systematic_variation);
		jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::NOMINAL);
		jer_sf_up   = resolution_sf.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::UP);
		jer_sf_down = resolution_sf.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::DOWN);


		// it's needed to control smearing calculation -- do it in according branches of genJet-matching bellow
		//fill_1d(string("control_jet_slimmedjet_jet_resolution"), 400, 0., 2., jet_resolution, weight);
		//fill_1d(string("control_jet_slimmedjet_jet_jer_scalefactor"), 400, 0., 2., jer_sf, weight);

		// matching to generation jet:
		//const reco::GenJet* genJet=jet.genJet();
		// the PAT doc describes it as "return the matched generated jet"
		// what's the matching procedure?
		// for now I'll do it manually in dR, as shown in Jet POG example
		// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
		const reco::GenJet* matched_genJet = nullptr;
		//double dR_max = 0.4/2; // 0.4 is the jet cone parameter of AK4 jets, which I use
		// moved it to parameters of the procedure
		for (unsigned int i=0; i<genJets.size(); i++)
			{
			reco::GenJet& genJet = genJets[i];
			double dR = reco::deltaR(jet, genJet);

			if (dR > dR_max) continue;

			double dPt = std::abs(genJet.pt() - jet.pt());
			double dPt_max_factor = 3*jet.pt(); // from twiki
			if (dPt > dPt_max_factor * jet_resolution) continue;

			matched_genJet = &genJet;
			}

		// calculate and apply the smearing factor to jet energy
		// with one of the two algorithms, if jet matched to genJet or not
		if (matched_genJet)
			{ // the scaling is ok
			/*
			fill_1d(string("control_jet_slimmedjet_matched_jet_resolution"), 400, 0., 2., jet_resolution, weight);
			fill_1d(string("control_jet_slimmedjet_matched_jet_jer_scalefactor"), 400, 0., 2., jer_sf, weight);
			*/

			double dPt = jet.pt() - matched_genJet->pt();
			//double genjetpt( genJet ? genJet->pt(): 0.);                    
			//std::vector<double> smearJER=utils::cmssw::smearJER(jet.pt(),jet.eta(),genjetpt);
			// using the local smear:
			//std::vector<double> JER_smearing_factor = smearJER(jet.pt(),jet.eta(),genjetpt);
			//double jer_smearing = JER_smearing_factor[0];
			jer_factor = TMath::Max(0., 1. + (jer_sf - 1) * dPt / jet.pt());
			jer_factor_up   = TMath::Max(0., 1. + (jer_sf_up   - 1) * dPt / jet.pt());
			jer_factor_down = TMath::Max(0., 1. + (jer_sf_down - 1) * dPt / jet.pt());
			jet.setP4(jet.p4()*jer_factor); // same as scaleEnergy in the Jet POG example
			// but they also do MIN_ENERGY thing
			// which is static constexpr const double MIN_JET_ENERGY = 1e-2;

			/*
			fill_1d(string("control_jet_slimmedjet_mc_jerSmearing_scaling"), 100, 0., 2., jer_factor, weight);
			if (record)
				{
				fill_2d(string("control_jet_slimmedjet_mc_jerSmearing_scaling_done"), 50, 0., 400., 50, -4., 4., jet.pt(), jet.eta(), weight);
				}
			*/
			}
		else
			{ // the smearing with gaussian randomizer
			//fill_1d(string("control_jet_slimmedjet_notmatched_jet_resolution"), 400, 0., 2., jet_resolution, weight);
			//fill_1d(string("control_jet_slimmedjet_notmatched_jet_jer_scalefactor"), 400, 0., 2., jer_sf, weight);

			// this is the example:
			//double sigma = jet_resolution * std::sqrt(jer_sf * jer_sf - 1);
			//double smearFactor = 1 + r3->Gaus(0, sigma);
			// this is the twiki:
			double smearFactor      = 1. + r3->Gaus(0, jet_resolution) * std::sqrt(TMath::Max(0., jer_sf*jer_sf - 1.));
			double smearFactor_up   = 1. + r3->Gaus(0, jet_resolution) * std::sqrt(TMath::Max(0., jer_sf_up*jer_sf_up - 1.));
			double smearFactor_down = 1. + r3->Gaus(0, jet_resolution) * std::sqrt(TMath::Max(0., jer_sf_down*jer_sf_down - 1.));
			jer_factor = TMath::Max(0., smearFactor);
			jer_factor_up   = TMath::Max(0., smearFactor_up);
			jer_factor_down = TMath::Max(0., smearFactor_down);

			// multiplying a Gaussian should = to multiplying the sigma
			jet.setP4(jet.p4()*jer_factor);

			/*
			fill_1d(string("control_jet_slimmedjet_mc_jerSmearing_stochastic_smearing"), 400, 0., 2., smearFactor, weight);
			if (record)
				{
				fill_2d(string("control_jet_slimmedjet_mc_jerSmearing_smearing_done"), 200, 0., 400., 100, -4., 4., jet.pt(), jet.eta(), weight);
				}
			*/
			}
		}
	jet.addUserFloat("jet_resolution", jet_resolution);
	jet.addUserFloat("jer_sf", jer_sf);
	jet.addUserFloat("jer_sf_up",   jer_sf_up);
	jet.addUserFloat("jer_sf_down", jer_sf_down);
	jet.addUserFloat("jer_factor", jer_factor);
	jet.addUserFloat("jer_factor_up", jer_factor_up);
	jet.addUserFloat("jer_factor_down", jer_factor_down);

	// FIXME: this is not to be re-set. Check that this is a desired non-feature.
	// i.e. check that the uncorrectedJet remains the same even when the corrected momentum is changed by this routine.
	//to get the raw jet again
	//selJets[ijet].setVal("torawsf",1./(newJECSF*newJERSF));

	// Add the jet momentum correction:
	LorentzVector jet_corr(0., 0., 0., 0.);
	jet_corr = jet.p4() - jet_initial_momentum;
	full_jet_corr += jet_corr;

	/*
	if (record)
		{
		fill_2d(string("control_jet_slimmedjet_jet_corr_pt_eta"), 200,   0., 400., 100, -4., 4., jet_corr.pt(), jet_corr.eta(), weight);
		fill_2d(string("control_jet_slimmedjet_jet_corr_pX_pY"),  100, -50.,  50., 100, -50., 50.,  jet_corr.X(), jet_corr.Y(), weight);
		fill_1d(string("control_jet_slimmedjet_jet_corr_pZ"),     100, -50.,  50., jet_corr.Z(), weight);
		// so, jes & jer shouldn't change the eta of the jet -- only the energy (pt)
		// record the correction per jet's eta and per original pt
		//fill_2d(string("control_jet_slimmedjet_full_jetcor_pt_per_jet_eta"), 400, 0., 400., 200, -4., 4., full_jet_corr.pt(), jet.eta(), weight);
		//fill_2d(string("control_jet_slimmedjet_full_jetcor_pt_per_jet_pt"), 400, 0., 400., 200, -4., 4.,  full_jet_corr.pt(), jet_initial_momentum.eta(), weight);
		}
	*/

	if(debug)
		{
		cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;
		cout << "-" << jet_initial_momentum.eta() << " " << jet_initial_momentum.pt() << " " << jet_initial_momentum.energy() << endl;
		}
	}


std::sort (selJets.begin(),  selJets.end(),  utils::sort_CandidatesByPt);

// ----------------------------------- here is the correctF jet correction point
// Propagate full_jet_corr to MET:
//met.setP4(met.p4() - full_jet_corr); // just return the full correction and propagate in place
// TODO: uncertainties?
// for leptons they are done as:
//met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
//met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
//met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
//met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction



if(debug) cout << "Jet Energy Corrections updated" << endl;


// TODO: should MET corrections be done here?
// METs with corrections
// not used now
//
//double met_pt_values[7];
//met_pt_values[0] = n_met.pt();
//met_pt_values[1] = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp).pt();
//met_pt_values[2] = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnDown).pt();
//met_pt_values[3] = mets[0].shiftedP4(pat::MET::METUncertainty::JetResUp).pt();
//met_pt_values[4] = mets[0].shiftedP4(pat::MET::METUncertainty::JetResDown).pt();
//met_pt_values[5] = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp).pt();
//met_pt_values[6] = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown).pt();


return 0;
}


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
	map<systematic_shift, LorentzVector>& full_jet_corr,
	map<systematic_shift, pat::JetCollection>& selJets,                          // output
	bool record, bool debug) // more output

{
// the PF ID (in Moriond17 recommended to be applied before corrections)
/*
JME::JetResolution resolution = JME::JetResolution(jetResolutionFileName);
JME::JetResolutionScaleFactor resolution_sf = JME::JetResolutionScaleFactor(jetResolutionSFFileName);
Variation m_systematic_variation = Variation::NOMINAL; // FIXME: it should be in some headers, included before... but remake it somehow
*/

for(size_t ijet=0; ijet<jets.size(); ijet++)
	{
	// the input jets here are, supposedly, slimmedJets from MINIAODs -- which are already corrected
	// -- apparently, the parameters that are used in PF ID calculation use the uncorrected values stored in the jets
	// so, everything is fine anyway
	pat::Jet& jet = jets[ijet];

	/*
	if (record)
		{
		fill_1d(string("control_jet_jets_forID_per_eta"), 100, -3., 3., jet.eta(), weight);
		fill_1d(string("control_jet_jets_forID_per_pt"),  100,  0., 400., jet.pt(), weight);
		}
	*/

	bool passID = passPFJetID(jetID, jet);

	// PUJet ID
	float PUJetID_descriminant = jet.userFloat("pileupJetId:fullDiscriminant");

	/*
	if (record)
		{
		fill_2d(string("control_jet_PUJetID_discriminant_vs_eta"), 100, -1.5, 1.5, 100, -3., 3., PUJetID_descriminant, jet.eta(), weight);
		}
	*/

	bool passPUJetID_Loose  = false;
	bool passPUJetID_Medium = false;
	bool passPUJetID_Tight  = false;
	// assume eta < 2.5 always
	if (jet.pt() > 30.)
		{
		passPUJetID_Tight  = PUJetID_descriminant >  0.86;
		passPUJetID_Medium = PUJetID_descriminant >  0.61;
		passPUJetID_Loose  = PUJetID_descriminant > -0.89;
		}
	else if (jet.pt() > 20.)
		{
		passPUJetID_Tight  = PUJetID_descriminant >  0.69; 
		passPUJetID_Medium = PUJetID_descriminant >  0.18; 
		passPUJetID_Loose  = PUJetID_descriminant > -0.97; 
		}
	else if (jet.pt() > 10.)
		{
		passPUJetID_Tight  = PUJetID_descriminant >  0.69; 
		passPUJetID_Medium = PUJetID_descriminant >  0.18; 
		passPUJetID_Loose  = PUJetID_descriminant > -0.97; 
		}
	else
		{
		passPUJetID_Tight  = PUJetID_descriminant >  0.69; 
		passPUJetID_Medium = PUJetID_descriminant >  0.18; 
		passPUJetID_Loose  = PUJetID_descriminant > -0.97; 
		}

	/*
	if (record)
		{
		if (passPUJetID_Loose)
			{
			fill_1d(string("control_jet_jets_passPUJetID_Loose_per_eta"), 100, -3., 3., jet.eta(), weight);
			fill_1d(string("control_jet_jets_passPUJetID_Loose_per_pt"),  100,  0., 400., jet.pt(), weight);
			}
		if (passPUJetID_Medium)
			{
			fill_1d(string("control_jet_jets_passPUJetID_Medium_per_eta"), 100, -3., 3., jet.eta(), weight);
			fill_1d(string("control_jet_jets_passPUJetID_Medium_per_pt"),  100,  0., 400., jet.pt(), weight);
			}
		if (passPUJetID_Tight)
			{
			fill_1d(string("control_jet_jets_passPUJetID_Tight_per_eta"), 100, -3., 3., jet.eta(), weight);
			fill_1d(string("control_jet_jets_passPUJetID_Tight_per_pt"),  100,  0., 400., jet.pt(), weight);
			}
		}
	*/

	bool passPUJetID        = false;
	if      (jetPUID == LoosePUJET)
		passPUJetID = passPUJetID_Loose;
	else if (jetPUID == MediumPUJET)
		passPUJetID = passPUJetID_Medium;
	else if (jetPUID == TightPUJET)
		passPUJetID = passPUJetID_Tight;

	if (with_PUID)
		passID &= (passPUJetID);

	if (passID)
		{
		for ( const auto s : jetSystematics )
			{
			selJets[s].push_back(jet);
			}

		/*
		if (record)
			{
			fill_2d(string("control_jet_jetsIDed_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("control_jet_jetsIDed_phi"), 128, -3.2, 3.2, jet.phi(), weight);
			}
		*/
		}
	}

// The jet corrections

// v6, adding jet corrections and b-tagging
//LorentzVector full_jet_corr(0., 0., 0., 0.);

for ( const auto s : jetSystematics )
	{
	pat::JetCollection& jet_col = selJets[s];
	full_jet_corr[s] = LorentzVector(0,0,0,0);
	for(size_t ijet=0; ijet<jet_col.size(); ijet++)
		{
		// TODO: so does this mean "in place"?
		pat::Jet& jet = jet_col[ijet];

		// for MC smearing
		//TRandom *r3 = new TRandom3();

		/*
		if (record && s == SYS_NOMINAL)
			{
			//fill_2d(string("slimmedjet_pt_eta"), 400, 0., 400., 200, -4., 4., jet.pt(), jet.eta(), weight);
			fill_2d(string("control_jet_slimmedjet_pt_eta"), 200, 0., 400., 100, -4., 4., jet.pt(), jet.eta(), weight);
			fill_1d(string("control_jet_slimmedjet_phi"), 128, -3.2, 3.2, jet.phi(), weight);
			}
		*/

		LorentzVector jet_initial_momentum = jet.p4();

		if(debug) cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;

		// no-corrections jet
		LorentzVector rawJet = jet.correctedP4("Uncorrected");

		/*
		if (record && s == SYS_NOMINAL)
			fill_2d(string("control_jet_slimmedjet_uncorrected_pt_eta"), 200, 0., 400., 100, -4., 4., rawJet.pt(), rawJet.eta(), weight);
		*/

		if(debug) cout << rawJet.eta() << " " << rawJet.pt() << " " << rawJet.energy() << endl;

		//double toRawSF=jet.correctedJet("Uncorrected").pt()/jet.pt();
		//LorentzVector rawJet(jet*toRawSF);
		// FactorizedCorrector with all jet correction files
		jesCor->setJetEta(rawJet.eta());
		jesCor->setJetPt(rawJet.pt());
		jesCor->setJetA(jet.jetArea());
		jesCor->setRho(rho);
		//jesCor->setNPV(nGoodPV); // not used in PU Jet ID example, shouldn't matter
		float jes_correction = jesCor->getCorrection();
		LorentzVector jesCorJet = (rawJet*jes_correction);

		// jet energy scale is corrected
		// now do the uncertainty calculation:
		//std::vector<float> jes_shifted_pts = smearJES(jesCorJet.pt(), jesCorJet.eta(), totalJESUnc);
		totalJESUnc->setJetEta(jesCorJet.eta());
		totalJESUnc->setJetPt(jesCorJet.pt()); // should be the corrected jet pt
		float relShift = fabs(totalJESUnc->getUncertainty(true));

		/*
		// up, down
		if (s == SYS_JES_UP)
			{
			fill_1d(string("control_jet_slimmedjet_jescorrection_") + systematic_shift_names[s], 200, 0., 2., jes_correction*(1+relShift), weight);
			jet.setP4(jesCorJet*(1 + relShift));
			}
		else if (s == SYS_JES_DOWN)
			{
			fill_1d(string("control_jet_slimmedjet_jescorrection_") + systematic_shift_names[s], 200, 0., 2., jes_correction*(1-relShift), weight);
			jet.setP4(jesCorJet*(1 - relShift));
			}
		else
			{
			if (s == SYS_NOMINAL)
				fill_1d(string("control_jet_slimmedjet_jescorrection_") + systematic_shift_names[s], 200, 0., 2., jes_correction, weight);
			jet.setP4(jesCorJet);
			}

		if (record && s == SYS_NOMINAL)
			fill_2d(string("control_jet_slimmedjet_jescor_pt_eta"), 200, 0., 400., 100, -4., 4., jet.pt(), jet.eta(), weight);
		*/

		if(debug) cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;

		//smear Jet Energy Resolution (JER)
		//https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
		// if a jet matches well to a genJet -- it can be just scaled and that's fine
		// if id doesn't -- one needs to smear it with a randomizer according to the resolution and other stuff

		// here is the matching of the jet:
		if(isMC)
			{

			/*
			if (record && s == SYS_NOMINAL)
			       {
			       fill_2d(string("control_jet_slimmedjet_mc_jerSmearing_before"), 200, 0., 400., 100, -4., 4., jet.pt(), jet.eta(), weight);
			       }
			*/
			// the JER SF and resolution for the jet:
			//std::vector<double> jer_sf_pair = JER_SF(jet.eta());
			//double jer_sf = jer_sf_pair[0];
			//double jer_resolution = jer_sf_pair[1]; // TODO: not sure about this -- is the table the same as what their tool returns?
			// getting it with the tool from files:
			double jet_resolution = resolution.getResolution({{JME::Binning::JetPt, jet.pt()}, {JME::Binning::JetEta, jet.eta()}, {JME::Binning::Rho, rho}});

			// the JER systematic is shift in JER scale factor
			// according to
			// https://gitlab.cern.ch/ttH/reference/blob/master/definitions/Moriond17.md#7-systematic-uncertainties
			// and other places
			double jer_sf;
			if (s == SYS_JER_UP)
				jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::UP);
			else if (s == SYS_JER_DOWN)
				jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::DOWN);
			else
				jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::NOMINAL);

			/*
			 * from CondFormats/JetMETObjects/interface/JetResolutionObject.h
			 *
			 * enum class Variation {
			 *   NOMINAL = 0,
			 *   DOWN = 1,
			 *   UP = 2
			 *};
			*/

			// it's needed to control smearing calculation -- do it in according branches of genJet-matching bellow
			//fill_1d(string("control_jet_slimmedjet_jet_resolution"), 400, 0., 2., jet_resolution, weight);
			//fill_1d(string("control_jet_slimmedjet_jet_jer_scalefactor"), 400, 0., 2., jer_sf, weight);

			// matching to generation jet:
			//const reco::GenJet* genJet=jet.genJet();
			// the PAT doc describes it as "return the matched generated jet"
			// what's the matching procedure?
			// for now I'll do it manually in dR, as shown in Jet POG example
			// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
			const reco::GenJet* matched_genJet = nullptr;
			//double dR_max = 0.4/2; // 0.4 is the jet cone parameter of AK4 jets, which I use
			// moved it to parameters of the procedure
			for (unsigned int i=0; i<genJets.size(); i++)
				{
				reco::GenJet& genJet = genJets[i];
				double dR = reco::deltaR(jet, genJet);

				if (dR > dR_max) continue;

				double dPt = std::abs(genJet.pt() - jet.pt());
				double dPt_max_factor = 3*jet.pt(); // from twiki
				if (dPt > dPt_max_factor * jet_resolution) continue;

				matched_genJet = &genJet;
				}

			// calculate and apply the smearing factor to jet energy
			// with one of the two algorithms, if jet matched to genJet or not
			if (matched_genJet)
				{ // the scaling is ok
				/*
				if (s == SYS_NOMINAL)
					{
					fill_1d(string("control_jet_slimmedjet_matched_jet_resolution"), 400, 0., 2., jet_resolution, weight);
					fill_1d(string("control_jet_slimmedjet_matched_jet_jer_scalefactor"), 400, 0., 2., jer_sf, weight);
					}
				*/

				double dPt = jet.pt() - matched_genJet->pt();
				//double genjetpt( genJet ? genJet->pt(): 0.);                    
				//std::vector<double> smearJER=utils::cmssw::smearJER(jet.pt(),jet.eta(),genjetpt);
				// using the local smear:
				//std::vector<double> JER_smearing_factor = smearJER(jet.pt(),jet.eta(),genjetpt);
				//double jer_smearing = JER_smearing_factor[0];
				double jer_smearing = TMath::Max(0., 1. + (jer_sf - 1) * dPt / jet.pt());
				jet.setP4(jet.p4()*jer_smearing); // same as scaleEnergy in the Jet POG example
				// but they also do MIN_ENERGY thing
				// which is static constexpr const double MIN_JET_ENERGY = 1e-2;

				/*
				if (s == SYS_NOMINAL)
					fill_1d(string("control_jet_slimmedjet_mc_jerSmearing_scaling"), 400, 0., 2., jer_smearing, weight);

				if (record && s == SYS_NOMINAL)
					{
					fill_2d(string("control_jet_slimmedjet_mc_jerSmearing_scaling_done"), 200, 0., 400., 100, -4., 4., jet.pt(), jet.eta(), weight);
					}
				*/
				}
			else
				{ // the smearing with gaussian randomizer
				/*
				if (s == SYS_NOMINAL)
					{
					fill_1d(string("control_jet_slimmedjet_notmatched_jet_resolution"), 400, 0., 2., jet_resolution, weight);
					fill_1d(string("control_jet_slimmedjet_notmatched_jet_jer_scalefactor"), 400, 0., 2., jer_sf, weight);
					}
				*/

				// this is the example:
				//double sigma = jet_resolution * std::sqrt(jer_sf * jer_sf - 1);
				//double smearFactor = 1 + r3->Gaus(0, sigma);
				// this is the twiki:
				double smearFactor = 1. + r3->Gaus(0, jet_resolution) * std::sqrt(TMath::Max(0., jer_sf*jer_sf - 1.));
				// multiplying a Gaussian should = to multiplying the sigma
				jet.setP4(jet.p4()*TMath::Max(0., smearFactor));

				/*
				if (s == SYS_NOMINAL)
					fill_1d(string("control_jet_slimmedjet_mc_jerSmearing_stochastic_smearing"), 400, 0., 2., smearFactor, weight);
				if (record && s == SYS_NOMINAL)
					{
					fill_2d(string("control_jet_slimmedjet_mc_jerSmearing_smearing_done"), 200, 0., 400., 100, -4., 4., jet.pt(), jet.eta(), weight);
					}
				*/
				}
			}

		// FIXME: this is not to be re-set. Check that this is a desired non-feature.
		// i.e. check that the uncorrectedJet remains the same even when the corrected momentum is changed by this routine.
		//to get the raw jet again
		//selJets[ijet].setVal("torawsf",1./(newJECSF*newJERSF));

		// Add the jet momentum correction:
		LorentzVector jet_corr(0., 0., 0., 0.);
		jet_corr = jet.p4() - jet_initial_momentum;
		full_jet_corr[s] += jet_corr;

		/*
		if (record && s == SYS_NOMINAL)
			{
			fill_2d(string("control_jet_slimmedjet_jet_corr_pt_eta"), 200,   0., 400., 100, -4., 4., jet_corr.pt(), jet_corr.eta(), weight);
			fill_2d(string("control_jet_slimmedjet_jet_corr_pX_pY"),  100, -50.,  50., 100, -50., 50.,  jet_corr.X(), jet_corr.Y(), weight);
			fill_1d(string("control_jet_slimmedjet_jet_corr_pZ"),     100, -50.,  50., jet_corr.Z(), weight);
			// so, jes & jer shouldn't change the eta of the jet -- only the energy (pt)
			// record the correction per jet's eta and per original pt
			//fill_2d(string("control_jet_slimmedjet_full_jetcor_pt_per_jet_eta"), 400, 0., 400., 200, -4., 4., full_jet_corr.pt(), jet.eta(), weight);
			//fill_2d(string("control_jet_slimmedjet_full_jetcor_pt_per_jet_pt"), 400, 0., 400., 200, -4., 4.,  full_jet_corr.pt(), jet_initial_momentum.eta(), weight);
			}
		*/

		if(debug)
			{
			cout << jet.eta() << " " << jet.pt() << " " << jet.energy() << endl;
			cout << "-" << jet_initial_momentum.eta() << " " << jet_initial_momentum.pt() << " " << jet_initial_momentum.energy() << endl;
			}

		}
	std::sort (jet_col.begin(),  jet_col.end(),  utils::sort_CandidatesByPt);
	}



// ----------------------------------- here is the correctF jet correction point
// Propagate full_jet_corr to MET:
//met.setP4(met.p4() - full_jet_corr); // just return the full correction and propagate in place
// TODO: uncertainties?
// for leptons they are done as:
//met.setUncShift(met.px() - muDiff.px()*0.01, met.py() - muDiff.py()*0.01, met.sumEt() - muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnUp);   //assume 1% uncertainty on muon rochester
//met.setUncShift(met.px() + muDiff.px()*0.01, met.py() + muDiff.py()*0.01, met.sumEt() + muDiff.pt()*0.01, pat::MET::METUncertainty::MuonEnDown); //assume 1% uncertainty on muon rochester
//met.setUncShift(met.px() - elDiff.px()*0.01, met.py() - elDiff.py()*0.01, met.sumEt() - elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnUp);   //assume 1% uncertainty on electron scale correction
//met.setUncShift(met.px() + elDiff.px()*0.01, met.py() + elDiff.py()*0.01, met.sumEt() + elDiff.pt()*0.01, pat::MET::METUncertainty::ElectronEnDown); //assume 1% uncertainty on electron scale correction



if(debug) cout << "Jet Energy Corrections updated" << endl;


// TODO: should MET corrections be done here?
// METs with corrections
// not used now
//
//double met_pt_values[7];
//met_pt_values[0] = n_met.pt();
//met_pt_values[1] = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnUp).pt();
//met_pt_values[2] = mets[0].shiftedP4(pat::MET::METUncertainty::JetEnDown).pt();
//met_pt_values[3] = mets[0].shiftedP4(pat::MET::METUncertainty::JetResUp).pt();
//met_pt_values[4] = mets[0].shiftedP4(pat::MET::METUncertainty::JetResDown).pt();
//met_pt_values[5] = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnUp).pt();
//met_pt_values[6] = mets[0].shiftedP4(pat::MET::METUncertainty::UnclusteredEnDown).pt();


return 0;
}


int processJets_Kinematics(pat::JetCollection& jets, // input
	//bool isMC,
	double weight,
	double pt_cut, double eta_cut,
	pat::JetCollection& selJets,                 // output
	bool record, bool debug) // more output

{

// ------------------------------- JETS SELECTION
// the kinematics only, now (at Moriond17) the IDs are recommended to be applied before corrections

//pat::JetCollection selJets;
//pat::JetCollection selJets20GeV, selJets30GeV;
// TODO: do all jet selection right here
// now selBJets are not used anywhere
// selJets pass cross-cleaning with taus later
// and b-tagging again
double mindphijmet (9999.);
for (unsigned int count_ided_jets = 0, ijet = 0; ijet < jets.size(); ++ijet)
	{
	pat::Jet& jet = jets[ijet];

	/*
	if (record)
		{
		fill_2d(string("control_jet_jetscorrected_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), jet.eta(), weight);
		fill_1d(string("control_jet_jetscorrected_phi"), 128, -3.2, 3.2, jet.phi(), weight);
		}
	*/

	// TODO: what do we do here exactly?
	// a loose selection on jets, and then tighten it later?
	// if (jet.pt() < 15 || fabs (jet.eta()) > 3.0) continue;
	// Was 4.7 in eta. Tightened for computing time. 3.0 ensures that we don't cut associations with leptons (0.4 from 2.4)

	//mc truth for this jet
	//const reco::GenJet * genJet = jet.genJet();
	//TString jetType (genJet && genJet->pt() > 0 ? "truejetsid" : "pujetsid");
	// TODO: this mctruth for jets it is never used in the code

	//jet id (done before corrections)
	//bool passPFloose = passPFJetID(jetID, jet); 
	// bool passPFloose = passPFJetID("Tight", jet); 
	//if (label=="Tight")
	// FIXME: check when pileup ID will come out

	// Jet Kinematics
	double eta = jet.eta();
	double pt  = jet.pt();
	//bool passKino = pt > 30. && fabs(eta) < 2.4;

	if ((fabs(eta) < eta_cut) && pt > pt_cut)
		{
		//if (pt > 30.)
			{
			// trying new fake rate application:
			// if there are only 2 30GeV jets -- they are not considered for faking taus
			// TODO: and what to do with 20GeV jets? in selection and in fake rates?
			//  --- for now I will just add 20-30 GeV jets for fake rates..
			//      one should actually try without them
			//selJets30GeV.push_back(jet);
			// drop the 20-30GeV jets for now
			selJets.push_back(jet);

			/*
			if (record)
				{
				fill_2d(string("control_jet_selJets_pt_eta"), 100, 0., 500., 200, -3., 3., jet.pt(), jet.eta(), weight);
				fill_1d(string("control_jet_selJets_phi"), 64, -3.2, 3.2, jet.phi(), weight);
				}
			*/
			// double dphijmet = fabs (deltaPhi (n_met.phi(), jet.phi()));
			//double dphijmet = fabs (deltaPhi (met.phi(), jet.phi()));
			//if (dphijmet < mindphijmet) mindphijmet = dphijmet;
			// FIXME: mindphijmet is not used anywhere now
			}
		//else if (pt > 20.)
			//{
			//selJets20GeV.push_back(jet);
			//}
		}
	}

return 0;
}

