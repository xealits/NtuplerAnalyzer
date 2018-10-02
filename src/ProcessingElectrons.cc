#include "UserCode/NtuplerAnalyzer/interface/ProcessingElectrons.h"
//#include "UserCode/NtuplerAnalyzer/interface/recordFuncs.h"



int processElectrons_ID_ISO_Kinematics(pat::ElectronCollection& electrons, reco::Vertex goodPV, double rho, double weight, // input
	patUtils::llvvElecId::ElecId el_ID, patUtils::llvvElecId::ElecId veto_el_ID,                                       // config/cuts
	patUtils::llvvElecIso::ElecIso el_ISO, patUtils::llvvElecIso::ElecIso veto_el_ISO,
	double pt_cut, double eta_cut, double veto_pt_cut, double veto_eta_cut,
	// output
	pat::ElectronCollection& selElectrons, pat::ElectronCollection& selElectrons_allIso,
	LorentzVector& elDiff,
	unsigned int& nVetoE_Iso, unsigned int& nVetoE_all, // all vetos
	bool record, bool debug) // more output
{

for(unsigned int count_idiso_electrons = 0, n=0; n<electrons.size (); ++n)
	{
	pat::Electron& electron = electrons[n];

	bool 
		passKin(true),     passId(true),     passIso(true),
		passVetoKin(true), passVetoId(true), passVetoIso(true);
	bool passSigma(false), passSigmaVeto(false);
	bool passImpactParameter(false), passImpactParameterVeto(true);
	// from passId( pat::Photon .. ) of PatUtils:
	//bool elevto = photon.hasPixelSeed();  //LQ  REACTIVATED FOR TIGHT ID, OTHERWISE MANY ELECtRONS pass the photon Id
	// and then, in Tight ID, they include:
	// && !elevto 
	//
	// So, it would be nice to add it to Tight Electron ID:
	//bool hasPixelSeed = electron.hasPixelSeed();
	// but electrons don't have this method
	// will have to cross-clean with photons or etc

	/* did it for excess of 1 highly weighted QCD in e-tau, it didn't change anything
	// removing all electrons close to tight Photons
	// actually, people do it the other way around, testing v9.5
	double minDRlg(9999.);
	for(size_t i=0; i<selPhotons.size(); i++)
		{
		minDRlg = TMath::Min(minDRlg, deltaR(electron.p4(), selPhotons[i].p4()));
		}
	if(minDRlg<0.1) continue;
	*/

	int lid(electron.pdgId()); // should always be 11

	//apply electron corrections
	/* no lepcorrs in 13.6
	if(abs(lid)==11)
		{
		elDiff -= electron.p4();
		ElectronEnCorrector.calibrate(electron, ev.eventAuxiliary().run(), edm::StreamID::invalidStreamID()); 
		//electron = patUtils::GenericLepton(electron.el); //recreate the generic electron to be sure that the p4 is ok
		elDiff += electron.p4();
		}
	*/

	//no need for charge info any longer
	//lid = abs(lid);
	//TString lepStr(lid == 13 ? "mu" : "e");
	// should always be 11!
			
	// no need to mess with photon ID // //veto nearby photon (loose electrons are many times photons...)
	// no need to mess with photon ID // double minDRlg(9999.);
	// no need to mess with photon ID // for(size_t ipho=0; ipho<selPhotons.size(); ipho++)
	// no need to mess with photon ID //   minDRlg=TMath::Min(minDRlg,deltaR(leptons[n].p4(),selPhotons[ipho].p4()));
	// no need to mess with photon ID // if(minDRlg<0.1) continue;


	// ------------------------- electron IDs
	//Cut based identification
	//https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
	// full5x5_sigmaIetaIeta is forgotten in our passId for electrons
	// getting high MC/Data in electrons now -- maybe due to photons,
	// checking sigmaIetaIeta, then photon rejection (hasPixelSeed()) and cross-cleaning with photons
	float eta = std::abs(electron.superCluster()->position().eta());
	float sigmaIetaIeta = electron.full5x5_sigmaIetaIeta();

	// std impact parameter -- the one suggested here
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
	//    barrel    endcap
	// d0   0.05      0.10
	// dz   0.10      0.20
	bool passStdImpactParameter = false;
	if (eta <= 1.479) // barrel, newer selection is precise?
		{
		passSigma =     sigmaIetaIeta < 0.00998; // Tight WP
		passSigmaVeto = sigmaIetaIeta < 0.0115;  // Veto WP
		passStdImpactParameter = electron.gsfTrack()->dxy() < 0.05 && electron.gsfTrack()->dz() < 0.1;
		}
	else if (eta > 1.479) // endcap
		{
		passSigma =     sigmaIetaIeta < 0.0292; // Tight WP
		passSigmaVeto = sigmaIetaIeta < 0.037;  // Veto WP
		passStdImpactParameter = electron.gsfTrack()->dxy() < 0.1 && electron.gsfTrack()->dz() < 0.2;
		}
	// Ichecked this sigma cut with EGamma ID page
	// link is in AN....

	//passImpactParameter = electron.dB() < 0.02;
	// what units is this? in the PAT example on top they say "we use < 0.02cm",
	// in recommendations for muons it is < 0.2 
	// and say "The 2 mm cut preserves efficiency for muons from decays of b and c hadrons"
	//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPATExampleTopQuarks
	//https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
	//
	// the electron ID from parUtils has the impact parameter commented out:
	//                  //dxy             < 0.05     &&
	//                  //dz              < 0.10     &&
	// -- dxy should be the same as dB
	// ..will have to test 0.05 and 0.02 if 0.2 doesn't work well..
	// it is suggested < 0.02 in el ID tables (link in AN...)
	// to study that I do a rough cut here and store dB, dxy dz to ntuple
	//passImpactParameter = electron.dB() < 0.2;
	passImpactParameter = passStdImpactParameter;
	passImpactParameterVeto = passImpactParameter; // same as nominal

	// to include into veto counters
	bool passStdImpactParameterVeto = passStdImpactParameter;

	passId     = patUtils::passId(electron, goodPV, el_ID,      patUtils::CutVersion::Moriond17Cut) && passSigma     && passImpactParameter;
	passVetoId = patUtils::passId(electron, goodPV, veto_el_ID, patUtils::CutVersion::Moriond17Cut) && passSigmaVeto && passImpactParameterVeto;

	// ------------------------- electron isolation

	passIso     = patUtils::passIso(electron, el_ISO,      patUtils::CutVersion::Moriond17Cut, rho);
	passVetoIso = patUtils::passIso(electron, veto_el_ISO, patUtils::CutVersion::Moriond17Cut, rho);


	// ---------------------------- Electron Kinematics
	//double leta(fabs(lid==11 ? lepton.el.superCluster()->eta() : lepton.eta()));
	double leta( electron.superCluster()->eta() );

	// ---------------------- Main lepton kin
	if(electron.pt() < pt_cut)                 passKin = false;
	if(leta > eta_cut)                         passKin = false;
	// also barrel-endcap window (these guys do it: https://gitlab.cern.ch/ttH/reference/blob/master/definitions/Moriond17.md#22-electron)
	if(leta > 1.4442 && leta < 1.5660)     passKin = false; // Crack veto

	// ---------------------- Veto lepton kin
	if (electron.pt () < veto_pt_cut)            passVetoKin = false;
	if (leta > veto_eta_cut)                     passVetoKin = false;
	if (leta > 1.4442 && leta < 1.5660) passVetoKin = false; // Crack veto


	//if (passKin     && passId     && passIso)
	if (passKin && passId) // saving all iso tight leptons for anti-iso QCD region
		{
		selElectrons_allIso.push_back(electron);
		if (passIso)
			selElectrons.push_back(electron);
		}
	else
		{
		// impact is on in IDs by default
		//if(passVetoKin && passVetoId && passVetoIso && passStdImpactParameter) nVetoE_IsoImp++;
		if(passVetoKin && passVetoId && passVetoIso) nVetoE_Iso++;
		if(passVetoKin && passVetoId) nVetoE_all++;
		}

	}

// TODO: there should be no need to sort selected electrons here again -- they are in order of Pt
std::sort(selElectrons.begin(),   selElectrons.end(),   utils::sort_CandidatesByPt);
// std::sort(selLeptons.begin(),   selLeptons.end(),   utils::sort_CandidatesByPt);
// std::sort(selLeptons_nocor.begin(),   selLeptons_nocor.end(),   utils::sort_CandidatesByPt);

return 0;
}

// old compatibility call
int processElectrons_ID_ISO_Kinematics(pat::ElectronCollection& electrons, reco::Vertex goodPV, double rho, double weight, // input
	patUtils::llvvElecId::ElecId el_ID, patUtils::llvvElecId::ElecId veto_el_ID,                                       // config/cuts
	patUtils::llvvElecIso::ElecIso el_ISO, patUtils::llvvElecIso::ElecIso veto_el_ISO,
	double pt_cut, double eta_cut, double veto_pt_cut, double veto_eta_cut,
	// output
	pat::ElectronCollection& selElectrons, LorentzVector& elDiff,
	unsigned int& nVetoE_Iso, unsigned int& nVetoE_all, // all vetos
	bool record, bool debug) // more output

{
pat::ElectronCollection selElectrons_allIso;
unsigned int nVetoE_IsoImp = 0;
 processElectrons_ID_ISO_Kinematics(electrons, goodPV, rho, weight, // input
	el_ID, veto_el_ID,                                       // config/cuts
	el_ISO, veto_el_ISO,
	pt_cut, eta_cut, veto_pt_cut, veto_eta_cut,
	// output
	selElectrons, selElectrons_allIso,
	elDiff,
	nVetoE_Iso, nVetoE_all, // all vetos
	record, debug); // more output
return 0;
}

/*
 * loop over all HLT objects, find if anything matches to given lepton
 */
bool processElectron_matchesHLTs(
	pat::Electron& electron,
	vector<pat::TriggerObjectStandAlone>& trig_objs,    // input: trigger objects to match against (so, these should match HLT of interest)
	float min_dR
	)
{
bool matched = false;
for (size_t i = 0; i < trig_objs.size(); i++)
	{
	pat::TriggerObjectStandAlone& obj = trig_objs[i];
	if (reco::deltaR(electron, obj) < min_dR)
		{
		matched = true;
		break;
		}
	}
return matched;
}

/*
 */
int processElectrons_MatchHLT(
	pat::ElectronCollection& electrons,
	vector<pat::TriggerObjectStandAlone>& trig_objs,    // input: trigger objects to match against (so, these should match HLT of interest)
	float min_dR,
	pat::ElectronCollection& electrons_matched
	)
{
// loop over all electrons
// for each electron check if it matches in dR (< min_dR) to any of given trigger objects
// save it if yes
// return amount of discarded leptons

unsigned int nDiscarded = 0;
for(unsigned int n=0; n<electrons.size(); ++n)
	{
	pat::Electron& electron = electrons[n];
	bool matched = false;
	for (size_t i = 0; i < trig_objs.size(); i++)
		{
		pat::TriggerObjectStandAlone& obj = trig_objs[i];
		if (reco::deltaR(electron, obj) < min_dR)
			{
			matched = true;
			electrons_matched.push_back(electron);
			break;
			}
		}
	if (!matched) nDiscarded += 1;
	}

return nDiscarded;
}

float relIso(pat::Electron& el, double rho)
{
float relIso = 0.0; 

//https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
float  chIso   = el.pfIsolationVariables().sumChargedHadronPt;
float  nhIso   = el.pfIsolationVariables().sumNeutralHadronEt;
float  gIso    = el.pfIsolationVariables().sumPhotonEt;
  
if (rho == 0)
	{
	float  puchIso = el.pfIsolationVariables().sumPUPt; 
	relIso  = (chIso + TMath::Max(0.,nhIso+gIso-0.5*puchIso)) / el.pt();
	}
else
	{
	float effArea = utils::cmssw::getEffectiveArea(11,el.superCluster()->eta());
	relIso  = (chIso + TMath::Max(0.,nhIso+gIso-rho*effArea)) / el.pt();
	}

return relIso;
}


