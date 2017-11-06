#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include <Math/VectorUtil.h>

#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

// the exact LorentzVector declaration
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

TString recoil_corrections_data_file("/HTT-utilities/RecoilCorrections/data/TypeIPFMET_2016BCD.root");
//TString recoil_corrections_data_file("/HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root");
RecoilCorrector* recoilPFMetCorrector = new RecoilCorrector(recoil_corrections_data_file);

float met_pt_recoilcor(float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets)
	{
	float NT_pfmetcorr_ex, NT_pfmetcorr_ey;
	recoilPFMetCorrector->CorrectByMeanResolution(
		met_px, // uncorrected type I pf met px (float)
		met_py, // uncorrected type I pf met py (float)
		gen_genPx, // generator Z/W/Higgs px (float)
		gen_genPy, // generator Z/W/Higgs py (float)
		gen_visPx, // generator visible Z/W/Higgs px (float)
		gen_visPy, // generator visible Z/W/Higgs py (float)
		njets,  // number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
		NT_pfmetcorr_ex, // corrected type I pf met px (float)
		NT_pfmetcorr_ey  // corrected type I pf met py (float)
		);
	return sqrt(NT_pfmetcorr_ey*NT_pfmetcorr_ey + NT_pfmetcorr_ex*NT_pfmetcorr_ex);
	}

float met_pt_recoilcor_x(float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets)
	{
	float NT_pfmetcorr_ex, NT_pfmetcorr_ey;
	recoilPFMetCorrector->CorrectByMeanResolution(
		met_px, // uncorrected type I pf met px (float)
		met_py, // uncorrected type I pf met py (float)
		gen_genPx, // generator Z/W/Higgs px (float)
		gen_genPy, // generator Z/W/Higgs py (float)
		gen_visPx, // generator visible Z/W/Higgs px (float)
		gen_visPy, // generator visible Z/W/Higgs py (float)
		njets,  // number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
		NT_pfmetcorr_ex, // corrected type I pf met px (float)
		NT_pfmetcorr_ey  // corrected type I pf met py (float)
		);
	return NT_pfmetcorr_ex;
	}

float met_pt_recoilcor_y(float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets)
	{
	float NT_pfmetcorr_ex, NT_pfmetcorr_ey;
	recoilPFMetCorrector->CorrectByMeanResolution(
		met_px, // uncorrected type I pf met px (float)
		met_py, // uncorrected type I pf met py (float)
		gen_genPx, // generator Z/W/Higgs px (float)
		gen_genPy, // generator Z/W/Higgs py (float)
		gen_visPx, // generator visible Z/W/Higgs px (float)
		gen_visPy, // generator visible Z/W/Higgs py (float)
		njets,  // number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
		NT_pfmetcorr_ex, // corrected type I pf met px (float)
		NT_pfmetcorr_ey  // corrected type I pf met py (float)
		);
	return NT_pfmetcorr_ey;
	}

float MTlep_met_pt_recoilcor(float lep_px, float lep_py,
	float met_px, float met_py,
	float gen_genPx, float gen_genPy, float gen_visPx, float gen_visPy,
	int njets)
	{
	float NT_pfmetcorr_ex, NT_pfmetcorr_ey;
	recoilPFMetCorrector->CorrectByMeanResolution(
		met_px, // uncorrected type I pf met px (float)
		met_py, // uncorrected type I pf met py (float)
		gen_genPx, // generator Z/W/Higgs px (float)
		gen_genPy, // generator Z/W/Higgs py (float)
		gen_visPx, // generator visible Z/W/Higgs px (float)
		gen_visPy, // generator visible Z/W/Higgs py (float)
		njets,  // number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
		NT_pfmetcorr_ex, // corrected type I pf met px (float)
		NT_pfmetcorr_ey  // corrected type I pf met py (float)
		);
	//return sqrt(NT_pfmetcorr_ey*NT_pfmetcorr_ey + NT_pfmetcorr_ex*NT_pfmetcorr_ex);
	double v1v2 = sqrt((lep_px*lep_px + lep_py*lep_py)*(NT_pfmetcorr_ey*NT_pfmetcorr_ey + NT_pfmetcorr_ex*NT_pfmetcorr_ex));
	return sqrt(2*(v1v2 - (lep_px*NT_pfmetcorr_ex + lep_py*NT_pfmetcorr_ey)));
	}

