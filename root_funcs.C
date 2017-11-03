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

const double b_Medium_WP = 0.8484;

// the exact LorentzVector declaration
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

float MTlep_c(float lep_px, float lep_py,
	float met_px, float met_py
	)
	{
	double v1v2 = sqrt((lep_px*lep_px + lep_py*lep_py)*(met_px*met_px + met_py*met_py));
	return sqrt(2*(v1v2 - (lep_px*met_px + lep_py*met_py)));
	}

double transverse_mass(LorentzVector v1, LorentzVector v2)                                                                                                                       
        {
        return sqrt(2*v1.pt()*v2.pt()*(1 - cos(v1.phi() - v2.phi())));
        }


