// -*- C++ -*-
//
// Package:    UserCode/NtuplerAnalyzer
// Class:      NtuplerAnalyzer
// 
/**\class NtuplerAnalyzer NtuplerAnalyzer.cc UserCode/NtuplerAnalyzer/plugins/NtuplerAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Oleksii Toldaiev
//         Created:  Thu, 13 Jul 2017 00:08:32 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// for TFileService -- to get to the output file
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// the vertex fitter for tau secondary vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"//VertexCollection

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "PhysicsTools/HepMCCandAlgos/interface/PDFWeightsHelper.h"

#include "TRandom3.h"

// MET recoil corrections for DY and WJets, from higgs->tautau group
// usage: https://github.com/CMS-HTT/RecoilCorrections/blob/master/instructions.txt
//#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
// do correction off-line in processing

// lepton ID/Iso prescriptions
#include "UserCode/llvv_fwk/interface/PatUtils.h"
#include "UserCode/llvv_fwk/interface/MacroUtils.h"
// lumiUtils
#include "UserCode/llvv_fwk/interface/LumiUtils.h"
// couple functions processing leptons
#include "UserCode/ttbar-leptons-80X/interface/ProcessingMuons.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingElectrons.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingTaus.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingJets.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingBJets.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingHLT.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingGenParticles.h"
#include "UserCode/ttbar-leptons-80X/interface/ProcessingDRCleaning.h"

static double pu_vector_NOMINAL[] = {0, 0.360609416811339, 0.910848525427002, 1.20629960507795, 0.965997726573782, 1.10708082813183, 1.14843491548622, 0.786526251164482, 0.490577792661333, 0.740680941110478,
0.884048630953726, 0.964813189764159, 1.07045369167689, 1.12497267309738, 1.17367530613108, 1.20239808206413, 1.20815108390021, 1.20049333094509, 1.18284686347315, 1.14408796655615,
1.0962284704313, 1.06549162803223, 1.05151011089581, 1.05159666626121, 1.05064452078328, 1.0491726301522, 1.05772537082991, 1.07279673875566, 1.0837536468865, 1.09536667397119,
1.10934472980173, 1.09375894592864, 1.08263679568271, 1.04289345879947, 0.9851490341672, 0.909983816540809, 0.821346330143864, 0.71704523475871, 0.609800913869359, 0.502935245638477,
0.405579825620816, 0.309696044611377, 0.228191137503131, 0.163380359253309, 0.113368437957202, 0.0772279997453792, 0.0508111733313502, 0.0319007262683943, 0.0200879459309245, 0.0122753366005436,
0.00739933885813127, 0.00437426967257811, 0.00260473545284139, 0.00157047254226743, 0.000969500595715493, 0.000733193118123283, 0.000669817107713128, 0.000728548958604492, 0.000934559691182011, 0.00133719688378802,
0.00186652283903214, 0.00314422244976771, 0.00406954793369611, 0.00467888840511915, 0.00505224284441512, 0.00562827194936864, 0.0055889504870752, 0.00522867039470319, 0.00450752163476433, 0.00395300774604375,
0.00330577167682956, 0.00308353042577215, 0.00277846504893301, 0.00223943190687725, 0.00196650068765464, 0.00184742734258922, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0};



// try simpler first
struct gen_matching {
	//Int_t sum;
	Int_t closest;
	Float_t dR;
};

struct gen_matching match_to_gen(const LorentzVector& p4,
	vector<LorentzVector>& gen_leps,
	vector<LorentzVector>& gen_taus,
	vector<LorentzVector>& gen_tau3ch,
	vector<LorentzVector>& gen_w_prods,
	vector<LorentzVector>& gen_b_prods)
{
Float_t min_dR = 9999.;
Float_t dR_cut = 0.4;
Int_t   min_id = 0;
edm::LogInfo ("Demo") << "act gen sizes " << gen_leps.size() << gen_taus.size() << gen_tau3ch.size() << gen_w_prods.size() << gen_b_prods.size();
edm::LogInfo ("Demo") << "act obj p4    " << p4.pt() << p4.eta() << p4.phi();
edm::LogInfo ("Demo") << "act gen pts  " << (gen_leps.size()>0 ? gen_leps[0].pt() : 0)
	<< (gen_taus.size()>0? gen_taus[0].pt():0)
	<< (gen_tau3ch.size()>0? gen_tau3ch[0].pt():0)
	<< (gen_w_prods.size()>0? gen_w_prods[0].pt():0)
	<< (gen_b_prods.size()>0? gen_b_prods[0].pt():0);
edm::LogInfo ("Demo") << "act matches";

// I need in-place quick loop with list of all inputs here, how to do it in C? -- do it later
int gen_id = 1;
for (unsigned int i = 0; i<gen_leps.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_leps[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		min_id = gen_id;
		}
	}

gen_id++;
for (unsigned int i = 0; i<gen_taus.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_taus[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		min_id = gen_id;
		}
	}

gen_id++;
for (unsigned int i = 0; i<gen_tau3ch.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_tau3ch[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		min_id = gen_id;
		}
	}

gen_id++;
for (unsigned int i = 0; i<gen_w_prods.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_w_prods[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		min_id = gen_id;
		}
	}

gen_id++;
for (unsigned int i = 0; i<gen_b_prods.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_b_prods[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		min_id = gen_id;
		}
	}

struct gen_matching match = {.closest=min_id, .dR=min_dR};
return match;
}

struct sv_pair {
	double flightLength;
	double flightLengthSignificance;
};

TRandom3 *r3 = new TRandom3();

typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> ROOT_TTree_vector3D;

struct sv_pair geometrical_SV(
	ROOT_TTree_vector3D& b_1, ROOT_TTree_vector3D& tr1,
	ROOT_TTree_vector3D& b_2, ROOT_TTree_vector3D& tr2,
	ROOT_TTree_vector3D& b_3, ROOT_TTree_vector3D& tr3
	)
	{
	//Float_t tracker_error = 0.002; // approximately systematic error on positions
	// it will cancel out with weights

	TVector3 b_vec1, b_vec2, b_vec3;
	//b_vec1.SetXYZ(-b1x, -b1y, -b1z); // 100% known that z here has giant error -- need to do something with it
	//b_vec2.SetXYZ(-b2x, -b2y, -b2z);
	//b_vec3.SetXYZ(-b3x, -b3y, -b3z);
	b_vec1.SetXYZ(b_1.X(), b_1.Y(), b_1.Z());
	b_vec2.SetXYZ(b_2.X(), b_2.Y(), b_2.Z());
	b_vec3.SetXYZ(b_3.X(), b_3.Y(), b_3.Z());

	// I need just the direction of tracks for geometry
	// thus making copy
	TVector3 t1, t2, t3;
	t1.SetXYZ(tr1.X(), tr1.Y(), tr1.Z());
	t2.SetXYZ(tr2.X(), tr2.Y(), tr2.Z());
	t3.SetXYZ(tr3.X(), tr3.Y(), tr3.Z());
	//t1.SetPtEtaPhi(v1pt, v1eta, v1phi);
	//t2.SetPtEtaPhi(v2pt, v2eta, v2phi);
	//t3.SetPtEtaPhi(v3pt, v3eta, v3phi);

	// weighted bis direction -- used in simple b SV (Friday result)
	TVector3 t_sum = t1 + t2 + t3;

	// root throws warning "zero vector can't be streched"
	// crab jobs crash with it
	// protective programming follows
	struct sv_pair sv_zeros = {0., 0.};
	if (t_sum.Mag() == 0)
		return sv_zeros;
	t_sum.SetMag(1);

	// after establishing direction of tau
	// tracks are only geometrical lines
	if (t1.Mag() == 0 || t2.Mag() == 0 || t3.Mag() == 0)
		return sv_zeros;
	t1.SetMag(1);
	t2.SetMag(1);
	t3.SetMag(1);

	//TVector3 tau = t1+t2+t3;
	//TVector3 tau;
	//tau.SetPtEtaPhi(taupt, taueta, tauphi);
	// tests show that tau direction and the sum are practically the same

	// find the "optimal direction"
	// -- direction of minimal angles betwee tracks and b-s in perpendicular plane

	// 2)
	// just shift bis direction randomly in max phi max theta deviations
	// choose best position, i.e. max sum b-track angles in transverse plane
	// thus, no Z changes
	double max_angle_sum = 0;
	TVector3 max_average = t_sum; // initial best direction is bis

	// find max phi and theta dev around bis
	double max_dPhi = 0, max_dTheta = 0;

	double dPhi = abs(t_sum.Phi() - t1.Phi());
	if (dPhi > max_dPhi) max_dPhi = dPhi;
	dPhi = abs(t_sum.Phi() - t2.Phi());
	if (dPhi > max_dPhi) max_dPhi = dPhi;
	dPhi = abs(t_sum.Phi() - t3.Phi());
	if (dPhi > max_dPhi) max_dPhi = dPhi;

	double dTheta = abs(t_sum.Theta() - t1.Theta());
	if (dTheta > max_dTheta) max_dTheta = dTheta;
	dTheta = abs(t_sum.Theta() - t2.Theta());
	if (dTheta > max_dTheta) max_dTheta = dTheta;
	dTheta = abs(t_sum.Theta() - t3.Theta());
	if (dTheta > max_dTheta) max_dTheta = dTheta;

	for (unsigned int i = 0; i<1000; i++)
		{
		//// uniform search around bis dir +- max dphi
		//double dPhi_shift   = max_dPhi * r3->Uniform() * 2 - max_dPhi;
		//double dTheta_shift = max_dTheta * r3->Uniform() * 2 - max_dTheta;
		//TVector3 direction = t_sum;
		//direction.SetPhi(t_sum.Phi() + dPhi_shift);
		//direction.SetTheta(t_sum.Theta() + dTheta_shift);

		// Gaussian + Markov walk from bis dir
		double dPhi_shift   = r3->Gaus(0, max_dPhi);
		double dTheta_shift = r3->Gaus(0, max_dTheta);
		// shift around current best (in principle I should also reduce sigma..)
		TVector3 direction = max_average;
		direction.SetPhi(max_average.Phi() + dPhi_shift);
		direction.SetTheta(max_average.Theta() + dTheta_shift);

		if (direction.Mag() == 0)
			return sv_zeros;
		direction.SetMag(1); // just in case

		// and to the direction
		// find perpendicular b-s
		TVector3 b_long1 = direction * (b_vec1.Dot(direction));
		TVector3 b_perp1 = b_vec1 - b_long1;
		TVector3 b_long2 = direction * (b_vec2.Dot(direction));
		TVector3 b_perp2 = b_vec2 - b_long2;
		TVector3 b_long3 = direction * (b_vec3.Dot(direction));
		TVector3 b_perp3 = b_vec3 - b_long3;

		// perpendicular parts of tracks
		TVector3 t1_long = direction * (t1.Dot(direction));
		TVector3 t1_perp = t1 - t1_long;
		TVector3 t2_long = direction * (t2.Dot(direction));
		TVector3 t2_perp = t2 - t2_long;
		TVector3 t3_long = direction * (t3.Dot(direction));
		TVector3 t3_perp = t3 - t3_long;

		double angle_sum = b_perp1.Angle(t1_perp) + b_perp2.Angle(t2_perp) + b_perp3.Angle(t3_perp);
		if (angle_sum > max_angle_sum)
			{
			max_angle_sum = angle_sum;
			max_average = direction;
			}
		}

	// and to optimal direction
	// find perpendicular b-s
	TVector3 b_long1 = max_average * (b_vec1.Dot(max_average));
	TVector3 b_perp1 = b_vec1 - b_long1;
	TVector3 b_long2 = max_average * (b_vec2.Dot(max_average));
	TVector3 b_perp2 = b_vec2 - b_long2;
	TVector3 b_long3 = max_average * (b_vec3.Dot(max_average));
	TVector3 b_perp3 = b_vec3 - b_long3;

	// perpendicular parts of tracks
	TVector3 t1_long = max_average * (t1.Dot(max_average));
	TVector3 t1_perp = t1 - t1_long;
	TVector3 t2_long = max_average * (t2.Dot(max_average));
	TVector3 t2_perp = t2 - t2_long;
	TVector3 t3_long = max_average * (t3.Dot(max_average));
	TVector3 t3_perp = t3 - t3_long;

	// project found b-s to perp tracks
	// in principle it should not be needed, since the direction is found to fit them together well
	// but let's try to get to simple SV best result
	if (t1_perp.Mag() == 0 || t2_perp.Mag() == 0 || t3_perp.Mag() == 0)
		return sv_zeros;
	t1_perp.SetMag(1);
	t2_perp.SetMag(1);
	t3_perp.SetMag(1);

	TVector3 b_long_perp1 = t1_perp * (b_perp1.Dot(t1_perp));
	TVector3 b_long_perp2 = t2_perp * (b_perp2.Dot(t2_perp));
	TVector3 b_long_perp3 = t3_perp * (b_perp3.Dot(t3_perp));

	// [let's try without these for now]

	/*
	TVector3 b_long_perp1 = b_perp1;
	TVector3 b_long_perp2 = b_perp2;
	TVector3 b_long_perp3 = b_perp3;
	*/


	// perpendiculars to bis direction, for reference
	// find perpendicular b-s
	TVector3 b_bis_long1 = t_sum * (b_vec1.Dot(t_sum));
	TVector3 b_bis_perp1 = b_vec1 - b_bis_long1;
	TVector3 b_bis_long2 = t_sum * (b_vec2.Dot(t_sum));
	TVector3 b_bis_perp2 = b_vec2 - b_bis_long2;
	TVector3 b_bis_long3 = t_sum * (b_vec3.Dot(t_sum));
	TVector3 b_bis_perp3 = b_vec3 - b_bis_long3;

	// perpendicular parts of tracks
	TVector3 t1_bis_long = t_sum * (t1.Dot(t_sum));
	TVector3 t1_bis_perp = t1 - t1_bis_long;
	TVector3 t2_bis_long = t_sum * (t2.Dot(t_sum));
	TVector3 t2_bis_perp = t2 - t2_bis_long;
	TVector3 t3_bis_long = t_sum * (t3.Dot(t_sum));
	TVector3 t3_bis_perp = t3 - t3_bis_long;

	// in the perp plane find b long to tracks
	// -- nope, no additional correction to b-s

	// the best point calculation
	// with just transverse b-s
	//TVector3 dV = t1 - t2;
	//TVector3 dB = b_perp1 - b_perp2;
	//double x12 = dV.Dot(dB) / dV.Mag2();
	//dV = t2 - t3;
	//dB = b_perp2 - b_perp3;
	//double x23 = dV.Dot(dB) / dV.Mag2();
	//dV = t3 - t1;
	//dB = b_perp3 - b_perp1;
	//double x31 = dV.Dot(dB) / dV.Mag2();

	// the best point calculation with projected b-s
	TVector3 dV = t1 - t2;
	TVector3 dB1 = b_long_perp1 - b_long_perp2;
	double x12 = - dV.Dot(dB1) / dV.Mag2();

	dV = t2 - t3;
	TVector3 dB2 = b_long_perp2 - b_long_perp3;
	double x23 = - dV.Dot(dB2) / dV.Mag2();

	dV = t3 - t1;
	TVector3 dB3 = b_long_perp3 - b_long_perp1;
	double x31 = - dV.Dot(dB3) / dV.Mag2();

	TVector3 bp12 = b_long_perp1 + x12 * t1;
	TVector3 bp21 = b_long_perp2 + x12 * t2;
	TVector3 bp_1 = 0.5*(bp12 + bp21);

	TVector3 bp23 = b_long_perp2 + x23 * t2;
	TVector3 bp32 = b_long_perp3 + x23 * t3;
	TVector3 bp_2 = 0.5*(bp23 + bp32);

	TVector3 bp31 = b_long_perp3 + x31 * t3;
	TVector3 bp13 = b_long_perp1 + x31 * t1;
	TVector3 bp_3 = 0.5*(bp31 + bp13);

	TVector3 bp_average = 0.3333*(bp_1 + bp_2 + bp_3);
	TVector3 bp_dev1 = bp_1 - bp_average;
	TVector3 bp_dev2 = bp_2 - bp_average;
	TVector3 bp_dev3 = bp_3 - bp_average;

	/*
	// and systematic error of tracker
	double syst12 = tracker_error / t1.Angle(t2); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst23 = tracker_error / t2.Angle(t3); // technically / Sin (or Tan), but Sin = Angle with these angles
	double syst31 = tracker_error / t3.Angle(t1); // technically / Sin (or Tan), but Sin = Angle with these angles
	//double syst = pow(syst12, 2) + pow(syst23, 2) + pow(syst31, 2);
	double syst12_weight = 1/syst12;
	double syst23_weight = 1/syst23;
	double syst31_weight = 1/syst31;

	// not weighted averages
	double x_average = (x12 + x23 + x31) / 3;
	double x_deviation = (pow(x12 - x_average, 2) + pow(x23 - x_average, 2) + pow(x31 - x_average, 2)) * 0.3333;
	//double x_dev_syst = x_deviation + syst;

	// weighted average with tracker errors
	//double x_average = (x12*syst12_weight + x23*syst23_weight + x31*syst31_weight) / (syst12_weight + syst23_weight + syst31_weight);
	//double x_deviation = (syst12_weight*pow(x12 - x_average, 2) + syst23_weight*pow(x23 - x_average, 2) + syst31_weight*pow(x31 - x_average, 2))/(2*(syst12_weight + syst23_weight + syst31_weight)/3);
	*/

	double convergence_factor = 1;

	//double triang_a = 0, triang_b = 0;
	//const double tan60 = 1.732, sin60 = 0.866;

	double conv1 = (bp12 - bp21).Mag();
	double conv2 = (bp23 - bp32).Mag();
	double conv3 = (bp31 - bp13).Mag();
	double conv_frac1 = conv1/dB1.Mag();
	double conv_frac2 = conv2/dB2.Mag();
	double conv_frac3 = conv3/dB3.Mag();
	double conv_frac_sum = conv_frac1 + conv_frac2 + conv_frac3;
	//double conv_frac_averaged = sqrt(pow(conv_frac1, 2) + pow(conv_frac2, 2) + pow(conv_frac3, 2));

	// SV out of all penalties
	double flightLength = 0, flightLengthSignificance = 0;
	// by relative convergence volume
	convergence_factor *= 1 / (1 + conv_frac1 * conv_frac2 * conv_frac3 / 0.027); // 0.027 = 0.3*0.3*0.3 -- when fractions are equal
	// by fraction sum-s, extracting correlation of divergences
	convergence_factor *= 1 / (1 + (conv_frac1/conv_frac_sum) * (conv_frac2/conv_frac_sum) * (conv_frac3/conv_frac_sum));

	// sign of flight
	if (bp_average.Dot(t_sum) > 0)
		{
		flightLength = bp_average.Mag();
		flightLengthSignificance = flightLength * convergence_factor / sqrt(pow(bp_dev1.Mag(), 2) + pow(bp_dev2.Mag(), 2) + pow(bp_dev3.Mag(), 2));
		}
	else
		{
		flightLength = - bp_average.Mag();
		flightLengthSignificance = flightLength * convergence_factor / sqrt(pow(bp_dev1.Mag(), 2) + pow(bp_dev2.Mag(), 2) + pow(bp_dev3.Mag(), 2));
		}

	struct sv_pair SV = {.flightLength = flightLength, .flightLengthSignificance = flightLengthSignificance};
	return SV;
	}

// ptocedures for initializations
namespace utils
	{
	namespace cmssw
		{
		// TODO: it is the same jetCorrector as in MacroUtils, only Fall_ prefix is set
		// Fall15_25nsV2_
		FactorizedJetCorrector* getJetCorrector(TString baseDir, TString pf, bool isMC)
			{
			gSystem->ExpandPathName(baseDir);
			//TString pf(isMC ? "MC" : "DATA");
			// TString pf("Fall15_25nsV2_");
			//pf += (isMC ? "MC" : "DATA");

			//order matters: L1 -> L2 -> L3 (-> Residuals)
			std::vector<std::string> jetCorFiles;
			std::cout << baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt" << std::endl;
			jetCorFiles.push_back((baseDir+"/"+pf+"_L1FastJet_AK4PFchs.txt").Data());
			jetCorFiles.push_back((baseDir+"/"+pf+"_L2Relative_AK4PFchs.txt").Data());
			jetCorFiles.push_back((baseDir+"/"+pf+"_L3Absolute_AK4PFchs.txt").Data());
			if(!isMC) jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
			// now there is a practically empty file Fall15_25nsV2_MC_L2L3Residual_AK4PFchs.txt
			// adding the run on it anyway
			//jetCorFiles.push_back((baseDir+"/"+pf+"_L2L3Residual_AK4PFchs.txt").Data());
			// it is dummy/empty file for MC and apparently is is not used
			// but in v13.1 it seemed to influence selection a bit
			// adding it for v13.4 -- will test later without it
			// and removing in 13.7 test -- compare with 13.4 & 13.4_repeat

			//init the parameters for correction
			std::vector<JetCorrectorParameters> corSteps;
			for(size_t i=0; i<jetCorFiles.size(); i++) corSteps.push_back(JetCorrectorParameters(jetCorFiles[i]));
			//return the corrector
			return new FactorizedJetCorrector(corSteps);
			}
		}
	}


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class NtuplerAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit NtuplerAnalyzer(const edm::ParameterSet&);
      ~NtuplerAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      unsigned int minTracks_;

	/* the new, never documented in Workbook stuff with tikens
	 * delivered by this guy:
	 * https://twiki.cern.ch/twiki/bin/view/Main/HanDle
	 */
	/*declare in the class */
	//edm::EDGetTokenT<reco::TrackCollection> tracks_;
	edm::EDGetTokenT<pat::MuonCollection> muons_;
	edm::EDGetTokenT<pat::ElectronCollection> electrons_;
	edm::EDGetTokenT<pat::TauCollection> taus_;
	edm::EDGetTokenT<reco::VertexCollection> vtx_;
	edm::EDGetTokenT<double> rho_;
	edm::EDGetTokenT<edm::TriggerResults> trigResults_, trigResultsRECO_, trigResultsPAT_;
	//triggerObjectsHandle.getByLabel(ev, "selectedPatTrigger");
	//iEvent.getByToken(triggerObjects_, triggerObjectsHandle);
	edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> triggerObjects_;
	edm::EDGetTokenT<LHEEventProduct> lheEPToken_;
	edm::EDGetTokenT< std::vector < PileupSummaryInfo > > puInfo_;
	edm::EDGetTokenT< std::vector < PileupSummaryInfo > > puInfo2_;
	edm::EDGetTokenT<reco::GenParticleCollection> genParticle_;
	edm::EDGetTokenT<GenEventInfoProduct> evt_;
	PDFWeightsHelper pdfweightshelper_;

	edm::EDGetTokenT<edm::View<pat::PackedCandidate>> tracks_;
	edm::EDGetTokenT<reco::BeamSpot> beamSpot_;

	edm::EDGetTokenT<pat::METCollection> mets_slimmedMETs_, mets_slimmedMETsMuEGClean_, mets_uncorrected_;

	edm::EDGetTokenT<vector<reco::GenJet>> genJets_;
	//genJetsHandle.getByLabel(ev, "slimmedGenJets");
	edm::EDGetTokenT<pat::JetCollection> jets_;
	//jetsHandle.getByLabel(ev, "slimmedJets");
	edm::EDGetTokenT<std::vector<reco::GenJet>  > genJetsToken_; // there is a twist why

	//edm::EDGetTokenT<edm::ValueMap<float> > petersonFragToken_;
	edm::EDGetTokenT<edm::ValueMap<float> > upFragToken_, centralFragToken_, downFragToken_, PetersonFragToken_, semilepbrUpToken_, semilepbrDownToken_;

	//RecoilCorrector* recoilPFMetCorrector;
	//TH2D* zPtMass_histo;

	bool record_tauID, record_tauIDantiIso, record_bPreselection, record_MonitorHLT, record_ElMu, record_Dilep, record_jets;

	TString dtag;
	bool isMC, aMCatNLO, isWJets, isDY, isTT;
	bool isLocal;
	bool withHLT;
	string  HLT_source,
		muHLT_MC1  , muHLT_MC2  ,
		muHLT_Data1, muHLT_Data2,
		elHLT_Data , elHLT_MC,
		lepMonitorHLT;

	jet_id    jetID;
	pu_jet_id jetPUID;

	TString jecDir;
	TString TjetResolutionFileName;
	TString TjetResolutionSFFileName;

	FactorizedJetCorrector *jesCor;
	JetCorrectionUncertainty *totalJESUnc;
	JME::JetResolution jet_resolution_in_pt;
	JME::JetResolutionScaleFactor jet_resolution_sf_per_eta;

	double  el_kino_cuts_pt, el_kino_cuts_eta, el_veto_kino_cuts_pt, el_veto_kino_cuts_eta,
		mu_kino_cuts_pt, mu_kino_cuts_eta, mu_veto_kino_cuts_pt, mu_veto_kino_cuts_eta,
		tau_kino_cuts_pt, tau_kino_cuts_eta;
	double jet_kino_cuts_pt, jet_kino_cuts_eta;
	double btag_threshold;

	lumiUtils::GoodLumiFilter goodLumiFilter;
	std::vector < std::string > urls; // = runProcess.getUntrackedParameter < std::vector < std::string > >("input");
	TString outUrl;

	edm::EDGetTokenT<bool> BadChCandFilterToken_;
	edm::EDGetTokenT<bool> BadPFMuonFilterToken_;

	// random numbers for corrections & uncertainties
	TRandom3 *r3;

	// gotta be a new option for definition of the interface
	//#include "ntupleOutput_leps.h"
	// so, I need to
	//  - declare the interface as class atributes
	//  - then make the TTree and connect branches in constructor
	//  - and reset/fill stuff in analyze
	// let's do it first manually for p4-s of leptons
	TTree* NT_output_ttree; 
	TH1D *event_counter, *weight_counter; 
	/*
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > NT_lep_p4;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* pt_lep_p4; // yep, vectors of complex objects require additional persistent pointers
	*/
	#define NTUPLE_INTERFACE_CLASS_DECLARE
	#include "UserCode/NtuplerAnalyzer/interface/ntupleOutput.h"
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
NtuplerAnalyzer::NtuplerAnalyzer(const edm::ParameterSet& iConfig) :
record_tauID         (iConfig.getParameter<bool>("record_tauID"))         ,
record_tauIDantiIso  (iConfig.getParameter<bool>("record_tauIDantiIso"))  ,
record_bPreselection (iConfig.getParameter<bool>("record_bPreselection")) ,
record_MonitorHLT    (iConfig.getParameter<bool>("record_MonitorHLT"))    ,
record_ElMu          (iConfig.getParameter<bool>("record_ElMu"))          ,
record_Dilep         (iConfig.getParameter<bool>("record_Dilep"))         ,
record_jets          (iConfig.getParameter<bool>("record_jets"))         ,
dtag       (iConfig.getParameter<std::string>("dtag")),
isMC       (iConfig.getParameter<bool>("isMC")),
isLocal    (iConfig.getParameter<bool>("isLocal")),
withHLT    (iConfig.getParameter<bool>("withHLT")),
HLT_source (iConfig.getParameter<std::string>("HLT_source")),
muHLT_MC1  (iConfig.getParameter<std::string>("muHLT_MC1")),
muHLT_MC2  (iConfig.getParameter<std::string>("muHLT_MC2")),
muHLT_Data1(iConfig.getParameter<std::string>("muHLT_Data1")),
muHLT_Data2(iConfig.getParameter<std::string>("muHLT_Data2")),
elHLT_Data (iConfig.getParameter<std::string>("elHLT_Data")),
elHLT_MC   (iConfig.getParameter<std::string>("elHLT_MC")),
lepMonitorHLT   (iConfig.getParameter<std::string>("lepMonitorHLT")),
jecDir     (iConfig.getParameter<std::string>("jecDir")),
TjetResolutionFileName     (iConfig.getParameter<std::string>("resolutionFile")),
TjetResolutionSFFileName   (iConfig.getParameter<std::string>("scaleFactorFile")),
el_kino_cuts_pt    (iConfig.getParameter<double>("el_kino_cuts_pt")),
el_kino_cuts_eta   (iConfig.getParameter<double>("el_kino_cuts_eta")),
el_veto_kino_cuts_pt    (iConfig.getParameter<double>("el_veto_kino_cuts_pt")),
el_veto_kino_cuts_eta   (iConfig.getParameter<double>("el_veto_kino_cuts_eta")),
mu_kino_cuts_pt    (iConfig.getParameter<double>("mu_kino_cuts_pt")),
mu_kino_cuts_eta   (iConfig.getParameter<double>("mu_kino_cuts_eta")),
mu_veto_kino_cuts_pt    (iConfig.getParameter<double>("mu_veto_kino_cuts_pt")),
mu_veto_kino_cuts_eta   (iConfig.getParameter<double>("mu_veto_kino_cuts_eta")),
tau_kino_cuts_pt    (iConfig.getParameter<double>("tau_kino_cuts_pt")),
tau_kino_cuts_eta   (iConfig.getParameter<double>("tau_kino_cuts_eta")),
jet_kino_cuts_pt    (iConfig.getParameter<double>("jet_kino_cuts_pt")),
jet_kino_cuts_eta   (iConfig.getParameter<double>("jet_kino_cuts_eta")),
btag_threshold   (iConfig.getParameter<double>("btag_threshold")),
goodLumiFilter   (iConfig.getUntrackedParameter<std::vector<edm::LuminosityBlockRange>>("lumisToProcess", std::vector<edm::LuminosityBlockRange>())),
urls   (iConfig.getUntrackedParameter <std::vector <std::string> >("input")),
outUrl (iConfig.getParameter<std::string>("outfile"))
//BadChCandFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"))),
//BadPFMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuonFilter")))
// here it breaks with
// MissingParameter: Parameter 'BadChargedCandidateFilter' not found.

{
	r3 = new TRandom3();

	/* define in constructor via call to consumes (magic thingy) */
	//tracks_    = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
	muons_     = consumes<pat::MuonCollection>    (edm::InputTag("slimmedMuons"));
	electrons_ = consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"));
	taus_ = consumes<pat::TauCollection>(edm::InputTag("slimmedTaus"));
	vtx_ = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));
	rho_ = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
	// declare consuming the HLT to be able to get triggers in the following
	if (withHLT)
		trigResults_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","", HLT_source));
	//trigResults2_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT2"));
	trigResultsRECO_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","RECO"));
	trigResultsPAT_     = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","PAT"));
	triggerObjects_ = consumes<vector<pat::TriggerObjectStandAlone>>(edm::InputTag("selectedPatTrigger"));
	//fwlite::Handle<double> rhoHandle;
	//rhoHandle.getByLabel(ev, "fixedGridRhoFastjetAll");
	lheEPToken_ = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
	puInfo_  = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));
	puInfo2_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"));
	genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));

	evt_  = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

	edm::FileInPath mc2hessianCSV = iConfig.getParameter<edm::FileInPath>("mc2hessianCSV");
	//pdfweightshelper_.Init(nPdfWeights_,nPdfEigWeights_,mc2hessianCSV);
	pdfweightshelper_.Init(100, 60, mc2hessianCSV);

	tracks_ = consumes<edm::View<pat::PackedCandidate>> (edm::InputTag("packedPFCandidates"));
	beamSpot_ = consumes<reco::BeamSpot> (edm::InputTag("offlineBeamSpot"));

	mets_slimmedMETs_ = consumes<pat::METCollection>(edm::InputTag("slimmedMETs"));
	mets_slimmedMETsMuEGClean_ = consumes<pat::METCollection>(edm::InputTag("slimmedMETsMuEGClean"));
	mets_uncorrected_ = consumes<pat::METCollection>(edm::InputTag("slimmedMETsUncorrected"));

	genJets_ = consumes<vector<reco::GenJet>>(edm::InputTag("slimmedGenJets"));
	genJetsToken_ = consumes<std::vector<reco::GenJet> >(edm::InputTag("particleLevel:jets")); // not sure if it is different from previous one
	jets_    = consumes<pat::JetCollection>(edm::InputTag("slimmedJets"));

	//petersonFragToken_(consumes<edm::ValueMap<float> >(edm::InputTag("bfragWgtProducer:PetersonFrag"))),
	upFragToken_        = consumes<edm::ValueMap<float> >(edm::InputTag("bfragWgtProducer:upFrag"));
	centralFragToken_   = consumes<edm::ValueMap<float> >(edm::InputTag("bfragWgtProducer:centralFrag"));
	downFragToken_      = consumes<edm::ValueMap<float> >(edm::InputTag("bfragWgtProducer:downFrag"));
	PetersonFragToken_  = consumes<edm::ValueMap<float> >(edm::InputTag("bfragWgtProducer:PetersonFrag"));
	semilepbrUpToken_   = consumes<edm::ValueMap<float> >(edm::InputTag("bfragWgtProducer:semilepbrUp"));
	semilepbrDownToken_ = consumes<edm::ValueMap<float> >(edm::InputTag("bfragWgtProducer:semilepbrDown"));


	//BadChCandFilterToken_ (consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter")));
	//BadPFMuonFilterToken_ (consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuonFilter")));
	// these break with no match for call to '(edm::EDGetTokenT<bool>) (edm::EDGetTokenT<bool>)'
	//BadChCandFilterToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"));
	//BadPFMuonFilterToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuonFilter"));
	// still MissingParameter: Parameter 'BadChargedCandidateFilter' not found.
	// this thing worked, not sure is it is correct:
	//BadChCandFilterToken_ = consumes<bool>(edm::InputTag("BadChargedCandidate"));
	//BadPFMuonFilterToken_ = consumes<bool>(edm::InputTag("BadPFMuon"));
	// these are not found
	BadPFMuonFilterToken_ = consumes<bool>(edm::InputTag("BadPFMuonFilter"));
	BadChCandFilterToken_ = consumes<bool>(edm::InputTag("BadChargedCandidateFilter"));
	/* try one of these strings:
	 * "BadParticleFilter",
	 *  PFCandidates  = cms.InputTag("particleFlow"),   # Collection to test
	 *  muons  = cms.InputTag("muons"),   # Collection to test
	 *  taggingMode   = cms.bool(False),
	 *  filterType  =cms.string("BadChargedCandidate"
	 *  
	 *  and for muons:
		BadPFMuonFilter = cms.EDFilter(
		    "BadParticleFilter",
		    PFCandidates  = cms.InputTag("particleFlow"),   # Collection to test
		    muons  = cms.InputTag("muons"),   # Collection to test 
		    taggingMode   = cms.bool(False),
		    filterType  =cms.string("BadPFMuon"),
		    maxDR         = cms.double(0.001),              # Maximum DR between reco::muon->innerTrack and pfCandidate 
		...
		BadChargedCandidateFilter = cms.EDFilter(
		    "BadParticleFilter",
		    PFCandidates  = cms.InputTag("particleFlow"),   # Collection to test
		    muons  = cms.InputTag("muons"),   # Collection to test
		    taggingMode   = cms.bool(False),
		    filterType  =cms.string("BadChargedCandidate"),
		    maxDR         = cms.double(0.00001),              # Maximum DR between reco::muon->innerTrack and pfCandidate 
		    minPtDiffRel = cms.double(0.00001),               # lower threshold on difference between pt of reco::muon->innerTrack and pfCandidate
		    minMuonTrackRelErr = cms.double(2.0),          # minimum ptError/pt on muon best track
		    innerTrackRelErr   = cms.double(1.0),          # minimum relPtErr on innerTrack
		    minMuonPt     = cms.double(100.0),               # minimum muon pt 
		    segmentCompatibility = cms.double(0.3),        # compatibility between the inner track and the segments in the muon spectrometer
		)
	 *
	 * -- there is no filter = True option!!
	 */

	// dtag configs
	bool period_BCD = !isMC && (dtag.Contains("2016B") || dtag.Contains("2016C") || dtag.Contains("2016D"));
	bool period_EF  = !isMC && (dtag.Contains("2016E") || dtag.Contains("2016F"));
	bool period_G   = !isMC && (dtag.Contains("2016G"));
	bool period_H   = !isMC && (dtag.Contains("2016H"));

	aMCatNLO = dtag.Contains("amcatnlo");
	isWJets = dtag.Contains("WJet") || dtag.Contains("W0Jet") || dtag.Contains("W1Jet") || dtag.Contains("W2Jet") || dtag.Contains("W3Jet") || dtag.Contains("W4Jet");
	isDY = dtag.Contains("DYJet");
	isTT = dtag.Contains("TT");

	/* do it offline
	// recoil corrector
	if (isDY || isWJets)
		{
		//TString recoil_corrections_data_file("${CMSSW_BASE}/src/HTT-utilities/RecoilCorrections/data/TypeIPFMET_2016BCD.root");
		TString recoil_corrections_data_file("/HTT-utilities/RecoilCorrections/data/TypeIPFMET_2016BCD.root");
		gSystem->ExpandPathName(recoil_corrections_data_file);
		recoilPFMetCorrector = new RecoilCorrector(recoil_corrections_data_file);
		}
	*/

	/*
	 * let's find the weight in processing
	 * also do the same in met corrector
	 */
	//if (isDY)
	//	{
	//	TString zPtMass_filename("${CMSSW_BASE}/src/HTT-utilities/RecoilCorrections/data/TypeIPFMET_2016BCD.root");
	//	gSystem->ExpandPathName(zPtMass_filename);
	//	zPtMass_histo = ;
	//	}


	// jet IDs, corrections, resolutions etc
	jetID = LooseJET; // TODO: move to Conf
	jetPUID = LoosePUJET;

	// JEC, JES, JER
	gSystem->ExpandPathName (jecDir);
	// v1
	// getJetCorrector(TString baseDir, TString pf, bool isMC)
	//https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite
	
	// in 2016 the corrections for data are per-period:
	// Summer16_23Sep2016BCDV4_DATA_        Summer16_23Sep2016EFV4_DATA_        Summer16_23Sep2016GV4_DATA_        Summer16_23Sep2016HV4_DATA_
	TString jet_corr_files;
	if (isMC)
		jet_corr_files = "/Summer16_23Sep2016V4_MC";
	else if (period_BCD)
		jet_corr_files = "/Summer16_23Sep2016BCDV4_DATA";
	else if (period_EF)
		jet_corr_files = "/Summer16_23Sep2016EFV4_DATA";
	else if (period_G)
		jet_corr_files = "/Summer16_23Sep2016GV4_DATA";
	else if (period_H)
		jet_corr_files = "/Summer16_23Sep2016HV4_DATA";
	jesCor = utils::cmssw::getJetCorrector (jecDir, jet_corr_files, isMC);
	totalJESUnc = new JetCorrectionUncertainty ((jecDir + jet_corr_files + "_Uncertainty_AK4PFchs.txt").Data());

	// resolution and scale-factors for the systematics
	gSystem->ExpandPathName(TjetResolutionFileName);
	gSystem->ExpandPathName(TjetResolutionSFFileName);

	string jetResolutionFileName   (TjetResolutionFileName);
	string jetResolutionSFFileName (TjetResolutionSFFileName);
	// <---- ROOT & CMSSW are best friends
	jet_resolution_in_pt = JME::JetResolution(jetResolutionFileName);
	jet_resolution_sf_per_eta = JME::JetResolutionScaleFactor(jetResolutionSFFileName);

	edm::Service<TFileService> fs;
	// init ttree
	NT_output_ttree = fs->make<TTree>("reduced_ttree", "TTree with reduced event data");
	event_counter  = fs->make<TH1D>( "events_counter"  , "pass category", 100,  0, 100);
	weight_counter = fs->make<TH1D>( "weight_counter"  , "pass category", 100,  0, 100); // for control of effect from PU reweighting, aMCatNLO gen -1 weights, top pt
	// in principle it should be orthogonal to the ntuple preselection
	// and be possible to do it after the preselection
	// N_presel / weighted N_presel = N_all / weighted N_all

	// connect the branch
	/*
	// set the additional pointer:
	pt_lep_p4 = &NT_lep_p4;
	NT_output_ttree->Branch("lep_p4", "vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &pt_lep_p4);
	// should be ok
	*/

	// connect the branch with macro:
	#undef NTUPLE_INTERFACE_CLASS_DECLARE
	#define NTUPLE_INTERFACE_CLASS_INITIALIZE
	#include "UserCode/NtuplerAnalyzer/interface/ntupleOutput.h"

	//now do what ever initialization is needed
	usesResource("TFileService");

	/* output via EDM stuff
	 * does it work ok with skipping events?
	 * test the straight forward way, if doesn't work -- leave it, save in TTree manualy
	 *
	 * error: 'produces' was not declared in this scope
	 * -- so it's really just for producers
	 */
	//produces<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >>("leps_p4");

}


NtuplerAnalyzer::~NtuplerAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
NtuplerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	LogInfo ("Demo") << "entered event " << iEvent.eventAuxiliary().run() << ' ' << iEvent.eventAuxiliary().luminosityBlock();

	// reset the output objects with macro
	#undef NTUPLE_INTERFACE_CLASS_DECLARE
	#undef NTUPLE_INTERFACE_CLASS_INITIALIZE
	#define NTUPLE_INTERFACE_CLASS_RESET
	#include "UserCode/NtuplerAnalyzer/interface/ntupleOutput.h"
	// if output contains stand-alone objects (not vector of LorentxVector-s, but just 1 LorentzVector, like MET or something)
	// you have to reset them yourself, since each object might have its' own method
	NT_met_init.SetXYZT(0,0,0,0);
	NT_met_uncorrected.SetXYZT(0,0,0,0);
	NT_met_corrected.SetXYZT(0,0,0,0);
	NT_met_slimmedMets.SetXYZT(0,0,0,0);
	NT_met_slimmedMETsMuEGClean.SetXYZT(0,0,0,0);

	NT_jets_full_correction.SetXYZT(0,0,0,0);

	NT_gen_t_w1_final_p4.SetXYZT(0,0,0,0);
	NT_gen_t_w2_final_p4.SetXYZT(0,0,0,0);
	NT_gen_t_b_final_p4.SetXYZT(0,0,0,0);
	NT_gen_tb_w1_final_p4.SetXYZT(0,0,0,0);
	NT_gen_tb_w2_final_p4.SetXYZT(0,0,0,0);
	NT_gen_tb_b_final_p4.SetXYZT(0,0,0,0);

	math::Error<3>::type pvCov;
	pvCov(0,0) = 999;
	pvCov(1,1) = 999;
	pvCov(2,2) = 999;
	NT_PV_cov = pvCov;


	unsigned int event_checkpoint = 0;
	event_counter->Fill(event_checkpoint++);
	double weight = 1;

	NT_indexevents = iEvent.id().event();
	NT_runNumber   = iEvent.id().run();
	NT_lumi        = iEvent.luminosityBlock();


	// for matching reco to gen I need
	// leptons (11, 13), taus (15) and tau od DM10 (3ch) from hard processes or from tt decay in case of tt
	// and visible products of W and b in case of tt
	// then a reco object will be matched to gen collections
	// the match defined by not overlapping sum:
	// lep tau tau3ch b    W
	// 1   2   4      8   16
	//vector<const reco::Candidate*> gen_leps, gen_taus, gen_tau3ch, gen_w_prods, gen_b_prods;
	vector<LorentzVector> gen_leps, gen_taus, gen_tau3ch, gen_w_prods, gen_b_prods;

	// couple things for MC:
	//  - gen aMCatNLO
	//  - gen nvtx
	//  - gen NUP (NUmber of Particles? needed for WNJets)
	//  - top pt-s, channel ID and other processing gen particles (save LorentzVector of generated taus, or non-neutrino part of generated tau)
	if(isMC)
		{
		LogInfo ("Demo") << "Processing MC";
		double weight_TopPT = 1, weight_Gen = 1;
		// ----------------------- aMCatNLO -1 weights
		//fwlite::Handle<GenEventInfoProduct> evt;
		//evt.getByLabel(ev, "generator");
		edm::Handle<GenEventInfoProduct> evt;
		iEvent.getByToken(evt_, evt);
		if(evt.isValid())
			{
			//if (debug) cout << "evt is valid, evt->weight() = " << evt->weight() << "\n";
			LogInfo("Demo") << "evt is valid, evt->weight() = " << evt->weight();
			weight_Gen = (evt->weight() > 0 ) ? 1. : -1. ;
			NT_aMCatNLO_weight = evt->weight();
			}

		// ----------------------- gen nvtx
		int ngenITpu = 0;
		edm::Handle < std::vector<PileupSummaryInfo>> puInfoH;
		//puInfoH.getByLabel (ev, "slimmedAddPileupInfo");
		iEvent.getByToken(puInfo_, puInfoH);
		if (!puInfoH.isValid())
			{
			//puInfoH.getByLabel( ev, "addPileupInfo" );
			iEvent.getByToken(puInfo2_, puInfoH);
			if (!puInfoH.isValid()) {printf("collection PileupSummaryInfo with name slimmedAddPileupInfo or addPileupInfo does not exist\n"); exit(0);}
			}
		// so here we have valid puInfoH
		// otherwise exit was called
		for (std::vector < PileupSummaryInfo >::const_iterator it = puInfoH->begin (); it != puInfoH->end (); it++)
			{
			//if (it->getBunchCrossing () == 0) ngenITpu += it->getPU_NumInteractions ();
			// guys and Mara use getTrueNumInteractions :
			if (it->getBunchCrossing () == 0) ngenITpu += it->getTrueNumInteractions();
			}
		NT_nvtx_gen = ngenITpu;
		LogInfo("Demo") << "nvtx gen is stored";

		// ----------------------- gen NUP and n of PUP
		edm::Handle < LHEEventProduct > lheEPHandle;
		//lheEPHandle.getByLabel (ev, "externalLHEProducer");
		iEvent.getByToken(lheEPToken_, lheEPHandle);
		if (isMC && lheEPHandle.isValid())
			{
			NT_gen_NUP = lheEPHandle->hepeup().NUP;
			NT_gen_n_PUP = lheEPHandle->hepeup().PUP.size();
			}

		LogInfo("Demo") << "NUP, PUP";

		if (isTT)
			{
			// I don't see this nominal weight and others in ST_tW powheg-pythia CUETP8M1 sample
			// but DY, WJets, QCD and s/t-channel all pass well
			// before figuring put what's going on I lave the waits only for TT

			// the howto is from

			// what he says about
			// > For such samples (and good practice in general),
			// > must
			// > use weight from
			// > GenEventInfoProduct::weight() to fill histograms, compute yields, train
			// > MVA’s, etc
			// I save as NT_aMCatNLO_weight
			//double nomlheweight = lheEPHandle->weights()[0].wgt; // the norm weight
			// it's not the norm weight
			// and on 1/100th of TT it always = 1.
			double nomlheweight    = lheEPHandle->originalXWGTUP();
			NT_gen_weight_norm = nomlheweight;
			LogInfo("Demo") << "PDFs, alphaS and nominal weight";

			// scale weights
			// got howto from
			// https://indico.cern.ch/event/494682/contributions/1172505/attachments/1223578/1800218/mcaod-Feb15-2016.pdf
			// printouts of weight labels should be in the git
			double  muf_nom_mur_nom_weight  = (fabs(lheEPHandle->weights()[0].wgt))/(fabs(nomlheweight));
			double   muf_up_mur_nom_weight  = (fabs(lheEPHandle->weights()[1].wgt))/(fabs(nomlheweight));
			double muf_down_mur_nom_weight  = (fabs(lheEPHandle->weights()[2].wgt))/(fabs(nomlheweight));
			double  muf_nom_mur_up_weight   = (fabs(lheEPHandle->weights()[3].wgt))/(fabs(nomlheweight));
			double   muf_up_mur_up_weight   = (fabs(lheEPHandle->weights()[4].wgt))/(fabs(nomlheweight));
			double muf_down_mur_up_weight   = (fabs(lheEPHandle->weights()[5].wgt))/(fabs(nomlheweight));
			double  muf_nom_mur_down_weight = (fabs(lheEPHandle->weights()[6].wgt))/(fabs(nomlheweight));
			double   muf_up_mur_down_weight = (fabs(lheEPHandle->weights()[7].wgt))/(fabs(nomlheweight));
			double muf_down_mur_down_weight = (fabs(lheEPHandle->weights()[8].wgt))/(fabs(nomlheweight));

			NT_gen_weights_renorm_fact.push_back(muf_down_mur_down_weight);
			NT_gen_weights_renorm_fact.push_back(muf_down_mur_nom_weight );
			NT_gen_weights_renorm_fact.push_back(muf_down_mur_up_weight  );
			NT_gen_weights_renorm_fact.push_back(muf_nom_mur_down_weight );
			NT_gen_weights_renorm_fact.push_back(muf_nom_mur_nom_weight  );
			NT_gen_weights_renorm_fact.push_back(muf_nom_mur_up_weight   );
			NT_gen_weights_renorm_fact.push_back(muf_up_mur_down_weight  );
			NT_gen_weights_renorm_fact.push_back(muf_up_mur_nom_weight   );
			NT_gen_weights_renorm_fact.push_back(muf_up_mur_up_weight    );
			LogInfo("Demo") << "scale weights";


			// and PDF uncertainties
			// also from https://indico.cern.ch/event/494682/contributions/1172505/attachments/1223578/1800218/mcaod-Feb15-2016.pdf
			// the labels in printouts correspond to one of NN
			// from the presentation:
			// > LHAID 260000, pdf uncertainties:  260001-260100,
			// > as uncertainties 265000,266000
			// issues:
			// 1) not sure if they do linearly map to the weights() array (assume they do)
			// 2) 13100 -- what are these?
			// from twiki:
			// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#PDF_uncertainties
			// talks only about NNPDF3
			// but there is an error in NNPDF3 which interferes with converting them to hessian form
			// and makes the systematic variations inconvenient,
			// which was noted by Till Arndt, and also that CT14 PDFs are stored in hessian form
			// https://indico.cern.ch/event/560833/contributions/2268181/attachments/1319955/1982168/Arndt_PDF.pdf
			// ...
			// did I just assume that the following weights are PDF CT14?
			// they have IDs [3001, 3060] and PDF set label [13100, 13156]+{13164, 13166, 11000}
			// looks like it could be 57 PDFs, 2 alphaS variations and nominal weight?
			// ...
			// aha!
			// no, I didn't guess, I noticed this bit here:
			// https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW#How_to_use_weights
			// You can see that the CT10 weight is the last in the truncated list (ID: 3001) and that there are other 109 weights before.
			// So in this case, CT10 will be weight number 109 (numbering starts from 0 as in C++) and ID = "3001". 
			// and the alphaS were assumed,
			// but I didn't save the last one 3060 with label 11000
			// they also suggest this:
			// > In the event loop, do: 
			// > int whichWeight = XXX;
			// > theWeight *= EvtHandle->weights()[whichWeight].wgt/EvtHandle->originalXWGTUP(); 
			// -- that's way the envelope weights are divided by it
			// and in PDF example it was also used?..
			// aha!
			// these are PDF label IDs:
			// https://lhapdf.hepforge.org/pdfsets.html
			// -- the guess on alphaS was right!
			//
			// last question:
			// should I divide them by that nominal weight?
			// aha! the pdf example is
			// PhysicsTools/HepMCCandAlgos/plugins/PDFWeightsTest.cc
			// -- I saved it in info-examples-notes/
			// it divides the weights by the value os weights()[0]
			// which is strange
			// -- let's not divide them now and see later?
			// the twiki page does suggest to divide everything though...
			// -- keeping the division after all, seems like a good guess

			//get the original mc replica weights
			//unsigned int nNNPPDF3Weights_ = 100, pdfWeightOffset_ = 9, nPdfEigWeights_ = 60; // as in example
			// -- for NNPDF3
			unsigned int nPdfWeights_ = 57, pdfWeightOffset_ = 9 + 100 + 2; // already, nPdfEigWeights_ = 60; // as in example
			// -- for CT pdf
			// pdf offset skips the renorm/fact weights
			// there are 100 replicas of NNPDF3.0
			// convert them into hessian matrix represenation of a 60-dim function (which covers all the correlations/deviations of the pdf uncretainty)
			// varying the hessian parameters we get systematic on the shape and everything
			//std::vector<double> inpdfweights(nPdfWeights_);
			for (unsigned int ipdf=0; ipdf<nPdfWeights_; ++ipdf)
				{
				unsigned int iwgt = ipdf + pdfWeightOffset_;
				//this is the weight to be used for evaluating uncertainties with mc replica weights
				//pdfweights_[ipdf] = lheInfo->weights()[iwgt].wgt*weight_/nomlheweight;
				//this is the raw weight to be fed to the mc2hessian convertor
				double wgtval = lheEPHandle->weights()[iwgt].wgt;
				//inpdfweights[ipdf] = wgtval;
				// for CT pdf hessians are stored
				NT_gen_weights_pdf_hessians.push_back(wgtval/nomlheweight);
				//NT_gen_weights_pdf_hessians.push_back(wgtval); // raw value
				}

			/* the CT pdfs are already hessian
			std::vector<double> outpdfweights(nPdfEigWeights_);
			//do the actual conversion, where the nominal lhe weight is needed as the reference point for the linearization
			pdfweightshelper_.DoMC2Hessian(nomlheweight,inpdfweights.data(),outpdfweights.data());

			for (unsigned int iwgt=0; iwgt<nPdfEigWeights_; ++iwgt)
				{
				double wgtval = outpdfweights[iwgt];
				//the is the weight to be used for evaluating uncertainties with hessian weights
				//pdfeigweights_[iwgt] = wgtval*weight_/nomlheweight;
				//NT_pdf_hessians.push_back(wgtval*NT_aMCatNLO_weight/nomlheweight);
				// it seem like an outdated formula to multiply by this weight
				// I'll add the multiplication in the processing if needed
				NT_gen_weights_pdf_hessians.push_back(wgtval/nomlheweight);
				}    
			*/


			// and alpha strong variation
			NT_gen_weight_alphas_1 = lheEPHandle->weights()[pdfWeightOffset_+ nPdfWeights_ + 0].wgt;
			NT_gen_weight_alphas_2 = lheEPHandle->weights()[pdfWeightOffset_+ nPdfWeights_ + 1].wgt;
			// some additional weight with label 11000, stored within the same index range 3001-3060
			NT_gen_weight_too      = lheEPHandle->weights()[pdfWeightOffset_+ nPdfWeights_ + 3].wgt;

			LogInfo("Demo") << "MC systematic weights";
			}

		// fragmentation and decay tables (of b->hadron) systematics
		edm::Handle<std::vector<reco::GenJet>> genJets2;
		iEvent.getByToken( genJetsToken_, genJets2); // https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer
		// https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer/blob/master/plugins/BFragmentationWeightProducer.cc
		// TODO: check the genJet info checks

		edm::Handle<edm::ValueMap<float> > centralFrag;
		iEvent.getByToken(centralFragToken_, centralFrag);
		edm::Handle<edm::ValueMap<float> > upFrag;
		iEvent.getByToken(upFragToken_, upFrag);
		edm::Handle<edm::ValueMap<float> > downFrag;
		iEvent.getByToken(downFragToken_, downFrag);
		edm::Handle<edm::ValueMap<float> > petersonFrag;
		iEvent.getByToken(PetersonFragToken_, petersonFrag);

		edm::Handle<edm::ValueMap<float> > semilepbrUp;
		iEvent.getByToken(semilepbrUpToken_, semilepbrUp);
		edm::Handle<edm::ValueMap<float> > semilepbrDown;
		iEvent.getByToken(semilepbrDownToken_, semilepbrDown);

		double weight_upFrag = 1, weight_centralFrag = 1, weight_downFrag = 1, weight_PetersonFrag = 1, weight_semilepbrUp = 1, weight_semilepbrDown = 1;
		for(auto genJet=genJets2->begin(); genJet!=genJets2->end(); ++genJet)
			{
			edm::Ref<std::vector<reco::GenJet> > genJetRef(genJets2, genJet - genJets2->begin()); // this looks really weird
			//cout << "pt=" << genJet->pt() << " id=" << genJet->pdgId() << " petersonFragWeight=" << (*petersonFrag)[genJetRef] << endl;
			//...
			weight_centralFrag   *= (*centralFrag)[genJetRef];
			weight_upFrag        *= (*upFrag)[genJetRef];
			weight_downFrag      *= (*downFrag)[genJetRef];
			weight_PetersonFrag  *= (*petersonFrag)[genJetRef];
			weight_semilepbrUp   *= (*semilepbrUp)[genJetRef];
			weight_semilepbrDown *= (*semilepbrDown)[genJetRef];
			}
		NT_gen_weight_centralFrag   = weight_centralFrag  ;
		NT_gen_weight_FragUp        = weight_upFrag       ;
		NT_gen_weight_FragDown      = weight_downFrag     ;
		NT_gen_weight_PetersonFrag  = weight_PetersonFrag ;
		NT_gen_weight_semilepbrUp   = weight_semilepbrUp  ;
		NT_gen_weight_semilepbrDown = weight_semilepbrDown;


		LogInfo("Demo") << "MC systematic weights for jet fragmentation";

		LogInfo ("Demo") << "Processing MC, gen particles";
		// ----------------------- GENERATED PARTICLES
		// parse gen particles tree and get top pt-s and channel
		// channels are needed for:
		// TTbar, Single-top, DY
		// for now implement thed only TTbar and Single-Top
		reco::GenParticleCollection gen;
		edm::Handle<reco::GenParticleCollection> genHandle;
		//genHandle.getByLabel(ev, "prunedGenParticles");
		iEvent.getByToken(genParticle_, genHandle);

		NT_gen_genPx = 0, NT_gen_genPy = 0, NT_gen_visPx = 0, NT_gen_visPy = 0;
		if(genHandle.isValid())
			{
			gen = *genHandle;
			// For reference, some PDG IDs:
			// QUARKS
			// d  1
			// u  2
			// s  3
			// c  4
			// b  5
			// t  6
			// b' 7
			// t' 8
			// g 21
			// gamma 22
			// Z     23
			// W     24
			// h     25
			// e, ve     11, 12
			// mu, vmu   13, 14
			// tau, vtau 15, 16

        		vector<const reco::Candidate*> t_b_parts, tb_b_parts, t_W1_parts, tb_W1_parts, t_W2_parts, tb_W2_parts;
			LogInfo ("Demo") << "Processing MC, gen particles, t decays and taus";
			NT_gen_pythia8_prompt_leptons_N = 0;
			NT_gen_N_wdecays = 0;
			NT_gen_N_zdecays = 0;
			LorentzVector genMomentum(0,0,0,0);
			NT_genPt = 0;
			NT_genMass = 0;
			//NT_zPtWeight = 1;
			for(size_t i = 0; i < gen.size(); ++ i)	
				{
				const reco::GenParticle & p = gen[i];
				//const reco::Candidate* p_cand = GenParticle inherits from candidate, let's try to just use it
				int id = p.pdgId();
				unsigned int a_id = abs(id);
				int st = p.status(); // TODO: check what is status in decat simulation (pythia for our TTbar set)
				int n_daughters = p.numberOfDaughters();

				// leptons from hard processes
				if ((a_id >= 11 && a_id <= 16 && p.fromHardProcessFinalState()) ||
					//(p.isDirectHardProcessTauDecayProduct()))
					(p.isDirectHardProcessTauDecayProductFinalState())) // same stuff
					{
					// saving the hard leptons from here if it is not TT
					if (!isTT)
						{
						if (a_id == 11 || a_id == 13)
							save_final_cands(&p, gen_leps);
						// if it is tau -- check if it is DM10+, i.e. decay to 3 charged particles
						if (a_id == 15)
							{
							int tau_id = simple_tau_decay_id(&p);
							// = 11, 13 for leptons and 20 + 5*(Nch-1) + Npi0 for hadrons
							// 20 + 10
							if (abs(tau_id) >= 30)
								save_final_cands(&p, gen_tau3ch);
							else
								save_final_cands(&p, gen_taus);
							}
						}
					// Save parameters for recoil corrections
					// relevant for DY and WJets
					if (isDY || isWJets)
						{
						if (isDY) genMomentum += p.p4();

						NT_gen_genPx += p.p4().Px();
						NT_gen_genPy += p.p4().Py();

						if ( !(a_id == 12 || a_id == 14 || a_id == 16) )
							{
							NT_gen_visPx += p.p4().Px();
							NT_gen_visPy += p.p4().Py();
							}
						}
					}

				if (abs(id) == 6 && n_daughters == 2)
					// if it is a t quark
					// it is a decay vertex of t to something
					// (could use p.isLastCopy())
					{
					// calculate top_pt weights:
					weight_TopPT *= TMath::Sqrt(top_pT_SF(p.pt()));
					// find the W decay channel in this top
					unsigned int d0_id = abs(p.daughter(0)->pdgId());
					unsigned int d1_id = abs(p.daughter(1)->pdgId());
					int W_num = d0_id == 24 ? 0 : (d1_id == 24 ? 1 : -1) ;
					if (W_num < 0) continue;
					const reco::Candidate * W = p.daughter( W_num );
					const reco::Candidate * b = p.daughter( 1 - W_num );
					const reco::Candidate * W_final = find_W_decay(W);
					int decay_id = 1;
					// = id of lepton or 1 for quarks
					for (unsigned int d_i = 0; d_i < 2; d_i++)
						{
						int d_i_pdgId = W_final->daughter(d_i)->pdgId();
						if (fabs(d_i_pdgId) == 11 || fabs(d_i_pdgId) == 13)
							{
							decay_id = d_i_pdgId;
							save_final_cands(W_final->daughter(d_i), gen_leps);
							}
						if (fabs(d_i_pdgId) == 15)
							{
							int tau_id = simple_tau_decay_id(W_final->daughter(d_i));
							// = 11, 13 for leptons and 20 + 5*(Nch-1) + Npi0 for hadrons
							// 20 + 10
							if (abs(tau_id) >= 30)
								save_final_cands(W_final->daughter(d_i), gen_tau3ch);
							else
								save_final_cands(W_final->daughter(d_i), gen_taus);
							decay_id = d_i_pdgId * abs(tau_id); // to keep the sign of the tau take it in abs
							}
						}
					// if W is not leptonic, and decay id is still = 1
					if (decay_id == 1)
						{
						save_final_cands(W_final, gen_w_prods);
						}
					save_final_cands(b, gen_b_prods);

					// save stuff, according to top Id, also save top p_T
					if (id>0)
						{
						NT_gen_t_pt  = p.pt();
						NT_gen_t_w_decay_id = decay_id;
						save_final_states(W_final->daughter(0), NT_gen_t_w1_final_p4s, NT_gen_t_w1_final_pdgIds, NT_gen_t_w1_final_statuses, t_W1_parts);
						save_final_states(W_final->daughter(1), NT_gen_t_w2_final_p4s, NT_gen_t_w2_final_pdgIds, NT_gen_t_w2_final_statuses, t_W2_parts);
						save_final_states(b, NT_gen_t_b_final_p4s, NT_gen_t_b_final_pdgIds, NT_gen_t_b_final_statuses, t_b_parts);
						for (unsigned int i=0; i<NT_gen_t_w1_final_p4s.size(); i++)
							{
							// skip neutrinos
							if (abs(NT_gen_t_w1_final_pdgIds[i]) == 12 || abs(NT_gen_t_w1_final_pdgIds[i]) == 14 || abs(NT_gen_t_w1_final_pdgIds[i]) == 16) continue;
							NT_gen_t_w1_final_p4 += NT_gen_t_w1_final_p4s[i];
							}
						for (unsigned int i=0; i<NT_gen_t_w2_final_p4s.size(); i++)
							{
							// skip neutrinos
							if (abs(NT_gen_t_w2_final_pdgIds[i]) == 12 || abs(NT_gen_t_w2_final_pdgIds[i]) == 14 || abs(NT_gen_t_w2_final_pdgIds[i]) == 16) continue;
							NT_gen_t_w2_final_p4 += NT_gen_t_w2_final_p4s[i];
							}
						for (unsigned int i=0; i<NT_gen_t_b_final_p4s.size(); i++)
							{
							// skip neutrinos
							if (abs(NT_gen_t_b_final_pdgIds[i]) == 12 || abs(NT_gen_t_b_final_pdgIds[i]) == 14 || abs(NT_gen_t_b_final_pdgIds[i]) == 16) continue;
							NT_gen_t_b_final_p4 += NT_gen_t_b_final_p4s[i];
							}
						}

					else
						{
						NT_gen_tb_pt = p.pt();
						NT_gen_tb_w_decay_id = decay_id;
						save_final_states(W_final->daughter(0), NT_gen_tb_w1_final_p4s, NT_gen_tb_w1_final_pdgIds, NT_gen_tb_w1_final_statuses, tb_W1_parts);
						save_final_states(W_final->daughter(1), NT_gen_tb_w2_final_p4s, NT_gen_tb_w2_final_pdgIds, NT_gen_tb_w2_final_statuses, tb_W2_parts);
						save_final_states(b, NT_gen_tb_b_final_p4s, NT_gen_tb_b_final_pdgIds, NT_gen_tb_b_final_statuses, tb_b_parts);
						for (unsigned int i=0; i<NT_gen_tb_w1_final_p4s.size(); i++)
							{
							// skip neutrinos
							if (abs(NT_gen_tb_w1_final_pdgIds[i]) == 12 || abs(NT_gen_tb_w1_final_pdgIds[i]) == 14 || abs(NT_gen_tb_w1_final_pdgIds[i]) == 16) continue;
							NT_gen_tb_w1_final_p4 += NT_gen_tb_w1_final_p4s[i];
							}
						for (unsigned int i=0; i<NT_gen_tb_w2_final_p4s.size(); i++)
							{
							// skip neutrinos
							if (abs(NT_gen_tb_w2_final_pdgIds[i]) == 12 || abs(NT_gen_tb_w2_final_pdgIds[i]) == 14 || abs(NT_gen_tb_w2_final_pdgIds[i]) == 16) continue;
							NT_gen_tb_w2_final_p4 += NT_gen_tb_w2_final_p4s[i];
							}
						for (unsigned int i=0; i<NT_gen_tb_b_final_p4s.size(); i++)
							{
							// skip neutrinos
							if (abs(NT_gen_tb_b_final_pdgIds[i]) == 12 || abs(NT_gen_tb_b_final_pdgIds[i]) == 14 || abs(NT_gen_tb_b_final_pdgIds[i]) == 16) continue;
							NT_gen_tb_b_final_p4 += NT_gen_tb_b_final_p4s[i];
							}
						}
					}

				// if it is a tau -- save the non-neutrino part to the output
				//  the status is 1 or 2
				//  1. final state, not decays, so it should never happen for tau
				//  2. decayed or fragmented -- the case for tau
				if (abs(id) == 15 && st == 1)
					NT_gen_tau_p4.push_back(p.p4()); 
				else if (abs(id) == 15 && st == 2)
					{
					// it's a final state tau
					// select its' daughters, skipping neutrinos
					// add their momenta -- use the sum as a visible_gen_tau
					LorentzVector vis_ds(0,0,0,0);
					LogInfo ("Demo") << "N tau daughters: " << n_daughters;
					for (int j = 0; j < n_daughters; ++j)
						{
						const reco::Candidate * d = p.daughter(j);
						unsigned int d_id = abs(d->pdgId());
						LogInfo ("Demo") << j << " tau daughter ID = " << d->pdgId();
						if (d_id == 12 || d_id == 14 || d_id == 16) continue;
						vis_ds += d->p4();
						}
					NT_gen_tau_p4.push_back(vis_ds); 
					}

				// Prompt leptons ID for Z->LL
				// pythia 8 stores prompt particles with status 21-29 ("hardest subprocess", PYTHIA 8 Worksheet for tutorial at ASP 2012 Summer School)
				// -- checked DYJetsToLL -- there is no Z (pdgId 23) particles, but prompt leptons work fine
				// thus save N prompt leptons in the process
				// and their ID
				if ((abs(id) == 11 || abs(id) == 13 || abs(id) == 15) && st > 20 && st < 30)
					{
					NT_gen_pythia8_prompt_leptons_N += 1;
					int gen_prompt_lepton_ID = id;
					if (abs(id) == 15)
						gen_prompt_lepton_ID *= simple_tau_decay_id(&p);
					NT_gen_pythia8_prompt_leptons_IDs.push_back(gen_prompt_lepton_ID);
					LogInfo ("Demo") << "Found (pythia8) prompt lepton: " << id << ' ' << gen_prompt_lepton_ID << ' ' << NT_gen_pythia8_prompt_leptons_N;
					}

				// but madgraph DY (50-Inf, i.e. the main one) has the Z-s............
				if (abs(id) == 23 && n_daughters == 2)
					{
					NT_gen_N_zdecays += 1;
					int d0_id = p.daughter(0)->pdgId();
					int d1_id = p.daughter(1)->pdgId();
					int a_d0_id = abs(d0_id);
					int a_d1_id = abs(d1_id);
					int lep_daughter = (a_d0_id == 11 || a_d0_id == 13 || a_d0_id == 15 ? 0 : (a_d1_id == 11 || a_d1_id == 13 || a_d1_id == 15 ? 1 : -1));
					if (lep_daughter >= 0)
						{
						if (a_d0_id == 15) d0_id *= simple_tau_decay_id(p.daughter(0));
						if (a_d1_id == 15) d1_id *= simple_tau_decay_id(p.daughter(1));
						}
					NT_gen_zdecays_IDs.push_back(d0_id);
					NT_gen_zdecays_IDs.push_back(d1_id);
					}

				// and W->Lnu processes, apparently prompt leptons don't work there -- sometimes lepton goes directly to final state
				// search for first W decay -- it's supposedly prompt
				// could merge this with t decay procedure, or search from the final state leptons..
				if (abs(id) == 24 && n_daughters == 2)
					{
					int wdecay_id = 1;
					int d0_id = abs(p.daughter(0)->pdgId());
					int d1_id = abs(p.daughter(1)->pdgId());
					int lep_daughter = (d0_id == 11 || d0_id == 13 || d0_id == 15 ? 0 : (d1_id == 11 || d1_id == 13 || d1_id == 15 ? 1 : -1));
					if (lep_daughter >= 0)
						{
						wdecay_id = p.daughter(lep_daughter)->pdgId();
						if (abs(wdecay_id) == 15)
							wdecay_id *= simple_tau_decay_id(p.daughter(lep_daughter));
						}

					NT_gen_N_wdecays += 1;
					NT_gen_wdecays_IDs.push_back(wdecay_id);
					}
				}

			if (isDY)
				{
				NT_genPt = genMomentum.Pt();
				NT_genMass = genMomentum.M();
				// weight in processing
				//NT_zPtWeight = zPtMass_histo->GetBinContent(zPtMass_histo->GetXaxis()->FindBin(NT_genMass), zPtMass_histo->GetYaxis()->FindBin(NT_genPt));
				}

			LogInfo ("Demo") << "Found: t decay = " << NT_gen_t_w_decay_id << " ; tb decay = " << NT_gen_tb_w_decay_id;
			}

		weight = pu_vector_NOMINAL[NT_nvtx_gen] * (aMCatNLO? weight_Gen : 1) * weight_TopPT;
		}
	weight_counter->Fill(event_checkpoint, weight);
	weight_counter->Fill(++event_checkpoint, weight * NT_aMCatNLO_weight);
	// in old times (2015) there was a recommendation to actually use this as weight -- checking this

	if (isMC)
		{
		LogInfo ("Demo") << "gen sizes " << gen_leps.size() << gen_taus.size() << gen_tau3ch.size() << gen_w_prods.size() << gen_b_prods.size();
		//LogInfo ("Demo") << "gen PDGs  " << (gen_leps.size()>0 ? gen_leps[0].pdgId() : 0)
		//	<< (gen_taus.size()>0?    gen_taus[0].pdgId():0)
		//	<< (gen_tau3ch.size()>0?  gen_tau3ch[0].pdgId():0)
		//	<< (gen_w_prods.size()>0? gen_w_prods[0].pdgId():0)
		//	<< (gen_b_prods.size()>0? gen_b_prods[0].pdgId():0);
		LogInfo ("Demo") << "gen pts  " << (gen_leps.size()>0 ? gen_leps[0].pt() : 0)
			<< (gen_taus.size()>0?    gen_taus[0].pt():0)
			<< (gen_tau3ch.size()>0?  gen_tau3ch[0].pt():0)
			<< (gen_w_prods.size()>0? gen_w_prods[0].pt():0)
			<< (gen_b_prods.size()>0? gen_b_prods[0].pt():0);
		LogInfo ("Demo") << "gen etas  " << (gen_leps.size()>0 ? gen_leps[0].eta() : 0)
			<< (gen_taus.size()>0?    gen_taus[0].eta():0)
			<< (gen_tau3ch.size()>0?  gen_tau3ch[0].eta():0)
			<< (gen_w_prods.size()>0? gen_w_prods[0].eta():0)
			<< (gen_b_prods.size()>0? gen_b_prods[0].eta():0);
		LogInfo ("Demo") << "gen phis  " << (gen_leps.size()>0 ? gen_leps[0].phi() : 0)
			<< (gen_taus.size()>0?    gen_taus[0].phi():0)
			<< (gen_tau3ch.size()>0?  gen_tau3ch[0].phi():0)
			<< (gen_w_prods.size()>0? gen_w_prods[0].phi():0)
			<< (gen_b_prods.size()>0? gen_b_prods[0].phi():0);
		}

	//Handle<reco::TrackCollection> tracks;
	//iEvent.getByToken( tracks_, tracks );
	//LogInfo("Demo") << "number of tracks "<<tracks->size();
	/*
	if( minTracks_ <= tracks->size() ) {
	   LogInfo("Demo") << "number of tracks "<<tracks->size();
	}
	*/

	// ------------------------------------------------- Apply MET FILTERS

	/* it seems the tagging mode is broken in these filters,
	 * thus trying to see if they work in just filtering mode
	 * --- the tests didn't work, all events in selected data files passed,
	 *     a couple of ~100k events files were used
	 */
	edm::Handle<bool> ifilterbadChCand;
	edm::Handle<bool> ifilterbadPFMuon;

	iEvent.getByToken(BadChCandFilterToken_, ifilterbadChCand);
	iEvent.getByToken(BadPFMuonFilterToken_, ifilterbadPFMuon);
	//bool  filterbadChCandidate = *ifilterbadChCand;
	//bool filterbadPFMuon = *ifilterbadPFMuon;
	NT_METfilterbadChCand = *ifilterbadChCand;
	NT_METfilterbadPFMuon = *ifilterbadPFMuon;

	// in https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2017
	// they say the bool is false if rejected by event

	//if (!(filterbadChCandidate && filterbadPFMuon)) return;
	//if ((filterbadChCandidate && filterbadPFMuon)) return;

	/*
	 * MET filters are data-only thing -- remove events before passing and counting lumi, since MC is then normalized to data lumi
	 * thus after passing lumi data and MC should only have the same cuts
	 *
	 * info on MET filters and their presence in MINIAOD:
	 *   https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#ETmiss_filters
	 *   https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2017
	 *     filter                      location                     data	MC(fullSim)	MC(fastSim)	comment
	 *     primary vertex filter       available in miniAOD v2	DONE	suggested	suggested	 
	 *     beam halo filter            available in miniAOD v2	DONE	suggested	not suggested	Beam Halo Presentation
	 *     HBHE noise filter	   available in miniAOD v2	DONE	suggested	suggested	HCAL DPG Presentation
	 *     HBHEiso noise filter	   available in miniAOD v2	DONE	suggested	suggested	same as above
	 *     ECAL TP filter              available in miniAOD v2	DONE	suggested	suggested	ECAL DPG Presentation
	 *     Bad PF Muon Filter          to be run on the fly 	DONE	suggested	suggested	PPD presentation
	 *     Bad Charged Hadron Filter   to be run on the fly 	DONE	suggested	suggested	PPD presentation
	 *     ee badSC noise filter       available in miniAOD v2	DONE	not suggested	not suggested
	 *   https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
	 *   https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes#MET_Recipes
	 *
	 *   https://twiki.cern.ch/twiki/bin/view/CMS/MissingET
	 *   and their hypernews:
	 *   https://hypernews.cern.ch/HyperNews/CMS/get/met.html
	 */

	//edm::TriggerResultsByName recoFilters = iEvent.triggerResultsByName("RECO"); //is present only if PAT (and miniAOD) is not run simultaniously with RECO
	edm::TriggerResultsByName patFilters  = iEvent.triggerResultsByName("PAT"); //is present only if PAT (and miniAOD) is not run simultaniously with RECO
	//if(!isMC && !metFilters.isValid()){metFilters = iEvent.triggerResultsByName("PAT");} //if not present, then it's part of RECO
	//if(!isMC && !metFilters.isValid()){       
	//	LogInfo("Demo") << "TriggerResultsByName for MET filters is not found in the process, as a consequence the MET filter is disabled for this event";
	//	return;
	//	}

	//if (!isMC && metFilters.isValid())
	// apparently MET POG suggests trying and looking at filters in MC
	// and they are there, in patFilters
	if (patFilters.isValid())
		{
		// event is good if all filters ar true
		NT_filters_hbhe             = utils::passTriggerPatterns(patFilters, "Flag_HBHENoiseFilter*", "Flag_HBHENoiseIsoFilter*");
		NT_filters_ecalDeadCellTrig = utils::passTriggerPatterns(patFilters, "Flag_EcalDeadCellTriggerPrimitiveFilter*");
		NT_filters_good_vertices    = utils::passTriggerPatterns(patFilters, "Flag_goodVertices");
		NT_filters_eebad            = utils::passTriggerPatterns(patFilters, "Flag_eeBadScFilter");
		NT_filters_halo             = utils::passTriggerPatterns(patFilters, "Flag_globalTightHalo2016Filter");
		NT_filters_halo_super       = utils::passTriggerPatterns(patFilters, "Flag_globalSuperTightHalo2016Filter");
		// 2016 thing: bad muons
		//bool flag_noBadMuons = utils::passTriggerPatterns(patFilters, "Flag_noBadMuons"); // <---- the bad muons are done on the fly with cfg.py thingy
		//bool flag_duplicateMuons = utils::passTriggerPatterns(patFilters, "Flag_duplicateMuons");
		// from
		// https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes#Event_flags
		// Three flags are saved in the event:
		//    Flag_badMuons: the event contained at least one PF muon of pT > 20 GeV that is flagged as bad
		//    Flag_duplicateMuons: the event contained at least one PF muon of pT > 20 GeV that is flagged as duplicate
		//    Flag_noBadMuons: the event does not contain any PF muon of pT > 20 GeV flagged as bad or duplicate (i.e. the event is safe)

		// --- thus the Flag_noBadMuons should be enough
		///
		// --- there is no such flag in Feb03 rereco
		//     they also don't exist in 07Aug rereco,
		//     but the Flag_BadPFMuon and the ChHadron do exist
		// --- they are in PAT filters, not in RECO
		//     apparently PAT should be used first -- it seems to have all needed flags
		//     and RECO if something is missing...
		//     also PAT is filled with stuff for MC
		//     the "on the fly" stuff I didn't manage to get working with their receipts
		//     and their test also doesn't filter anything in 100k events of SingleMuon RunB
		//     ...

		//if (! (filters1 & good_vertices & eebad & halo & flag_noBadMuons)) return;
		// save for study
		NT_pass_basic_METfilters = NT_filters_hbhe & NT_filters_ecalDeadCellTrig & NT_filters_good_vertices & NT_filters_eebad & NT_filters_halo;
		// these Flag_noBadMuons/Flag_duplicateMuons are MET flags (the issue with bad muons in 2016),
		// they are true if the MET got corrected and event is fine

		// 
		// add: BadChHadron and BadPFMuon -- it seems their name should be Flag_BadChHadron etc
		//
		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2017
		// Bad PF Muon Filter	to be run on the fly	
		// -- "run on the fly", no this flag in the data itself
		//
		// but at the same time:
		//
		// Note that with the in the re-miniaod you will have (will rerun as pointed out below for) the following flags for the "bad muon" events:
		//    Bad PF Muon Filter
		//    Bad Charged Hadrons
		//    Flag_badMuons -> New Giovanni's Filter that the MET is corrected for (flag is set to true if the MET got corrected)
		//    Flag_duplicateMuons -> New Giovanni's Filter that the MET is corrected for (flag is set to true if the MET got corrected)

		// aha https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#ETmiss_filters
		// "Note that many of the current recommended filters can be accessed directly from Miniaod using the flag stored in the TriggerResults,
		//  with the exception of Bad Charged Hadron and Bad Muon Filters."
		// --- so, 2 vs 1 that there should be no Flags for these two in MINIAOD
		//  they should be run on the fly
		///


		//
		// MET POG gives some names to their filters instead of givin the name in code
		// apparently the actual name in the code can be found at:
		// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
		//
		// and there is no BadChHandron
		// the closes to their names are:
		// BadChargedCandidateFilter BadPFMuonFilter
		//
		// -- need to print out what actually is in 03Feb ReReco & ask on hypernews.
		//
		//  found these:
		//  root [7] metFilters.triggerNames()
		//  (const std::vector<std::string> &)
		//  { "Flag_duplicateMuons", "Flag_badMuons", "Flag_noBadMuons",
		//    "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter",
		//    "Flag_CSCTightHaloFilter", "Flag_CSCTightHaloTrkMuUnvetoFilter", "Flag_CSCTightHalo2015Filter",
		//    "Flag_globalTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter",
		//    "Flag_HcalStripHaloFilter", "Flag_hcalLaserEventFilter",
		//    "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellBoundaryEnergyFilter",
		//    "Flag_goodVertices",
		//    "Flag_eeBadScFilter",
		//    "Flag_ecalLaserCorrFilter",
		//    "Flag_trkPOGFilters",
		//    "Flag_chargedHadronTrackResolutionFilter",
		//    "Flag_muonBadTrackFilter",
		//    "Flag_trkPOG_manystripclus53X", "Flag_trkPOG_toomanystripclus53X", "Flag_trkPOG_logErrorTooManyClusters",
		//    "Flag_METFilters" }
		///
		//
		// it seems the track of these two filters goes to:
		// https://indico.cern.ch/event/591506/contributions/2387636/attachments/1381281/2099935/2016_12_01_MET_Scanning_Report_PPD.pdf
		// https://twiki.cern.ch/twiki/bin/view/CMS/MissingETScanners#More_info_on_filter_bad_ChargedC
		// and back to
		// https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
		// -- Flag_BadChargedCandidateFilter
		// and Flag_BadPFMuonFilter
		}

	LogInfo ("Demo") << "passed MET filters";

	// PASS LUMI
	// done with some trick in crab/cmsRun python config
	if (!isMC && isLocal)
		{
		if (!goodLumiFilter.isGoodLumi(iEvent.eventAuxiliary().run(), iEvent.eventAuxiliary().luminosityBlock())) return; 
		}
	event_counter->Fill(event_checkpoint++);
	weight_counter->Fill(event_checkpoint, weight);


	edm::Handle<double> rhoHandle;
	iEvent.getByToken(rho_, rhoHandle);
	if(rhoHandle.isValid() ) NT_fixedGridRhoFastjetAll = *rhoHandle;
	//NT_fixedGridRhoFastjetAll = NT_fixedGridRhoFastjetAll;

	Handle<pat::MuonCollection> muons_h;
	iEvent.getByToken( muons_, muons_h );
	pat::MuonCollection muons = *muons_h;
	Handle<pat::ElectronCollection> electrons_h;
	iEvent.getByToken( electrons_, electrons_h );
	pat::ElectronCollection electrons = *electrons_h;
	Handle<pat::TauCollection> taus_h;
	iEvent.getByToken( taus_, taus_h );
	pat::TauCollection taus = *taus_h;

	//LogInfo("Demo") << "number of muons "<< muons.size();


	// HLT TRIGGER
	// for 2015 noHLT MC (like QCD)
	// in this case all triggers pass, but leptons never match to trigger objects
	bool lepMonitorTrigger = !withHLT;
	bool eTrigger = !withHLT;
	bool muTrigger1 = !withHLT, muTrigger2 = !withHLT, muTrigger = !withHLT;
	bool jetsHLT140 = !withHLT, jetsHLT400 = !withHLT, jetsHLT = !withHLT;

	// TriggerNames for TriggerObjects --------------------
	edm::Handle<edm::TriggerResults> trigResults; //our trigger result object
	//edm::InputTag * trigResultsTag; // the tag object, trigResults are extracted from the event via this tag

	string matched_elTriggerName("");
	string matched_muTriggerName1("");
	string matched_muTriggerName2("");
	//edm::TriggerResultsByName tr = ev.triggerResultsByName ("HLT");

	if (withHLT)
		{
		edm::TriggerResultsByName tr = iEvent.triggerResultsByName (HLT_source);
		//if (!tr.isValid ()){
			// HLT2 was a quirk of Spring16 MC campaigns (noHLT/reHLT/withHLT thing)
			// need to compare 2016-2015 in tau SV
			// using HLT2 as backup (DY50 with HLT is present only in reHLT campaign, wich has this HLT2 path)
			//tr = iEvent.triggerResultsByName ("HLT2");
			// it crashes weirdly with "no consumes" complaints
			//}

		if (!tr.isValid ()){
			cout << HLT_source << " is NOT valid!" << endl;
			return;
			}
		else
			{
			LogInfo("Demo") << "Trigger HLT is valid";
			//trigResultsTag = new edm::InputTag("TriggerResults","","HLT"); //make sure have correct process on MC
			// pass trigger
			// using this:
			//   bool passTriggerPatternsAndGetName(edm::TriggerResultsByName& tr, std::string& pathName, std::string pattern)
			// -- pathName is the matched part of the trigger name (as I got it)
			//    it is passed to select trigger objects
			lepMonitorTrigger   = utils::passTriggerPatterns(tr, lepMonitorHLT);
			eTrigger   = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_elTriggerName, elHLT_MC)   : utils::passTriggerPatternsAndGetName(tr, matched_elTriggerName, elHLT_Data));
			muTrigger1 = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName1, muHLT_MC1) : utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName1, muHLT_Data1));
			muTrigger2 = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName2, muHLT_MC2) : utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName2, muHLT_Data2));
			muTrigger = muTrigger1 || muTrigger2;

			jetsHLT140 = utils::passTriggerPatterns(tr, "HLT_PFJet140_v*");
		jetsHLT400 = utils::passTriggerPatterns(tr, "HLT_PFJet400_v*");
		}

		if (record_jets)
			{
			jetsHLT = jetsHLT140 || jetsHLT400;
			}
		}

	if (!(eTrigger || muTrigger || lepMonitorTrigger || jetsHLT)) return; // orthogonalization is done afterwards
	event_counter ->Fill(event_checkpoint++);
	weight_counter->Fill(event_checkpoint, weight);

	NT_HLT_el = eTrigger;
	NT_HLT_mu = muTrigger;
	NT_HLT_lepMonitor = lepMonitorTrigger;
	NT_HLT_jets140 = jetsHLT140;
	NT_HLT_jets400 = jetsHLT400;

	LogInfo ("Demo") << "passed HLT " << eTrigger << ' ' << muTrigger << '(' << muTrigger1 << ',' << muTrigger2 << ')' << ';' << matched_elTriggerName << ' ' << matched_muTriggerName1 << ',' << matched_muTriggerName2;

	// HLT matching
	// objects of our triggers
	vector<pat::TriggerObjectStandAlone> el_trig_objs;
	vector<pat::TriggerObjectStandAlone> mu_trig_objs, mu_trig_objs2;

	if (withHLT)
		{
		// names for trigger bits
		//edm::EDGetTokenT<edm::TriggerResults> trigResults_ = consumes<edm::TriggerResults>(trigResultsTag);
		//ev.getByLabel(*trigResultsTag, trigResults);
		iEvent.getByToken( trigResults_, trigResults );
		const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);

		//fwlite::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjectsHandle;
		edm::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjectsHandle;
		//triggerObjectsHandle.getByLabel(ev, "selectedPatTrigger");
		iEvent.getByToken(triggerObjects_, triggerObjectsHandle);
		if (!triggerObjectsHandle.isValid())
			{
			LogInfo("Demo") << "!triggerObjectsHandle.isValid()";
			return;
			}
		LogInfo ("Demo") << "got trigger objects";
		vector<pat::TriggerObjectStandAlone> trig_objs = *triggerObjectsHandle;

		if (eTrigger)
			{
			Processing_selectHLTobjects(trig_objs, trigNames, el_trig_objs, matched_elTriggerName);
			}
		if (muTrigger)
			{
			Processing_selectHLTobjects(trig_objs, trigNames, mu_trig_objs,  matched_muTriggerName1);
			Processing_selectHLTobjects(trig_objs, trigNames, mu_trig_objs2, matched_muTriggerName2);
			// vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
			mu_trig_objs.insert(mu_trig_objs.end(), mu_trig_objs2.begin(), mu_trig_objs2.end());
			}
		}

	LogInfo ("Demo") << "our trigger objects: " << el_trig_objs.size() << ',' << mu_trig_objs.size();

	// PRIMARY VERTEX
	reco::VertexCollection vtx;
	edm::Handle<reco::VertexCollection> vtxHandle;
	iEvent.getByToken(vtx_, vtxHandle);
	if(vtxHandle.isValid() ) vtx = *vtxHandle;
	NT_nvtx = vtx.size();

	reco::Vertex goodPV;                                                                                                                                                                                               
	unsigned int nGoodPV(0);
	for(size_t ivtx=0; ivtx<vtx.size(); ++ivtx)
		{
		//if(utils::isGoodVertex(vtx[ivtx]))
		// directly from rumors
		// some use:
		// * at least 4 degrees of freedome (ndof) (>=4) (!)
		// * Rho < 2 (impact parameter to the beam spot
		// * z < 24
		bool its_good = (!vtx[ivtx].isFake()) && vtx[ivtx].ndof() > 4 && abs(vtx[ivtx].z()) < 24 && abs(vtx[ivtx].position().Rho()) < 2;
		// it should be equivalent to patUtils procedure
		// only they use reverse: ! > 24 etc -- but without the >=, thus there is 1 bit of discrepancy
		if (its_good)
			{
			if(nGoodPV==0) goodPV=vtx[ivtx];
			nGoodPV++;
			// save info on the vertex
			NT_PV_x.push_back(vtx[ivtx].x());
			NT_PV_y.push_back(vtx[ivtx].y());
			NT_PV_z.push_back(vtx[ivtx].z());
			NT_PV_x_err.push_back(vtx[ivtx].xError());
			NT_PV_y_err.push_back(vtx[ivtx].yError());
			NT_PV_z_err.push_back(vtx[ivtx].zError());
			}
		}

	//weight = 1; // reset weights?

	// MUONS
	LorentzVector muDiff(0., 0., 0., 0.);
	unsigned int nVetoMu_Iso = 0, nVetoMu_all = 0;
	//pat::MuonCollection selIDMuons, selMuons;
	pat::MuonCollection selMuons;
	processMuons_ID_ISO_Kinematics(muons, goodPV, weight, patUtils::llvvMuonId::StdTight, patUtils::llvvMuonId::StdLoose, patUtils::llvvMuonIso::Tight, patUtils::llvvMuonIso::Loose,               
		mu_kino_cuts_pt, mu_kino_cuts_eta, mu_veto_kino_cuts_pt, mu_veto_kino_cuts_eta, selMuons, muDiff, nVetoMu_Iso, nVetoMu_all, false, false);

	//nVetoMu += processMuons_MatchHLT(selIDMuons, mu_trig_objs, 0.4, selMuons);

	// ELECTRONS
	//pat::ElectronCollection selIDElectrons, selElectrons;
	pat::ElectronCollection selElectrons;
	unsigned int nVetoE_IsoImp = 0, nVetoE_Iso = 0, nVetoE_all = 0;
	LorentzVector elDiff(0., 0., 0., 0.);
	processElectrons_ID_ISO_Kinematics(electrons, goodPV, NT_fixedGridRhoFastjetAll, weight, patUtils::llvvElecId::Tight, patUtils::llvvElecId::Loose, patUtils::llvvElecIso::Tight, patUtils::llvvElecIso::Loose,
		el_kino_cuts_pt, el_kino_cuts_eta, el_veto_kino_cuts_pt, el_veto_kino_cuts_eta, selElectrons, elDiff, nVetoE_IsoImp, nVetoE_Iso, nVetoE_all, false, false);

	//nVetoE += processElectrons_MatchHLT(selIDElectrons, el_trig_objs, 0.4, selElectrons);

	std::vector<patUtils::GenericLepton> selLeptons;
	for(size_t l=0; l<selElectrons.size(); ++l) selLeptons.push_back(patUtils::GenericLepton (selElectrons[l] ));
	for(size_t l=0; l<selMuons.size(); ++l)     selLeptons.push_back(patUtils::GenericLepton (selMuons[l]     ));
	std::sort(selLeptons.begin(), selLeptons.end(), utils::sort_CandidatesByPt);

	//LogInfo ("Demo") << "selected leptons: " << '(' << selIDElectrons.size() << ',' << selIDMuons.size() << ')' <<  selLeptons.size() << ' ' << nVetoE << ',' << nVetoMu;
	LogInfo ("Demo") << "selected leptons: " << '(' << selElectrons.size() << ',' << selMuons.size() << ')' <<  selLeptons.size() << ' ' << nVetoE_IsoImp << ',' << nVetoMu_Iso;

	//bool clean_lep_conditions = nVetoE==0 && nVetoMu==0 && nGoodPV != 0; // veto on std iso veto leptons
	//bool clean_lep_conditions = nVetoE_all==0 && nVetoMu_all==0 && nGoodPV != 0; // veto on all iso veto leptons
	bool clean_lep_conditions = nGoodPV != 0; // just good PV, the loosest req,save bit if no veto leps
	if (!(clean_lep_conditions && selLeptons.size() > 0 && selLeptons.size() < 3)) return;
	// exit now to reduce computation -- all record schemes have this requirement

	event_counter ->Fill(event_checkpoint++);
	weight_counter->Fill(event_checkpoint, weight);

	LogInfo ("Demo") << "passed lepton conditions ";

	bool leps_passed_relIso = true;
	for(size_t l=0; l<selMuons.size(); ++l)
		{
		NT_lep_p4.push_back(selMuons[l].p4());
		NT_lep_id.push_back(selMuons[l].pdgId());
		// mu_trig_objs or el_trig_objs
		NT_lep_matched_HLT.push_back(processMuon_matchesHLTs(selMuons[l], mu_trig_objs, 0.4));
		NT_lep_dz  .push_back(selMuons[l].muonBestTrack()->dz (goodPV.position()));
		NT_lep_dxy .push_back(selMuons[l].muonBestTrack()->dxy (goodPV.position()));
		NT_lep_dB.push_back(selMuons[l].dB());
		float rel_iso = relIso(selMuons[l], NT_fixedGridRhoFastjetAll);
		NT_lep_relIso.push_back(rel_iso);
		// using old procedures for now
		bool passIso = patUtils::passIso(selMuons[l], patUtils::llvvMuonIso::Tight, patUtils::CutVersion::Moriond17Cut);
		leps_passed_relIso &= passIso;
		if (isMC)
			{
			struct gen_matching match = match_to_gen(selMuons[l].p4(), gen_leps, gen_taus, gen_tau3ch, gen_w_prods, gen_b_prods);
			NT_lep_matching_gen   .push_back(match.closest);
			NT_lep_matching_gen_dR.push_back(match.dR);
			}
		}

	for(size_t l=0; l<selElectrons.size(); ++l)
		{
		NT_lep_p4.push_back(selElectrons[l].p4());
		NT_lep_id.push_back(selElectrons[l].pdgId());
		// mu_trig_objs or el_trig_objs
		NT_lep_matched_HLT.push_back(processElectron_matchesHLTs(selElectrons[l], el_trig_objs, 0.4));
		NT_lep_dz  .push_back(selElectrons[l].gsfTrack()->dz (goodPV.position()));
		NT_lep_dxy .push_back(selElectrons[l].gsfTrack()->dxy (goodPV.position()));
		NT_lep_dB.push_back(selElectrons[l].dB());
		float rel_iso = relIso(selElectrons[l], NT_fixedGridRhoFastjetAll);
		NT_lep_relIso.push_back(rel_iso);
		//bool passIso = patUtils::passIso(selMuons[l], el_ISO, patUtils::CutVersion::Moriond17Cut, rho);
		bool passIso = patUtils::passIso(selElectrons[l], patUtils::llvvElecIso::Tight, patUtils::CutVersion::Moriond17Cut, NT_fixedGridRhoFastjetAll);
		leps_passed_relIso &= passIso;
		if (isMC)
			{
			struct gen_matching match = match_to_gen(selElectrons[l].p4(), gen_leps, gen_taus, gen_tau3ch, gen_w_prods, gen_b_prods);
			NT_lep_matching_gen   .push_back(match.closest);
			NT_lep_matching_gen_dR.push_back(match.dR);
			}
		}

	NT_nleps = selLeptons.size();
	// for control
	NT_nleps_veto_el_isoimp = nVetoE_IsoImp;
	NT_nleps_veto_el_iso = nVetoE_Iso;
	NT_nleps_veto_el_all = nVetoE_all;
	NT_nleps_veto_mu_iso = nVetoMu_Iso;
	NT_nleps_veto_mu_all = nVetoMu_all;
	NT_no_std_veto_leps  = nVetoE_IsoImp == 0 && nVetoMu_Iso == 0;
	NT_no_iso_veto_leps  = nVetoE_Iso == 0 && nVetoMu_Iso == 0;

	NT_leps_ID = 1;
	for (unsigned int i = 0; i<selLeptons.size(); i++)
		{
		NT_leps_ID *= selLeptons[i].pdgId();
		}
	//NT_leps_ID = NT_leps_ID;

	LogInfo ("Demo") << "saved leptons";



	// MET
	pat::METCollection mets_slimmedMETs, mets_slimmedMETsMuEGClean;
	edm::Handle<pat::METCollection> metsHandle_slimmedMETs;
	edm::Handle<pat::METCollection> metsHandle_slimmedMETsMuEGClean;
	iEvent.getByToken(mets_slimmedMETs_, metsHandle_slimmedMETs);
	iEvent.getByToken(mets_slimmedMETsMuEGClean_, metsHandle_slimmedMETsMuEGClean);

	// 2016: slimmedMETsMuEGClean are corrected by muons and electrons, only in Data!
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes
	// slimmedMets should be there for both Data and MC
	if(metsHandle_slimmedMETs.isValid() )
		{
		const pat::MET& MET = metsHandle_slimmedMETs->front();
		NT_met_slimmedMets = MET.p4();
		NT_met_init = MET.p4();
		}
	// LorentzVector met = mets_slimmedMETs[0].p4 ();

	// these are valid only for data
	if(metsHandle_slimmedMETsMuEGClean.isValid() )
		{
		const pat::MET& MET2 = metsHandle_slimmedMETsMuEGClean->front();
		NT_met_slimmedMETsMuEGClean = MET2.p4();
		}

	// also for control let's get uncorrected met and compare the two:
	if (!isMC) // sadly this exists only in latest ReReco data made with 8.0.26 CMSSW, not in Summer16 MC
		{
		pat::METCollection mets_uncorrected;
		edm::Handle<pat::METCollection> mets_uncorrectedHandle;
		//mets_uncorrectedHandle.getByLabel(ev, "slimmedMETsUncorrected");
		iEvent.getByToken( mets_uncorrected_, mets_uncorrectedHandle);
		if(mets_uncorrectedHandle.isValid() ) mets_uncorrected = *mets_uncorrectedHandle;
		pat::MET met_uncorrected = mets_uncorrected[0];
		NT_met_uncorrected = met_uncorrected.p4();
		}


	// JETS
	/* jets additionally need initialization of:
	 * genJets
	 * jesCor, totalJESUnc,
	 * pass jet ID, PU jet ID (with/without PU),
	 * systematic variation (NOMINAL, the variation factors are saved per jet for offline)
	 * jet resolution in pt, eta (?)
	 * kinematic cuts
	 */

	// jets
	pat::JetCollection jets;
	edm::Handle<pat::JetCollection>jetsHandle;
	//jetsHandle.getByLabel(ev, "slimmedJets");
	iEvent.getByToken(jets_, jetsHandle);
	if(jetsHandle.isValid() ) jets = *jetsHandle;

	// get genJets from the event
	std::vector<reco::GenJet> genJets;
	edm::Handle<std::vector<reco::GenJet>> genJetsHandle;
	//genJetsHandle.getByLabel(ev, "slimmedGenJets");
	iEvent.getByToken( genJets_, genJetsHandle); // twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#GenJets
	if (genJetsHandle.isValid() ) genJets = *genJetsHandle;

	LorentzVector full_jet_corr(0., 0., 0., 0.);
	//pat::JetCollection IDjets;
	//pat::JetCollection selJets;
	pat::JetCollection selJetsNoLep;
	//map<systematic_shift, pat::JetCollection> IDjets;
	// it's filled with jetSystematics by processJets_CorrectJES_SmearJERnJES_ID_ISO_with_systematics
	//string jetID("Loose");
	//string jetPUID("MediumPU");
	//Variation jet_m_systematic_variation = Variation::NOMINAL;

	//processJets_CorrectJES_SmearJERnJES_ID_ISO(jets, genJets, isMC, weight, NT_fixedGridRhoFastjetAll, nGoodPV, jesCor, totalJESUnc, 0.4/2,
	//	jet_resolution_in_pt, jet_resolution_sf_per_eta, jet_m_systematic_variation, jetID, jetPUID, /*with_PU*/ false, r3, full_jet_corr, IDjets, false, false);
	// selecting jets in the following loop
	// similar to https://github.com/LLRCMS/LLRHiggsTauTau/blob/b8bc9146cab462fdaf8e4161761b5e70a08f4a65/NtupleProducer/plugins/HTauTauNtuplizer.cc

	NT_nbjets = 0;
	NT_nallbjets = 0;
	NT_njets = 0;
	NT_nalljets = 0;
	string btagger_label("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	for (unsigned int ijet=0; ijet<jets.size(); ijet++)
		{
		pat::Jet& jet = jets[ijet];
		// the only requirements are pt 111 and no-leptons 222 in dR and Loose Jet ID 333
		// (it is not a subject of systematic variations - right?)
		// maybe later remove the Loose ID requirement in ntuple for further studies
		// but I've never used it yet -- thus applying it for speed & size

		if (jet.pt() < jet_kino_cuts_pt) continue; // 111
		// all eta pass -- forward jets too, usefull for WJets control region

		//crossClean_in_dR(selJets, selLeptons, 0.4, selJetsNoLep, weight, string("selJetsNoLep"), false, false);
	        bool overlapWithLepton(false);
		float min_dR = 0.4;
	        for(unsigned int l=0; l<(unsigned int)selLeptons.size();++l)
	                {
	                if (reco::deltaR(jet, selLeptons[l])<min_dR)
	                        { overlapWithLepton=true; break; }
	                }
	        if (overlapWithLepton) continue; // 222

		// PF jet ID
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

		double abseta = fabs(jet.eta());

		bool looseJetID = false;
		bool tightJetID = false;
		bool tightLepVetoJetID = false;

		// and latest (at Moriond17) stuff:
		// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
		// quote:
		// > Note: All fractions are calculated with the raw/uncorrected energy of the jet (only then they add up to unity). So the PF JetID has to be applied before the jet energy corrections.
		// --- khmm....
		// in MINIAOD the jets are already corrected  --> one gets the uncorrected jet and reapplies the corrections
		// > ... collection of AK4 jets slimmedJets, made from ak4PFJetsCHS ... "These have standard jet energy corrections applied (L1FastJet, L2, L3), and a pT cut at 10 GeV"
		// so now one need to get the uncorrected jet --> 
		// TODO <-- it doesn't make sense: corrected jet differs only in pt, which doesn't matter in these cuts

		if (abseta <= 2.7)
			{
			looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abseta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abseta>2.4);
			tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abseta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abseta>2.4);
			tightLepVetoJetID = ((NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abseta<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abseta>2.4) );
			}
		else if (abseta <= 3.0)
			{
			looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
			tightJetID = looseJetID;
			}
		else
			{
			looseJetID = (NEMF<0.90 && NumNeutralParticles>10);
			tightJetID = looseJetID;
			}

		if (!looseJetID) continue; // 333

		NT_jet_id.push_back(jet.pdgId());

		// the initial slimmedJet is saved here
		LorentzVector jet_initial_p4 = jet.p4();
		NT_jet_initial_p4.push_back(jet_initial_p4);

		//NT_jet_rad.push_back(jet_radius(jet));
		NT_jet_etaetaMoment.push_back(jet.etaetaMoment());
		NT_jet_phiphiMoment.push_back(jet.phiphiMoment());
		NT_jet_pu_discr.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
		//NT_jet_pu_discr_updated.push_back(ijet->hasUserFloat("pileupJetIdUpdated:fullDiscriminant") ? ijet->userFloat("pileupJetIdUpdated:fullDiscriminant") : -999);
		float b_discriminator = jet.bDiscriminator(btagger_label);
		NT_jet_b_discr.push_back(b_discriminator);
		if (b_discriminator > btag_threshold) NT_nallbjets += 1;

		NT_jet_hadronFlavour.push_back(jet.hadronFlavour());
		NT_jet_partonFlavour.push_back(jet.partonFlavour());

		int jetid = 0;
		if (looseJetID)
			{
			++jetid;
			if (abseta < jet_kino_cuts_eta)
				{
				NT_njets += 1;
				if (b_discriminator > btag_threshold) NT_nbjets += 1;
				}
			// counter of our old jets: cut on pt, on eta (small eta, central jet) and PFID Loose
			}
		if (tightJetID) ++jetid;
		if (tightLepVetoJetID) ++jetid;

		NT_jet_PFID.push_back(jetid);

		// just one more parameter for info
		NT_jetCharge.push_back (jet.jetCharge());

		// corrections
		NT_jet_area.push_back (jet.jetArea());

		float jecFactor = jet.jecFactor("Uncorrected") ;
	 	//float jetRawPt = jecFactor * jet.pt();
		//float jetRawPt2 = ijet->pt() / jecFactor; // this is wrong
		NT_jet_uncorrected_jecFactor.push_back (jecFactor);
		LorentzVector rawJet = jet.correctedP4("Uncorrected");
		NT_jet_uncorrected_p4.    push_back(rawJet);
		// TODO: compare the two offline

		// save FactorizedCorrector factor
		//double toRawSF=jet.correctedJet("Uncorrected").pt()/jet.pt();
		//LorentzVector rawJet(jet*toRawSF);
		// FactorizedJetCorrector *jesCor; with all jet correction files
		jesCor->setJetEta(rawJet.eta());
		jesCor->setJetPt(rawJet.pt());
		jesCor->setJetA(jet.jetArea());
		jesCor->setRho(NT_fixedGridRhoFastjetAll);
		//jesCor->setNPV(nGoodPV); // not used in PU Jet ID example, shouldn't matter
		float jes_correction = jesCor->getCorrection();
		//jet.setP4(rawJet*jes_correction);
		//LorentzVector jesCorJet = (rawJet*jes_correction);
		//jet.addUserFloat("jes_correction", jes_correction);
		// TODO: compare jet_p4 (default MiniAOD jet) and uncorrected_p4 * jes_correction <- re-corrected jet
		NT_jet_jes_recorrection.push_back(jes_correction);
		jet.setP4(jet.p4()*jes_correction);
		// default jets are fully corrected, the initial slimmedJet is saved
		// the raw is saved too

		// jet energy scale has uncertainty
		//totalJESUnc->setJetEta(jet.eta());
		//totalJESUnc->setJetPt(jet.pt()); // should be the corrected jet pt <- de hell this means? do it after MC SF?
		// this is rel shift to jes_cor
		// so it need the corrected jet, but with new jes_cor
		totalJESUnc->setJetEta(rawJet.eta());
		totalJESUnc->setJetPt(rawJet.pt() * jes_correction);
		float relShift = fabs(totalJESUnc->getUncertainty(true));
		// use it with rawJet*jes_correction*(1 +- relShift)
		// since all these corrections are multiplication of p4
		// I can do this shift whenever I want
		// uncertainty shift is saved only for the NOMINAL jet, which is default MiniAOD one now
		NT_jet_jes_uncertainty.push_back(relShift);

		float dR_max = 0.4/2;
		double jet_resolution = -1;
		double jer_sf = -1;
		double jer_sf_up = -1;
		double jer_sf_down = -1;
		// and the final factors from these SFs
		double jer_factor = -1, jer_factor_up = -1, jer_factor_down = -1;
		// here is the matching of the jet:
		if(isMC)
			{
			LorentzVector gen_jet_p4(0,0,0,0);
			Float_t genjet_pt = -1, genjet_dR = -1;
			Int_t genjet_i = -1;
			Bool_t genjet_matched = false;

			// the JER SF and resolution for the jet:
			//std::vector<double> jer_sf_pair = JER_SF(jet.eta());
			//double jer_sf = jer_sf_pair[0];
			//double jer_resolution = jer_sf_pair[1]; // TODO: not sure about this -- is the table the same as what their tool returns?
			// getting it with the tool from files:
			jet_resolution = jet_resolution_in_pt.getResolution({{JME::Binning::JetPt, jet.pt()},
				{JME::Binning::JetEta, jet.eta()},
				{JME::Binning::Rho, NT_fixedGridRhoFastjetAll}});
			//jer_sf = resolution_sf.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, m_systematic_variation);
			jer_sf      = jet_resolution_sf_per_eta.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::NOMINAL);
			jer_sf_up   = jet_resolution_sf_per_eta.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::UP);
			jer_sf_down = jet_resolution_sf_per_eta.getScaleFactor({{JME::Binning::JetEta, jet.eta()}}, Variation::DOWN);

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
				genjet_dR = reco::deltaR(jet, genJet);

				if (genjet_dR > dR_max) continue;
				genjet_i = i;
				genjet_pt = genJet.pt();
				genjet_matched = true;

				double dPt = std::abs(genJet.pt() - jet.pt());
				double dPt_max_factor = 3*jet.pt(); // from twiki
				if (dPt > dPt_max_factor * jet_resolution) continue;

				matched_genJet = &genJet;
				}

			// calculate and apply the smearing factor to jet energy
			// with one of the two algorithms, if jet matched to genJet or not
			if (matched_genJet)
				{ // the scaling is ok
				gen_jet_p4 = matched_genJet->p4();

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
				}
			else
				{ // the smearing with gaussian randomizer

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
				}

			// THUS MC was shifted by SF jer_factor -- it got corrected and the correction has to be propagated to MET
			full_jet_corr += jet.p4() - jet_initial_p4; // initial jet + this difference = corrected jet
	
			// gen jet info only for MC
			//NT_jet_matched_genjet_p4.push_back(gen_jet_p4);
			NT_genjet_matched  .push_back( genjet_matched );
			NT_genjet_pt       .push_back( genjet_pt      );
			NT_genjet_dR       .push_back( genjet_dR      );
			NT_genjet_i        .push_back( genjet_i       );

			LogInfo ("Demo") << "match gen to jets";
			LogInfo ("Demo") << "gen sizes " << gen_leps.size() << gen_taus.size() << gen_tau3ch.size() << gen_w_prods.size() << gen_b_prods.size();
			// match to GENPRODUCTS
			struct gen_matching match = match_to_gen(jet.p4(), gen_leps, gen_taus, gen_tau3ch, gen_w_prods, gen_b_prods);
			NT_jet_matching_gen   .push_back(match.closest);
			NT_jet_matching_gen_dR.push_back(match.dR);
			}

		// the default jet is fully recorrected
		// but the corrections can be repeated offline
		NT_jet_p4.push_back(jet.p4());

		NT_jet_resolution.push_back(jet_resolution);
		NT_jet_jer_sf.        push_back(jer_sf);
		NT_jet_jer_sf_up.     push_back(jer_sf_up);
		NT_jet_jer_sf_down.   push_back(jer_sf_down);

		NT_jet_jer_factor.       push_back(jer_factor);
		NT_jet_jer_factor_up.    push_back(jer_factor_up);
		NT_jet_jer_factor_down.  push_back(jer_factor_down);

		/* just a note:
		 * so, my jets are default MiniAOD jets + MC SF for energy resolution
		 * the energy resolution UP/DOWN is saved and the energy scale uncertainty
		 *
		 * also I save JES factors: for reapplied Factorized Corrector and for uncorrected jet
		 * in principle I can jet nominal jets with reapplied Factorized Corrector, just without systematic uncertainties
		 * I save slimmedMETs[0] as main MET and propagate jet MC SF to it as met_corrected
		 * it corresponds to nominal jets
		 * there are other control mets
		 *
		 * it corresponds to LLRHiggs mets and jets
		 * and it is simple separation of different corrections and mets
		 */

		// save this jet -- it's used in tau dR match to jets
		selJetsNoLep.push_back(jet);
		/* requirements for there jets: pt and no leptons in dR
		 * no PFID requirement and no eta
		 */
		}
	NT_nalljets  = selJetsNoLep.size();

	NT_jets_full_correction = full_jet_corr; // initial jets + this correction = corrected jets
	// ALSO MET
	LorentzVector MET_corrected = NT_met_init - full_jet_corr; // checked: this is really negligible correction
	//float met_corrected = MET_corrected.pt();
	NT_met_corrected = MET_corrected;
	// TODO: I don't correct NT_met_slimmedMets and the other -- the correction should be applied offline if needed
	// in principle NT_met_init = NT_met_slimmedMets -- so NT_met_corrected saves the applied correction

	LogInfo ("Demo") << "saved jets";

	// default jets are fully corrected, the correction is propagated to this met,
	// the initial slimmedMet is saved
	// the slimmedJet, which corresponds to this met is saved too and the raw of this jet is saved
	// therefore if needed to redo met cor for only subset of selected jets
	// take slimmed met and slimmed jets
	// for slimmed jets take rawjets and the factors
	// apply them, caclulate the difference to slimmed and propagate to met

	/* apply them offline
	// APPLY RECOILD CORRECTIONS TO MET
	NT_pfmetcorr_ex = 0;
	NT_pfmetcorr_ey = 0;
	if (isDY || isWJets)
		{
		recoilPFMetCorrector->CorrectByMeanResolution(
			NT_met_corrected.Px(), // uncorrected type I pf met px (float)
			NT_met_corrected.Py(), // uncorrected type I pf met py (float)
			NT_gen_genPx, // generator Z/W/Higgs px (float)
			NT_gen_genPy, // generator Z/W/Higgs py (float)
			NT_gen_visPx, // generator visible Z/W/Higgs px (float)
			NT_gen_visPy, // generator visible Z/W/Higgs py (float)
			NT_nalljets,  // number of jets (hadronic jet multiplicity) (int) <-- they use jets with pt>30... here it's the same, only pt requirement (20), no eta or PF ID
			NT_pfmetcorr_ex, // corrected type I pf met px (float)
			NT_pfmetcorr_ey  // corrected type I pf met py (float)
			);
		LogInfo("Demo") << "recoil-corrected MET = " << NT_pfmetcorr_ex << ' ' << NT_pfmetcorr_ey;
		}
	*/

	//pat::JetCollection selJets;
	//processJets_Kinematics(IDjets, /*bool isMC,*/ weight, jet_kino_cuts_pt, jet_kino_cuts_eta, selJets, false, false);

	//pat::JetCollection selJetsNoLep;
	//crossClean_in_dR(selJets, selLeptons, 0.4, selJetsNoLep, weight, string("selJetsNoLep"), false, false);
	// and these are output jets for NTuple
	// they pass ID, corrected with JEC (smeared JES for MC)
	// pass kinematic cuts (pt, eta)
	// and dR-cleaned from selected leptons

	//std::sort (selJetsNoLep.begin(),  selJetsNoLep.end(),  utils::sort_CandidatesByPt);
	// no need for sort


	/*
	 * TAUS
	 */
	//LogInfo("Demo") << "taus.size() = "<< taus.size();
	//string tau_Loose_ID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
	string tau_VLoose_ID ("byVLooseIsolationMVArun2v1DBoldDMwLT");
	string tau_Loose_ID  ("byLooseIsolationMVArun2v1DBoldDMwLT");
	string tau_Medium_ID ("byMediumIsolationMVArun2v1DBoldDMwLT");
	string tau_Tight_ID  ("byTightIsolationMVArun2v1DBoldDMwLT");
	string tau_VTight_ID ("byVTightIsolationMVArun2v1DBoldDMwLT");
	string tau_decayMode       ("decayModeFinding");
	string tau_againstMuon     ("againstMuonTight3");
	string tau_againstElectron ("againstElectronTightMVA6");

	pat::TauCollection IDtaus, selTaus;
	processTaus_ID    (taus,   weight, tau_decayMode, tau_againstMuon, tau_againstElectron, IDtaus, false, false);
	//processTaus_ID_ISO    (taus,   weight, tau_decayMode, tau_VLoose_ID, tau_againstMuon, tau_againstElectron, IDtaus, false, false);
	processTaus_Kinematics(IDtaus, weight, tau_kino_cuts_pt, tau_kino_cuts_eta, selTaus,      false, false);

	pat::TauCollection selTausNoLep;
	crossClean_in_dR(selTaus,       selLeptons, 0.4, selTausNoLep,        weight, string("selTausNoLep"),        false, false);

	// and these are the NT output taus
	std::sort (selTausNoLep.begin(),  selTausNoLep.end(),  utils::sort_CandidatesByPt);

	// PREPARE TRACKS AND STUFF for REFITTING TAU SV AND PV

	// BEAMSPOT
	//
	reco::BeamSpot beamSpot;
	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByToken(beamSpot_, beamSpotHandle);
	if (beamSpotHandle.isValid()) beamSpot = *beamSpotHandle;

	// TRACKS
	//
	edm::Handle<edm::View<pat::PackedCandidate>> tracksHandle;
	iEvent.getByToken(tracks_, tracksHandle);
	//if (tracksHandle.isValid()) tracks = *tracksHandle;
	const edm::View<pat::PackedCandidate>* track_cands = tracksHandle.product();

	// some more feature for tracks
	// the python conf load some module TransientTrackBuilder from cmssw package TrackingTools.TransientTrack
	// probably the module stores this stuff to the event..
	edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);

	// select tracks associated with PV
	reco::TrackCollection pvTracks;
	reco::TrackCollection allTracks; // for taus (with possible SV) (testing now)

	//TLorentzVector aTrack;
	for(size_t i=0; i<track_cands->size(); ++i)
		{
		if((*track_cands)[i].charge()==0 || (*track_cands)[i].vertexRef().isNull()) continue;
		if(!(*track_cands)[i].bestTrack()) continue;

		unsigned int key = (*track_cands)[i].vertexRef().key();
		int quality = (*track_cands)[i].pvAssociationQuality();

		// here I need to select "good" tracks
		// save them to all tracks
		// and if they belong to PV save them to pv tracks
		if (!(key!=0 ||
			(quality!=pat::PackedCandidate::UsedInFitTight
			 && quality!=pat::PackedCandidate::UsedInFitLoose)))// continue;
			{
			pvTracks.push_back(*((*track_cands)[i].bestTrack()));
			// allTracks.push_back(*((*track_cands)[i].bestTrack())); // test for HelixLine Momentum is zero
			}

		// TODO: add requirement of "goodness"?
		allTracks.push_back(*((*track_cands)[i].bestTrack()));
		}

	//NT_ntaus = 0;
	NT_ntaus = selTausNoLep.size(); // all tau before MVA anti-jet iso
	std::vector<double > tracksToBeRemoved_PV; // compared by Pt due to the conflict of comparing const and not const iterators
	// these tracks correspond to all 3pi tau tracks in the event
	// without considering the ID of these taus
	// whenever refit works for a tau the corresponding tracks are removed
	for(size_t i=0; i<selTausNoLep.size(); ++i)
		{
		pat::Tau& tau = selTausNoLep[i];

		Int_t IDlev = 0;
		if (tau.tauID(tau_VTight_ID)) IDlev = 5;
		else if (tau.tauID(tau_Tight_ID))  IDlev = 4;
		else if (tau.tauID(tau_Medium_ID)) IDlev = 3;
		else if (tau.tauID(tau_Loose_ID))  IDlev = 2;
		else if (tau.tauID(tau_VLoose_ID)) IDlev = 1;

		//if (IDlev > 0) NT_ntaus += 1;

		NT_tau_id.push_back(tau.pdgId());
		NT_tau_decayMode.push_back(tau.decayMode());
		NT_tau_p4.push_back(tau.p4());
		NT_tau_IDlev.push_back(IDlev);
		NT_tau_leading_track_pt.push_back(tau.userFloat("leading_track_pt"));
		NT_tau_leadChargedHadrCand_pt.push_back(tau.userFloat("leadChargedHadrCand_pt"));
		NT_tau_leadNeutralCand_pt.push_back(tau.userFloat("leadNeutralCand_pt"));
		NT_tau_leadCand_pt.push_back(tau.userFloat("leadCand_pt"));
		NT_tau_hasSecondaryVertex.push_back(tau.hasSecondaryVertex());
		//NT_tau_hcalEnergy = tau.hcalEnergy();
		//NT_tau_hcalEnergyLeadChargedHadrCand = tau.hcalEnergyLeadChargedHadrCand();

		if (isMC)
			{
			LogInfo ("Demo") << "gen match to tau";
			struct gen_matching match = match_to_gen(tau.p4(), gen_leps, gen_taus, gen_tau3ch, gen_w_prods, gen_b_prods);
			NT_tau_matching_gen   .push_back(match.closest);
			NT_tau_matching_gen_dR.push_back(match.dR);
			}

		if (tau.hasSecondaryVertex())
			{
			//const pat::tau::TauPFEssential::CovMatrix& flightCovMatr = taus[i].flightLengthCov();
			float x = tau.flightLength().x();
			float y = tau.flightLength().y();
			float z = tau.flightLength().z();
			NT_tau_flightLength.push_back(x*x + y*y + z*z);
			NT_tau_flightLengthSignificance.push_back(tau.flightLengthSig());
			/*
			LogInfo("Demo") << "flightLengthSig = "<< taus[i].flightLengthSig();
			LogInfo("Demo") << "flightLength    = "<< taus[i].flightLength().x() << ',' << taus[i].flightLength().y() << ',' << taus[i].flightLength().z();
			LogInfo("Demo") << "flightLength Covar = " << endl <<
				flightCovMatr(0,0) << ',' << flightCovMatr(0,1) << ',' << flightCovMatr(0,2) << endl << 
				flightCovMatr(1,0) << ',' << flightCovMatr(1,1) << ',' << flightCovMatr(2,2) << endl << 
				flightCovMatr(2,0) << ',' << flightCovMatr(2,1) << ',' << flightCovMatr(2,2) << endl;
			*/
			}
		else
			{
			NT_tau_flightLength.push_back(-111);
			NT_tau_flightLengthSignificance.push_back(-111);
			}

		// dR match to jet
		Int_t matched_jet_number = -1;
		for (unsigned int i=0; i<selJetsNoLep.size(); i++)
			{
			pat::Jet& jet = selJetsNoLep[i];
			if (reco::deltaR(jet, tau) < 0.4)
				{
				matched_jet_number = i;
				break;
				}
			}

		NT_tau_dR_matched_jet.push_back(matched_jet_number); // number of the jet in jet vectors, if no match = -1

		// Re-reconstructed Secondary Vertex (SV) for taus with 3 pions
		// for each tau save the bool if fit is ok
		// and parameters,
		// default value for the parameters of not good fit are 999.
		// Also save matching quality (= sum of dR of the matching to general tracks).
		bool fitOk_SV = false;  
		std::vector<double > tracksToBeRemoved; // compared by Pt due to the conflict of comparing const and not const iterators
		// if fit ok some tracks have to be removed from PV refit
		TransientVertex transVtx_SV;
		double matchingQuality(0);
		if (tau.decayMode() > 9)
			{
			/*
			 * 1) unpack things necessary for the fitter, tau candidate tracks, beamSpot & general tracks
			 * 2) match tracks it tau candidates
			 * 3) fit adaptivefitter to the matched tracks
			 * 4) save with "quality score" = sum of dR of the matching
			 */

			// PAT tau tracks
			reco::CandidatePtrVector sigCands = tau.signalChargedHadrCands();//signalCands();

			// rebuild tau tracks (supposed o be more precise)
			// match allTracks (not pvTracks) to tau tracks
			std::vector<reco::TransientTrack> transTracks_tau;  

			for (reco::CandidatePtrVector::const_iterator itr = sigCands.begin(); itr != sigCands.end(); ++itr)
				{
				double deR(999.); 
				double checkqual(0);
				reco::Track closestTrack;

				//for(auto iter: pvTracks)
				for(auto iter: allTracks)
					{
					if(std::find(tracksToBeRemoved.begin(), tracksToBeRemoved.end(), iter.pt())!=tracksToBeRemoved.end())
						continue;
					if( sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2))  < deR)
						{
						deR = sqrt(pow(iter.eta() - (*itr)->p4().eta(),2) + pow(iter.phi() - (*itr)->p4().phi(),2));
						checkqual=deR;
						closestTrack = iter;
						}
					}

				matchingQuality+=checkqual;
				tracksToBeRemoved.push_back(closestTrack.pt());
				transTracks_tau.push_back(transTrackBuilder->build(closestTrack));
				//cout<<"  closestTrackiter eta  :  "<<   closestTrack.eta() << "   phi   " << closestTrack.phi() << "    pt  "<< closestTrack.pt() <<endl;
				}

			// do the tau SV refit
			if (transTracks_tau.size() >= 2 )
				{
				AdaptiveVertexFitter avf;
				avf.setWeightThreshold(0.001); 
				try
					{
					transVtx_SV = avf.vertex(transTracks_tau, beamSpot);
					fitOk_SV = true; 
					}
				catch (...)
					{
					fitOk_SV = false; 
					std::cout<<"Vtx fit failed!"<<std::endl;
					}
				}

			fitOk_SV = fitOk_SV && transVtx_SV.isValid() && fabs(transVtx_SV.position().x())<1 && fabs(transVtx_SV.position().y())<1;
			}

		// save the vertex
		//NT_tau_SV_fit_isOk.push_back(fitOk_SV);
		if(fitOk_SV)
			{
			// store the index of refited taus to the all taus vectors
			NT_tau_refited_index.push_back(NT_tau_SV_fit_matchingQuality.size()); // the last number in current vectors of refited taus

			// now add the new tau to refitted taus
			// set the tracks corresponding to tau to be removed for PV refit
			// tracksToBeRemoved_PV
			for (unsigned int i=0; i<tracksToBeRemoved.size(); i++)
				tracksToBeRemoved_PV.push_back(tracksToBeRemoved[i]);

			///NOTE: we take original vertex z position, as this gives the best reults on CP
			///variables. To be understood; probable reason are missing tracks with Pt<0.95GeV
			NT_tau_SV_fit_matchingQuality.push_back(matchingQuality);
			NT_tau_SV_fit_x.push_back(transVtx_SV.position().x());
			NT_tau_SV_fit_y.push_back(transVtx_SV.position().y());
			NT_tau_SV_fit_z.push_back(transVtx_SV.position().z());
			reco::Vertex secondaryVertex = transVtx_SV;
			// covariance matrix
			math::Error<3>::type svCov;
			secondaryVertex.fill(svCov);
			NT_tau_SV_cov.push_back(svCov);

			// and save the tracks (for kinematic fits on them)
			//transTracks_tau // <- not sure about using these for momentum TODO: check
			// using Signal Candidates of tau
			reco::CandidatePtrVector sigCands = tau.signalChargedHadrCands();
			// for control number of signal candidates (should only be = 3, the 3pi decay):
			//NT_tau_SV_fit_ntracks.push_back(sigCands.size());
			// tests show all is ok

			bool i = true;
			for (reco::CandidatePtrVector::const_iterator itr_cand = sigCands.begin(); itr_cand != sigCands.end(); ++itr_cand)
				{
				// save as Same Sign track
				// (something is reverted in signs here)
				if ((*itr_cand)->charge() * tau.pdgId() > 0)
					{
					NT_tau_SV_fit_track_OS_p4.push_back((*itr_cand)->p4());
					}
				else if (i) // first OS track
					{
					NT_tau_SV_fit_track_SS1_p4.push_back((*itr_cand)->p4());
					i = false;
					}
				else // second OS track
					NT_tau_SV_fit_track_SS2_p4.push_back((*itr_cand)->p4());
				}

			// loop through tracks and save their impact parameters
			// and match quality dR
			double min_dR_os(99999.), min_dR_ss1(99999.), min_dR_ss2(99999.);
			int matched_track_OS = -1, matched_track_SS1 = -1, matched_track_SS2 = -1;
			for(size_t i=0; i<track_cands->size(); ++i)
				{
				// TODO: these requirements are probably the reasone some tracks are not found for tau sigCands
				if((*track_cands)[i].charge()==0 || (*track_cands)[i].vertexRef().isNull()) continue;
				if(!(*track_cands)[i].bestTrack()) continue;

				auto track = (*track_cands)[i].bestTrack();

				// find closest matches to general track
				double dR_os  = sqrt(pow(track->eta() - NT_tau_SV_fit_track_OS_p4 .back().eta(), 2) + pow(track->phi() - NT_tau_SV_fit_track_OS_p4 .back().phi(), 2));
				double dR_ss1 = sqrt(pow(track->eta() - NT_tau_SV_fit_track_SS1_p4.back().eta(), 2) + pow(track->phi() - NT_tau_SV_fit_track_SS1_p4.back().phi(), 2));
				double dR_ss2 = sqrt(pow(track->eta() - NT_tau_SV_fit_track_SS2_p4.back().eta(), 2) + pow(track->phi() - NT_tau_SV_fit_track_SS2_p4.back().phi(), 2));
				// there is no order in saving tracks
				// in principle two sigCands can match to the same track
				// but then dR should be large
				// and matchQ (matchQuality) will be large
				if (dR_os < min_dR_os)
					{
					min_dR_os = dR_os;
					matched_track_OS = i;
					}
				if (dR_ss1 < min_dR_ss1)
					{
					min_dR_ss1 = dR_ss1;
					matched_track_SS1 = i;
					}
				if (dR_ss2 < min_dR_ss2)
					{
					min_dR_ss2 = dR_ss2;
					matched_track_SS2 = i;
					}
				}

			// tracks are matched, save parameters

			// quality of match
			NT_tau_SV_fit_track_OS_matched_track_dR .push_back(min_dR_os);
			NT_tau_SV_fit_track_SS1_matched_track_dR.push_back(min_dR_ss1);
			NT_tau_SV_fit_track_SS2_matched_track_dR.push_back(min_dR_ss2);
			NT_tau_SV_fit_track_matchQ.push_back(min_dR_os + min_dR_ss1 + min_dR_ss2);

			// track parameters:
			//  - vector of impact parameter
			//  - vertex reference key (0 is the PV of the event)
			//  - vertex association quality, whether the track was used in the fit

			// OS
			int track_index;
			unsigned int key;
			int quality;
			//TVector3 impact;

			track_index = matched_track_OS;

			key = (*track_cands)[track_index].vertexRef().key();
			quality = (*track_cands)[track_index].pvAssociationQuality();

			auto ref_vertex = *((*track_cands)[track_index].vertexRef());
			auto closest_point = (*track_cands)[track_index].vertex();
			auto distance = closest_point - ref_vertex.position();
			// distance is of class ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>
			//impact.SetXYZ(distance.x(), distance.y(), distance.z());

			NT_tau_SV_fit_track_OS_matched_track_vtxkey.push_back(key);
			NT_tau_SV_fit_track_OS_matched_track_vtxQ.push_back(quality);
			NT_tau_SV_fit_track_OS_matched_track_b.push_back(distance);

			// SS1
			track_index = matched_track_SS1;

			key = (*track_cands)[track_index].vertexRef().key();
			quality = (*track_cands)[track_index].pvAssociationQuality();

			ref_vertex = *((*track_cands)[track_index].vertexRef());
			closest_point = (*track_cands)[track_index].vertex();
			distance = closest_point - ref_vertex.position();
			// distance is of class ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>
			//impact.SetXYZ(distance.x(), distance.y(), distance.z());

			NT_tau_SV_fit_track_SS1_matched_track_vtxkey.push_back(key);
			NT_tau_SV_fit_track_SS1_matched_track_vtxQ.push_back(quality);
			NT_tau_SV_fit_track_SS1_matched_track_b.push_back(distance);

			// SS2
			track_index = matched_track_SS2;

			key = (*track_cands)[track_index].vertexRef().key();
			quality = (*track_cands)[track_index].pvAssociationQuality();

			ref_vertex = *((*track_cands)[track_index].vertexRef());
			closest_point = (*track_cands)[track_index].vertex();
			distance = closest_point - ref_vertex.position();
			// distance is of class ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>
			//impact.SetXYZ(distance.x(), distance.y(), distance.z());

			NT_tau_SV_fit_track_SS2_matched_track_vtxkey.push_back(key);
			NT_tau_SV_fit_track_SS2_matched_track_vtxQ.push_back(quality);
			NT_tau_SV_fit_track_SS2_matched_track_b.push_back(distance);

			//TVector3 tr_ss2;
			//TVector3 tr_os ;
			//TVector3 tr_ss1;
			//tr_ss2.SetXYZ(NT_tau_SV_fit_track_SS2_p4.back().X(), NT_tau_SV_fit_track_SS2_p4.back().Y(), NT_tau_SV_fit_track_SS2_p4.back().Z());
			//tr_os .SetXYZ(NT_tau_SV_fit_track_OS_p4 .back().X(), NT_tau_SV_fit_track_OS_p4 .back().Y(), NT_tau_SV_fit_track_OS_p4 .back().Z());
			//tr_ss1.SetXYZ(NT_tau_SV_fit_track_SS1_p4.back().X(), NT_tau_SV_fit_track_SS1_p4.back().Y(), NT_tau_SV_fit_track_SS1_p4.back().Z());

			// let's save also Vector3 of tracks? to not have to convert everything every time
			NT_tau_SV_fit_track_OS_matched_track_p3 .push_back(NT_tau_SV_fit_track_OS_p4 .back().Vect());
			NT_tau_SV_fit_track_SS1_matched_track_p3.push_back(NT_tau_SV_fit_track_SS1_p4.back().Vect());
			NT_tau_SV_fit_track_SS2_matched_track_p3.push_back(NT_tau_SV_fit_track_SS2_p4.back().Vect());

			// geometrical SV
			struct sv_pair geom_SV = geometrical_SV(
				NT_tau_SV_fit_track_SS2_matched_track_b.back(), NT_tau_SV_fit_track_SS2_matched_track_p3.back(),
				NT_tau_SV_fit_track_OS_matched_track_b .back(), NT_tau_SV_fit_track_OS_matched_track_p3 .back(),
				NT_tau_SV_fit_track_SS1_matched_track_b.back(), NT_tau_SV_fit_track_SS1_matched_track_p3.back());
			NT_tau_SV_geom_flightLen.push_back(geom_SV.flightLength);
			NT_tau_SV_geom_flightLenSign.push_back(geom_SV.flightLengthSignificance);
			}
		else
			{
			// store default index of refited taus to the all taus vectors
			NT_tau_refited_index.push_back(-1);
			}
		}
	LogInfo ("Demo") << "saved taus";

	// if some tracks were removed (SV was refitted for a tau)
	// refit PV
	// build transient tracks of the rest of the tracks
	// for PV refitting
	if (tracksToBeRemoved_PV.size()>0)
		{
		TransientVertex transVtx_PV;
		std::vector<reco::TransientTrack> transTracks_nottau;

		for(auto iter: pvTracks)
			{
			if(std::find(tracksToBeRemoved_PV.begin(), tracksToBeRemoved_PV.end(), iter.pt())!=tracksToBeRemoved_PV.end())
				continue;
			transTracks_nottau.push_back(transTrackBuilder->build(iter));
			}

		bool fitOk_PV = false;  
		if (transTracks_nottau.size() >= 2 ) {
			AdaptiveVertexFitter avf;
			avf.setWeightThreshold(0.001); 
			try {
				transVtx_PV = avf.vertex(transTracks_nottau, beamSpot);
				fitOk_PV = true; 
				}
			catch (...) {
				fitOk_PV = false; 
				std::cout<<"Vtx fit failed!"<<std::endl;
				}
			}

		fitOk_PV = fitOk_PV && transVtx_PV.isValid() && fabs(transVtx_PV.position().x())<1 && fabs(transVtx_PV.position().y())<1;

		// save the vertex
		NT_PV_fit_isOk = fitOk_PV;
		if(fitOk_PV)
			{
			// save position and cov of vertex
			NT_PV_fit_x = transVtx_PV.position().x();
			NT_PV_fit_y = transVtx_PV.position().y();
			NT_PV_fit_z = transVtx_PV.position().z();
			reco::Vertex pvVertex = transVtx_PV;
			// covariance matrix
			//math::Error<3>::type pvCov;
			//pvVertex.fill(pvCov);
			//NT_PV_cov = pvCov;
			pvVertex.fill(NT_PV_cov);
			}
		}



	LogInfo ("Demo") << "all particles/objects are selected, nbjets = " << NT_nbjets;

	//bool record_ntuple = (isSingleMu || isSingleE || pass_dileptons) && NT_nbjets > 0 && NT_tau_IDlev.size() > 0; // at least 1 b jet and 1 loose tau
	bool pass_leptons = clean_lep_conditions && leps_passed_relIso && selLeptons.size() == 1; // tau_ID oriented scheme, Dileptons are separate
	bool pass_leptons_all_iso = clean_lep_conditions && selLeptons.size() > 0 && selLeptons.size() < 3;
	bool record_ntuple = false;

	if (record_tauID)
		{
		/*
		 * tau ID preselection
		 * about twice less events then in our preselection with b jets
		 *
		 * should contain good WJets control sample
		 */
		record_ntuple |= pass_leptons && NT_ntaus > 0;
		}
	if (record_tauIDantiIso)
		{
		record_ntuple |= pass_leptons_all_iso && NT_ntaus > 0;
		}
	if (record_bPreselection)
		{
		/* leptons and at least 1 b jet
		 * our old preselection
		 */
		record_ntuple |= pass_leptons && NT_nbjets > 0;
		}
	if (record_MonitorHLT)
		{
		/* the HLT efficiency study
		 * only few events in lepMonitorTrigger
		 */
		record_ntuple |= pass_leptons && lepMonitorTrigger;
		}
	if (record_ElMu)
		{
		// all el-mu events (it's mainly TTbar, so should be few)
		record_ntuple |= clean_lep_conditions && abs(NT_leps_ID) == 143;
		}
	if (record_Dilep)
		{
		/* just all dileptons -- lots of DY here
		 * control for lepton IDs
		 * about = to tau ID preselection
		 */
		record_ntuple |= clean_lep_conditions && selLeptons.size() == 2;
		}
	if (record_jets)
		{
		record_ntuple |= jetsHLT && selJetsNoLep.size() > 0; // these are all-eta jets with pt-cut, anti-lep dR and Loose ID..
		}

	if (record_ntuple)
		{
		LogInfo ("Demo") << "recording " << record_tauID << record_bPreselection << record_MonitorHLT << record_ElMu << record_Dilep;
		event_counter ->Fill(event_checkpoint++);
		weight_counter->Fill(event_checkpoint, weight);

		NT_output_ttree->Fill();
		// EDM output
		//   but it's under if! not every event will store stuff -- see if it plays well with rest of EDM system
		//iEvent.put(NT_lep_p4, "leps_p4");
		// some use unknown C++ feature std::move :
		//iEvent.put(std::move(obj), "objname");
		}



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
NtuplerAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
NtuplerAnalyzer::endJob() 
{
if (isLocal && !isMC){
	goodLumiFilter.FindLumiInFiles(urls); // urls! why are they here at all? no even 1 comment in that "Utilities!"
	goodLumiFilter.DumpToJson(((outUrl.ReplaceAll(".root",""))+".json").Data());
	}
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NtuplerAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NtuplerAnalyzer);
