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
#include "UserCode/NtuplerAnalyzer/interface/PatUtils.h"
#include "UserCode/NtuplerAnalyzer/interface/MacroUtils.h"
// lumiUtils
#include "UserCode/NtuplerAnalyzer/interface/LumiUtils.h"
// couple functions processing leptons
#include "UserCode/NtuplerAnalyzer/interface/ProcessingMuons.h"
#include "UserCode/NtuplerAnalyzer/interface/ProcessingElectrons.h"
#include "UserCode/NtuplerAnalyzer/interface/ProcessingTaus.h"
#include "UserCode/NtuplerAnalyzer/interface/ProcessingJets.h"
#include "UserCode/NtuplerAnalyzer/interface/ProcessingBJets.h"
#include "UserCode/NtuplerAnalyzer/interface/ProcessingHLT.h"
#include "UserCode/NtuplerAnalyzer/interface/ProcessingGenParticles.h"
#include "UserCode/NtuplerAnalyzer/interface/ProcessingDRCleaning.h"

#include "UserCode/NtuplerAnalyzer/interface/GenDistrs.h"

//#include "UserCode/NtuplerAnalyzer/interface/RoccoR.cc"
//#include "UserCode/NtuplerAnalyzer/interface/RoccoR.h"

//#include "UserCode/NtuplerAnalyzer/interface/handy.h"


/*
 * hadnling the per-year sources (plugins, functions etc)
 * Ntupler receives their output in some per-year-specific format and packs into the TTree as needed
 */

// defaults
//#ifndef BFRAG_PROC
//#define BFRAG_PROC 2016
//#endif
// no defaults is probably safer



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

#define MAX_NVTX 80

// 80 numbers
static double pu_vector_mu[MAX_NVTX] = {0, 0.360609416811339, 0.910848525427002, 1.20629960507795, 0.965997726573782, 1.10708082813183, 1.14843491548622, 0.786526251164482, 0.490577792661333, 0.740680941110478,
0.884048630953726, 0.964813189764159, 1.07045369167689, 1.12497267309738, 1.17367530613108, 1.20239808206413, 1.20815108390021, 1.20049333094509, 1.18284686347315, 1.14408796655615,
1.0962284704313, 1.06549162803223, 1.05151011089581, 1.05159666626121, 1.05064452078328, 1.0491726301522, 1.05772537082991, 1.07279673875566, 1.0837536468865, 1.09536667397119,
1.10934472980173, 1.09375894592864, 1.08263679568271, 1.04289345879947, 0.9851490341672, 0.909983816540809, 0.821346330143864, 0.71704523475871, 0.609800913869359, 0.502935245638477,
0.405579825620816, 0.309696044611377, 0.228191137503131, 0.163380359253309, 0.113368437957202, 0.0772279997453792, 0.0508111733313502, 0.0319007262683943, 0.0200879459309245, 0.0122753366005436,
0.00739933885813127, 0.00437426967257811, 0.00260473545284139, 0.00157047254226743, 0.000969500595715493, 0.000733193118123283, 0.000669817107713128, 0.000728548958604492, 0.000934559691182011, 0.00133719688378802,
0.00186652283903214, 0.00314422244976771, 0.00406954793369611, 0.00467888840511915, 0.00505224284441512, 0.00562827194936864, 0.0055889504870752, 0.00522867039470319, 0.00450752163476433, 0.00395300774604375,
0.00330577167682956, 0.00308353042577215, 0.00277846504893301, 0.00223943190687725, 0.00196650068765464, 0.00184742734258922,  0, 0, 0, 0};

static double pu_vector_mu_up[MAX_NVTX] = {0, 0.351377216124927, 0.717199649125846, 1.14121536968772, 0.84885826611733, 1.00700929402897, 1.03428595270903, 0.717444379696992, 0.344078389355127, 0.499570875027422,
0.606614916257104, 0.632584599390169, 0.731450949466174, 0.827511723989754, 0.910682115553867, 0.960170981598162, 0.988896170761361, 1.02468865580207, 1.05296667126403, 1.05112033565679,
1.0269129153969, 1.00548641752714, 0.998316130432865, 1.01492587998551, 1.03753749807849, 1.05742218946485, 1.08503978097083, 1.12134132247053, 1.15585339474274, 1.19214399856171,
1.23308400947467, 1.24528804633732, 1.26786364716917, 1.26101551498967, 1.23297806722714, 1.18042533075471, 1.10534683838101, 1.00275591661645, 0.889094305531985, 0.768791254270252,
0.655054015673457, 0.533361034358457, 0.423095146361996, 0.329177839117034, 0.250352385505809, 0.188377378855567, 0.137852651411779, 0.0968577167707531, 0.0686240187247059, 0.0473889635126706,
0.0323695027438475, 0.0216752397914914, 0.0145119352923332, 0.00961177893634792, 0.00615582219138384, 0.00430085627914427, 0.00305735512896403, 0.00223567790438986, 0.00189369737638594, 0.00199585978316291,
0.00236236592656064, 0.00372472999463276, 0.00474687312579969, 0.00549508151576102, 0.00603023110946686, 0.0068545111910253, 0.00695838760530896, 0.00666224781277046, 0.00588243140681038, 0.00528714370892014,
0.00453424615273565, 0.00433985030329723, 0.00401493171035719, 0.00332436608713241, 0.00300063798808221, 0.00289925128977536, 0, 0, 0, 0};

static double pu_vector_mu_down[MAX_NVTX] = {0, 0.37361294640242, 1.1627791004568, 1.26890787896295, 1.10266790442705, 1.23456697093644, 1.26278991594152, 0.909648777562084, 0.759569490571151, 1.09035651921682,
1.34530547603283, 1.48713160105, 1.52535976889483, 1.49730550773404, 1.49792998045778, 1.49767851097519, 1.44431045398336, 1.3681909492045, 1.29912252494785, 1.2274279217797,
1.16525969099909, 1.12531044676724, 1.09094501417685, 1.06405434433422, 1.03997120824565, 1.0185716022098, 1.00560949501652, 0.997570939806059, 0.985543761409897, 0.972557804582185,
0.957832827239337, 0.9139572640153, 0.872252387173971, 0.808388185417578, 0.733817960498049, 0.650440963845892, 0.561688505024782, 0.466564380334112, 0.374428618658619, 0.28845274688129,
0.214909665968644, 0.149991974352384, 0.100014138338029, 0.0642260884603397, 0.0396553405911344, 0.0238687936736627, 0.0137921542898078, 0.00756854010632403, 0.00415483516246187, 0.00221776872027937,
0.00118249725637452, 0.000641889697310868, 0.000383647166012176, 0.000273637590071334, 0.000242902582071058, 0.000291239677209452, 0.000394091114279828, 0.000542541231466254, 0.000771067920964491, 0.00113596447675107,
0.00158061353194779, 0.00261959852500539, 0.00331800452823827, 0.00372426930370732, 0.00392086545082614, 0.00425479965493548, 0.00411256966391362, 0.00374240422174387, 0.00313603438166934, 0.00267155793176928,
0.00216878198028599, 0.00196249821290853, 0.00171433839159669, 0.00133866519755926, 0.00113810604240254, 0.00103447940224886, 0, 0, 0, 0};

static double pu_vector_el[MAX_NVTX] = {0,
   0.413231   ,    1.01701    ,    1.19502     ,   0.883906   ,    1.05852    ,    1.11823    ,    0.789439   ,    0.515477   ,    0.81338    ,    0.990148   ,
   1.0919     ,    1.21784    ,    1.28268     ,   1.33936    ,    1.37267    ,    1.38001    ,    1.37224    ,    1.35253    ,    1.30805    ,    1.25303    ,
   1.21761    ,    1.20085    ,    1.1987      ,   1.19257    ,    1.1807     ,    1.17079    ,    1.15238    ,    1.10667    ,    1.03375    ,    0.935086   ,
   0.793376   ,    0.65125    ,    0.502727    ,   0.369298   ,    0.25859    ,    0.173207   ,    0.110361   ,    0.0677957  ,    0.0403186  ,    0.0236369  ,
   0.0133546  ,    0.00746494 ,    0.00417626  ,   0.00233773 ,    0.0013288  ,    0.000757718,    0.000432788,    0.000266239,    0.000177605,    0.000137241,
   0.000125696,    0.000137018,    0.000167806 ,   0.000215108,    0.000313214,    0.000464376,    0.000669855,    0.000981399,    0.00148275 ,    0.00211313 ,
   0.0035872  ,    0.00465614 ,    0.005359    ,   0.00578897 ,    0.00645001 ,    0.00640537 ,    0.00599263 ,    0.00516618 ,    0.00453067 ,    0.00378886 ,
   0.00353415 ,    0.00318451 ,    0.0025667   ,   0.00225388 ,    0.00211741 ,    0          ,    0          ,    0          ,    0};

static double pu_vector_el_up[MAX_NVTX] = {0,
   0.402668    ,   0.803377   ,    1.15963    ,    0.764147   ,    0.966328    ,   0.995159   ,    0.71563    ,    0.354304   ,   0.541943   ,    0.674778   ,
   0.713035    ,   0.830366   ,    0.942616   ,    1.03882    ,    1.09589     ,   1.12909    ,    1.17068    ,    1.20376    ,   1.20191    ,    1.17404    ,
   1.14928     ,   1.14083    ,    1.15906    ,    1.18279    ,    1.20082     ,   1.22277    ,    1.24559    ,    1.25129    ,   1.23607    ,    1.19534    ,
   1.09539     ,   0.978694   ,    0.82546    ,    0.662451   ,    0.50547     ,   0.367764   ,    0.25369    ,    0.168007   ,   0.10706    ,    0.0667404  ,
   0.0397874   ,   0.0233291  ,    0.0136533  ,    0.00799737 ,    0.00476279  ,   0.00284044 ,    0.00167744 ,    0.00103389 ,   0.000648432,    0.000427764,
   0.000303899 ,   0.000247672,    0.000236803,    0.000258026,    0.000345092 ,   0.000494341,    0.000708128,    0.00104444 ,   0.00159927 ,    0.00231779 ,
   0.00400894  ,   0.00530831 ,    0.00623822 ,    0.00688571 ,    0.00784455  ,   0.00797042 ,    0.00763388 ,    0.00674129 ,   0.00605947 ,    0.00519674 ,
   0.00497399  ,   0.00460162 ,    0.00381017 ,    0.00343914 ,    0.00332295  ,   0          ,    0          ,    0          ,   0};

static double pu_vector_el_down[MAX_NVTX] = {0,
   0.428107   ,   1.29388    ,    1.22078   ,    1.02596    ,    1.1769     ,    1.24377    ,    0.921862   ,    0.814769   ,    1.20901   ,    1.51527   ,
   1.68838    ,   1.73792    ,    1.70826   ,    1.70984    ,    1.71038    ,    1.65067    ,    1.56442    ,    1.48535    ,    1.40303   ,    1.33164   ,
   1.28514    ,   1.24342    ,    1.20714   ,    1.16839    ,    1.12262    ,    1.06993    ,    0.999693   ,    0.900043   ,    0.778486  ,    0.644942  ,
   0.497564   ,   0.37052    ,    0.259917  ,    0.174109   ,    0.111585   ,    0.0687061  ,    0.0404941  ,    0.0232033  ,    0.0129877 ,    0.00721863,
   0.00388086 ,   0.00206418 ,    0.00109735,    0.000585006,    0.000321242,    0.000184087,    0.000114353,    8.57172e-5 ,    7.68135e-5,    8.09537e-5,
   9.42381e-5 ,   0.000118387,    0.0001548 ,    0.000202628,    0.000294603,    0.000431519,    0.00061181 ,    0.000878668,    0.00129927,    0.0018102 ,
   0.00300153 ,   0.00380243 ,    0.0042683 ,    0.00449375 ,    0.00487654 ,    0.00471355 ,    0.0042893  ,    0.00359432 ,    0.00306198,    0.00248572,
   0.0022493  ,   0.00196487 ,    0.0015343 ,    0.00130443 ,    0.00118566 ,    0          ,    0          ,    0          ,    0};


//bool sort_CandidatesByPt(const pat::GenericParticle &a, const pat::GenericParticle &b)  { return a.pt()>b.pt(); }
bool sort_TausByIDByPt(const pat::Tau &a, const pat::Tau &b)
	{
	return (a.userInt("IDlev") != b.userInt("IDlev") ?
		a.userInt("IDlev") >  b.userInt("IDlev") :
		a.pt()>b.pt());
	}


// try simpler first
struct dR_matching {
	//Int_t sum;
	Bool_t  matched;
	Float_t dR;
};

struct dR_matching dR_match_to_HLTs(
        pat::Muon& lepton,
        vector<pat::TriggerObjectStandAlone>& trig_objs,    // input: trigger objects to match against (so, these should match HLT of interest)
        float dR_cut
        )
{
Float_t min_dR = 999.;
for (size_t i = 0; i < trig_objs.size(); i++)
        {
        pat::TriggerObjectStandAlone& obj = trig_objs[i];
	Float_t dR = reco::deltaR(lepton, obj);
        if (dR < min_dR)
                {
                min_dR = dR;
                }
        }
struct dR_matching ret = {.matched = min_dR < dR_cut, .dR = min_dR};
return ret;
}

struct dR_matching dR_match_to_HLTs(
        pat::Electron& lepton,
        vector<pat::TriggerObjectStandAlone>& trig_objs,    // input: trigger objects to match against (so, these should match HLT of interest)
        float dR_cut
        )
{
Float_t min_dR = 999.;
for (size_t i = 0; i < trig_objs.size(); i++)
        {
        pat::TriggerObjectStandAlone& obj = trig_objs[i];
	Float_t dR = reco::deltaR(lepton, obj);
        if (dR < min_dR)
                {
                min_dR = dR;
                }
        }
struct dR_matching ret = {.matched = min_dR < dR_cut, .dR = min_dR};
return ret;
}




const reco::Candidate* has_tau_mother(
	const reco::Candidate* part, unsigned int depth)

{
unsigned int pdgId = abs(part->pdgId());
if (pdgId == 15)
	return part;
else if (depth > 0)
	{
	// loop through mother
	const reco::Candidate* cand = NULL;
	for (unsigned int d_i=0; d_i < part->numberOfMothers(); d_i++)
		{
		const reco::Candidate * mother = part->mother(d_i);
		depth -= 1;
		cand = has_tau_mother(mother, depth);
		if (cand != NULL)
			return cand;
		}
	return cand;
	}
return NULL;
}

bool has_mother(const reco::Candidate& part, unsigned int target_pdgId)
{
unsigned int pdgId = abs(part.pdgId());
unsigned int n_mothers = part.numberOfMothers();
edm::LogInfo ("Demo") << "has_mother: " << target_pdgId << ' ' << pdgId << ' ' << n_mothers;

if (pdgId == target_pdgId)
	return true;
else if (n_mothers == 0)
	return false;
else
	{
	// check pdgIDs of all mothers
	for (unsigned int m_i=0; m_i<n_mothers; m_i++)
		{
		unsigned int m_pdg = abs(part.mother(m_i)->pdgId());
		edm::LogInfo ("Demo") << "has_mother: m " << m_pdg;
		if (m_pdg == target_pdgId) return true;
		}
	// if still no match -- follow the first mother
	return has_mother(*(part.mother(0)), target_pdgId);
	}

edm::LogInfo ("Demo") << "has_mother: impossible case, returning the default false";
return false;
}






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
	vector<LorentzVector>& gen_taulep,
	vector<LorentzVector>& gen_w_prods,
	vector<LorentzVector>& gen_b_prods,
	vector<int>& gid_leps,
	vector<int>& gid_taus,
	vector<int>& gid_tau3ch,
	vector<int>& gid_taulep,
	vector<int>& gid_w_prods,
	vector<int>& gid_b_prods)

{
Float_t min_dR = 9999.;
Float_t dR_cut = 0.4;
Int_t   min_id = 0;
//LorentzVector min_p4(0., 0., 0., 0.);

edm::LogInfo ("Demo") << "act gen sizes " << gen_leps.size() << gen_taus.size() << gen_tau3ch.size() << gen_w_prods.size() << gen_b_prods.size();
edm::LogInfo ("Demo") << "act obj p4    " << p4.pt() << p4.eta() << p4.phi();
edm::LogInfo ("Demo") << "act gen pts  " << (gen_leps.size()>0 ? gen_leps[0].pt() : 0)
	<< (gen_taus.size()>0? gen_taus[0].pt():0)
	<< (gen_tau3ch.size()>0? gen_tau3ch[0].pt():0)
	<< (gen_w_prods.size()>0? gen_w_prods[0].pt():0)
	<< (gen_b_prods.size()>0? gen_b_prods[0].pt():0);
edm::LogInfo ("Demo") << "act matches";

// I need in-place quick loop with list of all inputs here, how to do it in C? -- do it later

// so 0 is no match to anything
// numbers for different kinds of particles
// and the sign is for the pdgId sign of the provenance particle of these final states
// i.e. the sign of W, b, tau or leptons
int gen_id = 1;
for (unsigned int i = 0; i<gen_leps.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_leps[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		int id = gid_leps[i];
		min_id = (id > 0 ? gen_id : -gen_id);
		//min_p4 = gen_leps[i];
		}
	}

gen_id++;
for (unsigned int i = 0; i<gen_taus.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_taus[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		int id = gid_taus[i];
		min_id = (id > 0 ? gen_id : -gen_id);
		//min_p4 = gen_taus[i];
		}
	}

gen_id++;
for (unsigned int i = 0; i<gen_tau3ch.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_tau3ch[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		int id = gid_tau3ch[i];
		min_id = (id > 0 ? gen_id : -gen_id);
		//min_p4 = gen_tau3ch[i];
		}
	}

gen_id++;
for (unsigned int i = 0; i<gen_taulep.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_taulep[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		int id = gid_taulep[i];
		min_id = (id > 0 ? gen_id : -gen_id);
		//min_p4 = gen_taulep[i];
		}
	}

gen_id++;
for (unsigned int i = 0; i<gen_w_prods.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_w_prods[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		int id = gid_w_prods[i];
		min_id = (id > 0 ? gen_id : -gen_id);
		//min_p4 = gen_w_prods[i];
		}
	}

gen_id++;
for (unsigned int i = 0; i<gen_b_prods.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_b_prods[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		int id = gid_b_prods[i];
		min_id = (id > 0 ? gen_id : -gen_id);
		//min_p4 = gen_b_prods[i];
		}
	}

//struct gen_matching match = {.closest=min_id, .dR=min_dR, .p4=min_p4};
struct gen_matching match = {.closest=min_id, .dR=min_dR};
return match;
}

// genmatchig to flat gen collections
struct gen_matching_collection {
	//Int_t sum;
	Int_t index;
	Float_t dR;
};

struct gen_matching_collection match_to_gen_collection(const LorentzVector& p4,
	vector<LorentzVector>& gen_collection)

{
Float_t min_dR = 9999.;
Float_t dR_cut = 0.4;
Int_t   min_index = -99;

// I need in-place quick loop with list of all inputs here, how to do it in C? -- do it later

// so 0 is no match to anything
// numbers for different kinds of particles
// and the sign is for the pdgId sign of the provenance particle of these final states
// i.e. the sign of W, b, tau or leptons
for (unsigned int i = 0; i<gen_collection.size(); i++)
	{
	Float_t dR = reco::deltaR(p4, gen_collection[i]);
	if (dR < dR_cut && dR < min_dR)
		{
		min_dR = dR;
		min_index = i;
		}
	}

struct gen_matching_collection match = {.index=min_index, .dR=min_dR};
return match;
}


enum nWeights {
	NOM_EVENTS = 1,
	NOMINAL    = 2,
	MU_PU      = 3,
	MU_PUUp    = 4,
	MU_PUDown  = 5,
	El_PU      = 6,
	El_PUUp    = 7,
	El_PUDown  = 8,
	TOPPT      = 9,
	// renorm refact scales
	M_NOM    = 10,
	MrUp     = 11,
	MrDown   = 12,
	MfUp     = 13,
	MfDown   = 14,
	MfrUp    = 15,
	MfrDown  = 16,
	// fragmentatopn
	Central       = 21,
	FragUp        = 22,
	FragDown      = 23,
	SemilepBRUp   = 24,
	SemilepBRDown = 25,
	PetersonUp    = 26,
	// PDF and alphaS
	PDF_NOM       = 27,
	AlphaSUp      = 28,
	AlphaSDown    = 29,

	// PDFs
	PDFCT14n1     =  31,
	PDFCT14n2     =  32,
	PDFCT14n3     =  33,
	PDFCT14n4     =  34,
	PDFCT14n5     =  35,
	PDFCT14n6     =  36,
	PDFCT14n7     =  37,
	PDFCT14n8     =  38,
	PDFCT14n9     =  39,
	PDFCT14n10    =  40,
	PDFCT14n11    =  41,
	PDFCT14n12    =  42,
	PDFCT14n13    =  43,
	PDFCT14n14    =  44,
	PDFCT14n15    =  45,
	PDFCT14n16    =  46,
	PDFCT14n17    =  47,
	PDFCT14n18    =  48,
	PDFCT14n19    =  49,
	PDFCT14n20    =  50,
	PDFCT14n21    =  51,
	PDFCT14n22    =  52,
	PDFCT14n23    =  53,
	PDFCT14n24    =  54,
	PDFCT14n25    =  55,
	PDFCT14n26    =  56,
	PDFCT14n27    =  57,
	PDFCT14n28    =  58,
	PDFCT14n29    =  59,
	PDFCT14n30    =  60,
	PDFCT14n31    =  61,
	PDFCT14n32    =  62,
	PDFCT14n33    =  63,
	PDFCT14n34    =  64,
	PDFCT14n35    =  65,
	PDFCT14n36    =  66,
	PDFCT14n37    =  67,
	PDFCT14n38    =  68,
	PDFCT14n39    =  69,
	PDFCT14n40    =  70,
	PDFCT14n41    =  71,
	PDFCT14n42    =  72,
	PDFCT14n43    =  73,
	PDFCT14n44    =  74,
	PDFCT14n45    =  75,
	PDFCT14n46    =  76,
	PDFCT14n47    =  77,
	PDFCT14n48    =  78,
	PDFCT14n49    =  79,
	PDFCT14n50    =  80,
	PDFCT14n51    =  81,
	PDFCT14n52    =  82,
	PDFCT14n53    =  83,
	PDFCT14n54    =  84,
	PDFCT14n55    =  85,
	PDFCT14n56    =  86,
};

// substituting stuff with another enum
enum names_Mfr {
	MUf_nom_MUr_nom    = 0,
	MUf_up_MUr_nom     = 1,
	MUf_down_MUr_nom   = 2,
	MUf_nom_MUr_up     = 3,
	MUf_up_MUr_up      = 4,
	MUf_down_MUr_up    = 5,
	MUf_nom_MUr_down   = 6,
	MUf_up_MUr_down    = 7,
	MUf_down_MUr_down  = 8
};


TRandom3 *r3 = new TRandom3();

//RoccoR  rc("rcdata.2016.v3");
//RoccoR  rc("UserCode/NtuplerAnalyzer/interface/rcdata.2016.v3");

typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> ROOT_TTree_vector3D;


struct sv_pair {
	double flightLength;
	double flightLengthSignificance;
};

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

	//string ("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017");
	string tau_VLoose_ID, tau_Loose_ID , tau_Medium_ID, tau_Tight_ID , tau_VTight_ID;
	//string ("byVVTightIsolationMVArun2017v2DBoldDMwLT2017");
	//
	edm::EDGetTokenT<reco::VertexCollection> vtx_;
	edm::EDGetTokenT<double> rho_;
	edm::EDGetTokenT<edm::TriggerResults> trigResults_, trigResultsRECO_, trigResultsPAT_;

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
	edm::EDGetTokenT<std::vector<reco::GenJet>  > genLepsToken_, genJetsToken_, genFatJetsToken_; // there is a twist why

	//edm::EDGetTokenT<edm::ValueMap<float> > petersonFragToken_;
	edm::EDGetTokenT<edm::ValueMap<float> > upFragToken_, centralFragToken_, downFragToken_, PetersonFragToken_, semilepbrUpToken_, semilepbrDownToken_;

	//RecoilCorrector* recoilPFMetCorrector;
	//TH2D* zPtMass_histo;

	bool record_ElTau, record_MuTau, record_tauCands, record_tauID, record_tauIDantiIso, record_bPreselection, record_MonitorHLT, record_ElMu, record_Dilep, record_jets, record_signal, record_all;
	bool record_lepTauVL_b, record_ElMu_b;

	TString dtag;
	bool isMC, aMCatNLO, isWJets, isDY, isTT, isSingleTop, is2016legacy, is2017data, is2017legacy;
	bool isLocal;
	bool withHLT;
	string  HLT_source,
		low_pt_muHLT_MC1  , low_pt_muHLT_MC2  ,
		low_pt_muHLT_Data1, low_pt_muHLT_Data2,
		low_pt_elHLT_Data , low_pt_elHLT_MC,
		muHLT_MC1  , muHLT_MC2  ,
		muHLT_Data1, muHLT_Data2,
		elHLT_Data , elHLT_MC,
		lepMonitorHLT;

	string
		low_pt32_elHLT, low_pt28_150HT_elHLT, low_pt30_35PFJet_elHLT,
		elmuHLT_1, elmuHLT_2, elmuHLT_3, elmuHLT_4,
		elelHLT_1, elelHLT_2,
		eltauHLT, mutauHLT1, mutauHLT2;


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
	double btag_threshold, btag_threshold_Medium;

	lumiUtils::GoodLumiFilter goodLumiFilter;
	std::vector < std::string > urls; // = runProcess.getUntrackedParameter < std::vector < std::string > >("input");
	TString outUrl;

	//triggerObjectsHandle.getByLabel(ev, "selectedPatTrigger");
	//iEvent.getByToken(triggerObjects_, triggerObjectsHandle);
	edm::EDGetTokenT<vector<pat::TriggerObjectStandAlone>> triggerObjects_;
	edm::InputTag triggerObjects_InputTag;
	//edm::EDGetTokenT triggerObjects_; // they do this in the workbook and it does not compile
	// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Trigger
	// error: invalid use of template-name 'edm::EDGetTokenT' without an argument list

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
	TH1D *event_counter, *weight_counter, *systematic_weights; 

	/*
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > NT_lep_p4;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* pt_lep_p4; // yep, vectors of complex objects require additional persistent pointers
	*/
	#define NTUPLE_INTERFACE_CLASS_DECLARE
	#include "UserCode/NtuplerAnalyzer/interface/ntupleOutput.h"
	#undef NTUPLE_INTERFACE_CLASS_DECLARE
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
record_ElTau         (iConfig.getParameter<bool>("record_ElTau"))         ,
record_MuTau         (iConfig.getParameter<bool>("record_MuTau"))         ,
record_tauCands      (iConfig.getParameter<bool>("record_tauCands"))      ,
record_tauID         (iConfig.getParameter<bool>("record_tauID"))         ,
record_tauIDantiIso  (iConfig.getParameter<bool>("record_tauIDantiIso"))  ,
record_bPreselection (iConfig.getParameter<bool>("record_bPreselection")) ,
record_MonitorHLT    (iConfig.getParameter<bool>("record_MonitorHLT"))    ,
record_ElMu          (iConfig.getParameter<bool>("record_ElMu"))          ,
record_Dilep         (iConfig.getParameter<bool>("record_Dilep"))         ,
record_jets          (iConfig.getParameter<bool>("record_jets"))          ,
record_signal        (iConfig.getParameter<bool>("record_signal"))        ,
record_all           (iConfig.getParameter<bool>("record_all"))           ,
record_lepTauVL_b    (iConfig.getParameter<bool>("record_lepTauVL_b"))    ,
record_ElMu_b        (iConfig.getParameter<bool>("record_ElMu_b"))        ,
dtag       (iConfig.getParameter<std::string>("dtag")),
isMC       (iConfig.getParameter<bool>("isMC")),
is2016legacy       (iConfig.getParameter<bool>("is2016legacy")),
isLocal    (iConfig.getParameter<bool>("isLocal")),
withHLT    (iConfig.getParameter<bool>("withHLT")),
HLT_source (iConfig.getParameter<std::string>("HLT_source")),
low_pt_muHLT_MC1  (iConfig.getParameter<std::string>("low_pt_muHLT_MC1")),
low_pt_muHLT_MC2  (iConfig.getParameter<std::string>("low_pt_muHLT_MC2")),
low_pt_muHLT_Data1(iConfig.getParameter<std::string>("low_pt_muHLT_Data1")),
low_pt_muHLT_Data2(iConfig.getParameter<std::string>("low_pt_muHLT_Data2")),
low_pt_elHLT_Data (iConfig.getParameter<std::string>("low_pt_elHLT_Data")),
low_pt_elHLT_MC   (iConfig.getParameter<std::string>("low_pt_elHLT_MC")),
muHLT_MC1  (iConfig.getParameter<std::string>("muHLT_MC1")),
muHLT_MC2  (iConfig.getParameter<std::string>("muHLT_MC2")),
muHLT_Data1(iConfig.getParameter<std::string>("muHLT_Data1")),
muHLT_Data2(iConfig.getParameter<std::string>("muHLT_Data2")),
elHLT_Data (iConfig.getParameter<std::string>("elHLT_Data")),
elHLT_MC   (iConfig.getParameter<std::string>("elHLT_MC")),
lepMonitorHLT   (iConfig.getParameter<std::string>("lepMonitorHLT")),

low_pt32_elHLT          (iConfig.getParameter<std::string>("low_pt32_elHLT")),
low_pt28_150HT_elHLT    (iConfig.getParameter<std::string>("low_pt28_150HT_elHLT")),
low_pt30_35PFJet_elHLT  (iConfig.getParameter<std::string>("low_pt30_35PFJet_elHLT")),
elmuHLT_1               (iConfig.getParameter<std::string>("elmuHLT_1")),
elmuHLT_2               (iConfig.getParameter<std::string>("elmuHLT_2")),
elmuHLT_3               (iConfig.getParameter<std::string>("elmuHLT_3")),
elmuHLT_4               (iConfig.getParameter<std::string>("elmuHLT_4")),
elelHLT_1               (iConfig.getParameter<std::string>("elelHLT_1")),
elelHLT_2               (iConfig.getParameter<std::string>("elelHLT_2")),

eltauHLT               (iConfig.getParameter<std::string>("eltauHLT")),
mutauHLT1              (iConfig.getParameter<std::string>("mutauHLT1")),
mutauHLT2              (iConfig.getParameter<std::string>("mutauHLT2")),

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
btag_threshold          (iConfig.getParameter<double>("btag_threshold")),
btag_threshold_Medium   (iConfig.getParameter<double>("btag_threshold_Medium")),
goodLumiFilter   (iConfig.getUntrackedParameter<std::vector<edm::LuminosityBlockRange>>("lumisToProcess", std::vector<edm::LuminosityBlockRange>())),
urls   (iConfig.getUntrackedParameter <std::vector <std::string> >("input")),
outUrl (iConfig.getParameter<std::string>("outfile")),
//BadChCandFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"))),
//BadPFMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuonFilter")))
// here it breaks with
// MissingParameter: Parameter 'BadChargedCandidateFilter' not found.
//triggerObjects_ (consumes(iConfig.getParameter<edm::InputTag>("hlt_objects")) //"selectedPatTrigger"))
// error: no matching function for call to 'NtuplerAnalyzer::consumes(edm::InputTag)'
triggerObjects_InputTag (iConfig.getParameter<edm::InputTag>("hlt_objects"))

{
	edm::LogInfo ("Demo") << "building the object";

	r3 = new TRandom3();

	/* define in constructor via call to consumes (magic thingy) */
	//tracks_    = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
	muons_     = consumes<pat::MuonCollection>    (edm::InputTag("slimmedMuons"));
	electrons_ = consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"));
	//taus_ = consumes<pat::TauCollection>(edm::InputTag("slimmedTaus"));
	//taus_ = consumes<pat::TauCollection>(edm::InputTag("NewTauIDsEmbedded")); // the taus with embedded reprocessed 2017v2 tau IDs for 2016 legacy
	taus_ = consumes<pat::TauCollection>(edm::InputTag(iConfig.getParameter<std::string>("tau_objs_name")));

	//string ("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017");
	tau_VLoose_ID = iConfig.getParameter<std::string>("tau_VLoose_ID");
	tau_Loose_ID  = iConfig.getParameter<std::string>("tau_Loose_ID" );
	tau_Medium_ID = iConfig.getParameter<std::string>("tau_Medium_ID");
	tau_Tight_ID  = iConfig.getParameter<std::string>("tau_Tight_ID" );
	tau_VTight_ID = iConfig.getParameter<std::string>("tau_VTight_ID");
	//string ("byVVTightIsolationMVArun2017v2DBoldDMwLT2017");

	vtx_ = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));
	rho_ = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
	// declare consuming the HLT to be able to get triggers in the following
	trigResults_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","", HLT_source));
	edm::LogInfo ("Demo") << "HLT results";
	//trigResults2_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT2"));
	trigResultsRECO_    = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","RECO"));
	trigResultsPAT_     = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","PAT"));
	//triggerObjects_ = consumes<vector<pat::TriggerObjectStandAlone>>(edm::InputTag("selectedPatTrigger"));
	//triggerObjects_ = consumes<vector<pat::TriggerObjectStandAlone>>(edm::InputTag("slimmedPatTrigger")); // damn!!!!!!!!!!!!!
	triggerObjects_ = consumes<vector<pat::TriggerObjectStandAlone>>(triggerObjects_InputTag); // damn!!!!!!!!!!!!!

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
	genLepsToken_ = consumes<std::vector<reco::GenJet> >(edm::InputTag("particleLevel:leptons"));
	genJetsToken_ = consumes<std::vector<reco::GenJet> >(edm::InputTag("particleLevel:jets")); // not sure if it is different from previous one
	genFatJetsToken_ = consumes<std::vector<reco::GenJet> >(edm::InputTag("particleLevel:fatjets"));
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

	edm::LogInfo ("Demo") << "AAA";

	// dtag configs
	bool period_2016BCD = !isMC && (dtag.Contains("2016B") || dtag.Contains("2016C") || dtag.Contains("2016D"));
	bool period_2016EF  = !isMC && (dtag.Contains("2016E") || dtag.Contains("2016F"));
	bool period_2016G   = !isMC && (dtag.Contains("2016G"));
	bool period_2016H   = !isMC && (dtag.Contains("2016H"));
	is2017legacy = dtag.Contains("2017legacy");
	is2017data   = !isMC && (dtag.Contains("2017")); // TODO check 2017 data handling with dtag
	bool period_2017legacyB  = !isMC &&  dtag.Contains("2017B");
	bool period_2017legacyC  = !isMC &&  dtag.Contains("2017C");
	bool period_2017legacyDE = !isMC && (dtag.Contains("2017D") || dtag.Contains("2017E"));
	bool period_2017legacyF  = !isMC &&  dtag.Contains("2017F");

	aMCatNLO = dtag.Contains("amcatnlo");
	isWJets = dtag.Contains("WJet") || dtag.Contains("W0Jet") || dtag.Contains("W1Jet") || dtag.Contains("W2Jet") || dtag.Contains("W3Jet") || dtag.Contains("W4Jet");
	isDY = dtag.Contains("DYJet");
	isTT = dtag.Contains("TT");
	isSingleTop = dtag.Contains("SingleT") || dtag.Contains("tchannel") || dtag.Contains("schannel");

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
	if (is2016legacy)
		{
		// Summer16_07Aug2017BCD_V11_DATA.tar.gz  Summer16_07Aug2017EF_V11_DATA.tar.gz  Summer16_07Aug2017GH_V11_DATA.tar.gz  Summer16_07Aug2017_V11_MC.tar.gz
		if (isMC)
			jet_corr_files = "/Summer16_07Aug2017_V11_MC";
		else if (period_2016BCD)
			jet_corr_files = "/Summer16_07Aug2017BCD_V11_DATA";
		else if (period_2016EF)
			jet_corr_files = "/Summer16_07Aug2017EF_V11_DATA";
		else if (period_2016G || period_2016H)
			jet_corr_files = "/Summer16_07Aug2017GH_V11_DATA";
		//else jet_corr_files = "/Summer16_07Aug2017GH_V11_DATA";
		}

	else if (is2017legacy)
		{
		// Fall17_17Nov2017B_V32_DATA.tar.gz
		// Fall17_17Nov2017C_V32_DATA.tar.gz
		// Fall17_17Nov2017DE_V32_DATA.tar.gz
		// Fall17_17Nov2017F_V32_DATA.tar.gz
		// Fall17_17Nov2017_V32_MC.tar.gz
		if (isMC)
			jet_corr_files = "Fall17_17Nov2017_V32_MC";
		else if (period_2017legacyB)
			jet_corr_files = "Fall17_17Nov2017B_V32_DATA";
		else if (period_2017legacyC)
			jet_corr_files = "Fall17_17Nov2017C_V32_DATA";
		else if (period_2017legacyDE)
			jet_corr_files = "Fall17_17Nov2017DE_V32_DATA";
		else if (period_2017legacyF)
			jet_corr_files = "Fall17_17Nov2017F_V32_DATA";
		}

	else
		{
		// original 2016 rereco, Moriond17
		if (isMC)
			jet_corr_files = "/Summer16_23Sep2016V4_MC";
		else if (period_2016BCD)
			jet_corr_files = "/Summer16_23Sep2016BCDV4_DATA";
		else if (period_2016EF)
			jet_corr_files = "/Summer16_23Sep2016EFV4_DATA";
		else if (period_2016G)
			jet_corr_files = "/Summer16_23Sep2016GV4_DATA";
		else if (period_2016H)
			jet_corr_files = "/Summer16_23Sep2016HV4_DATA";
		else jet_corr_files = "/Summer16_23Sep2016HV4_DATA";
		}

	jesCor = utils::cmssw::getJetCorrector (jecDir, jet_corr_files, isMC);
	totalJESUnc = new JetCorrectionUncertainty ((jecDir + jet_corr_files + "_Uncertainty_AK4PFchs.txt").Data());

	edm::LogInfo ("Demo") << "BBB";

	// resolution and scale-factors for the systematics
	gSystem->ExpandPathName(TjetResolutionFileName);
	gSystem->ExpandPathName(TjetResolutionSFFileName);

	string jetResolutionFileName   (TjetResolutionFileName);
	string jetResolutionSFFileName (TjetResolutionSFFileName);
	// <---- ROOT & CMSSW are best friends
	jet_resolution_in_pt = JME::JetResolution(jetResolutionFileName);
	jet_resolution_sf_per_eta = JME::JetResolutionScaleFactor(jetResolutionSFFileName);

	edm::LogInfo ("Demo") << "CCC";

	edm::Service<TFileService> fs;
	// init ttree
	NT_output_ttree = fs->make<TTree>("reduced_ttree", "TTree with reduced event data");
	event_counter      = fs->make<TH1D>( "events_counter"     , "pass category", 100,  0, 100);
	weight_counter     = fs->make<TH1D>( "weight_counter"     , "pass category", 100,  0, 100);
	systematic_weights = fs->make<TH1D>( "systematic_weights" , "pass category", 200,  0, 200);

	/* We want to save gen level distributions
	 * for particular sub-processes, i.e. the final states,
	 * in each dtag, i.e. the hard process.
	 * Let's define the histograms per the leptonic part of the final state:
	 * elel, elmu, mumu, eltau1h, eltau3h, mutau1h, mutau3h, el and mu (for el+jets and mu+jets), tau1h and tau3h (for W+jets->tau+jets).
	 * In principle, the final state definition might include the hadronic part of the final state.
	 * But we do not need it now, for the current processes under study.
	 * We can add it later.

	 * This way we will cover all dtags (ttbar, DY, and W+jets) in one go:
	 * define an ID of the final state calculated in each dtag differently,
	 * save the corresponding histograms.
	 */

	GenDistrs_make_histos_in_FileService(fs);

	//ttbar_elmu_lep_pt = fs->make<TH1D>( "ttbar_lep_pt" , "ttbar_lep_pt", 200,  0, 200);
	// for control of effect from PU reweighting, aMCatNLO gen -1 weights, top pt

	/*
	 * add more bins for various systematic weights in case they are shape-systematics:
	 * PDF, fragmentation, QCD scale etc
	 * -- these per-event weights are normalized to some "central" weight
	 *    does it cover the event weight sum normalization?
	 *    anyway add these counters just in case
	 *
	 * I need nominal weight of all events
	 * and this same weight * various systematics
	 * 1) what's nominal?
	 * 2) which weight-based systematics are there?
	 *
	 * 1) nominal = PU * aMCatNLO-1 * z_mass_pt (DY)
	 *    pu for EL and MU distrs
	 *    -- z_mass_pt is out, hence there are 2 nominals, for MU and EL PU
	 *
	 * 2) systematics:
	 *    pu up/down
	 *    top pt
	 *    fragmentations
	 *    pdfs and alphaS (with constrain of 2)?
	 *    qcd renorm/refact
	 *
	 * it looks like weight_aMCatNLO * weight_PU * weight_th_sys
	 * they all are orthogonal
	 * therefore, just save them separately
	 * find the variation with respect to N of events
	 * compensate for it
	 *
	 * and hardcode the bin numbers for reference.. make enum for this?
	 */

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
	#define NTUPLE_INTERFACE_CLASS_INITIALIZE
	#include "UserCode/NtuplerAnalyzer/interface/ntupleOutput.h"
	#undef NTUPLE_INTERFACE_CLASS_INITIALIZE

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

	edm::LogInfo ("Demo") << "done";

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
	#define NTUPLE_INTERFACE_CLASS_RESET
	#include "UserCode/NtuplerAnalyzer/interface/ntupleOutput.h"
	#undef NTUPLE_INTERFACE_CLASS_RESET
	// if output contains stand-alone objects (not vector of LorentzVector-s, but just 1 LorentzVector, like MET or something)
	// you have to reset them yourself, since each object might have its' own method
	NT_met_init.SetXYZT(0,0,0,0);

	NT_met_init_shift_UnclusteredEnUp   .SetXYZT(0,0,0,0);
	NT_met_init_shift_UnclusteredEnDown .SetXYZT(0,0,0,0);
	NT_met_init_shift_JetEnUp           .SetXYZT(0,0,0,0);
	NT_met_init_shift_JetEnDown         .SetXYZT(0,0,0,0);
	NT_met_init_shift_JetResUp          .SetXYZT(0,0,0,0);
	NT_met_init_shift_JetResDown        .SetXYZT(0,0,0,0);
	NT_met_init_shift_MuonEnUp          .SetXYZT(0,0,0,0);
	NT_met_init_shift_MuonEnDown        .SetXYZT(0,0,0,0);
	NT_met_init_shift_ElectronEnUp      .SetXYZT(0,0,0,0);
	NT_met_init_shift_ElectronEnDown    .SetXYZT(0,0,0,0);
	NT_met_init_shift_TauEnUp           .SetXYZT(0,0,0,0);
	NT_met_init_shift_TauEnDown         .SetXYZT(0,0,0,0);

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

	NT_gen_decay_lep1_p4.SetXYZT(0,0,0,0);
	NT_gen_decay_lep2_p4.SetXYZT(0,0,0,0);

	NT_gen_decay_bjet1_p4.SetXYZT(0,0,0,0);
	NT_gen_decay_bjet2_p4.SetXYZT(0,0,0,0);

	NT_gen_decay_jet1_p4.SetXYZT(0,0,0,0);
	NT_gen_decay_jet2_p4.SetXYZT(0,0,0,0);

	NT_gen_decay_missing_p4.SetXYZT(0,0,0,0);

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

	/*
	 * weird cms software needs pat jets to get the jet flavour of reco::GenJet
	 */

	// jets
	pat::JetCollection jets;
	edm::Handle<pat::JetCollection>jetsHandle;
	//jetsHandle.getByLabel(ev, "slimmedJets");
	iEvent.getByToken(jets_, jetsHandle);
	if(jetsHandle.isValid() ) jets = *jetsHandle;

	// for matching reco to gen I need
	// leptons (11, 13), taus (15) and tau od DM10 (3ch) from hard processes or from tt decay in case of tt
	// and visible products of W and b in case of tt
	// then a reco object will be matched to gen collections
	// the match defined by not overlapping sum:
	// lep tau tau3ch b    W
	// 1   2   4      8   16
	//vector<const reco::Candidate*> gen_leps, gen_taus, gen_tau3ch, gen_w_prods, gen_b_prods;
	vector<LorentzVector> gen_leps, gen_taus, gen_tau3ch, gen_taulep, gen_w_prods, gen_b_prods, gen_tt_tau_prods;
	// parallel provenance ids
	// save the corresponding ID of the particle
	// -- to match b to W+ and b-bar to W- etc
	vector<int>           gid_leps, gid_taus, gid_tau3ch, gid_taulep, gid_w_prods, gid_b_prods;

	bool isTTSignal = false;

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
			NT_gen_weight_norm = nomlheweight; // this one also = 1
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

			// TODO: substitute this with another enum
			//NT_gen_weights_renorm_fact.push_back(muf_down_mur_down_weight); // 0
			//NT_gen_weights_renorm_fact.push_back(muf_down_mur_nom_weight ); // 1
			//NT_gen_weights_renorm_fact.push_back(muf_down_mur_up_weight  ); // 2
			//NT_gen_weights_renorm_fact.push_back(muf_nom_mur_down_weight ); // 3
			//NT_gen_weights_renorm_fact.push_back(muf_nom_mur_nom_weight  ); // 4
			//NT_gen_weights_renorm_fact.push_back(muf_nom_mur_up_weight   );
			//NT_gen_weights_renorm_fact.push_back(muf_up_mur_down_weight  );
			//NT_gen_weights_renorm_fact.push_back(muf_up_mur_nom_weight   );
			//NT_gen_weights_renorm_fact.push_back(muf_up_mur_up_weight    ); // 8

			//LogInfo("Demo") << "renorm refact len " << NT_gen_weights_renorm_fact.size();
			NT_gen_weights_renorm_fact.resize(9);
			//LogInfo("Demo") << "renorm refact len " << NT_gen_weights_renorm_fact.size();
			NT_gen_weights_renorm_fact[MUf_down_MUr_down ]  = muf_down_mur_down_weight ;
			NT_gen_weights_renorm_fact[MUf_down_MUr_nom  ]  = muf_down_mur_nom_weight  ;
			NT_gen_weights_renorm_fact[MUf_down_MUr_up   ]  = muf_down_mur_up_weight   ;
			NT_gen_weights_renorm_fact[MUf_nom_MUr_down  ]  = muf_nom_mur_down_weight  ;
			NT_gen_weights_renorm_fact[MUf_nom_MUr_nom   ]  = muf_nom_mur_nom_weight   ;
			NT_gen_weights_renorm_fact[MUf_nom_MUr_up    ]  = muf_nom_mur_up_weight    ;
			NT_gen_weights_renorm_fact[MUf_up_MUr_down   ]  = muf_up_mur_down_weight   ;
			NT_gen_weights_renorm_fact[MUf_up_MUr_nom    ]  = muf_up_mur_nom_weight    ;
			NT_gen_weights_renorm_fact[MUf_up_MUr_up     ]  = muf_up_mur_up_weight     ;
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
			// according to
			// https://lhapdf.hepforge.org/pdfsets
			// 13100 is CT14nlo and has 57 members
			// after these in the weights are (see the printout with names)
			// 13164 -- CT14nlo_as_0117
			// 13166 -- CT14nlo_as_0119
			// these are variation of alpha S
			// and after them is
			// 11000 -- CT10nlo with 53 parameters
			// but there are no parameters,
			// it is immediatly folowed by 25200,25201 etc

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

			// CT14 pdfs have 57 members: 0th nominal + 56 variations

			//get the original mc replica weights
			//unsigned int nNNPPDF3Weights_ = 100, pdfWeightOffset_ = 9, nPdfEigWeights_ = 60; // as in example
			// -- for NNPDF3
			// save nominal NNPDF30 weight
			//NT_gen_weight_pdf_nn30_nominal = lheEPHandle->weights()[9].wgt;
			// -- there is no nominal value 260000, only 260001+

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

		#if defined(BFRAG_PROC)
		#if BFRAG_PROC==2016
		// gen2 leptons for signal acceptance info
		edm::Handle<std::vector<reco::GenJet>> genLeps2;
		iEvent.getByToken( genLepsToken_, genLeps2);

		for(auto genLep=genLeps2->begin(); genLep!=genLeps2->end(); ++genLep)
			{
			// this jet should have some info on whether it's b jet
			NT_gen2_leptons_pdgId.push_back(genLep->pdgId());
			NT_gen2_leptons_p4.push_back(genLep->p4());
			}
		#endif
		#endif

		//// fatjets just in case
		//edm::Handle<std::vector<reco::GenJet>> genFatJets2;
		//iEvent.getByToken( genFatJetsToken_, genFatJets2);
		//for(auto genFatJet=genFatJets2->begin(); genFatJet!=genFatJets2->end(); ++genFatJet)
		//	{
		//	// this jet should have some info on whether it's b jet
		//	NT_gen2_fatjets_pdgId.push_back(genFatJet->pdgId());
		//	}
		#if defined(BFRAG_PROC)
		#if BFRAG_PROC==2016

		// fragmentation and decay tables (of b->hadron) systematics
		edm::Handle<std::vector<reco::GenJet>> genJets2;
		iEvent.getByToken( genJetsToken_, genJets2);
		// https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer
		// in particular:
		// > in the analyzer method get the genJets,
		// > the weights and loop over the jets to analyse them,
		// > and take the product of the jet weights as event weight
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
			// this jet should have some info on whether it's b jet
			NT_gen2_jets_pdgId.push_back(genJet->pdgId());
			NT_gen2_jets_p4.push_back(genJet->p4());

			// get matching to gen leptons
			Float_t dR_to_sel = 999.;
			for(size_t l=0; l<NT_gen2_leptons_p4.size(); ++l)
				{
				double dR = reco::deltaR(genJet->p4(), NT_gen2_leptons_p4[l]);
				if (dR < dR_to_sel)
					dR_to_sel = dR;
				}
			NT_gen2_jets_lep_dR.push_back(dR_to_sel);
			NT_gen2_jets_lep_dR_matched.push_back(dR_to_sel < 0.4);

			edm::Ref<std::vector<reco::GenJet> > genJetRef(genJets2, genJet - genJets2->begin()); // this looks really weird
			//cout << "pt=" << genJet->pt() << " id=" << genJet->pdgId() << " petersonFragWeight=" << (*petersonFrag)[genJetRef] << endl;
			//...
			//Check if the jet is b-jet
			//sadly, these genJets are reco::GenJet
			//they don't have partonFlavour, only pat jets have it
			//match within dR < 0.4 to pat jets, get the flavour
			float min_dR_to_patJet = 999., dR_max = 0.4;
			unsigned int min_dR_partonFlavour = 0;
			for(auto patJet=jets.begin(); patJet!=jets.end(); ++patJet)
				{
				// whatchout, the deltaR does not get through iterators..
				float dR_to_patJet = reco::deltaR(*patJet, *genJet);
				if (dR_to_patJet < dR_max && dR_to_patJet < min_dR_to_patJet)
					{
					min_dR_to_patJet = dR_to_patJet;
					min_dR_partonFlavour = abs(patJet->partonFlavour());
					// break; maybe?
					}
				}

			// these systematic weights are only for b-jets
			if (min_dR_partonFlavour != 5) continue;

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
		#endif

		LogInfo("Demo") << "MC systematic weights for jet fragmentation";
		#endif

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

				/* these don't work well -- brute force all no-daughter genparts
				*/
				// new homogeneous gen final states
				if (p.fromHardProcessFinalState() || p.isPromptFinalState() || p.isDirectHardProcessTauDecayProductFinalState())
					{
					LogInfo ("Demo") << "got hard process final state";
					NT_gen_final_p4              .push_back(p.p4());
					NT_gen_final_PromptFinal     .push_back(p.isPromptFinalState());
					NT_gen_final_fromHardFinal   .push_back(p.fromHardProcessFinalState());
					NT_gen_final_PromptTauDecay  .push_back(p.isDirectHardProcessTauDecayProductFinalState());
					NT_gen_final_pdgId           .push_back(p.pdgId());
					NT_gen_final_status          .push_back(p.status());
					NT_gen_final_ndaughters      .push_back(p.numberOfDaughters());
					NT_gen_final_chainId         .push_back(parse_chain_id(p));
					LogInfo ("Demo") << "saved hard process final state";
					}

				// prompt decayed are b-mesons and w-jet hadrons and tau stuff
				// save them separately
				/*
				if (p.isPromptDecayed())
					{
					NT_gen_PromptDecayed_p4      .push_back(p.p4());
					NT_gen_PromptDecayed_pdgId   .push_back(p.pdgId());
					NT_gen_PromptDecayed_status  .push_back(p.status());
					NT_gen_PromptDecayed_ndaughters  .push_back(p.numberOfDaughters());
					NT_gen_PromptDecayed_chainId .push_back(parse_chain_id(&p));
					// and save their final states
					save_final_states(&p, decay_final_states);
					for (size_t f_i=0; f_i < decay_final_states.size(); f_i++)
						{
						NT_gen_final_p4              .push_back(decay_final_states[f_i]->p4());
						// and cms awesomeness... the candidates don't have all the following methods
						// call for daughter returns candidates but GenParticles are needed
						///
						NT_gen_final_PromptFinal     .push_back(((const reco::GenParticle *) decay_final_states[f_i])->isPromptFinalState());
						NT_gen_final_fromHardFinal   .push_back(((const reco::GenParticle *) decay_final_states[f_i])->fromHardProcessFinalState());
						NT_gen_final_PromptTauDecay  .push_back(((const reco::GenParticle *) decay_final_states[f_i])->isDirectHardProcessTauDecayProductFinalState());
						NT_gen_final_pdgId           .push_back(((const reco::GenParticle *) decay_final_states[f_i])->pdgId());
						NT_gen_final_status          .push_back(((const reco::GenParticle *) decay_final_states[f_i])->status());
						NT_gen_final_ndaughters      .push_back(decay_final_states[f_i]->numberOfDaughters());
						NT_gen_final_chainId         .push_back(parse_chain_id(decay_final_states[f_i]));
						}
					}
				*/

				/* brute force does not work due to collor reconnection crap and it is not obvious that a simple hack will solve it, a complex hack is too messy
				if (p.numberOfDaughters() == 0)
					{
					Int_t chain_id = parse_chain_id(&p);
					if (chain_id != 0)
						{
						NT_gen_final_p4              .push_back(p.p4());
						NT_gen_final_PromptFinal     .push_back(p.isPromptFinalState());
						NT_gen_final_fromHardFinal   .push_back(p.fromHardProcessFinalState());
						NT_gen_final_PromptTauDecay  .push_back(p.isDirectHardProcessTauDecayProductFinalState());
						NT_gen_final_pdgId           .push_back(p.pdgId());
						NT_gen_final_status          .push_back(p.status());
						//NT_gen_final_ndaughters      .push_back(p.numberOfDaughters());
						NT_gen_final_chainId         .push_back(chain_id);
						}
					}
				*/

				// leptons from hard processes
				// basically, this is needed only for DY -- WJets do always have the W, which is caught later
				if ((a_id >= 11 && a_id <= 16 && p.fromHardProcessFinalState()) ||
					//(p.isDirectHardProcessTauDecayProduct()))
					(p.isDirectHardProcessTauDecayProductFinalState())) // same stuff
					{
					LogInfo ("Demo") << "got hard process lepton";
					// saving the hard leptons from here if it is not TT
					// apparently single top and WJets also have explicit W in the decay tree
					// the W is caught later
					if (!(isTT || isSingleTop || isWJets))
						{
						if (a_id == 11 || a_id == 13)
							{
							save_final_cands(&p, gen_leps, gid_leps, p.pdgId());
							NT_gen_match_lep_id.push_back(id);
							}
						// if it is tau -- check if it is DM10+, i.e. decay to 3 charged particles
						else if (a_id == 15)
							{
							int tau_id = simple_tau_decay_id(&p);
							// = 11, 13 for leptons and 20 + 5*(Nch-1) + Npi0 for hadrons
							// 20 + 10
							if (abs(tau_id) >= 30)
								{
								save_final_cands(&p, gen_tau3ch, gid_tau3ch, p.pdgId());
								NT_gen_match_tau3ch_id.push_back(id);
								}
							else if (abs(tau_id) >= 15)
								{
								save_final_cands(&p, gen_taus, gid_taus, p.pdgId());
								NT_gen_match_tau_id.push_back(id);
								}
							else
								{
								save_final_cands(&p, gen_taulep, gid_taulep, p.pdgId());
								NT_gen_match_tau_id.push_back(id);
								}
							}
						else
							{
							const reco::Candidate* tau_mother = has_tau_mother(&p, 5);
							if (tau_mother != NULL)
								{
								int tau_id = simple_tau_decay_id(tau_mother);
								if (abs(tau_id) >= 30)
									{
									save_final_cands(tau_mother, gen_tau3ch, gid_tau3ch, p.pdgId());
									NT_gen_match_tau3ch_id.push_back(id);
									}
								else if (abs(tau_id) >= 15)
									{
									save_final_cands(tau_mother, gen_taus, gid_taus, p.pdgId());
									NT_gen_match_tau_id.push_back(id);
									}
								else
									{
									save_final_cands(tau_mother, gen_taulep, gid_taulep, p.pdgId());
									NT_gen_match_tau_id.push_back(id);
									}
								}
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

					if (isDY && a_id == 15)
						{
						int tau_id = simple_tau_decay_id(&p);
						// save the gen tau (final in the tau-tau modeling chain?)
						NT_gen_tt_tau_orig_p4.push_back(p.p4());
						NT_gen_tt_tau_simpleID.push_back(tau_id);

						// save the visible part of the tt tau (all gen taus grab not-tt stuff
						LorentzVector gen_tt_tau_vis(0,0,0,0);
						sum_final_cands(&p, gen_tt_tau_prods, gen_tt_tau_vis, true);
						NT_gen_tt_tau_vis_p4.push_back(gen_tt_tau_vis);
						// save invis part (might not be present)
						LorentzVector gen_tt_tau_invis(0,0,0,0);
						sum_final_cands(&p, gen_tt_tau_prods, gen_tt_tau_invis, false);
						NT_gen_tt_tau_invis_p4.push_back(gen_tt_tau_invis);
						}
					LogInfo ("Demo") << "processed hard process lepton";
					}

				if (abs(id) == 6 && n_daughters == 2)
					// if it is a t quark
					// it is a decay vertex of t to something
					// (could use p.isLastCopy())
					{
					LogInfo ("Demo") << "got decayed t";
					// calculate top_pt weights:
					weight_TopPT *= TMath::Sqrt(top_pT_SF(p.pt()));
					LogInfo ("Demo") << "got decayed t pt";
					// find the W decay channel in this top
					unsigned int d0_id = abs(p.daughter(0)->pdgId());
					LogInfo ("Demo") << "got decayed daughter 0 pdgId = " << d0_id;
					unsigned int d1_id = abs(p.daughter(1)->pdgId());
					LogInfo ("Demo") << "got decayed daughter 1 pdgId = " << d1_id;
					int W_num = d0_id == 24 ? 0 : (d1_id == 24 ? 1 : -1) ;
					if (W_num < 0) continue;
					const reco::Candidate * W = p.daughter( W_num );
					const reco::Candidate * b = p.daughter( 1 - W_num );
					LogInfo ("Demo") << "got decayed fixed daughters";
					const reco::Candidate * W_final = find_W_decay(W);
					LogInfo ("Demo") << "t " << id << " decay W "  << W_final ->pdgId() << " final pointer " << W_final << " #daughters = " << W_final->numberOfDaughters();
					LogInfo ("Demo") << "t " << id << " decay b  " << b       ->pdgId() << " pointer "       << b       << " #daughters = " << b->numberOfDaughters();

					int decay_id = 1;
					// = id of lepton or 1 for quarks
					for (unsigned int d_i = 0; d_i < 2; d_i++)
						{
						int d_i_pdgId = W_final->daughter(d_i)->pdgId();
						LogInfo ("Demo") << "W boson daughter " << d_i << " " << d_i_pdgId;
						if (fabs(d_i_pdgId) == 11 || fabs(d_i_pdgId) == 13)
							{
							decay_id = d_i_pdgId;
							save_final_cands(W_final->daughter(d_i), gen_leps, gid_leps, d_i_pdgId);
							}
						if (fabs(d_i_pdgId) == 15)
							{
							int tau_id = simple_tau_decay_id(W_final->daughter(d_i));
							// save the gen tau (final in the tau-tau modeling chain?)
							NT_gen_tt_tau_orig_p4.push_back(W_final->daughter(d_i)->p4());
							NT_gen_tt_tau_simpleID.push_back(tau_id);

							// save the visible part of the tt tau (all gen taus grab not-tt stuff
							LorentzVector gen_tt_tau_vis(0,0,0,0);
							sum_final_cands(W_final->daughter(d_i), gen_tt_tau_prods, gen_tt_tau_vis, true);
							NT_gen_tt_tau_vis_p4.push_back(gen_tt_tau_vis);
							// save invis part (might not be present)
							LorentzVector gen_tt_tau_invis(0,0,0,0);
							sum_final_cands(W_final->daughter(d_i), gen_tt_tau_prods, gen_tt_tau_invis, false);
							NT_gen_tt_tau_invis_p4.push_back(gen_tt_tau_invis);

							// = 11, 13 for leptons and 20 + 5*(Nch-1) + Npi0 for hadrons
							// 20 + 10
							LogInfo ("Demo") << "tau decay w ID " << tau_id;
							if (abs(tau_id) >= 30)
								{
								save_final_cands(W_final->daughter(d_i), gen_tau3ch, gid_tau3ch, d_i_pdgId);
								LogInfo ("Demo") << "saved to 3h taus";
								}
							else if (abs(tau_id) >= 15)
								{
								save_final_cands(W_final->daughter(d_i), gen_taus, gid_taus, d_i_pdgId);
								LogInfo ("Demo") << "saved to  h taus";
								}
							else
								{
								save_final_cands(W_final->daughter(d_i), gen_taulep, gid_taulep, d_i_pdgId);
								LogInfo ("Demo") << "saved to l taus";
								}
							decay_id = d_i_pdgId * abs(tau_id); // to keep the sign of the tau take it in abs
							}
						}
					LogInfo ("Demo") << "processed leptonic W";
					// if W is not leptonic, and decay id is still = 1
					if (decay_id == 1)
						{
						LogInfo ("Demo") << "W was not leptonic";
						save_final_cands(W_final, gen_w_prods, gid_w_prods, W_final->pdgId());
						LogInfo ("Demo") << "processed hadronic W";
						//NT_gen_match_w_id.push_back(id); // TODO: remake the saving important gen particloes
						}
					save_final_cands(b, gen_b_prods, gid_b_prods, b->pdgId());
					LogInfo ("Demo") << "processed b";

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
					LogInfo ("Demo") << "processed t";
					}

				// if it is a tau -- save the non-neutrino part to the output
				//  the status is 1 or 2
				//  1. final state, not decays, so it should never happen for tau
				//  2. decayed or fragmented -- the case for tau
				if (abs(id) == 15)
					{
					LogInfo ("Demo") << "got a tau " << st;
					if (st == 1)
						NT_gen_tau_p4.push_back(p.p4()); 
					else if (st == 2)
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
					LogInfo ("Demo") << "processed tau";
					}

				// Prompt leptons ID for Z->LL
				// pythia 8 stores prompt particles with status 21-29 ("hardest subprocess", PYTHIA 8 Worksheet for tutorial at ASP 2012 Summer School)
				// -- checked DYJetsToLL -- there is no Z (pdgId 23) particles, but prompt leptons work fine
				// thus save N prompt leptons in the process
				// and their ID
				if ((abs(id) == 11 || abs(id) == 13 || abs(id) == 15) && st > 20 && st < 30)
					{
					LogInfo ("Demo") << "got prompt lepton";
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
					LogInfo ("Demo") << "got Z";
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
					// to suite single-top
					// skip the W from top decay -- it's picked up above
					LogInfo ("Demo") << "consider W with mothers: " << p.numberOfMothers();
					if (has_mother(p, 6)) continue;
					// and save the final states from this W
					// if it is "normal" final state daughters:
					// leptons, light jets or b, but not top
					// b is for the "exotic" single top s/t channels
					// it should catch tW correctly, other channels don't matter now

					int wdecay_id = 1;
					int d0_id = abs(p.daughter(0)->pdgId());
					int d1_id = abs(p.daughter(1)->pdgId());
					LogInfo ("Demo") << "processing standalone W " << d0_id << ' ' << d1_id;

					int lep_daughter = (d0_id == 11 || d0_id == 13 || d0_id == 15 ? 0 : (d1_id == 11 || d1_id == 13 || d1_id == 15 ? 1 : -1));
					if (lep_daughter >= 0)
						{
						wdecay_id = p.daughter(lep_daughter)->pdgId();
						if (abs(wdecay_id) == 15)
							{
							int tau_id = simple_tau_decay_id(p.daughter(lep_daughter));
							// = 11, 13 for leptons and 20 + 5*(Nch-1) + Npi0 for hadrons
							// 20 + 10
							if (abs(tau_id) >= 30)
								save_final_cands(p.daughter(lep_daughter), gen_tau3ch, gid_tau3ch, wdecay_id);
							else if (abs(tau_id) >= 15)
								save_final_cands(p.daughter(lep_daughter), gen_taus, gid_taus, wdecay_id);
							else
								save_final_cands(p.daughter(lep_daughter), gen_taulep, gid_taulep, wdecay_id);
							wdecay_id *= tau_id;
							}
						else
							{
							save_final_cands(p.daughter(lep_daughter), gen_leps, gid_leps, wdecay_id);
							}
						}
					// the t/s channel case of b-jet
					// in WZTo2LTo2Q
					// I get decay of:
					// processing standalone W 5 4
					// and it comes direclty from up-down collision in protons (m 1 m 2)
					// -- it is the W of the process and it decays into b+s
					// however, W cannot decay into 2 b-s, hence if else
					if (d0_id == 5)
						save_final_cands(p.daughter(0), gen_b_prods, gid_b_prods, p.daughter(0)->pdgId());
					else if (d1_id == 5)
						save_final_cands(p.daughter(1), gen_b_prods, gid_b_prods, p.daughter(1)->pdgId());

					// normal W->jets
					if (d0_id < 5)
						{
						save_final_cands(p.daughter(0), gen_w_prods, gid_w_prods, id);
						}
					if (d1_id < 5)
						{
						save_final_cands(p.daughter(1), gen_w_prods, gid_w_prods, id);
						}
					// TODO: check if it can go backwards in the t-channel tree?

					NT_gen_N_wdecays += 1;
					NT_gen_wdecays_IDs.push_back(wdecay_id);
					LogInfo ("Demo") << "done with W";
					}
				}

			if (isDY)
				{
				NT_genPt = genMomentum.Pt();
				NT_genMass = genMomentum.M();
				// weight in processing
				//NT_zPtWeight = zPtMass_histo->GetBinContent(zPtMass_histo->GetXaxis()->FindBin(NT_genMass), zPtMass_histo->GetYaxis()->FindBin(NT_genPt));
				}

			if (isTT)
				{
				// lepton+tauh decay
				isTTSignal |= (abs(NT_gen_t_w_decay_id) == 11 || abs(NT_gen_t_w_decay_id) == 13)   && (abs(NT_gen_tb_w_decay_id) > 15*15);
				isTTSignal |= (abs(NT_gen_tb_w_decay_id) == 11 || abs(NT_gen_tb_w_decay_id) == 13) && (abs(NT_gen_t_w_decay_id) > 15*15);
				isTTSignal |= record_all;

				// the general gen level IDs for ttbar
				NT_gen_decay_lep1_id = NT_gen_t_w_decay_id;
				NT_gen_decay_lep2_id = NT_gen_tb_w_decay_id;

				// 
				if (abs(NT_gen_decay_lep1_id) > 1) NT_gen_decay_lep1_p4 = (NT_gen_t_w1_final_p4 + NT_gen_t_w2_final_p4);
				if (abs(NT_gen_decay_lep1_id) > 1) NT_gen_decay_lep2_p4 = (NT_gen_tb_w1_final_p4 + NT_gen_tb_w2_final_p4);

				NT_gen_decay_bjet1_p4 = NT_gen_t_b_final_p4;
				NT_gen_decay_bjet2_p4 = NT_gen_tb_b_final_p4;
				// and no
				//NT_gen_decay_jet1_p4
				//NT_gen_decay_jet2_p4
				}
			LogInfo ("Demo") << "Found: t decay = " << NT_gen_t_w_decay_id << " ; tb decay = " << NT_gen_tb_w_decay_id << " is Sig " << isTTSignal;

			// save the gen-levela distributions
			// for the general final state IDs
			GenDistrs_record((struct GenDistrs_recorded_gen_objects) {
					{NT_gen_decay_lep1_id, NT_gen_decay_lep2_id},
					{&NT_gen_decay_lep1_p4, &NT_gen_decay_lep2_p4},
					//&NT_gen_decay_jet1_p4, &NT_gen_decay_jet2_p4,
					{&NT_gen_decay_bjet1_p4, &NT_gen_decay_bjet2_p4},
				});
			}

		if(NT_nvtx_gen <0 )        NT_nvtx_gen=0; // TODO check vertexes procedure
		if(NT_nvtx_gen >MAX_NVTX ) NT_nvtx_gen=MAX_NVTX;
		weight = pu_vector_NOMINAL[NT_nvtx_gen] * (aMCatNLO? weight_Gen : 1) * weight_TopPT;

		// the referense of all events
		// and the weights for the rate of each systematic
		// multiply the rates and compensate them in proc
		double nominal_syst_weight = (aMCatNLO? weight_Gen : 1);
		systematic_weights->Fill(NOM_EVENTS);
		systematic_weights->Fill(NOMINAL, nominal_syst_weight);
		if (NT_nvtx_gen < MAX_NVTX)
			{
			systematic_weights->Fill(MU_PU     , nominal_syst_weight * pu_vector_mu[NT_nvtx_gen]);
			systematic_weights->Fill(MU_PUUp   , nominal_syst_weight * pu_vector_mu_up[NT_nvtx_gen]);
			systematic_weights->Fill(MU_PUDown , nominal_syst_weight * pu_vector_mu_down[NT_nvtx_gen]);
			systematic_weights->Fill(El_PU     , nominal_syst_weight * pu_vector_el[NT_nvtx_gen]);
			systematic_weights->Fill(El_PUUp   , nominal_syst_weight * pu_vector_el_up[NT_nvtx_gen]);
			systematic_weights->Fill(El_PUDown , nominal_syst_weight * pu_vector_el_down[NT_nvtx_gen]);
			}
		// else the PU weight = 0

		if (isTT)
			{
			systematic_weights->Fill(TOPPT , nominal_syst_weight * weight_TopPT);

			// renorm refact scales
			systematic_weights->Fill(M_NOM   , nominal_syst_weight * NT_gen_weights_renorm_fact[MUf_nom_MUr_nom   ]);
			systematic_weights->Fill(MrUp    , nominal_syst_weight * NT_gen_weights_renorm_fact[MUf_nom_MUr_up    ]);
			systematic_weights->Fill(MrDown  , nominal_syst_weight * NT_gen_weights_renorm_fact[MUf_nom_MUr_down  ]);
			systematic_weights->Fill(MfUp    , nominal_syst_weight * NT_gen_weights_renorm_fact[  MUf_up_MUr_nom  ]);
			systematic_weights->Fill(MfDown  , nominal_syst_weight * NT_gen_weights_renorm_fact[MUf_down_MUr_nom  ]);
			systematic_weights->Fill(MfrUp   , nominal_syst_weight * NT_gen_weights_renorm_fact[MUf_up_MUr_up     ]);
			systematic_weights->Fill(MfrDown , nominal_syst_weight * NT_gen_weights_renorm_fact[MUf_down_MUr_down ]);
			// and these are fake
			//NT_gen_weights_renorm_fact[MUf_down_MUr_up   ]
			//NT_gen_weights_renorm_fact[MUf_up_MUr_down   ]

			// fragmentatopn
			systematic_weights->Fill(Central       , nominal_syst_weight * NT_gen_weight_centralFrag  );
			systematic_weights->Fill(FragUp        , nominal_syst_weight * NT_gen_weight_FragUp       );
			systematic_weights->Fill(FragDown      , nominal_syst_weight * NT_gen_weight_FragDown     );
			systematic_weights->Fill(PetersonUp    , nominal_syst_weight * NT_gen_weight_PetersonFrag );
			systematic_weights->Fill(SemilepBRUp   , nominal_syst_weight * NT_gen_weight_semilepbrUp  );
			systematic_weights->Fill(SemilepBRDown , nominal_syst_weight * NT_gen_weight_semilepbrDown);
			// PDF and alphaS
			systematic_weights->Fill(PDF_NOM       , nominal_syst_weight * NT_gen_weights_pdf_hessians[0]);
			systematic_weights->Fill(AlphaSUp      , nominal_syst_weight * NT_gen_weight_alphas_1);
			systematic_weights->Fill(AlphaSDown    , nominal_syst_weight * NT_gen_weight_alphas_2);

			// PDFs
			systematic_weights->Fill(PDFCT14n1     , nominal_syst_weight * NT_gen_weights_pdf_hessians[1]);
			systematic_weights->Fill(PDFCT14n2     , nominal_syst_weight * NT_gen_weights_pdf_hessians[2]);
			systematic_weights->Fill(PDFCT14n3     , nominal_syst_weight * NT_gen_weights_pdf_hessians[3]);
			systematic_weights->Fill(PDFCT14n4     , nominal_syst_weight * NT_gen_weights_pdf_hessians[4]);
			systematic_weights->Fill(PDFCT14n5     , nominal_syst_weight * NT_gen_weights_pdf_hessians[5]);
			systematic_weights->Fill(PDFCT14n6     , nominal_syst_weight * NT_gen_weights_pdf_hessians[6]);
			systematic_weights->Fill(PDFCT14n7     , nominal_syst_weight * NT_gen_weights_pdf_hessians[7]);
			systematic_weights->Fill(PDFCT14n8     , nominal_syst_weight * NT_gen_weights_pdf_hessians[8]);
			systematic_weights->Fill(PDFCT14n9     , nominal_syst_weight * NT_gen_weights_pdf_hessians[9]);
			systematic_weights->Fill(PDFCT14n10    , nominal_syst_weight * NT_gen_weights_pdf_hessians[10]);
			systematic_weights->Fill(PDFCT14n11    , nominal_syst_weight * NT_gen_weights_pdf_hessians[11]);
			systematic_weights->Fill(PDFCT14n12    , nominal_syst_weight * NT_gen_weights_pdf_hessians[12]);
			systematic_weights->Fill(PDFCT14n13    , nominal_syst_weight * NT_gen_weights_pdf_hessians[13]);
			systematic_weights->Fill(PDFCT14n14    , nominal_syst_weight * NT_gen_weights_pdf_hessians[14]);
			systematic_weights->Fill(PDFCT14n15    , nominal_syst_weight * NT_gen_weights_pdf_hessians[15]);
			systematic_weights->Fill(PDFCT14n16    , nominal_syst_weight * NT_gen_weights_pdf_hessians[16]);
			systematic_weights->Fill(PDFCT14n17    , nominal_syst_weight * NT_gen_weights_pdf_hessians[17]);
			systematic_weights->Fill(PDFCT14n18    , nominal_syst_weight * NT_gen_weights_pdf_hessians[18]);
			systematic_weights->Fill(PDFCT14n19    , nominal_syst_weight * NT_gen_weights_pdf_hessians[19]);
			systematic_weights->Fill(PDFCT14n20    , nominal_syst_weight * NT_gen_weights_pdf_hessians[20]);
			systematic_weights->Fill(PDFCT14n21    , nominal_syst_weight * NT_gen_weights_pdf_hessians[21]);
			systematic_weights->Fill(PDFCT14n22    , nominal_syst_weight * NT_gen_weights_pdf_hessians[22]);
			systematic_weights->Fill(PDFCT14n23    , nominal_syst_weight * NT_gen_weights_pdf_hessians[23]);
			systematic_weights->Fill(PDFCT14n24    , nominal_syst_weight * NT_gen_weights_pdf_hessians[24]);
			systematic_weights->Fill(PDFCT14n25    , nominal_syst_weight * NT_gen_weights_pdf_hessians[25]);
			systematic_weights->Fill(PDFCT14n26    , nominal_syst_weight * NT_gen_weights_pdf_hessians[26]);
			systematic_weights->Fill(PDFCT14n27    , nominal_syst_weight * NT_gen_weights_pdf_hessians[27]);
			systematic_weights->Fill(PDFCT14n28    , nominal_syst_weight * NT_gen_weights_pdf_hessians[28]);
			systematic_weights->Fill(PDFCT14n29    , nominal_syst_weight * NT_gen_weights_pdf_hessians[29]);
			systematic_weights->Fill(PDFCT14n30    , nominal_syst_weight * NT_gen_weights_pdf_hessians[30]);
			systematic_weights->Fill(PDFCT14n31    , nominal_syst_weight * NT_gen_weights_pdf_hessians[31]);
			systematic_weights->Fill(PDFCT14n32    , nominal_syst_weight * NT_gen_weights_pdf_hessians[32]);
			systematic_weights->Fill(PDFCT14n33    , nominal_syst_weight * NT_gen_weights_pdf_hessians[33]);
			systematic_weights->Fill(PDFCT14n34    , nominal_syst_weight * NT_gen_weights_pdf_hessians[34]);
			systematic_weights->Fill(PDFCT14n35    , nominal_syst_weight * NT_gen_weights_pdf_hessians[35]);
			systematic_weights->Fill(PDFCT14n36    , nominal_syst_weight * NT_gen_weights_pdf_hessians[36]);
			systematic_weights->Fill(PDFCT14n37    , nominal_syst_weight * NT_gen_weights_pdf_hessians[37]);
			systematic_weights->Fill(PDFCT14n38    , nominal_syst_weight * NT_gen_weights_pdf_hessians[38]);
			systematic_weights->Fill(PDFCT14n39    , nominal_syst_weight * NT_gen_weights_pdf_hessians[39]);
			systematic_weights->Fill(PDFCT14n40    , nominal_syst_weight * NT_gen_weights_pdf_hessians[40]);
			systematic_weights->Fill(PDFCT14n41    , nominal_syst_weight * NT_gen_weights_pdf_hessians[41]);
			systematic_weights->Fill(PDFCT14n42    , nominal_syst_weight * NT_gen_weights_pdf_hessians[42]);
			systematic_weights->Fill(PDFCT14n43    , nominal_syst_weight * NT_gen_weights_pdf_hessians[43]);
			systematic_weights->Fill(PDFCT14n44    , nominal_syst_weight * NT_gen_weights_pdf_hessians[44]);
			systematic_weights->Fill(PDFCT14n45    , nominal_syst_weight * NT_gen_weights_pdf_hessians[45]);
			systematic_weights->Fill(PDFCT14n46    , nominal_syst_weight * NT_gen_weights_pdf_hessians[46]);
			systematic_weights->Fill(PDFCT14n47    , nominal_syst_weight * NT_gen_weights_pdf_hessians[47]);
			systematic_weights->Fill(PDFCT14n48    , nominal_syst_weight * NT_gen_weights_pdf_hessians[48]);
			systematic_weights->Fill(PDFCT14n49    , nominal_syst_weight * NT_gen_weights_pdf_hessians[49]);
			systematic_weights->Fill(PDFCT14n50    , nominal_syst_weight * NT_gen_weights_pdf_hessians[50]);
			systematic_weights->Fill(PDFCT14n51    , nominal_syst_weight * NT_gen_weights_pdf_hessians[51]);
			systematic_weights->Fill(PDFCT14n52    , nominal_syst_weight * NT_gen_weights_pdf_hessians[52]);
			systematic_weights->Fill(PDFCT14n53    , nominal_syst_weight * NT_gen_weights_pdf_hessians[53]);
			systematic_weights->Fill(PDFCT14n54    , nominal_syst_weight * NT_gen_weights_pdf_hessians[54]);
			systematic_weights->Fill(PDFCT14n55    , nominal_syst_weight * NT_gen_weights_pdf_hessians[55]);
			systematic_weights->Fill(PDFCT14n56    , nominal_syst_weight * NT_gen_weights_pdf_hessians[56]);
			}
		}
	weight_counter->Fill(event_checkpoint, weight);
	weight_counter->Fill(++event_checkpoint, weight * NT_aMCatNLO_weight);
	// in old times (2015) there was a recommendation to actually use this as weight -- checking this

	if (isMC)
		{
		LogInfo ("Demo") << "gen sizes " << gen_leps.size() << ' ' << gen_taus.size() << ' ' << gen_tau3ch.size() << ' ' << gen_taulep.size() << ' ' << gen_w_prods.size() << ' ' << gen_b_prods.size();
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

	/* sting interface in 2016
	string filters_name = "";
	// 2016 data, Aug 2017 rereco
	if (is2016legacy)
		{
		filters_name = "RECO";
		}
	else
		{
		// 2016 data, Feb 2017 rereco
		filters_name = "PAT";
		}
	*/

	/* Token interface in 2016
	*/
	// get trigger results for all the tokens:
	//trigger result object for "HLT" triggers, for "RECO" and "PAT", these are 3 different sources of triggers
	edm::Handle<edm::TriggerResults> trigResults, trigResultsRECO, trigResultsPAT, trigResults_with_met_filters;
	//edm::InputTag * trigResultsTag; // the tag object, trigResults are extracted from the event via this tag

	iEvent.getByToken( trigResults_, trigResults );
	iEvent.getByToken( trigResultsRECO_, trigResultsRECO );
	if(!is2017data)
	{
		iEvent.getByToken( trigResultsPAT_, trigResultsPAT );
	}

	// 2016 data, Aug 2017 rereco
	if (is2016legacy || is2017data)
		{
		//is present only if PAT (and miniAOD) is not run simultaniously with RECO
		trigResults_with_met_filters = trigResultsRECO;
		}
	else
		{
		// 2016 data, Feb 2017 rereco, 2017legacy
		trigResults_with_met_filters = trigResultsPAT;
		}

	edm::TriggerResultsByName patFilters = iEvent.triggerResultsByName(*trigResults_with_met_filters);

	//if(!isMC && !metFilters.isValid()){metFilters = iEvent.triggerResultsByName("PAT");} //if not present, then it's part of RECO
	//if(!isMC && !metFilters.isValid()){       
	//	LogInfo("Demo") << "TriggerResultsByName for MET filters is not found in the process, as a consequence the MET filter is disabled for this event";
	//	return;
	//	}

	//if (!isMC && metFilters.isValid())
	// apparently MET POG suggests trying and looking at filters in MC
	// and they are there, in patFilters
	//if (!isMC && patFilters.isValid())
	if (!isMC)
		{
		// event is good if all filters ar true
		NT_filters_hbhe             = utils::passTriggerPatterns(patFilters, "Flag_HBHENoiseFilter*", "Flag_HBHENoiseIsoFilter*");
		NT_filters_ecalDeadCellTrig = utils::passTriggerPatterns(patFilters, "Flag_EcalDeadCellTriggerPrimitiveFilter*");
		NT_filters_good_vertices    = utils::passTriggerPatterns(patFilters, "Flag_goodVertices");
		NT_filters_eebad            = utils::passTriggerPatterns(patFilters, "Flag_eeBadScFilter");
		NT_filters_halo             = utils::passTriggerPatterns(patFilters, "Flag_globalTightHalo2016Filter");
		NT_filters_halo_super       = utils::passTriggerPatterns(patFilters, "Flag_globalSuperTightHalo2016Filter");

		// 2016 thing: bad muons
		if (is2016legacy)
			{
			NT_BadChargedCandidateFilter = utils::passTriggerPatterns(patFilters, "Flag_BadChargedCandidateFilter");
			NT_BadPFMuonFilter           = utils::passTriggerPatterns(patFilters, "Flag_BadPFMuonFilter");
			}
		else if(!is2017data) // TODO make bool for 2017 data
			{
			// TODO: check what this is? the 2016original?
			NT_filters_noBadMuons     = utils::passTriggerPatterns(patFilters, "Flag_noBadMuons"); // <---- the bad muons are done on the fly with cfg.py thingy
			NT_filters_duplicateMuons = utils::passTriggerPatterns(patFilters, "Flag_duplicateMuons");
			}
		//NT_BadChargedCandidateFilter = utils::passTriggerPatterns(patFilters, "Flag_BadChargedCandidateFilter");
		//NT_BadPFMuonFilter           = utils::passTriggerPatterns(patFilters, "Flag_BadPFMuonFilter");
		//NT_filters_noBadMuons        = utils::passTriggerPatterns(patFilters, "Flag_noBadMuons"); // <---- the bad muons are done on the fly with cfg.py thingy
		//NT_filters_duplicateMuons    = utils::passTriggerPatterns(patFilters, "Flag_duplicateMuons");

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
	bool low_pt_eTrigger = !withHLT;
	bool low_pt_muTrigger1 = !withHLT, low_pt_muTrigger2 = !withHLT, low_pt_muTrigger = !withHLT;
	bool jetsHLT140 = !withHLT, jetsHLT400 = !withHLT, jetsHLT = !withHLT;

	bool NT_HLT_el_low_pt32         = !withHLT;
	bool NT_HLT_el_low_pt28_150HT   = !withHLT;
	bool NT_HLT_el_low_pt30_35PFJet = !withHLT;
	bool NT_HLT_elmu_1 = !withHLT;
	bool NT_HLT_elmu_2 = !withHLT;
	bool NT_HLT_elmu_3 = !withHLT;
	bool NT_HLT_elmu_4 = !withHLT;
	bool NT_HLT_elel_1 = !withHLT;
	bool NT_HLT_elel_2 = !withHLT;

	// TriggerNames for TriggerObjects --------------------

	string matched_elTriggerName("");
	string matched_muTriggerName1("");
	string matched_muTriggerName2("");
	//edm::TriggerResultsByName tr = ev.triggerResultsByName ("HLT");

	// HLT matching
	// objects of our triggers
	vector<pat::TriggerObjectStandAlone> el_trig_objs;
	vector<pat::TriggerObjectStandAlone> mu_trig_objs, mu_trig_objs1, mu_trig_objs2;

	if (withHLT)
		{
		//edm::TriggerResultsByName tr = iEvent.triggerResultsByName (HLT_source);
		edm::TriggerResultsByName tr = iEvent.triggerResultsByName(*trigResults);

		//if (!tr.isValid ()){
			// HLT2 was a quirk of Spring16 MC campaigns (noHLT/reHLT/withHLT thing)
			// need to compare 2016-2015 in tau SV
			// using HLT2 as backup (DY50 with HLT is present only in reHLT campaign, wich has this HLT2 path)
			//tr = iEvent.triggerResultsByName ("HLT2");
			// it crashes weirdly with "no consumes" complaints
			//}

		if (!tr.isValid ()){
			//cout << HLT_source << " is NOT valid!" << endl;
			//return;
			throw cms::Exception("NoTrigger") << "!tr.isValid()" << "\n";
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
			eTrigger   = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_elTriggerName,  elHLT_MC)  : utils::passTriggerPatternsAndGetName(tr, matched_elTriggerName,  elHLT_Data));
			muTrigger1 = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName1, muHLT_MC1) : utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName1, muHLT_Data1));
			muTrigger2 = (isMC ? utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName2, muHLT_MC2) : utils::passTriggerPatternsAndGetName(tr, matched_muTriggerName2, muHLT_Data2));
			muTrigger = muTrigger1 || muTrigger2;

			low_pt_eTrigger   = (isMC ?
					utils::passTriggerPatterns(tr, low_pt_elHLT_MC)  :
					utils::passTriggerPatterns(tr, low_pt_elHLT_Data));
			low_pt_muTrigger1 = (isMC ?
					utils::passTriggerPatterns(tr, low_pt_muHLT_MC1) :
					utils::passTriggerPatterns(tr, low_pt_muHLT_Data1));
			low_pt_muTrigger2 = (isMC ?
					utils::passTriggerPatterns(tr, low_pt_muHLT_MC2) :
					utils::passTriggerPatterns(tr, low_pt_muHLT_Data2));
			low_pt_muTrigger = low_pt_muTrigger1 || low_pt_muTrigger2;

			NT_HLT_el_low_pt32         = utils::passTriggerPatterns(tr, low_pt32_elHLT);
			NT_HLT_el_low_pt28_150HT   = utils::passTriggerPatterns(tr, low_pt28_150HT_elHLT);
			NT_HLT_el_low_pt30_35PFJet = utils::passTriggerPatterns(tr, low_pt30_35PFJet_elHLT);
			NT_HLT_elmu_1              = utils::passTriggerPatterns(tr, elmuHLT_1);
			NT_HLT_elmu_2              = utils::passTriggerPatterns(tr, elmuHLT_2);
			NT_HLT_elmu_3              = utils::passTriggerPatterns(tr, elmuHLT_3);
			NT_HLT_elmu_4              = utils::passTriggerPatterns(tr, elmuHLT_4);
			NT_HLT_elel_1              = utils::passTriggerPatterns(tr, elelHLT_1);
			NT_HLT_elel_2              = utils::passTriggerPatterns(tr, elelHLT_2);

			NT_HLT_eltau              = utils::passTriggerPatterns(tr, eltauHLT);
			NT_HLT_mutau1             = utils::passTriggerPatterns(tr, mutauHLT1);
			NT_HLT_mutau2             = utils::passTriggerPatterns(tr, mutauHLT2);

			jetsHLT140 = utils::passTriggerPatterns(tr, "HLT_PFJet140_v*");
			jetsHLT400 = utils::passTriggerPatterns(tr, "HLT_PFJet400_v*");
			}

		if (record_jets)
			{
			jetsHLT = jetsHLT140 || jetsHLT400;
			}

		// get trigger objects for HLT matching

		// names for trigger bits
		//edm::EDGetTokenT<edm::TriggerResults> trigResults_ = consumes<edm::TriggerResults>(trigResultsTag);
		//ev.getByLabel(*trigResultsTag, trigResults);
		const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);

		//fwlite::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjectsHandle;
		edm::Handle<vector<pat::TriggerObjectStandAlone>> triggerObjectsHandle;
		//triggerObjectsHandle.getByLabel(ev, "selectedPatTrigger");
		iEvent.getByToken(triggerObjects_, triggerObjectsHandle);
		if (!triggerObjectsHandle.isValid())
			{
			//LogInfo("Demo") << "!triggerObjectsHandle.isValid()";
			//return;
			throw cms::Exception("NoTriggerObjects") << "!triggerObjectsHandle.isValid() -- no products under the given name " << triggerObjects_InputTag.label() << "\n";
			}
		LogInfo ("Demo") << "got trigger objects";
		vector<pat::TriggerObjectStandAlone> trig_objs = *triggerObjectsHandle;

		if (eTrigger)
			{
			Processing_selectHLTobjects(trig_objs, trigNames, el_trig_objs, matched_elTriggerName);
			}
		if (muTrigger)
			{
			if (muTrigger1)
				{ Processing_selectHLTobjects(trig_objs, trigNames, mu_trig_objs1, matched_muTriggerName1); }
			if (muTrigger2)
				{ Processing_selectHLTobjects(trig_objs, trigNames, mu_trig_objs2, matched_muTriggerName2); }
			// vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
			mu_trig_objs.insert(mu_trig_objs.end(), mu_trig_objs1.begin(), mu_trig_objs1.end());
			mu_trig_objs.insert(mu_trig_objs.end(), mu_trig_objs2.begin(), mu_trig_objs2.end());
			}

		}

	event_checkpoint++;

	event_checkpoint++;
	if (eTrigger)
		event_counter ->Fill(event_checkpoint);

	event_checkpoint++;
	if (muTrigger)
		event_counter ->Fill(event_checkpoint);

	event_checkpoint++;

	bool any_trigger = eTrigger || muTrigger || low_pt_eTrigger || low_pt_muTrigger || lepMonitorTrigger || jetsHLT ||
			NT_HLT_el_low_pt32         ||
			NT_HLT_el_low_pt28_150HT   ||
			NT_HLT_el_low_pt30_35PFJet ||
			NT_HLT_elmu_1              ||
			NT_HLT_elmu_2              ||
			NT_HLT_elmu_3              ||
			NT_HLT_elmu_4              ||
			NT_HLT_elel_1              ||
			NT_HLT_elel_2              ||
			NT_HLT_eltau               ||
			NT_HLT_mutau1              ||
			NT_HLT_mutau2              ||
			false;

	if (!record_all && !record_signal && !any_trigger) return; // orthogonalization is done afterwards
	event_counter ->Fill(event_checkpoint++);
	weight_counter->Fill(event_checkpoint, weight);

	NT_HLT_el = eTrigger;
	NT_HLT_mu = muTrigger;
	NT_HLT_el_low_pt = low_pt_eTrigger;
	NT_HLT_mu_low_pt = low_pt_muTrigger;
	NT_HLT_lepMonitor = lepMonitorTrigger;
	NT_HLT_jets140 = jetsHLT140;
	NT_HLT_jets400 = jetsHLT400;

	LogInfo ("Demo") << "passed HLT " << eTrigger << ' ' << muTrigger << '(' << muTrigger1 << ',' << muTrigger2 << ')' << ';' << matched_elTriggerName << ' ' << matched_muTriggerName1 << ',' << matched_muTriggerName2 << ' ' << el_trig_objs.size() << ' ' << mu_trig_objs.size() << '(' << mu_trig_objs1.size() << ',' << mu_trig_objs2.size() << ')';
	//LogInfo ("Demo") << "our trigger objects: " << el_trig_objs.size() << ' ' << mu_trig_objs.size();
	LogInfo ("Demo") << "passed low pt HLT " << low_pt_eTrigger << ' ' << low_pt_muTrigger ;
	LogInfo ("Demo") << "passed  other HLT "
			<< " " << NT_HLT_el_low_pt32
			<< " " << NT_HLT_el_low_pt28_150HT
			<< " " << NT_HLT_el_low_pt30_35PFJet
			<< " " << NT_HLT_elmu_1
			<< " " << NT_HLT_elmu_2
			<< " " << NT_HLT_elmu_3
			<< " " << NT_HLT_elmu_4
			<< " " << NT_HLT_elel_1
			<< " " << NT_HLT_elel_2 ;

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

	LogInfo ("Demo") << "selected vertex";
	//weight = 1; // reset weights?

	// MUONS
	LorentzVector muDiff(0., 0., 0., 0.);
	unsigned int nVetoMu_Iso = 0, nVetoMu_all = 0;
	//pat::MuonCollection selIDMuons, selMuons;
	pat::MuonCollection selMuons, selMuons_allIso;
	//processMuons_ID_ISO_Kinematics(muons, goodPV, weight, patUtils::llvvMuonId::StdTight, patUtils::llvvMuonId::StdLoose, patUtils::llvvMuonIso::Tight, patUtils::llvvMuonIso::Loose,               
	processMuons_ID_ISO_Kinematics(muons, goodPV, weight, patUtils::llvvMuonId::StdTight, patUtils::llvvMuonId::StdLoose, patUtils::llvvMuonIso::StdTight, patUtils::llvvMuonIso::StdLoose,
		mu_kino_cuts_pt, mu_kino_cuts_eta, mu_veto_kino_cuts_pt, mu_veto_kino_cuts_eta,
		selMuons, selMuons_allIso,
		muDiff, nVetoMu_Iso, nVetoMu_all, false, false);

	//nVetoMu += processMuons_MatchHLT(selIDMuons, mu_trig_objs, 0.4, selMuons);
	LogInfo ("Demo") << "passed muons" << selMuons.size() << ' ' << selMuons_allIso.size();

	// ELECTRONS
	//pat::ElectronCollection selIDElectrons, selElectrons;
	pat::ElectronCollection selElectrons, selElectrons_allIso;
	unsigned int nVetoE_Iso = 0, nVetoE_all = 0;
	LorentzVector elDiff(0., 0., 0., 0.);
	//processElectrons_ID_ISO_Kinematics(electrons, goodPV, NT_fixedGridRhoFastjetAll, weight, patUtils::llvvElecId::Tight, patUtils::llvvElecId::Loose, patUtils::llvvElecIso::Tight, patUtils::llvvElecIso::Loose,
	processElectrons_ID_ISO_Kinematics(electrons, goodPV, NT_fixedGridRhoFastjetAll, weight, patUtils::llvvElecId::StdTight, patUtils::llvvElecId::StdLoose, patUtils::llvvElecIso::StdTight, patUtils::llvvElecIso::StdLoose,
		el_kino_cuts_pt, el_kino_cuts_eta, el_veto_kino_cuts_pt, el_veto_kino_cuts_eta,
		selElectrons, selElectrons_allIso,
		elDiff, nVetoE_Iso, nVetoE_all, false, false);

	LogInfo ("Demo") << "passed electrons" << selElectrons.size() << ' ' << selElectrons_allIso.size();
	//nVetoE += processElectrons_MatchHLT(selIDElectrons, el_trig_objs, 0.4, selElectrons);

	std::vector<patUtils::GenericLepton> selLeptons;
	for(size_t l=0; l<selElectrons.size(); ++l) selLeptons.push_back(patUtils::GenericLepton (selElectrons[l] ));
	for(size_t l=0; l<selMuons.size(); ++l)     selLeptons.push_back(patUtils::GenericLepton (selMuons[l]     ));
	//std::sort(selLeptons.begin(), selLeptons.end(), utils::sort_CandidatesByPt);

	std::vector<patUtils::GenericLepton> selLeptons_allIso;
	for(size_t l=0; l<selElectrons_allIso.size(); ++l) selLeptons_allIso.push_back(patUtils::GenericLepton (selElectrons_allIso[l] ));
	for(size_t l=0; l<selMuons_allIso.size(); ++l)     selLeptons_allIso.push_back(patUtils::GenericLepton (selMuons_allIso[l]     ));

	//LogInfo ("Demo") << "selected leptons: " << '(' << selIDElectrons.size() << ',' << selIDMuons.size() << ')' <<  selLeptons.size() << ' ' << nVetoE << ',' << nVetoMu;
	LogInfo ("Demo") << "selected leptons: " << '(' << selElectrons.size() << ',' << selMuons.size() << ')' <<  selLeptons.size() << ' ' << nVetoE_Iso << ',' << nVetoMu_Iso;

	/*
	 * TAUS preliminary
	 */
	//string tau_Loose_ID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
	//string tau_VLoose_ID ("byVLooseIsolationMVArun2v1DBoldDMwLT");
	//string tau_Loose_ID  ("byLooseIsolationMVArun2v1DBoldDMwLT");
	//string tau_Medium_ID ("byMediumIsolationMVArun2v1DBoldDMwLT");
	//string tau_Tight_ID  ("byTightIsolationMVArun2v1DBoldDMwLT");
	//string tau_VTight_ID ("byVTightIsolationMVArun2v1DBoldDMwLT");
	string tau_decayMode       ("decayModeFinding");
	//string tau_decayMode       ("decayModeFindingOldDMs");
	string tau_againstMuon     ("againstMuonTight3");
	string tau_againstElectron ("againstElectronTightMVA6");

	//("byIsolationMVArun2017v2DBoldDMwLTraw2017");

	pat::TauCollection IDtaus, selTaus;
	processTaus_ID    (taus,   weight, tau_decayMode, tau_againstMuon, tau_againstElectron, IDtaus, false, false);
	//processTaus_ID_ISO    (taus,   weight, tau_decayMode, tau_VLoose_ID, tau_againstMuon, tau_againstElectron, IDtaus, false, false);
	processTaus_Kinematics(IDtaus, weight, tau_kino_cuts_pt, tau_kino_cuts_eta, selTaus,      false, false);

	//pat::TauCollection selTausNoLep;
	//crossClean_in_dR(selTaus,       selLeptons, 0.4, selTausNoLep,        weight, string("selTausNoLep"),        false, false);
	//NT_ntaus = 0;
	NT_ntaus = selTaus.size(); // all tau before MVA anti-jet iso
	NT_ntaus_idVL = 0;
	NT_ntaus_idM  = 0;
	NT_ntaus_idT  = 0;

	// and these are the NT output taus
	// taus are sorted py pt
	//std::sort (selTaus.begin(),  selTaus.end(),  utils::sort_CandidatesByPt);

	// sort taus by IDlev and pT
	// first add IDlev into the tau

	for(size_t i=0; i<selTaus.size(); ++i)
		{
		pat::Tau& tau = selTaus[i];

		Int_t IDlev = 0;
		if      (tau.tauID(tau_VTight_ID) >= 0.5) IDlev = 5;
		else if (tau.tauID(tau_Tight_ID)  >= 0.5) IDlev = 4;
		else if (tau.tauID(tau_Medium_ID) >= 0.5) IDlev = 3;
		else if (tau.tauID(tau_Loose_ID)  >= 0.5) IDlev = 2;
		else if (tau.tauID(tau_VLoose_ID) >= 0.5) IDlev = 1;
		tau.addUserInt("IDlev", IDlev);
		}
	std::sort (selTaus.begin(),  selTaus.end(),  sort_TausByIDByPt);

	//bool clean_lep_conditions = nVetoE==0 && nVetoMu==0 && nGoodPV != 0; // veto on std iso veto leptons
	//bool clean_lep_conditions = nVetoE_all==0 && nVetoMu_all==0 && nGoodPV != 0; // veto on all iso veto leptons
	bool clean_lep_conditions = nGoodPV != 0; // just good PV, the loosest req,save bit if no veto leps
	//if (!(clean_lep_conditions && ((selLeptons.size() > 0 && selLeptons.size() < 3 && nVetoE_Iso == 0 && nVetoMu_Iso == 0) || (selLeptons_allIso.size() == 1 && nVetoE_all == 0 && nVetoMu_all == 0)) )) return;
	if (!record_all && !record_signal && !clean_lep_conditions) return;
	if (!record_all && !record_signal && !((selLeptons.size() > 0 && selLeptons.size() < 3) || selLeptons_allIso.size() == 1 || selTaus.size() > 0)) return;
	// exit now to reduce computation -- all record schemes have this requirement

	event_counter ->Fill(event_checkpoint++);
	weight_counter->Fill(event_checkpoint, weight);

	LogInfo ("Demo") << "passed lepton conditions ";


	event_checkpoint++;

	event_checkpoint++;
	if (selElectrons.size() == 1)
		event_counter ->Fill(event_checkpoint);

	event_checkpoint++;
	if (selMuons.size() == 1)
		event_counter ->Fill(event_checkpoint);

	event_checkpoint++;
	if (selLeptons.size() == 1)
		event_counter ->Fill(event_checkpoint);


	bool no_veto = nVetoE_Iso == 0 && nVetoMu_Iso == 0;
	event_checkpoint++;
	if (selElectrons.size() == 1 && no_veto)
		event_counter ->Fill(event_checkpoint);

	event_checkpoint++;
	if (selMuons.size() == 1 && no_veto)
		event_counter ->Fill(event_checkpoint);

	event_checkpoint++;
	if (selLeptons.size() == 1 && no_veto)
		event_counter ->Fill(event_checkpoint);

	event_checkpoint++;

	for(size_t l=0; l<selMuons.size(); ++l)
		{
		LogInfo ("Demo") << "saving mu " << l;

		NT_lep_p4.push_back(selMuons[l].p4());
		NT_lep_id.push_back(selMuons[l].pdgId());
		// mu_trig_objs or el_trig_objs
		// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#2016_Data
		// > Offline muons selected for analysis should always be associated to the HLT objects that triggered the event. Recommended matching criterion: ΔR(HLT object, offline muon) < 0.1.
		struct dR_matching match_to_HLTs = dR_match_to_HLTs(selMuons[l], mu_trig_objs, 0.1);
		NT_lep_matched_HLT    .push_back(match_to_HLTs.matched);
		NT_lep_matched_HLT_dR .push_back(match_to_HLTs.dR);
		NT_lep_dz  .push_back(selMuons[l].muonBestTrack()->dz (goodPV.position()));
		NT_lep_dxy .push_back(selMuons[l].muonBestTrack()->dxy (goodPV.position()));
		NT_lep_dB.push_back(selMuons[l].dB());
		NT_lep_N_trackerLayersWithMeasurement .push_back(selMuons[l].innerTrack()->hitPattern().trackerLayersWithMeasurement());

		float rel_iso = relIso(selMuons[l], NT_fixedGridRhoFastjetAll);
		NT_lep_relIso.push_back(rel_iso);
		// using old procedures for now
		NT_lep_passIso.push_back(patUtils::passIso(selMuons[l], patUtils::llvvMuonIso::Tight, patUtils::CutVersion::Moriond17Cut));
		if (isMC)
			{
			struct gen_matching match = match_to_gen(selMuons[l].p4(), gen_leps, gen_taus, gen_tau3ch, gen_taulep, gen_w_prods, gen_b_prods, gid_leps, gid_taus, gid_tau3ch, gid_taulep, gid_w_prods, gid_b_prods);
			NT_lep_matching_gen   .push_back(match.closest);
			NT_lep_matching_gen_dR.push_back(match.dR);
			if (isDY || isTT)
				{
				struct gen_matching_collection match2 = match_to_gen_collection(selMuons[l].p4(), NT_gen_tt_tau_vis_p4);
				NT_lep_matching_gen_collection   .push_back(match2.index);
				NT_lep_matching_gen_collection_dR.push_back(match2.dR);
				}
			}
		//double dataSF = rc.kScaleDT(selMuons[l].charge(), selMuons[l].pt(), selMuons[l].eta(), selMuons[l].phi(), 0, 0);
		//NT_lep_correction.push_back(dataSF);

		NT_lep_energy_ScaleUp      .push_back(0.);
		NT_lep_energy_ScaleDown    .push_back(0.);
		NT_lep_energy_SmearUp      .push_back(0.);
		NT_lep_energy_SmearDown    .push_back(0.);
		NT_lep_energy_2016legacy_ScaleEtUp    .push_back(0.);
		NT_lep_energy_2016legacy_ScaleEtDown  .push_back(0.);
        	NT_lep_energy_parameter1_r9.push_back(0.);
		}

	LogInfo ("Demo") << "saved muons";

	for(size_t l=0; l<selElectrons.size(); ++l)
		{
		//NT_lep_p4.push_back(selElectrons[l].p4());
		// legacy 2016 nominal corrected electron energy:
		NT_lep_p4.push_back(selElectrons[l].p4() * selElectrons[l].userFloat("ecalTrkEnergyPostCorr") / selElectrons[l].energy());
		// systematic variations of the correction:
		// Scale and Smear -- saving them for both data and MC, but apply only in MC
		NT_lep_energy_ScaleUp   .push_back(selElectrons[l].userFloat("energyScaleUp"));
		NT_lep_energy_ScaleDown .push_back(selElectrons[l].userFloat("energyScaleDown"));
		NT_lep_energy_SmearUp   .push_back(selElectrons[l].userFloat("energySigmaUp"));
		NT_lep_energy_SmearDown .push_back(selElectrons[l].userFloat("energySigmaDown"));
		// legacy 2016 feature: energy invervion at 45GeV, correction with this additional variation
		if (is2016legacy) {
			NT_lep_energy_2016legacy_ScaleEtUp   .push_back(selElectrons[l].userFloat("energyScaleEtUp"));
			NT_lep_energy_2016legacy_ScaleEtDown .push_back(selElectrons[l].userFloat("energyScaleEtDown"));
			}
		else {
			NT_lep_energy_2016legacy_ScaleEtUp   .push_back(0.);
			NT_lep_energy_2016legacy_ScaleEtDown .push_back(0.);
			}

		NT_lep_id.push_back(selElectrons[l].pdgId());
		// mu_trig_objs or el_trig_objs
		// Propagating 0.1 -> 0.2 trigger match to el untill found the recommended value
		struct dR_matching match_to_HLTs = dR_match_to_HLTs(selElectrons[l], el_trig_objs, 0.2);
		NT_lep_matched_HLT    .push_back(match_to_HLTs.matched);
		NT_lep_matched_HLT_dR .push_back(match_to_HLTs.dR);
		NT_lep_dz  .push_back(selElectrons[l].gsfTrack()->dz (goodPV.position()));
		NT_lep_dxy .push_back(selElectrons[l].gsfTrack()->dxy (goodPV.position()));
		NT_lep_dB.push_back(selElectrons[l].dB());
		NT_lep_N_trackerLayersWithMeasurement .push_back(selElectrons[l].gsfTrack()->hitPattern().trackerLayersWithMeasurement());

		float rel_iso = relIso(selElectrons[l], NT_fixedGridRhoFastjetAll);
		NT_lep_relIso.push_back(rel_iso);
		//bool passIso = patUtils::passIso(selMuons[l], el_ISO, patUtils::CutVersion::Moriond17Cut, rho);
		NT_lep_passIso.push_back(patUtils::passIso(selElectrons[l], patUtils::llvvElecIso::Tight, patUtils::CutVersion::Moriond17Cut, NT_fixedGridRhoFastjetAll));
		if (isMC)
			{
			struct gen_matching match = match_to_gen(selElectrons[l].p4(), gen_leps, gen_taus, gen_tau3ch, gen_taulep, gen_w_prods, gen_b_prods, gid_leps, gid_taus, gid_tau3ch, gid_taulep, gid_w_prods, gid_b_prods);
			NT_lep_matching_gen   .push_back(match.closest);
			NT_lep_matching_gen_dR.push_back(match.dR);

			if (isDY || isTT)
				{
				struct gen_matching_collection match2 = match_to_gen_collection(selElectrons[l].p4(), NT_gen_tt_tau_vis_p4);
				NT_lep_matching_gen_collection   .push_back(match2.index);
				NT_lep_matching_gen_collection_dR.push_back(match2.dR);
				}
			}

        	NT_lep_energy_parameter1_r9.push_back(selElectrons[l].r9());
		}

	NT_nleps = selLeptons.size();

	LogInfo ("Demo") << "saved electrons";


	// for control
	//NT_nleps_veto_el_isoimp = nVetoE_IsoImp;
	NT_nleps_veto_el_iso = nVetoE_Iso;
	NT_nleps_veto_el_all = nVetoE_all;
	NT_nleps_veto_mu_iso = nVetoMu_Iso;
	NT_nleps_veto_mu_all = nVetoMu_all;
	//NT_no_std_veto_leps  = nVetoE_IsoImp == 0 && nVetoMu_Iso == 0;
	NT_no_iso_veto_leps  = nVetoE_Iso == 0 && nVetoMu_Iso == 0;

	NT_leps_ID = 1;
	for (unsigned int i = 0; i<selLeptons.size(); i++)
		{
		NT_leps_ID *= selLeptons[i].pdgId();
		}

	// all iso leptons for QCD study
	/* this selection won't reproduce the veto on relIso
	 * -- TODO: not sure if no-veto for QCD is fine
	 */
	NT_leps_ID_allIso = 1;

	for(size_t l=0; l<selElectrons_allIso.size(); ++l)
		{
		NT_leps_ID_allIso *= selElectrons_allIso[l].pdgId();

		NT_lep_alliso_p4.push_back(selElectrons_allIso[l].p4());
		NT_lep_alliso_id.push_back(selElectrons_allIso[l].pdgId());
		float rel_iso = relIso(selElectrons_allIso[l], NT_fixedGridRhoFastjetAll);
		NT_lep_alliso_relIso.push_back(rel_iso);
		struct dR_matching match_to_HLTs = dR_match_to_HLTs(selElectrons_allIso[l], el_trig_objs, 0.2);
		//NT_lep_alliso_matched_HLT.push_back(processElectron_matchesHLTs(selElectrons_allIso[l], el_trig_objs, 0.2));
		NT_lep_alliso_matched_HLT   .push_back(match_to_HLTs.matched);
		NT_lep_alliso_matched_HLT_dR.push_back(match_to_HLTs.dR);
		if (isMC)
			{
			struct gen_matching match = match_to_gen(selElectrons_allIso[l].p4(), gen_leps, gen_taus, gen_tau3ch, gen_taulep, gen_w_prods, gen_b_prods, gid_leps, gid_taus, gid_tau3ch, gid_taulep, gid_w_prods, gid_b_prods);
			NT_lep_alliso_matching_gen   .push_back(match.closest);
			NT_lep_alliso_matching_gen_dR.push_back(match.dR);
			}
		}

	for(size_t l=0; l<selMuons_allIso.size(); ++l)
		{
		NT_leps_ID_allIso *= selMuons_allIso[l].pdgId();

		NT_lep_alliso_p4.push_back(selMuons_allIso[l].p4());
		NT_lep_alliso_id.push_back(selMuons_allIso[l].pdgId());
		float rel_iso = relIso(selMuons_allIso[l], NT_fixedGridRhoFastjetAll);
		NT_lep_alliso_relIso.push_back(rel_iso);
		struct dR_matching match_to_HLTs = dR_match_to_HLTs(selMuons_allIso[l], mu_trig_objs, 0.1);
		//NT_lep_alliso_matched_HLT.push_back(processMuon_matchesHLTs(selMuons_allIso[l], mu_trig_objs, 0.1));
		NT_lep_alliso_matched_HLT   .push_back(match_to_HLTs.matched);
		NT_lep_alliso_matched_HLT_dR.push_back(match_to_HLTs.dR);
		if (isMC)
			{
			struct gen_matching match = match_to_gen(selMuons_allIso[l].p4(), gen_leps, gen_taus, gen_tau3ch, gen_taulep, gen_w_prods, gen_b_prods, gid_leps, gid_taus, gid_tau3ch, gid_taulep, gid_w_prods, gid_b_prods);
			NT_lep_alliso_matching_gen   .push_back(match.closest);
			NT_lep_alliso_matching_gen_dR.push_back(match.dR);
			}
		}


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

		// https://cmssdt.cern.ch/lxr/source/DataFormats/PatCandidates/interface/MET.h?%21v=CMSSW_9_4_0
		// systematic shifts
		/*
		 *       enum METUncertainty {
		 *     JetResUp=0, JetResDown=1, JetEnUp=2, JetEnDown=3,
		 *         MuonEnUp=4, MuonEnDown=5, ElectronEnUp=6, ElectronEnDown=7,
		 *     TauEnUp=8, TauEnDown=9, UnclusteredEnUp=10, UnclusteredEnDown=11,
		 *     PhotonEnUp=12, PhotonEnDown=13, NoShift=14, METUncertaintySize=15,
		 *     JetResUpSmear=16, JetResDownSmear=17, METFullUncertaintySize=18
		 *       };
		 */

		NT_met_init_shift_UnclusteredEnUp   = MET.shiftedP4(pat::MET::UnclusteredEnUp);
		NT_met_init_shift_UnclusteredEnDown = MET.shiftedP4(pat::MET::UnclusteredEnDown);
		NT_met_init_shift_JetEnUp           = MET.shiftedP4(pat::MET::JetEnUp);
		NT_met_init_shift_JetEnDown         = MET.shiftedP4(pat::MET::JetEnDown);
		NT_met_init_shift_JetResUp          = MET.shiftedP4(pat::MET::JetResUp);
		NT_met_init_shift_JetResDown        = MET.shiftedP4(pat::MET::JetResDown);
		NT_met_init_shift_MuonEnUp          = MET.shiftedP4(pat::MET::MuonEnUp);
		NT_met_init_shift_MuonEnDown        = MET.shiftedP4(pat::MET::MuonEnDown);
		NT_met_init_shift_ElectronEnUp      = MET.shiftedP4(pat::MET::ElectronEnUp);
		NT_met_init_shift_ElectronEnDown    = MET.shiftedP4(pat::MET::ElectronEnDown);
		NT_met_init_shift_TauEnUp           = MET.shiftedP4(pat::MET::TauEnUp);
		NT_met_init_shift_TauEnDown         = MET.shiftedP4(pat::MET::TauEnDown);

		/*
		NT_met_init_shift_JetResUpSmear     = MET.shiftedP4(pat::MET::JetResUpSmear);
		NT_met_init_shift_JetResDownSmear   = MET.shiftedP4(pat::MET::JetResDownSmear);
		*/

		LogInfo ("Demo") << "met slimmed";
		}
	// LorentzVector met = mets_slimmedMETs[0].p4 ();

	// these are valid only for data
	if(metsHandle_slimmedMETsMuEGClean.isValid() )
		{
		const pat::MET& MET2 = metsHandle_slimmedMETsMuEGClean->front();
		NT_met_slimmedMETsMuEGClean = MET2.p4();
		LogInfo ("Demo") << "met MuEgClean";
		}

	// also for control let's get uncorrected met and compare the two:
	if (!isMC) // sadly this exists only in latest ReReco data made with 8.0.26 CMSSW, not in Summer16 MC
		{
		pat::METCollection mets_uncorrected;
		edm::Handle<pat::METCollection> mets_uncorrectedHandle;
		//mets_uncorrectedHandle.getByLabel(ev, "slimmedMETsUncorrected");
		iEvent.getByToken( mets_uncorrected_, mets_uncorrectedHandle);
		if(mets_uncorrectedHandle.isValid() )
			{
			mets_uncorrected = *mets_uncorrectedHandle;
			pat::MET met_uncorrected = mets_uncorrected[0];
			NT_met_uncorrected = met_uncorrected.p4();
			LogInfo ("Demo") << "met uncor";
			}
		}


	LogInfo ("Demo") << "passed mets";

	// JETS
	/* jets additionally need initialization of:
	 * genJets
	 * jesCor, totalJESUnc,
	 * pass jet ID, PU jet ID (with/without PU),
	 * systematic variation (NOMINAL, the variation factors are saved per jet for offline)
	 * jet resolution in pt, eta (?)
	 * kinematic cuts
	 */

	// PAT jets are opened above
	// to have access to their partonFlavour/hadronFlavour
	// at processing MC weights

	// get genJets from the event
	std::vector<reco::GenJet> genJets;
	edm::Handle<std::vector<reco::GenJet>> genJetsHandle;
	//genJetsHandle.getByLabel(ev, "slimmedGenJets");
	iEvent.getByToken( genJets_, genJetsHandle); // twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#GenJets
	if (genJetsHandle.isValid() ) genJets = *genJetsHandle;

	LorentzVector full_jet_corr(0., 0., 0., 0.);
	//pat::JetCollection IDjets;
	pat::JetCollection selJets;
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
	NT_nbjets_noVLooseTau = 0;
	NT_nbjets_noMediumTau = 0;
	NT_nbjets_noTightTau  = 0;

	NT_nMbjets = 0;
	NT_nMbjets_noVLooseTau = 0;
	NT_nMbjets_noMediumTau = 0;
	NT_nMbjets_noTightTau  = 0;

	NT_nallbjets = 0;
	NT_njets = 0;
	NT_nalljets = 0;
	string btagger_label("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	string btagger_2016legacy_1("pfDeepCSVJetTags:probb");
	string btagger_2016legacy_2("pfDeepCSVJetTags:probbb");

	for (unsigned int ijet=0; ijet<jets.size(); ijet++)
		{
		pat::Jet& jet = jets[ijet];
		// the only requirements are pt 111 and no-leptons 222 in dR and Loose Jet ID 333
		// (it is not a subject of systematic variations - right?)
		// maybe later remove the Loose ID requirement in ntuple for further studies
		// but I've never used it yet -- thus applying it for speed & size

		if (jet.pt() < jet_kino_cuts_pt) continue; // 111
		// all eta pass -- forward jets too, usefull for WJets control region

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

		/*
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
		*/

		// for 2017 94X see:
		// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
		if (abseta <= 2.7)
			{
			// no Loose ID anymore!
			//looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abseta<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abseta>2.4);
			// efficiency > 99% everywhere
			// no CEMF<0.99
			tightJetID        =  (NHF<0.90 && NEMF<0.90 && NumConst>1) &&             ((abseta<=2.4 && CHF>0 && CHM>0) || abseta>2.4);
			tightLepVetoJetID = ((NHF<0.90 && NEMF<0.90 && NumConst>1  && MUF<0.8) && ((abseta<=2.4 && CHF>0 && CHM>0 && CEMF<0.80) || abseta>2.4));
			}
		else if (abseta <= 3.0)
			{
			//looseJetID = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2);
			tightJetID = (NEMF>0.02 && NEMF<0.99 && NumNeutralParticles>2);
			}
		else
			{
			//looseJetID = (NEMF<0.90 && NumNeutralParticles>10);
			tightJetID = (NEMF<0.90 && NHF>0.02 && NumNeutralParticles>10);
			}

		// no Loose ID anymore!
		// efficiency > 99% everywhere
		looseJetID = tightJetID;
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
		NT_jet_b_discr_original2016.push_back(b_discriminator);
		//if (b_discriminator > btag_threshold) NT_nallbjets += 1;

		// legacy 2016 rereco
		float b_discriminator_legacy2016_1 = jet.bDiscriminator(btagger_2016legacy_1);
		float b_discriminator_legacy2016_2 = jet.bDiscriminator(btagger_2016legacy_2);
		NT_jet_b_DeepCSV_legacy2016_1.push_back(b_discriminator_legacy2016_1);
		NT_jet_b_DeepCSV_legacy2016_2.push_back(b_discriminator_legacy2016_2);
		NT_jet_b_discr.push_back(b_discriminator_legacy2016_1 + b_discriminator_legacy2016_2);
		if ((b_discriminator_legacy2016_1 + b_discriminator_legacy2016_2) > btag_threshold) NT_nallbjets += 1;

		NT_jet_hadronFlavour.push_back(jet.hadronFlavour());
		NT_jet_partonFlavour.push_back(jet.partonFlavour());

		int jetid = 0;
		if (looseJetID)
			{
			++jetid;
			if (abseta < jet_kino_cuts_eta)
				{
				NT_njets += 1;
				if (b_discriminator > btag_threshold)
					{
					NT_nbjets += 1;
					if (b_discriminator > btag_threshold_Medium) {NT_nMbjets +=1;}

					// match to taus
					bool matched_VLoose = false, matched_Medium = false, matched_Tight = false;
					for(size_t i=0; i<selTaus.size(); ++i)
						{
						pat::Tau& tau = selTaus[i];
						Float_t tau_dR = reco::deltaR(jet, tau);
						if (tau_dR>0.4) continue;

						// a matched tau
						if (tau.tauID(tau_Tight_ID) >= 0.5)
							{
							matched_VLoose = true;
							matched_Medium = true;
							matched_Tight  = true;
							}
						else if (tau.tauID(tau_Medium_ID) >= 0.5)
							{
							matched_VLoose = true;
							matched_Medium = true;
							}
						else if (tau.tauID(tau_VLoose_ID) >= 0.5)
							{
							matched_VLoose = true;
							}
						}

					if (!matched_VLoose)
						{
						NT_nbjets_noVLooseTau += 1;
						NT_nbjets_noMediumTau += 1;
						NT_nbjets_noTightTau  += 1;
						if (b_discriminator > btag_threshold_Medium)
							{
							NT_nMbjets_noVLooseTau +=1;
							NT_nMbjets_noMediumTau +=1;
							NT_nMbjets_noTightTau  +=1;
							}
						}
					else if (!matched_Medium)
						{
						NT_nbjets_noMediumTau += 1;
						NT_nbjets_noTightTau  += 1;
						if (b_discriminator > btag_threshold_Medium)
							{
							NT_nMbjets_noMediumTau +=1;
							NT_nMbjets_noTightTau  +=1;
							}
						}
					else if (!matched_Tight )
						{
						NT_nbjets_noTightTau  += 1;
						if (b_discriminator > btag_threshold_Medium)
							{
							NT_nMbjets_noTightTau  +=1;
							}
						}

					}
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
		jet.setP4(rawJet*jes_correction);
		// default jets are fully corrected, the initial slimmedJet is saved
		// the raw is saved too
		/*
		 * I checked this re-correction compared to initial miniaod jets
		 * -- the are equal, the miniaod jets correction is up to date
		 * only the JER is missing in MC.
		 */

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
			// here I use the jes-corrected jet pt
			// TODO: maybe the uncorrected jet is needed? (difference is not big)

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
	
			// gen jet info only for MC
			//NT_jet_matched_genjet_p4.push_back(gen_jet_p4);
			NT_genjet_matched  .push_back( genjet_matched );
			NT_genjet_pt       .push_back( genjet_pt      );
			NT_genjet_dR       .push_back( genjet_dR      );
			NT_genjet_i        .push_back( genjet_i       );

			LogInfo ("Demo") << "match gen to jets";
			LogInfo ("Demo") << "gen sizes " << gen_leps.size() << gen_taus.size() << gen_tau3ch.size() << gen_w_prods.size() << gen_b_prods.size();
			// match to GENPRODUCTS
			struct gen_matching match = match_to_gen(jet.p4(), gen_leps, gen_taus, gen_tau3ch, gen_taulep, gen_w_prods, gen_b_prods, gid_leps, gid_taus, gid_tau3ch, gid_taulep, gid_w_prods, gid_b_prods);
			NT_jet_matching_gen   .push_back(match.closest);
			NT_jet_matching_gen_dR.push_back(match.dR);
			LogInfo ("Demo") << "matched gen";
			}

		// propagate the jet correction to met and whatnot
		full_jet_corr += jet.p4() - jet_initial_p4; // initial jet + this difference = corrected jet

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
		selJets.push_back(jet);
		/* requirements for there jets: pt and no leptons in dR
		 * no PFID requirement and no eta
		 */

		// match to sel leptons
		//crossClean_in_dR(selJets, selLeptons, 0.4, selJets, weight, string("selJets"), false, false);
	        //bool overlapWithLepton(false);
		//float min_dR = 0.4;
	        //for(unsigned int l=0; l<(unsigned int)selLeptons.size();++l)
	        //        {
	        //        if (reco::deltaR(jet, selLeptons[l])<min_dR)
	        //                { overlapWithLepton=true; break; }
	        //        }
	        //if (overlapWithLepton) continue; // 222
		Float_t dR_to_sel = 999.;
		for(size_t l=0; l<selLeptons.size(); ++l)
			{
			double dR = reco::deltaR(jet, selLeptons[l]);
			if (dR < dR_to_sel)
				dR_to_sel = dR;
			}
		NT_jet_matching_lep    .push_back(dR_to_sel < 0.4);
		NT_jet_matching_lep_dR .push_back(dR_to_sel);

		// match to all-iso leptons
		Float_t dR_to_alliso = 999.;
		for(size_t l=0; l<selLeptons_allIso.size(); ++l)
			{
			double dR = reco::deltaR(jet, selLeptons_allIso[l]);
			if (dR < dR_to_alliso)
				dR_to_alliso = dR;
			}
		NT_jet_matching_allIso_lep    .push_back(dR_to_alliso < 0.4);
		NT_jet_matching_allIso_lep_dR .push_back(dR_to_alliso);

		}
	NT_nalljets  = selJets.size();

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

	//pat::JetCollection selJets;
	//crossClean_in_dR(selJets, selLeptons, 0.4, selJets, weight, string("selJets"), false, false);
	// and these are output jets for NTuple
	// they pass ID, corrected with JEC (smeared JES for MC)
	// pass kinematic cuts (pt, eta)
	// and dR-cleaned from selected leptons

	//std::sort (selJets.begin(),  selJets.end(),  utils::sort_CandidatesByPt);
	// no need for sort


	/*
	 * TAUS final
	 */
	//LogInfo("Demo") << "taus.size() = "<< taus.size();
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

	std::vector<double > tracksToBeRemoved_PV; // compared by Pt due to the conflict of comparing const and not const iterators
	// these tracks correspond to all 3pi tau tracks in the event
	// without considering the ID of these taus
	// whenever refit works for a tau the corresponding tracks are removed
	for(size_t i=0; i<selTaus.size(); ++i)
		{
		pat::Tau& tau = selTaus[i];

		//Int_t IDlev = 0;
		//if (tau.tauID(tau_VTight_ID)) IDlev = 5;
		//else if (tau.tauID(tau_Tight_ID))  IDlev = 4;
		//else if (tau.tauID(tau_Medium_ID)) IDlev = 3;
		//else if (tau.tauID(tau_Loose_ID))  IDlev = 2;
		//else if (tau.tauID(tau_VLoose_ID)) IDlev = 1;
		Int_t IDlev = tau.userInt("IDlev");

		//if (IDlev < 0) continue; // save only candidates with minimal ID
		if (IDlev > 0) NT_ntaus_idVL += 1;
		if (IDlev > 2) NT_ntaus_idM  += 1;
		if (IDlev > 3) NT_ntaus_idT  += 1;

		NT_tau_id.push_back(tau.pdgId());
		NT_tau_decayMode.push_back(tau.decayMode());
		NT_tau_p4.push_back(tau.p4());
		NT_tau_IDlev.push_back(IDlev);
		NT_tau_IDmedium_discr.push_back(tau.tauID(tau_Medium_ID));
		NT_tau_leading_track_pt.push_back(tau.userFloat("leading_track_pt"));
		NT_tau_leadChargedHadrCand_pt.push_back(tau.userFloat("leadChargedHadrCand_pt"));
		NT_tau_leadNeutralCand_pt.push_back(tau.userFloat("leadNeutralCand_pt"));
		NT_tau_leadCand_pt.push_back(tau.userFloat("leadCand_pt"));
		NT_tau_hasSecondaryVertex.push_back(tau.hasSecondaryVertex());
		//NT_tau_hcalEnergy = tau.hcalEnergy();
		//NT_tau_hcalEnergyLeadChargedHadrCand = tau.hcalEnergyLeadChargedHadrCand();

		// PAT isolation parameters
		NT_tau_n_isolationCands              .push_back(tau.isolationCands              () .size());
		NT_tau_n_isolationChargedHadrCands   .push_back(tau.isolationChargedHadrCands   () .size());
		NT_tau_n_isolationGammaCands         .push_back(tau.isolationGammaCands         () .size());
		NT_tau_n_isolationNeutrHadrCands     .push_back(tau.isolationNeutrHadrCands     () .size());

		NT_tau_n_signalCands              .push_back(tau.signalCands              () .size());
		NT_tau_n_signalChargedHadrCands   .push_back(tau.signalChargedHadrCands   () .size());
		NT_tau_n_signalGammaCands         .push_back(tau.signalGammaCands         () .size());
		NT_tau_n_signalNeutrHadrCands     .push_back(tau.signalNeutrHadrCands     () .size());

		// save the p4 momenta for tau candidates with an ID
		if (IDlev > 0)
			{
			//std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >> cands_p4;
			for(const auto& cand: tau.isolationCands())
				{
				//cands_p4.push_back(cand->p4());
				NT_tau_all_isolationCands_tau_index  .push_back(i);
				NT_tau_all_isolationCands            .push_back(cand->p4());
				}

			//cands_p4.clear();
			for(const auto& cand: tau.isolationChargedHadrCands())
				{
				NT_tau_all_isolationCharged_tau_index  .push_back(i);
				NT_tau_all_isolationCharged_pdgId      .push_back(cand->pdgId());
				NT_tau_all_isolationCharged            .push_back(cand->p4());
				}

			for(const auto& cand: tau.isolationGammaCands())
				{
				NT_tau_all_isolationGamma_tau_index  .push_back(i);
				NT_tau_all_isolationGamma            .push_back(cand->p4());
				}

			for(const auto& cand: tau.isolationNeutrHadrCands())
				{
				NT_tau_all_isolationNeutr_tau_index  .push_back(i);
				NT_tau_all_isolationNeutr_pdgId      .push_back(cand->pdgId());
				NT_tau_all_isolationNeutr            .push_back(cand->p4());
				}

			// again save the p4 momenta
			for(const auto& cand: tau.signalCands())
				{
				//cands_p4.push_back(cand->p4());
				NT_tau_all_signalCands_tau_index  .push_back(i);
				NT_tau_all_signalCands            .push_back(cand->p4());
				}

			//cands_p4.clear();
			for(const auto& cand: tau.signalChargedHadrCands())
				{
				NT_tau_all_signalCharged_tau_index  .push_back(i);
				NT_tau_all_signalCharged_pdgId      .push_back(cand->pdgId());
				NT_tau_all_signalCharged            .push_back(cand->p4());
				}

			for(const auto& cand: tau.signalGammaCands())
				{
				NT_tau_all_signalGamma_tau_index  .push_back(i);
				NT_tau_all_signalGamma            .push_back(cand->p4());
				}

			for(const auto& cand: tau.signalNeutrHadrCands())
				{
				NT_tau_all_signalNeutr_tau_index  .push_back(i);
				NT_tau_all_signalNeutr_pdgId      .push_back(cand->pdgId());
				NT_tau_all_signalNeutr            .push_back(cand->p4());
				}
			}

		// these have 0 in the output
		//NT_tau_n_isolationPFCands            .push_back(tau.isolationPFCands            () .size());
		//NT_tau_n_isolationPFChargedHadrCands .push_back(tau.isolationPFChargedHadrCands () .size());
		//NT_tau_n_isolationPFGammaCands       .push_back(tau.isolationPFGammaCands       () .size());
		//NT_tau_n_isolationPFNeutrHadrCands   .push_back(tau.isolationPFNeutrHadrCands   () .size());

		// match to sel leptons
		Float_t dR_to_sel = 999.;
		for(size_t l=0; l<selLeptons.size(); ++l)
			{
			double dR = reco::deltaR(tau, selLeptons[l]);
			if (dR < dR_to_sel)
				dR_to_sel = dR;
			}
		NT_tau_matching_lep    .push_back(dR_to_sel < 0.4);
		NT_tau_matching_lep_dR .push_back(dR_to_sel);

		// match to all-iso leptons
		Float_t dR_to_alliso = 999.;
		for(size_t l=0; l<selLeptons_allIso.size(); ++l)
			{
			double dR = reco::deltaR(tau, selLeptons_allIso[l]);
			if (dR < dR_to_alliso)
				dR_to_alliso = dR;
			}
		NT_tau_matching_allIso_lep    .push_back(dR_to_alliso < 0.4);
		NT_tau_matching_allIso_lep_dR .push_back(dR_to_alliso);

		if (isMC)
			{
			LogInfo ("Demo") << "gen match to tau";
			struct gen_matching match = match_to_gen(tau.p4(), gen_leps, gen_taus, gen_tau3ch, gen_taulep, gen_w_prods, gen_b_prods, gid_leps, gid_taus, gid_tau3ch, gid_taulep, gid_w_prods, gid_b_prods);
			NT_tau_matching_gen   .push_back(match.closest);
			NT_tau_matching_gen_dR.push_back(match.dR);

			if (isDY || isTT)
				{
				struct gen_matching_collection match2 = match_to_gen_collection(tau.p4(), NT_gen_tt_tau_vis_p4);
				NT_tau_matching_gen_collection   .push_back(match2.index);
				NT_tau_matching_gen_collection_dR.push_back(match2.dR);
				}
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
		for (unsigned int i=0; i<selJets.size(); i++)
			{
			pat::Jet& jet = selJets[i];
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
			int  track_pdgId = (*track_cands)[track_index].pdgId();
			auto distance = closest_point - ref_vertex.position();
			// distance is of class ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>
			//impact.SetXYZ(distance.x(), distance.y(), distance.z());

			NT_tau_SV_fit_track_OS_matched_track_vtxkey.push_back(key);
			NT_tau_SV_fit_track_OS_matched_track_vtxQ.push_back(quality);
			NT_tau_SV_fit_track_OS_matched_track_b.push_back(distance);
			NT_tau_SV_fit_track_OS_matched_track_pdgId.push_back(track_pdgId);

			// SS1
			track_index = matched_track_SS1;

			key = (*track_cands)[track_index].vertexRef().key();
			quality = (*track_cands)[track_index].pvAssociationQuality();

			ref_vertex = *((*track_cands)[track_index].vertexRef());
			closest_point = (*track_cands)[track_index].vertex();
			distance = closest_point - ref_vertex.position();
			// distance is of class ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>
			//impact.SetXYZ(distance.x(), distance.y(), distance.z());
			track_pdgId = (*track_cands)[track_index].pdgId();

			NT_tau_SV_fit_track_SS1_matched_track_vtxkey.push_back(key);
			NT_tau_SV_fit_track_SS1_matched_track_vtxQ.push_back(quality);
			NT_tau_SV_fit_track_SS1_matched_track_b.push_back(distance);
			NT_tau_SV_fit_track_SS1_matched_track_pdgId.push_back(track_pdgId);

			// SS2
			track_index = matched_track_SS2;

			key = (*track_cands)[track_index].vertexRef().key();
			quality = (*track_cands)[track_index].pvAssociationQuality();

			ref_vertex = *((*track_cands)[track_index].vertexRef());
			closest_point = (*track_cands)[track_index].vertex();
			distance = closest_point - ref_vertex.position();
			// distance is of class ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>
			//impact.SetXYZ(distance.x(), distance.y(), distance.z());
			track_pdgId = (*track_cands)[track_index].pdgId();

			NT_tau_SV_fit_track_SS2_matched_track_vtxkey.push_back(key);
			NT_tau_SV_fit_track_SS2_matched_track_vtxQ.push_back(quality);
			NT_tau_SV_fit_track_SS2_matched_track_b.push_back(distance);
			NT_tau_SV_fit_track_SS2_matched_track_pdgId.push_back(track_pdgId);

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
	bool pass_leptons = clean_lep_conditions && selLeptons.size() > 0 && selLeptons.size() < 3; // tau_ID oriented scheme, Dileptons are separate
	bool pass_leptons_all_iso = clean_lep_conditions && selLeptons_allIso.size() > 0 && selLeptons_allIso.size() < 3;
	bool record_ntuple = false;

	bool record_tauID_cond = pass_leptons && NT_ntaus > 0 && selLeptons.size() == 1;

	event_checkpoint++;

	event_checkpoint++;
	if (selElectrons.size() == 1 && record_tauID_cond)
		event_counter ->Fill(event_checkpoint);

	event_checkpoint++;
	if (selMuons.size() == 1 && record_tauID_cond)
		event_counter ->Fill(event_checkpoint);

	event_checkpoint++;

	if (record_tauCands)
		{
		record_ntuple |= NT_ntaus > 0;
		}

	if (record_tauID)
		{
		/*
		 * tau ID preselection
		 * about twice less events then in our preselection with b jets
		 *
		 * should contain good WJets control sample
		 */
		record_ntuple |= record_tauID_cond;
		}

	if (record_lepTauVL_b)
		{
		record_ntuple |= record_tauID_cond && NT_ntaus_idVL > 0 && NT_nbjets_noVLooseTau > 0;
		}

	if (record_ElTau)
		{
		record_ntuple |= record_tauID_cond && selElectrons.size() == 1;
		}
	if (record_MuTau)
		{
		record_ntuple |= record_tauID_cond && selMuons.size() == 1;
		}

	if (record_tauIDantiIso)
		{
		record_ntuple |= pass_leptons_all_iso && NT_ntaus > 0 && selLeptons_allIso.size() == 1;
		}
	if (record_bPreselection)
		{
		/* leptons and at least 1 b jet
		 * our old preselection
		 */
		record_ntuple |= pass_leptons && NT_nbjets > 0 && selLeptons.size() == 1;
		}
	if (record_MonitorHLT)
		{
		/* the HLT efficiency study
		 * only few events in lepMonitorTrigger
		 */
		record_ntuple |= pass_leptons && lepMonitorTrigger && selLeptons.size() == 1;
		}

	if (record_ElMu_b)
		{
		record_ntuple |= pass_leptons && abs(NT_leps_ID) == 143 && NT_nbjets_noMediumTau > 0 && NT_njets > 1; // TODO: this is our ooold selection, it is stupid -- there is no need to ask for non-b-tagged jets
		}
	if (record_ElMu)
		{
		// all el-mu events (it's mainly TTbar, so should be few)
		record_ntuple |= pass_leptons && abs(NT_leps_ID) == 143;
		}
	if (record_Dilep)
		{
		/* just all dileptons -- lots of DY here
		 * control for lepton IDs
		 * about = to tau ID preselection
		 */
		record_ntuple |= pass_leptons && selLeptons.size() == 2;
		}

	if (record_jets)
		{
		record_ntuple |= jetsHLT && selJets.size() > 0; // these are all-eta jets with pt-cut, anti-lep dR and Loose ID..
		}
	if (record_signal && isTTSignal)
		{
		record_ntuple = true; // all events of signal
		}
	if (record_all)
		{
		record_ntuple = true; // all events
		}

	if (record_ntuple)
		{
		LogInfo ("Demo") << "recording " << NT_indexevents << ' ' << record_tauCands << record_tauID << record_tauIDantiIso << record_bPreselection << record_MonitorHLT << record_ElMu << record_Dilep << record_jets << record_signal << isTTSignal;
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
GenDistrs_cleanup();

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
