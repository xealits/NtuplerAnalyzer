
//#include "TMath.h" // LorentzVector?
#include "UserCode/NtuplerAnalyzer/interface/GenDistrs.h"

struct LorentzVector_pointer_pair sorted_byPt_LorentzVectors(const LorentzVector& v1, const LorentzVector& v2)
	{
	struct LorentzVector_pointer_pair res;

	auto pt1 = v1.pt();
	auto pt2 = v2.pt();
	auto leading_p    = pt1 > pt2 ? &v1 : &v2;
	auto subleading_p = pt1 < pt2 ? &v1 : &v2;
	res.first  = leading_p;
	res.second = subleading_p;

	return res;
	}

#define define_gen_level_pt(name)  name ##_pt  = fs->make<TH1D>( #name "_pt",  #name "_pt",  200,  0, 200)
#define define_gen_level_eta(name) name ##_eta = fs->make<TH1D>( #name "_eta", #name "_eta", 200,  -3.0, 3.0)

int genDistrs_make_histos_in_FileService(edm::Service<TFileService>& fs)
	{
	define_gen_level_pt(elel_el1);
	define_gen_level_pt(elel_el2);
	define_gen_level_pt(elel_b1);
	define_gen_level_pt(elel_b2);
	define_gen_level_eta(elel_el1);
	define_gen_level_eta(elel_el2);
	define_gen_level_eta(elel_b1);
	define_gen_level_eta(elel_b2);

	define_gen_level_pt(eltaul_el);
	define_gen_level_pt(eltaul_tau);
	define_gen_level_pt(eltaul_b1);
	define_gen_level_pt(eltaul_b2);
	define_gen_level_eta(eltaul_el);
	define_gen_level_eta(eltaul_tau);
	define_gen_level_eta(eltaul_b1);
	define_gen_level_eta(eltaul_b2);

	define_gen_level_pt(eltau1h_el);
	define_gen_level_pt(eltau1h_tau);
	define_gen_level_pt(eltau1h_b1);
	define_gen_level_pt(eltau1h_b2);
        define_gen_level_eta(eltau1h_el);
        define_gen_level_eta(eltau1h_tau);
        define_gen_level_eta(eltau1h_b1);
        define_gen_level_eta(eltau1h_b2);

        define_gen_level_pt(eltau3h_el);
        define_gen_level_pt(eltau3h_tau);
        define_gen_level_pt(eltau3h_b1);
        define_gen_level_pt(eltau3h_b2);
        define_gen_level_eta(eltau3h_el);
        define_gen_level_eta(eltau3h_tau);
        define_gen_level_eta(eltau3h_b1);
        define_gen_level_eta(eltau3h_b2);

	return 0; // Success
	}

