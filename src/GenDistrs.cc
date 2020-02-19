
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

