
//#include "TMath.h" // LorentzVector?
#include "UserCode/NtuplerAnalyzer/interface/GenDistrs.h"

struct LorentzVector_pointer_pair {
	const LorentzVector* first;
	const LorentzVector* second;
};

//struct LorentzVector_pointer_pair sorted_byPt_LorentzVectors(const LorentzVector& v1, const LorentzVector& v2);
//void genDistrs_fill_sorted_pair(const LorentzVector& v1, const LorentzVector& v2);

struct LorentzVector_pointer_pair sorted_byPt_LorentzVectors(const LorentzVector* v1, const LorentzVector* v2)
	{
	struct LorentzVector_pointer_pair res;

	auto pt1 = v1->pt();
	auto pt2 = v2->pt();
	auto leading_p    = pt1 > pt2 ? v1 : v2;
	auto subleading_p = pt1 < pt2 ? v1 : v2;
	res.first  = leading_p;
	res.second = subleading_p;

	return res;
	}

#define define_gen_level_pt(name)  name ##_pt  = fs->make<TH1D>( #name "_pt",  #name "_pt",  200,  0, 200)
#define define_gen_level_eta(name) name ##_eta = fs->make<TH1D>( #name "_eta", #name "_eta", 200,  -3.0, 3.0)

int GenDistrs_make_histos_in_FileService(edm::Service<TFileService>& fs)
	{
	define_gen_level_pt(elel_l1);
	define_gen_level_pt(elel_l2);
	define_gen_level_pt(elel_b1);
	define_gen_level_pt(elel_b2);
	define_gen_level_eta(elel_l1);
	define_gen_level_eta(elel_l2);
	define_gen_level_eta(elel_b1);
	define_gen_level_eta(elel_b2);

	define_gen_level_pt(eltaul_l1);
	define_gen_level_pt(eltaul_l2);
	define_gen_level_pt(eltaul_b1);
	define_gen_level_pt(eltaul_b2);
	define_gen_level_eta(eltaul_l1);
	define_gen_level_eta(eltaul_l2);
	define_gen_level_eta(eltaul_b1);
	define_gen_level_eta(eltaul_b2);

	define_gen_level_pt(taulel_l1);
	define_gen_level_pt(taulel_l2);
	define_gen_level_pt(taulel_b1);
	define_gen_level_pt(taulel_b2);
	define_gen_level_eta(taulel_l1);
	define_gen_level_eta(taulel_l2);
	define_gen_level_eta(taulel_b1);
	define_gen_level_eta(taulel_b2);

	define_gen_level_pt(eltau1h_l1);
	define_gen_level_pt(eltau1h_l2);
	define_gen_level_pt(eltau1h_b1);
	define_gen_level_pt(eltau1h_b2);
        define_gen_level_eta(eltau1h_l1);
        define_gen_level_eta(eltau1h_l2);
        define_gen_level_eta(eltau1h_b1);
        define_gen_level_eta(eltau1h_b2);

        define_gen_level_pt(eltau3h_l1);
        define_gen_level_pt(eltau3h_l2);
        define_gen_level_pt(eltau3h_b1);
        define_gen_level_pt(eltau3h_b2);
        define_gen_level_eta(eltau3h_l1);
        define_gen_level_eta(eltau3h_l2);
        define_gen_level_eta(eltau3h_b1);
        define_gen_level_eta(eltau3h_b2);

	return 0; // Success
	}

struct record_histos {
	TH1D* pts[2];
	TH1D* etas[2];
};
int record_a_pair_sorted_by_pt(const LorentzVector* v1, const LorentzVector* v2, struct record_histos histos)
	{
	// sort by pT
	struct LorentzVector_pointer_pair sorted_p4s = sorted_byPt_LorentzVectors(v1, v2);
	// save
	histos.pts[0]->Fill(sorted_p4s.first->pt());
	histos.pts[1]->Fill(sorted_p4s.second->pt());
	histos.etas[0]->Fill(sorted_p4s.first->eta());
	histos.etas[1]->Fill(sorted_p4s.second->eta());
	return 0; //Success
	}

/** sort by PDG ID

does not merge different tauh yet!
*/
int record_leptons_sorted_by_ID(const Int_t id1, const Int_t id2,
		const LorentzVector* v1, const LorentzVector* v2,
		struct record_histos histos)
	{
	//
	unsigned int lowest_id = id1 < id2 ? 0 : 1;
	auto another_id = 1 - lowest_id;
	histos.pts[lowest_id]  ->Fill(v1->pt());
	histos.pts[another_id] ->Fill(v2->pt());
	histos.etas[lowest_id] ->Fill(v1->eta());
	histos.etas[another_id]->Fill(v2->eta());
	return 0; //Success
	}

int GenDistrs_record(struct GenDistrs_general_gen_params ps)
	{
	auto gen_decay_lep1_id = ps.gen_decay_lep1_id;
	auto gen_decay_lep2_id = ps.gen_decay_lep2_id;
	auto NT_gen_decay_lep1_p4 = ps.NT_gen_decay_lep1_p4;
	auto NT_gen_decay_lep2_p4 = ps.NT_gen_decay_lep2_p4;
	auto NT_gen_decay_jet1_p4 = ps.NT_gen_decay_jet1_p4;
	auto NT_gen_decay_jet2_p4 = ps.NT_gen_decay_jet2_p4;
	auto NT_gen_decay_bjet1_p4 = ps.NT_gen_decay_bjet1_p4;
	auto NT_gen_decay_bjet2_p4 = ps.NT_gen_decay_bjet2_p4;

	if (abs(gen_decay_lep1_id) == 11 && abs(gen_decay_lep2_id) == 11)
		{
		// sort leptons by pT
		record_a_pair_sorted_by_pt(NT_gen_decay_lep1_p4, NT_gen_decay_lep2_p4,
			{{elel_l1_pt, elel_l2_pt}, {elel_l1_eta, elel_l2_eta}});
		// and b jets
		record_a_pair_sorted_by_pt(NT_gen_decay_bjet1_p4, NT_gen_decay_bjet2_p4,
			{{elel_b1_pt, elel_b2_pt}, {elel_b1_eta, elel_b2_eta}});
		}

	// el tau->lepton
	if (abs(gen_decay_lep1_id) == 11 && (abs(gen_decay_lep2_id) == 15*11 || abs(gen_decay_lep2_id) == 15*13))
		{
		//eltaul_l1_pt ->Fill(NT_gen_decay_lep1_p4->pt());
		//eltaul_l2_pt->Fill(NT_gen_decay_lep2_p4->pt());
		//eltaul_l1_eta ->Fill(NT_gen_decay_lep1_p4->eta());
		//eltaul_l2_eta->Fill(NT_gen_decay_lep2_p4->eta());
		record_leptons_sorted_by_ID(abs(gen_decay_lep1_id), abs(gen_decay_lep2_id),
			NT_gen_decay_lep1_p4, NT_gen_decay_lep2_p4,
			{{eltaul_l1_pt, eltaul_l2_pt}, {eltaul_l1_eta, eltaul_l2_eta}});

		// and b jets
		record_a_pair_sorted_by_pt(NT_gen_decay_bjet1_p4, NT_gen_decay_bjet2_p4,
			{{eltaul_b1_pt, eltaul_b2_pt}, {eltaul_b1_eta, eltaul_b2_eta}});
		}

	// to test that this works lets separately save the taul el case
	if (abs(gen_decay_lep2_id) == 11 && (abs(gen_decay_lep1_id) == 15*11 || abs(gen_decay_lep1_id) == 15*13))
		{
		record_leptons_sorted_by_ID(abs(gen_decay_lep1_id), abs(gen_decay_lep2_id),
			NT_gen_decay_lep1_p4, NT_gen_decay_lep2_p4,
			{{taulel_l1_pt, taulel_l2_pt}, {taulel_l1_eta, taulel_l2_eta}});

		// and b jets
		record_a_pair_sorted_by_pt(NT_gen_decay_bjet1_p4, NT_gen_decay_bjet2_p4,
			{{taulel_b1_pt, taulel_b2_pt}, {taulel_b1_eta, taulel_b2_eta}});
		}


	// el tau->1charged
	if (abs(gen_decay_lep1_id) == 11 && (abs(gen_decay_lep2_id) > 15*15 && abs(gen_decay_lep2_id) < 15*30))
		{
		eltau1h_l1_pt ->Fill(NT_gen_decay_lep1_p4->pt());
		eltau1h_l2_pt->Fill(NT_gen_decay_lep2_p4->pt());
		eltau1h_l1_eta ->Fill(NT_gen_decay_lep1_p4->eta());
		eltau1h_l2_eta->Fill(NT_gen_decay_lep2_p4->eta());

		// and b jets
		record_a_pair_sorted_by_pt(NT_gen_decay_bjet1_p4, NT_gen_decay_bjet2_p4,
			{{eltau1h_b1_pt, eltau1h_b2_pt}, {eltau1h_b1_eta, eltau1h_b2_eta}});
		}

	//el tau -> 3charged
	// practically the same as 3charged hadrons
	// TODO merge these procedures
	if (abs(gen_decay_lep1_id) == 11 && abs(gen_decay_lep2_id) >= 15*30)
		{
		eltau3h_l1_pt ->Fill(NT_gen_decay_lep1_p4->pt());
		eltau3h_l2_pt->Fill(NT_gen_decay_lep2_p4->pt());
		eltau3h_l1_eta ->Fill(NT_gen_decay_lep1_p4->eta());
		eltau3h_l2_eta->Fill(NT_gen_decay_lep2_p4->eta());

		// and b jets
		record_a_pair_sorted_by_pt(NT_gen_decay_bjet1_p4, NT_gen_decay_bjet2_p4,
			{{eltau3h_b1_pt, eltau3h_b2_pt}, {eltau3h_b1_eta, eltau3h_b2_eta}});
		}

	return 0;
	}
