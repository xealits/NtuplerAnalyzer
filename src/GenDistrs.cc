
//#include "TMath.h" // LorentzVector?
#include "UserCode/NtuplerAnalyzer/interface/GenDistrs.h"

struct lep_id_range {
	Int_t min;
	Int_t max;
};

struct dilep_id_ranges {
	struct lep_id_range r1;
	struct lep_id_range r2;
};

bool id_in_range(Int_t id, struct lep_id_range range)
	{
	Int_t absid = abs(id);
	return absid >= range.min && absid <= range.max;
	}

struct recorded_histos {
	TH1D* pts[2];
	TH1D* etas[2];
};

/** the definition of parameters recorded in each final state
*/
struct recorded_gen_histos {
	TH1D* h_lep_pts[2];
	TH1D* h_lep_etas[2];
	TH1D* h_bjet_pts[2];
	TH1D* h_bjet_etas[2];
};

class FinalStateProcess
{
	// process name: elel, eltau3h etc
	std::string proc_name;

	/* nope, let's set rigid structure
	// to allow for some flexibility
	proc_selection_func passes_function;
	histo_record_function record_function;
	*/
	struct dilep_id_ranges proc_selection_ranges;

	// the recorded gen-level parameters
	struct recorded_gen_histos params;
public:
	//FinalStateProcess(std::string proc_name, struct lep_id_range, struct lep_id_range);
	//FinalStateProcess(std::string proc_name, proc_selection_func, histo_record_function, edm::Service<TFileService>& fs);
	FinalStateProcess(std::string proc_name,
		struct dilep_id_ranges, // final state definition with leptons
		edm::Service<TFileService>& fs);

	// passes when both leptons fall within separate ranges, no matter which one where
	//bool passes(Int_t lep1_id, Int_t lep2_id)
	bool passes(struct GenDistrs_recorded_gen_objects objs)
		{
		UInt_t abs1 = abs(objs.lep_ids[0]);
		UInt_t abs2 = abs(objs.lep_ids[1]);
		return (id_in_range(abs1, proc_selection_ranges.r1) && id_in_range(abs2, proc_selection_ranges.r2)) ||
			(id_in_range(abs1, proc_selection_ranges.r2) && id_in_range(abs2, proc_selection_ranges.r1));
		}

	// record given gen-level objects: leptons, bjets, other jets etc
	// sort leptons by ID then by pT
	void record_histos(struct GenDistrs_recorded_gen_objects);
};

#define quick_fs_histo_pt(obj)  fs->make<TH1D>((proc_name + "_" #obj "_pt" ).c_str(), (proc_name + "_" #obj "_pt ").c_str(),  200,  0, 200)
#define quick_fs_histo_eta(obj) fs->make<TH1D>((proc_name + "_" #obj "_eta").c_str(), (proc_name + "_" #obj "_eta").c_str(),  200, -3.0, 3.0)

FinalStateProcess::FinalStateProcess(std::string name,
		struct dilep_id_ranges dilep_ranges,
		edm::Service<TFileService>& fs) :
proc_name(name),
proc_selection_ranges(dilep_ranges)
	{
	params.h_lep_pts[0]  = quick_fs_histo_pt(l0);
	params.h_lep_pts[1]  = quick_fs_histo_pt(l1);
	params.h_lep_etas[0] = quick_fs_histo_eta(l0);
	params.h_lep_etas[1] = quick_fs_histo_eta(l1);

	params.h_bjet_pts[0]  = quick_fs_histo_pt(b0);
	params.h_bjet_pts[1]  = quick_fs_histo_pt(b1);
	params.h_bjet_etas[0] = quick_fs_histo_eta(b0);
	params.h_bjet_etas[1] = quick_fs_histo_eta(b1);
	}

#define define_gen_level_pt(name)  name ##_pt  = fs->make<TH1D>( #name "_pt",  #name "_pt",  200,  0, 200)
#define define_gen_level_eta(name) name ##_eta = fs->make<TH1D>( #name "_eta", #name "_eta", 200,  -3.0, 3.0)

std::vector<FinalStateProcess*> all_final_states;

int GenDistrs_make_histos_in_FileService(edm::Service<TFileService>& fs)
	{
	/*
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
	*/
	//std::vector<FinalStateProcess*> all_final_states;

	FinalStateProcess* fstate;
	fstate = new FinalStateProcess("elel",   {{11,11}, {11,11}} , fs);
	all_final_states.push_back(fstate);
	fstate = new FinalStateProcess("eltaul", {{11,11}, {15*11,15*11}}, fs);
	all_final_states.push_back(fstate);
	//fstate = new FinalStateProcess("taulel", {{11,11}, {11,11}} );
	//all_final_states.push_back(fstate);
	fstate = new FinalStateProcess("eltau1h", {{11,11}, {15*15,15*29}}, fs);
	all_final_states.push_back(fstate);
	fstate = new FinalStateProcess("eltau3h", {{11,11}, {15*30,1000*1000}}, fs);
	all_final_states.push_back(fstate);

	return 0; // Success
	}

int GenDistrs_cleanup(void)
	{
	for (unsigned int i=0; i<all_final_states.size(); i++)
		{
		delete all_final_states[i]; all_final_states[i] = NULL;
		}
	return 0;
	}

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

int record_a_pair_sorted_by_pt(const LorentzVector* v1, const LorentzVector* v2, struct recorded_histos histos)
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
		struct recorded_histos histos)
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

//int GenDistrs_record(struct GenDistrs_general_gen_params ps)
int GenDistrs_record(struct GenDistrs_recorded_gen_objects objs)
	{
	/*
	auto gen_decay_lep1_id = ps.gen_decay_lep1_id;
	auto gen_decay_lep2_id = ps.gen_decay_lep2_id;
	auto NT_gen_decay_lep1_p4 = ps.NT_gen_decay_lep1_p4;
	auto NT_gen_decay_lep2_p4 = ps.NT_gen_decay_lep2_p4;
	auto NT_gen_decay_jet1_p4 = ps.NT_gen_decay_jet1_p4;
	auto NT_gen_decay_jet2_p4 = ps.NT_gen_decay_jet2_p4;
	auto NT_gen_decay_bjet1_p4 = ps.NT_gen_decay_bjet1_p4;
	auto NT_gen_decay_bjet2_p4 = ps.NT_gen_decay_bjet2_p4;
	*/

	for (unsigned int i=0; i<all_final_states.size(); i++)
		{
		if (all_final_states[i]->passes(objs))
			{
			all_final_states[i]->record_histos(objs);
			}
		}

	/*
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
	*/

	return 0;
	}

void FinalStateProcess::record_histos(struct GenDistrs_recorded_gen_objects gen_obj)
	{
	// sort leptons by IDs if they are different
	if (abs(gen_obj.lep_ids[0]) != abs(gen_obj.lep_ids[1]))
		{
		record_leptons_sorted_by_ID(abs(gen_obj.lep_ids[0]), abs(gen_obj.lep_ids[1]),
			gen_obj.leps[0], gen_obj.leps[1],
			//{.pts=params.h_lep_pts, .etas=params.h_lep_etas});
			{{params.h_lep_pts[0], params.h_lep_pts[1]}, {params.h_lep_etas[0], params.h_lep_etas[1]}});
		}
	else
	// otherwise by pT
		{
		record_a_pair_sorted_by_pt(
			gen_obj.leps[0], gen_obj.leps[1],
			{{params.h_lep_pts[0], params.h_lep_pts[1]}, {params.h_lep_etas[0], params.h_lep_etas[1]}});
			//{.pts=params.h_lep_pts, .etas=params.h_lep_etas});
		}
	// sort jets by pT
	record_a_pair_sorted_by_pt(gen_obj.bjets[0], gen_obj.bjets[1],
		{{params.h_bjet_pts[0], params.h_bjet_pts[1]}, {params.h_bjet_etas[0], params.h_bjet_etas[1]}});
		//{.pts=params.h_bjet_pts, .etas=params.h_bjet_etas});
	}
