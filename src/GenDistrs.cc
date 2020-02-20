
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

struct LorentzVector_pointer_pair {
	const LorentzVector* first;
	const LorentzVector* second;
};

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


class FinalStateProcess
{
	// process name: elel, eltau3h etc
	std::string proc_name;

	// nope, a rigid structure with inclusive ranges for abs(PDG ID)
	struct dilep_id_ranges proc_selection_ranges;

	// the recorded gen-level parameters
	struct recorded_gen_histos params;
public:
	FinalStateProcess(std::string proc_name,
		struct dilep_id_ranges, // final state definition with leptons
		edm::Service<TFileService>& fs);

	// passes when both leptons fall within separate ranges, no matter which one where
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

void FinalStateProcess::record_histos(struct GenDistrs_recorded_gen_objects gen_obj)
	{
	// sort leptons by IDs if they are different
	if (abs(gen_obj.lep_ids[0]) != abs(gen_obj.lep_ids[1]))
		{
		record_leptons_sorted_by_ID(abs(gen_obj.lep_ids[0]), abs(gen_obj.lep_ids[1]),
			gen_obj.leps[0], gen_obj.leps[1],
			{{params.h_lep_pts[0], params.h_lep_pts[1]}, {params.h_lep_etas[0], params.h_lep_etas[1]}});
		}
	else
	// otherwise by pT
		{
		record_a_pair_sorted_by_pt(
			gen_obj.leps[0], gen_obj.leps[1],
			{{params.h_lep_pts[0], params.h_lep_pts[1]}, {params.h_lep_etas[0], params.h_lep_etas[1]}});
		}
	// sort jets by pT
	record_a_pair_sorted_by_pt(gen_obj.bjets[0], gen_obj.bjets[1],
		{{params.h_bjet_pts[0], params.h_bjet_pts[1]}, {params.h_bjet_etas[0], params.h_bjet_etas[1]}});
	}

std::vector<FinalStateProcess*> all_final_states;

int GenDistrs_make_histos_in_FileService(edm::Service<TFileService>& fs)
	{
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

int GenDistrs_record(struct GenDistrs_recorded_gen_objects objs)
	{
	for (unsigned int i=0; i<all_final_states.size(); i++)
		{
		if (all_final_states[i]->passes(objs))
			{
			all_final_states[i]->record_histos(objs);
			}
		}

	return 0;
	}

