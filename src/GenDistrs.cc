
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
	TH1D* h_lep_control_ids[2];
	TH1D* h_lep_pts[2];
	TH2D* h_lep_pt_pt;
	TH1D* h_lep_etas[2];
	TH2D* h_lep_pts_etas[2];
	TH1D* h_jet_pts[2];
	TH1D* h_jet_etas[2];
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

#define quick_fs_histo_control_id(obj)     fs->make<TH1D>((proc_name + "_" #obj "_id" ).c_str(), (proc_name + "_" #obj "_id ").c_str(),  200,  0, 200)
#define quick_fs_histo_pt(obj)     fs->make<TH1D>((proc_name + "_" #obj "_pt" ).c_str(), (proc_name + "_" #obj "_pt ").c_str(),  200,  0, 200)
#define quick_fs_histo_eta(obj)    fs->make<TH1D>((proc_name + "_" #obj "_eta").c_str(), (proc_name + "_" #obj "_eta").c_str(),  200, -3.0, 3.0)
#define quick_fs_histo_pt_eta(obj) fs->make<TH2D>((proc_name + "_" #obj "_pt_eta").c_str(), (proc_name + "_" #obj "_pt_eta").c_str(), 40, 0, 200,  40, -3.0, 3.0)
#define quick_fs_histo_pt_pt(obj)  fs->make<TH2D>((proc_name + "_" #obj "_pt_pt").c_str(), (proc_name + "_" #obj "_pt_pt").c_str(), 40, 0, 200,  40, 0, 200)

FinalStateProcess::FinalStateProcess(std::string name,
		struct dilep_id_ranges dilep_ranges,
		edm::Service<TFileService>& fs) :
proc_name(name),
proc_selection_ranges(dilep_ranges)
	{
	params.h_lep_control_ids[0]  = quick_fs_histo_control_id  (l0);
	params.h_lep_control_ids[1]  = quick_fs_histo_control_id  (l1);
	params.h_lep_pts[0]  = quick_fs_histo_pt  (l0);
	params.h_lep_pts[1]  = quick_fs_histo_pt  (l1);
	params.h_lep_pt_pt   = quick_fs_histo_pt_pt (l0l1);
	params.h_lep_etas[0] = quick_fs_histo_eta (l0);
	params.h_lep_etas[1] = quick_fs_histo_eta (l1);
	params.h_lep_pts_etas[0] = quick_fs_histo_pt_eta (l0);
	params.h_lep_pts_etas[1] = quick_fs_histo_pt_eta (l1);

	params.h_jet_pts[0]  = quick_fs_histo_pt  (j0);
	params.h_jet_pts[1]  = quick_fs_histo_pt  (j1);
	params.h_jet_etas[0] = quick_fs_histo_eta (j0);
	params.h_jet_etas[1] = quick_fs_histo_eta (j1);

	params.h_bjet_pts[0]  = quick_fs_histo_pt  (b0);
	params.h_bjet_pts[1]  = quick_fs_histo_pt  (b1);
	params.h_bjet_etas[0] = quick_fs_histo_eta (b0);
	params.h_bjet_etas[1] = quick_fs_histo_eta (b1);
	}

void FinalStateProcess::record_histos(struct GenDistrs_recorded_gen_objects gen_obj)
	{
	unsigned int leading_lep, sublead_lep;
	// sort leptons by IDs if they are different
	if (abs(gen_obj.lep_ids[0]) != abs(gen_obj.lep_ids[1]))
		{
		//record_leptons_sorted_by_ID(abs(gen_obj.lep_ids[0]), abs(gen_obj.lep_ids[1]),
		//	gen_obj.leps[0], gen_obj.leps[1],
		//	{{params.h_lep_pts[0], params.h_lep_pts[1]}, {params.h_lep_etas[0], params.h_lep_etas[1]}});
		//histos.pts[lowest_id]  ->Fill(v1->pt());
		//histos.pts[another_id] ->Fill(v2->pt());

		unsigned int lowest_id = abs(gen_obj.lep_ids[0]) < abs(gen_obj.lep_ids[1]) ? 0 : 1;
		auto another_id = 1 - lowest_id;
		// many variables with names to describe the semantics of what is done
		leading_lep = lowest_id;
		sublead_lep = another_id;
		}

	else
	// otherwise by pT
		{
		//record_a_pair_sorted_by_pt(
		//	gen_obj.leps[0], gen_obj.leps[1],
		//	{{params.h_lep_pts[0], params.h_lep_pts[1]}, {params.h_lep_etas[0], params.h_lep_etas[1]}});

		double pts[2] = {gen_obj.leps[0]->pt(), gen_obj.leps[1]->pt()};

		unsigned int larger_pt = pts[0] > pts[1] ? 0 : 1;
		auto another_pt = 1 - larger_pt;
		leading_lep = larger_pt;
		sublead_lep = another_pt;
		}

	params.h_lep_control_ids[leading_lep] ->Fill(gen_obj.lep_ids[0]);
	params.h_lep_control_ids[sublead_lep] ->Fill(gen_obj.lep_ids[1]);
	params.h_lep_pts[leading_lep] ->Fill(gen_obj.leps[0]->pt());
	params.h_lep_pts[sublead_lep] ->Fill(gen_obj.leps[1]->pt());
	params.h_lep_pt_pt ->Fill(gen_obj.leps[0]->pt(), gen_obj.leps[1]->pt());
	params.h_lep_etas[leading_lep] ->Fill(gen_obj.leps[0]->eta());
	params.h_lep_etas[sublead_lep]->Fill(gen_obj.leps[1]->eta());
	params.h_lep_pts_etas[leading_lep]  ->Fill(gen_obj.leps[0]->pt(), gen_obj.leps[0]->eta());
	params.h_lep_pts_etas[sublead_lep] ->Fill(gen_obj.leps[1]->pt(), gen_obj.leps[1]->eta());

	// sort b jets by pT
	record_a_pair_sorted_by_pt(gen_obj.bjets[0], gen_obj.bjets[1],
		{{params.h_bjet_pts[0], params.h_bjet_pts[1]}, {params.h_bjet_etas[0], params.h_bjet_etas[1]}});
	// same for other jets
	record_a_pair_sorted_by_pt(gen_obj.jets[0], gen_obj.jets[1],
		{{params.h_jet_pts[0], params.h_jet_pts[1]}, {params.h_jet_etas[0], params.h_jet_etas[1]}});
	}

std::vector<FinalStateProcess*> all_final_states;

int GenDistrs_make_histos_in_FileService(edm::Service<TFileService>& fs)
	{
	FinalStateProcess* fstate;
	fstate = new FinalStateProcess("elel",   {{11,11}, {11,11}} , fs);
	all_final_states.push_back(fstate);
	fstate = new FinalStateProcess("mumu",   {{13,13}, {13,13}} , fs);
	all_final_states.push_back(fstate);
	fstate = new FinalStateProcess("elmu",   {{11,11}, {13,13}} , fs);
	all_final_states.push_back(fstate);

	fstate = new FinalStateProcess("eltaul", {{11,11}, {15*11,15*11}}, fs);
	all_final_states.push_back(fstate);
	//fstate = new FinalStateProcess("taulel", {{11,11}, {11,11}} );
	//all_final_states.push_back(fstate);

	fstate = new FinalStateProcess("eltau1h", {{11,11}, {15*15,15*29}}, fs);
	all_final_states.push_back(fstate);
	fstate = new FinalStateProcess("eltau3h", {{11,11}, {15*30,1000*1000}}, fs);
	all_final_states.push_back(fstate);
	fstate = new FinalStateProcess("mutau1h", {{13,13}, {15*15,15*29}}, fs);
	all_final_states.push_back(fstate);
	fstate = new FinalStateProcess("mutau3h", {{13,13}, {15*30,1000*1000}}, fs);
	all_final_states.push_back(fstate);

	// lepton+jets, in ttbar and wjets
	fstate = new FinalStateProcess("el",   {{1,1}, {11,11}} , fs);
	all_final_states.push_back(fstate);
	fstate = new FinalStateProcess("mu",   {{1,1}, {13,13}} , fs);
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

