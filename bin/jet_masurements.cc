#include <iostream>                                                                                                                   

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TH1F.h"                                                                                                                      
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TNtuple.h"

#include "TMath.h"
#include <Math/Expression.h>
#include <Math/BinaryOperators.h>
#include <Math/VectorUtil.h>


#include <map>
#include <string>
#include <vector>

//#include "dtag_xsecs.h"
//#define INPUT_DTAGS_START 9

#define N_ITEMS(array) sizeof(array)/sizeof(array[0])

using namespace std;

float b_tag_wp_medium = 0.8484;

//static float bin_edges_pt  [] = { 0, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 65, 70, 75, 80, 100, 150, 500 };
static float bin_edges_pt  [] = { 0, 20, 500 };
static float bin_edges_eta [] = { -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5 };

size_t n_bin_edges_pt  = N_ITEMS(bin_edges_pt)  - 1;
size_t n_bin_edges_eta = N_ITEMS(bin_edges_eta) - 1;


/*
 * measurement histograms
 */

TH2D* btag_b_hadronFlavour_candidates = new TH2D( "btag_b_hadronFlavour_candidates"           , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_eta, bin_edges_eta);
TH2D* btag_b_hadronFlavour_candidates_tagged = new TH2D( "btag_b_hadronFlavour_candidates_tagged"    , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_eta, bin_edges_eta);
TH2D* btag_c_hadronFlavour_candidates = new TH2D( "btag_c_hadronFlavour_candidates"           , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_eta, bin_edges_eta);
TH2D* btag_c_hadronFlavour_candidates_tagged = new TH2D( "btag_c_hadronFlavour_candidates_tagged"    , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_eta, bin_edges_eta);
TH2D* btag_udsg_hadronFlavour_candidates = new TH2D( "btag_udsg_hadronFlavour_candidates"        , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_eta, bin_edges_eta);
TH2D* btag_udsg_hadronFlavour_candidates_tagged = new TH2D( "btag_udsg_hadronFlavour_candidates_tagged" , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_eta, bin_edges_eta);

TH2D* jes_cor_2d_gen_reco_barrel_b = new TH2D( "jes_cor_2d_gen_reco_barrel_b" , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_pt, bin_edges_pt);
TH2D* jes_cor_2d_gen_reco_endcap_b = new TH2D( "jes_cor_2d_gen_reco_endcap_b" , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_pt, bin_edges_pt);

TH2D* jes_cor_2d_gen_reco_barrel_c = new TH2D( "jes_cor_2d_gen_reco_barrel_c" , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_pt, bin_edges_pt);
TH2D* jes_cor_2d_gen_reco_endcap_c = new TH2D( "jes_cor_2d_gen_reco_endcap_c" , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_pt, bin_edges_pt);

TH2D* jes_cor_2d_gen_reco_barrel_udsg = new TH2D( "jes_cor_2d_gen_reco_barrel_udsg" , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_pt, bin_edges_pt);
TH2D* jes_cor_2d_gen_reco_endcap_udsg = new TH2D( "jes_cor_2d_gen_reco_endcap_udsg" , "", n_bin_edges_pt,  bin_edges_pt, n_bin_edges_pt, bin_edges_pt);

TH1D* jes_cor_1d_gen_reco_barrel_b = new TH1D( "jes_cor_1d_gen_reco_barrel_b" , "", n_bin_edges_pt,  bin_edges_pt);
TH1D* jes_cor_1d_gen_reco_endcap_b = new TH1D( "jes_cor_1d_gen_reco_endcap_b" , "", n_bin_edges_pt,  bin_edges_pt);
TH1D* jes_cor_1d_gen_reco_barrel_b_unit = new TH1D( "jes_cor_1d_gen_reco_barrel_b_unit" , "", n_bin_edges_pt,  bin_edges_pt);
TH1D* jes_cor_1d_gen_reco_endcap_b_unit = new TH1D( "jes_cor_1d_gen_reco_endcap_b_unit" , "", n_bin_edges_pt,  bin_edges_pt);

TH1D* jes_cor_1d_gen_reco_barrel_c = new TH1D( "jes_cor_1d_gen_reco_barrel_c" , "", n_bin_edges_pt,  bin_edges_pt);
TH1D* jes_cor_1d_gen_reco_endcap_c = new TH1D( "jes_cor_1d_gen_reco_endcap_c" , "", n_bin_edges_pt,  bin_edges_pt);
TH1D* jes_cor_1d_gen_reco_barrel_c_unit = new TH1D( "jes_cor_1d_gen_reco_barrel_c_unit" , "", n_bin_edges_pt,  bin_edges_pt);
TH1D* jes_cor_1d_gen_reco_endcap_c_unit = new TH1D( "jes_cor_1d_gen_reco_endcap_c_unit" , "", n_bin_edges_pt,  bin_edges_pt);

TH1D* jes_cor_1d_gen_reco_barrel_udsg = new TH1D( "jes_cor_1d_gen_reco_barrel_udsg" , "", n_bin_edges_pt,  bin_edges_pt);
TH1D* jes_cor_1d_gen_reco_endcap_udsg = new TH1D( "jes_cor_1d_gen_reco_endcap_udsg" , "", n_bin_edges_pt,  bin_edges_pt);
TH1D* jes_cor_1d_gen_reco_barrel_udsg_unit = new TH1D( "jes_cor_1d_gen_reco_barrel_udsg_unit" , "", n_bin_edges_pt,  bin_edges_pt);
TH1D* jes_cor_1d_gen_reco_endcap_udsg_unit = new TH1D( "jes_cor_1d_gen_reco_endcap_udsg_unit" , "", n_bin_edges_pt,  bin_edges_pt);


//int stacked_histo_distr (int argc, char *argv[])
int main (int argc, char *argv[])
{
/*
char usage_string[128] = "[--verbose] [--normalize] lumi distr projection rebin_factor x_axis_min_range x_axis_max_range name_tag dir dtags";
if (argc < 8)
    {
    std::cout << "Usage : " << argv[0] << usage_string << std::endl;
    return 1;
    }
*/

char usage_string[128] = " [--[test][verbose]] output_file input_file [input_file..]";
if (argc < 1 + 2)
    {
    std::cout << "Usage : " << argv[0] << usage_string << std::endl;
    return 1;
    }

gROOT->Reset();

//// printout the arguments
//printf("argc %5d\n", argc);
//for (int i=0; i<argc; i++)
//    printf("%5d %s\n", i, argv[i]);

unsigned int input_starts = 1;

// check for options
string inp1(argv[1]);
bool run_test = false;
bool be_verbose = false;

if (inp1.find("--") != std::string::npos)
	// then it's an option string
	{
	input_starts += 1;
	if (inp1.find("test")    != std::string::npos) run_test = true;
	if (inp1.find("verbose") != std::string::npos) be_verbose = true;
	}

TString output_file(argv[input_starts + 0]);
//TString inp_file(argv[input_starts + 1]);

for (int inp_f = input_starts + 1; inp_f < argc; inp_f++)
{
if (be_verbose) printf("inp_f %5d %s\n", inp_f, argv[inp_f]);
TString inp_file(argv[inp_f]);

TFile* tfile = TFile::Open(inp_file);

if (tfile == NULL)
	{
	std::cout << "No input file: " << inp_file << std::endl;
	return 2;
	}

TTree* NT_output_ttree = (TTree*) tfile->Get("ntupler/reduced_ttree");

// connect the branch with macro:
#define NTUPLE_INTERFACE_OPEN
#include "UserCode/NtuplerAnalyzer/interface/noSV_ntupleOutput.h"
#undef NTUPLE_INTERFACE_OPEN

unsigned long int n_entries = NT_output_ttree->GetEntries();
for (unsigned long int iev = 0; iev < n_entries; iev++)
	{
	NT_output_ttree->GetEntry(iev);

	if (run_test && iev > 10) break;

	//// tests
	//cout << iev << ": " << NT_leps_ID ; //<< " " << NT_jet_p4->size();
	//for (unsigned int i=0; i<NT_lep_relIso->size(); i++)
	//	{
	//	cout << " " << ((*NT_lep_relIso)[i]);
	//	}
	//cout << " jets ";
	//for (unsigned int i=0; i<NT_jet_p4->size(); i++)
	//	{
	//	cout << " " << ((*NT_jet_initial_p4)[i].pt());
	//	}
	//cout << " objects ";
	//cout << " " << ((*NT_met_init).pt()) << endl;
	//cout << endl;

	// loop over jets in the event
	for (unsigned int i=0; i<NT_jet_p4->size(); i++)
		{
		//cout << " " << ((*NT_jet_initial_p4)[i].pt());
		auto p4   = (*NT_jet_initial_p4)[i];

		// jet conditions
		auto pfid = (*NT_jet_PFID)[i];
		if (pfid < 1 or abs(p4.eta()) > 2.5 or p4.pt() < 30.) continue; // Loose PFID and eta

		if ((*NT_jet_matching_lep)[i]) continue;
		// ignoring match to tau

		auto jet_b_discr = (*NT_jet_b_discr)[i];
		bool b_tagged_medium = jet_b_discr > b_tag_wp_medium;
		bool b_tagged = b_tagged_medium;
		auto HF = (*NT_jet_hadronFlavour)[i];

		if (HF == 5)
			{
			btag_b_hadronFlavour_candidates     ->Fill(p4.pt(), p4.eta());
			if (b_tagged) btag_b_hadronFlavour_candidates_tagged->Fill(p4.pt(), p4.eta());
			}

		else if (HF == 4)
			{
			btag_c_hadronFlavour_candidates     ->Fill(p4.pt(), p4.eta());
			if (b_tagged) btag_c_hadronFlavour_candidates_tagged ->Fill(p4.pt(), p4.eta());
			}

		else
			{
			btag_udsg_hadronFlavour_candidates  ->Fill(p4.pt(), p4.eta());
			if (b_tagged) btag_udsg_hadronFlavour_candidates_tagged ->Fill(p4.pt(), p4.eta());
			}

		// gen-match info for JES corrections
		if ((*NT_genjet_matched)[i])
			{
			auto pt_gen  = (*NT_genjet_pt)[i];
			auto pt_reco = p4.pt();

			if (abs(p4.eta()) < 1.5)
			// barrel
				{
				if (HF == 5)
					{
					jes_cor_2d_gen_reco_barrel_b->Fill(pt_gen, pt_reco);
					jes_cor_1d_gen_reco_barrel_b->Fill(pt_gen, pt_reco/pt_gen);
					jes_cor_1d_gen_reco_barrel_b_unit->Fill(pt_gen);
					}
				else if (HF == 4)
					{
					jes_cor_2d_gen_reco_barrel_c->Fill(pt_gen, pt_reco);
					jes_cor_1d_gen_reco_barrel_c->Fill(pt_gen, pt_reco/pt_gen);
					jes_cor_1d_gen_reco_barrel_c_unit->Fill(pt_gen);
					}
				else
					{
					jes_cor_2d_gen_reco_barrel_udsg->Fill(pt_gen, pt_reco);
					jes_cor_1d_gen_reco_barrel_udsg->Fill(pt_gen, pt_reco/pt_gen);
					jes_cor_1d_gen_reco_barrel_udsg_unit->Fill(pt_gen);
					}
				}

			else
				{
				// endcap
				if (HF == 5)
					{
					jes_cor_2d_gen_reco_endcap_b->Fill(pt_gen, pt_reco);
					jes_cor_1d_gen_reco_endcap_b->Fill(pt_gen, pt_reco/pt_gen);
					jes_cor_1d_gen_reco_endcap_b_unit->Fill(pt_gen);
					}
				else if (HF == 4)
					{
					jes_cor_2d_gen_reco_endcap_c->Fill(pt_gen, pt_reco);
					jes_cor_1d_gen_reco_endcap_c->Fill(pt_gen, pt_reco/pt_gen);
					jes_cor_1d_gen_reco_endcap_c_unit->Fill(pt_gen);
					}
				else
					{
					jes_cor_2d_gen_reco_endcap_udsg->Fill(pt_gen, pt_reco);
					jes_cor_1d_gen_reco_endcap_udsg->Fill(pt_gen, pt_reco/pt_gen);
					jes_cor_1d_gen_reco_endcap_udsg_unit->Fill(pt_gen);
					}
				}
			}

		}
	}

tfile->Close();
}


TFile* fout = new TFile(output_file, "RECREATE");
fout->Write();

//for name, histo in histos.items():
//    histo.Write()

btag_b_hadronFlavour_candidates ->Write();
btag_b_hadronFlavour_candidates_tagged ->Write();
btag_c_hadronFlavour_candidates ->Write();
btag_c_hadronFlavour_candidates_tagged ->Write();
btag_udsg_hadronFlavour_candidates ->Write();
btag_udsg_hadronFlavour_candidates_tagged ->Write();

jes_cor_2d_gen_reco_barrel_b ->Write();
jes_cor_2d_gen_reco_endcap_b ->Write();

jes_cor_2d_gen_reco_barrel_c ->Write();
jes_cor_2d_gen_reco_endcap_c ->Write();

jes_cor_2d_gen_reco_barrel_udsg ->Write();
jes_cor_2d_gen_reco_endcap_udsg ->Write();

jes_cor_1d_gen_reco_barrel_b ->Write();
jes_cor_1d_gen_reco_endcap_b ->Write();
jes_cor_1d_gen_reco_barrel_b_unit ->Write();
jes_cor_1d_gen_reco_endcap_b_unit ->Write();

jes_cor_1d_gen_reco_barrel_c ->Write();
jes_cor_1d_gen_reco_endcap_c ->Write();
jes_cor_1d_gen_reco_barrel_c_unit ->Write();
jes_cor_1d_gen_reco_endcap_c_unit ->Write();

jes_cor_1d_gen_reco_barrel_udsg ->Write();
jes_cor_1d_gen_reco_endcap_udsg ->Write();
jes_cor_1d_gen_reco_barrel_udsg_unit ->Write();
jes_cor_1d_gen_reco_endcap_udsg_unit ->Write();

fout->Write();
fout->Close();

return 0;
}

