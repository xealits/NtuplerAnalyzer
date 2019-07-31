#include <iostream>
#include <stdio.h>

void accumulate_b_effs (TString acum_file, TString new_file)
{

	/*
	char *filename = "./qcd_EM_Pt30to80_0.root";
	TFile *f = TFile::Open( filename );

	TTree *t = (TTree *) f->Get( "myEventSelector/data" );
	MiniEvent_t ev;
	*/

	TFile *f_new  = (TFile*) new TFile( new_file  );
	TFile *f_acum = (TFile*) new TFile( acum_file, "UPDATE" );

	TString distr_name("btag_udsg_hadronFlavour_candidates_tagged");
	if (!f_new->GetListOfKeys()->Contains(distr_name))
		{
		cerr << distr_name << " is not found in " << new_file << endl;
		return;
		}

	TH2D * b_candidates    = (TH2D*) f_new->Get("btag_b_hadronFlavour_candidates");
	TH2D * b_tagged        = (TH2D*) f_new->Get("btag_b_hadronFlavour_candidates_tagged");
	TH2D * c_candidates    = (TH2D*) f_new->Get("btag_c_hadronFlavour_candidates");
	TH2D * c_tagged        = (TH2D*) f_new->Get("btag_c_hadronFlavour_candidates_tagged");
	TH2D * udsg_candidates = (TH2D*) f_new->Get("btag_udsg_hadronFlavour_candidates");
	TH2D * udsg_tagged     = (TH2D*) f_new->Get("btag_udsg_hadronFlavour_candidates_tagged");

	b_candidates    -> SetName(new_file + "_btag_b_hadronFlavour_candidates");
	b_tagged        -> SetName(new_file + "_btag_b_hadronFlavour_candidates_tagged");
	c_candidates    -> SetName(new_file + "_btag_c_hadronFlavour_candidates");
	c_tagged        -> SetName(new_file + "_btag_c_hadronFlavour_candidates_tagged");
	udsg_candidates -> SetName(new_file + "_btag_udsg_hadronFlavour_candidates");
	udsg_tagged     -> SetName(new_file + "_btag_udsg_hadronFlavour_candidates_tagged");

	f_acum->Write();

	b_candidates    -> Write();
	b_tagged        -> Write();
	c_candidates    -> Write();
	c_tagged        -> Write();
	udsg_candidates -> Write();
	udsg_tagged     -> Write();

	f_acum->Close();
	f_new->Close();
}

