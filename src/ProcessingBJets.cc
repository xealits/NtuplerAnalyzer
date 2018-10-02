#include "UserCode/NtuplerAnalyzer/interface/ProcessingBJets.h"
//#include "UserCode/NtuplerAnalyzer/interface/recordFuncs.h"


/* defined in header
struct bTaggingEfficiencyHistograms {
	TH2F* b_alljet   ;
	TH2F* b_tagged   ;
	TH2F* c_alljet   ;
	TH2F* c_tagged   ;
	TH2F* udsg_alljet;
	TH2F* udsg_tagged;
	};
*/

double bTagging_b_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta)
	{
	Int_t global_bin_id = bEffs.b_alljet->FindBin(pt, eta); // the binning should be the same in both histograms
        Double_t N_alljets  = bEffs.b_alljet->GetBinContent(global_bin_id);
	Double_t N_tagged   = bEffs.b_tagged->GetBinContent(global_bin_id);
	if (N_alljets < 0.001 || N_tagged < 0.001)
		{
		//fill_2d(string("btag_eff_retrieved0_flavour_b"), 250, 0., 500., 200, -4., 4., pt, eta, 1);
		//fill_btag_efficiency(string("btag_eff_retrieved0_flavour_b"), pt, eta, 1);
		return 0; // whatch out -- equality of floats
		}

	return N_tagged/N_alljets;
	}

double bTagging_c_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta)
	{
	Int_t global_bin_id = bEffs.c_alljet->FindBin(pt, eta); // the binning should be the same in both histograms
        Double_t N_alljets  = bEffs.c_alljet->GetBinContent(global_bin_id);
	Double_t N_tagged   = bEffs.c_tagged->GetBinContent(global_bin_id);
	if (N_alljets < 0.001 || N_tagged < 0.001)
		{
		//fill_btag_efficiency(string("btag_eff_retrieved0_flavour_c"), pt, eta, 1);
		return 0;
		}
	return N_tagged/N_alljets;
	}

double bTagging_udsg_jet_efficiency(struct bTaggingEfficiencyHistograms& bEffs, double& pt, double& eta)
	{
	Int_t global_bin_id = bEffs.udsg_alljet->FindBin(pt, eta); // the binning should be the same in both histograms
        Double_t N_alljets  = bEffs.udsg_alljet->GetBinContent(global_bin_id);
	Double_t N_tagged   = bEffs.udsg_tagged->GetBinContent(global_bin_id);
	if (N_alljets < 0.001 || N_tagged < 0.001)
		{
		//fill_btag_efficiency(string("btag_eff_retrieved0_flavour_udsg"), pt, eta, 1);
		return 0;
		}
	return N_tagged/N_alljets;
	}




int processBJets_BTag(pat::JetCollection& jets, bool isMC, double& weight, double& bTaggingSF_eventWeight, // input
	BTagCalibrationReader& btagCal, // BTagSFUtil& btsfutil, old b-tag SF weighting, done with bEffs now
	struct bTaggingEfficiencyHistograms& bEffs,
	string& b_tagger_label, float b_tag_WP,
	pat::JetCollection& selBJets,                          // output
	bool record, bool debug) // more output

{
// for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)

for (size_t ijet = 0; ijet < jets.size(); ++ijet)
	{
	pat::Jet& jet = jets[ijet];

	double eta=jet.eta();
	double pt=jet.pt();

	float b_discriminator = jet.bDiscriminator(b_tagger_label);
	//fill_1d(string("btag_discriminator"), 200, -1.0, 1.0, b_discriminator, weight);

	bool hasCSVtag(b_discriminator > b_tag_WP);
	//bool raw_CSV_tag = hasCSVtag;
	//bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.935);
	bool hasCSVtag_BTagUp(false), hasCSVtag_BTagDown(false);

	/*
	if (record)
		{
		fill_2d(string("btag_b_tagging_candidate_jets_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
		fill_1d(string("btag_b_tagging_candidate_jets_phi"), 128, -3.2, 3.2, jet.phi(), weight);
		if (hasCSVtag) fill_2d(string("btag_b_tagging_candidate_jets_pt_eta_tagged"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
		}
	*/

	//update according to the SF measured by BTV
	// new fency procedure with CSV files
	// 80X calibrators in btagCal
	// usage:
	// btagCal.eval_auto_bounds("central", BTagEntry::FLAV_B, b_jet.eta(), b_jet.pt());
	// -- one calibrator for central value, and up/down
	// thus the specification of the value to callibrate,
	// instead of different callibrators
	if(isMC){
		// int flavId=jet.partonFlavour();
		int flavId=jet.hadronFlavour();
		// also: patJet->genParton().pdgId()
		// fill_btag_eff(string("mc_all_b_tagging_candidate_jets_pt_eta"), jet.pt(), eta, weight);

		double sf = 1.0, eff = 1.0;
		/* If the jet is tagged -- weight *= SF of the jet
		 * if not weight *= (1 - eff*SF)/(1 - eff)
		 */

		if (abs(flavId)==5) {
			// get SF for the jet
			sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_B, eta, pt, 0.);

			// get eff for the jet
			eff = bTagging_b_jet_efficiency(bEffs, pt, eta);
			}
		else if(abs(flavId)==4) {
			sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_C, eta, pt, 0.);

			eff = bTagging_c_jet_efficiency(bEffs, pt, eta);
			}
		else {
			sf = btagCal.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, eta, pt, 0.);

			eff = bTagging_udsg_jet_efficiency(bEffs, pt, eta);
			}

		double jet_weight_factor = 1;
		if (hasCSVtag) // a tagged jet
			{
			jet_weight_factor = sf;
			}
		else // not tagged
			{
			// truncate efficiency to 0 and 0.99
			eff = (eff < 0. ? 0. : (eff > 0.99 ? 0.99 : eff));
			jet_weight_factor = (1 - sf*eff) / (1 - eff);
			}

		bTaggingSF_eventWeight *= jet_weight_factor;
		}

	// now record the btagging discriminators with the reweighted weight
	double full_weight = weight * bTaggingSF_eventWeight;
	if (isMC)
		{
		int flavId=jet.hadronFlavour();
		}

	if(hasCSVtag || hasCSVtag_BTagUp || hasCSVtag_BTagDown)
		{
		selBJets.push_back(jet);
		}
	}


return 0;
}



int processBJets_BTag_for_sys(pat::JetCollection& jets, bool isMC, double& weight, double& bTaggingSF_eventWeight, // input
	BTagCalibrationReader& btagCal, // BTagSFUtil& btsfutil, old b-tag SF weighting, done with bEffs now
	struct bTaggingEfficiencyHistograms& bEffs,
	string& b_tagger_label, float b_tag_WP,
	pat::JetCollection& selBJets,                          // output
	const string& btag_sys_point,
	bool record, bool record_control, bool debug) // more output

{
// for (size_t ijet = 0; ijet < selJetsNoLep.size(); ++ijet)

for (size_t ijet = 0; ijet < jets.size(); ++ijet)
	{
	pat::Jet& jet = jets[ijet];

	double eta=jet.eta();
	double pt=jet.pt();

	float b_discriminator = jet.bDiscriminator(b_tagger_label);
	if (record_control) fill_1d(string("btag_discriminator"), 200, -1.0, 1.0, b_discriminator, weight);

	bool hasCSVtag(b_discriminator > b_tag_WP);
	//bool raw_CSV_tag = hasCSVtag;
	//bool hasCSVtag(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.935);
	bool hasCSVtag_BTagUp(false), hasCSVtag_BTagDown(false);

	if (record)
		{
		fill_2d(string("btag_b_tagging_candidate_jets_pt_eta"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
		fill_1d(string("btag_b_tagging_candidate_jets_phi"), 128, -3.2, 3.2, jet.phi(), weight);
		if (hasCSVtag) fill_2d(string("btag_b_tagging_candidate_jets_pt_eta_tagged"), 250, 0., 500., 200, -4., 4., jet.pt(), eta, weight);
		}

	//update according to the SF measured by BTV
	// new fency procedure with CSV files
	// 80X calibrators in btagCal
	// usage:
	// btagCal.eval_auto_bounds("central", BTagEntry::FLAV_B, b_jet.eta(), b_jet.pt());
	// -- one calibrator for central value, and up/down
	// thus the specification of the value to callibrate,
	// instead of different callibrators
	if(isMC){
		// int flavId=jet.partonFlavour();
		int flavId=jet.hadronFlavour();
		// also: patJet->genParton().pdgId()
		// fill_btag_eff(string("mc_all_b_tagging_candidate_jets_pt_eta"), jet.pt(), eta, weight);

		double sf = 1.0, eff = 1.0;
		/* If the jet is tagged -- weight *= SF of the jet
		 * if not weight *= (1 - eff*SF)/(1 - eff)
		 */

		if (abs(flavId)==5) {
			if (record)
				{
				/* recording the jets, per-flavour, candidates and tagged --- for the efficiency measurement and for control
				 * NOTICE: the recorded jets are weighted with the current genWeight, pile-up etc
				 */
				//fill_2d(string("btag_b_hadronFlavour_candidates"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				fill_1d(string("btag_discriminator_b_hadronFlavour"), 200, -1.0, 1.0, b_discriminator, weight);
				fill_btag_efficiency(string("btag_b_hadronFlavour_candidates"), pt, eta, weight); // shouldn't weight be = 1?
				if (hasCSVtag)
					fill_btag_efficiency(string("btag_b_hadronFlavour_candidates_tagged"), pt, eta, weight); // shouldn't weight be = 1?
					//fill_2d(string("btag_b_hadronFlavour_candidates_tagged"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				}

			// get SF for the jet
			sf = btagCal.eval_auto_bounds(btag_sys_point, BTagEntry::FLAV_B, eta, pt, 0.);
			if (record) fill_1d(string("btag_sf_flavour_b"), 100, 0., 2.,   sf, weight);

			// get eff for the jet
			eff = bTagging_b_jet_efficiency(bEffs, pt, eta);
			if (record) fill_1d(string("btag_eff_flavour_b"), 100, 0., 2.,   eff, weight);
			}
		else if(abs(flavId)==4) {

			if (record)
				{ 
				//fill_2d(string("btag_c_hadronFlavour_candidates"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				fill_1d(string("btag_discriminator_c_hadronFlavour"), 200, -1.0, 1.0, b_discriminator, weight);
				fill_btag_efficiency(string("btag_c_hadronFlavour_candidates"), pt, eta, weight);
				if (hasCSVtag)
					fill_btag_efficiency(string("btag_c_hadronFlavour_candidates_tagged"), pt, eta, weight);
					//fill_2d(string("btag_c_hadronFlavour_candidates_tagged"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				}

			sf = btagCal.eval_auto_bounds(btag_sys_point, BTagEntry::FLAV_C, eta, pt, 0.);
			if (record) fill_1d(string("btag_sf_flavour_c"), 100, 0., 2.,   sf, weight);

			eff = bTagging_c_jet_efficiency(bEffs, pt, eta);
			if (record) fill_1d(string("btag_eff_flavour_c"), 100, 0., 2.,   eff, weight);
			}
		else {

			if (record)
				{ 
				//fill_2d(string("btag_udsg_hadronFlavour_candidates"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				fill_1d(string("btag_discriminator_udsg_hadronFlavour"), 200, -1.0, 1.0, b_discriminator, weight);
				fill_btag_efficiency(string("btag_udsg_hadronFlavour_candidates"), pt, eta, weight);
				if (hasCSVtag)
					fill_btag_efficiency(string("btag_udsg_hadronFlavour_candidates_tagged"), pt, eta, weight);
					//fill_2d(string("btag_udsg_hadronFlavour_candidates_tagged"), 250, 0., 500., 200, -4., 4., pt, eta, weight);
				}

			sf = btagCal.eval_auto_bounds(btag_sys_point, BTagEntry::FLAV_UDSG, eta, pt, 0.);
			if (record) fill_1d(string("btag_sf_flavour_udsg"), 100, 0., 2.,   sf, weight);

			eff = bTagging_udsg_jet_efficiency(bEffs, pt, eta);
			if (record) fill_1d(string("btag_eff_flavour_udsg"), 100, 0., 2.,   eff, weight);
			}

		double jet_weight_factor = 1;
		if (hasCSVtag) // a tagged jet
			{
			jet_weight_factor = sf;
			if (record_control)
				{
				if (abs(flavId)==5)
					fill_1d(string("btag_jetweight_flavour_b_tagged"), 100, 0., 2., jet_weight_factor, 1);
				else if (abs(flavId)==4)
					fill_1d(string("btag_jetweight_flavour_c_tagged"), 100, 0., 2., jet_weight_factor, 1);
				else
					fill_1d(string("btag_jetweight_flavour_udsg_tagged"), 100, 0., 2., jet_weight_factor, 1);
				}
			}
		else // not tagged
			{
			// truncate efficiency to 0 and 0.99
			eff = (eff < 0. ? 0. : (eff > 0.99 ? 0.99 : eff));
			jet_weight_factor = (1 - sf*eff) / (1 - eff);
			if (record_control)
				{
				if (abs(flavId)==5)
					fill_1d(string("btag_jetweight_flavour_b_notag"), 100, 0., 2., jet_weight_factor, 1);
				else if (abs(flavId)==4)
					fill_1d(string("btag_jetweight_flavour_c_notag"), 100, 0., 2., jet_weight_factor, 1);
				else
					fill_1d(string("btag_jetweight_flavour_udsg_notag"), 100, 0., 2., jet_weight_factor, 1);
				}
			}

		bTaggingSF_eventWeight *= jet_weight_factor;
		}

	// now record the btagging discriminators with the reweighted weight
	double full_weight = weight * bTaggingSF_eventWeight;
	if (record_control) fill_1d(string("btag_discriminator_weighted"), 200, -1.0, 1.0, b_discriminator, full_weight);
	if (isMC && record_control)
		{
		int flavId=jet.hadronFlavour();
		if (abs(flavId)==5)
			fill_1d(string("btag_discriminator_b_hadronFlavour_weighted"), 200, -1.0, 1.0, b_discriminator, full_weight);
		else if (abs(flavId)==4)
			fill_1d(string("btag_discriminator_c_hadronFlavour_weighted"), 200, -1.0, 1.0, b_discriminator, full_weight);
		else
			fill_1d(string("btag_discriminator_udsg_hadronFlavour_weighted"), 200, -1.0, 1.0, b_discriminator, full_weight);
		}

	if(hasCSVtag || hasCSVtag_BTagUp || hasCSVtag_BTagDown)
		{
		selBJets.push_back(jet);
		if (record)
			{
			fill_2d(string("btag_control_jet_selBJets_pt_eta"), 250, 0., 500., 200, -4., 4., pt, jet.eta(), weight);
			fill_1d(string("btag_control_jet_selBJets_phi"), 128, -3.2, 3.2, jet.phi(), weight);
			}
		}
	}


return 0;
}


