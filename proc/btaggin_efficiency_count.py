import argparse
import logging
import ctypes
from os.path import isfile, basename


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "count the jets for b tagging efficiency histogram in pT-eta",
    epilog = """Example:\npython btag_eff_test.root gstore_outdirs/v31v27/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v31v27_MC2016_Summer16_TTJets_powheg/180522_225221/0000/*root"""
    )

parser.add_argument("output_file", type=str, default="output.root", help="filename for output")
parser.add_argument("--histo-prefix", type=str, default="", help="to match the old naming of efficiency histos: start with dtag_")

parser.add_argument("--test",   action='store_true', help="test run on 10 events")
parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")

parser.add_argument('input_files', nargs='+', help="""the files to sum up, passed by shell, as:

gstore_outdirs/v31v27/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v31v27_MC2016_Summer16_TTJets_powheg/180522_225221/0000/*root""")


args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)


logging.info("import ROOT")

import ROOT
from ROOT import TH1D, TH2D, TFile, TTree, gROOT
from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
from ROOT import TLegend
#from ROOT import THStack

logging.info("done")

from plotting_root import rgb

gROOT.SetBatch()


assert not isfile(args.output_file)

# the counting histos for bSF efficiencies

bin_edges_pt  = [ 0, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 65, 70, 75, 80, 100, 150, 500 ]
bin_edges_eta = [ -2.5, -2.4, -2.3, -2.2, -2.1, -2., -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5 ]

def make_root_histo_custom_bins(bin_edges):
    #bin_edges = [float(b) for b in args.custom_range.split(',')]
    n_bin_edges = len(bin_edges)
    n_bins = n_bin_edges - 1
    logging.debug("making %d bins" % n_bins)

    # first method
    root_bin_edges = (ctypes.c_double * n_bin_edges)(*bin_edges)
    # use as:
    #temp_output_histo = TH1D("histotemp", "", n_bins, root_bin_edges) # root commands can access it by the name
    return (n_bins, root_bin_edges)

n_bin_pt,  root_bin_pt  = make_root_histo_custom_bins(bin_edges_pt)
n_bin_eta, root_bin_eta = make_root_histo_custom_bins(bin_edges_eta)

histos = {}
def histo(*inparams):
    TH_type = inparams[0]
    in_name = inparams[1]
    name = args.histo_prefix + in_name
    the_params = inparams[2:]
    histos[in_name] = TH_type(*((name,) + the_params))

histo(TH2D, "btag_b_hadronFlavour_candidates"           , "", n_bin_pt,  root_bin_pt, n_bin_eta, root_bin_eta)
histo(TH2D, "btag_b_hadronFlavour_candidates_tagged"    , "", n_bin_pt,  root_bin_pt, n_bin_eta, root_bin_eta)
histo(TH2D, "btag_c_hadronFlavour_candidates"           , "", n_bin_pt,  root_bin_pt, n_bin_eta, root_bin_eta)
histo(TH2D, "btag_c_hadronFlavour_candidates_tagged"    , "", n_bin_pt,  root_bin_pt, n_bin_eta, root_bin_eta)
histo(TH2D, "btag_udsg_hadronFlavour_candidates"        , "", n_bin_pt,  root_bin_pt, n_bin_eta, root_bin_eta)
histo(TH2D, "btag_udsg_hadronFlavour_candidates_tagged" , "", n_bin_pt,  root_bin_pt, n_bin_eta, root_bin_eta)


histo(TH2D, "jes_cor_2d_gen_reco_barrel_b" , "", n_bin_pt,  root_bin_pt, n_bin_pt, root_bin_pt)
histo(TH2D, "jes_cor_2d_gen_reco_endcap_b" , "", n_bin_pt,  root_bin_pt, n_bin_pt, root_bin_pt)

histo(TH2D, "jes_cor_2d_gen_reco_barrel_c" , "", n_bin_pt,  root_bin_pt, n_bin_pt, root_bin_pt)
histo(TH2D, "jes_cor_2d_gen_reco_endcap_c" , "", n_bin_pt,  root_bin_pt, n_bin_pt, root_bin_pt)

histo(TH2D, "jes_cor_2d_gen_reco_barrel_udsg" , "", n_bin_pt,  root_bin_pt, n_bin_pt, root_bin_pt)
histo(TH2D, "jes_cor_2d_gen_reco_endcap_udsg" , "", n_bin_pt,  root_bin_pt, n_bin_pt, root_bin_pt)

histo(TH1D, "jes_cor_1d_gen_reco_barrel_b" , "", n_bin_pt,  root_bin_pt)
histo(TH1D, "jes_cor_1d_gen_reco_endcap_b" , "", n_bin_pt,  root_bin_pt)
histo(TH1D, "jes_cor_1d_gen_reco_barrel_b_unit" , "", n_bin_pt,  root_bin_pt)
histo(TH1D, "jes_cor_1d_gen_reco_endcap_b_unit" , "", n_bin_pt,  root_bin_pt)

histo(TH1D, "jes_cor_1d_gen_reco_barrel_c" , "", n_bin_pt,  root_bin_pt)
histo(TH1D, "jes_cor_1d_gen_reco_endcap_c" , "", n_bin_pt,  root_bin_pt)
histo(TH1D, "jes_cor_1d_gen_reco_barrel_c_unit" , "", n_bin_pt,  root_bin_pt)
histo(TH1D, "jes_cor_1d_gen_reco_endcap_c_unit" , "", n_bin_pt,  root_bin_pt)

histo(TH1D, "jes_cor_1d_gen_reco_barrel_udsg" , "", n_bin_pt,  root_bin_pt)
histo(TH1D, "jes_cor_1d_gen_reco_endcap_udsg" , "", n_bin_pt,  root_bin_pt)
histo(TH1D, "jes_cor_1d_gen_reco_barrel_udsg_unit" , "", n_bin_pt,  root_bin_pt)
histo(TH1D, "jes_cor_1d_gen_reco_endcap_udsg_unit" , "", n_bin_pt,  root_bin_pt)


# loop over all events in all files
# counting the non-lep matched jets
# per their hadron flavour
# b-tagged and all
b_tag_wp_medium = 0.8484 # Medium WP

for filename in args.input_files:
    if not isfile(filename):
        logging.info("missing: " + filename)
        continue

    logging.debug(filename)

    tfile = TFile(filename)
    ttree = tfile.Get("ntupler/reduced_ttree")

    for iev, ev in enumerate(ttree):
        if args.test and iev > 10: break
        # only elmu selection
        #pass_elmu_id = ev.leps_ID == -11*13 and ev.HLT_mu and ev.no_iso_veto_leps and \
        #    (ev.lep_matched_HLT[0] if abs(ev.lep_id[0]) == 13 else ev.lep_matched_HLT[1]) and \
        #    (ev.lep_p4[0].pt() > 30 and abs(ev.lep_p4[0].eta()) < 2.4) and \
        #    (ev.lep_p4[1].pt() > 30 and abs(ev.lep_p4[1].eta()) < 2.4)
        #if not pass_elmu_id: continue
        # all 
        # loop through jets of main selection

        for i in xrange(ev.jet_p4.size()):
            p4   = ev.jet_initial_p4[i]

            # jet conditions
            pfid = ev.jet_PFID[i]
            if pfid < 1 or abs(p4.eta()) > 2.5 or p4.pt() < 30.: continue # Loose PFID and eta

            if ev.jet_matching_lep[i]: continue
            # ignoring match to tau

            jet_b_discr = ev.jet_b_discr[i]
            b_tagged = b_tagged_medium = jet_b_discr > b_tag_wp_medium
            HF = ev.jet_hadronFlavour[i]

            if HF == 5:
                histos['btag_b_hadronFlavour_candidates']     .Fill(p4.pt(), p4.eta())
                if b_tagged: histos['btag_b_hadronFlavour_candidates_tagged'] .Fill(p4.pt(), p4.eta())

            elif HF == 4:
                histos['btag_c_hadronFlavour_candidates']     .Fill(p4.pt(), p4.eta())
                if b_tagged: histos['btag_c_hadronFlavour_candidates_tagged'] .Fill(p4.pt(), p4.eta())

            else:
                histos['btag_udsg_hadronFlavour_candidates']     .Fill(p4.pt(), p4.eta())
                if b_tagged: histos['btag_udsg_hadronFlavour_candidates_tagged'] .Fill(p4.pt(), p4.eta())

            # gen-match info for JES corrections
            if ev.genjet_matched[i]:
                pt_gen  = ev.genjet_pt[i]
                pt_reco = p4.pt()

                if abs(p4.eta()) < 1.5:
                    # barrel
                    if HF == 5:
                        histos["jes_cor_2d_gen_reco_barrel_b"].Fill(pt_gen, pt_reco)
                        histos["jes_cor_1d_gen_reco_barrel_b"].Fill(pt_gen, pt_reco/pt_gen)
                        histos["jes_cor_1d_gen_reco_barrel_b_unit"].Fill(pt_gen)
                    elif HF == 4:
                        histos["jes_cor_2d_gen_reco_barrel_c"].Fill(pt_gen, pt_reco)
                        histos["jes_cor_1d_gen_reco_barrel_c"].Fill(pt_gen, pt_reco/pt_gen)
                        histos["jes_cor_1d_gen_reco_barrel_c_unit"].Fill(pt_gen)
                    else:
                        histos["jes_cor_2d_gen_reco_barrel_udsg"].Fill(pt_gen, pt_reco)
                        histos["jes_cor_1d_gen_reco_barrel_udsg"].Fill(pt_gen, pt_reco/pt_gen)
                        histos["jes_cor_1d_gen_reco_barrel_udsg_unit"].Fill(pt_gen)

                else:
                    # endcap
                    if HF == 5:
                        histos["jes_cor_2d_gen_reco_endcap_b"].Fill(pt_gen, pt_reco)
                        histos["jes_cor_1d_gen_reco_endcap_b"].Fill(pt_gen, pt_reco/pt_gen)
                        histos["jes_cor_1d_gen_reco_endcap_b_unit"].Fill(pt_gen)
                    elif HF == 4:
                        histos["jes_cor_2d_gen_reco_endcap_c"].Fill(pt_gen, pt_reco)
                        histos["jes_cor_1d_gen_reco_endcap_c"].Fill(pt_gen, pt_reco/pt_gen)
                        histos["jes_cor_1d_gen_reco_endcap_c_unit"].Fill(pt_gen)
                    else:
                        histos["jes_cor_2d_gen_reco_endcap_udsg"].Fill(pt_gen, pt_reco)
                        histos["jes_cor_1d_gen_reco_endcap_udsg"].Fill(pt_gen, pt_reco/pt_gen)
                        histos["jes_cor_1d_gen_reco_endcap_udsg_unit"].Fill(pt_gen)

    tfile.Close()

fout = TFile(args.output_file, "RECREATE")
fout.Write()

for name, histo in histos.items():
    histo.Write()

fout.Write()
fout.Close()

