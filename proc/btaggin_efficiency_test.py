import argparse
import logging
import ctypes
from os.path import isfile, basename




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "sumup TTree Draw",
    epilog = """Example:\npython btaggin_efficiency_test.py btag_eff_test.root btag_effs_elmu_test1.root gstore_outdirs/v31v27/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v31v27_MC2016_Summer16_TTJets_powheg/180522_225221/0000/*root"""
    )

parser.add_argument("output_file", type=str, default="output.root", help="filename for output")
parser.add_argument("eff_file",    type=str, default="output.root", help="file with b effs")

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



input_files = args.input_files

assert isfile(args.eff_file)
eff_f = TFile(args.eff_file)

efficiency_pt_eta_b        = eff_f.Get("MC2016_Summer16_TTJets_powheg_btag_b_hadronFlavour_candidates")
efficiency_pt_eta_b_tag    = eff_f.Get("MC2016_Summer16_TTJets_powheg_btag_b_hadronFlavour_candidates_tagged")
efficiency_pt_eta_c        = eff_f.Get("MC2016_Summer16_TTJets_powheg_btag_c_hadronFlavour_candidates")
efficiency_pt_eta_c_tag    = eff_f.Get("MC2016_Summer16_TTJets_powheg_btag_c_hadronFlavour_candidates_tagged")
efficiency_pt_eta_udsg     = eff_f.Get("MC2016_Summer16_TTJets_powheg_btag_udsg_hadronFlavour_candidates")
efficiency_pt_eta_udsg_tag = eff_f.Get("MC2016_Summer16_TTJets_powheg_btag_udsg_hadronFlavour_candidates_tagged")

efficiency_pt_eta_b_tag        .Divide( efficiency_pt_eta_b    )
efficiency_pt_eta_c_tag        .Divide( efficiency_pt_eta_c    )
efficiency_pt_eta_udsg_tag     .Divide( efficiency_pt_eta_udsg )

import support_btagging_sf
from support_btagging_sf import calc_btag_sf_weight

support_btagging_sf.bEff_histo_b    = efficiency_pt_eta_b_tag
support_btagging_sf.bEff_histo_c    = efficiency_pt_eta_c_tag
support_btagging_sf.bEff_histo_udsg = efficiency_pt_eta_udsg_tag

b_tag_wp_medium = 0.8 #0.8484

histos = {}
def histo(*args):
    histos[args[0]] = TH1D(*args)

histo("bMjet_pt_events",  "", 50, 0., 200.)
histo("bMjet_pt_njets",   "", 50, 0., 200.)
histo("bMjet_pt_nbtags",  "", 50, 0., 200.)
histo("bMjet_pt_flav_nb",    "", 50, 0., 200.)
histo("bMjet_pt_flav_nc",    "", 50, 0., 200.)
histo("bMjet_pt_flav_nudsg", "", 50, 0., 200.)

histo("bMjet_pt_flav_sf_b",    "", 50, 0., 200.)
histo("bMjet_pt_flav_sf_c",    "", 50, 0., 200.)
histo("bMjet_pt_flav_sf_udsg", "", 50, 0., 200.)

histo("bMjet_pt_flav_eff_b",    "", 50, 0., 200.)
histo("bMjet_pt_flav_eff_c",    "", 50, 0., 200.)
histo("bMjet_pt_flav_eff_udsg", "", 50, 0., 200.)

histo("bMjet_pt_flav_weight_b",    "", 50, 0., 200.)
histo("bMjet_pt_flav_weight_c",    "", 50, 0., 200.)
histo("bMjet_pt_flav_weight_udsg", "", 50, 0., 200.)

histo("bMjet_pt_bSF_weight", "", 50, 0., 200.)
histo("bMjet_pt_bSF_b"     , "", 50, 0., 200.)
histo("bMjet_pt_bSF_c"     , "", 50, 0., 200.)
histo("bMjet_pt_bSF_udsg"  , "", 50, 0., 200.)

histo("weight_all_simple",  "", 100, 0., 2.)
histo("weight_all",  "", 100, 0., 2.)
histo("weight_b",    "", 100, 0., 2.)
histo("weight_c",    "", 100, 0., 2.)
histo("weight_udsg", "", 100, 0., 2.)

histo("bMjet0_flavour", "", 10, 0., 10.)
histo("bMjet1_flavour", "", 10, 0., 10.)
histo("bRjet0_flavour", "", 10, 0., 10.)


# test b-tagging weights effect on bMjet pt in elmu close
# for this:
# - per bMjet pt N jets, N b-tagged jets and diff flavours
# - for different flavours the weight of the jet
# - weight of event (in simple elmu)
# - weight per bMjet pt and what it consists of?
#   weight = product of each jet weights
#   split in flavours
#   N jets, N each flavour, weight from each flavour
#   which one is the bMjet?
# effect of efficiency on it?
#   need for each flavour per bMjet pt SF and eff
for filename in input_files:
    if not isfile(filename):
        logging.info("missing: " + filename)
        continue

    logging.debug(filename)

    tfile = TFile(filename)
    ttree = tfile.Get("ntupler/reduced_ttree")

    for iev, ev in enumerate(ttree):
        #if iev > 10: break
        # only elmu selection
        pass_elmu_id = ev.leps_ID == -11*13 and ev.HLT_mu and ev.no_iso_veto_leps and \
            (ev.lep_matched_HLT[0] if abs(ev.lep_id[0]) == 13 else ev.lep_matched_HLT[1]) and \
            (ev.lep_p4[0].pt() > 30 and abs(ev.lep_p4[0].eta()) < 2.4) and \
            (ev.lep_p4[1].pt() > 30 and abs(ev.lep_p4[1].eta()) < 2.4)

        if not pass_elmu_id: continue

        n_bjets = 0
        n_jets  = 0
        n_b_jets    = 0
        n_c_jets    = 0
        n_udsg_jets = 0

        flavours_tagged     = []
        flavours_tagged_not = []

        btagged_pts     = []

        bweight_all  = 1.
        bweight_b    = 1.
        bweight_c    = 1.
        bweight_udsg = 1.

        flav_sum_weight_b    = 0.
        flav_sum_weight_c    = 0.
        flav_sum_weight_udsg = 0.

        flav_sum_sf_b    = 0.
        flav_sum_sf_c    = 0.
        flav_sum_sf_udsg = 0.

        flav_sum_eff_b    = 0.
        flav_sum_eff_c    = 0.
        flav_sum_eff_udsg = 0.

        # loop through jets of main selection
        for i in xrange(ev.jet_p4.size()):

            p4   = ev.jet_initial_p4[i]
            pfid = ev.jet_PFID[i]
            if pfid < 1 or abs(p4.eta()) > 2.5 or p4.pt() < 30.: continue # Loose PFID and eta

            if ev.jet_matching_lep[i]: continue
            # ignoring match to tau

            n_jets  += 1

            jet_b_discr = ev.jet_b_discr[i]
            b_tagged = b_tagged_medium = jet_b_discr > b_tag_wp_medium
            HF = ev.jet_hadronFlavour[i]

            jet_weight_bSF, b_sf, b_eff = calc_btag_sf_weight(b_tagged_medium, HF, p4.pt(), p4.eta())

            if b_tagged:
                n_bjets += 1
                flavours_tagged.append(HF)
                btagged_pts.append(p4.pt())
            else:
                flavours_tagged_not.append(HF)

            bweight_all *= jet_weight_bSF

            if HF == 5:
                bweight_b *= jet_weight_bSF
                flav_sum_weight_b += jet_weight_bSF
                flav_sum_sf_b     += b_sf
                flav_sum_eff_b    += b_eff
                n_b_jets += 1
            elif HF == 4:
                bweight_c *= jet_weight_bSF
                flav_sum_weight_c += jet_weight_bSF
                flav_sum_sf_c     += b_sf
                flav_sum_eff_c    += b_eff
                n_c_jets += 1
            else:
                bweight_udsg *= jet_weight_bSF
                flav_sum_weight_udsg += jet_weight_bSF
                flav_sum_sf_udsg     += b_sf
                flav_sum_eff_udsg    += b_eff
                n_udsg_jets += 1

        histos['weight_all_simple'].Fill(bweight_all)

        if not (n_jets > 1 and n_bjets > 0): continue
        # elmu close

        histos['weight_all'].Fill(bweight_all)
        histos['weight_b'].Fill(bweight_b)
        histos['weight_c'].Fill(bweight_c)
        histos['weight_udsg'].Fill(bweight_udsg)

        histos["bMjet0_flavour"].Fill(flavours_tagged[0])
        if len(flavours_tagged) > 1:
            histos["bMjet1_flavour"].Fill(flavours_tagged[1])
        if len(flavours_tagged_not) > 0:
            histos["bRjet0_flavour"].Fill(flavours_tagged_not[0])

        bMjet_pt = btagged_pts[0]

        histos["bMjet_pt_events"] .Fill(bMjet_pt)
        histos["bMjet_pt_njets"]  .Fill(bMjet_pt, n_jets)
        histos["bMjet_pt_nbtags"] .Fill(bMjet_pt, n_bjets)

        histos["bMjet_pt_bSF_weight"] .Fill(bMjet_pt, bweight_all)
        histos["bMjet_pt_bSF_b"] .Fill(bMjet_pt, bweight_b)
        histos["bMjet_pt_bSF_c"] .Fill(bMjet_pt, bweight_c)
        histos["bMjet_pt_bSF_udsg"] .Fill(bMjet_pt, bweight_udsg)

        histos["bMjet_pt_flav_nb"]   .Fill(bMjet_pt, n_b_jets)
        histos["bMjet_pt_flav_nc"]   .Fill(bMjet_pt, n_c_jets)
        histos["bMjet_pt_flav_nudsg"].Fill(bMjet_pt, n_udsg_jets)

        histos["bMjet_pt_flav_sf_b"]   .Fill(bMjet_pt, flav_sum_sf_b)
        histos["bMjet_pt_flav_sf_c"]   .Fill(bMjet_pt, flav_sum_sf_c)
        histos["bMjet_pt_flav_sf_udsg"].Fill(bMjet_pt, flav_sum_sf_udsg)

        histos["bMjet_pt_flav_eff_b"]   .Fill(bMjet_pt, flav_sum_eff_b)
        histos["bMjet_pt_flav_eff_c"]   .Fill(bMjet_pt, flav_sum_eff_c)
        histos["bMjet_pt_flav_eff_udsg"].Fill(bMjet_pt, flav_sum_eff_udsg)

        histos["bMjet_pt_flav_weight_b"]   .Fill(bMjet_pt, flav_sum_weight_b)
        histos["bMjet_pt_flav_weight_c"]   .Fill(bMjet_pt, flav_sum_weight_c)
        histos["bMjet_pt_flav_weight_udsg"].Fill(bMjet_pt, flav_sum_weight_udsg)


    tfile.Close()

fout = TFile(args.output_file, "RECREATE")
fout.Write()

for name, histo in histos.items():
    histo.Write()

fout.Write()
fout.Close()

