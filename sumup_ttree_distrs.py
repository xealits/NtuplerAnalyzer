import argparse
import logging
from os.path import isfile

logging.basicConfig(level=logging.DEBUG)

import ROOT
from ROOT import TH1D, TFile, TTree, gROOT
from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
from ROOT import TLegend
#from ROOT import THStack

gROOT.SetBatch()



file_pattern = "/gstore/t3cms/store/user/otoldaie/v19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v19_MC2016_Summer16_TTJets_powheg/180226_022336/0000/MC2016_Summer16_TTJets_powheg_%d.root"

tau_all    = TH1D("os_tauall", "", 44, -2, 20)
tau_b      = TH1D("os_taub", "", 44, -2, 20)
tau_w      = TH1D("os_tauw", "", 44, -2, 20)
tau_w_c    = TH1D("os_tauw_c", "", 44, -2, 20)
tau_w_notc = TH1D("os_tauw_notc", "", 44, -2, 20)

ss_tau_all    = TH1D("ss_tauall", "", 44, -2, 20)
ss_tau_b      = TH1D("ss_taub", "", 44, -2, 20)
ss_tau_w      = TH1D("ss_tauw", "", 44, -2, 20)
ss_tau_w_c    = TH1D("ss_tauw_c", "", 44, -2, 20)
ss_tau_w_notc = TH1D("ss_tauw_notc", "", 44, -2, 20)

for filename in (file_pattern % i for i in range(100)):
    if not isfile(filename):
        continue

    logging.info(filename)

    tfile = TFile(filename)
    ttree = tfile.Get("ntupler/reduced_ttree")

    # OS
    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_id[0] * lep_id[0] < 0")
    tau_all.Add(ttree.GetHistogram())

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_id[0] * lep_id[0] < 0 && tau_matching_gen[0] == 5")
    tau_b.Add(ttree.GetHistogram())

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_id[0] * lep_id[0] < 0 && tau_matching_gen[0] == 4")
    tau_w.Add(ttree.GetHistogram())

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_id[0] * lep_id[0] < 0 && tau_matching_gen[0] == 4 && tau_dR_matched_jet[0] >= 0 && abs(jet_hadronFlavour[tau_dR_matched_jet[0]]) == 4")
    tau_w_c.Add(ttree.GetHistogram())

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_id[0] * lep_id[0] < 0 && tau_matching_gen[0] == 4 && tau_dR_matched_jet[0] >= 0 && abs(jet_hadronFlavour[tau_dR_matched_jet[0]]) != 4")
    tau_w_notc.Add(ttree.GetHistogram())

    # SS
    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_id[0] * lep_id[0] > 0")
    ss_tau_all.Add(ttree.GetHistogram())

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_id[0] * lep_id[0] > 0 && tau_matching_gen[0] == 5")
    ss_tau_b.Add(ttree.GetHistogram())

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_id[0] * lep_id[0] > 0 && tau_matching_gen[0] == 4")
    ss_tau_w.Add(ttree.GetHistogram())

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_id[0] * lep_id[0] > 0 && tau_matching_gen[0] == 4 && tau_dR_matched_jet[0] >= 0 && abs(jet_hadronFlavour[tau_dR_matched_jet[0]]) == 4")
    ss_tau_w_c.Add(ttree.GetHistogram())

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_id[0] * lep_id[0] > 0 && tau_matching_gen[0] == 4 && tau_dR_matched_jet[0] >= 0 && abs(jet_hadronFlavour[tau_dR_matched_jet[0]]) != 4")
    ss_tau_w_notc.Add(ttree.GetHistogram())

    tfile.Close()

tau_b      .SetLineColor(kGreen)
tau_w      .SetLineColor(kRed)
tau_w_c    .SetLineColor(kYellow)
tau_w_notc .SetLineColor(kYellow)
tau_w_notc .SetLineStyle(7)

ss_tau_b      .SetLineColor(kGreen)
ss_tau_w      .SetLineColor(kRed)
ss_tau_w_c    .SetLineColor(kYellow)
ss_tau_w_notc .SetLineColor(kYellow)
ss_tau_w_notc .SetLineStyle(7)

leg = TLegend(0.5, 0.5, 0.9, 0.9)
leg.SetName("legOS")

leg.AddEntry(tau_all    , "tt->lj, OS, all", 'L')
leg.AddEntry(tau_b      , "tt->lj, OS, b prod", 'L')
leg.AddEntry(tau_w      , "tt->lj, OS, w prod", 'L')
leg.AddEntry(tau_w_c    , "tt->lj, OS, w, C flav", 'L')
leg.AddEntry(tau_w_notc , "tt->lj, OS, w, not C", 'L')

leg2 = TLegend(0.5, 0.5, 0.9, 0.9)
leg2.SetName("legSS")

leg2.AddEntry(tau_all    , "tt->lj, SS, all", 'L')
leg2.AddEntry(tau_b      , "tt->lj, SS, b prod", 'L')
leg2.AddEntry(tau_w      , "tt->lj, SS, w prod", 'L')
leg2.AddEntry(tau_w_c    , "tt->lj, SS, w, C flav", 'L')
leg2.AddEntry(tau_w_notc , "tt->lj, SS, w, not C", 'L')

fout = TFile("sumup.root", "RECREATE")
fout.Write()

tau_all    .Write()
tau_b      .Write()
tau_w      .Write()
tau_w_c    .Write()
tau_w_notc .Write()

ss_tau_all    .Write()
ss_tau_b      .Write()
ss_tau_w      .Write()
ss_tau_w_c    .Write()
ss_tau_w_notc .Write()

leg  .Write()
leg2 .Write()

fout.Write()
fout.Close()

