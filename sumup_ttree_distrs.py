import argparse
import logging
from os.path import isfile

logging.basicConfig(level=logging.DEBUG)

import ROOT
from ROOT import THStack, TLegend, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor



def get_distrs(ttree)
    tau_b = ttree.GetHistogram()
    tau_b = 
    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_matching_gen[0] == 5")

file_pattern = "/gstore/t3cms/store/user/otoldaie/v19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v19_MC2016_Summer16_TTJets_powheg/180226_022336/0000/MC2016_Summer16_TTJets_powheg_%d.root"

tau_b = TH1D("taub", "", 44, -2, 20)
tau_w = TH1D("tauw", "", 44, -2, 20)
tau_w_c    = TH1D("tauw_c", "", 44, -2, 20)
tau_w_notc = TH1D("tauw_notc", "", 44, -2, 20)

for filename in (file_pattern % i for i in range(100)):
    if not isfile(filename):
        continue

    tfile = TFile(filename)
    ttree = tfile.Get("ntupler/reduced_ttree")

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_matching_gen[0] == 5")
    tau_b.Add(ttree.GetHistogram())

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_matching_gen[0] == 4")
    tau_w.Add(ttree.GetHistogram())

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_matching_gen[0] == 4 && tau_dR_matched_jet[0] >= 0 && abs(jet_hadronFlavour[tau_dR_matched_jet[0]]) == 4")
    tau_w_c.Add(ttree.GetHistogram())

    ttree.Draw("tau_SV_geom_flightLenSign[tau_refited_index[0]]>>h(44,-2,20)", "(abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 13 || abs(gen_t_w_decay_id * gen_tb_w_decay_id) == 11) && tau_refited_index[0] >= 0 && tau_IDlev[0] > 2 && tau_matching_gen[0] == 4 && tau_dR_matched_jet[0] >= 0 && abs(jet_hadronFlavour[tau_dR_matched_jet[0]]) != 4")
    tau_w_notc.Add(ttree.GetHistogram())

    tfile.Close()

fout = TFile("sumup.root", "RECREATE")
fout.Write()

tau_b      .Write()
tau_w      .Write()
tau_w_c    .Write()
tau_w_notc .Write()

fout.Write()
fout.Close()

