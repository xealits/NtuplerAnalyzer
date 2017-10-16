import logging
from numpy import array, sqrt
from numpy.linalg import norm
from sys import argv


gen_match = '-g' in argv
logging.basicConfig(level=logging.DEBUG)


# ROOT is heavy, init it after comline
from ROOT import TFile, TTree, TH1D, TH2D, gROOT, gSystem, TCanvas
gROOT.SetBatch(True)

if '-w' in argv:
    filename, nick = "/eos/user/o/otoldaie/ttbar-leptons-80X_data/v12.2_merged-sets/MC2016_Summer16_WJets_madgraph.root", 'wjets'
else:
    filename, nick = "/eos/user/o/otoldaie/ttbar-leptons-80X_data/v12.2_merged-sets/MC2016_Summer16_TTJets_powheg.root" , 'tt'

isTT = 'TT' in filename

in_file = TFile(filename, "read")
ntuple  = in_file.Get("ntupler/reduced_ttree")

logging.debug("opened file and ntuple")

out_filename = "test_full_loop_on_taus_%s.root" % (nick if not gen_match else "genMatch_" + nick)
logging.debug("will write to " + out_filename)
out_file = TFile(out_filename, "recreate")
h_refit_SV_cov00 = TH1D("refit_SV_cov00_%s" % nick, "", 100, 0, 0.00001)
h_refit_SV_cov01 = TH1D("refit_SV_cov01_%s" % nick, "", 100, 0, 0.00001)
h_refit_SV_cov02 = TH1D("refit_SV_cov02_%s" % nick, "", 100, 0, 0.00001)
h_refit_SV_cov10 = TH1D("refit_SV_cov10_%s" % nick, "", 100, 0, 0.00001)
h_refit_SV_cov11 = TH1D("refit_SV_cov11_%s" % nick, "", 100, 0, 0.00001)
h_refit_SV_cov12 = TH1D("refit_SV_cov12_%s" % nick, "", 100, 0, 0.00001)
h_refit_SV_cov20 = TH1D("refit_SV_cov20_%s" % nick, "", 100, 0, 0.00001)
h_refit_SV_cov21 = TH1D("refit_SV_cov21_%s" % nick, "", 100, 0, 0.00001)
h_refit_SV_cov22 = TH1D("refit_SV_cov22_%s" % nick, "", 100, 0, 0.00001)

h_refit_SV_cov00.SetDirectory(out_file)
h_refit_SV_cov01.SetDirectory(out_file)
h_refit_SV_cov02.SetDirectory(out_file)
h_refit_SV_cov11.SetDirectory(out_file)
h_refit_SV_cov11.SetDirectory(out_file)
h_refit_SV_cov12.SetDirectory(out_file)
h_refit_SV_cov22.SetDirectory(out_file)
h_refit_SV_cov21.SetDirectory(out_file)
h_refit_SV_cov22.SetDirectory(out_file)

h_refit_flightLen_to_other_PV = TH1D("refit_flightLen_to_other_PV_%s" % nick, "", 100, 0, 0.5)

h_refit_flightLen     = TH1D("refit_flightLen_%s" % nick, "", 100, 0, 0.1)
h_refit_flightSign    = TH1D("refit_flightSign_%s" % nick, "", 100, 0, 50)
h_refit_flightLen_lt  = TH1D("refit_flightLen_lt", "", 100, 0, 0.1)
h_refit_flightLen_lj  = TH1D("refit_flightLen_lj", "", 100, 0, 0.1)
h_refit_flightSign_lt = TH1D("refit_flightSign_lt", "", 100, 0, 50)
h_refit_flightSign_lj = TH1D("refit_flightSign_lj", "", 100, 0, 50)

h_pat_flightLen = TH1D("pat_flightLen_%s" % nick, "", 100, 0, 0.1)
h_pat_flightLen_lt = TH1D("pat_flightLen_lt", "", 100, 0, 0.1)
h_pat_flightLen_lj = TH1D("pat_flightLen_lj", "", 100, 0, 0.1)
h_pat_flightSign = TH1D("pat_flightSign_%s" % nick, "", 100, 0, 50)
h_pat_flightSign_lt = TH1D("pat_flightSign_lt", "", 100, 0, 50)
h_pat_flightSign_lj = TH1D("pat_flightSign_lj", "", 100, 0, 50)

h_pat_bTag_flightSign    = TH2D("pat_bTag_flightSign_%s" % nick, "", 100, 0, 50, 100, 0, 1)
h_pat_bTag_flightSign_lt = TH2D("pat_bTag_flightSign_lt",        "", 100, 0, 50, 100, 0, 1)
h_pat_bTag_flightSign_lj = TH2D("pat_bTag_flightSign_lj",        "", 100, 0, 50, 100, 0, 1)
#tau_dR_matched_jet

h_pat_bTag_refit_flightSign    = TH2D("pat_bTag_refit_flightSign_%s" % nick, "", 100, 0, 50, 100, 0, 1)
h_pat_bTag_refit_flightSign_lt = TH2D("pat_bTag_refit_flightSign_lt",        "", 100, 0, 50, 100, 0, 1)
h_pat_bTag_refit_flightSign_lj = TH2D("pat_bTag_refit_flightSign_lj",        "", 100, 0, 50, 100, 0, 1)
#tau_dR_matched_jet

h_refit_flightSign_lt_Btag = TH1D("refit_flightSign_lt_Btag", "", 100, 0, 50)
h_refit_flightSign_lj_Btag = TH1D("refit_flightSign_lj_Btag", "", 100, 0, 50)
h_pat_flightSign_lt_Btag   = TH1D("pat_flightSign_lt_Btag",   "", 100, 0, 50)
h_pat_flightSign_lj_Btag   = TH1D("pat_flightSign_lj_Btag",   "", 100, 0, 50)

h_refit_flightSign_lt_noBtag = TH1D("refit_flightSign_lt_noBtag", "", 100, 0, 50)
h_refit_flightSign_lj_noBtag = TH1D("refit_flightSign_lj_noBtag", "", 100, 0, 50)
h_pat_flightSign_lt_noBtag   = TH1D("pat_flightSign_lt_noBtag",   "", 100, 0, 50)
h_pat_flightSign_lj_noBtag   = TH1D("pat_flightSign_lj_noBtag",   "", 100, 0, 50)

h_n_goodPV = TH1D("n_goodPV_%s" % nick, "", 50, 0, 50)

logging.debug("set histograms")

N_events = 100000 # ntuple.GetEntries()
n_events = 0

# looping:
for i, event in enumerate(ntuple):
    if event.tau_SV_fit_isOk.size() < 1 or not event.tau_SV_fit_isOk[0]: continue
    n_events += 1
    if n_events > N_events: break

    # match the 0 tau to any of gen taus
    if gen_match and not any(sqrt((gen.Eta() - event.tau_p4[0].Eta())**2 + (gen.Phi() - event.tau_p4[0].Phi())**2) < 0.1 for gen in event.gen_tau_p4): continue
    # ref:
    # http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_5_3_9/doc/html/d5/d6b/DataFormats_2Math_2interface_2deltaR_8h_source.html#l00019

    h_n_goodPV.Fill(event.PV_x.size())
    h_refit_SV_cov00.Fill(event.tau_SV_cov[0](0, 0))
    h_refit_SV_cov01.Fill(event.tau_SV_cov[0](0, 1))
    h_refit_SV_cov02.Fill(event.tau_SV_cov[0](0, 2))
    h_refit_SV_cov10.Fill(event.tau_SV_cov[0](1, 1))
    h_refit_SV_cov11.Fill(event.tau_SV_cov[0](1, 1))
    h_refit_SV_cov12.Fill(event.tau_SV_cov[0](1, 2))
    h_refit_SV_cov20.Fill(event.tau_SV_cov[0](2, 2))
    h_refit_SV_cov21.Fill(event.tau_SV_cov[0](2, 1))
    h_refit_SV_cov22.Fill(event.tau_SV_cov[0](2, 2))
    pv = array((event.PV_x[0], event.PV_y[0], event.PV_z[0]))
    pv_err = array((event.PV_x_err[0], event.PV_y_err[0], event.PV_z_err[0]))
    sv = array((event.tau_SV_fit_x[0], event.tau_SV_fit_y[0], event.tau_SV_fit_z[0]))
    sv_err = array((event.tau_SV_cov[0](0,0), event.tau_SV_cov[0](1,1), event.tau_SV_cov[0](2,2)))
    flight_vec = pv - sv
    flight_len = norm(flight_vec)
    flight_vec_err = sqrt(sv_err**2 + pv_err**2)
    flight_sign = norm(flight_vec / flight_vec_err)
    pat_flightLen  = event.tau_flightLength[0]
    pat_flightSign = event.tau_flightLengthSignificance[0]

    if event.PV_x.size() > 1:
      for i in range(1, event.PV_x.size()):
        pv_i = array((event.PV_x[i], event.PV_y[0], event.PV_z[i]))
        flight_len_i = norm(pv_i - sv)
        h_refit_flightLen_to_other_PV.Fill(flight_len_i)


    #h_pat_bTag_flightSign_lt = TH2D("pat_bTag_flightSign_lt",        "", 100, 0, 1, 100, 0, 50)
    #h_pat_bTag_flightSign_lj = TH2D("pat_bTag_flightSign_lj",        "", 100, 0, 1, 100, 0, 50)
    ##tau_dR_matched_jet
    h_pat_flightLen .Fill(pat_flightLen)
    h_pat_flightSign.Fill(pat_flightSign)
    h_refit_flightLen .Fill(flight_len)
    h_refit_flightSign.Fill(flight_sign)
    i_matched_tau_jet = event.tau_dR_matched_jet[0]
    if i_matched_tau_jet > -1:
        h_pat_bTag_flightSign.Fill(pat_flightSign, event.jet_b_discr[i_matched_tau_jet])
        h_pat_bTag_refit_flightSign.Fill(flight_sign, event.jet_b_discr[i_matched_tau_jet])

    if isTT:
        if abs(event.gen_t_w_decay_id * event.gen_tb_w_decay_id) == 13: # lj
            h_refit_flightLen_lj .Fill(flight_len)
            h_refit_flightSign_lj.Fill(flight_sign)
            h_pat_flightLen_lj .Fill(pat_flightLen)
            h_pat_flightSign_lj.Fill(pat_flightSign)
            if i_matched_tau_jet > -1:
                h_pat_bTag_flightSign_lj.Fill(pat_flightSign, event.jet_b_discr[i_matched_tau_jet])
                h_pat_bTag_refit_flightSign_lj.Fill(flight_sign, event.jet_b_discr[i_matched_tau_jet])
                if event.jet_b_discr[i_matched_tau_jet] < 0.85:
                    h_refit_flightSign_lj_noBtag .Fill(flight_sign)
                    h_pat_flightSign_lj_noBtag   .Fill(pat_flightSign)
                else:
                    h_refit_flightSign_lj_Btag .Fill(flight_sign)
                    h_pat_flightSign_lj_Btag   .Fill(pat_flightSign)

        if (abs(event.gen_t_w_decay_id) > 15*15 and abs(event.gen_tb_w_decay_id) == 13) or (abs(event.gen_t_w_decay_id) == 13 and abs(event.gen_tb_w_decay_id) > 15*15): # lt
            h_refit_flightLen_lt .Fill(flight_len)
            h_refit_flightSign_lt.Fill(flight_sign)
            h_pat_flightLen_lt .Fill(pat_flightLen)
            h_pat_flightSign_lt.Fill(pat_flightSign)
            if i_matched_tau_jet > -1:
                h_pat_bTag_flightSign_lt.Fill(pat_flightSign, event.jet_b_discr[i_matched_tau_jet])
                h_pat_bTag_refit_flightSign_lt.Fill(flight_sign, event.jet_b_discr[i_matched_tau_jet])
                if event.jet_b_discr[i_matched_tau_jet] < 0.85:
                    h_refit_flightSign_lt_noBtag .Fill(flight_sign)
                    h_pat_flightSign_lt_noBtag   .Fill(pat_flightSign)
                else:
                    h_refit_flightSign_lt_Btag .Fill(flight_sign)
                    h_pat_flightSign_lt_Btag   .Fill(pat_flightSign)



logging.debug("writing output")

out_file.Write()

