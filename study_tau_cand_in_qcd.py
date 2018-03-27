import argparse
import logging
from os.path import isfile




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "study tau cands in OS/SS iso/antiiso MC QCD",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

#parser.add_argument('draw_com', type=str, help='Draw("??")')
#parser.add_argument('--cond-com', type=str, default="", help='Draw("", "??")')
#parser.add_argument('--histo-name',  type=str, default="out_histo", help='the ROOTName for output')
#parser.add_argument('--histo-range', type=str, default=None, help='optionally set the range')
#parser.add_argument('--histo-color', type=str, default=None, help='optional rgb color, like `255,255,255`')

#parser.add_argument("--per-weight",  action='store_true', help="normalize by event weight of datasets")

parser.add_argument("output", type=str, default="output.root", help="filename for output")
parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")

parser.add_argument('input_files', nargs='+', help="""the files to sum up, passed by shell,
as:

/gstore/t3cms/store/user/otoldaie/v19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v19_MC2016_Summer16_TTJets_powheg/180226_022336/0000/MC2016_Summer16_TTJets_powheg_*.root""")


args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)


logging.info("import ROOT")

import ROOT
from ROOT import TH1D, TH2D, TFile, TTree, gROOT
from ROOT import TMath, TLorentzVector
from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
from ROOT import TLegend
#from ROOT import THStack

logging.info("done")

from plotting_root import rgb

gROOT.SetBatch()



input_files = args.input_files



histos = {}

for isolation in (True, False):
    for os_sign in (True, False):
        channel = ('iso' if isolation else 'antiiso') + '_' + ('OS' if os_sign else 'SS')
        histos[(isolation, os_sign)] = {
          'tau_pt':       TH1D("tau_pt_%s" % channel,  "", 50, 0, 200),
          'lep_tau_mass': TH1D("lep_tau_mass_%s" % channel,  "", 50, 0, 200),
          'lep_tau_Mt':   TH1D("lep_tau_Mt_%s" % channel,  "", 50, 0, 400),
          'lep_tau_dPhi': TH1D("lep_tau_dPhi_%s"% channel,  "", 50, -3.5, 3.5),
          'lep_tau_dR':   TH1D("lep_tau_dR_%s"% channel,  "", 50, 0, 5),
          }

    channel = 'iso' if isolation else 'antiiso'
    histos[isolation] = {
          'first_ss_os':    TH1D("first_ss_os_%s" % channel,  "", 3, 0, 3),
          'ntaus_os':       TH1D("ntaus_os_%s" % channel,  "", 10, 0, 10),
          'ntaus_ss':       TH1D("ntaus_ss_%s" % channel,  "", 10, 0, 10),
          'ntaus_os_ss':    TH2D("ntaus_os_ss_%s" % channel,  "", 10, 0, 10, 10, 0, 10),
      }

def calc_n_save_params(event):
    # condition for roughly alliso presel
    event_pass = abs(event.leps_ID_allIso) == 13 and event.nleps_veto_mu_all == 0 and event.nleps_veto_el_all == 0 # or abs(event.leps_ID_allIso) == 11
    #&& tau_matching_allIso_lep_dR[0]>0.4 && lep_alliso_id[0] < 0.15
    if not event_pass: return None

    if abs(event.leps_ID_allIso) == 13:
        iso_lep = event.lep_alliso_relIso[0] < 0.15

    #lep_tlor = TLorentzVector(event.lep_alliso_p4[0])
    #tj_p4 = TLorentzVector(p4.X(), p4.Y(), p4.Z(), p4.T())
    lep_tlor = TLorentzVector(event.lep_alliso_p4[0].X(), event.lep_alliso_p4[0].Y(), event.lep_alliso_p4[0].Z(), event.lep_alliso_p4[0].T())

    ntaus_os = 0
    ntaus_ss = 0

    first_done = False

    for tau_id, tau_p4, lep_match in zip(event.tau_id, event.tau_p4, event.tau_matching_allIso_lep):
        # don't match lep
        if lep_match: continue

        is_os = tau_id * event.lep_alliso_id[0] < 0
        if is_os:
            ntaus_os += 1
        else:
            ntaus_ss += 1

        if not first_done:
            histos[iso_lep]['first_ss_os'].Fill(is_os)
            first_done = True

        histos[(iso_lep, is_os)]["tau_pt"].Fill(tau_p4.pt())
        lep_tau = event.lep_alliso_p4[0] + tau_p4
        histos[(iso_lep, is_os)]["lep_tau_mass"].Fill(lep_tau.mass())
        histos[(iso_lep, is_os)]["lep_tau_Mt"].Fill(lep_tau.Mt())
        lep_minus_tau = event.lep_alliso_p4[0] - tau_p4
        histos[(iso_lep, is_os)]["lep_tau_dPhi"].Fill(lep_minus_tau.Phi())
        #tau_tlor = TLorentzVector(tau_p4)
        tau_tlor = TLorentzVector(tau_p4.X(), tau_p4.Y(), tau_p4.Z(), tau_p4.T())
        dR = lep_tlor.DeltaR(tau_tlor)
        histos[(iso_lep, is_os)]["lep_tau_dR"].Fill(dR)

    histos[iso_lep]["ntaus_os"].Fill(ntaus_os)
    histos[iso_lep]["ntaus_ss"].Fill(ntaus_ss)
    histos[iso_lep]["ntaus_os_ss"].Fill(ntaus_os, ntaus_ss)


for filename in input_files:
    if not isfile(filename):
        logging.info("missing: " + filename)
        continue

    logging.debug(filename)

    tfile = TFile(filename)
    ttree = tfile.Get("ntupler/reduced_ttree")

    # Draw the file and sum up
    # 
    # TOFIX: without explicit range the histos won't add up
    #ttree.Draw(draw_command, args.cond_com)

    for iev, ev in enumerate(ttree):
        calc_n_save_params(ev)

    tfile.Close()

fout = TFile(args.output, "RECREATE")
fout.Write()

for isolation in (True, False):
    for os_sign in (True, False):
        for h in histos[(isolation, os_sign)].values():
            h.Write()
    for h in histos[isolation].values():
        h.Write()

fout.Write()
fout.Close()

