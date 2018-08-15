import argparse
import logging
from os.path import isfile, basename
from math import sqrt

#parser = argparse.ArgumentParser(
#    formatter_class = argparse.RawDescriptionHelpFormatter,
#    description = "sumup 1 hist in root files",
#    epilog = """Example:\npython sumup_histo.py histoname gstore_outdirs/v28/SingleMuon/Ntupler_v28_Data13TeV_SingleMuon2016*/1*/*/*root"""
#    )
#
#parser.add_argument('histo_name', type=str, help='full name of TH from root of the files, like ntupler/weight_counter')
#
#parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")
#parser.add_argument("--output", type=str, default="output.root", help="filename for output")
#
#parser.add_argument('input_files', nargs='+', help="""the files to sum up, passed by shell, as:
#
#/gstore/t3cms/store/user/otoldaie/v19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v19_MC2016_Summer16_TTJets_powheg/180226_022336/0000/MC2016_Summer16_TTJets_powheg_*.root""")


logging.basicConfig(level=logging.INFO)

logging.info("import ROOT")

import ROOT
from ROOT import TH1D, TFile, TTree, gROOT
from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
from ROOT import TLegend
#from ROOT import THStack

logging.info("done")




nom_file  = "merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg.root"
sys_files = [
"merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_CUETP8M2T4down.root",
"merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_CUETP8M2T4up.root",
"merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_fsrdown.root",
"merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_fsrup.root",
"merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_hdampDOWN.root",
"merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_hdampUP.root",
"merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_isrdown.root",
"merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_isrup.root",
]


#  KEY: TH1D	all_ev_eltau;1	
#  KEY: TH1D	cut_ev_eltau;1	
#  KEY: TH1D	all_ev_mutau;1	
#  KEY: TH1D	cut_ev_mutau;1	

tfile = TFile(nom_file)

histo_all_ev_eltau = tfile.Get("all_ev_eltau")
histo_cut_ev_eltau = tfile.Get("cut_ev_eltau")
histo_all_ev_mutau = tfile.Get("all_ev_mutau")
histo_cut_ev_mutau = tfile.Get("cut_ev_mutau")


histo_both_all = histo_all_ev_eltau + histo_all_ev_mutau
histo_both_cut = histo_cut_ev_eltau + histo_cut_ev_mutau


# the nominal acceptance and weight variations
histo_both_cut.Divide(histo_both_all)
histo_cut_ev_mutau.Divide(histo_all_ev_mutau)
histo_cut_ev_eltau.Divide(histo_all_ev_eltau)

acc_nom_el, acc_nom_mu, acc_nom_both = histo_cut_ev_eltau.GetBinContent(1), histo_cut_ev_mutau.GetBinContent(1), histo_both_cut.GetBinContent(1)
print [acc_nom_el, acc_nom_mu, acc_nom_both]

variations = [0, 0, 0]

# TODO: weight variations
for var_index in range(2,100):
    acc_var_el, acc_var_mu, acc_var_both = histo_cut_ev_eltau.GetBinContent(var_index), histo_cut_ev_mutau.GetBinContent(var_index), histo_both_cut.GetBinContent(var_index)
    if acc_var_el < 0.01: break

    var_el   = (acc_var_el   - acc_nom_el)**2
    var_mu   = (acc_var_mu   - acc_nom_mu)**2
    var_both = (acc_var_both - acc_nom_both)**2

    #print var_el, var_mu, var_both

    variations[0] += var_el
    variations[1] += var_mu
    variations[2] += var_both

tfile.Close()

# dataset variations
for sys_filename in sys_files:
    sys_file = TFile(sys_filename)

    histo_all_ev_eltau = sys_file.Get("all_ev_eltau")
    histo_cut_ev_eltau = sys_file.Get("cut_ev_eltau")
    histo_all_ev_mutau = sys_file.Get("all_ev_mutau")
    histo_cut_ev_mutau = sys_file.Get("cut_ev_mutau")

    histo_both_all = histo_all_ev_eltau + histo_all_ev_mutau
    histo_both_cut = histo_cut_ev_eltau + histo_cut_ev_mutau

    histo_both_cut.Divide(histo_both_all)
    histo_cut_ev_mutau.Divide(histo_all_ev_mutau)
    histo_cut_ev_eltau.Divide(histo_all_ev_eltau)

    acc_sys_el, acc_sys_mu, acc_sys_both = histo_cut_ev_eltau.GetBinContent(1), histo_cut_ev_mutau.GetBinContent(1), histo_both_cut.GetBinContent(1)

    var_el   = (acc_sys_el   - acc_nom_el)**2
    var_mu   = (acc_sys_mu   - acc_nom_mu)**2
    var_both = (acc_sys_both - acc_nom_both)**2

    variations[0] += var_el
    variations[1] += var_mu
    variations[2] += var_both

    logging.info("%-80s %f %f %f" % (sys_filename, acc_sys_el, acc_sys_mu, acc_sys_both))

    sys_file.Close()

print [sqrt(v) for v in variations]

