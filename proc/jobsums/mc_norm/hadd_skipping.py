import argparse
import logging
from os.path import isfile, basename
from math import sqrt


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "hadd histos for os/ss clusre test, skipping the bins with no content",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument("--outfile", type=str, default='hadd_skipped.root', help="the filename for output")
parser.add_argument("--debug", action='store_true', help="DEBUG level of logging")
parser.add_argument('input_files', nargs='+', help="the files to hadd up")


args = parser.parse_args()


if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)


logging.debug("import ROOT")

import ROOT
from ROOT import TH1D, TFile, TTree, gROOT
from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
from ROOT import TLegend
#from ROOT import THStack

logging.debug("done")



first_file = "MC2016_Summer16_QCD_EMEnriched_Pt-20to30.root"

files_to_hadd = [
"MC2016_Summer16_QCD_EMEnriched_Pt-30to50.root",
"MC2016_Summer16_QCD_EMEnriched_Pt-50to80.root",
"MC2016_Summer16_QCD_EMEnriched_Pt-80to120.root",
"MC2016_Summer16_QCD_EMEnriched_Pt-120to170.root",
"MC2016_Summer16_QCD_EMEnriched_Pt-170to300.root",
"MC2016_Summer16_QCD_EMEnriched_Pt-300toInf.root",
]

first_file = args.input_files[0]
files_to_hadd = args.input_files[1:]

"MC2016_Summer16_QCD_EMEnriched_Pt-30to40.root", "MC2016_Summer16_QCD_EMEnriched_Pt-30toInf.root", "MC2016_Summer16_QCD_EMEnriched_Pt-40toInf.root"



#TH1D* h_el    = (TH1D*) alliso_el_sel   ->Get("qcd/NOMINAL/alliso_el_sel_qcd_NOMINAL_relIso")
#TH1D* h_el_ss = (TH1D*) alliso_el_sel_ss->Get("qcd/NOMINAL/alliso_el_sel_ss_qcd_NOMINAL_relIso")

first_tfile = TFile(first_file)

h_el    = first_tfile.Get("alliso_el_sel/qcd/NOMINAL/alliso_el_sel_qcd_NOMINAL_relIso").Clone()
h_el_ss = first_tfile.Get("alliso_el_sel_ss/qcd/NOMINAL/alliso_el_sel_ss_qcd_NOMINAL_relIso").Clone()

h_el    .SetName("h_el")
h_el_ss .SetName("h_el_ss")

h_el    .Reset()
h_el_ss .Reset()

for fname in files_to_hadd:
    logging.info(fname)
    tfile = TFile(fname)

    hadd_el    = tfile.Get("alliso_el_sel/qcd/NOMINAL/alliso_el_sel_qcd_NOMINAL_relIso")
    hadd_el_ss = tfile.Get("alliso_el_sel_ss/qcd/NOMINAL/alliso_el_sel_ss_qcd_NOMINAL_relIso")

    # clear bins if there is no content in one of the hists
    for bini in range(hadd_el.GetSize()):
        # how to handle zero?
        #if hadd_el.GetBinContent(bini) == 0 or hadd_el_ss.GetBinContent(bini) == 0:
        cont_os = hadd_el    .GetBinContent(bini)
        cont_ss = hadd_el_ss .GetBinContent(bini)
        e = 0.1
        if (cont_os < e * cont_ss) or (cont_ss < e * cont_os):
            hadd_el   .SetBinContent(bini, 0)
            hadd_el_ss.SetBinContent(bini, 0)

    h_el    .Add(hadd_el   )
    h_el_ss .Add(hadd_el_ss)



# write out
outfile = TFile(args.outfile, 'RECREATE')
outfile.cd()
outfile.Write()

h_el    .Write()
h_el_ss .Write()

outfile.Write()
outfile.Close()


