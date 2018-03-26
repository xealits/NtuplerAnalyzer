import argparse
import logging
from os.path import isfile




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "study PDF weights etc",
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
from ROOT import TH1D, TFile, TTree, gROOT
from ROOT import TMath
from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
from ROOT import TLegend
#from ROOT import THStack

logging.info("done")

from plotting_root import rgb

gROOT.SetBatch()



input_files = args.input_files



histo_pdf_nom  = TH1D("pdf_nom",  "", 50, 0, 2)
histo_pdf_up   = TH1D("pdf_up",   "", 50, 0, 2)
histo_pdf_down = TH1D("pdf_down", "", 50, 0, 2)
histo_alpha_1  = TH1D("alpha_1",  "", 50, 0, 2)
histo_alpha_2  = TH1D("alpha_2",  "", 50, 0, 2)

def calc_weights(event):
    pdf_nom = event.gen_weights_pdf_hessians[0]
    if not pdf_nom > 0:
        print event.indexevents, event.runNumber, event.lumi
        return 0, 0, 0, 0, 0
    alpha1  = event.gen_weight_alphas_1 / pdf_nom
    alpha2  = event.gen_weight_alphas_2 / pdf_nom
    pdf_up   = 0
    pdf_down = 0
    for pdf_w in event.gen_weights_pdf_hessians[1:]:
        pdf_w /= pdf_nom
        if pdf_w > 1.:
            pdf_up += pdf_w*pdf_w
        else:
            pdf_down += pdf_w*pdf_w
    return pdf_nom, alpha1, alpha2, pdf_up, pdf_down


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
        pdf_nom, alpha1, alpha2, pdf_up, pdf_down = calc_weights(ev)

        histo_pdf_nom  .Fill(pdf_nom)
        histo_pdf_up   .Fill(pdf_up)
        histo_pdf_down .Fill(pdf_down)
        histo_alpha_1  .Fill(alpha1)
        histo_alpha_2  .Fill(alpha2)

    tfile.Close()

fout = TFile(args.output, "RECREATE")
fout.Write()

histo_pdf_nom  .Write()
histo_pdf_up   .Write()
histo_pdf_down .Write()
histo_alpha_1  .Write()
histo_alpha_2  .Write()

fout.Write()
fout.Close()

