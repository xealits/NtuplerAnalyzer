import argparse
import logging
from os.path import isfile, basename




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "sumup 1 hist in root files",
    epilog = """Example:\npython sumup_histo.py histoname gstore_outdirs/v28/SingleMuon/Ntupler_v28_Data13TeV_SingleMuon2016*/1*/*/*root"""
    )

parser.add_argument('histo_name', type=str, help='full name of TH from root of the files, like ntupler/weight_counter')

parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")
parser.add_argument("--output", type=str, default="output.root", help="filename for output")

parser.add_argument('input_files', nargs='+', help="""the files to sum up, passed by shell, as:

/gstore/t3cms/store/user/otoldaie/v19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v19_MC2016_Summer16_TTJets_powheg/180226_022336/0000/MC2016_Summer16_TTJets_powheg_*.root""")


args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)


logging.info("import ROOT")

import ROOT
from ROOT import TH1D, TFile, TTree, gROOT
from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
from ROOT import TLegend
#from ROOT import THStack

logging.info("done")

from plotting_root import rgb

gROOT.SetBatch()



input_files = args.input_files


out_histo = None

for filename in input_files:
    if not isfile(filename):
        logging.info("missing: " + filename)
        continue

    logging.debug(filename)

    tfile  = TFile(filename)
    thisto = tfile.Get(args.histo_name)

    if not out_histo:
        out_histo = thisto.Clone()
        out_histo.Clone()
        out_histo.SetDirectory(0)
    else:
        out_histo.Add(thisto)

    tfile.Close()

fout = TFile(args.output, "RECREATE")
fout.Write()

out_histo.Write()

fout.Write()
fout.Close()

