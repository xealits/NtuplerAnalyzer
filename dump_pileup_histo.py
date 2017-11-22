import argparse
import logging


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = 'save the "pileup" TH1D as csv or .C macro to get the numbers in text',
    )

parser.add_argument("pu_file",    help='file with the TH1D "pileup"')
parser.add_argument("--th-name", type=str, default='pileup', help='name of the TH1D to dump ("pileup" by default)')
parser.add_argument("--debug", action='store_true', help="DEBUG logging threshold")
parser.add_argument("--info",  action='store_true', help="INFO logging threshold")
parser.add_argument("--raw",   action='store_true', help="dump .C macro")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
elif args.info:
    logging.basicConfig(level=logging.INFO)
else:
    logging.basicConfig(level=logging.WARNING)



logging.info("import ROOT")

import ROOT
from ROOT import TFile

pu_file = TFile(args.pu_file)

histo = pu_file.Get(args.th_name)

if args.raw:
    histo.SaveAs(args.pu_file.replace('.root', '_%s.C' % args.th_name))
else:
    print "bin_n,pu,pu_err"
    for bin_n in range(histo.GetXaxis().GetNbins()):
        print "%d,%f,%f" % (bin_n, histo.GetBinContent(bin_n), histo.GetBinError(bin_n))

