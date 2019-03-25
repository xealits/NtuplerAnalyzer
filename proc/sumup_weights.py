from datetime import datetime
import argparse
import logging
from os.path import isfile, basename
from sys import exit


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "hadd the specified histos in the files",
    epilog = """Example:
time python sumup_weights.py --overwrite --debug temp/test_weight.root lstore_outdirs/v37/test1/MC2016_Summer16_TTJets_powheg/MC2016_Summer16_TTJets_powheg_*.root
python sumup_weights.py --overwrite --debug temp/test_weight.root --target-histos ntupler/weight_counter,ntupler/events_counter,ntupler/systematic_weights gstore_outdirs/v37/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v37_MC2016_Summer16_TTJets_powheg/190207_134619/0000/MC2016_Summer16_TTJets_powheg_*.root
"""
    )

parser.add_argument('out_file',   help='output file')

parser.add_argument('-t', '--target-histos', type=str, default='weight_counter,events_counter,systematic_weights,control_counters', help="the histograms to hadd `histo1,path/to/h2,etc`")

'''
  KEY: TH1D	events_counter;1	pass category
  KEY: TH1D	weight_counter;1	pass category
  KEY: TH1D	systematic_weights;1	pass category
  KEY: TH1D	control_counters;1	
'''

parser.add_argument("--overwrite", action='store_true', help="overwrite output")
parser.add_argument("--debug", action='store_true', help="set logging DEBUG")

parser.add_argument('inp_files', nargs='+', help='the files with weight histos to sumup')


args = parser.parse_args()

if isfile(args.out_file) and not args.overwrite:
    print 'output file exists: %s' % args.out_file
    exit(0)

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

logging.debug('importing ROOT')
preroot_time = datetime.now()
import ROOT
from ROOT import TFile,TDirectory
from ROOT import gDirectory
logging.debug('done')

start_time = datetime.now()

# the buffer histogram for adding all the weights in the files
# it is written to the output at the end
#weight_counter = None
hadd_histos = {}
hadd_histos_names = args.target_histos.split(',')

for filename in args.inp_files:
    if not isfile(filename):
        print 'the input file does not exist or is not a file: %s' % filename
        continue

    tfile = TFile(filename)

    for hname in hadd_histos_names:
        hname_h = tfile.Get(hname)
        assert bool(hname_h)
        if hname in hadd_histos:
            hadd_histos[hname].Add(hname_h)
        else:
            hname_h_0 = hname_h.Clone()
            hname_h_0.SetDirectory(0)
            hname_h_0.SetName(hname.split('/')[-1])
            hadd_histos[hname] = hname_h_0

fout = TFile(args.out_file, "RECREATE")
fout.Write()

# save weight counter etc in the top directory of the output file
fout.cd()

#weight_counter.SetName('weight_counter')
#weight_counter.Write()
for histo in hadd_histos.values():
    histo.Write()

fout.Write()
fout.Close()

