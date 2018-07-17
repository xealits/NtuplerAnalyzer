from os.path import isfile
import argparse
import logging



parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "rebin histograms to 1 bin",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument('distr', type=str, default='Mt_lep_met_f', help="distribution name to rebin")
parser.add_argument('files', type=str, default=None, help="list of filesnames to process, without .root separated by coma, if not given the default batch is used")
parser.add_argument("--debug", action="store_true", help="debug logging")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)


logging.info("import ROOT")

import ROOT
from ROOT import TFile
#from ROOT.TObject import kOverwrite


files = args.files.split(',') if args.files else files_default

logging.info(repr(files))

for filename in files:
    # check if file exists
    if not isfile(filename):
        logging.warning("file doesn't exist: %s" % filename)
        continue

    f = TFile(filename, "UPDATE")

    logging.info("%s" % (filename))

    for chan in list(f.GetListOfKeys()):
        channel = chan.GetName()
        if channel in ('weight_counter', 'events_counter', 'control_counters', 'systematic_weights'):
            continue
        #try:
        for proc in list(chan.ReadObj().GetListOfKeys()):
            #if process == 'mu_presel':
                #logging.debug(str(proc))
            nick = proc.GetName()
            for sys in list(proc.ReadObj().GetListOfKeys()):
                sys_name = sys.GetName()

                histoname = "{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}".format(chan=channel, proc=nick, sys=sys_name, distr=args.distr)
                histo = f.Get(histoname)
                if not histo:
                    continue
                logging.debug("%s" % (histoname))
                N_bins_to_merge = histo.GetNbinsX()
                histo.Rebin(N_bins_to_merge)

    f.Write("", TFile.kOverwrite)
    f.Close()


