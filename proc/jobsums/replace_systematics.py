import argparse
import logging
from os.path import isfile

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "replace requested systematic with another one in place",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument("data_file",  help="Data file name")
parser.add_argument("-c", "--channel",  type=str, default='ctr_old_mu_sel', help="selection channel")
parser.add_argument("-p", "--process",  type=str, default='qcd', help="the target process (usually it is done to qcd)")
parser.add_argument("-d", "--distr-name",  type=str, default='Mt_lep_met', help='the distribution')
parser.add_argument("-s", "--systematics",  type=str, default='FSRUp,NOMINAL', help="systematic variations for substitution: TARGET,REPLACEMENT (FSRUp,NOMINAL)")

args = parser.parse_args()

assert isfile(args.data_file)

logging.info("import ROOT")

import ROOT
from ROOT import TFile, THStack, TLegend, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
#import plotting_root


target_sys, replacement_sys = args.systematics.split(',')
channel = args.channel
process = args.process
distr   = args.distr_name

f = TFile(args.data_file, "UPDATE")
#f.cd()

'''
channel = f.Get(args.channel)
channel.cd()

process = f.Get(args.process)
process.cd()

sys_t
'''

f.ls()
f.Get(channel).ls()
f.Get(channel + '/' + process).ls()
f.Get(channel + '/' + process + '/' + target_sys).ls()
f.Get(channel + '/' + process + '/' + replacement_sys).ls()
# ROOT could not get the following target/replacement histograms
# this .ls() sequence was for debuggiing
# but now it works
# ROOT craft

histo_target      = f.Get('%s/%s/%s/%s' % (channel, process, target_sys, '_'.join([channel, process, target_sys, distr])))
histo_replacement = f.Get('%s/%s/%s/%s' % (channel, process, replacement_sys, '_'.join([channel, process, replacement_sys, distr])))

histo_target      .Print()
histo_replacement .Print()

# loop through each bin and substitute value
for bini in range(histo_target.GetSize()):
    repl_content = histo_replacement.GetBinContent (bini)
    repl_error   = histo_replacement.GetBinError   (bini)
    histo_target.SetBinContent (bini, repl_content)
    histo_target.SetBinError   (bini, repl_error)

f.Write("", TFile.kOverwrite)
f.Close()


