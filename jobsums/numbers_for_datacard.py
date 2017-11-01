'''
bin         mu_lj     mu_lj     mu_lj     mu_lj         mu_lj       mu_lj       mu_lj    mu_lj        mu_lj     mu_lj     mu_lj
process     dy_tautau dy_other  wjets     s_top_mutau   s_top_other s_top_lj    dibosons tt_taultauh  tt_lj     tt_other  tt_mutau
process     10        9         8         7             6           5           4        3            2         1         0
rate        117.997   10.96     417.34    169.45        29.42       205.302     25.73    228.08       3523.18   328.59    3750.46     ------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

processes_mu = ['dy_tautau', 'dy_other',  'wjets',     's_top_mutau',   's_top_other', 's_top_lj',    'dibosons', 'tt_taultauh',  'tt_lj',     'tt_other',  'tt_mutau']
processes_mu_id = {'dy_tautau':10, 'dy_other': 9,  'wjets': 8,     's_top_mutau': 7,   's_top_other': 6, 's_top_lj': 5,    'dibosons': 4, 'tt_taultauh': 3,  'tt_lj': 2,     'tt_other': 1,  'tt_mutau':0}

processes_el = ['dy_tautau', 'dy_other',  'wjets',     's_top_eltau',   's_top_other', 's_top_lj',    'dibosons', 'tt_taultauh',  'tt_lj',     'tt_other',  'tt_eltau']
processes_el_id = {'dy_tautau':10, 'dy_other': 9,  'wjets': 8,     's_top_eltau': 7,   's_top_other': 6, 's_top_lj': 5,    'dibosons': 4, 'tt_taultauh': 3,  'tt_lj': 2,     'tt_other': 1,  'tt_eltau':0}


import argparse


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "extract observed and expected values in el_lj/el_lj_out or mu_lj/mu_lj_out channels",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument("data_file",    help="full Data+MC file name")
parser.add_argument("-c", "--channel",  type=str, help="if given then get the numbers from here")
parser.add_argument("-m", "--mu",  action = "store_true", help="muon channels: mu_lj, mu_lj_out")
parser.add_argument("-e", "--el",  action = "store_true", help="electron channels: el_lj, el_lj_out")
parser.add_argument("-v", "--verbose",  action = "store_true", help="log debug output")


args = parser.parse_args()

import logging

if args.verbose:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

logging.info("import ROOT")

import ROOT
from ROOT import TFile, THStack, TLegend, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
import plotting_root


if args.channel:
    channels = [args.channel]
    processes = processes_mu
    processes_id = processes_mu_id
elif args.mu:
    channels = ['mu_lj', 'mu_lj_out']
    processes = processes_mu
    processes_id = processes_mu_id
elif args.el:
    channels = ['el_lj', 'el_lj_out']
    processes = processes_el
    processes_id = processes_el_id
else:
    channels = []

fdata = TFile(args.data_file)

sys_name = 'NOMINAL'
distr_name = 'Mt_lep_met'

for channel in channels:
    chan = fdata.Get(channel)
    #if chan.GetName() == 'weight_counter' or chan.GetName() == 'events_counter':
    #    continue
    #if chan.GetName() != channel and chan.GetName() != shape_channel:
    #    continue

    if not chan:
        logging.warning('no channel %s' % channel)
        continue

    histos = []

    #processes_keys = list(chan.ReadObj().GetListOfKeys())
    #processes_keys = list(chan.GetListOfKeys())
    #sorted_pkeys = sorted(processes_keys, key=lambda pkey: pkey.GetName() if pkey.GetName() not in ('qcd') else 'z_' + pkey.GetName())
    #for process in sorted_pkeys:
    for process in processes:
       # if the channel with shapes for dy and wjets are given
       # don't use other channels
       #if shape_channel and nick in ('dy', 'wjets') and chan.GetName() != shape_channel:
       #    continue

       histo_name = '_'.join([channel, process, sys_name, distr_name])

       logging.debug(histo_name)
       logging.debug('%s/%s/%s' % (process, sys_name, histo_name))

       histo = chan.Get('%s/%s/%s' % (process, sys_name, histo_name))
       logging.debug(histo.Integral())

       histos.append(histo)


    print 'bin           ' + ''.join('%-15s' % channel for _ in processes)
    print 'process       ' + ''.join('%-15s' % proc for proc in processes)
    print 'process       ' + ''.join('%-15d' % processes_id[proc] for proc in processes)
    print 'rate          ' + ''.join('%-15.3f' % histo.Integral() for histo in histos)

    histo_name = '_'.join([channel, 'data', sys_name, distr_name])
    data_histo = fdata.Get('%s/%s/%s/%s' % (channel, 'data', sys_name, histo_name))
    print 'obs %f' % data_histo.Integral()

