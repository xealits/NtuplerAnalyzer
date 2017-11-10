'''
bin         mu_lj     mu_lj     mu_lj     mu_lj         mu_lj       mu_lj       mu_lj    mu_lj        mu_lj     mu_lj     mu_lj
process     dy_tautau dy_other  wjets     s_top_mutau   s_top_other s_top_lj    dibosons tt_taultauh  tt_lj     tt_other  tt_mutau
process     10        9         8         7             6           5           4        3            2         1         0
rate        117.997   10.96     417.34    169.45        29.42       205.302     25.73    228.08       3523.18   328.59    3750.46     ------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

all_known_sorted_processes = [ 'dy_tautau', 'dy_other',  'wjets', 'dibosons',  's_top_mutau',   's_top_other', 's_top_lj',    'tt_taultauh',  'tt_lj',     'tt_other',  'tt_mutau', 'tt_eltau']

processes_mu = ['dy_tautau', 'dy_other',  'wjets',     's_top_mutau',   's_top_other', 's_top_lj',    'dibosons', 'tt_taultauh',  'tt_lj',     'tt_other',  'tt_mutau']
#processes_mu = set(['dy_tautau', 'dy_other',  'wjets',     's_top_mutau',   's_top_other', 's_top_lj',    'dibosons', 'tt_taultauh',  'tt_lj',     'tt_other',  'tt_mutau'])
processes_mu_id = {'dy_tautau':10, 'dy_other': 9,  'wjets': 8,     's_top_mutau': 7,   's_top_other': 6, 's_top_lj': 5,    'dibosons': 4, 'tt_taultauh': 3,  'tt_lj': 2,     'tt_other': 1,  'tt_mutau':0}

processes_el = ['dy_tautau', 'dy_other',  'wjets',     's_top_eltau',   's_top_other', 's_top_lj',    'dibosons', 'tt_taultauh',  'tt_lj',     'tt_other',  'tt_eltau']
#processes_el = set(['dy_tautau', 'dy_other',  'wjets',     's_top_eltau',   's_top_other', 's_top_lj',    'dibosons', 'tt_taultauh',  'tt_lj',     'tt_other',  'tt_eltau'])
processes_el_id = {'dy_tautau':10, 'dy_other': 9,  'wjets': 8,     's_top_eltau': 7,   's_top_other': 6, 's_top_lj': 5,    'dibosons': 4, 'tt_taultauh': 3,  'tt_lj': 2,     'tt_other': 1,  'tt_eltau':0}


import argparse


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "extract observed and expected values in el_lj/el_lj_out or mu_lj/mu_lj_out channels",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument("data_file",    help="full Data+MC file name")
parser.add_argument("-c", "--channels",  type=str, help="if given then get the numbers from here, enter multiple channels with 'mu_presel,mu_sel'")
parser.add_argument("-d", "--distr",    type=str, default='Mt_lep_met', help="fitted distribution")
parser.add_argument("-m", "--mu",  action = "store_true", help="muon channels: mu_lj, mu_lj_out")
parser.add_argument("-e", "--el",  action = "store_true", help="electron channels: el_lj, el_lj_out")
parser.add_argument("-v", "--verbose",  action = "store_true", help="log debug output")
#parser.add_argument("-y", "--event-yields", type=str, default='text', help="output in the form of event yield table (set -y latex for latex table output)")
parser.add_argument("-y", "--event-yields", action='store_true', help="output in the form of event yield table (set -y latex for latex table output)")


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


'''
if args.channel and 'mu' in args.channel:
    channels = [args.channel]
    processes = processes_mu
    processes_id = processes_mu_id
elif args.channel and 'el' in args.channel:
    channels = [args.channel]
    processes = processes_el
    processes_id = processes_el_id
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
'''

if not args.channels and args.mu:
    channels = ['mu_lj', 'mu_lj_out']
elif not args.channels and args.el:
    channels = ['el_lj', 'el_lj_out']
elif args.channels:
    channels = args.channels.split(',')

fdata = TFile(args.data_file)

sys_name = 'NOMINAL'
distr_name = args.distr

#if args.event_yields:
#proc_yields = [] #[('channel', [(proc_nick, yield)...]) ...]
proc_yields = {} # proc_nick: {channel: yield}

data_yields = {} # channel: yield

for channel in channels:
    chan = fdata.Get(channel)
    #if chan.GetName() == 'weight_counter' or chan.GetName() == 'events_counter':
    #    continue
    #if chan.GetName() != channel and chan.GetName() != shape_channel:
    #    continue
    logging.info(channel)

    if not chan:
        logging.warning('no channel %s' % channel)
        continue

    histos = []

    if 'mu' in channel:
        processes = processes_mu
        processes_id = processes_mu_id
    else:
        processes = processes_el
        processes_id = processes_el_id

    #proc_yields.append((channel, []))

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

       #proc_yields[-1][1].append((process, histo.Integral()))
       proc_yields.setdefault(process, {})[channel] = histo.Integral()

    histo_name = '_'.join([channel, 'data', sys_name, distr_name])
    full_path = '%s/%s/%s/%s' % (channel, 'data', sys_name, histo_name)
    data_histo = fdata.Get(full_path)
    data_yields[channel] = data_histo.Integral()

    if not args.event_yields:
        print 'bin           ' + ''.join('%-15s' % channel for _ in processes)
        print 'process       ' + ''.join('%-15s' % proc for proc in processes)
        print 'process       ' + ''.join('%-15d' % processes_id[proc] for proc in processes)
        print 'rate          ' + ''.join('%-15.3f' % histo.Integral() for histo in histos)

        print full_path
        print 'obs %f' % data_histo.Integral()
        print 'mc sum = %f' % sum(h.Integral() for h in histos)


def string_yield(integral):
    if integral is None:
        return '%10s' % ''
    else:
        return '%10.1f' % integral

separator = ' $ & $ ' if args.event_yields == 'latex' else ''

if args.event_yields:
    print '%15s' % 'process' + separator.join(['%10s' % channel for channel in channels])
    #for process, chan_d in proc_yields.items():
    mc_sums = {} # channel: sum of integrals
    for process in all_known_sorted_processes:
        if process in proc_yields:
            chan_d = proc_yields[process]
            print '%15s' % process + separator.join([string_yield(chan_d.get(channel)) for channel in channels])
            # also calc sum
            for channel in channels:
                mc_sums.setdefault(channel, 0)
                mc_sums[channel] += chan_d.get(channel, 0)
        else:
            continue

    print '%15s' % 'mc_sum' + separator.join(['%10.1f' % mc_sums[channel] for channel in channels])
    print '%15s' % 'data'   + separator.join(['%10.1f' % data_yields[channel] for channel in channels])



