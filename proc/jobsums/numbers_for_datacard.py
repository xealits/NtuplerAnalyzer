'''
bin         mu_lj     mu_lj     mu_lj     mu_lj         mu_lj       mu_lj       mu_lj    mu_lj        mu_lj     mu_lj     mu_lj
process     dy_tautau dy_other  wjets     s_top_mutau   s_top_other s_top_lj    dibosons tt_taultauh  tt_lj     tt_other  tt_mutau
process     10        9         8         7             6           5           4        3            2         1         0
rate        117.997   10.96     417.34    169.45        29.42       205.302     25.73    228.08       3523.18   328.59    3750.46     ------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

all_known_sorted_processes = [
    'tt_mutau', 'tt_eltau',
    'tt_elmu',  'tt_ltaul', 'tt_taueltaumu',
    'tt_taultauh',  'tt_lj', 'tt_taulj', 'tt_other',
    's_top_elmu',
    's_top_eltau',  's_top_mutau',   's_top_other', 's_top_lj',
    'dy_tautau', 'dy_other',
    'wjets_other', 'wjets_taul', # 'wjets_tauh',
    'dibosons',
    'qcd']

# tt_taulj
# wjets_taul
# wjets_tauh

processes_mu = ['dy_tautau', 'dy_other',  'wjets_other', 'wjets_taul', 'wjets_tauh',     's_top_mutau',   's_top_other', 's_top_lj',    'dibosons', 'tt_taultauh',  'tt_lj',     'tt_other',  'tt_mutau', 'tt_taulj', 'qcd']
#processes_mu = set(['dy_tautau', 'dy_other',  'wjets_other',     's_top_mutau',   's_top_other', 's_top_lj',    'dibosons', 'tt_taultauh',  'tt_lj',     'tt_other',  'tt_mutau'])
processes_mu_id = {'qcd': 14, 'dy_tautau': 13, 'dy_other': 12,  'wjets_other': 11,  'wjets_taul': 10,  'wjets_tauh': 9,     's_top_mutau': 8,   's_top_other': 7, 's_top_lj': 6,    'dibosons': 5,     'tt_other': 4, 'tt_taultauh': 3, 'tt_taulj': 2,  'tt_lj': 1,  'tt_mutau':0}

processes_el = ['dy_tautau', 'dy_other',  'wjets_other', 'wjets_taul', 'wjets_tauh',     's_top_eltau',   's_top_other', 's_top_lj',    'dibosons', 'tt_taultauh',  'tt_lj',     'tt_other',  'tt_eltau', 'tt_taulj', 'qcd']
#processes_el = set(['dy_tautau', 'dy_other',  'wjets_other',     's_top_eltau',   's_top_other', 's_top_lj',    'dibosons', 'tt_taultauh',  'tt_lj',     'tt_other',  'tt_eltau'])
processes_el_id = {'qcd': 14, 'dy_tautau': 13, 'dy_other': 12,  'wjets_other': 11,  'wjets_taul': 10,  'wjets_tauh': 9,     's_top_eltau': 8,   's_top_other': 7, 's_top_lj': 6,    'dibosons': 5,     'tt_other': 4, 'tt_taultauh': 3, 'tt_taulj': 2,  'tt_lj': 1,  'tt_eltau':0}

processes_id = {'qcd': 18, 'dy_tautau': 17, 'dy_other': 16,
   'wjets_other': 15,  'wjets_taul': 14,  'wjets_tauh': 13,
   's_top_eltau': 12, 's_top_mutau': 11, 's_top_elmu': 10,
   's_top_other': 9, 's_top_lj': 8,
   'dibosons': 7,
   'tt_other': 6, 'tt_taultauh': 5, 'tt_taulj': 4,  'tt_lj': 3,
   'tt_ltaul': 2, 'tt_taueltaumu': 1,
   'tt_eltau': 0, 'tt_mutau':0, 'tt_elmu':0}

"""
\text{dy}\rightarrow{\tau\tau}
\text{dy}\rightarrow{other}
\text{wjets}
\text{dibosons}
\text{top}^{single}\rightarrow{other}
\text{top}^{single}\rightarrow{\ell j}
\text{top}^{single}\rightarrow{\mu\tau}
t\bar{t}\rightarrow{other}
t\bar{t}\rightarrow{\tau_{l}\tau_{h}}
t\bar{t}\rightarrow{\ell j}
t\bar{t}\rightarrow{\mu\tau}
\text{mc}_{\text{sum}}
\text{data}
"""

process_latex_strings = {
'dy_tautau' : "$\\text{DY}\\rightarrow{\\tau_{\\ell}\\tau_{h}}+ {\\rm jets}$" ,
'dy_other'  : "$\\text{DY}\\rightarrow{\\rm other}$"                       ,

    'wjets_other':     "$\\text{W+jets} \\rightarrow\\ell\\nu+ {\\rm jets}$"         ,
    'wjets_taul':      "$\\text{W+jets} \\rightarrow\\tau_{\\ell}\\nu + {\\rm jets}$" ,
    'wjets_tauh':      "$\\text{W+jets} \\rightarrow\\tau_{h}\\nu + {\\rm jets}$"    ,

    's_top_eltau':   "$\\text{single top}\\rightarrow{{\\rm \\ell} \\tau}+X$" ,
    's_top_mutau':   "$\\text{single top}\\rightarrow{{\\rm \\ell} \\tau}+X$" ,
    's_top_other':   "$\\text{single top}\\rightarrow{\\rm other}$"         ,
    's_top_lj':      "$\\text{single top}\\rightarrow{\\ell j}+X$"          ,

    'tt_taultauh':  "$\\ttbar\\rightarrow{\\tau_{\\ell}\\tau_{h} \\nu\\nu b\\bar{b} }$",
    'tt_taulj':     "$\\ttbar\\rightarrow{\\tau_{\\ell}\\nu jj b\\bar{b}}$",
    'tt_lj':     "$\\ttbar\\rightarrow{\\ell\\nu jj b\\bar{b}}$",
    'tt_other':  "$\\ttbar\\rightarrow{\\rm other}$",
    'tt_mutau':  "$\\ttbar\\rightarrow{\\ell\\tau}\\nu\\nu b\\bar{b}$ ($\\ell=e,\\mu$)",
    'tt_eltau':  "$\\ttbar\\rightarrow{\\ell\\tau}\\nu\\nu b\\bar{b}$ ($\\ell=e,\\mu$)",

    'dibosons':                       '$\\text{dibosons}$',
    'qcd':                            '$\\text{QCD}$'
}



"""
'dy_tautau' : "$\text{DY}\rightarrow{\tau_{\ell}\tau_{h}}+ {\rm jets}$" ,
'dy_other'  : "$\text{DY}\rightarrow{\rm other}$"                       ,

    'wjets_other':     "$\text{W+jets} \rightarrow\ell\nu+ {\rm jets}$"         ,
    'wjets_taul':      "$\text{W+jets} \rightarrow\tau_{\ell}\nu + {\rm jets}$" ,
    'wjets_tauh':      "$\text{W+jets} \rightarrow\tau_{h}\nu + {\rm jets}$"    ,

    's_top_eltau':   "$\text{single top}\rightarrow{{\rm \ell} \tau}+X$" ,
    's_top_mutau':   "$\text{single top}\rightarrow{{\rm \ell} \tau}+X$" ,
    's_top_other':   "$\text{single top}\rightarrow{\rm other}$"         ,
    's_top_lj':      "$\text{single top}\rightarrow{\ell j}+X$"          ,

    'tt_taultauh':  "$\ttbar\rightarrow{\tau_{\ell}\tau_{h} \nu\nu b\bar{b} }$",
    'tt_taulj':     "$\ttbar\rightarrow{\tau_{\ell}\nu jj b\bar{b}}$",
    'tt_lj':     "$\ttbar\rightarrow{\ell\nu jj b\bar{b}}$",
    'tt_other':  "$\ttbar\rightarrow{\rm other}$",
    'tt_mutau':  "$\ttbar\rightarrow{\ell\tau}\nu\nu b\bar{b}$ ($\ell=e,\mu$)",
    'tt_eltau':  "$\ttbar\rightarrow{\ell\tau}\nu\nu b\bar{b}$ ($\ell=e,\mu$)",

    'dibosons':                       '$\\text{dibosons}$',
    'qcd':                            '$\\text{QCD}$'
}
"""


import argparse


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "extract observed and expected values in el_lj/el_lj_out or mu_lj/mu_lj_out channels",
    epilog = """Example:
python numbers_for_datacard.py histosel-roots/v25p3_yields/histosel_v25p3_yields.root -y -d yield -c ctr_old_mu_presel,ctr_old_mu_sel,ctr_old_mu_sel_lj,ctr_old_mu_sel_ljout --latex"""
    )

parser.add_argument("data_file",    help="full Data+MC file name")
parser.add_argument("-c", "--channels",  type=str, help="if given then get the numbers from here, enter multiple channels with 'mu_presel,mu_sel'")
parser.add_argument("-d", "--distr",    type=str, default='Mt_lep_met', help="fitted distribution")
parser.add_argument("-s", "--sys",    type=str, default='NOMINAL', help="systematic shift to consider")
parser.add_argument("-m", "--mu",  action = "store_true", help="muon channels: mu_lj, mu_lj_out")
parser.add_argument("-e", "--el",  action = "store_true", help="electron channels: el_lj, el_lj_out")
parser.add_argument("--debug",     action = "store_true", help="log debug output")
parser.add_argument("--no-sums",   action = "store_false", default=True, help="don't collect MC sum")
#parser.add_argument("-y", "--event-yields", type=str, default='text', help="output in the form of event yield table (set -y latex for latex table output)")
parser.add_argument("-y", "--event-yields", action='store_true', help="output in the form of event yield table (set -y latex for latex table output)")
parser.add_argument("-r", "--ratios",       action='store_true', help="output ratios to data")

parser.add_argument("--skip-procs",    type=str, default='', help="skip these processes")

parser.add_argument("--wjets",    type=float, help="factor of wjets")

parser.add_argument("--range-min",    type=float, help="calculate from")
parser.add_argument("--range-max",    type=float, help="calculate up to")

parser.add_argument("--latex",       action='store_true', help="output in latex")


args = parser.parse_args()

import logging

if args.debug:
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

logging.debug("data file " + args.data_file)

fdata = TFile(args.data_file)

sys_name   = args.sys
distr_name = args.distr

logging.debug("%s %s" % (sys_name, distr_name))

#if args.event_yields:
#proc_yields = [] #[('channel', [(proc_nick, yield)...]) ...]
proc_yields = {} # proc_nick: {channel: yield}

data_yields = {} # channel: yield

def range_integral(histo):
    if args.range_max or args.range_min:
        integral = 0.
        uncertainty = 0.
        for bini in range(histo.GetSize()):
            bin_x   = histo.GetBinCenter(bini)
            if args.range_max and bin_x > args.range_max:
                break
            if args.range_min and bin_x < args.range_min:
                continue
            integral    += histo.GetBinContent(bini)
            uncertainty += histo.GetBinError(bini)
    else:
        integral = histo.Integral()
        uncertainty = 0.
        for bini in range(histo.GetSize()):
            bin_x   = histo.GetBinCenter(bini)
            uncertainty += histo.GetBinError(bini)

    return integral, uncertainty

# sumhistos in the channel
mc_sums_sumhisto = {}

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

    # get sum from the sum histo
    # sums_NOMINAL/mc_sum2_NOMINAL
    if args.no_sums:
        sum_histo = chan.Get("sums_NOMINAL/mc_sum2_NOMINAL")
        mc_sums_sumhisto[channel] = range_integral(sum_histo)

    histos = []
    processes = []

    #if 'mu' in channel:
    #    processes = processes_mu
    #    processes_id = processes_mu_id
    #else:
    #    processes = processes_el
    #    processes_id = processes_el_id

    #proc_yields.append((channel, []))

    #processes_keys = list(chan.ReadObj().GetListOfKeys())
    #processes_keys = list(chan.GetListOfKeys())
    #sorted_pkeys = sorted(processes_keys, key=lambda pkey: pkey.GetName() if pkey.GetName() not in ('qcd') else 'z_' + pkey.GetName())
    #for process in sorted_pkeys:
    #for process in processes:
    for process in (p.GetName() for p in list(chan.GetListOfKeys())):
       # if the channel with shapes for dy and wjets are given
       # don't use other channels
       #if shape_channel and nick in ('dy', 'wjets') and chan.GetName() != shape_channel:
       #    continue
       logging.debug(process)

       if process not in all_known_sorted_processes:
           continue

       if args.skip_procs and process in args.skip_procs:
           continue

       histo_name = '_'.join([channel, process, sys_name, distr_name])

       logging.debug(histo_name)
       histo_name = '%s/%s/%s' % (process, sys_name, histo_name)
       logging.debug(histo_name)

       histo = chan.Get(histo_name)
       if not histo:
           logging.info('no %s' % histo_name)
           continue
       logging.debug(histo.Integral())
       processes.append(process)

       if 'wjet' in process and args.wjets:
           histo.Scale(args.wjets)

       histos.append(histo)

       #proc_yields[-1][1].append((process, histo.Integral()))
       proc_yields.setdefault(process, {})[channel] = range_integral(histo)

    histo_name = '_'.join([channel, 'data', 'NOMINAL', distr_name])
    full_path = '%s/%s/%s/%s' % (channel, 'data', 'NOMINAL', histo_name)
    data_histo = fdata.Get(full_path)

    logging.debug(full_path)
    data_yields[channel] = range_integral(data_histo)

    if not args.event_yields:
        processes.sort(key=lambda k: processes_id[k])
        print 'bin           ' + ''.join('%-25s' % channel for _ in processes)
        print 'process       ' + ''.join('%-25s' % proc for proc in processes)
        print 'process       ' + ''.join('%-25d' % processes_id[proc] for proc in processes)
        print 'rate          ' + ''.join('%-25.3f' % (proc_yields[proc][channel][0] if proc in proc_yields else 0) for proc in processes)

        print full_path
        #print 'obs %f' % data_histo.Integral()
        #print 'mc sum = %f' % sum(h.Integral() for h in histos)
        print 'obs %f'      % data_yields[channel][0]
        print 'mc sum = %f' % sum(proc[channel][0] for proc in proc_yields.values())


proc_s = '%40s'
item_s = '%10'

def string_yield(integral):
    if integral is None:
        return (item_s + 's') % ''
    else:
        if args.ratios:
            return (item_s + '.3f') % integral
        else:
            return (item_s + '.0f') % integral

separator = ' & ' if args.latex else ''
line_end = ' \\\\' if args.latex else ''

if args.event_yields:
    print proc_s % 'process' + separator + separator.join([(item_s + 's') % channel for channel in channels]) + line_end
    #for process, chan_d in proc_yields.items():
    mc_sums   = {} # channel: sum of integrals

    for process in all_known_sorted_processes:
        if process in proc_yields:
            chan_d = proc_yields[process]
            line = proc_s % (process_latex_strings[process] if args.latex else process)
            # check the process in each requested channel
            # process = row
            # channel = column
            for channel in channels:
                #
                if not channel in chan_d: continue
                integral  , uncertainty      = chan_d.get(channel)
                data_yield, uncertainty_data = data_yields[channel]
                if integral and args.ratios:
                    integral = integral / data_yield
                line += separator + '$' + string_yield(integral) + ' \\pm ' + string_yield(uncertainty) + ' $'
            print line + line_end
            # also calc sum
            for channel in channels:
                mc_sums.setdefault(channel, [0, 0])
                mc_sums[channel][0] += chan_d.get(channel, (0, 0))[0]
                mc_sums[channel][1] += chan_d.get(channel, (0, 0))[1]
        else:
            continue


    if args.no_sums:
        if args.ratios:
            print proc_s % 'mc_sum_ratio' + separator + separator.join([(item_s + '.3f') % (mc_sums[channel][0]/data_yields[channel][0]) for channel in channels]) + line_end
        else:
            #print proc_s % 'mc_sum' + separator + separator.join([('$' + item_s + '.1f' + ' \\pm ' + item_s + '.1f $') % tuple(mc_sums[channel]) for channel in channels]) + line_end
            print proc_s % 'mc_sum' + separator + separator.join([('$' + item_s + '.0f' + ' \\pm ' + item_s + '.0f $') % tuple(mc_sums_sumhisto[channel]) for channel in channels]) + line_end

    print proc_s % 'data' + separator + separator.join([('$' + item_s + '.0f' + ' \\pm ' + item_s + '.0f $') % data_yields[channel] for channel in channels]) + line_end



