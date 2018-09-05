from os.path import isfile
import argparse


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "find the systematic variations of the yield in a channel for given process",
    epilog = """Example:
python yield_variations.py mc_v25v26pR5_updowns_with_data.root -c ctr_old_mu_sel_lj -p tt_taultauh,tt_lj,tt_other,tt_mutau,tt_taulj -s FSR,ISR,HDAMP,TuneCUETP8M2T4"""
    )

parser.add_argument("data_file",    help="full Data+MC file name")
parser.add_argument("-c", "--channels",  type=str, help="the channels to look in")
parser.add_argument("-p", "--processes", type=str, default='tt_taultauh,tt_lj,tt_other,tt_mutau,tt_taulj', help="processes to consider")
parser.add_argument("-d", "--distr",     type=str, default='yield', help="the distribution (yield by default)")
parser.add_argument("-s", "--sys",       type=str, default='FSR,ISR,HDAMP,TuneCUETP8M2T4', help="systematics to consider")

parser.add_argument("--debug",      action = "store_true", help="log debug output")
parser.add_argument("--to-average", action = "store_true", help="ratio to average of systematics, instead of the nominal tt")


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

assert isfile(args.data_file)

fdata = TFile(args.data_file)

channels  = []
processes = []
for channel in args.channels.split(','):
    for process in args.processes.split(','):
        channels.append(channel)
        processes.append(process)
print ' '*20 + '  ' + '  '.join(('%20s' % chan for chan in channels))
print ' '*20 + '  ' + '  '.join(('%20s' % proc for proc in processes))

for sys in args.sys.split(','):
    sys_variations = []
    for channel in args.channels.split(','):
        for process in args.processes.split(','):
            #print process

            sysname_up = "{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}".format(chan=channel, proc=process, sys=sys+'Up',   distr=args.distr)
            #print sysname_up
            systematic_up   = fdata.Get(sysname_up)
            systematic_down = fdata.Get("{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}".format(chan=channel, proc=process, sys=sys+'Down', distr=args.distr))
            if not systematic_up:
                print "no " + sys + 'Up'
                sys_up_yield   = 0.
            else:
                sys_up_yield   = systematic_up.Integral()

            if not systematic_down:
                print "no " + sys + 'Down'
                sys_down_yield = 0.
            else:
                sys_down_yield = systematic_down.Integral()

            if args.to_average:
                reference = (sys_up_yield + sys_down_yield) * 0.5
            else:
                nominal = fdata.Get("{chan}/{proc}/NOMINAL/{chan}_{proc}_NOMINAL_{distr}".format(chan=channel, proc=process, distr=args.distr))
                reference = nominal.Integral()

            logging.debug('%s %s %s %f %f %f' % (channel, process, sys, sys_up_yield, sys_down_yield, reference))
            if reference == 0:
                logging.warning('reference == 0 in %s %s %s' % (channel, process, sys))
                continue

            var_up   = sys_up_yield   / reference
            var_down = sys_down_yield / reference
            sys_variations.append((var_up, var_down))

    print ('%20s  ' % sys) + '  '.join(('%20s' % ('%.3f/%.3f' % (var_down, var_up)) for var_up, var_down in sys_variations))


