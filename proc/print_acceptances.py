import argparse
import logging
import os
from os.path import isfile, basename
import ctypes
from array import array
from sys import exit
from collections import OrderedDict



parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "signal acceptance printout",
    epilog = """Example:
python print_acceptances.py temp/NtuplerAnalyzer_test_METfiltersOFF_TTJets2_signal_recordAll.root
"""
    )

#parser.add_argument("input_filename",  type=str, help="the path tot he file with acceptances")
parser.add_argument("--sys-weights", action='store_true', help="print weight-based systematics too")
parser.add_argument("--yields",      action='store_true', help="print number of events processed")
#parser.add_argument("--ratio",       action='store_true', help="calculate ratio to the first value")
parser.add_argument("--all-procs",   action='store_true', help="print for all processes")
parser.add_argument("--with-stat",   action='store_true', help="print stat error too")
parser.add_argument("--debug",       action='store_true', help="DEBUG level of logging")

parser.add_argument('input_files', nargs="+", help="""files with acceptances""")


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

if args.all_procs:
    results = {
    "mutau"    : OrderedDict(),
    "eltau"    : OrderedDict(),
    "muj"      : OrderedDict(),
    "elj"      : OrderedDict(),
    "taumutauh": OrderedDict(),
    "taueltauh": OrderedDict(),
    "taumuj"   : OrderedDict(),
    "tauelj"   : OrderedDict(),
    "other"    : OrderedDict()}
else:
    results = {
    "mutau"    : OrderedDict(),
    "eltau"    : OrderedDict()}

# nicknames of the systematic files
filenames = []

# weight names to indexes
pdf_names = ['PDFCTn%s' % str(n) for n in range(1, 57)]
sys_weght_name_ind = {n+1: name for n, name in enumerate(['nom', 'scale_renorm_Up', 'scale_renorm_Down', 'scale_refrag_Up', 'scale_refrag_Down', 'scale_comb_Up', 'scale_comb_Down', 'alpha_Up', 'alpha_Down'] + pdf_names)}

range_length = (len(sys_weght_name_ind) + 1) if args.sys_weights else 2

# nicknames of the nominal and systematic files
nominal_file_n = 0

for file_i, input_filename in enumerate(args.input_files):
    if not isfile(input_filename):
        logging.info("missing: " + input_filename)
        exit(1)

    logging.debug(input_filename)

    if   'CUETP8M2T4down' in input_filename:
        file_nickname = 'tune_Down'
    elif 'CUETP8M2T4up'   in input_filename:
        file_nickname = 'tune_Up'
    elif 'hdampUP'   in input_filename:
        file_nickname = 'hdampUP'
    elif 'hdampDOWN' in input_filename:
        file_nickname = 'hdampDOWN'
    elif 'fsrup'   in input_filename:
        file_nickname = 'fsrup'
    elif 'fsrdown' in input_filename:
        file_nickname = 'fsrdown'
    elif 'isrup'   in input_filename:
        file_nickname = 'isrup'
    elif 'isrdown' in input_filename:
        file_nickname = 'isrdown'
    else:
        file_nickname = 'nominal_' + str(nominal_file_n)
        nominal_file_n += 1

    tfile = TFile(input_filename)

    all_histos = [tfile.Get("all_ev_" + name) for name in results]
    cut_histos = [tfile.Get("cut_ev_" + name) for name in results]

    sum_all_histos = all_histos[0].Clone()
    sum_cut_histos = cut_histos[0].Clone()

    sum_all_histos.SetName("sum_all_histos_%s" % file_nickname)
    sum_cut_histos.SetName("sum_cut_histos_%s" % file_nickname)

    # calculate the sums of all processes
    for histo_cut, histo_all in zip(cut_histos[1:], all_histos[1:]):
        logging.debug('adding the sum')
        sum_all_histos += histo_all
        sum_cut_histos += histo_cut

    sum_cut_histos.Divide(sum_all_histos)

    # divide for acceptance and save all processes
    for name, numbers in results.items():
        histo_all = tfile.Get("all_ev_" + name)
        histo_cut = tfile.Get("cut_ev_" + name)
        histo_cut.Divide(histo_all)

        if args.yields:
            numbers[(file_nickname, 'event_yield')] = (histo_all.GetBinContent(1), histo_all.GetBinError(1))

        for i in range(1,range_length):
            weight_name = sys_weght_name_ind[i]
            numbers[(file_nickname, weight_name)] = (histo_cut.GetBinContent(i), histo_cut.GetBinError(i))

    # save the acceptance of the sum of the processes (usefull only for signal processes!)
    results['sum'] = OrderedDict()
    if args.yields:
        results['sum'][(file_nickname, 'event_yield')] = (sum_all_histos.GetBinContent(1), sum_all_histos.GetBinError(1))

    for i in range(1,range_length):
        weight_name = sys_weght_name_ind[i]
        results['sum'][(file_nickname, weight_name)]   = (sum_cut_histos.GetBinContent(i), sum_cut_histos.GetBinError(i))

## print stuff in rows
#for name, numbers in results.items():
#    #
#    if args.with_stat:
#        print "%20s " % name + ' '*5, ' '.join("%10.4f +- %6.4f" % (val, unc) for val, unc in numbers.values())
#    else:
#        print "%20s " % name + ' '*5, ' '.join("%10.4f" % val for val, _ in numbers.values())

# print stuff in columns per process

# revert it for syst per row
per_sys = OrderedDict()
for pn, numbers in results.items():
    for sys_name, val_unc in numbers.items():
        per_sys.setdefault(sys_name, OrderedDict())[pn] = val_unc

proc_names = results.keys()

print ' '*21 + ' '*5 + ' '.join(['%10s' % p for p in proc_names])
for name, numbers in per_sys.items():
    #
    if args.with_stat:
        print ("%10s %-20s " % name) + ' '*5 + ' '.join(["%10.4f +- %6.4f" % numbers[pn] for pn in proc_names])
    else:
        #print str(numbers)
        #print [numbers[pn] for pn in proc_names]
        #print [numbers[pn][0] for pn in proc_names]
        #print [val for pn in proc_names for val in numbers[pn]]
        ##print [val for pn in proc_names for val, _ in numbers[pn]]
        ##print ["%10.4f" % val for pn in proc_names for val, _ in numbers[pn]]
        #print ["%10.4f" % numbers[pn][0] for pn in proc_names]
        #print ' '.join(["%10.4f" % numbers[pn][0] for pn in proc_names])
        #print name
        #print ("%20s " % name)
        #print ("%20s " % name) + ' '*5
        print ("%10s %-20s " % name) + ' '*5 + ' '.join(["%10.4f" % numbers[pn][0] for pn in proc_names])

## TODO: what is this ratio for? it does not work now
#if args.ratio:
#    for name, numbers in results.items():
#        print "%20s " % name + ' '*5, numbers[0][0] / (sum(numbers[1:][0])/(len(numbers) - 1))

