import argparse
import logging
import os
from os.path import isfile, basename
import ctypes
from array import array
from sys import exit
from collections import OrderedDict
from math import sqrt



parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "signal acceptance printout",
    epilog = """Example:
python print_acceptances.py temp/NtuplerAnalyzer_test_METfiltersOFF_TTJets2_signal_recordAll.root
python print_acceptances.py --sys-weights lstore_outdirs/merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg*.root --devs merge
python print_acceptances.py --sys-weights lstore_outdirs/merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg*.root --devs merge_files_rel --with-stat --how-print cols
python print_acceptances.py --sys-weights lstore_outdirs/merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg*.root --devs merge_files
"""
    )

parser.add_argument("--debug",       action='store_true', help="DEBUG level of logging")

#parser.add_argument("input_filename",  type=str, help="the path tot he file with acceptances")
parser.add_argument("--sys-weights", action='store_true', help="print weight-based systematics too")
parser.add_argument("--yields",      action='store_true', help="print number of events processed")
#parser.add_argument("--ratio",       action='store_true', help="calculate ratio to the first value")
parser.add_argument("--all-procs",   action='store_true', help="print for all processes")
parser.add_argument("--with-stat",   action='store_true', help="print stat error too")
parser.add_argument("--devs",        type=str, help="calculate the relative systematic deviations")
parser.add_argument("--how-print",   type=str, default='rows', help="the default output is in rows, if any argument is given print columns")

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

sum_channel = OrderedDict()

# nicknames of the systematic files
filenames = []

# weight names to indexes
systs_pdfs = ['PDFCTn%s' % str(n) for n in range(1, 57)]
systs_scales = ['scale_renorm_Up', 'scale_renorm_Down', 'scale_refrag_Up', 'scale_refrag_Down', 'scale_comb_Up', 'scale_comb_Down']
systs_scales_renorm = ['scale_renorm_Up', 'scale_renorm_Down']
systs_scales_refact = ['scale_refrag_Up', 'scale_refrag_Down']
systs_scales_both   = ['scale_comb_Up',   'scale_comb_Down']

systs_alphas = ['alpha_Up', 'alpha_Down']

#systs_used_groups = [('scales', systs_scales), ('alphas', systs_alphas), ('pdfs', systs_pdfs)]
#systs_used_groups = [('scales_renorm', systs_scales_renorm), ('scales_refact', systs_scales_refact), ('scales_both', systs_scales_both), ('alphas', systs_alphas), ('pdfs', systs_pdfs)]
systs_used_groups = [('scales_renorm', systs_scales_renorm), ('scales_refact', systs_scales_refact), ('scales_both', systs_scales_both), ('pdfs', systs_pdfs + systs_alphas)]

sys_weght_name_ind = {n+1: name for n, name in enumerate(['nom'] + [sys_name for _, snames in systs_used_groups for sys_name in snames])}
range_length = (len(sys_weght_name_ind) + 1) if args.sys_weights else 2

# nicknames of the nominal and systematic files
nominal_file_n = 0

systematic_filenames = {'CUETP8M2T4down': 'tune_Down',
'CUETP8M2T4up': 'tune_Up',
'hdampUP':   'hdamp_Up',
'hdampDOWN': 'hdamp_Down',
'fsrup':   'fsr_Up',
'fsrdown': 'fsr_Down',
'isrup':   'isr_Up',
'isrdown': 'isr_Down'}


for file_i, input_filename in enumerate(args.input_files):
    if not isfile(input_filename):
        logging.info("missing: " + input_filename)
        exit(1)

    logging.debug(input_filename)

    file_nickname = None
    for fname in systematic_filenames:
        if fname in input_filename:
            file_nickname = systematic_filenames[fname]
            break

    if not file_nickname:
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
    for name, file_numbers in results.items():
        histo_all = tfile.Get("all_ev_" + name)
        histo_cut = tfile.Get("cut_ev_" + name)
        histo_cut.Divide(histo_all)

        file_numbers[file_nickname] = OrderedDict()

        if args.yields:
            file_numbers[file_nickname]['event_yield'] = (histo_all.GetBinContent(1), histo_all.GetBinError(1))

        for i in range(1,range_length):
            weight_name = sys_weght_name_ind[i]
            file_numbers[file_nickname][weight_name] = (histo_cut.GetBinContent(i), histo_cut.GetBinError(i))

    # save the acceptance of the sum of the processes (usefull only for signal processes!)
    sum_channel[file_nickname] = OrderedDict()
    if args.yields:
        sum_channel[file_nickname]['event_yield'] = (sum_all_histos.GetBinContent(1), sum_all_histos.GetBinError(1))

    for i in range(1,range_length):
        weight_name = sys_weght_name_ind[i]
        sum_channel[file_nickname][weight_name]   = (sum_cut_histos.GetBinContent(i), sum_cut_histos.GetBinError(i))

results['sum'] = sum_channel

# calculate relative deviations of systematics from nominals
# the statistical uncertainty of deviations in each nominal file does not matter
# it's relevant only for the difference different files (different events)
if args.devs:
    for channel, file_results in results.items():
        for file_nickname, numbers in file_results.items():
            nominal, _ = numbers['nom']

            # when merging systematic groups, construct new dict with results
            if 'merge' in args.devs:
                update_numbers = OrderedDict()
                update_numbers['nom'] = numbers['nom']
                for sys_group_name, sys_group in systs_used_groups:
                    if any(sys_name in numbers for sys_name in sys_group):
                        if 'rel' in args.devs:
                            devs = [(sys_val-nominal)/nominal for sname, (sys_val, _) in numbers.items() if sname in sys_group]
                        else:
                            devs = [sys_val-nominal for sname, (sys_val, _) in numbers.items() if sname in sys_group]
                        logging.debug(devs)
                        sqrt_sum_dev = sqrt(sum(dev**2 for dev in devs))
                        logging.debug(sqrt_sum_dev)
                        update_numbers[sys_group_name] = sqrt_sum_dev, 0

                file_results[file_nickname] = update_numbers

            else:
                for sys_name, (sys_val, sys_err) in numbers.items():
                    if sys_name == 'nom': continue
                    if 'rel' in args.devs:
                        numbers[sys_name] = ((sys_val-nominal)/nominal, sys_err/nominal)
                    else:
                        numbers[sys_name] = (sys_val-nominal, sys_err)

if 'merge_files' in args.devs:
    for channel, file_results in results.items():
        # get the nominal file:
        nominal, _ = file_results['nominal_0']['nom']

        # for all other files get and merge the Up-Down nominal deviations
        merged_files = OrderedDict()
        for sys_name, sys_group_names in [('tune', ('tune_Up', 'tune_Down')), ('hdamp', ('hdamp_Up', 'hdamp_Down')), ('tune', ('tune_Up', 'tune_Down')), ('isr', ('isr_Up', 'isr_Down')), ('fsr', ('fsr_Up', 'fsr_Down'))]:
            sys_group = [file_results[fn]['nom'] for fn in sys_group_names if fn in sys_group_names]
            if 'rel' in args.devs:
                sys_dev  = sqrt(sum(((nominal - sys_val) / nominal)**2 for sys_val, _ in sys_group))
                sys_stat = sqrt(sum((sys_sta / nominal)**2 for _, sys_sta in sys_group))/len(sys_group)
            else:
                sys_dev  = sqrt(sum((nominal - sys_val)**2 for sys_val, _ in sys_group))
                sys_stat = sqrt(sum((sys_sta)**2 for _, sys_sta in sys_group))/len(sys_group)
            merged_files.setdefault(sys_name, OrderedDict())['nom'] = sys_dev, sys_stat

        # pop the separate up-downs deviations from non-nominal files:
        for fn in file_results:
            if fn == 'nominal_0': continue
            file_results.pop(fn)

        # update with the merged devs
        file_results.update(merged_files)

## print stuff in rows
if args.how_print == 'rows':
    # get the sysnames from any first channel saved in results
    sys_names = [' %20s' % ' '.join((fname, sysname)) for fname, fsysts in results.values()[0].items() for sysname in fsysts]
    print "%20s " % 'channel' + ' '*5, ' '.join(sys_names)
    for channel, file_numbers in results.items():
        numbers = [v for f in file_numbers.values() for v in f.values()]
        #
        if args.with_stat:
            print "%20s " % channel + ' '*5, ' '.join(" %10.4f +- %6.4f" % (val, unc) for val, unc in numbers)
        else:
            print "%20s " % channel + ' '*5, ' '.join(" %20.4f"          % val for val, _ in numbers)

# print stuff in columns per process
else:
    # revert it for syst per row
    per_sys = OrderedDict()
    for pn, numbers in results.items():
        for file_nickname in numbers:
          for sys_name, val_unc in numbers[file_nickname].items():
            per_sys.setdefault((file_nickname, sys_name), OrderedDict())[pn] = val_unc

    proc_names = results.keys()
    logging.debug(proc_names)

    print ' '*21 + ' '*5 + ' '.join(['%10s' % p for p in proc_names])
    for name, numbers in per_sys.items():
        #
        logging.debug(name)
        logging.debug(repr(numbers))
        if args.with_stat:
            print ("%10s %-20s " % name) + ' '*5 + ' '.join(["%10.4f +- %6.4f" % numbers[pn] for pn in proc_names])
        else:
            print ("%10s %-20s " % name) + ' '*5 + ' '.join(["%10.4f" % numbers[pn][0] for pn in proc_names])

## TODO: what is this ratio for? it does not work now
#if args.ratio:
#    for name, numbers in results.items():
#        print "%20s " % name + ' '*5, numbers[0][0] / (sum(numbers[1:][0])/(len(numbers) - 1))

