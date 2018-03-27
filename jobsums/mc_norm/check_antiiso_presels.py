import argparse
import logging
from os.path import isfile




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "check alliso presel contributions of qcds",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)


files = [
'MC2016_Summer16_QCD_HT-100-200.root',
'MC2016_Summer16_QCD_HT-200-300.root',
'MC2016_Summer16_QCD_HT-300-500.root',
'MC2016_Summer16_QCD_HT-500-700.root',
'MC2016_Summer16_QCD_HT-700-1000.root',
'MC2016_Summer16_QCD_HT-1000-1500.root',
'MC2016_Summer16_QCD_HT-1500-2000.root',
'MC2016_Summer16_QCD_HT-2000-Inf.root',
]





logging.info("import ROOT")

import ROOT
from ROOT import TH1D, TH2D, TFile, TTree, gROOT
from ROOT import TMath, TLorentzVector
from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
from ROOT import TLegend
#from ROOT import THStack


channel = 'ctr_old_mu_presel_alliso'
proc = 'qcd'
sys = 'NOMINAL'
distr = 'lep_relIso'

print "%60s    OS    SS" % "proc_file"
integral_os_all = 0
integral_ss_all = 0
for fin in files:
    tfile = TFile(fin)
    integral_os = tfile.Get("%s/%s/%s/" % (channel, proc, sys) + "%s_%s_%s_%s" % (channel, proc, sys, distr)).Integral()
    integral_ss = tfile.Get("%s/%s/%s/" % (channel + '_ss', proc, sys) + "%s_%s_%s_%s" % (channel + '_ss', proc, sys, distr)).Integral()
    integral_os_all += integral_os
    integral_ss_all += integral_ss
    print "%60s %f %f %f" % (fin, integral_os, integral_ss, integral_os/integral_ss)
print "%10s %f %f %f" % ('sum', integral_os_all, integral_ss_all, integral_os_all/ integral_ss_all)

print "%60s    OS    SS" % "per bin proc_file"
integral_os_isobin  = 0
integral_ss_isobin  = 0
integral_os_antibin = 0
integral_ss_antibin = 0
per_bin_anti = {}
for fin in files:
    tfile = TFile(fin)
    histo_os = tfile.Get("%s/%s/%s/" % (channel, proc, sys) + "%s_%s_%s_%s" % (channel, proc, sys, distr))
    histo_ss = tfile.Get("%s/%s/%s/" % (channel + '_ss', proc, sys) + "%s_%s_%s_%s" % (channel + '_ss', proc, sys, distr))

    isobin = 1
    integral_iso_os = histo_os.GetBinContent(isobin)
    integral_iso_ss = histo_ss.GetBinContent(isobin)

    integral_anti_os = 0.
    integral_anti_ss = 0.
    per_bin_anti_per_file = {}
    for bini in range(isobin, histo_os.GetSize()):
        reliso_yeild_os = histo_os.GetBinContent(bini)
        reliso_yeild_ss = histo_ss.GetBinContent(bini)

        if bini > isobin:
            integral_anti_os += reliso_yeild_os
            integral_anti_ss += reliso_yeild_ss
        #if integral_anti_ss > 0:
        #per_bin_anti.append((bini, integral_anti_os/ integral_anti_ss))
        if bini in per_bin_anti:
            per_bin_anti[bini][0] += reliso_yeild_os
            per_bin_anti[bini][1] += reliso_yeild_ss
        else:
            per_bin_anti[bini] = [reliso_yeild_os, reliso_yeild_ss]
        per_bin_anti_per_file[bini] = [reliso_yeild_os, reliso_yeild_ss]

    integral_os_isobin  += integral_iso_os
    integral_ss_isobin  += integral_iso_ss
    integral_os_antibin += integral_anti_os
    integral_ss_antibin += integral_anti_ss
    print "%60s  %f %f %f  |  %f %f %f" % (fin, integral_iso_os, integral_iso_ss, integral_iso_os/integral_iso_ss if integral_iso_ss > 0 else integral_iso_os,  integral_anti_os, integral_anti_ss, integral_anti_os/integral_anti_ss if integral_anti_ss > 0 else integral_anti_os)

    for bini, (os, ss) in per_bin_anti_per_file.items():
        print bini, os, ss, os/ss if ss > 0 else '-- %f' % os

print "%10s  %f %f %f   |   %f %f %f" % ('sum', integral_os_isobin, integral_ss_isobin, integral_os_isobin/ integral_ss_isobin, integral_os_antibin, integral_ss_antibin, integral_os_antibin/ integral_ss_antibin)

for bini, (os, ss) in per_bin_anti.items():
    print bini, os, ss, os/ss if ss > 0 else '-- %f' % os

print "iso"
channel = 'ctr_old_mu_presel'
integral_os_iso = 0
integral_ss_iso = 0
for fin in files:
    tfile = TFile(fin)
    integral_os = tfile.Get("%s/%s/%s/" % (channel, proc, sys) + "%s_%s_%s_%s" % (channel, proc, sys, distr)).Integral()
    integral_ss = tfile.Get("%s/%s/%s/" % (channel + '_ss', proc, sys) + "%s_%s_%s_%s" % (channel + '_ss', proc, sys, distr)).Integral()
    integral_os_iso += integral_os
    integral_ss_iso += integral_ss
    print "%60s %f %f %f" % (fin, integral_os, integral_ss, integral_os/integral_ss)
print "%10s %f %f %f" % ('sum', integral_os_iso, integral_ss_iso, integral_os_iso/ integral_ss_iso)


