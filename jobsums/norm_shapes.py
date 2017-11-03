import argparse
import logging
logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "stack histos",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument("mc_file",    help="MC file name")
#parser.add_argument("-x", "--shape", type=str, default='', help="selection for shape distributions of wjets and dy")
#parser.add_argument("-m", "--muon",     action = "store_true", help="don't save the root file, plot stacked histogram in png")
#parser.add_argument("-e", "--electron", action = "store_true", help="don't save the root file, make ratio plot (in addition to stack or alone)")

'''
so I need to take histo scale it and save in another directory inside a file, overwriting the histo there.
'''

args = parser.parse_args()

logging.info("import ROOT")

import ROOT
from ROOT import TFile, THStack, TLegend, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
import plotting_root


infile = TFile(args.mc_file, 'UPDATE')

channels = ['mu_lj', 'mu_lj_out', 'el_lj', 'el_lj_out']
shapes   = ['mu_presel', 'mu_presel', 'el_presel', 'el_presel']

processes = ['dy_tautau', 'dy_other', 'wjets']

distr_names = ['Mt_lep_met', 'Mt_lep_met_f']

infile.cd()

for chname, shname in zip(channels, shapes):

    for proc in processes:
        systs_dir = infile.Get('%s/%s' % (chname, proc)).GetListOfKeys()

        for sys in systs_dir:
          for distr_name in distr_names:
            sys_name = sys.GetName()
            histo_name = '_'.join([chname, proc, sys_name, distr_name])
            histo_shape_name = '_'.join([shname, proc, sys_name, distr_name])

            logging.info(histo_name)

            histo       = infile.Get('%s/%s/%s/%s' % (chname, proc, sys_name, histo_name))
            histo_shape = infile.Get('%s/%s/%s/%s' % (shname, proc, sys_name, histo_shape_name))

            histo_shape.Scale(histo.Integral() / histo_shape.Integral())

            #histo = histo_shape.Clone()
            #histo.SetName(histo_name)

            #histo_dir = infile.Get('%s/%s/%s' % (chname, proc, sys_name))
            #histo_dir.cd()
            #histo.SetDirectory(histo_dir)
            # ROOT is nonsense

            sum_histo_err, sum_shape_err = 0., 0.
            for i in range(histo.GetSize()):
                sum_histo_err += histo.GetBinError(i)
                sum_shape_err += histo_shape.GetBinError(i)

            error_ratio = sum_histo_err / sum_shape_err

            for i in range(histo.GetSize()):
                histo_err = histo.GetBinError(i)
                shape_err = histo_shape.GetBinError(i)
                histo.SetBinContent(i, histo_shape.GetBinContent(i))
                #if shape_err > histo_err:
                #    histo.SetBinError(i, shape_err)
                histo.SetBinError(i, shape_err * error_ratio)


infile.Write()
infile.Close()

