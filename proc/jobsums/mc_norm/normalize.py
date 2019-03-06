from os.path import isfile
import argparse
import logging
logging.basicConfig(level=logging.DEBUG)



parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "normalize histograms to the cross sections of processes",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument('-f', '--files', type=str, default=None, help="list of filesnames to process, without .root separated by coma, if not given the default batch is used")
parser.add_argument('--per-weight', action='store_true', help='normalize to event weight')
parser.add_argument('--per-usual-weight', action='store_true', help='normalize to usual event weight, recorded in v37 test1')

args = parser.parse_args()



logging.info("import ROOT")

import ROOT
from ROOT import TFile
#from ROOT.TObject import kOverwrite

from dtag_xsecs import dtag_xsecs, usual_gen_lumi_weights

files_default = [
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo",
"MC2016_Summer16_DYJetsToLL_50toInf_madgraph",
"MC2016_Summer16_QCD_HT-100-200",
"MC2016_Summer16_QCD_HT-1000-1500",
"MC2016_Summer16_QCD_HT-1500-2000",
"MC2016_Summer16_QCD_HT-200-300",
"MC2016_Summer16_QCD_HT-2000-Inf",
"MC2016_Summer16_QCD_HT-300-500",
"MC2016_Summer16_QCD_HT-500-700",
"MC2016_Summer16_QCD_HT-700-1000",
"MC2016_Summer16_SingleT_tW_5FS_powheg",
"MC2016_Summer16_SingleTbar_tW_5FS_powheg",
"MC2016_Summer16_TTJets_powheg",
"MC2016_Summer16_TTJets_powheg_FSRDown",
"MC2016_Summer16_TTJets_powheg_FSRUp",
"MC2016_Summer16_TTJets_powheg_ISRDown",
"MC2016_Summer16_TTJets_powheg_ISRUp",
"MC2016_Summer16_W1Jets_madgraph",
"MC2016_Summer16_W2Jets_madgraph",
"MC2016_Summer16_W3Jets_madgraph",
"MC2016_Summer16_W4Jets_madgraph",
"MC2016_Summer16_WJets_amcatnlo",
"MC2016_Summer16_WJets_madgraph",
"MC2016_Summer16_WWTo2L2Nu_powheg",
"MC2016_Summer16_WWToLNuQQ_powheg",
"MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin",
"MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin",
"MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin",
"MC2016_Summer16_WZTo3LNu_powheg",
"MC2016_Summer16_ZZTo2L2Nu_powheg",
"MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin",
"MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo",
"MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg",
"MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg"]


files = args.files.split(',') if args.files else files_default

print files

for fname in files:
    # check if file exists
    filename = fname + '.root'
    if not isfile(filename):
        logging.warning("file doesn't exist: %s" % filename)
        continue

    f = TFile(filename, "UPDATE")
    #f = TFile(fname + '.root')
    weight_counter = f.Get('weight_counter')
    #weight_counter->GetBinContent(2)
    if args.per_weight:
        event_weight = weight_counter.GetBinContent(2)
    elif args.per_usual_weight:
        if fname in usual_gen_lumi_weights:
            event_weight = usual_gen_lumi_weights[fname]
        else:
            print "no usual gen lumi weight for %s" % fname
            event_weight = 1.0
    else:
        event_weight = 1.0

    # to include the "shapeness" of various th weight-base systematics
    # since v25 now event_weight = nominal_w * PU_rate * sys_rate
    # nominal_w is just amcatanlo
    # PU rate and other sys rates are calculated as PU_weight / nominal_events <- IT DOES NOT INCLUDE -1 aMCatNLO!
    # ok, in v25p1 leave at as is
    #event_weight = 
    scale = dtag_xsecs[fname] / event_weight


    logging.info("%s %f %f" % (fname, scale, dtag_xsecs[fname]))

    for proc in list(f.GetListOfKeys()):
        process = proc.GetName()
        if process in ('weight_counter', 'events_counter', 'control_counters', 'systematic_weights'):
            continue
        #try:
        for chan in list(proc.ReadObj().GetListOfKeys()):
            #if process == 'mu_presel':
                #logging.debug(str(chan))
            nick = chan.GetName()
            for sys in list(chan.ReadObj().GetListOfKeys()):
                sys_name = sys.GetName()

                if 'alliso' in process and sys_name != 'NOMINAL':
                #if sys_name != 'NOMINAL':
                    continue # it takes forever

                for histo_key in list(sys.ReadObj().GetListOfKeys()):
                    h = histo_key.ReadObj()
                    integral_init = h.Integral()

                    h.Sumw2(ROOT.kTRUE) # to correctly save the errors after scaling
                    # it complains with
                    # TH1D::Sumw2:0: RuntimeWarning: Sum of squares of weights structure already created
                    # but apparently it fixes stuff in wjets
                    # -- filter the ROOT's stderr on commandline
                    # check if root is loosing the errors somewhere here
                    if process == 'mu_sel' and sys.GetName() == 'NOMINAL' and histo_key.GetName() == '_'.join([process, nick, sys_name, 'Mt_lep_met']):
                        err_init = 0.
                        sum_of_weights_init = h.GetSumOfWeights()
                        for i in range(h.GetSize()):
                            err_init += h.GetBinError(i)

                    h.Scale(scale)
                    #h.Print()
                    if process == 'mu_sel' and sys.GetName() == 'NOMINAL' and histo_key.GetName() == '_'.join([process, nick, sys_name, 'Mt_lep_met']):
                        err_scaled = 0.
                        for i in range(h.GetSize()):
                            err_scaled += h.GetBinError(i)
                        #print proc.GetName(), chan.GetName(), sys.GetName(), histo_key.GetName()
                        #logging.debug("%s %s %10.f %5.5f %5.5f" % (h.GetName(), nick, event_weight, integral_init, h.Integral()))
                        logging.debug("%s\t%s %2.6f %3.2f %10.f  %3.3f / %5.5f  %3.3f / %5.5f %3.3f %3.3f %3.3f %3.3f, %5.5f  %5.5f" % (fname, nick, scale, (scale*35000), event_weight, err_init, integral_init, err_scaled, h.Integral(), err_init/integral_init if integral_init > 0 else err_init, err_scaled/h.Integral() if h.Integral() > 0 else err_scaled, h.GetSumOfWeights(), sum_of_weights_init, integral_init, h.Integral()))

                    # root Update doesn't actually update but makes new clone of the object...
                    #h.Write(h.GetName(), TFile.kOverwrite)
                    # now these copies are at top of the file...
        #except:
            #print "failed"

    f.Write("", TFile.kOverwrite)
    f.Close()


