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

args = parser.parse_args()



logging.info("import ROOT")

import ROOT
from ROOT import TFile
#from ROOT.TObject import kOverwrite


W_lep_br = 0.108;
W_qar_br = 0.676;

W_lep_br2 = W_lep_br*W_lep_br;
W_qar_br2 = W_qar_br*W_qar_br;

br_tau_electron = 0.1785;
br_tau_muon     = 0.1736;
br_tau_lepton   = br_tau_electron + br_tau_muon;
br_tau_hadronic = 1 - br_tau_lepton;

ttbar_xsec = 831.76;

xsecs = {
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo": 18610,
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo_": 18610,
"MC2016_Summer16_DYJetsToLL_50toInf_madgraph": 5765.4, # but McM shows 4970
#"MC2016_Summer16_DYJetsToLL_50toInf_madgraph": 6025.2, # KIT TauTau uses this one FIXME: the DY x-section was updated by about 5%
# link: https://github.com/cms-analysis/HiggsAnalysis-KITHiggsToTauTau/blob/704c2adb0148f62580c68d60dc309075639ff6e1/data/ArtusConfigs/Includes/settingsCrossSection.json
# -- it's DY amcatnlo
# Spring15
# commit https://github.com/cms-analysis/HiggsAnalysis-KITHiggsToTauTau/commit/bbd27735feecd7f56e8c390e25301fe0b1c16bb2
# updated GluGluJets
# tag SMHTT2016 CMSSW_7_4_7
# FIXME: "the DY x-section was updated by about 5%" -- so, need to check what's the official value after all, confirm 5765
#"MC2016_Summer16_DYJetsToLL_50toInf_madgraph_": 6025.2,
"MC2016_Summer16_DYJetsToLL_50toInf_madgraph_": 4970,

"MC2016_Summer16_QCD_HT-100-200":   27540000    * 0.0723985,
"MC2016_Summer16_QCD_HT-100-200_":  27540000    * 0.0723985,
"MC2016_Summer16_QCD_HT-200-300":    1717000    * 0.0967776,
"MC2016_Summer16_QCD_HT-200-300_":   1717000    * 0.0967776,
"MC2016_Summer16_QCD_HT-300-500":     351300    * 0.107775,
"MC2016_Summer16_QCD_HT-300-500_":    351300    * 0.107775,
"MC2016_Summer16_QCD_HT-500-700":      31630    * 0.141555,
"MC2016_Summer16_QCD_HT-500-700_":     31630    * 0.141555,
"MC2016_Summer16_QCD_HT-700-1000":      6802    * 0.1437,
"MC2016_Summer16_QCD_HT-700-1000_":     6802    * 0.1437,
"MC2016_Summer16_QCD_HT-1000-1500":     1206    * 0.160749,
"MC2016_Summer16_QCD_HT-1000-1500_":    1206    * 0.160749,
"MC2016_Summer16_QCD_HT-1500-2000":      120.4  * 0.141555,
"MC2016_Summer16_QCD_HT-1500-2000_":     120.4  * 0.141555,
"MC2016_Summer16_QCD_HT-2000-Inf":        25.25 * 0.135489,
"MC2016_Summer16_QCD_HT-2000-Inf_":       25.25 * 0.135489,

"MC2016_Summer16_SingleT_tW_5FS_powheg_":    35.6,
"MC2016_Summer16_SingleT_tW_5FS_powheg":    35.6,
"MC2016_Summer16_SingleTbar_tW_5FS_powheg_": 35.6,
"MC2016_Summer16_SingleTbar_tW_5FS_powheg": 35.6,

'TT_TuneCUETP8M2T4_13TeV-powheg-fsrdown-pythia8'     : ttbar_xsec,
'TT_TuneCUETP8M2T4_13TeV-powheg-fsrup-pythia8'       : ttbar_xsec,
'TT_TuneCUETP8M2T4_13TeV-powheg-isrdown-pythia8'     : ttbar_xsec,
'TT_TuneCUETP8M2T4_13TeV-powheg-isrup-pythia8'       : ttbar_xsec,
'TT_TuneCUETP8M2T4down_13TeV-powheg-pythia8'         : ttbar_xsec,
'TT_TuneCUETP8M2T4up_13TeV-powheg-pythia8'           : ttbar_xsec,
'TT_hdampDOWN_TuneCUETP8M2T4_13TeV-powheg-pythia8'   : ttbar_xsec,
'TT_hdampUP_TuneCUETP8M2T4_13TeV-powheg-pythia8'     : ttbar_xsec,

'MC2016_Summer16_TTJets_powheg_isrup'           : ttbar_xsec,
'MC2016_Summer16_TTJets_powheg_isrdown'         : ttbar_xsec,
'MC2016_Summer16_TTJets_powheg_hdampUP'         : ttbar_xsec,
'MC2016_Summer16_TTJets_powheg_hdampDOWN'       : ttbar_xsec,
'MC2016_Summer16_TTJets_powheg_fsrup'           : ttbar_xsec,
'MC2016_Summer16_TTJets_powheg_fsrdown'         : ttbar_xsec,
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4up'    : ttbar_xsec,
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4down'  : ttbar_xsec,

"MC2016_Summer16_TTJets_powheg_FSRDown": ttbar_xsec,
"MC2016_Summer16_TTJets_powheg_FSRUp":   ttbar_xsec,
"MC2016_Summer16_TTJets_powheg_ISRDown": ttbar_xsec,
"MC2016_Summer16_TTJets_powheg_ISRUp":   ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg"     : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_add" : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_1"   : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_2"   : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_3"   : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_4"   : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_5"   : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_2_aelmtuu" :  ttbar_xsec * W_lep_br2 * 2 * br_tau_lepton ,
 "MC2016_Summer16_TTJets_powheg_2_elmu"    :  ttbar_xsec * W_lep_br2 * 2 ,
 "MC2016_Summer16_TTJets_powheg_2_other"   :  ttbar_xsec * (1 - 2*W_lep_br2 - 2*W_lep_br2*br_tau_lepton) ,
 "MC2016_Summer16_TTJets_powheg_aattuu"      : ttbar_xsec * W_lep_br2 ,
 "MC2016_Summer16_TTJets_powheg_aeltu"       : ttbar_xsec * W_lep_br2 * 2 ,
 "MC2016_Summer16_TTJets_powheg_amtuu"       : ttbar_xsec * W_lep_br2 * 2 ,
 "MC2016_Summer16_TTJets_powheg_aqtu"        : ttbar_xsec * W_lep_br*W_qar_br * 2,
 "MC2016_Summer16_TTJets_powheg_aaelmttuuu":  ttbar_xsec * W_lep_br2 *2 * br_tau_muon     * br_tau_electron,
 "MC2016_Summer16_TTJets_powheg_aaeellttuu":  ttbar_xsec * W_lep_br2    * br_tau_electron * br_tau_electron ,
 "MC2016_Summer16_TTJets_powheg_aaehlttuu" :  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_electron ,
 "MC2016_Summer16_TTJets_powheg_aammttuuuu":  ttbar_xsec * W_lep_br2    * br_tau_muon     * br_tau_muon ,
 "MC2016_Summer16_TTJets_powheg_aahmttuuu" :  ttbar_xsec * W_lep_br2 *2 * br_tau_hadronic * br_tau_muon ,
 "MC2016_Summer16_TTJets_powheg_aahhttuu"  :  ttbar_xsec * W_lep_br2    * br_tau_hadronic * br_tau_hadronic,
 "MC2016_Summer16_TTJets_powheg_aelmtuu"   :  ttbar_xsec * W_lep_br2         * 2 * br_tau_lepton ,
 "MC2016_Summer16_TTJets_powheg_aeelltu"   :  ttbar_xsec * W_lep_br2         * 2 * br_tau_electron ,
 "MC2016_Summer16_TTJets_powheg_ammtuuu"   :  ttbar_xsec * W_lep_br2         * 2 * br_tau_muon ,
 "MC2016_Summer16_TTJets_powheg_aelqtu"    :  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_electron ,
 "MC2016_Summer16_TTJets_powheg_amqtuu"    :  ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_muon ,
 "MC2016_Summer16_TTJets_powheg_aehltu"     : ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic,
 "MC2016_Summer16_TTJets_powheg_ahmtuu"     : ttbar_xsec * W_lep_br2 * 2         * br_tau_hadronic,
 "MC2016_Summer16_TTJets_powheg_ahqtu"      : ttbar_xsec * W_lep_br*W_qar_br * 2 * br_tau_hadronic,
 "MC2016_Summer16_TTJets_powheg_eell"      : ttbar_xsec * W_lep_br2 ,
 "MC2016_Summer16_TTJets_powheg_elmu"      : ttbar_xsec * W_lep_br2 * 2 ,
 "MC2016_Summer16_TTJets_powheg_elq"       : ttbar_xsec * W_lep_br*W_qar_br * 2,
 "MC2016_Summer16_TTJets_powheg_mmuu"      : ttbar_xsec * W_lep_br2 ,
 "MC2016_Summer16_TTJets_powheg_mqu"       : ttbar_xsec * W_lep_br*W_qar_br * 2,
 "MC2016_Summer16_TTJets_powheg_qq"        : ttbar_xsec * W_qar_br2 ,
 "MC2016_Summer16_TTJets_powheg_aattuu"      : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aeltu"       : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_amtuu"       : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aqtu"        : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aaelmttuuu":  ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aaeellttuu":  ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aaehlttuu" :  ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aammttuuuu":  ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aahmttuuu" :  ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aahhttuu"  :  ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aelmtuu"   :  ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aeelltu"   :  ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_ammtuuu"   :  ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aelqtu"    :  ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_amqtuu"    :  ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aehltu"     : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_ahmtuu"     : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_ahqtu"      : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_eell"      : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_elmu"      : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_elq"       : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_mmuu"      : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_mqu"       : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_qq"        : ttbar_xsec,
 "MC2016_Summer16_TTJets_powheg_aattuu"      : ttbar_xsec * W_lep_br2 ,
 "MC2016_Summer16_TTJets_powheg_aeltu"       : ttbar_xsec * W_lep_br2 * 2 ,
 "MC2016_Summer16_TTJets_powheg_amtuu"       : ttbar_xsec * W_lep_br2 * 2 ,
 "MC2016_Summer16_TTJets_powheg_aqtu"        : ttbar_xsec * W_lep_br*W_qar_br ,
 "MC2016_Summer16_TTJets_powheg_eell"        : ttbar_xsec * W_lep_br2 ,
 "MC2016_Summer16_TTJets_powheg_elmu"        : ttbar_xsec * W_lep_br2 * 2 ,
 "MC2016_Summer16_TTJets_powheg_elq"         : ttbar_xsec * W_lep_br*W_qar_br ,
 "MC2016_Summer16_TTJets_powheg_mmuu"        : ttbar_xsec * W_lep_br2 ,
 "MC2016_Summer16_TTJets_powheg_mqu"         : ttbar_xsec * W_lep_br*W_qar_br ,
 "MC2016_Summer16_TTJets_powheg_qq"          : ttbar_xsec * W_qar_br2 ,

"MC2016_Summer16_W1Jets_madgraph_perc": 8921.44,
"MC2016_Summer16_W1Jets_madgraph_": 9493,
"MC2016_Summer16_W2Jets_madgraph_perc": 2838.64,
"MC2016_Summer16_W2Jets_madgraph_": 3120,
"MC2016_Summer16_W3Jets_madgraph_perc":  861.73,
"MC2016_Summer16_W3Jets_madgraph_": 942.3,
"MC2016_Summer16_W4Jets_madgraph_perc":  354.83,
"MC2016_Summer16_W4Jets_madgraph_": 524.2,

"MC2016_Summer16_WJets_amcatnlo":      61526.7, #49143.25,# norm study in madgraph gives about 49k x-sec # 61526.7 - 9493 - 3120 - 942.3 - 524.2, # = 47447.2
"MC2016_Summer16_WJets_amcatnlo_full": 61526.7,
"MC2016_Summer16_WJets_amcatnlo_":     61526.7,
"MC2016_Summer16_WJets_madgraph_075":  50690 - 9493 - 3120 - 942.3 - 524.2,
"MC2016_Summer16_WJets_madgraph_perc": 37713.4 ,

## TODO: need to replay this study..
## o.65 factor
#"MC2016_Summer16_WJets_madgraph":  61526.7 - 9493 - 3120 - 942.3 - 524.2, # 50690 - 9493 - 3120 - 942.3 - 524.2},
"MC2016_Summer16_WJets_madgraph":  50690 - 9493 - 3120 - 942.3 - 524.2, # 50690 - 9493 - 3120 - 942.3 - 524.2},
"MC2016_Summer16_W0Jets_madgraph":  50690, # with the NUP
"MC2016_Summer16_W1Jets_madgraph":  9493,
"MC2016_Summer16_W2Jets_madgraph":  3120,
"MC2016_Summer16_W3Jets_madgraph":   942.3,
"MC2016_Summer16_W4Jets_madgraph":   524.2,

# with the MiT factors
# https://github.com/MiT-HEP/ChargedHiggs/blob/8927a1101b87bd0e4bb8059a0cdd4e957eef08e7/aux/DYNjets.txt
# -- also had it in some report and something similar in a phd
### WNjets MG ## xSec=61526.7
# 0.557053 0.299967 0.0981022 0.0296045 0.0152728 0 0 0 0
# 61526.7 .* [0.557053 0.299967 0.0981022 0.0296045 0.0152728] = 34273.6  18456.0  6035.9  1821.47  939.685

#"MC2016_Summer16_WJets_madgraph":  34273.6,
#"MC2016_Summer16_W0Jets_madgraph": 34273.6,
#"MC2016_Summer16_W1Jets_madgraph": 18456.0,
#"MC2016_Summer16_W2Jets_madgraph":  6035.9,
#"MC2016_Summer16_W3Jets_madgraph":  1821.47,
#"MC2016_Summer16_W4Jets_madgraph":   939.685,

 "MC2016_Summer16_WWTo2L2Nu_powheg"               :  12.178  ,
 "MC2016_Summer16_WWToLNuQQ_powheg"               :  49.997  ,
 "MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin"   : 10.71  ,
 "MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin"     : 3.033  ,
 "MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin"      : 5.595  ,
 "MC2016_Summer16_WZTo3LNu_powheg"                : 4.42965  ,
 "MC2016_Summer16_ZZTo2L2Nu_powheg"               : 0.564  ,
 "MC2016_Summer16_ZZTo2L2Nu_powheg"               : 1.256  ,
 "MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin"      : 3.22   ,
"MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo_": 10.11, #3.36},
"MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo": 10.11, #3.36},
"MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg_": 80.95, #70.69/2},
"MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg": 80.95, #70.69/2},
"MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg_": 136.02, #70.69/2},
"MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg": 136.02, #70.69/2},
}

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
    else:
        event_weight = 1.0

    # to include the "shapeness" of various th weight-base systematics
    # since v25 now event_weight = nominal_w * PU_rate * sys_rate
    # nominal_w is just amcatanlo
    # PU rate and other sys rates are calculated as PU_weight / nominal_events <- IT DOES NOT INCLUDE -1 aMCatNLO!
    # ok, in v25p1 leave at as is
    #event_weight = 
    scale = xsecs[fname] / event_weight


    logging.info("%s %f %f" % (fname, scale, xsecs[fname]))

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


