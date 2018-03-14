from os.path import isfile
import argparse
import logging
logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "normalize histograms to luminosity",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument('-l', '--lumi', type=float, default=None, help="luminosity in pb to use instead of defaults")
parser.add_argument('-m', '--mu', action = "store_true", help="use SingleMuon full json lumi = 35867.060")
parser.add_argument('-e', '--el', action = "store_true", help="use SingleElectron reduced json lumi = 31341.774")
parser.add_argument('-f', '--files', type=str, default=None, help="list of filesnames to process, without .root separated by coma, if not given the default batch is used")

args = parser.parse_args()

control_channel = 'mu_lj'

if args.mu:
    lumi = 35867.060
    control_channel = 'mu_lj'
elif args.el:
    lumi = 31341.774
    control_channel = 'el_lj'

if args.lumi:
    lumi = args.lumi


logging.info("lumi = %f" % lumi)
logging.info("control_channel = %s" % control_channel)

logging.info("import ROOT")

import ROOT
from ROOT import TFile


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
#"MC2016_Summer16_DYJetsToLL_50toInf_madgraph_": 6025.2, // FIXME: the DY x-section was updated by about 5%
"MC2016_Summer16_DYJetsToLL_50toInf_madgraph_": 4970, # FIXME: the DY x-section was updated by about 5%
"MC2016_Summer16_QCD_HT-100-200":  27540000,
"MC2016_Summer16_QCD_HT-100-200_":  27540000,
"MC2016_Summer16_QCD_HT-1000-1500":  1206,
"MC2016_Summer16_QCD_HT-1000-1500_":  1206,
"MC2016_Summer16_QCD_HT-1500-2000":  120.4,
"MC2016_Summer16_QCD_HT-1500-2000_":  120.4,
"MC2016_Summer16_QCD_HT-200-300":  1717000,
"MC2016_Summer16_QCD_HT-200-300_":  1717000,
"MC2016_Summer16_QCD_HT-2000-Inf":  25.25,
"MC2016_Summer16_QCD_HT-2000-Inf_":  25.25,
"MC2016_Summer16_QCD_HT-300-500":  351300,
"MC2016_Summer16_QCD_HT-300-500_":  351300,
"MC2016_Summer16_QCD_HT-500-700":  31630,
"MC2016_Summer16_QCD_HT-500-700_":  31630,
"MC2016_Summer16_QCD_HT-700-1000":  6802,
"MC2016_Summer16_QCD_HT-700-1000_":  6802,
"MC2016_Summer16_SingleT_tW_5FS_powheg_":    35.6,
"MC2016_Summer16_SingleT_tW_5FS_powheg":    35.6,
"MC2016_Summer16_SingleTbar_tW_5FS_powheg_": 35.6,
"MC2016_Summer16_SingleTbar_tW_5FS_powheg": 35.6,
 "MC2016_Summer16_TTJets_powheg"      : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_1"      : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_2"      : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_3"      : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_4"      : ttbar_xsec ,
 "MC2016_Summer16_TTJets_powheg_5"      : ttbar_xsec ,
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
"MC2016_Summer16_W1Jets_madgraph": 0.65 *  9493,
"MC2016_Summer16_W1Jets_madgraph_": 9493,
"MC2016_Summer16_W2Jets_madgraph_perc": 2838.64,
"MC2016_Summer16_W2Jets_madgraph": 0.65 *  3120,
"MC2016_Summer16_W2Jets_madgraph_": 3120,
"MC2016_Summer16_W3Jets_madgraph_perc":  861.73,
"MC2016_Summer16_W3Jets_madgraph": 0.65 *  942.3,
"MC2016_Summer16_W3Jets_madgraph_": 942.3,
"MC2016_Summer16_W4Jets_madgraph_perc":  354.83,
"MC2016_Summer16_W4Jets_madgraph": 0.65 *  524.2,
"MC2016_Summer16_W4Jets_madgraph_": 524.2,
"MC2016_Summer16_WJets_amcatnlo": 61526.7 - 9493 - 3120 - 942.3 - 524.2,
"MC2016_Summer16_WJets_amcatnlo_full": 61526.7,
"MC2016_Summer16_WJets_amcatnlo_": 61526.7,
"MC2016_Summer16_WJets_madgraph_075": 50690 - 9493 - 3120 - 942.3 - 524.2,
"MC2016_Summer16_WJets_madgraph_perc": 37713.4 ,
"MC2016_Summer16_WJets_madgraph":  0.65 * 61526.7, # 50690 - 9493 - 3120 - 942.3 - 524.2},
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
"MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg",
# also overflows in data:
"SingleMuon",
"SingleElectron",
]

files = args.files.split(',') if args.files else files_default

print files



for fname in files:
    # check if file exists
    filename = fname + '.root'
    if not isfile(filename):
        logging.warning("file doesn't exist: %s" % filename)
        continue

    isMC = 'MC' in fname or fname[:3] == 'TT_'
    logging.info("isMC %d" % isMC)

    f = TFile(fname + '.root', "UPDATE")
    #weight_counter = f.Get('weight_counter')
    #weight_counter->GetBinContent(2)
    #scale = xsecs[fname] / weight_counter.GetBinContent(2)
    # norm to lumi mu+tau
    scale = lumi

    logging.info("%s %f" % (fname, scale))

    for proc in list(f.GetListOfKeys()):
        process = proc.GetName()
        if process in ('weight_counter', 'events_counter', 'control_counters'):
            continue

        # mu tracking factor
        mu_factor = 1.
        #if 'mu' in process:
        #    mu_factor = 0.975
        #    #logging.debug('mu channel, multiplying by %f for mu tracking' % mu_factor)

        #try:
        for chan in list(proc.ReadObj().GetListOfKeys()):
            nick = chan.GetName()
            if nick == 'none': #'tt_mutau':
                tauIDSF_factor = 0.90
            elif nick in ('tt_mutau3ch', 'tt_eltau3ch', 'tt_mutau', 'tt_eltau', 'tt_taultauh', 'dy_tautau', 's_top_eltau', 's_top_mutau', 'dy_tautau'):
                tauIDSF_factor = 0.95
            elif any(ch in nick for ch in ('ctr_old_mu_sel', 'ctr_old_mu_sel', 'ctr_old_el_sel', 'optel_tight_2L1M', 'optmu_tight_2L1M')):
                tauIDSF_factor = 0.95
            else:
                tauIDSF_factor = 1.

            for sys in list(chan.ReadObj().GetListOfKeys()):
                sys_name = sys.GetName()

                #logging.debug(sys_name)

                if 'alliso' in process and sys_name != 'NOMINAL':
                #if sys_name != 'NOMINAL':
                    continue # it takes forever

                # there is this 3% of error in PU right now.. correcting ad-hoc, TODO: re-proc jobs with new PU (frist get the new PU, then re-proc.. and add PU-effect-meter distr and other meters on Mt)

                #cor = 0.965 if 'mu_' in process else 1. # muons and electrons are a bit different eras........
                cor = 1.
                if 'PUUp' in sys_name:
                    pu_factor = cor * 1. / 0.9979 # 0.97 # 1./ 0.9979
                elif 'PUDown' in sys_name:
                    pu_factor = cor * 1. / 1.0485 # 1.17 # 1./ 1.485
                else:
                    pu_factor = cor * 1. / 1.022  # 1.06 # 1./ 1.02135 an 1/1.014 with weight counter..

                for histo_key in list(sys.ReadObj().GetListOfKeys()):
                    h = histo_key.ReadObj()
                    h.ClearUnderflowAndOverflow() # this might not work as expected! -- it leaves N entries the same
                    #h.Print()
                    #h.Sumw2() # to keep errors correctly
                    # it still complains that the structure is there

                    if process == control_channel and sys_name == 'NOMINAL' and histo_key.GetName() == '_'.join([process, nick, sys_name, 'Mt_lep_met']):
                        integral_init = h.Integral()
                        err_init = 0.
                        for i in range(h.GetSize()):
                            err_init += h.GetBinError(i)

                    h.Sumw2(ROOT.kTRUE) # to correctly save the errors after scaling
                    # it will raise Warnin for histograms which already have this -- have to filter it offline
                    if isMC:
                        h.Scale(scale * tauIDSF_factor * pu_factor * mu_factor)
                        if h.Integral() < 0:
                            h.Scale(-1.)

                    #h.Print()
                    if process == control_channel and sys_name == 'NOMINAL' and histo_key.GetName() == '_'.join([process, nick, sys_name, 'Mt_lep_met']):
                        err_scaled = 0.
                        for i in range(h.GetSize()):
                            err_scaled += h.GetBinError(i)
                        logging.debug("%s\t%s %10.f %5.5f %5.5f %3.3f %3.3f" % (fname, nick, scale, err_init/integral_init if integral_init > 0 else err_init, err_scaled/h.Integral() if h.Integral() > 0 else err_scaled, integral_init, h.Integral()))

        #except AttributeError as e:
            #continue #print e
            #print e

    f.Write()
    f.Close()

