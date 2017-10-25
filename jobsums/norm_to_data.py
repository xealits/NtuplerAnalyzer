import logging
logging.basicConfig(level=logging.DEBUG)

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


fname = 'mc2_norm'

f = TFile(fname + '.root', "UPDATE")
#weight_counter = f.Get('weight_counter')
#weight_counter->GetBinContent(2)
#scale = xsecs[fname] / weight_counter.GetBinContent(2)
# norm to lumi mu+tau
scale = 35600.

logging.info("%s %f" % (fname, scale))

for proc in list(f.GetListOfKeys()):
    try:
      for chan in list(proc.ReadObj().GetListOfKeys()):
        for sys in list(chan.ReadObj().GetListOfKeys()):
            for histo_key in list(sys.ReadObj().GetListOfKeys()):
                h = histo_key.ReadObj()
                #h.Print()
                h.Scale(scale)
                #h.Print()
    except:
        pass

f.Write('', TFile.kOverwrite)
f.Close()


