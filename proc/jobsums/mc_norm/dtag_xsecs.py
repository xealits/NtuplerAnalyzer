
W_lep_br = 0.108;
W_qar_br = 0.676;

W_lep_br2 = W_lep_br*W_lep_br;
W_qar_br2 = W_qar_br*W_qar_br;

br_tau_electron = 0.1785;
br_tau_muon     = 0.1736;
br_tau_lepton   = br_tau_electron + br_tau_muon;
br_tau_hadronic = 1 - br_tau_lepton;

ttbar_xsec = 831.76;

dtag_xsecs = {
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo":  18610,
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo_": 18610,
"MC2016_Summer16_DYJetsToLL_50toInf_madgraph": 6225.42, #5765.4, # but McM shows 4970
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

"MC2016_Summer16_QCD_HT-100-200":   27540000    ,#* 0.0723985,
"MC2016_Summer16_QCD_HT-100-200_":  27540000    ,#* 0.0723985,
"MC2016_Summer16_QCD_HT-200-300":    1717000    ,#* 0.0967776,
"MC2016_Summer16_QCD_HT-200-300_":   1717000    ,#* 0.0967776,
"MC2016_Summer16_QCD_HT-300-500":     351300    ,#* 0.107775,
"MC2016_Summer16_QCD_HT-300-500_":    351300    ,#* 0.107775,
"MC2016_Summer16_QCD_HT-500-700":      31630    ,#* 0.141555,
"MC2016_Summer16_QCD_HT-500-700_":     31630    ,#* 0.141555,
"MC2016_Summer16_QCD_HT-700-1000":      6802    ,#* 0.1437,
"MC2016_Summer16_QCD_HT-700-1000_":     6802    ,#* 0.1437,
"MC2016_Summer16_QCD_HT-1000-1500":     1206    ,#* 0.160749,
"MC2016_Summer16_QCD_HT-1000-1500_":    1206    ,#* 0.160749,
"MC2016_Summer16_QCD_HT-1500-2000":      120.4  ,#* 0.141555,
"MC2016_Summer16_QCD_HT-1500-2000_":     120.4  ,#* 0.141555,
"MC2016_Summer16_QCD_HT-2000-Inf":        25.25 ,#* 0.135489,
"MC2016_Summer16_QCD_HT-2000-Inf_":       25.25 ,#* 0.135489,

"MC2016_Summer16_QCD_MuEnriched_Pt5_1000toInf"      : 10.4305   , # /  0.15544,
"MC2016_Summer16_QCD_MuEnriched_Pt5_1000toInf_ext1" : 10.4305   , # /  0.15544,
"MC2016_Summer16_QCD_MuEnriched_Pt5_15to20"         : 1273190000, # /  0.003,
"MC2016_Summer16_QCD_MuEnriched_Pt5_20to30"         : 558528000 , # /  0.0053,
"MC2016_Summer16_QCD_MuEnriched_Pt5_30to50"         : 139803000 , # /  0.01182,
"MC2016_Summer16_QCD_MuEnriched_Pt5_50to80"         :  19222500 , # /  0.02276,
"MC2016_Summer16_QCD_MuEnriched_Pt5_80to120"        :   2758420 , # /  0.03844,
"MC2016_Summer16_QCD_MuEnriched_Pt5_80to120_ext1"   :   2758420 , # /  0.03844,
"MC2016_Summer16_QCD_MuEnriched_Pt5_120to170"       : 469797    , # /  0.05362,
"MC2016_Summer16_QCD_MuEnriched_Pt5_170to300"       : 117989    , # /  0.07335,
"MC2016_Summer16_QCD_MuEnriched_Pt5_170to300_ext1"  : 117989    , # /  0.07335,
"MC2016_Summer16_QCD_MuEnriched_Pt5_300to470"       : 7820.25   , # /  0.10196,
"MC2016_Summer16_QCD_MuEnriched_Pt5_300to470_ext1"  : 7820.25   , # /  0.10196,
"MC2016_Summer16_QCD_MuEnriched_Pt5_300to470_ext2"  : 7820.25   , # /  0.10196,
"MC2016_Summer16_QCD_MuEnriched_Pt5_470to600"       : 645.528   , # /  0.12242,
"MC2016_Summer16_QCD_MuEnriched_Pt5_470to600_ext1"  : 645.528   , # /  0.12242,
"MC2016_Summer16_QCD_MuEnriched_Pt5_470to600_ext2"  : 645.528   , # /  0.12242,
"MC2016_Summer16_QCD_MuEnriched_Pt5_600to800"       : 187.109   , # /  0.13412,
"MC2016_Summer16_QCD_MuEnriched_Pt5_600to800_ext1"  : 187.109   , # /  0.13412,
"MC2016_Summer16_QCD_MuEnriched_Pt5_800to1000"      : 32.3486   , # /  0.14552,
"MC2016_Summer16_QCD_MuEnriched_Pt5_800to1000_ext1" : 32.3486   , # /  0.14552,
"MC2016_Summer16_QCD_MuEnriched_Pt5_800to1000_ext2" : 32.3486   , # /  0.14552,
"MC2016_Summer16_QCD_MuEnriched_Pt15_20toInf"       : 720648000 , # /  0.00042,

"MC2016_Summer16_QCD_EMEnriched_Pt-20to30"        : 557600000 * 0.0096  ,
"MC2016_Summer16_QCD_EMEnriched_Pt-30to50"        : 136000000 * 0.073   ,
"MC2016_Summer16_QCD_EMEnriched_Pt-30to50_ext1"   : 136000000 * 0.073   ,
"MC2016_Summer16_QCD_EMEnriched_Pt-50to80"        : 19800000  * 0.146   ,
"MC2016_Summer16_QCD_EMEnriched_Pt-50to80_ext1"   : 19800000  * 0.146   ,
"MC2016_Summer16_QCD_EMEnriched_Pt-80to120"       : 2800000   * 0.151   ,
"MC2016_Summer16_QCD_EMEnriched_Pt-80to120_ext1"  : 2800000   * 0.151   ,
"MC2016_Summer16_QCD_EMEnriched_Pt-120to170"      : 477000    * 0.162   ,
"MC2016_Summer16_QCD_EMEnriched_Pt-120to170_ext1" : 477000    * 0.162   ,
"MC2016_Summer16_QCD_EMEnriched_Pt-170to300"      : 114000    * 0.165   ,
"MC2016_Summer16_QCD_EMEnriched_Pt-300toInf"      : 9000      * 0.15    ,
"MC2016_Summer16_QCD_EMEnriched_Pt-30toInf"       : 162060000 * 0.0016  ,
"MC2016_Summer16_QCD_EMEnriched_Pt-30to40"        : 108000000 * 0.000225,
"MC2016_Summer16_QCD_EMEnriched_Pt-40toInf"       : 54120000  * 0.002   ,

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
"MC2016_Summer16_WJets_madgraph_reduced":   52940 - 9493 - 3120 - 942.3 - 524.2, # 50690 - 9493 - 3120 - 942.3 - 524.2},
"MC2016_Summer16_WJets_madgraph":   52940  , #50690  ,
"MC2016_Summer16_W0Jets_madgraph":  52940  , #50690  , # with the NUP
"MC2016_Summer16_W1Jets_madgraph":  9493   ,
"MC2016_Summer16_W2Jets_madgraph":  3120   ,
"MC2016_Summer16_W3Jets_madgraph":   942.3 ,
"MC2016_Summer16_W4Jets_madgraph":   524.2 ,

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

 "MC2016_Summer16_ZZTo2L2Nu_powheg"               : 0.564  ,

 "MC2016_Summer16_WWTo2L2Nu_powheg"               : 12.178  ,
 "MC2016_Summer16_WWToLNuQQ_powheg"               : 49.997  ,
 "MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin"   : 10.71  ,
 "MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin"     : 3.033  ,
 "MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin"      : 5.595  ,
 "MC2016_Summer16_WZTo3LNu_powheg"                : 4.42965  ,
 "MC2016_Summer16_ZZTo2L2Nu_powheg"               : 1.256  ,
 "MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin"      : 3.22   ,
"MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo_": 10.11, #3.36},
"MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo": 10.11, #3.36},
"MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg_": 80.95, #70.69/2},
"MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg": 80.95, #70.69/2},
"MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg_": 136.02, #70.69/2},
"MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg": 136.02, #70.69/2},
}

usual_gen_lumi_weights = {
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo" : 92667256.031132,
"MC2016_Summer16_DYJetsToLL_50toInf_madgraph" : 147827597.073751,

"MC2016_Summer16_QCD_HT-100-200" : 80668282.010224,
"MC2016_Summer16_QCD_HT-1000-1500" : 4831794.914995,
"MC2016_Summer16_QCD_HT-1500-2000" : 4025637.189128,
"MC2016_Summer16_QCD_HT-200-300" : 18982179.906991,
"MC2016_Summer16_QCD_HT-2000-Inf" : 2019120.954068,
"MC2016_Summer16_QCD_HT-300-500" : 17273854.173847,
"MC2016_Summer16_QCD_HT-500-700" : 19192242.333873,
"MC2016_Summer16_QCD_HT-700-1000" : 15846448.822448,

#"MC2016_Summer16_W1Jets_madgraph" : 42268723.655958,
#"MC2016_Summer16_W2Jets_madgraph" : 58623419.131360,
#"MC2016_Summer16_W3Jets_madgraph" : 57773140.789413,
#"MC2016_Summer16_W4Jets_madgraph" : 30412479.957944,

"MC2016_Summer16_W1Jets_madgraph" : 27024349.212653,
"MC2016_Summer16_W3Jets_madgraph" : 38799730.358834,
"MC2016_Summer16_W4Jets_madgraph" : 19897709.356113,
"MC2016_Summer16_W2Jets_madgraph" : 37639183.614431,
#"MC2016_Summer16_WJets_madgraph"  : 85737026.803850,
#"MC2016_Summer16_W0Jets_madgraph" : 85737026.803850,
"MC2016_Summer16_WJets_madgraph"  : 51952262.172543,
"MC2016_Summer16_W0Jets_madgraph" : 51952262.172543,

#"MC2016_Summer16_WJets_amcatnlo"  : 143451839.458732,
"MC2016_Summer16_WJets_amcatnlo"  : 151108470.560272,

"MC2016_Summer16_WWTo2L2Nu_powheg" : 2026759.028423,
"MC2016_Summer16_WWToLNuQQ_powheg" : 2026281.926199,
"MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin" : 14254071.531880,
"MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin" : 955200.742387,
"MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin" : 16101693.833498,
"MC2016_Summer16_WZTo3LNu_powheg" : 2020951.120024,
"MC2016_Summer16_ZZTo2L2Nu_powheg" : 8965538.308207,
"MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin" : 9822865.311869,

"MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo" : 637697.524850,
"MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg" : 39844542.965987,
"MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg" : 68953757.295413,
"MC2016_Summer16_SingleT_tW_5FS_powheg" : 7079946.591460,
"MC2016_Summer16_SingleTbar_tW_5FS_powheg" : 7060819.557982,

"MC2016_Summer16_TTJets_powheg" : 75829144.704647,
}

