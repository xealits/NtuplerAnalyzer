import os
import argparse
import logging


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "just run sumup draw with standard histo names",
    epilog = """Example:\npython per_dtag_script.py MC2016_Summer16_DYJetsToLL_50toInf_madgraph"""
    )

parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")
parser.add_argument('input_dtags', nargs='+', help="the dtags to process")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

dtags = {
'MC2016_Summer16_DYJetsToLL_10to50_amcatnlo'   : ('dy', 18610),
'MC2016_Summer16_DYJetsToLL_50toInf_madgraph'  : ('dy', 5765.4),
'MC2016_Summer16_QCD_HT-100-200'    : ('qcd', 27540000    * 0.0723985),
'MC2016_Summer16_QCD_HT-200-300'    : ('qcd',  1717000    * 0.0967776),
'MC2016_Summer16_QCD_HT-300-500'    : ('qcd',   351300    * 0.107775 ),
'MC2016_Summer16_QCD_HT-500-700'    : ('qcd',    31630    * 0.141555 ),
'MC2016_Summer16_QCD_HT-700-1000'   : ('qcd',     6802    * 0.1437   ),
'MC2016_Summer16_QCD_HT-1000-1500'  : ('qcd',     1206    * 0.160749 ),
'MC2016_Summer16_QCD_HT-1500-2000'  : ('qcd',      120.4  * 0.141555 ),
'MC2016_Summer16_QCD_HT-2000-Inf'   : ('qcd',       25.25 * 0.135489 ),
'MC2016_Summer16_QCD_MuEnriched_Pt15_20toInf'  : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_1000toInf' : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_120to170'  : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_15to20'    : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_170to300'  : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_20to30'    : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_300to470'  : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_30to50'    : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_470to600'  : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_50to80'    : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_600to800'  : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_800to1000' : ('qcd_muenriched', 1.),
'MC2016_Summer16_QCD_MuEnriched_Pt5_80to120'   : ('qcd_muenriched', 1.),

'MC2016_Summer16_QCD_EMEnriched_Pt-120to170.root' : ('qcd_emenriched', 1.),
'MC2016_Summer16_QCD_EMEnriched_Pt-170to300.root' : ('qcd_emenriched', 1.),
'MC2016_Summer16_QCD_EMEnriched_Pt-20to30.root'   : ('qcd_emenriched', 1.),
'MC2016_Summer16_QCD_EMEnriched_Pt-300toInf.root' : ('qcd_emenriched', 1.),
'MC2016_Summer16_QCD_EMEnriched_Pt-30to40.root'   : ('qcd_emenriched', 1.),
'MC2016_Summer16_QCD_EMEnriched_Pt-30to50.root'   : ('qcd_emenriched', 1.),
'MC2016_Summer16_QCD_EMEnriched_Pt-30toInf.root'  : ('qcd_emenriched', 1.),
'MC2016_Summer16_QCD_EMEnriched_Pt-40toInf.root'  : ('qcd_emenriched', 1.),
'MC2016_Summer16_QCD_EMEnriched_Pt-50to80.root'   : ('qcd_emenriched', 1.),
'MC2016_Summer16_QCD_EMEnriched_Pt-80to120.root'  : ('qcd_emenriched', 1.),

'MC2016_Summer16_SingleT_tW_5FS_powheg'                     : ('single_top', 35.6),
'MC2016_Summer16_SingleTbar_tW_5FS_powheg'                  : ('single_top', 35.6),
'MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo'      : ('single_top', 10.11),
'MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg' : ('single_top', 80.95),
'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg'     : ('single_top', 136.02),
'MC2016_Summer16_TTJets_powheg'  : ('ttbar', 831.76),

#'MC2016_Summer16_TTJets_powheg_CUETP8M2T4down' : genprocs_mutau_tt,
#'MC2016_Summer16_TTJets_powheg_CUETP8M2T4up'   : genprocs_mutau_tt,
#'MC2016_Summer16_TTJets_powheg_fsrdown'        : genprocs_mutau_tt,
#'MC2016_Summer16_TTJets_powheg_fsrup'          : genprocs_mutau_tt,
#'MC2016_Summer16_TTJets_powheg_hdampDOWN'      : genprocs_mutau_tt,
#'MC2016_Summer16_TTJets_powheg_hdampUP'        : genprocs_mutau_tt,
#'MC2016_Summer16_TTJets_powheg_isrdown'        : genprocs_mutau_tt,
#'MC2016_Summer16_TTJets_powheg_isrup'          : genprocs_mutau_tt,

'MC2016_Summer16_WWTo2L2Nu_powheg'             : ('dibosons',  12.178 ),
'MC2016_Summer16_WWToLNuQQ_powheg'             : ('dibosons',  49.997 ),
'MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin' : ('dibosons',  10.71  ),
'MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin'   : ('dibosons',  3.033  ),
'MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin'    : ('dibosons',  5.595  ),
'MC2016_Summer16_WZTo3LNu_powheg'              : ('dibosons',  4.42965),
'MC2016_Summer16_ZZTo2L2Nu_powheg'             : ('dibosons',  1.256  ),
'MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin'    : ('dibosons',  3.22   ),

'MC2016_Summer16_WJets_amcatnlo'   : ('wjets', 61526.7),
'MC2016_Summer16_WJets_madgraph'   : ('wjets', 50690 ),
'MC2016_Summer16_W1Jets_madgraph'  : ('wjets', 9493  ),
'MC2016_Summer16_W2Jets_madgraph'  : ('wjets', 3120  ),
'MC2016_Summer16_W3Jets_madgraph'  : ('wjets',  942.3),
'MC2016_Summer16_W4Jets_madgraph'  : ('wjets',  524.2),

'SingleMuon'     : ('data', 1.),
'SingleElectron' : ('data', 1.),
}

std_dtags = [
'MC2016_Summer16_QCD_HT-100-200',
'MC2016_Summer16_QCD_HT-200-300',
'MC2016_Summer16_QCD_HT-300-500',
'MC2016_Summer16_QCD_HT-500-700',
'MC2016_Summer16_QCD_HT-700-1000',
'MC2016_Summer16_QCD_HT-1000-1500',
'MC2016_Summer16_QCD_HT-1500-2000',
'MC2016_Summer16_QCD_HT-2000-Inf',
'MC2016_Summer16_WWTo2L2Nu_powheg',
'MC2016_Summer16_WWToLNuQQ_powheg',
'MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin',
'MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin',
'MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin',
'MC2016_Summer16_WZTo3LNu_powheg',
'MC2016_Summer16_ZZTo2L2Nu_powheg',
'MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin',
'MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo',
'MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg',
'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg',
'MC2016_Summer16_DYJetsToLL_10to50_amcatnlo',
'MC2016_Summer16_DYJetsToLL_50toInf_madgraph',
'MC2016_Summer16_SingleT_tW_5FS_powheg',
'MC2016_Summer16_SingleTbar_tW_5FS_powheg',
'MC2016_Summer16_W1Jets_madgraph',
'MC2016_Summer16_W2Jets_madgraph',
'MC2016_Summer16_W3Jets_madgraph',
'MC2016_Summer16_W4Jets_madgraph',
'MC2016_Summer16_WJets_amcatnlo',
'MC2016_Summer16_WJets_madgraph',
'SingleElectron',
'SingleMuon',
]

chan = 'dy_mutau_sel'
sys = 'NOMINAL'
distr = 'geom_tau_sv_sign'

muon_chan = False
channel, pt_cut = ('abs(leps_ID) == 13 && HLT_mu', 26) if muon_chan else ('abs(leps_ID) == 13 && HLT_mu', 29)

com = """python sumup_ttree_draw.py "tau_SV_geom_flightLenSign[tau_refited_index[0]]" --histo-range 21,-1,20 --cond-com "%s && lep_matched_HLT[0] && tau_IDlev[0] > 2 && tau_refited_index[0] > -1 && lep_p4[0].pt() > %d && abs(lep_p4[0].eta()) < 2.4 && tau_p4[0].pt() > 20 && abs(tau_p4[0].eta()) < 2.4 && nbjets < 2 && lep_id[0]*tau_id[0] < 0"     {options}        --save-weight --histo-name "{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}" --output {outfile} gstore_outdirs/v40/*/Ntupler_v40_{dtag}*/*/0000/*root """ % (channel, pt_cut)

# with Mt < 40
com = """python sumup_ttree_draw.py "tau_SV_geom_flightLenSign[tau_refited_index[0]]" --histo-range 21,-1,20 --cond-com "%s && lep_matched_HLT[0] && tau_IDlev[0] > 2 && tau_refited_index[0] > -1 && lep_p4[0].pt() > %d && abs(lep_p4[0].eta()) < 2.4 && tau_p4[0].pt() > 20 && abs(tau_p4[0].eta()) < 2.4 && nbjets < 2 && lep_id[0]*tau_id[0] < 0 && sqrt(2*(sqrt((met_init.Px()*met_init.Px() + met_init.Py()*met_init.Py())*(lep_p4[0].Px()*lep_p4[0].Px() + lep_p4[0].Py()*lep_p4[0].Py())) - (met_init.Px()*lep_p4[0].Px() + met_init.Py()*lep_p4[0].Py()))) < 40."     {options}        --save-weight --histo-name "{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}" --output {outfile} gstore_outdirs/v40/*/Ntupler_v40_{dtag}*/*/0000/*root """ % (channel, pt_cut)

#for dtag in std_dtags:
for dtag in args.input_dtags:
    proc, xsec = dtags[dtag]
    options = ('--per-weight --scale %f' % xsec) if 'MC' in dtag else ''
    proc_command = com.format(chan=chan, sys=sys, distr=distr, options=options, outfile = 'quick-test/v40/' + dtag + '.root', dtag=dtag, proc=proc)

    os.system(proc_command)




