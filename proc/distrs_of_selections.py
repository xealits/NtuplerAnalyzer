import logging
import argparse


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "generate jobs for distrs of selections",
    epilog = "Example:\n$ python distrs_of_selections.py lstore_outdirs/merge-sets/v36/test1/ jobsums/distrs/v36test1_distrs1/all/"
    )

parser.add_argument('inp_dir',    help="input directory, lstore or merge")
parser.add_argument('out_dir',    help="output directory, like jobsums/distrs/v36test1_distrs1/all/")
parser.add_argument('--set-sys',     type=str, help="set systematics for all samples")
parser.add_argument('--set-samples', type=str, help="set samples to process")
parser.add_argument('--select-channels', default='all_leptau_joined', type=str, help="channels to select (with hard-coded distrs if not overwritten)")
parser.add_argument('--set-distrs',  type=str, help='overwrite distrs "draw_command,distr_name[-dr2,nm2]"')
parser.add_argument('--jobs-dir',     type=str, default='batch_jobs', help="the directory for jobs scripts")

parser.add_argument("-d", "--debug",    action='store_true', help="DEBUG level of logging")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)


distr_ranges = {'Mt_lep_met_c': '--custom-range 0,20,40,60,80,100,130,160,200,250',
    'Mt_lep_met_c2': '--custom-range 0,20,40,60,80,100,120,140,170,200,250,500',
    'Mt_lep_met_f':  '--histo-range 20,0,250',
    'met_f':         '--histo-range 20,0,300',
    'met_c':         '--custom-range 0,20,40,60,80,100,120,140,200,500',
    'dilep_mass':    '--histo-range 100,0,400',
    'lep_pt':        '--histo-range 20,0,150',
    'lep_eta':       '--histo-range 26,-2.6,2.6',
    'tau_sv_sign':   '--histo-range 42,-1,20',
    'tau_pt':        '--histo-range 20,0,100',
    'tau_eta':       '--histo-range 26,-2.6,2.6',
    'yield':         '--histo-range 3,0.0,3.0'
}

distrs_leptonic = [('std_mt_vars', 'Mt_lep_met_c'), ('std_mt_vars', 'Mt_lep_met_c2'), ('event_met_lep_mt', 'Mt_lep_met_f'), ('event_dilep_mass', 'dilep_mass'), ('event_leptons[0].pt()', 'lep_pt')]
distrs_tauonic_std  = [('std_mt_vars', 'Mt_lep_met_c'), ('std_mt_vars', 'Mt_lep_met_c2'), ('event_met_lep_mt', 'Mt_lep_met_f'), ('event_dilep_mass', 'dilep_mass'), ('event_leptons[0].pt()', 'lep_pt'),
        ('event_taus[0].pt()', 'tau_pt'), ('event_taus[0].eta()', 'tau_eta'),
        ('event_taus_sv_sign[0]', 'tau_sv_sign')]

distrs_tauonic_2  = [
        ('event_leptons[0].eta()', 'lep_eta'),
        ('std_met_vars', 'met_c'), ('std_met_vars', 'met_f')]

distrs_tauonic_fit  = [('std_mt_vars', 'Mt_lep_met_c'), ('std_mt_vars', 'Mt_lep_met_c2'), ('event_met_lep_mt', 'Mt_lep_met_f')]

distrs_tauonic  = distrs_tauonic_std + distrs_tauonic_2
distrs_tauonic  = distrs_tauonic_std

distrs_for_leptau = distrs_tauonic_fit

# channels and conditions
channels_defs = {
'tt_dileptons'    : [(['tt_elmu'], '', distrs_leptonic)],
'tt_leptauSV'     : [(['el_selSV', 'el_selSVVloose', 'el_selSV_ss', 'el_selSVVloose_ss', 'mu_selSV', 'mu_selSVVloose', 'mu_selSV_ss', 'mu_selSVVloose_ss'], '', distrs_tauonic)],

'tt_leptau'       : [(['el_sel',       'el_sel_ss',       'mu_sel',       'mu_sel_ss'], '', distrs_for_leptau)],
'tt_leptau_lj'    : [(['el_sel_lj',    'el_sel_lj_ss',    'mu_sel_lj',    'mu_sel_lj_ss'], '', distrs_for_leptau)],
'tt_leptau_ljout' : [(['el_sel_ljout', 'el_sel_ljout_ss', 'mu_sel_ljout', 'mu_sel_ljout_ss'], '', distrs_for_leptau)],

'tt_leptau_Vloose'       : [(['el_selVloose',       'el_selVloose_ss',       'mu_selVloose',       'mu_selVloose_ss'], '', distrs_for_leptau)],
'tt_leptau_Vloose_lj'    : [(['el_selVloose_lj',    'el_selVloose_lj_ss',    'mu_selVloose_lj',    'mu_selVloose_lj_ss'], '', distrs_for_leptau)],
'tt_leptau_Vloose_ljout' : [(['el_selVloose_ljout', 'el_selVloose_ljout_ss', 'mu_selVloose_ljout', 'mu_selVloose_ljout_ss'], '', distrs_for_leptau)],

'tt_leptau_Tight'       : [(['el_selTight',       'el_selTight_ss',       'mu_selTight',       'mu_selTight_ss'], '', distrs_for_leptau)],
'tt_leptau_Tight_lj'    : [(['el_selTight_lj',    'el_selTight_lj_ss',    'mu_selTight_lj',    'mu_selTight_lj_ss'], '', distrs_for_leptau)],
'tt_leptau_Tight_ljout' : [(['el_selTight_ljout', 'el_selTight_ljout_ss', 'mu_selTight_ljout', 'mu_selTight_ljout_ss'], '', distrs_for_leptau)],

'dy_dileptons'    : [(['dy_mumu',  'dy_elel'],  '--cond-com "std_mt_vars < 40."', distrs_leptonic)],
'dy_leptau'       : [(['dy_mutau', 'dy_eltau'], '--cond-com "event_taus_sv_sign[0] > 2.5 && std_mt_vars < 40."', distrs_tauonic)],
}

channels_defs['sv_test']    = channels_defs['tt_dileptons'] + channels_defs['tt_leptau'] + channels_defs['tt_leptauSV'] + channels_defs['dy_dileptons'] + channels_defs['dy_leptau']
channels_defs['std_leptau'] = channels_defs['tt_leptau']    + channels_defs['tt_leptau_lj'] + channels_defs['tt_leptau_ljout']
channels_defs['all_leptau'] = channels_defs['tt_leptau'] + channels_defs['tt_leptau_lj'] + channels_defs['tt_leptau_ljout'] + \
    channels_defs['tt_leptau_Tight']  + channels_defs['tt_leptau_Tight_lj']  + channels_defs['tt_leptau_Tight_ljout'] + \
    channels_defs['tt_leptau_Vloose'] + channels_defs['tt_leptau_Vloose_lj'] + channels_defs['tt_leptau_Vloose_ljout']
select_channels_all_leptau = channels_defs['all_leptau']
# simple group
channels_defs['all_leptau_joined'] = [([','.join(chans)], cond, distrs) for chans, cond, distrs  in select_channels_all_leptau]
select_channels_all_leptau_joined  = channels_defs['all_leptau_joined']

select_sparse_channels = channels_defs ['all_leptau']
select_joined_channels = channels_defs ['all_leptau_joined']

select_channels = channels_defs[args.select_channels]

# dtags and systematics
sample_info = {
'data': (['SingleMuon', 'SingleElectron'], ['nom']),
'tt'  : (['MC2016_Summer16_TTJets_powheg'],  ["nom,common", "obj", "tt_weights", "tt_hard", "tt_pdf1", "tt_pdf10", "tt_pdf20", "tt_pdf30", "tt_pdf40", "tt_pdf50,tt_alpha"]), #select_sparse_channels
'other_mc': (['MC2016_Summer16_DYJetsToLL_10to50_amcatnlo',
'MC2016_Summer16_DYJetsToLL_50toInf_madgraph',
'MC2016_Summer16_SingleT_tW_5FS_powheg',
'MC2016_Summer16_SingleTbar_tW_5FS_powheg',
'MC2016_Summer16_W1Jets_madgraph',
'MC2016_Summer16_W2Jets_madgraph',
'MC2016_Summer16_W3Jets_madgraph',
'MC2016_Summer16_W4Jets_madgraph',
'MC2016_Summer16_WJets_madgraph',
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
'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg'], ['nom', 'common', 'obj']),

'qcd_mc': (['MC2016_Summer16_QCD_HT-100-200',
'MC2016_Summer16_QCD_HT-200-300',
'MC2016_Summer16_QCD_HT-300-500',
'MC2016_Summer16_QCD_HT-500-700',
'MC2016_Summer16_QCD_HT-700-1000',
'MC2016_Summer16_QCD_HT-1000-1500',
'MC2016_Summer16_QCD_HT-1500-2000',
'MC2016_Summer16_QCD_HT-2000-Inf'], ['nom'])
}


intro = """
#!/bin/sh
export X509_USER_PROXY=/home/t3cms/olek/x509_proxy
export SCRAM_ARCH=slc6_amd64_gcc530
export BUILD_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
export CMS_PATH=$VO_CMS_SW_DIR
cd /lstore/cms/olek/CMSSW_8_0_26_patch1/src/
cmsenv
cd UserCode/NtuplerAnalyzer/proc/

{draw_commands}
"""

draw_command_template = """python -W ignore sumup_ttree_draw.py --cut-w0jets "{draw}"           --ttree ttree_out {histo_range} --output {output_dir}/{dtag}_{chan}_{sys}_{distr_name}.root  --std-histos --histo-name {chan}/std_procs/{sys}/{distr_name} {cond} --per-weight --try-xsec --save-weight {merge_dir}/*{dtag}*.root"""

#batch_jobs/j_${dtag}_${chans}_${systs}
#merge_dir ='v27/dilep1'
merge_dir = 'lstore_outdirs/merge-sets/v25/resub1'
merge_dir = 'lstore_outdirs/merge-sets/v36/test1'
merge_dir = args.inp_dir

output_dir = 'temp/'
output_dir = 'quick-test/v27-dilep2'
output_dir = 'quick-test/v25v26-resub2_data_resub2'
output_dir = 'quick-test/v25v26-resub2'
output_dir = 'quick-test/v25v26-resub3_new_data'
output_dir = 'quick-test/v36-test1'
output_dir = args.out_dir


samples = [(['MC2016_Summer16_W2Jets_madgraph'], ['nom', 'common', 'obj'])]
samples = [(['MC2016_Summer16_W3Jets_madgraph'], ['nom'])]
samples = [(['MC2016_Summer16_TTJets_powheg'],   ['nom', 'common', 'obj', "tt_weights", "tt_hard", "tt_pdf1", "tt_pdf10", "tt_pdf20", "tt_pdf30", "tt_pdf40", "tt_pdf50,tt_alpha"])]
samples = [sample_info['qcd_mc']]
samples = [sample_info['data']]
samples = [sample_info['data'], sample_info['tt'], sample_info['other_mc'], sample_info['qcd_mc']]

if args.set_samples:
    samples_list = args.set_samples.split(',')
    samples = [sample_info.get(samp) in samples_list] # TODO: implement separate dtags

# set systematics
if args.set_sys:
    syst_list = args.set_sys.split(',')
    samples = [(dtags, syst_list) for dtags, _ in samples]

logging.info("merge_dir = %s  -->  output_dir = %s" % (merge_dir, output_dir))

for dtags, systs in samples:
    for dtag in dtags:
        for sys in systs:
            for chans, cond, distrs  in select_channels:
                if args.set_distrs:
                    distrs = [distr.split(',') for distr in args.set_distrs.split('-')]
                for chan in chans:
                    with open(args.jobs_dir + '/j_{dtag}_{chan}_{sys}'.format(dtag=dtag, chan=chan, sys=sys), 'w') as f:
                        draw_coms = '\n'.join([draw_command_template.format(draw=draw, distr_name=name, dtag=dtag, chan=chan, output_dir=output_dir, merge_dir=merge_dir, sys=sys, cond=cond, histo_range=distr_ranges[name]) for draw, name in distrs])
                        f.write(intro.format(draw_commands=draw_coms) + '\n')


