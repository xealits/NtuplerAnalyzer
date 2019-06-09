import argparse
import logging
from yaml import load

import glob
from os import listdir, mkdir, makedirs
from os.path import isdir, isfile, join, getmtime


'''
Input:
* ntupler run version
* processing version
* dtags for the job.

Make the list of all jobs,
split into sublists for requested queues,
save for reference in files, in per-job directories nodeN/queueM/,
so the result is a directory with N subdirectories for computing nodes and M queues on that node.
The queues are submitted over ssh in shell or make.

example of 1 job in a queue:
python channel_distrs_full_loop.py -l logss /lstore/cms/olek/outdirs/v17/1/MC2016_Summer16_DYJetsToLL_10to50_amcatnlo/ -i /gstore/t3cms/store/user/otoldaie/v17/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/Ntupler_v17_MC2016_Summer16_DYJetsToLL_10to50_amcatnlo/*/0000/MC2016_Summer16_DYJetsToLL_10to50_amcatnlo_1.root   || true

with input params:
python channel_distrs_full_loop.py -l logss /lstore/cms/olek/outdirs/vNTUPLER/vPROC/DTAG/ -i /gstore/t3cms/store/user/otoldaie/vNTUPLER/DSET-FIRST-NAME/Ntupler_vNTUPLER_DTAG/*/0000/DTAG_1.root   || true


-- TODO: the 0000 might be different, and what about first name dataset?
dtag is used for extended datasets -- the first name is the same for them, therefore on practice dtag = first name
'''

parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "generate job queues",
    epilog = "Example:\n$ python genjobs.py v17 oldrun5 --dtags std"
    )

parser.add_argument("vNtupler",    help="Ntupler output version")
parser.add_argument("vProc",       help="Processing version")
parser.add_argument("chan_def",    help="definition of selection channels")
parser.add_argument("-q", "--queue-dir",      type=str, default='./queue_dir',
        help="the directory for storing the generated job queues (default ./queue_dir)")
parser.add_argument("-n", "--ntupler-dir",    type=str, default='./gstore_outdirs',
        help="the directory with Ntupler outputs (default ./gstore_outdirs)")
parser.add_argument("-p", "--processing-dir", type=str, default='./lstore_outdirs/',
        help="the directory with Processed output (default ./lstore_outdirs/)")
parser.add_argument("-c", "--cmssw-project-dir", type=str, default='/lstore/cms/olek/CMSSW_8_0_26_patch1/src/',
        help="the directory with CMSSW/src to find crrections etc (default '/lstore/cms/olek/CMSSW_8_0_26_patch1/src/')")

parser.add_argument("--jobs-name-tag", type=str, default='job', help='change the name of the job script (default is "job", use for job management)')

parser.add_argument("--dtags", type=str, default='std',
        help='dtags or groups of dtags to submit, separated by coma like "std,updowns" (default std)')

parser.add_argument("--without-dtags", type=str, default='',
        help='dtags or groups of dtags to remove from the submission list (default "")')

parser.add_argument("--scheme", type=str, help="the scheme of queue as 5,5,0,15,15 for 1,2,3,4,5 nodes")

parser.add_argument("--submit", type=str, default='online', help="the type of the jobs (online by default, other option is 'queue')")
parser.add_argument("--mem-size",  type=str, default='1G', help="make the queue jobs with given memory size (1G default)")

parser.add_argument("--metmuegclean", type=str, default='true', help="use slimmedMETsMuEGClean MET for data")
parser.add_argument("--options", type=str, help="more options of stage2")

parser.add_argument("--do-W-stitching", action='store_true', help="turn ON skipping NUP events of inclusive sample")
parser.add_argument("--all-jets",       action='store_true', help="propagate corrections to met from all selected jets, including lep and tau matched")
parser.add_argument("--without-bSF",    action='store_true', help="don't apply b tagging SF")
parser.add_argument("--old-loop",       action='store_true', help="run on old full selection")
parser.add_argument("--acceptance-study", action='store_true', help="run the acceptance study")

parser.add_argument("-d", "--debug",    action='store_true', help="DEBUG level of logging")


args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

logging.debug("parsed args: %s" % repr(args))

assert all(isdir(d) for d in [args.queue_dir, args.ntupler_dir, args.processing_dir])


# dsets info is needed to get the DSET-FIRST-NAME of the requested dtag
# (TODO: probably I can substitute the dtag with the firstname altogether
#  -- no! I do need to distinguish ext datasets for running separately on them -- FIRST-NAME does not distinguish)
dsets_info_filename = "./dsets_info.yaml"
with open(dsets_info_filename) as f:
    dsets_info = load(f)
logging.debug("loaded dsets info from %s" % repr(dsets_info_filename))

def get_dset(dtag):
    '''get_dset(dtag)

    the info is stored as '/dset/full/name': {..., 'dtag'}

    returns first-name of the dataset corresponding to this dtag
    '''

    for dset, info in dsets_info.items():
        if info['dtag'] == dtag:
            return dset.split('/')[1]

    return None


dtags_std_mc_min_no_wjets_noqcd = [
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
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo",
"MC2016_Summer16_DYJetsToLL_50toInf_madgraph",
"MC2016_Summer16_SingleT_tW_5FS_powheg",
"MC2016_Summer16_SingleTbar_tW_5FS_powheg",
"MC2016_Summer16_TTJets_powheg",
]

dtags_mc_dibosons = [
"MC2016_Summer16_WWTo2L2Nu_powheg",
"MC2016_Summer16_WWToLNuQQ_powheg",
"MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin",
"MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin",
"MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin",
"MC2016_Summer16_WZTo3LNu_powheg",
"MC2016_Summer16_ZZTo2L2Nu_powheg",
"MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin",
]

dtags_std_mc_qcd = [
"MC2016_Summer16_QCD_HT-100-200",
"MC2016_Summer16_QCD_HT-1000-1500",
"MC2016_Summer16_QCD_HT-1500-2000",
"MC2016_Summer16_QCD_HT-200-300",
"MC2016_Summer16_QCD_HT-2000-Inf",
"MC2016_Summer16_QCD_HT-300-500",
"MC2016_Summer16_QCD_HT-500-700",
"MC2016_Summer16_QCD_HT-700-1000",
]

dtags_mc_qcd_ext = [
"MC2016_Summer16_QCD_HT-50-100",
"MC2016_Summer16_QCD_HT-200-300_ext1",
"MC2016_Summer16_QCD_HT-300-500_ext1",
"MC2016_Summer16_QCD_HT-500-700_ext1",
"MC2016_Summer16_QCD_HT-700-1000_ext1",
"MC2016_Summer16_QCD_HT-1000-1500_ext1",
"MC2016_Summer16_QCD_HT-1500-2000_ext1",
"MC2016_Summer16_QCD_HT-2000-Inf_ext1",
]

dtags_mc_qcd_mu_enriched = [
"MC2016_Summer16_QCD_MuEnriched_Pt5_1000toInf",
"MC2016_Summer16_QCD_MuEnriched_Pt5_120to170",
"MC2016_Summer16_QCD_MuEnriched_Pt5_15to20",
"MC2016_Summer16_QCD_MuEnriched_Pt5_170to300",
"MC2016_Summer16_QCD_MuEnriched_Pt5_20to30",
"MC2016_Summer16_QCD_MuEnriched_Pt15_20toInf",
"MC2016_Summer16_QCD_MuEnriched_Pt5_300to470",
"MC2016_Summer16_QCD_MuEnriched_Pt5_30to50",
"MC2016_Summer16_QCD_MuEnriched_Pt5_470to600",
"MC2016_Summer16_QCD_MuEnriched_Pt5_50to80",
"MC2016_Summer16_QCD_MuEnriched_Pt5_600to800",
"MC2016_Summer16_QCD_MuEnriched_Pt5_800to1000",
"MC2016_Summer16_QCD_MuEnriched_Pt5_80to120",
]

dtags_mc_qcd_em_enriched = [
"MC2016_Summer16_QCD_EMEnriched_Pt-20to30",
"MC2016_Summer16_QCD_EMEnriched_Pt-30to50",
"MC2016_Summer16_QCD_EMEnriched_Pt-50to80",
"MC2016_Summer16_QCD_EMEnriched_Pt-80to120",
"MC2016_Summer16_QCD_EMEnriched_Pt-120to170",
"MC2016_Summer16_QCD_EMEnriched_Pt-170to300",
"MC2016_Summer16_QCD_EMEnriched_Pt-300toInf",
"MC2016_Summer16_QCD_EMEnriched_Pt-30toInf",
"MC2016_Summer16_QCD_EMEnriched_Pt-30to40",
"MC2016_Summer16_QCD_EMEnriched_Pt-40toInf",
]

dtags_mc_qcd_em_enriched_ext = [
"MC2016_Summer16_QCD_EMEnriched_Pt-30to50_ext1",
"MC2016_Summer16_QCD_EMEnriched_Pt-50to80_ext1",
"MC2016_Summer16_QCD_EMEnriched_Pt-80to120_ext1",
"MC2016_Summer16_QCD_EMEnriched_Pt-120to170_ext1",
]




dtags_std_mc_min_no_wjets = dtags_std_mc_min_no_wjets_noqcd + dtags_std_mc_qcd

dtags_wjets_madgraph = [
"MC2016_Summer16_W1Jets_madgraph",
"MC2016_Summer16_W2Jets_madgraph",
"MC2016_Summer16_W3Jets_madgraph",
"MC2016_Summer16_W4Jets_madgraph",
"MC2016_Summer16_WJets_madgraph",
]

dtags_std_mc_min = dtags_std_mc_min_no_wjets + dtags_wjets_madgraph

dtags_std_mc_ext_wjets_madgraph = [
"MC2016_Summer16_WJets_madgraph_ext2_v1",
"MC2016_Summer16_W4Jets_madgraph_ext1",
"MC2016_Summer16_W4Jets_madgraph_ext2",
"MC2016_Summer16_W3Jets_madgraph_ext1",
"MC2016_Summer16_W2Jets_madgraph_ext1",
]

dtags_std_mc_dy = [
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo",
"MC2016_Summer16_DYJetsToLL_50toInf_madgraph",
]

dtags_std_mc_ext_dy = [
"MC2016_Summer16_DYJetsToLL_50toInf_madgraph_ext2_v1",
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo_v1_ext1",
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo_v2",
]

dtags_std_mc_ext = dtags_std_mc_ext_wjets_madgraph + dtags_std_mc_ext_dy

dtags_std_mc     = dtags_std_mc_min + dtags_std_mc_ext

dtags_wjets      = dtags_wjets_madgraph + dtags_std_mc_ext_wjets_madgraph


dtags_wjets_amcatnlo = [
"MC2016_Summer16_WJets_amcatnlo",
"MC2016_Summer16_WJets_amcatnlo_ext2_v2",
]

dtags_std_mc_amcatnlo = dtags_std_mc_min_no_wjets + dtags_std_mc_ext_dy + dtags_wjets_amcatnlo

dtags_singletop_mc = [
"MC2016_Summer16_SingleT_tW_5FS_powheg",
"MC2016_Summer16_SingleTbar_tW_5FS_powheg",
"MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo",
"MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg",
"MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg",
]

# the dtags used in Ntupler runs up to v17 have names different to systematics
# I probably won't fix it -- just will switch to using FIRST-NAME and drop dtag
dtags_updowns_old = [
"TT_TuneCUETP8M2T4_13TeV-powheg-fsrdown-pythia8",
"TT_TuneCUETP8M2T4_13TeV-powheg-fsrup-pythia8  ",
"TT_TuneCUETP8M2T4_13TeV-powheg-isrdown-pythia8",
"TT_TuneCUETP8M2T4_13TeV-powheg-isrup-pythia8  ",
"TT_TuneCUETP8M2T4down_13TeV-powheg-pythia8",
"TT_TuneCUETP8M2T4up_13TeV-powheg-pythia8",
"TT_hdampDOWN_TuneCUETP8M2T4_13TeV-powheg-pythia8",
"TT_hdampUP_TuneCUETP8M2T4_13TeV-powheg-pythia8",
]

dtags_updowns = [
"MC2016_Summer16_TTJets_powheg_fsrdown",
"MC2016_Summer16_TTJets_powheg_fsrup",
"MC2016_Summer16_TTJets_powheg_isrdown",
"MC2016_Summer16_TTJets_powheg_isrup",
"MC2016_Summer16_TTJets_powheg_CUETP8M2T4down",
"MC2016_Summer16_TTJets_powheg_CUETP8M2T4up",
"MC2016_Summer16_TTJets_powheg_hdampDOWN",
"MC2016_Summer16_TTJets_powheg_hdampUP",
]

dtags_updowns_extra = [
"MC2016_Summer16_TTJets_powheg_fsrdown_ext1",
"MC2016_Summer16_TTJets_powheg_fsrdown_ext2",
"MC2016_Summer16_TTJets_powheg_fsrup_ext1",
"MC2016_Summer16_TTJets_powheg_fsrup_ext2",
"MC2016_Summer16_TTJets_powheg_isrdown_ext1",
"MC2016_Summer16_TTJets_powheg_isrdown_ext2",
"MC2016_Summer16_TTJets_powheg_isrup_ext2",
"MC2016_Summer16_TTJets_powheg_CUETP8M2T4down_ext1",
"MC2016_Summer16_TTJets_powheg_CUETP8M2T4up_ext1",
"MC2016_Summer16_TTJets_powheg_hdampDOWN_ext1",
"MC2016_Summer16_TTJets_powheg_hdampUP_ext1",
]

dtags_updowns_all_fsr = [
"MC2016_Summer16_TTJets_powheg_fsrdown",
"MC2016_Summer16_TTJets_powheg_fsrup",
"MC2016_Summer16_TTJets_powheg_fsrdown_ext1",
"MC2016_Summer16_TTJets_powheg_fsrdown_ext2",
"MC2016_Summer16_TTJets_powheg_fsrup_ext1",
"MC2016_Summer16_TTJets_powheg_fsrup_ext2",
]

# also
"MC2016_Summer16_TTJets_amcatnlo_backup",
"MC2016_Summer16_TTJets_powheg_herwigpp",
"MC2016_Summer16_TTJets_powheg_GluonMoveCRTune",
"MC2016_Summer16_TTJets_powheg_QCDCRTune",
"MC2016_Summer16_TTJets_amcatnlo"
"MC2016_Summer16_TTJets_powheg_QCDCRTune_ext1"
"MC2016_Summer16_TTJets_powheg_herwigpp_ext2"
"MC2016_Summer16_TTJets_powheg_herwigpp_ext3"


dtags_std_data_mu = [
"Data13TeV_SingleMuon2016B_03Feb2017_ver2",
"Data13TeV_SingleMuon2016C_03Feb2017_v1",
"Data13TeV_SingleMuon2016D_03Feb2017_v1",
"Data13TeV_SingleMuon2016E_03Feb2017_v1",
"Data13TeV_SingleMuon2016F_03Feb2017_v1",
"Data13TeV_SingleMuon2016G_03Feb2017_v1",
"Data13TeV_SingleMuon2016H_03Feb2017_ver2",
"Data13TeV_SingleMuon2016H_03Feb2017_ver3",
]

dtags_data_mu_aug = [
'Data13TeV_SingleMuon2016B_07Aug17_ver2',
'Data13TeV_SingleMuon2016C_07Aug17_v1',
'Data13TeV_SingleMuon2016D_07Aug17_v1',
'Data13TeV_SingleMuon2016E_07Aug17_v1',
'Data13TeV_SingleMuon2016F_07Aug17_v1',
'Data13TeV_SingleMuon2016G_07Aug17_v1',
'Data13TeV_SingleMuon2016H_07Aug17_v1',
]

dtags_std_data_el = [
"Data13TeV_SingleElectron2016B_03Feb2017_ver2",
"Data13TeV_SingleElectron2016C_03Feb2017_v1",
"Data13TeV_SingleElectron2016D_03Feb2017_v1",
"Data13TeV_SingleElectron2016E_03Feb2017_v1",
"Data13TeV_SingleElectron2016F_03Feb2017_v1",
"Data13TeV_SingleElectron2016G_03Feb2017_v1",
"Data13TeV_SingleElectron2016H_03Feb2017_ver2",
"Data13TeV_SingleElectron2016H_03Feb2017_ver3",
]

dtags_std_data = dtags_std_data_mu + dtags_std_data_el

dtags_std          = dtags_std_mc          + dtags_std_data
dtags_std_amcatnlo = dtags_std_mc_amcatnlo + dtags_std_data
dtags_std_min      = dtags_std_mc_min      + dtags_std_data


dtag_groups = {'std': dtags_std, 'std_amcatnlo': dtags_std_amcatnlo,
    'std_min': dtags_std_min, 'std_data': dtags_std_data, 'data_el': dtags_std_data_el, 'data_mu': dtags_std_data_mu,
    'std_mc': dtags_std_mc,
    'dy_std': dtags_std_mc_dy,
    'dy_ext': dtags_std_mc_ext_dy,
    'std_qcd': dtags_std_mc_qcd,
    'qcd_ext': dtags_mc_qcd_ext,
    'qcd_muenriched': dtags_mc_qcd_mu_enriched,
    'qcd_emenriched': dtags_mc_qcd_em_enriched,
    'qcd_emenriched_ext': dtags_mc_qcd_em_enriched_ext,
    'updowns': dtags_updowns, 'updowns_extra': dtags_updowns_extra,
    'updowns_allfsr': dtags_updowns_all_fsr,
    'singletop': dtags_singletop_mc,
    'dibosons':  dtags_mc_dibosons,
    'wjets': dtags_wjets,
    'wjets_amcatnlo': dtags_wjets_amcatnlo,
    'data_mu_aug': dtags_data_mu_aug,
}


# parse the requested dtags and groups
requested_dtags = set()

logging.debug("got %d %s dtags" % (len(args.dtags.split(',')), args.dtags))
#for req in (args.dtags.split(',') if ',' in args.dtags else [args.dtags,]):
for req in args.dtags.split(','):
    logging.debug("%s" % req)
    if req in dtag_groups:
        requested_dtags.update(dtag_groups[req])
    else:
        requested_dtags.update([req])

logging.debug("got %d %s subtract dtags" % (len(args.without_dtags.split(',')), args.without_dtags))
for req in args.without_dtags.split(','):
    logging.debug("%s" % req)
    if req in dtag_groups:
        requested_dtags.difference_update(dtag_groups[req])
    else:
        requested_dtags.difference_update([req])

logging.info("got %d requested dtags" % len(requested_dtags))
for dtag in requested_dtags:
    logging.debug(dtag)


# make the jobs
job_dir_template = args.ntupler_dir + "/{vntupler}/{dset}/Ntupler_{vntupler}_{dtag}/" #*/0000/"
# list this dir to get the available files
add_options = ""
if args.do_W_stitching:
    add_options += "--do-W-stitching "
if args.all_jets:
    add_options += "--all-jets "
if args.without_bSF:
    add_options += "--without-bSF "
if args.old_loop:
    add_options += "--old-loop "
if args.metmuegclean:
    add_options += "--metmuegclean %s " % args.metmuegclean
if args.options:
    add_options += "--options %s " % args.options

if args.acceptance_study:
    job_template = "python signal_acceptance.py " + args.processing_dir + "/{vntupler}/{vproc}/{dtag}/ {job_file}   || true"
    # time python signal_acceptance.py out_signal_accept_v34.root gstore_outdirs/v34/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v34_MC2016_Summer16_TTJets_powheg/180717_001103/0000/MC2016_Summer16_TTJets_powheg_90.root

elif args.submit == 'online':
    job_template = "python channel_distrs_full_loop.py " + add_options + " -l logss " + args.processing_dir + "/{vntupler}/{vproc}/{dtag}/ {chans} -i {job_file}   || true"
elif 'queue' in args.submit:
    job_template = "python channel_distrs_full_loop.py " + add_options + " -l logss " + args.processing_dir + "/{vntupler}/{vproc}/{dtag}/ {chans} -i {job_file}"
else:
    raise ValueError('unknown type of jobs submition "%s"' % args.submit)


# make just flat list of jobs commands
job_commands = []
for dtag in requested_dtags:
    logging.debug(dtag)
    # get first-name of the dataset corresponding to this dtag
    dset = get_dset(dtag)

    if not dset:
        logging.error("dset is not found for %s" % dtag)
        continue

    job_dirname = job_dir_template.format(vntupler=args.vNtupler, dset=dset, dtag=dtag)

    if isdir(job_dirname):
        logging.debug(job_dirname)
    else:
        logging.error('no directory ' + job_dirname)
        continue

    # there might be several submisions
    # pick up the last one

    job_submisions = [join(job_dirname, d) for d in listdir(job_dirname)]

    if len(job_submisions) != 1:
        logging.warning('not 1 dir: %d picking up the last one %s' % (len(job_submisions), dtag))
    job_dirname = job_submisions[0]
    logging.debug(job_dirname)

    job_dirname = max(job_submisions, key=getmtime)

    # need to check if I should use all subdirs here
    # then I'll just glob everything, including .root
    if len(listdir(job_dirname)) > 1:
        logging.error("many out dirs: %s" % job_dirname)

    job_dirname = job_dirname + '/0000/'

    # available .root files:
    rootfiles = [f for f in glob.glob(job_dirname + '/*.root') if isfile(f)]
    logging.debug("%d .root files" % len(rootfiles))

    for job_file in rootfiles:
        job_commands.append(job_template.format(vntupler=args.vNtupler, vproc=args.vProc, dtag=dtag, job_file=job_file, chans=args.chan_def))

logging.info("made %d job commands" % len(job_commands))

# split the jobs into queues according to the scheme
# at the moment standard scheme only:
standard_scheme = {1:5, 2:5, 3:5, 4:13, 5:13}
scheme_one3 = {1:6, 2:6, 3:1, 4:14, 5:14}
scheme_two3 = {1:6, 2:6, 3:2, 4:14, 5:14}
scheme_no3 =  {1:6, 2:6, 3:0, 4:14, 5:14}
scheme_no3_loaded =  {1:7, 2:7, 3:0, 4:15, 5:15}
scheme_no3no5 =  {1:7, 2:7, 3:0, 4:15, 5:0}
scheme = scheme_no3_loaded # standard_scheme

if args.scheme:
    custom_scheme = [int(i) for i in args.scheme.split(',')]
    scheme = {i:custom_scheme[i-1] for i in range(1,6)}

logging.info("scheme: %s" % str(scheme))

n_queues = sum(scheme.values())

#>>> def chunkify(lst,n):
#...     return [lst[i::n] for i in xrange(n)]

# just keep all jobs for queue submition
# split into queues for online submition
queues = [job_commands] if args.submit == 'queue' else [job_commands[i::n_queues] for i in xrange(n_queues)]

logging.info("made %d queues for %s jobs" % (n_queues, args.submit))


# make the job directory of the processing and write the queues
proc_queues_dir = join(args.queue_dir, args.vNtupler, args.vProc)
if not isdir(proc_queues_dir):
    makedirs(proc_queues_dir)

# if the jobs is queue make just 1 directory "queue" with all the jobs to submit
if 'queue' in args.submit:
    # TODO append the shell template to the job
    # operation: in `proc/` run `source queue_dir/v40/u4test/jobs_dir/job_1`
    if 'queue_online' == args.submit:
        job_template = """{{job}}\n"""
    else:
        job_template = """#!/bin/sh
pwd
export X509_USER_PROXY={X509_USER_PROXY}
export SCRAM_ARCH={SCRAM_ARCH}
export BUILD_ARCH={SCRAM_ARCH}
export VO_CMS_SW_DIR={VO_CMS_SW_DIR}
source $VO_CMS_SW_DIR/cmsset_default.sh
export CMS_PATH=$VO_CMS_SW_DIR
cd {project_dir}
cmsenv
cd UserCode/NtuplerAnalyzer/proc/
{{job}}
"""
    '''eval `scramv1 runtime -sh`
    cd -
    ulimit -c 0;'''

    from os import environ

    vars_for_the_job = dict(environ)
    # after boot.tcsh everything is in the vars
    project_dir = args.cmssw_project_dir
    ntupler_proc_dir = 'UserCode/NtuplerAnalyzer/proc/'
    vars_for_the_job.update(project_dir=project_dir)
    job_template = job_template.format(**vars_for_the_job)

    # make the jobs dir
    jobs_dir = join(project_dir, ntupler_proc_dir, proc_queues_dir, 'jobs_dir')
    if not isdir(jobs_dir):
        mkdir(jobs_dir)

    # write all the job files
    job_filenames = []
    for i, a_job in enumerate(job_commands):
        job_name = '/%s_%05d' % (args.jobs_name_tag, i)
        job_filename = jobs_dir + job_name
        job_filenames.append(job_filename)
        with open(job_filename, 'w') as f:
            f.write(job_template.format(job=a_job))

    #make the queue submition file
    submition_file = jobs_dir + '/submit_%s' % args.jobs_name_tag
    with open(submition_file, 'w') as f:
        f.write('\n'.join("""qsub -l h_vmem={mem_size} '{jobsh}' """.format(mem_size=args.mem_size, jobsh=j_fname) for j_fname in job_filenames) + '\n')

    print 'submit list in:'
    print submition_file

elif args.submit == 'online':
    com_file_template = """
cd
source bootup.tcsh
cmsenv
cd UserCode/NtuplerAnalyzer/proc

{queue_commands}

bash
"""

    queue_command_template = "bash {queue_dir}/{queue_name}  &"

    # write the queue files for each node
    for nod, n_queues in scheme.items():
        nod_queues_dir = join(proc_queues_dir, str(nod))
        if not isdir(nod_queues_dir):
            mkdir(nod_queues_dir)
        nod_queues, queues = queues[:n_queues], queues[n_queues:]

        queue_commands = []
        for i, nod_queue in enumerate(nod_queues):
            queue_name = '/q%d' % i
            queue_filename = nod_queues_dir + queue_name
            queue_commands.append(queue_command_template.format(queue_dir=nod_queues_dir, queue_name=queue_name))
            with open(queue_filename, 'w') as f:
                f.write('\n'.join(nod_queue) + '\n')

        # write the command file for the node
        with open(nod_queues_dir + '/com', 'w') as f:
            f.write(com_file_template.format(queue_commands = '\n'.join(queue_commands)))

else:
    raise ValueError('unknown type of jobs submition "%s"' % args.submit)




