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
parser.add_argument("-n", "--ntupler-dir",    type=str, default='/gstore/t3cms/store/user/otoldaie/',
        help="the directory with Ntupler outputs (default /gstore/t3cms/store/user/otoldaie/)")
parser.add_argument("-p", "--processing-dir", type=str, default='/lstore/cms/olek/outdirs/',
        help="the directory with Processed output (default /lstore/cms/olek/outdirs/)")

parser.add_argument("--dtags", type=str, default='std',
        help='dtags or groups of dtags to submit, separated by coma like "std,updowns" (default std)')

parser.add_argument("--without-dtags", type=str, default='',
        help='dtags or groups of dtags to remove from the submission list (default "")')

parser.add_argument("--scheme", type=str, help="the scheme of queue as 5,5,0,15,15 for 1,2,3,4,5 nodes")

parser.add_argument("--do-W-stitching", action='store_true', help="turn ON skipping NUP events of inclusive sample")
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


dtags_std_mc = [
"MC2016_Summer16_QCD_HT-100-200",
"MC2016_Summer16_QCD_HT-1000-1500",
"MC2016_Summer16_QCD_HT-1500-2000",
"MC2016_Summer16_QCD_HT-200-300",
"MC2016_Summer16_QCD_HT-2000-Inf",
"MC2016_Summer16_QCD_HT-300-500",
"MC2016_Summer16_QCD_HT-500-700",
"MC2016_Summer16_QCD_HT-700-1000",
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
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo_v1_ext1",
"MC2016_Summer16_DYJetsToLL_10to50_amcatnlo_v2",
"MC2016_Summer16_DYJetsToLL_50toInf_madgraph",
"MC2016_Summer16_DYJetsToLL_50toInf_madgraph_ext2_v1",
"MC2016_Summer16_SingleT_tW_5FS_powheg",
"MC2016_Summer16_SingleTbar_tW_5FS_powheg",
"MC2016_Summer16_TTJets_powheg",
"MC2016_Summer16_W1Jets_madgraph",
"MC2016_Summer16_W2Jets_madgraph",
"MC2016_Summer16_W2Jets_madgraph_ext1",
"MC2016_Summer16_W3Jets_madgraph",
"MC2016_Summer16_W3Jets_madgraph_ext1",
"MC2016_Summer16_W4Jets_madgraph",
"MC2016_Summer16_W4Jets_madgraph_ext1",
"MC2016_Summer16_W4Jets_madgraph_ext2",
"MC2016_Summer16_WJets_madgraph",
"MC2016_Summer16_WJets_madgraph_ext2_v1",
]

dtags_wjets = [
"MC2016_Summer16_W1Jets_madgraph",
"MC2016_Summer16_W2Jets_madgraph",
"MC2016_Summer16_W2Jets_madgraph_ext1",
"MC2016_Summer16_W3Jets_madgraph",
"MC2016_Summer16_W3Jets_madgraph_ext1",
"MC2016_Summer16_W4Jets_madgraph",
"MC2016_Summer16_W4Jets_madgraph_ext1",
"MC2016_Summer16_W4Jets_madgraph_ext2",
"MC2016_Summer16_WJets_madgraph",
"MC2016_Summer16_WJets_madgraph_ext2_v1",
]

dtags_wjets_amcatnlo = [
"MC2016_Summer16_WJets_amcatnlo",
"MC2016_Summer16_WJets_amcatnlo_ext2_v2",
]

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

# also
"MC2016_Summer16_TTJets_amcatnlo_backup",
"MC2016_Summer16_TTJets_powheg_herwigpp",
"MC2016_Summer16_TTJets_powheg_GluonMoveCRTune",
"MC2016_Summer16_TTJets_powheg_QCDCRTune",


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


dtags_std = dtags_std_mc + dtags_std_data


dtag_groups = {'std': dtags_std, 'std_data': dtags_std_data, 'data_el': dtags_std_data_el, 'data_mu': dtags_std_data_mu,
    'std_mc': dtags_std_mc, 'updowns': dtags_updowns,
    'singletop': dtags_singletop_mc,
    'wjets': dtags_wjets,
    'wjets_amcatnlo': dtags_wjets_amcatnlo,
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


# make the jobs
job_dir_template = args.ntupler_dir + "/{vntupler}/{dset}/Ntupler_{vntupler}_{dtag}/" #*/0000/"
# list this dir to get the available files
if args.do_W_stitching:
    job_template = "python channel_distrs_full_loop.py --do-W-stitching -l logss " + args.processing_dir + "/{vntupler}/{vproc}/{dtag}/ {chans} -i {job_file}   || true"
else:
    job_template = "python channel_distrs_full_loop.py -l logss " + args.processing_dir + "/{vntupler}/{vproc}/{dtag}/ {chans} -i {job_file}   || true"

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
    logging.debug(job_dirname)

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

queues = [job_commands[i::n_queues] for i in xrange(n_queues)]

logging.info("made %d queues" % n_queues)


# make the job directory of the processing and write the queues
proc_queues_dir = join(args.queue_dir, args.vNtupler, args.vProc)
if not isdir(proc_queues_dir):
    makedirs(proc_queues_dir)

com_file_template = """
cd
source bootup.tcsh
cmsenv
cd UserCode/NtuplerAnalyzer

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




