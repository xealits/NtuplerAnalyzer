import argparse
import logging
from yaml import load

from os import listdir, mkdir
from os.path import isdir, isfile, join

logging.basicConfig(level=logging.INFO)


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
parser.add_argument("-q", "--queue-dir",      type=str, default='./queue_dir',
        help="the directory for storing the generated job queues (default ./queue_dir)")
parser.add_argument("-n", "--ntupler-dir",    type=str, default='/gstore/t3cms/store/user/otoldaie/',
        help="the directory with Ntupler outputs (default /gstore/t3cms/store/user/otoldaie/)")
parser.add_argument("-p", "--processing-dir", type=str, default='/lstore/cms/olek/outdirs/',
        help="the directory with Processed output (default /lstore/cms/olek/outdirs/)")

parser.add_argument("--dtags", type=str, default='std',
        help='dtags or groups of dtags to submit, separated by coma like "std,updowns" (default std)')

# example of true-storing options
#parser.add_argument("-p", "--plot",  action='store_true', help="don't save the root file, plot stacked histogram in png")


args = parser.parse_args()
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
            return info['dtag'].split('/')[1]

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
"MC2016_Summer16_W1Jets_madgraph",
"MC2016_Summer16_W2Jets_madgraph",
"MC2016_Summer16_W2Jets_madgraph_ext1",
"MC2016_Summer16_W3Jets_madgraph",
"MC2016_Summer16_W3Jets_madgraph_ext1",
"MC2016_Summer16_W4Jets_madgraph",
"MC2016_Summer16_W4Jets_madgraph_ext1",
"MC2016_Summer16_W4Jets_madgraph_ext2",
"MC2016_Summer16_WJets_amcatnlo",
"MC2016_Summer16_WJets_amcatnlo_ext2_v2",
"MC2016_Summer16_WJets_madgraph",
"MC2016_Summer16_WJets_madgraph_ext2_v1",
"MC2016_Summer16_TTJets_powheg",
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



dtags_std_data = [
"Data13TeV_SingleMuon2016B_03Feb2017_ver2",
"Data13TeV_SingleMuon2016C_03Feb2017_v1",
"Data13TeV_SingleMuon2016D_03Feb2017_v1",
"Data13TeV_SingleMuon2016E_03Feb2017_v1",
"Data13TeV_SingleMuon2016F_03Feb2017_v1",
"Data13TeV_SingleMuon2016G_03Feb2017_v1",
"Data13TeV_SingleMuon2016H_03Feb2017_ver2",
"Data13TeV_SingleMuon2016H_03Feb2017_ver3",
"Data13TeV_SingleElectron2016B_03Feb2017_ver2",
"Data13TeV_SingleElectron2016C_03Feb2017_v1",
"Data13TeV_SingleElectron2016D_03Feb2017_v1",
"Data13TeV_SingleElectron2016E_03Feb2017_v1",
"Data13TeV_SingleElectron2016F_03Feb2017_v1",
"Data13TeV_SingleElectron2016G_03Feb2017_v1",
"Data13TeV_SingleElectron2016H_03Feb2017_ver2",
"Data13TeV_SingleElectron2016H_03Feb2017_ver3",
]



dtags_std = dtags_std_mc + dtags_std_data


dtag_groups = {'std': dtags_std, 'std_data': dtags_std_data, 'std_mc': dtags_std_mc, 'updowns': dtags_updowns_old}


# parse the requested dtags and groups
requested_dtags = set()
for req in args.dtags.split(','):
    if req in dtag_groups:
        requested_dtags.update(dtag_groups[req])
    else:
        requested_dtags.update(req)

logging.info("got %d requested dtags" % len(requested_dtags))


# make the jobs
job_dir_template = args.ntupler_dir + "/{vntupler}/{dset}/Ntupler_{vntupler}_{dtag}/*/0000/"
# list this dir to get the available files
job_template = "python channel_distrs_full_loop.py -l logss " + args.processing_dir + "/{vntupler}/{vproc}/{dtag}/ -i {job_file}   || true"

# make just flat list of jobs commands
job_commands = []
for dtag in requested_dtags:
    # get first-name of the dataset corresponding to this dtag
    dset = get_dset(dtag)

    if not dset:
        logging.error("dset is not found for %s" % dtag)
        continue

    job_dirname = job_dir_template.format(vntupler=args.vNtupler, dset=dset, dtag=dtag)

    # available .root files:
    rootfiles = [f for f in listdir(job_dirname) if isfile(join(job_dirname, f)) and f[-5:] == '.root']

    for job_file in rootfiles:
        job_commands.append(job_template.format(vntupler=args.vNtupler, vproc=args.vProc, dtag=dtag, job_file=job_file))

logging.info("made %d job commands" % len(job_commands))

# split the jobs into queues according to the scheme
# at the moment standard scheme only:
standard_scheme = {1:5, 2:5, 3:5, 4:13, 5:13}
scheme = standard_scheme

n_queues = sum(standard_scheme.values())

#>>> def chunkify(lst,n):
#...     return [lst[i::n] for i in xrange(n)]

queues = [job_commands[i::n_queues] for i in xrange(n_queues)]

logging.info("made %d queues" % n_queues)


# make the job directory of the processing and write the queues
proc_queues_dir = join(args.queue_dir, args.vProc)
mkdir(proc_queues_dir)

for nod, n_queues in scheme.items():
    nod_queues_dir = join(proc_queues_dir, nod)
    mkdir(nod_queues_dir)
    nod_queues, queues = queues[:n_queues], queues[n_queues:]

    for i, nod_queue in enumerate(nod_queues):
        queue_filename = nod_queues_dir + '/q%d' % i
        with open(queue_filename, 'w') as f:
            f.write('\n'.join(nod_queue))





