import logging
import os
from os import listdir, makedirs
from os.path import isfile, isdir, join, exists
import argparse


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "generate loop selection jobs",
    epilog = """Example:\npython loop_selections.py"""
    )

parser.add_argument('-d', '--debug',  action='store_true', help="DEBUG level of logging")
parser.add_argument('--no-queue',     action='store_true', help="do not prepare a queue job")
parser.add_argument('--only-nominal', action='store_true', help="only nominal syst")

parser.add_argument('--job-dir', type=str, default='batch_jobs/', help='set a custom directory for job files')
parser.add_argument('--chan-groups', type=str, help='set the channel groups for jobs')

parser.add_argument('--chans', type=str, help='set the channels to select')
parser.add_argument('--sys',   type=str, help='set the systematics to run on')
parser.add_argument('--distributions', type=str, help='set the distributions to produce')

parser.add_argument('--nt',     type=str, default='v37',    help='NT of the run')
parser.add_argument('--proc',   type=str, default='test1',  help='proc')
parser.add_argument('--distrs', type=str, default='bunch2', help='distrs bunch')

parser.add_argument('--options', type=str, default='', help='options of the sumup_loop.py script')
parser.add_argument('--dtags',         type=str, help='set which dtags to process')
parser.add_argument('--without-dtags', type=str, help='exclude these dtags from processing')

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

nt         = args.nt
proc       = args.proc
distrs_run = args.distrs

template_out_dir = 'jobsums/distrs/{nt}_{proc}_{distrs_run}/all/{dtag}/'

#job_dir = 'batch_jobs/'

if args.no_queue:
    template_queue_job = "%s"
else:
    template_queue_job = """
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

%s
"""


template_command = "time python sumup_loop.py {options} {out_file} {inp_file} {definitions}"

from std_defs import sample_info, channels_distrs, extend_full_sys_list

# construct {dtag: sys}
standard_dtags = {}
standard_all_dtags = {}
for nickname, (dtags, systs) in sample_info.items():
    s = systs
    standard_all_dtags.update({dtag: s for dtag in dtags})
    if 'tt_syst' in nickname: continue
    standard_dtags.update({dtag: s for dtag in dtags})


# get the requested dtags or use the standard ones
requested_dtags = {}
if args.dtags:
    user_dtags = args.dtags.split(',')
    requested_dtags.update({dtag: standard_all_dtags[dtag] for dtag in user_dtags})
else:
    requested_dtags = standard_dtags

if args.without_dtags:
    for dtag in args.without_dtags.split(','):
        requested_dtags.pop(dtag)

logging.debug(requested_dtags)

# construct [(channels, distrs)] list
#requested_chan_defs = []

all_std_chan_groups = ['tt_dileptons',
                     'fit_tt_leptau', 'fit_tt_leptau_lj', 'fit_tt_leptau_ljout',
                     'fit_tt_leptau_Vloose', 'fit_tt_leptau_Vloose_lj', 'fit_tt_leptau_Vloose_ljout',
                     'fit_tt_leptau_Tight', 'fit_tt_leptau_Tight_lj', 'fit_tt_leptau_Tight_ljout',
                     'tt_leptau', 'tt_leptau_lj', 'tt_leptau_ljout',
                     'tt_leptau_Vloose', 'tt_leptau_Vloose_lj', 'tt_leptau_Vloose_ljout',
                     'tt_leptau_Tight', 'tt_leptau_Tight_lj', 'tt_leptau_Tight_ljout',
                     'dy_dileptons', 'wjets']

now_chan_groups = ['tt_dileptons', 'wjets', 'dy_dileptons', 'fit_tt_leptau', 'fit_tt_leptau_lj', 'fit_tt_leptau_ljout', 'fit_tt_leptau_Vloose', 'fit_tt_leptau_Vloose_lj', 'fit_tt_leptau_Vloose_ljout']

if args.chan_groups:
    now_chan_groups = args.chan_groups.split(',')

# find files of each dtag, construct the job command
for dtag, systs in requested_dtags.items():
  logging.debug('%s %s' % (dtag, repr(systs)))

  # unpack the nicknamed systematics
  systs = extend_full_sys_list(systs)
  if args.sys:
      systs = args.sys.split(',')

  # create output directory for the whole dtag
  out_dtag = template_out_dir.format(nt=nt, proc=proc, distrs_run=distrs_run, dtag=dtag) # + '/' + chan_group

  # prepare the loop plotting jobs on the stage2 files
  # the directory with the input
  mypath = 'lstore_outdirs/%s/%s/%s' % (nt, proc, dtag)

  if not isdir(mypath):
      continue

  # get the filenames of input
  onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and '.root' in f]

  if not onlyfiles:
      continue

  all_definitions = []

  for chan_group in now_chan_groups:

      chans, distrs, allowed_systs = channels_distrs[chan_group]
      allowed_systs = extend_full_sys_list(allowed_systs)
      if args.only_nominal:
          allowed_systs = ['NOMINAL']
      elif args.sys:
          allowed_systs = args.sys.split(',')

      if args.distributions:
          distrs = args.distributions.split(',')

      if args.chans:
          chans = args.chans.split(',')

      # prepare pre-definitions, only systematics are missing
      #template_def = '{ch}/std/{{systs}}/{d}'.format(ch=','.join(chans), d=','.join(distrs))
      #requested_chan_defs.append(template_def)
      requested_chan_defs = ['{ch}/std/{{systs}}/{d}'.format(ch=ch, d=','.join(distrs)) for ch in chans]

      pre_definitions = ' '.join(requested_chan_defs)

      # make the definitions for the dtag with the systs
      #definitions = pre_definitions.format(systs=systs)
      definitions = pre_definitions.format(systs=','.join([s for s in systs if s in allowed_systs]))
      logging.debug('%s' % definitions)
      all_definitions.append((chan_group, definitions))

      # create output directory for this group of channels
      out_dir = out_dtag + '/' + chan_group
      if not os.path.exists(out_dir):
          os.makedirs(out_dir)

  # make a job for each input filename
  for fname in onlyfiles:
      inp_file = mypath  + '/' + fname
      job_script = args.job_dir + '/' + '_'.join(('j', fname)) #, chan_group))
      with open(job_script, 'w') as f:
          #
          commands = []
          for chan_group, defs in all_definitions:
              out_file = out_dtag + '/' + chan_group + '/' + fname
              commands.append(template_command.format(inp_file=inp_file, out_file=out_file, definitions=defs, options=args.options))
              #command = '\n'.join(template_command.format(inp_file=inp_file, out_file=out_file, definitions=a_def) for a_def in defs)

          all_commands = '\n'.join(commands)
          f.write(template_queue_job % all_commands + '\n')

