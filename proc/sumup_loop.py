from datetime import datetime
import argparse
import logging
from os.path import isfile, basename
from sys import exit


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "run event loop recording specified standard distributions",
    epilog = """Example:\npython sumup_loop.py --debug temp/loop_tt_95.root lstore_outdirs/v37/test1/MC2016_Summer16_TTJets_powheg/MC2016_Summer16_TTJets_powheg_95.root mu_sel/std/nom/Mt_lep_met_f"""
    )

parser.add_argument('out_file',   help='output file')
parser.add_argument('inp_file',   help='input file')
parser.add_argument('--ttree',    type=str, default='ttree_out', help='path to the TTree in the file')

parser.add_argument("--wjets-nup-cut", action='store_true', help="do cut WJets for NJets stitching")

parser.add_argument("--save-weight", action='store_true', help="save the weight counters in the output file")
parser.add_argument("--per-weight",  action='store_true', help="normalize by event weight of datasets")
parser.add_argument('--try-xsec',    action='store_true', help="scale histos by xsec")

parser.add_argument('--test',        type=int, help="test run on N events")
parser.add_argument('--time',        action='store_true', help="print timing counters")
parser.add_argument('--overwrite',   action='store_true', help="overwrite output if needed")
parser.add_argument('-d', '--debug', action='store_true', help="DEBUG level of logging")

parser.add_argument('loop_definitions', nargs='+', help='standard chan/proc/sys/distr and shortcuts like: chan/all, chan/all/nom/Mt_lep_met_f, chan/incl/nom/all')

args = parser.parse_args()

if isfile(args.out_file) and not args.overwrite:
    print 'output file exists: %s' % args.out_file
    exit(0)

assert     isfile(args.inp_file)

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

logging.debug('importing ROOT')
preroot_time = datetime.now()
import ROOT
from ROOT import TFile,TDirectory
from ROOT import gDirectory
logging.debug('done')

start_time = datetime.now()
if args.time:
    print "root import time %s" % str(start_time - preroot_time)

"""
define each part separately:
chan
gen procs
sys weights (and file-based?)
distrs with ranges

then definitions:
dtag in filename + chan = the set of possible {gen_proc, sys, distrs}

define all histos before the loop, save in a nested dict with definition names
histos[chan][proc][sys][distr]

and save defs as:
defs [chan name] = (chan def_funcs[sys], [procs], [systs], [distrs])
procs  = [(name, proc def_func)]
systs  = [(name, weight_func)]
distrs = [(name, distr_funcs[sys], range?)]

--- an event can fall into multiple chan and proc defs,
    the distr param and sys must be calculated once
"""

import gen_proc_defs

# find the dtag of the files
dtag = None
for dt in gen_proc_defs.dtags_procs.keys():
    if dt in args.inp_file:
        dtag = dt
        break

assert dtag

dtag_proc_defs = gen_proc_defs.dtags_procs[dtag]

logging.debug(dtag_proc_defs)

# general process of this dtag:
general_process = dtag_proc_defs['el'][0] + '_other'
# it will be chosen if no other definitions pass


import std_defs

def proc_def_func(ids):
    return lambda event: event.gen_proc_id in ids

histos = {}
defined_channels = {}


fout = TFile(args.out_file, "RECREATE")
fout.Write()
output_dirs = []

# parse the task, and create the histos for output
# marked the procedure for performance profiling
def parse_and_create_histos():
  for definition in args.loop_definitions:
    # chan/proc/sys/distr
    dchan, dproc, dsys, ddistr = definition.split('/')
    # each is a coma-separated def

    systematics = []
    # [(sys_name, sys_weight function)]
    # std systematics per dtag
    if dsys == 'std':
        logging.debug('making dsys=std')
        for dtags, systs in std_defs.std_dtag_systs.values():
            if dtag in dtags:
                systematics = systs
                break
        if not systematics:
            raise ValueError('do not know std systematics for dtag %s in %s' % (dtag, args.inp_file))

    # manual systematics from comline
    for sys in dsys.split(','):
        logging.debug('making sys %s' % sys)

        # fill in the named blocks of systematics
        if sys in std_defs.named_systs2_weights_all:
            for k, v in std_defs.named_systs2_weights_all[sys].items():
                systematics.append((k, v))

        # fill in the individual systematics
        elif sys in std_defs.systs2_weights_all:
            systematics.append((sys, std_defs.systs2_weights_all[sys]))
        else:
            print "SKIPPING UNKNOWN SYSTEMATIC %s" % sys

    logging.debug(systematics)

    distrs = []
    # manual distrs
    for dname in ddistr.split(','):
        logging.debug('making distr %s' % dname)
        distr_def, distr_range = std_defs.distr_defs[dname]
        distrs.append((dname, distr_def, distr_range))

    for chan in dchan.split(','):
        logging.debug('making chan %s' % chan)
        # channels selection functions, including the possible object systematic selections
        #chan_def = {}
        #sel_func, sel_stage_dict_funcs = std_defs.all_channels_ev_loop[chan]
        #for sysname in sel_stage_dict_funcs:
        #    chan_def[sysname] = lambda ev: sel_func(sel_stage_dict_funcs[sysname](ev), ev)
        chan_def = std_defs.all_channels_ev_loop[chan]
        # so, I'll recalculate selection stage for all systematics
        # TODO; I need to explicitly separate object systematics, which require all these recalculations

        # process definitions in this channel for this dtag
        # gen procs are only standard for now
        # TODO: add inclusive processes
        # the processes are defined for different flavours
        if dproc == 'std':
            if 'el_' in chan[:3]:
                proc_defs = dtag_proc_defs['el']
            else:
                proc_defs = dtag_proc_defs.get('mu', dtag_proc_defs['el'])
        elif dproc == 'el':
            proc_defs = dtag_proc_defs['el']
        elif dproc == 'mu':
            proc_defs = dtag_proc_defs['mu']
        else:
            raise ValueError('process definition %s is unknown, in %s' % (dproc, definition))

        # proc defs -> into a flat list [(name, function(ev))]
        main_name, subprocs = proc_defs
        proc_def = [(main_name + '_' + subname, proc_def_func(ids)) for subname, ids in subprocs]
        proc_def.append((main_name + '_other', lambda _: False))

        # channel definitions
        defined_channels[chan] = (chan_def, proc_def, systematics, distrs)

        # trying to output quickly
        #output_channel_dir = TDirectory(chan, chan)
        fout.cd()
        output_channel_dir = fout.mkdir(chan)
        output_dirs.append(output_channel_dir)

        # and create histograms
        for proc, _ in proc_def:
          output_channel_dir.cd()
          proc_dir = output_channel_dir.mkdir(proc)
          for sys, _ in systematics:
            proc_dir.cd()
            sys_dir = proc_dir.mkdir(sys) # + '/')
            for distr, _, distr_range in distrs:
                histo_name = '_'.join([chan, proc, sys, distr])
                out_histo = std_defs.make_histo(histo_name, distr_range)
                out_histo.SetDirectory(sys_dir)
                histos[(chan, proc, sys, distr)] = out_histo

parse_and_create_histos()

# ----------------

logging.debug(args.inp_file)
tfile = TFile(args.inp_file)
ttree = tfile.Get(args.ttree)

logging.debug('TTree %s' % args.ttree)
logging.debug('with %d entries' % ttree.GetEntries())

# separated the procedure for profilimng performance
def record(histos, def_tuple, param, weight):
    histos[def_tuple].Fill(param, weight)

weight_counter = None

# marked the procedure for profiling
def loop_ttree():
  for iev, ev in enumerate(ttree):
    if args.test and iev > args.test:
        logging.debug('breaking the event loop for test')
        break

    if args.wjets_nup_cut and 'WJets_madgraph' in dtag and ev.nup > 5:
        continue

    # no functools.lru_cache in python 2.7...
    event_weights = {}
    event_distrs  = {}
    #event_gen_procs  = {}

    # test if event passes channel selections
    for chan_name, (chan_def, chan_procs, chan_systs, chan_distrs) in defined_channels.items():
        #if not chan_def(ev): continue
        # channel definition can vary with object systematics
        # --- these systematics have separate selection_stage parameters in main selection
        #     to implement the movement of events between selection categories
        # precalculate the chan def of nominal systematic:
        chan_func, chan_stage_at_sys = chan_def
        chan_nom_stage = chan_stage_at_sys['NOMINAL'](ev)
        chan_pass_nominal = chan_func(chan_nom_stage, ev)
        logging.debug("%s %s %d" % (chan_name, repr(chan_pass_nominal), chan_nom_stage))
        #if not chan_pass_nominal: continue

        # define event's gen proc ID within this channel
        event_chan_proc = general_process
        for proc_name, proc_def in chan_procs:
            if not proc_def(ev): continue
            event_chan_proc = proc_name
            break

        # calculate requested distr params and systematic weights
        # if they are not already calculated
        for sname, sweight_func in chan_systs:
            # pass channel definition at this systematic
            if sname in chan_stage_at_sys:
                sys_stage = chan_stage_at_sys[sname](ev)
                if not chan_func(sys_stage, ev): continue
            elif not chan_pass_nominal: continue
            #if not chan_def.get(sname, chan_def['NOMINAL'])(ev): continue
            logging.debug(chan_pass_nominal)

            # get systematic weight
            # hopefully getting hashed item is faster than recomputing the weight (TODO: probably it's not?)
            if sname in event_weights:
                weight = event_weights[sname]
            else:
                weight = sweight_func(ev)

            # calculate parameters for the requested distrs
            for dname, distr_func, _ in chan_distrs:
                if dname in event_distrs:
                    param = event_distrs[dname]
                else:
                    param = distr_func.get(sname, distr_func['NOMINAL'])(ev)

                # record along the way
                #histos[(chan_name, event_chan_proc, sname, dname)].Fill(param, weight)
                record(histos, (chan_name, event_chan_proc, sname, dname), param, weight)
                # 3 loops deep calculation

    if args.per_weight or args.save_weight:
        wcounter = tfile.Get('ntupler/weight_counter')
        if not wcounter:
            # try top level
            wcounter = tfile.Get('weight_counter')
        assert bool(wcounter)
        if not weight_counter:
            weight_counter = wcounter
            weight_counter.SetDirectory(0)
            weight_counter.SetName("sumup_weight_counter")
        else:
            weight_counter.Add(wcounter)

loop_ttree()

logging.debug('finished the event loop')

if args.time:
    loop_time = str(datetime.now() - start_time)
    print "loop finished, time elapsed %s" % loop_time


# scale histos as needed

# by MC event weight
if args.per_weight and dtag not in ('data', 'SingleMuon', 'SingleElectron'):
    #weight_counter = tfile.Get('ntupler/weight_counter')
    #out_histo.Scale(1./weight_counter.GetBinContent(2))
    for histo in output_histos.values():
        histo.Scale(1./weight_counter.GetBinContent(2))

# by xsec
if args.try_xsec and dtag not in ('data', 'SingleMuon', 'SingleElectron'):
    from per_dtag_script import dtags
    _, xsec = dtags.get(dtag, (None, 1.))
    #out_histo.Scale(xsec)
    for histo in histos.values():
        histo.Scale(xsec)

logging.debug('scaled all')

# write histos out

#fout = TFile(args.out_file, "RECREATE")
#fout.Write()

#fout.cd()
#top
#print 'AAAA', gDirectory.pwd()
#
#for chan_dir in output_dirs:
#    #print chan_dir.GetDirectory()
#    print type(chan_dir)
#    print chan_dir.GetName()
#    #chan_dir.Write()
#    gDirectory.WriteObject(chan_dir, chan_dir.GetName())

'''
for path_tuple, histo in histos.items():
    logging.debug(repr(path_tuple))
    histo_path = '/'.join(path_tuple[:-1]) # skip the distr name -- histo name is used there
    # check if this directory already exists in the file (a feature for future)
    #out_dir_name = ''.join(part + '/' for part in histo_path)

    #if histo_path and fout.Get(histo_path):
    #    logging.debug('found  ' + histo_path)
    #    out_dir = fout.Get(histo_path)
    #else:
    #    logging.debug('making ' + histo_path)
    #    # iteratively create each directory fout -> part0 -> part1 -> ...
    #    # somehow root did not work for creating them in 1 go
    #    out_dir = fout
    #    for directory in histo_path.split('/'):
    #        logging.debug('making ' + directory)
    #        nested_dir = out_dir.Get(directory) if out_dir.Get(directory) else out_dir.mkdir(directory)
    #        nested_dir.cd()
    #        out_dir = nested_dir

    # new directory making
    #fout.cd()
    fout.mkdir(histo_path)
    #fout.cd(histo_path)

    #histo.SetDirectory(out_dir)
    #histo.SetDirectory(gDirectory) # suggestion on root forums, does not work
    histo.SetDirectory(fout.Get(histo_path))
    histo.Write()
'''

# save weight counters etc in the top directory in the file
if args.save_weight:
    fout.cd()
    weight_counter.SetName('weight_counter')
    weight_counter.Write()

fout.Write()
fout.Close()

if args.time:
    work_time = str(datetime.now() - start_time)
    print "output written, time elapsed %s" % work_time

