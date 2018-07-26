import argparse
import logging
from os.path import isfile, basename
import os
#from subprocess import call




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "plot selections for all known MC dtags: channel/proc/systematic",
    epilog = """Example:\npython selections.py "met_init.pt()" --histo-range 200,0,200 --histo-name data_met_init --output data_met_init.root gstore_outdirs/v28/SingleMuon/Ntupler_v28_Data13TeV_SingleMuon2016*/1*/*/*root
python sumup_ttree_distrs.py "event_leptons[0].pt()" --ttree ttree_out --cond-com "selection_stage == 5" --histo-range 50,0,200 --output test1_lep_pt.root --histo-name foo/bar/test1_lep_pt --save-weight MC2016_Summer16_TTJets_powheg_test1.root """
    )

parser.add_argument('draw_com', type=str, help='Draw("??") definition of the distribution')
parser.add_argument('--cond-com', type=str, default="", help='Draw("", "??") definition of the selected channel')
parser.add_argument('--weight',   type=str, default="", help='Draw("", "?? * ()") definition of the weight for the selection')
parser.add_argument('--ttree',    type=str, default="ntupler/reduced_ttree", help='path to the TTree in the file')
parser.add_argument('--distr',  type=str, default="lep_pt", help='the name of the distribution')
parser.add_argument('--sys',    type=str, default="NOMINAL", help='the systematic')
parser.add_argument('--chan',   type=str, default="ctr_old_mu_sel", help='the name of the channel')
parser.add_argument('--histo-range',  type=str, default=None, help='optionally set the range')
parser.add_argument('--custom-range', type=str, default=None, help='optionally set the custom range')
#parser.add_argument('--histo-color', type=str, default=None, help='optional rgb color, like `255,255,255`')

parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")
#parser.add_argument("--output", type=str, default="output.root", help="filename for output")

#parser.add_argument("--save-weight", action='store_true', help="save the weight counters in the output file")
#parser.add_argument("--per-weight",  action='store_true', help="normalize by event weight of datasets")
#parser.add_argument("--scan",        action='store_true', help="also scan events ant print out")
#parser.add_argument("--get-maximum", action='store_true', help="just find maxima on the draw_com")

parser.add_argument('input_files', nargs='+', help="""the files to sum up, passed by shell, as: merge-sets/v25/pStage2Test2/MC*root""")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

from gen_proc_defs import *

genprocs_mutau_dy    = ('dy',    [('tautau', [genproc_dy_tautau])])
genprocs_mutau_qcd   = ('qcd',   [])
genprocs_mutau_s_top = ('s_top', [('mutau', [genproc_stop_mu]),  ('lj', [genproc_stop_lj])])
genprocs_mutau_tt    = ('tt',    [('mutau', [genproc_tt_mutau, genproc_tt_mutau3ch]),
                                  ('lj', [genproc_tt_lj, genproc_tt_ljb, genproc_tt_ljw, genproc_tt_ljo]),
                                  ('taultauh', [genproc_tt_taultauh]),
                                  ('taulj', [genproc_tt_taulj])])
genprocs_mutau_wjets  = ('wjets', [('taul', [genproc_wjets_taul]), ('tauh', [genproc_wjets_tauh])])
genprocs_mutau_dibosons = ('dibosons', [])

dtags = {
'MC2016_Summer16_DYJetsToLL_10to50_amcatnlo'   : genprocs_mutau_dy,
'MC2016_Summer16_DYJetsToLL_50toInf_madgraph'  : genprocs_mutau_dy,
'MC2016_Summer16_QCD_HT-100-200'    : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_HT-1000-1500'  : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_HT-1500-2000'  : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_HT-200-300'    : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_HT-2000-Inf'   : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_HT-300-500'    : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_HT-500-700'    : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_HT-700-1000'   : genprocs_mutau_qcd,
'MC2016_Summer16_SingleT_tW_5FS_powheg'                     : genprocs_mutau_s_top,
'MC2016_Summer16_SingleTbar_tW_5FS_powheg'                  : genprocs_mutau_s_top,
'MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo'      : genprocs_mutau_s_top,
'MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg' : genprocs_mutau_s_top,
'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg'     : genprocs_mutau_s_top,
'MC2016_Summer16_TTJets_powheg'  : genprocs_mutau_tt,
'MC2016_Summer16_WWTo2L2Nu_powheg'             : genprocs_mutau_dibosons,
'MC2016_Summer16_WWToLNuQQ_powheg'             : genprocs_mutau_dibosons,
'MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin' : genprocs_mutau_dibosons,
'MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin'   : genprocs_mutau_dibosons,
'MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin'    : genprocs_mutau_dibosons,
'MC2016_Summer16_WZTo3LNu_powheg'              : genprocs_mutau_dibosons,
'MC2016_Summer16_ZZTo2L2Nu_powheg'             : genprocs_mutau_dibosons,
'MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin'    : genprocs_mutau_dibosons,
'MC2016_Summer16_W1Jets_madgraph'  : genprocs_mutau_wjets,
'MC2016_Summer16_W2Jets_madgraph'  : genprocs_mutau_wjets,
'MC2016_Summer16_W3Jets_madgraph'  : genprocs_mutau_wjets,
'MC2016_Summer16_W4Jets_madgraph'  : genprocs_mutau_wjets,
'MC2016_Summer16_WJets_madgraph'   : genprocs_mutau_wjets,
'MC2016_Summer16_WJets_amcatnlo'   : genprocs_mutau_wjets,
}


if args.custom_range:
    command = """python sumup_ttree_distrs.py "{draw_com}" --ttree ttree_out """ + "--custom-range {}".format(args.custom_range) + """ --cond-com "{selection}" --output test1_lep_pt_{dtag}_{proc}.root  --histo-name {chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}  --save-weight --per-weight {dtag_file}"""
else:
    command = """python sumup_ttree_distrs.py "{draw_com}" --ttree ttree_out """ + "--histo-range {}".format(args.histo_range) + """ --cond-com "{selection}" --output test1_lep_pt_{dtag}_{proc}.root  --histo-name {chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}  --save-weight --per-weight {dtag_file}"""

for input_file in args.input_files:
    matching_dtags = [dtag for dtag in dtags.keys() if dtag in input_file]
    if not matching_dtags:
        logging.warning("missing dtag for: " + input_file)
        continue

    assert len(matching_dtags) == 1
    assert isfile(input_file) # checks if a file (not directory or anything else) with this name exists

    dtag = matching_dtags[0]

    main_name, proc_defs = dtags[dtag]
    included_ids = []

    # if there are processes defined for this dtag
    # add the "other" process, which will automatically include all not included proc ids
    if proc_defs:
        proc_defs.append(('other', []))
    # otherwise, like in case of qcd of dibosons, all events are treated as 1 process

    for proc_name, proc_ids in proc_defs:
        # check that new ids have not been already processed
        assert not any(new_id in included_ids for new_id in proc_ids)
        # now save the ids of this process in the included
        for an_id in proc_ids:
            included_ids.append(an_id)

        # join the ids of this process in 1 string selection command
        # example: (gen_proc_id == 11 || gen_proc_id == 12 || gen_proc_id == 21)
        if not proc_name == 'other':
            proc_selection =  '(%s)' % (" || ".join('gen_proc_id == %d' % an_id for an_id in proc_ids))
        else:
            proc_selection = '!(%s)' % (" || ".join('gen_proc_id == %d' % an_id for an_id in included_ids))

        full_selection = args.cond_com + ' && ' + proc_selection
        if args.weight:
            full_selection = "%s * (%s)" % (args.weight, full_selection)

        proc_command = command.format(draw_com = args.draw_com,
                                      selection = full_selection,
                                      dtag = dtag,
                                      chan  = args.chan,
                                      proc  = main_name + '_' + proc_name,
                                      sys   = args.sys,
                                      distr = args.distr,
                                      dtag_file = input_file)

        logging.debug(proc_command)

        # launch the command
        #command_list = [operand for operand in proc_command.split()]
        #logging.debug(repr(command_list))
        #call(command_list)
        # ok, this is too fancy

        os.system(proc_command)

