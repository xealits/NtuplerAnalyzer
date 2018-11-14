import argparse
import logging
from os.path import isfile, basename
from sys import exit
import ctypes
from sys import exit




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "sumup TTree Draw",
    epilog = """Example:\npython sumup_ttree_draw.py "met_init.pt()" --histo-range 200,0,200 --histo-name data_met_init --output data_met_init.root gstore_outdirs/v28/SingleMuon/Ntupler_v28_Data13TeV_SingleMuon2016*/1*/*/*root
python sumup_ttree_draw.py "event_leptons[0].pt()" --ttree ttree_out --cond-com "selection_stage == 5" --histo-range 50,0,200 --output test1_lep_pt.root --histo-name foo/bar/test1_lep_pt --save-weight MC2016_Summer16_TTJets_powheg_test1.root """
    )

parser.add_argument('draw_com', type=str, help='Draw("??")')
parser.add_argument('--cond-com', type=str, default="", help='Draw("", "??")')
parser.add_argument('--ttree',    type=str, default="ntupler/reduced_ttree", help='path to the TTree in the file')
parser.add_argument('--histo-name',  type=str, default="out_histo", help='the ROOTName for output')
parser.add_argument('--std-histo-name',  type=str, default="distr", help='construct standard name for histogram: take "path" part from --histo-name and "distr" part from here, attach as chan_proc_sys_distr')
parser.add_argument('--histo-color', type=str, default=None, help='optional rgb color, like `255,255,255`')
parser.add_argument('--histo-range',  type=str, default=None, help='optionally set the range')
parser.add_argument('--custom-range', type=str, default=None, help='optionally set the range in custom bins')
parser.add_argument("--try-xsec",  action='store_true', help="try to normalize to cross-section, if the output filename contains dtag")

parser.add_argument("--std-histos", action='store_true', help="parse standard histo names, support standard channels, processes and systematics")
parser.add_argument("--el-procs",   action='store_true', help="resolving the std procs pick up the el-tau selection")
parser.add_argument("--cut-w0jets",   action='store_true', help="resolving the std procs pick up the el-tau selection")

parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")
parser.add_argument("--test",   action='store_true', help="test run, don't draw and don't write")
parser.add_argument("--output", type=str, default="output.root", help="filename for output")
parser.add_argument("--force",  action='store_true', help="force rewriting the output file if it exists")

parser.add_argument("--save-weight", action='store_true', help="save the weight counters in the output file")
parser.add_argument("--per-weight",  action='store_true', help="normalize by event weight of datasets")
parser.add_argument("--scale", type=float, help="the value to scale histo")
parser.add_argument("--scan",        action='store_true', help="also scan events ant print out")
parser.add_argument("--get-maximum", action='store_true', help="just find maxima on the draw_com")

parser.add_argument('input_files', nargs='+', help="""the files to sum up, passed by shell, as:

/gstore/t3cms/store/user/otoldaie/v19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v19_MC2016_Summer16_TTJets_powheg/180226_022336/0000/MC2016_Summer16_TTJets_powheg_*.root""")


args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

if isfile(args.output) and not args.force:
    print "the file already exists: %s" % args.output
    exit(0)

# [/]histo_path/histo_name
histo_name = args.histo_name.split('/')[-1]
histo_path = args.histo_name.split('/')[:-1]
# skip empty parts in the path (sequences of ////)
histo_path = [part for part in histo_path if part]

distr_name = histo_name
if args.std_histos:
    # 4 strings for names and they are not empty
    assert len(args.histo_name.split('/')) == 4 and all(args.histo_name.split('/'))
    channels, procs, systs, distr_name = args.histo_name.split('/')
    #histo_path = [chan, proc, sys]
    # later proc and sys are expanded to standard definitions
    # if proc == 'std_procs': ...
    # if sys  == 'std_sys': ,,,
    # there will be many histo_path-s for each definition




logging.info("import ROOT")

import ROOT
from ROOT import TH1D, TFile, TTree, gROOT
from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
from ROOT import TLegend
#from ROOT import THStack

logging.info("done")

from plotting_root import rgb

gROOT.SetBatch()



input_files = args.input_files


'''
out_histo  = TH1D("os_tauall", "", 44, -2, 20)
tau_b      = TH1D("os_taub", "", 44, -2, 20)
tau_w      = TH1D("os_tauw", "", 44, -2, 20)
tau_w_c    = TH1D("os_tauw_c", "", 44, -2, 20)
tau_w_notc = TH1D("os_tauw_notc", "", 44, -2, 20)

ss_tau_all    = TH1D("ss_tauall", "", 44, -2, 20)
ss_tau_b      = TH1D("ss_taub", "", 44, -2, 20)
ss_tau_w      = TH1D("ss_tauw", "", 44, -2, 20)
ss_tau_w_c    = TH1D("ss_tauw_c", "", 44, -2, 20)
ss_tau_w_notc = TH1D("ss_tauw_notc", "", 44, -2, 20)
'''


if args.histo_color:
    logging.debug("color: " + args.histo_color)
else:
    logging.debug("no color")

'''
root [4] double pt_bins[] = {20, 30, 40, 50, 70, 90, 120, 150, 200, 300, 2000};
root [5] 
root [5] sizeof
ROOT_prompt_5:2:1: error: expected expression
;
^
root [6] sizeof(pt_bins)
(unsigned long) 88
root [7] 
root [7] sizeof(double)
(unsigned long) 8
root [8] 
root [8] int pt_bins_n = 10
(int) 10
root [9] 
root [9] 
root [9] TH1D* hpts = new TH1D("hpts", "", pt_bins_n, pt_bins)
(TH1D *) 0x4db6750
root [10] hpts
(TH1D *) 0x4db6750
root [11] 
root [11] 
root [11] ttree_out->Draw("event_leptons[0].pt()>>hpts")
'''

if args.try_xsec:
    from per_dtag_script import dtags

# Draw command
# 
temp_output_histo = None # histo-template for custom bins

draw_command = '{met_lep_mt_var}' if args.draw_com == 'std_mt_vars' else args.draw_com
if args.custom_range:
    bin_edges = [float(b) for b in args.custom_range.split(',')]
    n_bin_edges = len(bin_edges)
    n_bins = n_bin_edges - 1
    logging.debug("making %d bins from %s" % (n_bins, args.custom_range))

    # first method
    root_bin_edges = (ctypes.c_double * n_bin_edges)(*bin_edges)
    temp_output_histo = TH1D("histotemp", "", n_bins, root_bin_edges) # root commands can access it by the name

    draw_command = draw_command + '>>' + "histotemp"
elif args.histo_range:
    draw_command = draw_command + ">>h(%s)" % args.histo_range
else:
    # TOFIX: without explicit range the histos won't add up
    draw_command = draw_command

logging.debug("draw: " + draw_command)
logging.debug("cond: " + args.cond_com)
if temp_output_histo:
    logging.debug("temp: " + temp_output_histo.GetName())
    #temp_output_histo.SetDirectory(0)
    # and some preparation for root stuff
    ROOT.gROOT.ProcessLine('TFile* intfile;')
    ROOT.gROOT.ProcessLine('TTree* inttree;')

out_histo = None
weight_counter = None



'''
contruct the output histos
for standard outputs
definitions of:
* sub-processes
* main selections
* systematics

condition string is defined by chan, proc and systematic like:
(chan && proc) * sys_weight

and I use the chan in 1 of systematics: LEP el or LEP mu


for a given dtag (a file or a bunch of them)
all chan/proc/sys apply independently
-- a list of (name, 'def_string') for each


if the std-histos option is given
the chan, proc, sys are parsed for names of standard defs or special names as 'std_procs' or 'all'

the goal is construction of all condition strings and loop over them for each file
the condition string consists of (<args.cond> && <def from std procs> && <std selections>) * <std weight>

-- std conditions have a loop of processes
   and of object systematics
   std weights loop weight systs

the loop is done over processes and systematics
for systematics there is an if for selection (objects) or weights
'''


all_std_channels = {
'mu_sel':          '{selection_stage}==107',
'mu_sel_ss':       '{selection_stage}==106',
'el_sel':          '{selection_stage}==117',
'el_sel_ss':       '{selection_stage}==116',
'mu_selVloose':    '({selection_stage}==107 || {selection_stage}==105)',
'mu_selVloose_ss': '({selection_stage}==106 || {selection_stage}==104)',
'el_selVloose':    '({selection_stage}==117 || {selection_stage}==115)',
'el_selVloose_ss': '({selection_stage}==116 || {selection_stage}==114)',

'mu_sel_ljout':          '{selection_stage}==107 && event_jets_lj_var >  60.',
'mu_sel_ljout_ss':       '{selection_stage}==106 && event_jets_lj_var >  60.',
'el_sel_ljout':          '{selection_stage}==117 && event_jets_lj_var >  60.',
'el_sel_ljout_ss':       '{selection_stage}==116 && event_jets_lj_var >  60.',
'mu_selVloose_ljout':    '({selection_stage}==107 || {selection_stage}==105) && event_jets_lj_var >  60.',
'mu_selVloose_ljout_ss': '({selection_stage}==106 || {selection_stage}==104) && event_jets_lj_var >  60.',
'el_selVloose_ljout':    '({selection_stage}==117 || {selection_stage}==115) && event_jets_lj_var >  60.',
'el_selVloose_ljout_ss': '({selection_stage}==116 || {selection_stage}==114) && event_jets_lj_var >  60.',

'mu_sel_lj':          '{selection_stage}==107 && event_jets_lj_var <= 60.',
'mu_sel_lj_ss':       '{selection_stage}==106 && event_jets_lj_var <= 60.',
'el_sel_lj':          '{selection_stage}==117 && event_jets_lj_var <= 60.',
'el_sel_lj_ss':       '{selection_stage}==116 && event_jets_lj_var <= 60.',
'mu_selVloose_lj':    '({selection_stage}==107 || {selection_stage}==105) && event_jets_lj_var <= 60.',
'mu_selVloose_lj_ss': '({selection_stage}==106 || {selection_stage}==104) && event_jets_lj_var <= 60.',
'el_selVloose_lj':    '({selection_stage}==117 || {selection_stage}==115) && event_jets_lj_var <= 60.',
'el_selVloose_lj_ss': '({selection_stage}==116 || {selection_stage}==114) && event_jets_lj_var <= 60.',
}

# new selection stages
all_std_channels = {
'mu_sel':          ('{selection_stage}==  7', 'selection_stage'),
'mu_sel_ss':       ('{selection_stage}==  6', 'selection_stage'),
'el_sel':          ('{selection_stage}== 17', 'selection_stage'),
'el_sel_ss':       ('{selection_stage}== 16', 'selection_stage'),
'mu_selVloose':    ('({selection_stage}==  7 || {selection_stage}==  5)', 'selection_stage'),
'mu_selVloose_ss': ('({selection_stage}==  6 || {selection_stage}==  4)', 'selection_stage'),
'el_selVloose':    ('({selection_stage}== 17 || {selection_stage}== 15)', 'selection_stage'),
'el_selVloose_ss': ('({selection_stage}== 16 || {selection_stage}== 14)', 'selection_stage'),

'mu_sel_ljout':          ('{selection_stage}==  7 && event_jets_lj_var >  60.', 'selection_stage'),
'mu_sel_ljout_ss':       ('{selection_stage}==  6 && event_jets_lj_var >  60.', 'selection_stage'),
'el_sel_ljout':          ('{selection_stage}== 17 && event_jets_lj_var >  60.', 'selection_stage'),
'el_sel_ljout_ss':       ('{selection_stage}== 16 && event_jets_lj_var >  60.', 'selection_stage'),
'mu_selVloose_ljout':    ('({selection_stage}==  7 || {selection_stage}==  5) && event_jets_lj_var >  60.', 'selection_stage'),
'mu_selVloose_ljout_ss': ('({selection_stage}==  6 || {selection_stage}==  4) && event_jets_lj_var >  60.', 'selection_stage'),
'el_selVloose_ljout':    ('({selection_stage}== 17 || {selection_stage}== 15) && event_jets_lj_var >  60.', 'selection_stage'),
'el_selVloose_ljout_ss': ('({selection_stage}== 16 || {selection_stage}== 14) && event_jets_lj_var >  60.', 'selection_stage'),

'mu_sel_lj':          ('{selection_stage}==  7 && event_jets_lj_var <= 60.', 'selection_stage'),
'mu_sel_lj_ss':       ('{selection_stage}==  6 && event_jets_lj_var <= 60.', 'selection_stage'),
'el_sel_lj':          ('{selection_stage}== 17 && event_jets_lj_var <= 60.', 'selection_stage'),
'el_sel_lj_ss':       ('{selection_stage}== 16 && event_jets_lj_var <= 60.', 'selection_stage'),
'mu_selVloose_lj':    ('({selection_stage}==  7 || {selection_stage}==  5) && event_jets_lj_var <= 60.', 'selection_stage'),
'mu_selVloose_lj_ss': ('({selection_stage}==  6 || {selection_stage}==  4) && event_jets_lj_var <= 60.', 'selection_stage'),
'el_selVloose_lj':    ('({selection_stage}== 17 || {selection_stage}== 15) && event_jets_lj_var <= 60.', 'selection_stage'),
'el_selVloose_lj_ss': ('({selection_stage}== 16 || {selection_stage}== 14) && event_jets_lj_var <= 60.', 'selection_stage'),

# additional channels
'dy_mutau': ('({selection_stage}== 102 || {selection_stage}== 103)', 'selection_stage_dy'),
'dy_eltau': ('({selection_stage}== 112 || {selection_stage}== 113)', 'selection_stage_dy'),
'dy_mumu':  ('({selection_stage}== 102 || {selection_stage}== 103 || {selection_stage}== 105)', 'selection_stage_dy_mumu'),
'dy_elel':  ('({selection_stage}== 112 || {selection_stage}== 113 || {selection_stage}== 115)', 'selection_stage_dy_mumu'),

'tt_elmu':  ('({selection_stage}== 205)', 'selection_stage_em'),
}


if args.std_histos:
    if   channels == 'all':
         channels = ','.join(all_std_channels.keys())
    elif channels == 'mu_std':
         channels = 'mu_sel_ljout,mu_sel_ljout_ss,mu_selVloose_ljout,mu_selVloose_ljout_ss,mu_sel_lj,mu_sel_lj_ss,mu_selVloose_lj,mu_selVloose_lj_ss'
    elif channels == 'el_std':
         channels = 'el_sel_ljout,el_sel_ljout_ss,el_selVloose_ljout,el_selVloose_ljout_ss,el_sel_lj,el_sel_lj_ss,el_selVloose_lj,el_selVloose_lj_ss'
    elif channels == 'all_std':
         channels = 'mu_sel_ljout,mu_sel_ljout_ss,mu_selVloose_ljout,mu_selVloose_ljout_ss,mu_sel_lj,mu_sel_lj_ss,mu_selVloose_lj,mu_selVloose_lj_ss,el_sel_ljout,el_sel_ljout_ss,el_selVloose_ljout,el_selVloose_ljout_ss,el_sel_lj,el_sel_lj_ss,el_selVloose_lj,el_selVloose_lj_ss'

# object systematics change the selection
#"selection_stage_JERDown ${cond}" 'JERDown': event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root
#"selection_stage_JERUp   ${cond}" 'JERUp'  : event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root
#"selection_stage_JESDown ${cond}" 'JESDown': event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root
#"selection_stage_JESUp   ${cond}" 'JESUp'  : event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root
#"selection_stage_TESDown ${cond}" 'TESDown': event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root
#"selection_stage_TESUp   ${cond}" 'TESUp'  : event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root

systs_objects = {
'JERDown':  'selection_stage_JERDown',
'JERUp'  :  'selection_stage_JERUp'  ,
'JESDown':  'selection_stage_JESDown',
'JESUp'  :  'selection_stage_JESUp'  ,
'TESDown':  'selection_stage_TESDown',
'TESUp'  :  'selection_stage_TESUp'  ,
}

systs_objects_mt_variation = {
'JERUp'    : 'event_met_lep_mt_JERUp',   
'JERDown'  : 'event_met_lep_mt_JERDown', 
'JESUp'    : 'event_met_lep_mt_JESUp',   
'JESDown'  : 'event_met_lep_mt_JESDown', 
'TESUp'    : 'event_met_lep_mt_TESUp',   
'TESDown'  : 'event_met_lep_mt_TESDown', 
}

# check is a syst name is in objects
#    if so -> change the selection and pick up nominal weight
#    else  -> check it's in the weights

systs_weights_nominal = {
'NOMINAL': 'event_weight',
}

systs_weights_common = {
'PUUp'   : "event_weight*event_weight_PUUp/event_weight_PU"    ,
'PUDown' : "event_weight*event_weight_PUDown/event_weight_PU"    ,
'bSFUp'  : "event_weight*event_weight_bSFUp/event_weight_bSF"  ,
'bSFDown': "event_weight*event_weight_bSFDown/event_weight_bSF",

'LEPelIDUp'   : "event_weight*event_weight_LEPidUp"   ,
'LEPelIDDown' : "event_weight*event_weight_LEPidDown" ,
'LEPelTRGUp'  : "event_weight*event_weight_LEPtrgUp"  ,
'LEPelTRGDown': "event_weight*event_weight_LEPtrgDown",
'LEPmuIDUp'   : "event_weight*event_weight_LEPidUp"   ,
'LEPmuIDDown' : "event_weight*event_weight_LEPidDown" ,
'LEPmuTRGUp'  : "event_weight*event_weight_LEPtrgUp"  ,
'LEPmuTRGDown': "event_weight*event_weight_LEPtrgDown",
}

systs_weights_tt = {
'TOPPTDown'     : 'event_weight'                           ,
'TOPPTUp'       : 'event_weight*event_weight_toppt'        ,
'FragUp'        : 'event_weight*event_weight_FragUp'       ,
'FragDown'      : 'event_weight*event_weight_FragDown'     ,
'SemilepBRUp'   : 'event_weight*event_weight_SemilepBRUp'  ,
'SemilepBRDown' : 'event_weight*event_weight_SemilepBRDown',
'PetersonUp'    : 'event_weight*event_weight_PetersonUp'   ,
'PetersonDown'  : 'event_weight'                           ,
}

systs_weights_tt_hard = {
"MrUp"    : "event_weight*event_weight_me_f_rUp",
"MrDown"  : "event_weight*event_weight_me_f_rDn",
"MfUp"    : "event_weight*event_weight_me_fUp_r",
"MfDown"  : "event_weight*event_weight_me_fDn_r",
"MfrUp"   : "event_weight*event_weight_me_frUp" ,
"MfrDown" : "event_weight*event_weight_me_frDn" ,
}

systs_weights_tt_alpha = {
'AlphaSUp'  : "event_weight*event_weight_AlphaS_up",
'AlphaSDown': "event_weight*event_weight_AlphaS_dn",
}

systs_weights_tt_pdf1 = {
'PDFCT14n1Up'     : "event_weight*event_weight_pdf[0]" ,
'PDFCT14n2Up'     : "event_weight*event_weight_pdf[1]" ,
'PDFCT14n3Up'     : "event_weight*event_weight_pdf[2]" ,
'PDFCT14n4Up'     : "event_weight*event_weight_pdf[3]" ,
'PDFCT14n5Up'     : "event_weight*event_weight_pdf[4]" ,
'PDFCT14n6Up'     : "event_weight*event_weight_pdf[5]" ,
'PDFCT14n7Up'     : "event_weight*event_weight_pdf[6]" ,
'PDFCT14n8Up'     : "event_weight*event_weight_pdf[7]" ,
'PDFCT14n9Up'     : "event_weight*event_weight_pdf[8]" ,
'PDFCT14n10Up'    : "event_weight*event_weight_pdf[9]" ,
}

systs_weights_tt_pdf10 = {
'PDFCT14n11Up'    : "event_weight*event_weight_pdf[10]",
'PDFCT14n12Up'    : "event_weight*event_weight_pdf[11]",
'PDFCT14n13Up'    : "event_weight*event_weight_pdf[12]",
'PDFCT14n14Up'    : "event_weight*event_weight_pdf[13]",
'PDFCT14n15Up'    : "event_weight*event_weight_pdf[14]",
'PDFCT14n16Up'    : "event_weight*event_weight_pdf[15]",
'PDFCT14n17Up'    : "event_weight*event_weight_pdf[16]",
'PDFCT14n18Up'    : "event_weight*event_weight_pdf[17]",
'PDFCT14n19Up'    : "event_weight*event_weight_pdf[18]",
'PDFCT14n20Up'    : "event_weight*event_weight_pdf[19]",
}

systs_weights_tt_pdf20 = {
'PDFCT14n21Up'    : "event_weight*event_weight_pdf[20]",
'PDFCT14n22Up'    : "event_weight*event_weight_pdf[21]",
'PDFCT14n23Up'    : "event_weight*event_weight_pdf[22]",
'PDFCT14n24Up'    : "event_weight*event_weight_pdf[23]",
'PDFCT14n25Up'    : "event_weight*event_weight_pdf[24]",
'PDFCT14n26Up'    : "event_weight*event_weight_pdf[25]",
'PDFCT14n27Up'    : "event_weight*event_weight_pdf[26]",
'PDFCT14n28Up'    : "event_weight*event_weight_pdf[27]",
'PDFCT14n29Up'    : "event_weight*event_weight_pdf[28]",
'PDFCT14n30Up'    : "event_weight*event_weight_pdf[29]",
}

systs_weights_tt_pdf30 = {
'PDFCT14n31Up'    : "event_weight*event_weight_pdf[30]",
'PDFCT14n32Up'    : "event_weight*event_weight_pdf[31]",
'PDFCT14n33Up'    : "event_weight*event_weight_pdf[32]",
'PDFCT14n34Up'    : "event_weight*event_weight_pdf[33]",
'PDFCT14n35Up'    : "event_weight*event_weight_pdf[34]",
'PDFCT14n36Up'    : "event_weight*event_weight_pdf[35]",
'PDFCT14n37Up'    : "event_weight*event_weight_pdf[36]",
'PDFCT14n38Up'    : "event_weight*event_weight_pdf[37]",
'PDFCT14n39Up'    : "event_weight*event_weight_pdf[38]",
'PDFCT14n40Up'    : "event_weight*event_weight_pdf[39]",
}

systs_weights_tt_pdf40 = {
'PDFCT14n41Up'    : "event_weight*event_weight_pdf[40]",
'PDFCT14n42Up'    : "event_weight*event_weight_pdf[41]",
'PDFCT14n43Up'    : "event_weight*event_weight_pdf[42]",
'PDFCT14n44Up'    : "event_weight*event_weight_pdf[43]",
'PDFCT14n45Up'    : "event_weight*event_weight_pdf[44]",
'PDFCT14n46Up'    : "event_weight*event_weight_pdf[45]",
'PDFCT14n47Up'    : "event_weight*event_weight_pdf[46]",
'PDFCT14n48Up'    : "event_weight*event_weight_pdf[47]",
'PDFCT14n49Up'    : "event_weight*event_weight_pdf[48]",
'PDFCT14n50Up'    : "event_weight*event_weight_pdf[49]",
}

systs_weights_tt_pdf50 = {
'PDFCT14n51Up'    : "event_weight*event_weight_pdf[50]",
'PDFCT14n52Up'    : "event_weight*event_weight_pdf[51]",
'PDFCT14n53Up'    : "event_weight*event_weight_pdf[52]",
'PDFCT14n54Up'    : "event_weight*event_weight_pdf[53]",
'PDFCT14n55Up'    : "event_weight*event_weight_pdf[54]",
'PDFCT14n56Up'    : "event_weight*event_weight_pdf[55]",
}


systs_weights_tt_alpha_pdf = {
'AlphaSUp'  : "event_weight*event_weight_AlphaS_up",
'AlphaSDown': "event_weight*event_weight_AlphaS_dn",

'PDFCT14n1Up'     : "event_weight*event_weight_pdf[0]" ,
'PDFCT14n2Up'     : "event_weight*event_weight_pdf[1]" ,
'PDFCT14n3Up'     : "event_weight*event_weight_pdf[2]" ,
'PDFCT14n4Up'     : "event_weight*event_weight_pdf[3]" ,
'PDFCT14n5Up'     : "event_weight*event_weight_pdf[4]" ,
'PDFCT14n6Up'     : "event_weight*event_weight_pdf[5]" ,
'PDFCT14n7Up'     : "event_weight*event_weight_pdf[6]" ,
'PDFCT14n8Up'     : "event_weight*event_weight_pdf[7]" ,
'PDFCT14n9Up'     : "event_weight*event_weight_pdf[8]" ,
'PDFCT14n10Up'    : "event_weight*event_weight_pdf[9]" ,
'PDFCT14n11Up'    : "event_weight*event_weight_pdf[10]",
'PDFCT14n12Up'    : "event_weight*event_weight_pdf[11]",
'PDFCT14n13Up'    : "event_weight*event_weight_pdf[12]",
'PDFCT14n14Up'    : "event_weight*event_weight_pdf[13]",
'PDFCT14n15Up'    : "event_weight*event_weight_pdf[14]",
'PDFCT14n16Up'    : "event_weight*event_weight_pdf[15]",
'PDFCT14n17Up'    : "event_weight*event_weight_pdf[16]",
'PDFCT14n18Up'    : "event_weight*event_weight_pdf[17]",
'PDFCT14n19Up'    : "event_weight*event_weight_pdf[18]",
'PDFCT14n20Up'    : "event_weight*event_weight_pdf[19]",
'PDFCT14n21Up'    : "event_weight*event_weight_pdf[20]",
'PDFCT14n22Up'    : "event_weight*event_weight_pdf[21]",
'PDFCT14n23Up'    : "event_weight*event_weight_pdf[22]",
'PDFCT14n24Up'    : "event_weight*event_weight_pdf[23]",
'PDFCT14n25Up'    : "event_weight*event_weight_pdf[24]",
'PDFCT14n26Up'    : "event_weight*event_weight_pdf[25]",
'PDFCT14n27Up'    : "event_weight*event_weight_pdf[26]",
'PDFCT14n28Up'    : "event_weight*event_weight_pdf[27]",
'PDFCT14n29Up'    : "event_weight*event_weight_pdf[28]",
'PDFCT14n30Up'    : "event_weight*event_weight_pdf[29]",
'PDFCT14n31Up'    : "event_weight*event_weight_pdf[30]",
'PDFCT14n32Up'    : "event_weight*event_weight_pdf[31]",
'PDFCT14n33Up'    : "event_weight*event_weight_pdf[32]",
'PDFCT14n34Up'    : "event_weight*event_weight_pdf[33]",
'PDFCT14n35Up'    : "event_weight*event_weight_pdf[34]",
'PDFCT14n36Up'    : "event_weight*event_weight_pdf[35]",
'PDFCT14n37Up'    : "event_weight*event_weight_pdf[36]",
'PDFCT14n38Up'    : "event_weight*event_weight_pdf[37]",
'PDFCT14n39Up'    : "event_weight*event_weight_pdf[38]",
'PDFCT14n40Up'    : "event_weight*event_weight_pdf[39]",
'PDFCT14n41Up'    : "event_weight*event_weight_pdf[40]",
'PDFCT14n42Up'    : "event_weight*event_weight_pdf[41]",
'PDFCT14n43Up'    : "event_weight*event_weight_pdf[42]",
'PDFCT14n44Up'    : "event_weight*event_weight_pdf[43]",
'PDFCT14n45Up'    : "event_weight*event_weight_pdf[44]",
'PDFCT14n46Up'    : "event_weight*event_weight_pdf[45]",
'PDFCT14n47Up'    : "event_weight*event_weight_pdf[46]",
'PDFCT14n48Up'    : "event_weight*event_weight_pdf[47]",
'PDFCT14n49Up'    : "event_weight*event_weight_pdf[48]",
'PDFCT14n50Up'    : "event_weight*event_weight_pdf[49]",
'PDFCT14n51Up'    : "event_weight*event_weight_pdf[50]",
'PDFCT14n52Up'    : "event_weight*event_weight_pdf[51]",
'PDFCT14n53Up'    : "event_weight*event_weight_pdf[52]",
'PDFCT14n54Up'    : "event_weight*event_weight_pdf[53]",
'PDFCT14n55Up'    : "event_weight*event_weight_pdf[54]",
'PDFCT14n56Up'    : "event_weight*event_weight_pdf[55]",
}

systs_weights_tt_updowns = {
'TuneCUETP8M2T4Down' : 'event_weight',
'TuneCUETP8M2T4Up'   : 'event_weight',
'FSRDown'            : 'event_weight',
'FSRUp'              : 'event_weight',
'HDAMPDown'          : 'event_weight',
'HDAMPUp'            : 'event_weight',
'ISRDown'            : 'event_weight',
'ISRUp'              : 'event_weight',
}

updown_dtags_to_sys = {
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4down' : 'TuneCUETP8M2T4Down',
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4up'   : 'TuneCUETP8M2T4Up',
'MC2016_Summer16_TTJets_powheg_fsrdown'        : 'FSRDown',
'MC2016_Summer16_TTJets_powheg_fsrup'          : 'FSRUp',
'MC2016_Summer16_TTJets_powheg_hdampDOWN'      : 'HDAMPDown',
'MC2016_Summer16_TTJets_powheg_hdampUP'        : 'HDAMPUp',
'MC2016_Summer16_TTJets_powheg_isrdown'        : 'ISRDown',
'MC2016_Summer16_TTJets_powheg_isrup'          : 'ISRUp',
}

named_systs_weights_all = {'nom': systs_weights_nominal,
'common': systs_weights_common,
'tt_weights': systs_weights_tt,
'tt_hard':  systs_weights_tt_hard,
'tt_pdf':   systs_weights_tt_alpha_pdf,
'tt_alpha': systs_weights_tt_alpha,
'tt_pdf1':  systs_weights_tt_pdf1,
'tt_pdf10': systs_weights_tt_pdf10,
'tt_pdf20': systs_weights_tt_pdf20,
'tt_pdf30': systs_weights_tt_pdf30,
'tt_pdf40': systs_weights_tt_pdf40,
'tt_pdf50': systs_weights_tt_pdf50,
}


systs_weights_all = {}
for s_d in named_systs_weights_all.values():
    systs_weights_all.update(s_d)

systs_weights_all.update(systs_weights_tt_updowns)

#systs = std_systematics_per_dtag[dtag]
usual_mc_systs = 'nom,common,obj'
std_systematics_per_dtag = {
'MC2016_Summer16_QCD_HT-100-200'    : 'nom',
'MC2016_Summer16_QCD_HT-1000-1500'  : 'nom',
'MC2016_Summer16_QCD_HT-1500-2000'  : 'nom',
'MC2016_Summer16_QCD_HT-200-300'    : 'nom',
'MC2016_Summer16_QCD_HT-2000-Inf'   : 'nom',
'MC2016_Summer16_QCD_HT-300-500'    : 'nom',
'MC2016_Summer16_QCD_HT-500-700'    : 'nom',
'MC2016_Summer16_QCD_HT-700-1000'   : 'nom',

'MC2016_Summer16_SingleT_tW_5FS_powheg'                     : usual_mc_systs,
'MC2016_Summer16_SingleTbar_tW_5FS_powheg'                  : usual_mc_systs,
'MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo'      : usual_mc_systs,
'MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg' : usual_mc_systs,
'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg'     : usual_mc_systs,

'MC2016_Summer16_TTJets_powheg'  : 'nom,common,obj,tt_weights,tt_hard,tt_pdf',
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4down' : 'NOMINAL',
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4up'   : 'NOMINAL',
'MC2016_Summer16_TTJets_powheg_fsrdown'        : 'NOMINAL',
'MC2016_Summer16_TTJets_powheg_fsrup'          : 'NOMINAL',
'MC2016_Summer16_TTJets_powheg_hdampDOWN'      : 'NOMINAL',
'MC2016_Summer16_TTJets_powheg_hdampUP'        : 'NOMINAL',
'MC2016_Summer16_TTJets_powheg_isrdown'        : 'NOMINAL',
'MC2016_Summer16_TTJets_powheg_isrup'          : 'NOMINAL',

'MC2016_Summer16_WWTo2L2Nu_powheg'             : usual_mc_systs,
'MC2016_Summer16_WWToLNuQQ_powheg'             : usual_mc_systs,
'MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin' : usual_mc_systs,
'MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin'   : usual_mc_systs,
'MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin'    : usual_mc_systs,
'MC2016_Summer16_WZTo3LNu_powheg'              : usual_mc_systs,
'MC2016_Summer16_ZZTo2L2Nu_powheg'             : usual_mc_systs,
'MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin'    : usual_mc_systs,

'MC2016_Summer16_W1Jets_madgraph'  : usual_mc_systs,
'MC2016_Summer16_W2Jets_madgraph'  : usual_mc_systs,
'MC2016_Summer16_W3Jets_madgraph'  : usual_mc_systs,
'MC2016_Summer16_W4Jets_madgraph'  : usual_mc_systs,
'MC2016_Summer16_WJets_madgraph'   : usual_mc_systs,
'MC2016_Summer16_WJets_amcatnlo'   : usual_mc_systs,

'MC2016_Summer16_DYJetsToLL_10to50_amcatnlo'   : usual_mc_systs,
'MC2016_Summer16_DYJetsToLL_50toInf_madgraph'  : usual_mc_systs,

'SingleMuon'     : 'NOMINAL',
'SingleElectron' : 'NOMINAL',
}


# processes

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

genprocs_eltau_s_top = ('s_top', [('eltau', [genproc_stop_el]),  ('lj', [genproc_stop_lj])])
genprocs_eltau_tt    = ('tt',    [('eltau', [genproc_tt_eltau, genproc_tt_eltau3ch]),
                                  ('lj', [genproc_tt_lj, genproc_tt_ljb, genproc_tt_ljw, genproc_tt_ljo]),
                                  ('taultauh', [genproc_tt_taultauh]),
                                  ('taulj', [genproc_tt_taulj])])

dtags_procs = {
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

'MC2016_Summer16_QCD_MuEnriched_Pt15_20toInf'  : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_1000toInf' : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_120to170'  : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_15to20'    : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_170to300'  : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_20to30'    : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_300to470'  : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_30to50'    : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_470to600'  : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_50to80'    : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_600to800'  : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_800to1000' : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_MuEnriched_Pt5_80to120'   : genprocs_mutau_qcd,

'MC2016_Summer16_QCD_EMEnriched_Pt-120to170.root' : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-170to300.root' : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-20to30.root'   : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-300toInf.root' : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-30to40.root'   : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-30to50.root'   : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-30toInf.root'  : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-40toInf.root'  : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-50to80.root'   : genprocs_mutau_qcd,
'MC2016_Summer16_QCD_EMEnriched_Pt-80to120.root'  : genprocs_mutau_qcd,

'MC2016_Summer16_SingleT_tW_5FS_powheg'                     : genprocs_mutau_s_top,
'MC2016_Summer16_SingleTbar_tW_5FS_powheg'                  : genprocs_mutau_s_top,
'MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo'      : genprocs_mutau_s_top,
'MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg' : genprocs_mutau_s_top,
'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg'     : genprocs_mutau_s_top,

'MC2016_Summer16_TTJets_powheg'  : genprocs_mutau_tt,
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4down' : genprocs_mutau_tt,
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4up'   : genprocs_mutau_tt,
'MC2016_Summer16_TTJets_powheg_fsrdown'        : genprocs_mutau_tt,
'MC2016_Summer16_TTJets_powheg_fsrup'          : genprocs_mutau_tt,
'MC2016_Summer16_TTJets_powheg_hdampDOWN'      : genprocs_mutau_tt,
'MC2016_Summer16_TTJets_powheg_hdampUP'        : genprocs_mutau_tt,
'MC2016_Summer16_TTJets_powheg_isrdown'        : genprocs_mutau_tt,
'MC2016_Summer16_TTJets_powheg_isrup'          : genprocs_mutau_tt,

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

'SingleMuon'     : ('data', []),
'SingleElectron' : ('data', []),
}

dtags_procs_el = {
'MC2016_Summer16_TTJets_powheg'               : genprocs_eltau_tt,
'MC2016_Summer16_SingleT_tW_5FS_powheg'       : genprocs_eltau_s_top,
'MC2016_Summer16_SingleTbar_tW_5FS_powheg'    : genprocs_eltau_s_top,
'MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo' : genprocs_eltau_s_top,
'MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg' : genprocs_eltau_s_top,
'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg' : genprocs_eltau_s_top,
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4down':  genprocs_eltau_tt,
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4up'  :  genprocs_eltau_tt,
'MC2016_Summer16_TTJets_powheg_fsrdown'       :  genprocs_eltau_tt,
'MC2016_Summer16_TTJets_powheg_fsrup'         :  genprocs_eltau_tt,
'MC2016_Summer16_TTJets_powheg_hdampDOWN'     :  genprocs_eltau_tt,
'MC2016_Summer16_TTJets_powheg_hdampUP'       :  genprocs_eltau_tt,
'MC2016_Summer16_TTJets_powheg_isrdown'       :  genprocs_eltau_tt,
'MC2016_Summer16_TTJets_powheg_isrup'         :  genprocs_eltau_tt,
}



dtag = None
if args.std_histos:
    # set processes
    if procs == 'std_procs':
        # find dtag and procs by the first input files

        # first check TT systematics
        for d in updown_dtags_to_sys:
            if d in input_files[0]:
                dtag = d
                break

        if not dtag:
          for d in dtags_procs:
            if d in input_files[0]:
                dtag = d
                break

        if not dtag:
            raise ValueError("for std_procs could not find the dtag for %s" % input_files[0])

    logging.debug("dtag, %s" % dtag)
    #procs = dtags_procs[dtag]

    # unfold systematics and weights
    # used per_dtag for standard systematics
    if systs == 'per_dtag':
        for d in std_systematics_per_dtag:
            if d in dtag:
                systs = std_systematics_per_dtag[dtag]
                break
        # error if still per_dtag
        if systs == 'per_dtag':
            raise ValueError('per_dtag, not found dtag for %s' % dtag)

    systs_in = systs.split(',')

    # if it is TT systematic dataset -- substitute the syst correctly
    if dtag in updown_dtags_to_sys:
        systs_in = [updown_dtags_to_sys[dtag]]

    systs = {}
    for sys in systs_in:
        if sys in named_systs_weights_all:
            systs.update(named_systs_weights_all[sys])
        elif sys == 'obj':
            systs.update(systs_objects)
        elif sys in systs_weights_all:
            systs[sys] = systs_weights_all[sys]
        elif sys in systs_objects:
            systs[sys] = systs_objects[sys]
        else:
            print "skipping %s" % sys

    isWJetsInclusive = dtag == "MC2016_Summer16_WJets_madgraph"

logging.debug('dtag = %s' % dtag)

if args.get_maximum:
   maximum = -1111111.


# draw and add to the output histos
output_histos = {} # (chan, proc, sys): histo
def draw_and_save(ttree, def_tuple, draw_command, condition_string, test):
    logging.debug(','.join(def_tuple))
    logging.debug('Draw(%s, %s)' % (draw_command, condition_string))

    if test:
        return None

    ttree.Draw(draw_command, condition_string)

    if temp_output_histo:
        # if there is temp histo then the histo was written there
        logging.debug(temp_output_histo.Integral())
        logging.debug(ROOT.histotemp.Integral())
        histo = temp_output_histo
    else:
        histo = ttree.GetHistogram()
        histo.SetDirectory(0)

    # handle errors
    histo.Sumw2()

    if def_tuple not in output_histos:
        out_histo = histo.Clone()
        #out_histo.SetName('_'.join(*def_tuple, distr_name))
        out_histo.SetName('_'.join([str(i) for i in def_tuple] + [str(distr_name)]))

        out_histo.SetDirectory(0)
        out_histo.Sumw2()
        output_histos[def_tuple] = out_histo
    else:
        output_histos[def_tuple].Add(histo)


for filename in input_files:
    if not isfile(filename):
        logging.info("missing: " + filename)
        continue

    logging.debug(filename)

    tfile = TFile(filename)
    ttree = tfile.Get(args.ttree)
    if temp_output_histo:
        temp_output_histo.SetDirectory(tfile)
        # in ROOT the TTree.Draw command "sees" only histograms in current working directory
        # after opening new file ROOT moves working directory into that file
        # therefore move temp tehre too
        # probably it moving current directory to some 0 directory could work
        # and as probable is the possibility to bring more troubles that way

    # Draw the file and sum up
    # 
    # TOFIX: without explicit range the histos won't add up, need to pick up the range of the first histo and propagate to the following?

    if args.std_histos:
        # TODO loop here
        #condition_strings = args.cond_com

        for chan in channels.split(','):
          # the channels defines the formula for the selection
          # which needs the object-based selection index

          if 'el_sel' in chan and dtag in dtags_procs_el:
              main_name, proc_defs = dtags_procs_el[dtag]
          else:
              main_name, proc_defs = dtags_procs[dtag]

          logging.debug("channel procs, %s %s" % (main_name, repr(proc_defs)))

          if proc_defs:
              proc_defs.append(('other', []))
              # the 'other' procs pick upp the not included ids
          #else:
          #    # it's 1 inclusive process
          #    proc_defs.append((None, []))

          included_ids = []
          for proc_name, proc_ids in (proc_defs if proc_defs else [(None, [])]):
            logging.debug(repr(proc_name))
            # check that new ids have not been already processed
            assert not any(new_id in included_ids for new_id in proc_ids)
            # now save the ids of this process in the included
            for an_id in proc_ids:
                included_ids.append(an_id)

            # join the ids of this process in 1 string selection command
            # example: (gen_proc_id == 11 || gen_proc_id == 12 || gen_proc_id == 21)
            if proc_name is None:
                proc_selection = None
            elif not proc_name == 'other':
                proc_selection =  '(%s)' % (" || ".join('gen_proc_id == %d' % an_id for an_id in proc_ids))
            else:
                proc_selection = '!(%s)' % (" || ".join('gen_proc_id == %d' % an_id for an_id in included_ids))

            if proc_selection is None:
                conditions = []
            else:
                conditions = [proc_selection]

            logging.debug("conditions %s" % repr(conditions))

            if args.cond_com:
                conditions.append(args.cond_com)

            if args.cut_w0jets and isWJetsInclusive:
                conditions.append('nup < 6')
            #if args.weight:
            #    condition_string = "%s * (%s)" % (args.weight, condition_string)

            logging.debug("conditions %s" % repr(conditions))
            #
            for sys_name in systs:
                # if it is object systematics -> change the selection and multiply nominal weight
                # else don't change the selection
                # = get the selection index name and the syst weight
                draw_command_final = draw_command
                ## it is the object systematic -- the selection index changes and the met variation, which propagates to precomputed mt
                ##selection_stage = systs_objects[sys_name]
                ##selection_stage = 'selection_stage'
                #if args.draw_com == 'std_mt_vars':
                #    draw_command_final = draw_command.format(met_lep_mt_var='event_met_lep_mt')
                #    logging.debug('substituted draw command to default %s at %s sys' % (draw_command_final, sys_name))

                if args.draw_com == 'std_mt_vars':
                    if sys_name in systs_objects:
                        draw_command_final = draw_command.format(met_lep_mt_var=systs_objects_mt_variation[sys_name])
                        logging.debug('substituted draw command to %s at %s sys' % (draw_command_final, sys_name))
                    else:
                        draw_command_final = draw_command.format(met_lep_mt_var='event_met_lep_mt')
                        logging.debug('substituted draw command to default %s at %s sys' % (draw_command_final, sys_name))

                final_cond = conditions[:]
                if chan in all_std_channels:
                    selection_stage = all_std_channels[chan][1]
                    final_cond.append(all_std_channels[chan][0].format(selection_stage=(selection_stage + '_' + sys_name if sys_name in systs_objects else selection_stage)))

                if dtag != 'data':
                    sys_weight = systs_weights_all.get(sys_name, systs_weights_nominal['NOMINAL'])
                    final_cond = '(%s) * %s' % (' && '.join(final_cond), sys_weight)
                else:
                    final_cond = '(%s)' % (' && '.join(final_cond))

                draw_and_save(ttree, (chan, (main_name + '_' + proc_name) if proc_name else main_name, sys_name), draw_command_final, final_cond, args.test)

    else:
        draw_and_save(ttree, tuple(histo_path), draw_command, args.cond_com, args.test)

    if args.test:
        continue

    if args.get_maximum:
        m = ttree.GetMaximum(args.draw_com)
        print "%30s %f" % (basename(filename), m)
        maximum = max(m, maximum)

    if args.scan:
        print filename
        ttree.Scan(args.draw_com, args.cond_com)

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

    # on closing the file all objects in it might be deleted
    # therefore move away the temp histo
    if temp_output_histo:
        temp_output_histo.SetDirectory(0)

    tfile.Close()

if args.test:
    exit(0)

# scale MC
if args.per_weight and dtag not in ('data', 'SingleMuon', 'SingleElectron'):
    #weight_counter = tfile.Get('ntupler/weight_counter')
    #out_histo.Scale(1./weight_counter.GetBinContent(2))
    for histo in output_histos.values():
        histo.Scale(1./weight_counter.GetBinContent(2))


# scale MC
if args.try_xsec and dtag not in ('data', 'SingleMuon', 'SingleElectron'):
    _, xsec = dtags.get(dtag, (None, 1.))
    #out_histo.Scale(xsec)
    for histo in output_histos.values():
        histo.Scale(xsec)

if args.scale:
    for histo in output_histos.values():
        histo.Scale(args.scale)

if args.histo_color:
    for histo in output_histos.values():
        histo.SetLineColor(rgb(*(int(i) for i in args.histo_color.split(','))))


if args.get_maximum:
   print "Max:", maximum

''' no legend for now
leg = TLegend(0.5, 0.5, 0.9, 0.9)
leg.SetName("legOS")

leg.AddEntry(tau_all    , "tt->lj, OS, all", 'L')
leg.AddEntry(tau_b      , "tt->lj, OS, b prod", 'L')
leg.AddEntry(tau_w      , "tt->lj, OS, w prod", 'L')
leg.AddEntry(tau_w_c    , "tt->lj, OS, w, C flav", 'L')
leg.AddEntry(tau_w_notc , "tt->lj, OS, w, not C", 'L')

leg2 = TLegend(0.5, 0.5, 0.9, 0.9)
leg2.SetName("legSS")

leg2.AddEntry(tau_all    , "tt->lj, SS, all", 'L')
leg2.AddEntry(tau_b      , "tt->lj, SS, b prod", 'L')
leg2.AddEntry(tau_w      , "tt->lj, SS, w prod", 'L')
leg2.AddEntry(tau_w_c    , "tt->lj, SS, w, C flav", 'L')
leg2.AddEntry(tau_w_notc , "tt->lj, SS, w, not C", 'L')
'''

fout = TFile(args.output, "RECREATE")
fout.Write()

fout.cd()

# corresponding to old protocol everything is saved as
# channel/process/systematic/channel_process_systematic_distr
# -- the final part is the name of the histogram, which must be unique in root
#    the rest is for convenience
# in principle, root handles many same name histos in different directories, but can run into some weirdness on practice
# that's why this protocol is still kept

# if trying the xsec -- plug the nickname too and construct standard name if needed
# -- no, now the xsec is just a scale factor
#if args.try_xsec and len(histo_path) > 1:
#    histo_path[1] = nickname
#    logging.debug("new path: " + '/'.join(histo_path))
#    if args.std_histo_name:
#        #
#        distr = args.std_histo_name
#        histo_name = "{path_part}_{distr}".format(path_part='_'.join(histo_path), distr=distr)
#        logging.debug("new histo name: " + histo_name)
#        out_histo.SetName(histo_name)

for path_tuple, histo in output_histos.items():
    histo_path = '/'.join(path_tuple)
    # check if this directory already exists in the file (a feature for future)
    #out_dir_name = ''.join(part + '/' for part in histo_path)

    if histo_path and fout.Get(histo_path):
        logging.debug('found  ' + histo_path)
        out_dir = fout.Get(histo_path)
    else:
        logging.debug('making ' + histo_path)
        # iteratively create each directory fout -> part0 -> part1 -> ...
        # somehow root did not work for creating them in 1 go
        out_dir = fout
        for directory in histo_path.split('/'):
            logging.debug('making ' + directory)
            nested_dir = out_dir.Get(directory) if out_dir.Get(directory) else out_dir.mkdir(directory)
            nested_dir.cd()
            out_dir = nested_dir

    histo.SetDirectory(out_dir)
    histo.Write()

# save weight counters etc in the top directory in the file
if args.save_weight:
    fout.cd()
    weight_counter.SetName('weight_counter')
    weight_counter.Write()

fout.Write()
fout.Close()

