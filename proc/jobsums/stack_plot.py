import argparse
import logging
from os.path import isfile
from math import sqrt, isnan
from sys import exit

logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "stack histos",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument("mc_file",    help="MC file name")
parser.add_argument("data_file",  help="Data file name")
parser.add_argument("-c", "--channel",  type=str, default='mu_sel', help="selection channels of events to consider (can be a list of channels for shape evolution study)")
parser.add_argument("-s", "--systematic",  type=str, default='NOMINAL', help="systematic variation")
parser.add_argument("-d", "--distr-name",  type=str, default='Mt_lep_met', help='recorded distribution or ratio "distr1/distr2" for fake rates (or sys effect)')
parser.add_argument("-x", "--shape",  type=str, default='', help="selection for shape distributions of wjets and dy")
parser.add_argument("-p", "--plot",  action='store_true', help="don't save the root file, plot stacked histogram in png")
parser.add_argument("-r", "--ratio", action='store_true', help="don't save the root file, make ratio plot (in addition to stack or alone)")
parser.add_argument("-l", "--logy", action='store_true', help="set logarithmic scale for Y axis of stack plot")
parser.add_argument("--normalize", action='store_true', help="normalize the integral to data")
parser.add_argument("-o", "--output-directory", type=str, default='', help="optional output directory")

parser.add_argument("--uncertainty-systematic",  type=str, help="add systematic variations to the hs sum uncertainty")

parser.add_argument("--rebin", type=int, help="rebin the histograms")
parser.add_argument("--infinite-bin-errors", type=float, help="if MC bins get infinite errors (like after rebining) use the sqrt(content)*by this factor (choose from the most relevant MC)")

parser.add_argument("--ratio-range", type=float, default=0.5, help="range of ratio plot (1-range 1+range)")

parser.add_argument("--lumi", type=float, help="use to skip the final normalizing step")

parser.add_argument("--y-max",   type=float, help="set the maximum on Y axis")
parser.add_argument("--x-range", type=str,   help="set the range on X axis")
parser.add_argument("--y-range", type=str,   help="set the range on Y axis `0,1.5`")

#parser.add_argument("--Vloose-shapes",   action='store_true', help="get TT proc shapes from selVloose")
parser.add_argument("--merge-procs", type=str,   help="merge some procs with the format <a_proc>,<another>,<add them to this one>.<and this proc goes into>,<another target> and rename target procs with <old_name:new_name>")

parser.add_argument("--sort-dy",   action='store_true', help="sort DY up")

parser.add_argument("--no-data",      action='store_true', help="don't get data")
parser.add_argument("--no-data-plot", action='store_true', help="don't draw data")

parser.add_argument("--bin-norm",  action='store_true', help="normalize per bin width")

parser.add_argument("--qcd",   type=float, default=0., help="get QCD from corresponding _ss channel and transfer it with the given factor (try --qcd 1)")
parser.add_argument("--qcd-factor",   type=float, default=1., help="factor for MC QCD")
parser.add_argument("--osss",     action='store_true', help="plot the ratio in OS/SS of data-other bkcg")
parser.add_argument("--osss-mc",  action='store_true', help="plot the ratio in OS/SS of MC QCD")
parser.add_argument("--osss-pl",  action='store_true', help="plot the data-other bkcg in the OS and SS selections")
parser.add_argument("--no-label-osss",  action='store_true', help="no anti-iso label")
parser.add_argument("--vert-lines", type=str, help="add vertical lines to the plot at 'x1[,x2]' positions")

parser.add_argument("--nulify-qcd-bins", type=float, help="nulify the negative bins of data - other MC subtraction")

parser.add_argument("--wjets", type=float, default=0., help="scale factor for wjets (from the control region)")
parser.add_argument("--factor-dy", type=float, help="scale factor for dy (from the control region)")
parser.add_argument("--factor-procs", type=str, help="scale factor for given procs p1,p2,p3:factor -- for effect of rate systematics, tau ID and mis-ID")
parser.add_argument("--factor-xsec",  type=str, help="guess x-sec/gen_lumi scale based on usual gen_lumi")

parser.add_argument("--factor-everything", type=float, help="in v37 test13 bunchCONTROLMET everything is processed 9 times for some reason -- factor to scale that out")

parser.add_argument("--factor-rate-systematic", type=str, help="scale the procs by the factor to simulate the rate syst (tau ID, mis-ID): tt_eltau,tt_taultauh,s_top_eltau,s_top_other,dy_tautau:0.95 tt_lj,tt_taulj,tt_other:0.5")
parser.add_argument("--rename-systematic",  type=str, help="save under different name, use for the rate systematics")

parser.add_argument("--draw-overflows", type=float, default=0., help="draw the overflow bin with the specified width")

parser.add_argument("--fake-rate", action='store_true', help='(ad-hocish) the fake rate in MC and data, i.e. ratio of whatever two distributions given as "nominator/denominator"')

parser.add_argument("--cumulative", action='store_true', help="plot stack of cumulative distribution (from right)")
parser.add_argument("--cumulative-fractions", action='store_true', help="cumulative distributions, each bin normalized to 1")

#parser.add_argument("-f", "--form-shapes", type=str, default='usual', help="plot the shapes of distributions normalized to 1")
parser.add_argument("-f", "--form-shapes", action='store_true', help="plot the shapes of distributions normalized to 1")
#parser.add_argument("-e", "--shape-evolution", type=str, default='usual', help="plot the same distr in different channels normalized to 1")
parser.add_argument("--processes", type=str, default='all', help="set processes to consider (all by default)")

parser.add_argument("--data-nick", type=str, default='data', help="data nick")

parser.add_argument("--exp-legend",   action='store_true', help="experimentary legend drawing")
parser.add_argument("--legend-shift", type=float,          help="shift legend on X axis")
parser.add_argument("--skip-legend",  action='store_true', help="don't plot the legend (full view of distribution)")

parser.add_argument("--skip-QCD", action='store_true', help="skip MC QCD")

parser.add_argument("--drop-bin", type=float, help="drop values in the bin")

parser.add_argument("--lumi-label", type=float, help="set lumi label on the plot")
parser.add_argument("--title-x", type=str, default="", help="set title of X axis on the plot")
parser.add_argument("--title-y", type=str, default="", help="set title of Y axis on the plot")
parser.add_argument("--title",   type=str, default="default", help="set title of the plot")

parser.add_argument("--output-name",   type=str, help="name for the output file")
parser.add_argument("--output-suffix", type=str, default='', help="additional suffix in output name")

parser.add_argument("--overwrite",   action='store_true', help="overwrite the output even if it exists")

args = parser.parse_args()

assert isfile(args.mc_file)
assert isfile(args.data_file)

logging.info("import ROOT")

import ROOT
from ROOT import TFile, THStack, TLegend, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TLine
import plotting_root
from draw_overflows import DrawOverflow

#channel = "mu_presel"
#channel = "mu_sel"
#channel = "mu_lj"
#channel = "mu_lj_out"

#sys_name = "NOMINAL"
#distr_name = 'mu_presel_NOMINAL_wjets_Mt_lep_met'

#distr_name = 'Mt_lep_met'
#'met'
#'Mt_lep_met'
##'Mt_lep_met_d'
#'Mt_tau_met'
#'njets'
#'nbjets'
#'dijet_trijet_mass'
#distr_name = 'dijet_trijet_mass'
if args.fake_rate:
    assert '/' in args.distr_name

distr_names = args.distr_name.split('/')

if '/' in distr_names:
    logging.info("getting %d distrs from %s" % (len(distr_names), args.distr_name))

processes_requirement = None
#if args.processes == 'all':
#    processes_requirement = None
#else:
#    processes_requirement = args.processes.split(',')


'''
if args.shape_evolution:
    if args.shape_evolution == 'usual':
        channels = ['wjets', 'tt_mutau', 'tt_eltau', 'tt_lj']
    else:
        channels = args.shape_evolution.split(',')
else:
'''
channels = args.channel.split(',')

logging.info("channels, processes requirement = %s, %s" % (str(channels), processes_requirement))

fdata = TFile(args.data_file)

# only nominal in data
#histo_data = fdata.Get(channel + '/data/' + sys_name + '/' + '_'.join([channel, sys_name, 'data', distr_name]))

# TODO: remove ad-hoc of PU tests
# Mt_lep_met_f_w_mu_trk_b
#histos_data_distrs = [(distr_name, [(histo, 'data', channel)]]
histos_data_distrs    = []
histos_data_distrs_ss = []
histos_data_per_distr    = []
histos_data_per_distr_ss = []

data_nick = args.data_nick
qcd_nick = 'qcd'

if not args.no_data:
  for distr_name in distr_names:
    data_distr_name = distr_name[:]
    if data_distr_name in ('Mt_lep_met_f_w_mu_trk_b', 'Mt_lep_met_f_w_mu_trk_h'):
        data_distr_name = 'Mt_lep_met_f'
    data_distr_name.replace('_w_pu_sum', '_w_pu').replace('_w_pu_b', '_w_pu').replace('_w_pu_h2', '_w_pu').replace('_pu_h2', '').replace('_pu_b', '').replace('_pu_sum', '')
    # get overflows if requested
    if args.draw_overflows:
        #h_overflows = draw_overflows.DrawOverflow(histo, args.draw_overflows)
        histos_data_distrs.append((data_distr_name, [(DrawOverflow(fdata.Get(channel + '/' + data_nick + '/NOMINAL/' + '_'.join([channel, data_nick, 'NOMINAL', data_distr_name])), args.draw_overflows), data_nick, channel) for channel in channels]))
    else:
        histos_data_distrs.append((data_distr_name, [(fdata.Get(channel + '/' + data_nick + '/NOMINAL/' + '_'.join([channel, data_nick, 'NOMINAL', data_distr_name])), data_nick, channel) for channel in channels]))

    if args.qcd > 0. or args.osss or args.osss_mc or args.osss_pl:
        if args.draw_overflows:
            histos_data_distrs_ss.append((data_distr_name, [(DrawOverflow(fdata.Get(channel + '_ss' + '/' + data_nick + '/NOMINAL/' + '_'.join([channel + '_ss', data_nick, 'NOMINAL', data_distr_name])), args.draw_overflows), data_nick, channel) for channel in channels]))
        else:
            histos_data_distrs_ss.append((data_distr_name, [(fdata.Get(channel + '_ss' + '/' + data_nick + '/NOMINAL/' + '_'.join([channel + '_ss', data_nick, 'NOMINAL', data_distr_name])), data_nick, channel) for channel in channels]))

  if args.cumulative or args.cumulative_fractions:
    histos_data_per_distr    = [(name, [(histo.GetCumulative(False), n, c) for histo, n, c in distrs]) for name, distrs in histos_data_distrs]
    histos_data_per_distr_ss = [(name, [(histo.GetCumulative(False), n, c) for histo, n, c in distrs]) for name, distrs in histos_data_distrs_ss]
  else:
    histos_data_per_distr    = histos_data_distrs
    histos_data_per_distr_ss = histos_data_distrs_ss

  if args.factor_everything:
    for _, histos  in histos_data_per_distr + histos_data_per_distr_ss:
      for h, _, _ in histos:
        h.Scale(args.factor_everything)


# normalize data distrs to different bin width:
# or rebin
if not args.no_data and (args.bin_norm or args.rebin):
  for _, distrs in histos_data_per_distr + histos_data_per_distr_ss:
      for histo, _, _ in distrs:
        if args.rebin:
            histo.Rebin(args.rebin)
        if args.bin_norm:
          for bini in range(histo.GetSize()):
              content = histo.GetBinContent(bini)
              error   = histo.GetBinError(bini)
              width   = histo.GetXaxis().GetBinUpEdge(bini) - histo.GetXaxis().GetBinLowEdge(bini)
              histo.SetBinContent(bini, content/width)
              histo.SetBinError(bini, error/width)

logging.info("# data histograms = %d" % len(histos_data_per_distr))

factor_rate_systematic = None
if args.factor_rate_systematic:
    # tt_eltau,tt_taultauh,s_top_eltau,s_top_other,dy_tautau:0.95
    procs, factor = args.factor_rate_systematic.split(':')
    factor_rate_systematic = (procs.split(','), float(factor))

syst_factors = {
              "FragUp":             0.958852, #1.001123,
            "FragDown":             1.02936,  #0.999304,
       #  "SemilepBRUp":             1.008565,
       #"SemilepBRDown":             0.987050,
       #   "PetersonUp":             1.000032,
         "SemilepBRUp":             1.002650,
       "SemilepBRDown":             1.035812,
          "PetersonUp":             0.990578,
# it seems all b-frag related are wrong -- probably I got the weights before the fix there
                "MrUp":             0.896601,
              "MrDown":             1.114154,
                "MfUp":             0.980610,
              "MfDown":             1.025511,
               "MfrUp":             0.874602,
             "MfrDown":             1.135832,
          "AlphaSDown":             1.015573,
            "AlphaSUp":             0.984924,
         "PDFCT14n1Up":             0.999214,
        "PDFCT14n10Up":             0.992985,
        "PDFCT14n11Up":             1.011667,
        "PDFCT14n12Up":             0.991012,
        "PDFCT14n13Up":             1.011832,
        "PDFCT14n14Up":             0.994285,
        "PDFCT14n15Up":             1.020625,
        "PDFCT14n16Up":             0.985324,
        "PDFCT14n17Up":             0.986681,
        "PDFCT14n18Up":             1.015958,
        "PDFCT14n19Up":             0.999024,
         "PDFCT14n2Up":             0.999610,
        "PDFCT14n20Up":             1.001383,
        "PDFCT14n21Up":             1.000811,
        "PDFCT14n22Up":             0.998574,
        "PDFCT14n23Up":             1.004126,
        "PDFCT14n24Up":             0.996106,
        "PDFCT14n25Up":             1.012630,
        "PDFCT14n26Up":             0.987756,
        "PDFCT14n27Up":             1.003978,
        "PDFCT14n28Up":             0.994568,
        "PDFCT14n29Up":             1.002681,
         "PDFCT14n3Up":             0.991550,
        "PDFCT14n30Up":             0.999111,
        "PDFCT14n31Up":             1.001791,
        "PDFCT14n32Up":             0.997694,
        "PDFCT14n33Up":             0.996328,
        "PDFCT14n34Up":             1.008113,
        "PDFCT14n35Up":             0.985943,
        "PDFCT14n36Up":             1.013638,
        "PDFCT14n37Up":             0.993136,
        "PDFCT14n38Up":             1.010604,
        "PDFCT14n39Up":             0.995664,
         "PDFCT14n4Up":             0.995349,
        "PDFCT14n40Up":             1.004863,
        "PDFCT14n41Up":             0.999538,
        "PDFCT14n42Up":             1.000247,
        "PDFCT14n43Up":             0.997833,
        "PDFCT14n44Up":             1.005385,
        "PDFCT14n45Up":             1.004722,
        "PDFCT14n46Up":             1.001524,
        "PDFCT14n47Up":             0.999380,
        "PDFCT14n48Up":             0.993950,
        "PDFCT14n49Up":             0.993164,
         "PDFCT14n5Up":             0.994886,
        "PDFCT14n50Up":             1.003200,
        "PDFCT14n51Up":             1.001966,
        "PDFCT14n52Up":             1.000614,
        "PDFCT14n53Up":             0.968178,
        "PDFCT14n54Up":             1.010604,
        "PDFCT14n55Up":             0.992717,
        "PDFCT14n56Up":             1.012645,
         "PDFCT14n6Up":             1.004233,
         "PDFCT14n7Up":             1.001432,
         "PDFCT14n8Up":             0.992782,
         "PDFCT14n9Up":             1.008439,
}

# the acceptance factors for updowns
#INFO:root:merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_CUETP8M2T4down.root              0.168221 0.175395 0.171805
#INFO:root:merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_CUETP8M2T4up.root                0.168400 0.175170 0.171786
#INFO:root:merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_fsrdown.root                     0.170905 0.177813 0.174358
#INFO:root:merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_fsrup.root                       0.164538 0.171178 0.167857
#INFO:root:merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_hdampDOWN.root                   0.166712 0.173022 0.169864
#INFO:root:merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_hdampUP.root                     0.168795 0.175712 0.172251
#INFO:root:merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_isrdown.root                     0.167539 0.175365 0.171449
#INFO:root:merge-sets/v34/t1/MC2016_Summer16_TTJets_powheg_isrup.root                       0.168896 0.176164 0.172531

#[0.16871733092913158, 0.17559440519640376, 0.17215718299241814]
acceptance_el = 0.168717
acceptance_mu = 0.175594

syst_factors_updowns_el = {
   'TuneCUETP8M2T4Down': acceptance_el / 0.168221,
   'TuneCUETP8M2T4Up':   acceptance_el / 0.168400,
   'FSRDown':            acceptance_el / 0.170905,
   'FSRUp':              acceptance_el / 0.164538,
   'HDAMPDown':          acceptance_el / 0.166712,
   'HDAMPUp':            acceptance_el / 0.168795,
   'ISRDown':            acceptance_el / 0.167539,
   'ISRUp':              acceptance_el / 0.168896,
  }

syst_factors_updowns_mu = {
   'TuneCUETP8M2T4Down': acceptance_mu / 0.175395,
   'TuneCUETP8M2T4Up':   acceptance_mu / 0.175170,
   'FSRDown':            acceptance_mu / 0.177813,
   'FSRUp':              acceptance_mu / 0.171178,
   'HDAMPDown':          acceptance_mu / 0.173022,
   'HDAMPUp':            acceptance_mu / 0.175712,
   'ISRDown':            acceptance_mu / 0.175365,
   'ISRUp':              acceptance_mu / 0.176164,
  }

syst_factors_updowns = {
'tt_eltau':
  syst_factors_updowns_el,
'tt_mutau':
  syst_factors_updowns_mu
}

data_nicks = ('data', 'data_other')

def get_histos(infile, channels, shape_channel, sys_name, distr_name, skip_QCD=args.skip_QCD):
    """get_histos(infile)

    the file contains usual structure
    <channels>/<processes>/<systematics>/<histograms>

    returns [(histo, nick)]
    """
    used_histos = [] # (histo, nick)

    for channel in channels:
        chan = infile.Get(channel)
        if not chan:
            logging.info('no channel %s' % channel)
            continue
        processes_keys = list(chan.GetListOfKeys())
        #sorted_pkeys = sorted(processes_keys, key=lambda pkey: pkey.GetName() if pkey.GetName() not in ('qcd') else 'z_' + pkey.GetName())
        #for process in sorted_pkeys:
        for process in processes_keys:
           nick = process.GetName()
           # pick up the qcd nick from MC
           if 'qcd' in nick:
               qcd_nick = nick

           if skip_QCD and 'qcd' in nick:
               continue
           if processes_requirement and nick not in processes_requirement:
               continue
           #logging.info(nick)

           if nick in data_nicks:
               continue

           # handlinn=g the tt systematics:
           # all other channels get their nominal distr
           # but it must be renamed
           fixed_sys_name = sys_name
           if any(sn in sys_name for sn in ['TOPPT', 'FSR', 'ISR', 'HDAMP', 'TuneCUETP8M2T4', 'QCDbasedCRTune', 'GluonMoveCRTune']):
               fixed_sys_name = sys_name if 'tt' in nick else 'NOMINAL'

           # TODO: PU test hack, remove it
           if nick not in ('wjets', 'dy_other', 'dy_tautau') and distr_name in ('lep_pt_pu_h2', 'lep_pt_pu_b', 'lep_pt_pu_sum',
                'met_pu_h2', 'met_pu_b', 'met_pu_sum',
                'Mt_lep_met_f_pu_h2', 'Mt_lep_met_f_pu_b', 'Mt_lep_met_f_pu_sum',):
               histo_name = '_'.join([channel, nick, fixed_sys_name, distr_name.replace('_pu_h2', '').replace('_pu_b', '').replace('_pu_sum', '')])
           else:
               histo_name = '_'.join([channel, nick, fixed_sys_name, distr_name])
           logging.debug("%s" % histo_name)

           h_init = process.ReadObj().Get(fixed_sys_name + '/' + histo_name)

           if not h_init:
               logging.debug("absent %s" % fixed_sys_name)
               fixed_sys_name = 'NOMINAL'
               histo_name = '_'.join([channel, nick, fixed_sys_name, distr_name])
               h_init = process.ReadObj().Get(fixed_sys_name + '/' + histo_name)

           if not h_init:
               logging.debug("absent NOMINAL %s" % histo_name)
               continue

           # just take shapes for all channels
           if shape_channel: # and nick in ('dy_other', 'dy_tautau', 'wjets'):
               histo_name = '_'.join([shape_channel, nick, fixed_sys_name, distr_name])
               h_shape_path = shape_channel + '/' + nick + '/' + fixed_sys_name + '/' + histo_name
               logging.info("shape from %s" % h_shape_path)

               h_shape = infile.Get(h_shape_path)
               histo = h_shape.Clone()
               if h_shape.Integral() == 0:
                   print "both = 0: %f %f" % (h_init.Integral(), h_shape.Integral())
               else:
                   histo.Scale(h_init.Integral() / h_shape.Integral())
           else:
               histo = h_init

           # rebin if asked
           if args.rebin:
               histo.Rebin(args.rebin)

           # the fast normalization
           if args.lumi:
               tauIDSF_factor = 1.
               #if nick == 'none': #'tt_mutau':
               #    tauIDSF_factor = 0.90
               ##elif nick in ('tt_mutau3ch', 'tt_eltau3ch', 'tt_mutau', 'tt_eltau', 'tt_taultauh', 'dy_tautau', 's_top_eltau', 's_top_mutau', 'dy_tautau'):
               ##    tauIDSF_factor = 0.95
               #elif any(ch in nick for ch in ('ctr_old_mu_sel', 'ctr_old_mu_sel', 'ctr_old_el_sel', 'optel_tight_2L1M', 'optmu_tight_2L1M')):
               #    tauIDSF_factor = 0.95
               #else:
               #    tauIDSF_factor = 1.
               tau_selection = any(ch in channel for ch in ('mu_sel', 'el_sel', 'mu_sel_ss', 'el_sel_ss', 'mu_selSV', 'el_selSV', 'mu_selSV_ss', 'el_selSV_ss', 'dy_mutau', 'ctr_old_mu_sel', 'ctr_old_mu_sel', 'ctr_old_el_sel', 'optel_tight_2L1M', 'optmu_tight_2L1M'))
               tau_process   = nick in ('tt_mutau3ch', 'tt_eltau3ch', 'tt_mutau', 'tt_eltau', 'tt_taultauh', 'dy_tautau', 's_top_eltau', 's_top_mutau')
               if tau_selection and tau_process:
                   tauIDSF_factor = 0.95

               pu_factor = 1.
               #if 'PUUp' in fixed_sys_name:
               #    pu_factor = 1. / (1.01995 if '_el_' in channel or 'el_sel' in channel else 0.9979) # 0.97 # 1./ 0.9979
               #elif 'PUDown' in fixed_sys_name:
               #    pu_factor = 1. / (1.07766 if '_el_' in channel or 'el_sel' in channel else 1.0485) # 1.17 # 1./ 1.485
               #else:
               #    pu_factor = 1. / (1.04678 if '_el_' in channel or 'el_sel' in channel else 1.022)  # 1.06 # 1./ 1.02135 an 1/1.014 with weight counter..
               #    #pu_factor = 1. / (1.04678 if '_el_' in channel else 1.014)
               if 'PUUp' in fixed_sys_name:
                   #pu_factor = 1. / (1.00892 if '_el_' in channel or 'el_sel' in channel else 0.994846) # 0.97 # 1./ 0.9979
                   # el PU is the same as mu, since v37 and Golden json for el
                   pu_factor = 1. / 0.994846
                   # extra norm correction in v37 test13 due to new PU calculation (which does not change NOMINAL norm)
                   # found a miscalculation in PU Up and Down weight -- maybe no need for norm corrections now
                   #pu_factor *= 0.95
                   ## damn! still deviation from old normalizations -- need to return to the old PU, calculated at ntupling
                   #if channel == 'el_sel_lj':
                   #    pu_factor *= 0.972207
                   #elif channel == 'el_sel_ljout':
                   #    pu_factor *= 0.985877
                   #elif channel == 'mu_sel_lj':
                   #    pu_factor *= 0.987258
                   #elif channel == 'mu_sel_ljout':
                   #    pu_factor *= 0.994168
               elif 'PUDown' in fixed_sys_name:
                   #pu_factor = 1. / (1.05718 if '_el_' in channel or 'el_sel' in channel else 1.03657) # 1.17 # 1./ 1.485
                   pu_factor = 1. / 1.03657
                   # extra norm correction in v37 test13
                   #pu_factor *= 1.05
                   ##
                   #if channel == 'el_sel_lj':
                   #    pu_factor *= 1.00701
                   #elif channel == 'el_sel_ljout':
                   #    pu_factor *= 1.00723
                   #elif channel == 'mu_sel_lj':
                   #    pu_factor *= 1.01054
                   #elif channel == 'mu_sel_ljout':
                   #    pu_factor *= 1.01047
               else:
                   #pu_factor = 1. / (1.03106 if '_el_' in channel or 'el_sel' in channel else 1.01388)  # 1.06 # 1./ 1.02135 an 1/1.014 with weight counter..
                   pu_factor = 1. / 1.01388
                   #pu_factor = 1. / (1.04678 if '_el_' in channel else 1.014)

               th_factor = 1.
               if nick[:3] == 'tt_':
                   # watch this tricky bit: Mfr come from v26p1
                   # and something was not completely done there (which was fixed in p2 of v25)
                   # hence there is a factor of difference, hopefuly the shape is not that much affected
                   # 0.958
                   #if 'MrUp' in fixed_sys_name:
                   #    th_factor = 1. / 0.8966
                   #elif 'MrDown' in fixed_sys_name:
                   #    th_factor = 1. / 1.1139
                   #elif 'MfUp' in fixed_sys_name:
                   #    th_factor = 1. / 0.9805
                   #elif 'MfDown' in fixed_sys_name:
                   #    th_factor = 1. / 1.0257
                   #elif 'MfrUp' in fixed_sys_name:
                   #    th_factor = 1. / 0.8747
                   #elif 'MfrDown' in fixed_sys_name:
                   #    th_factor = 1. / 1.1358
                   #elif 'AlphaSUp' in fixed_sys_name:
                   #    th_factor = 1. / 0.98
                   #elif 'AlphaSDown' in fixed_sys_name:
                   #    th_factor = 1. / 1.015
                   #elif 'SemilepBRUp' in fixed_sys_name:
                   #    th_factor = 1. / 1.0087
                   #elif 'SemilepBRDown' in fixed_sys_name:
                   #    th_factor = 1. / 0.987
                   th_factor = 1. / syst_factors.get(fixed_sys_name, 1.)

               ## TODO: remove this, just testing the bSF normalization from simle elmu to close elmu
               #th_factor *= 1. / 0.984960

               # updowns factor:
               if nick[:3] == 'tt_':
                   #th_factor *= syst_factors_updowns.get('nick', {}).get(fixed_sys_name, 1.)
                   if 'el_sel' in channel:
                       th_factor *= syst_factors_updowns_el.get(fixed_sys_name, 1.)
                   else:
                       th_factor *= syst_factors_updowns_mu.get(fixed_sys_name, 1.)

               final_factor = args.lumi * tauIDSF_factor * pu_factor * th_factor
               init_integral = histo.Integral()
               histo.Scale(final_factor)
               logging.info("final factor %20s %20f   %5f %5f %5f     %20f = %10.2f / %10f" % (nick, args.lumi, tauIDSF_factor, pu_factor, th_factor, final_factor, histo.Integral(), init_integral))

           # wjets normalization
           if args.wjets > 0 and 'wjets' in nick:
               histo.Scale(args.wjets)

           # dy normalization
           if args.factor_dy and 'dy' in nick:
               histo.Scale(args.factor_dy)

           if args.factor_procs and nick in args.factor_procs:
               _, factor = args.factor_procs.split(':')
               histo.Scale(float(factor))

           if factor_rate_systematic and nick in factor_rate_systematic[0]:
               histo.Scale(factor_rate_systematic[1])

           if args.factor_everything:
               histo.Scale(args.factor_everything)

           # mc qcd normalization
           if nick in ('qcd', 'qcd_other'):
               histo.Scale(args.qcd_factor)

           # rename the histo for the correct systematic name in case of TT systematics
           histo_name = '_'.join([channel, nick, sys_name, distr_name])
           histo.SetName(histo_name)

           #histo = histo_key.ReadObj()
           logging.info("%s   %s   %x = %f %f,  %f / %f" % (histo.GetName(), histo_name, histo_name == '_'.join([channel, nick, fixed_sys_name, distr_name]), histo.GetEntries(), histo.Integral(), histo.GetBinError(2), histo.GetBinContent(2)))

           # get overflows if requested
           if args.draw_overflows:
               h_overflows = DrawOverflow(histo, args.draw_overflows)
               histo = h_overflows

           # normalize to different bin width:
           if args.bin_norm:
             for bini in range(histo.GetSize()):
                content = histo.GetBinContent(bini)
                error   = histo.GetBinError(bini)
                width   = histo.GetXaxis().GetBinUpEdge(bini) - histo.GetXaxis().GetBinLowEdge(bini)
                histo.SetBinContent(bini, content/width)
                histo.SetBinError(bini, error/width)

           # drop bin
           if args.drop_bin:
               bn_to_drop = histo.FindBin(args.drop_bin)
               histo.SetBinContent (bn_to_drop, 0)
               histo.SetBinError   (bn_to_drop, 0)

           if args.cumulative or args.cumulative_fractions:
               used_histos.append((histo.GetCumulative(False), nick, channel)) # hopefully root wont screw this up
           else:
               used_histos.append((histo, nick, channel)) # hopefully root wont screw this up

    return used_histos


f = TFile(args.mc_file)

def get_histos_with_data_qcd(sys_name):
    used_histos_per_distr = [(distr_name, get_histos(f, channels, args.shape, sys_name, distr_name, skip_QCD = args.skip_QCD)) for distr_name in distr_names]

    # USING NOMINAL QCD:
    # ss MC for QCD is taken at NOMINAL, data is always nominal
    # then SS MC is scaled
    # and substituted in OS histos

    # for ss channels I need data - MC sum
    # and I substitute QCD with this difference everywhere
    hs_sums2_ss = [] # these will be sums of no-QCD MC in SS region

    # sum all MC except qcd
    def datadriven_qcd(channel, distr, sys):
        logging.warning("datadriven qcd")
        # get SS data and MC
        used_histos_per_distr_ss = get_histos(f, [channel + '_ss'], args.shape, sys, distr)
        data_hist = fdata.Get(channel + '_ss' + '/' + data_nick + '/NOMINAL/' + '_'.join([channel + '_ss', data_nick, 'NOMINAL', distr]))

        if args.draw_overflows:
            h_overflows = DrawOverflow(data_hist, args.draw_overflows)
            data_hist = h_overflows

        # find QCD shape in SS
        mc_hist = used_histos_per_distr_ss[0][0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
        mc_hist.SetName('mc_sum2_ss')
        for h, nick, channel in used_histos_per_distr_ss[1:]:
            if 'qcd' in nick: continue # don't include QCD MC -- its sum of no-QCD MC
            mc_hist.Add(h)

        logging.warning("qcd hist  %d  %d" % (data_hist.GetEntries(), mc_hist.GetEntries()))
        qcd_hist = data_hist - mc_hist
        # nulify negative bins
        for bini in range(qcd_hist.GetSize()):
            if qcd_hist.GetBinContent(bini) < 0:
                qcd_hist.SetBinContent(bini, 0)

        # scale by the given factor and return
        qcd_hist.Scale(args.qcd)

        return qcd_hist


    if args.qcd > 0. or args.osss or args.osss_mc or args.osss_pl:
        # get all the same distributions for _ss channel
        used_histos_per_distr_ss = [(distr_name, get_histos(f, [c + '_ss' for c in channels], args.shape, sys_name, distr_name)) for distr_name in distr_names]
        #used_histos_per_distr_ss = [(distr_name, get_histos(f, [c + '_ss' for c in channels], args.shape, 'NOMINAL', distr_name)) for distr_name in distr_names]

        for distr, histos in used_histos_per_distr_ss:
            # loop through normal list
            hs_sum2 = histos[0][0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
            hs_sum2.SetName('mc_sum2_ss')

            for h, nick, channel in histos[1:]: # I sum different channels?... the old hacked not fixed yet
                if 'qcd' in nick: continue # don't include QCD MC -- its sum of no-QCD MC
                hs_sum2.Add(h)
            channel = histos[0][2]
            logging.debug("all channels are the same: %r" % all(ch == channel for _, _, ch in histos))

            hs_sum2.SetFillStyle(3004);
            hs_sum2.SetFillColor(1);
            hs_sum2.SetMarkerStyle(1)
            hs_sum2.SetMarkerColor(0)

            # TODO: it turns messy, need to fix the channels
            hs_sums2_ss.append((distr, [(hs_sum2, "mc_sum", channel)]))

    # subtract Data - MC in SS, get the qcd distr
    if args.qcd > 0.:
        # find difference and substitute QCD
        # histos_data_per_distr_ss = (distr, hists = [(histo, nick=data, channel)...])
        # I need the same structure for the sums
        qcd_hists = {} # distr, channel: qcd hist
        for (distr, mc_sums), (_, datas) in zip(hs_sums2_ss, histos_data_per_distr_ss):
            for (mc_hist, _, channel), (data_hist, _, _) in zip(mc_sums, datas):
                logging.debug('calulating qcd  %d  %d' % (data_hist.GetSize(), mc_hist.GetSize()))

                qcd_hist = data_hist - mc_hist
                logging.debug('data qcd got in channel, Integral = %f, difference %f' % (qcd_hist.Integral(), data_hist.Integral() - mc_hist.Integral()))

                # negative qcd bins are equalized to zero:
                '''
                for (Int_t i=0; i<=histo->GetSize(); i++)
                        {
                        //yAxis->GetBinLowEdge(3)
                        double content = histo->GetBinContent(i);
                        double error   = histo->GetBinError(i);
                        double width   = histo->GetXaxis()->GetBinUpEdge(i) - histo->GetXaxis()->GetBinLowEdge(i);
                        histo->SetBinContent(i, content/width);
                        histo->SetBinError(i, error/width);
                        }
                '''

                # here I found nominal "yield" datadriven qcd
                # in case a shape channel is given -- find the qcd there and normalize to this yield
                if args.shape:
                    #shape_qcd = datadriven_qcd(args.shape, distr, 'NOMINAL') # nominal only for now
                    shape_qcd = datadriven_qcd(args.shape, distr, sys_name)
                    shape_qcd.Scale(qcd_hist.Integral() / shape_qcd.Integral())
                    qcd_hist = shape_qcd

                    logging.debug('data qcd got shape from %s, Integral = %f' % (args.shape, qcd_hist.Integral()))

                for bini in range(qcd_hist.GetSize()):
                    if qcd_hist.GetBinContent(bini) < 0:
                        qcd_hist.SetBinContent(bini, 0)

                #datadriven_qcd_name = "qcd"
                # optmu_tight_2L1M_wjets_NOMINAL_Mt_lep_met
                # I remove _ss from channel name
                # qcd_nick is defined from the MC qcd name or it is 'qcd'
                datadriven_qcd_name = '%s_%s_%s_%s' % (channel[:-3], 'qcd', args.systematic, distr)
                qcd_hist.SetName(datadriven_qcd_name)
                qcd_hist.Scale(args.qcd)
                qcd_hists[(distr, channel)] = qcd_hist
                logging.debug('data qcd for (%s, %s) Integral = %f' % (distr, channel, qcd_hist.Integral()))



    hs_sums2 = []

    hs_sums_noqcd = []

    # calculate MC sum and substitute data driven qcd
    for distr, histos in used_histos_per_distr:
        # loop through normal list
        hs_sum2 = histos[0][0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
        hs_sum2.SetName('mc_sum2')
        hs_sum2.Reset()

        hs_sum_noqcd = histos[0][0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
        hs_sum_noqcd.SetName('mc_sum_noqcd')
        hs_sum_noqcd.Reset()

        # substitute QCD
        for i, (h, nick, channel) in enumerate(histos[0:]):
            if args.qcd > 0. and 'qcd' in nick:
                # then take the qcd from the differences in ss region
                h = qcd_hists[(distr, channel + '_ss')]
                # and substitute it in the list (if it won't brak everything)
                histos[i] = (h, nick, channel) # h is new here!
            hs_sum2.Add(h)
            if 'qcd' not in nick:
                hs_sum_noqcd.Add(h)

        hs_sum2.SetFillStyle(3004);
        hs_sum2.SetFillColor(1);
        hs_sum2.SetMarkerStyle(1)
        hs_sum2.SetMarkerColor(0)

        hs_sums2.append(hs_sum2)
        hs_sums_noqcd.append(hs_sum_noqcd)

    return used_histos_per_distr, hs_sums2, hs_sums_noqcd, hs_sums2_ss


sys_name = args.systematic
used_histos_per_distr, hs_sums2, hs_sums_noqcd, hs_sums2_ss = get_histos_with_data_qcd(sys_name)

# here all processes, including data-driven qcd, are finally calculated

# for systematic uncertainties the sums of MC are needed at different systematics, with qcd recalculated
# I need to get MC sums and propagate the uncertainty to the one which goes into the plot
# the sum on the plot is taken from hs_sum1.Draw("same e2") from hs_sum1 = hs_sums1[0]
# which is summed from hs stacks
# but it actually equals to hs_sum2 = hs_sums2[0]

if args.uncertainty_systematic:
    logging.debug("adding systematic uncertainties %s" % args.uncertainty_systematic)
    nominal_sum = hs_sums2[0]
    for sys_var in args.uncertainty_systematic.split(','):
        _, hs_sums2_Up,   _, _ = get_histos_with_data_qcd(sys_var + 'Up')
        _, hs_sums2_Down, _, _ = get_histos_with_data_qcd(sys_var + 'Down')
        hs_sum2_Up   = hs_sums2_Up[0]
        hs_sum2_Down = hs_sums2_Down[0]
        # for each bin add in quadrature the deviation from nominal
        for bini in range(nominal_sum.GetSize()):
            err = nominal_sum.GetBinError(bini)
            dev_up   = abs(nominal_sum.GetBinContent(bini) - hs_sum2_Up.GetBinContent(bini))
            dev_down = abs(nominal_sum.GetBinContent(bini) - hs_sum2_Down.GetBinContent(bini))
            #dev = 0.5 * (dev_up + dev_down)
            dev = max(dev_up, dev_down)
            nominal_sum.SetBinError(bini, sqrt(err**2 + dev**2))


# merge some of the processes if asked
#used_histos_per_distr = [(distr_name, get_histos(f, channels, args.shape, sys_name, distr_name, skip_QCD = args.skip_QCD)) for distr_name in distr_names]
# [(distr_name, [(histo, nick, channel)])]
if args.merge_procs:
    procs_to_merge = {} # proc: target proc
    new_names = {}
    for procs_bunch in [bunch.split(',') for bunch in args.merge_procs.split('.') if ',' in bunch]:
        merged_procs, target_proc = procs_bunch[:-1], procs_bunch[-1]
        if ':' in target_proc:
            target_proc, new_name = target_proc.split(':')
            new_names[target_proc] = new_name
        logging.debug('merging %s to %s' % (repr(merged_procs), target_proc))
        for proc in merged_procs:
            procs_to_merge[proc] = target_proc

    merged_histos = []
    for distr, histos in used_histos_per_distr:
        histos_per_nick = dict([(n,h) for h,n,_ in histos])
        channel = histos[0][2]
        # FIXME it will break if there are multiple channels
        # add up the merged procs to the target
        # and pop out the merged procs
        for merged_proc, target in procs_to_merge.items():
            histos_per_nick[target].Add(histos_per_nick[merged_proc])
        for merged_proc in procs_to_merge.keys():
            histos_per_nick.pop(merged_proc)

        # the remaining histos are the targets or the untouched processes
        # rename them with the new names
        for target_nick, histo in histos_per_nick.items():
            if target_nick not in new_names: continue
            new_nick = new_names[target_nick]
            old_name = histos_per_nick[target_nick].GetName()
            new_name = old_name.replace(target_nick, new_nick)
            histos_per_nick[target_nick].SetName(new_name)

        merged_histos.append((distr, [(h, new_names.get(n, n), channel) for n,h in histos_per_nick.items()]))

    used_histos_per_distr = merged_histos

hs_sum2 = hs_sums2[0]
if args.qcd > 0. or args.osss or args.osss_mc or args.osss_pl:
    hs_sum_noqcd    = hs_sums_noqcd[0]
    #hs_sum_noqcd_ss = hs_sums2_ss[0]
    #hs_sums2_ss.append((distr, [(hs_sum2, "mc_sum", channel)]))
    hs_sum_noqcd_ss = hs_sums2_ss[0][1][0][0]


# normalize MC processes and data to MC sum bin-by-bin
if args.cumulative_fractions:
    for (name, histos), mc_sum in zip(used_histos_per_distr, hs_sums2):
        for bin_n in range(mc_sum.GetXaxis().GetNbins() + 1):
            mc_sum_bin = mc_sum.GetBinContent(bin_n)
            if mc_sum_bin != 0:
                # all the mc processes
                for h, nick, channel in histos:
                    mc_bin = h.GetBinContent(bin_n)
                    ratio = mc_bin/mc_sum_bin
                    h.SetBinContent(bin_n, ratio)
                    h.SetBinError(bin_n, h.GetBinError(bin_n)/mc_sum_bin)
                # data
                for distr_name, histos_per_channel in histos_data_per_distr:
                    for h, _, _ in histos_per_channel:
                        data_bin = h.GetBinContent(bin_n)
                        ratio = data_bin/mc_sum_bin
                        h.SetBinContent(bin_n, ratio)
                        h.SetBinError(bin_n, h.GetBinError(bin_n)/mc_sum_bin)
            # the sum
            mc_sum.SetBinContent(bin_n, 1.)
            if mc_sum_bin != 0:
                mc_sum.SetBinError(bin_n, mc_sum.GetBinError(bin_n)/mc_sum_bin)

# ??
histos_data_sums_per_distr = []
for distr_name, histos_per_channel in histos_data_per_distr:
    histos_data_sum = histos_per_channel[0][0].Clone()
    histos_data_sum.Sumw2(ROOT.kTRUE) # are errors saved here?
    for h, _, _ in histos_per_channel[1:]:
        histos_data_sum.Add(h)
    histos_data_sums_per_distr.append(histos_data_sum)

histos_data_sum = histos_data_sums_per_distr[0] # [:1]

if args.normalize:
    ratio = histos_data_sum.Integral() / hs_sum2.Integral()
    hs_sum2.Scale(ratio)
    for name, histos in used_histos_per_distr:
        for h, nick, channel in histos:
            h.Scale(ratio)

'''
done more generally now

if args.form_shapes:
    if args.form_shapes == 'usual':
        shape_nicks = ['wjets', 'tt_mutau', 'tt_eltau', 'tt_lj']
    else:
        shape_nicks = args.form_shapes.split(',')
    used_histos = [i for i in used_histos if i[1] in shape_nicks]
'''

# set data histo styles and
# add data to the legend
#for histo_data in histos_data:
#    histo_data.SetMarkerStyle(21)
#    leg.AddEntry(histo_data, "data", "e1 p")
# no, add just the data_sum entry
if args.exp_legend:
    shift = args.legend_shift if args.legend_shift else 0.
    leg = TLegend(0.8 - shift, 0.55, 1.   - shift, 0.92)
else:
    shift = args.legend_shift if args.legend_shift else 0.
    leg = TLegend(0.7 - shift, 0.55, 0.89 - shift, 0.92)

# data is first in the legend
if not args.fake_rate and not args.skip_legend and not (args.no_data or args.no_data_plot):
    histos_data_sum.SetMarkerStyle(21)
    leg.AddEntry(histos_data_sum, "data", "lep")

# add the legend entry for MC sum error band
leg.AddEntry(hs_sum2, "MC sum", 'f')

# get MC stack and legend for it
#hs, leg = plotting_root.stack_n_legend(used_histos)
if args.legend_shift:
    hs_legs_per_distr = [(distr_name, plotting_root.stack_n_legend(histos, args.legend_shift, args.exp_legend, sort_dy=args.sort_dy, leg=leg)) for distr_name, histos in used_histos_per_distr]
else:
    hs_legs_per_distr = [(distr_name, plotting_root.stack_n_legend(histos, exp_legend=args.exp_legend, sort_dy=args.sort_dy, leg=leg)) for distr_name, histos in used_histos_per_distr]

# sum of MC to get the sum of errors

# loop through TList
hs_sums1 = []

for distr_name, (hs, _) in hs_legs_per_distr:
    logging.info("summing hs_sum1 for %s" % distr_name)
    histos = hs.GetHists() # TList
    hs_sum1 = histos[0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
    hs_sum1.SetName('mc_sum1')

    for h in histos[1:]:
        hs_sum1.Add(h)

    hs_sum1.SetFillStyle(3004);
    hs_sum1.SetFillColor(1);
    hs_sum1.SetMarkerStyle(1)
    hs_sum1.SetMarkerColor(0)

    hs_sums1.append(hs_sum1)

'''
Plot ratio of first two distrs for fake rate.
Otherwise plot just the 0th distr.
TODO: somehow improve/simplify/make useful this system.
'''

logging.info("len hs_sums1 %d" % len(hs_sums1))
hs_sum1 = hs_sums1[0]
if len(histos_data_sums_per_distr) > 1:
    histos_data_sum2 = histos_data_sums_per_distr[1]
    hs_sum1_2 = hs_sums1[1]

hs, leg =  hs_legs_per_distr[0][1]
used_histos = used_histos_per_distr[0][1]
#histos_data_per_distr[0][0][0].Print()

logging.info("data       = %15f %15f" % (0. if args.no_data else histos_data_sum.Integral(), sum(histos_data_sum.GetBinError(bini) for bini in range(histos_data_sum.GetSize()))))
logging.info("mc sum     = %15f %15f" % (hs_sum1.Integral(), hs_sum2.Integral()))

logging.info("mc sum unc = %15f %15f" % (sum(hs_sum1.GetBinError(bini) for bini in range(hs_sum1.GetSize())), sum(hs_sum2.GetBinError(bini) for bini in range(hs_sum2.GetSize()))))

if args.infinite_bin_errors:
  for bini in range(hs_sum2.GetSize()):
    if isnan(hs_sum2.GetBinError(bini)):
        content = hs_sum2.GetBinContent(bini)
        if content < 0:
            hs_sum2.SetBinContent(bini, 0)
            hs_sum2.SetBinError(bini, 0)
        else:
            hs_sum2.SetBinError(bini, 3*sqrt(hs_sum2.GetBinContent(bini)))

logging.info("mc sum unc = %15f %15f" % (sum(hs_sum1.GetBinError(bini) for bini in range(hs_sum1.GetSize())), sum(hs_sum2.GetBinError(bini) for bini in range(hs_sum2.GetSize()))))

out_dir = args.output_directory + '/' if args.output_directory else './'

if not args.no_data:
    histos_data = histos_data_per_distr[0][1]

    histos_data[0][0].GetXaxis().SetLabelFont(63)
    histos_data[0][0].GetXaxis().SetLabelSize(14) # labels will be 14 pixels
    histos_data[0][0].GetYaxis().SetLabelFont(63)
    histos_data[0][0].GetYaxis().SetLabelSize(14) # labels will be 14 pixels


if not args.plot and not args.ratio and not args.osss_pl:
    '''
    For given channel, systematic and distribution
    save the selected distributions, a THStack of them,
    their sum, data distribution and legend.

    Save as usual in TDirs:
    channel/
         proc/
            sys/
                distrs

    also save MC sums:
    cannel/
        sums/
	    sys/
	        distrs
    '''

    options = ""
    options += "_data-qcd" if args.qcd > 0. else ""
    options += ("_shape-%s" % args.shape) if args.shape else ""

    if args.rename_systematic:
        sys_name = args.rename_systematic

    filename = out_dir + args.mc_file.split('.root')[0].replace('/', ',') + '_HISTOSEL_%s_%s_%s%s.root' % (distr_name, channel, sys_name, options)
    logging.info('saving root %s' % filename)

    if not args.overwrite and isfile(filename):
        print filename, "exists, nothing written"
        exit(0)

    fout = TFile(filename, 'RECREATE')
    #fout = TFile(filename, 'CREATE')
    fout.cd()

    # create the channel/sys/ TDir
    #out_dir = fout.mkdir('%s/%s/%s' % (distr_name, channel, sys_name))
    # this doesn't work
    #out_dir = fout.mkdir(distr_name)
    #out_dir.cd()
    chan_dir = fout.mkdir(channel)
    chan_dir.cd()
    #out_dir = out_dir.mkdir(sys_name)
    #out_dir.cd()

    #histo_data.Write() #'data')
    for h, nick, _ in histos_data + used_histos:
        logging.info(nick)
        proc_dir = chan_dir.mkdir(nick)
        proc_dir.cd()
        syst_dir = proc_dir.mkdir(sys_name)
        syst_dir.cd()

        h.SetDirectory(syst_dir)
        h.Write()

    chan_dir.cd()
    #syst_dir = proc_dir.mkdir(sys_name) # object of this name already exists...
    # now this is why Unix is such a great system: it does have hierarchical structures, global vs local scope, modularity
    # ROOT was made in 90th and it does not have these features
    sums_dir = chan_dir.mkdir("sums_" + sys_name)
    sums_dir.cd()
    #syst_dir = proc_dir.mkdir(sys_name)
    #syst_dir.cd()
    #hs      .SetName (hs     .GetName() + '_' + sys_name)
    hs_sum1 .SetName (hs_sum1.GetName() + '_' + sys_name)
    hs_sum2 .SetName (hs_sum2.GetName() + '_' + sys_name)
    #hs      .SetDirectory (sums_dir) # 'THStack' object has no attribute 'SetDirectory'
    # least expected nonuniformity of interface
    hs_sum1 .SetDirectory (sums_dir)
    hs_sum2 .SetDirectory (sums_dir)
    #hs      .Write()
    hs_sum1 .Write()
    hs_sum2 .Write()

    ## AttributeError: 'THStack' object has no attribute 'SetDirectory'
    ## therefore no THStack in the output
    #for stuff in [hs_sum1, hs_sum2]:
    #    stuff.SetDirectory(out_dir)
    #    stuff.Write()
    ## also
    ## AttributeError: 'TLegend' object has no attribute 'SetDirectory'
    ##leg.SetDirectory(out_dir)
    ##leg.Write('leg')

    fout.Write()
    fout.Close()

    #from ROOT import gROOT
    #gROOT.Reset()
    logging.info("no segfault here")

elif args.form_shapes:
    from ROOT import gStyle, gROOT, TCanvas, TPad
    gROOT.SetBatch()
    cst = TCanvas("cst","stacked hists",10,10,700,700)
    gStyle.SetOptStat(0)
    pad = TPad("pad","This is pad", 0., 0.,  1., 1.)
    pad.Draw()

    # normalize histograms to 1
    for histo_data, _, _ in histos_data:
        histo_data.Sumw2(ROOT.kTRUE) # to correctly save the errors after scaling
        histo_data.Scale(1/histo_data.Integral()) # but errors are not scaled here?

    for histo, nick, _ in used_histos:
        #if nick not in args.processes:
            #continue
        histo.Scale(1/histo.Integral())

    gStyle.SetHatchesSpacing(2)

    # plot stuff
    pad.cd()

    # find max Y for plot
    max_y = max([h.GetMaximum() for h, _, _ in used_histos + histos_data])
    if args.y_max:
        max_y = args.y_max

    if args.y_range:
        min_y, max_y = [float(x) for x in args.y_range.split(',')]
        histos_data[0][0].SetMinimum(min_y)

    if not args.no_data:
        histos_data[0][0].SetMaximum(max_y * 1.1)
        histos_data[0][0].SetXTitle(distr_name)
        histos_data[0][0].Draw('e1 p')

    histos_loop = [h_record for h_record in histos_data[1:] + used_histos if h_record[1] in args.processes]
    for i, (histo, nick, _) in enumerate(histos_loop):
        if nick not in args.processes:
            continue
        #histo.SetFillColor( # it shouldn't be needed with hist drawing option
        # nope, it's needed...
        histo.SetLineColor(plotting_root.nick_colour[nick])
        histo.SetLineWidth(4)
        histo.SetMaximum(max_y * 1.1)
        histo.SetXTitle(distr_name)

        if i < 2:
            histo.SetFillStyle([3354, 3345][i])
        else:
            histo.SetFillColorAlpha(0, 0.0)
        print nick
        histo.Draw("hist" if args.no_data and i<=len(histos_data[1:]) else "e same")

    # draw last histogram with markers to overlay over the hist
    last_histo, _, _ = histos_loop[-1]
    last_histo.Draw("e same")
    if not args.no_data:
        histos_data[0][0].Draw('same e1 p')

    if not args.skip_legend:
        leg.Draw("same")

    data_driv_qcd = "_data-qcd" if args.qcd > 0. else ""

    #cst.SaveAs(out_dir + "shapes%s%s%s.png" %('_' + args.mc_file.split('.')[0], '_' + args.processes, '_data-qcd' if args.qcd > 0. else ''))
    fname = '_'.join((args.mc_file.replace('/', ',').split('.root')[0], args.data_file.replace('/', ',').split('.root')[0], distr_name, channel, sys_name))
    cst.SaveAs(out_dir + fname + '_shapes_%d-channels_%s-processes%s.png' % (len(channels), '-'.join(args.processes.split(',')), data_driv_qcd))

elif args.osss or args.osss_mc or args.osss_pl:
    from ROOT import gStyle, gROOT, TCanvas, TPad
    gROOT.SetBatch()
    gStyle.SetOptStat(0)
    cst = TCanvas("cst","stacked hists",10,10,700,700)

    # prepare output target
    outname = out_dir + '_'.join((args.mc_file.replace('/', ',').split('.root')[0], args.data_file.replace('/', ',').split('.root')[0], distr_name, channel, sys_name))
    if args.osss:
        outname += '_osss-factor'
    elif args.osss_mc:
        outname += '_osss-mc-factor'
    elif args.osss_pl:
        outname += '_osss-parts'

    if args.logy:
        pad1.SetLogy()
        outname += '_logy'

    # the no-qcd MC sum is in
    # hs_sum_noqcd and hs_sum_noqcd_ss
    # the data is in
    # histos_data_per_distr and histos_data_per_distr_ss
    #histos_data_distrs = [(distr_name, [(histo, 'data', channel)]]

    if args.osss or args.osss_pl:
        histo_diff_os = histos_data_per_distr   [0][1][0][0] - hs_sum_noqcd
        histo_diff_ss = histos_data_per_distr_ss[0][1][0][0] - hs_sum_noqcd_ss
    elif args.osss_mc:
        #used_histos_per_distr = [(distr_name, get_histos(f, channels, args.shape, sys_name, distr_name, skip_QCD = args.skip_QCD)) for distr_name in distr_names]
        # [(distr_name, [h, nick, channel in histos])]
        histo_diff_os = [h for h, nick, _ in used_histos_per_distr   [0][1] if nick == 'qcd'][0]
        histo_diff_ss = [h for h, nick, _ in used_histos_per_distr_ss[0][1] if nick == 'qcd'][0]

    # nulify negative bins
    if args.nulify_qcd_bins is not None:
        for histo in (histo_diff_os, histo_diff_ss):
          for bini in range(histo.GetSize()):
            if histo.GetBinContent(bini) < 0.:
                #print histo.GetBinContent(bini)
                histo.SetBinContent(bini, args.nulify_qcd_bins)
            if isnan(histo.GetBinError(bini)):
                #print histo.GetBinContent(bini)
                histo.SetBinError(bini, sqrt(abs(histo.GetBinContent(bini))))
    # for some reason the plot still looks like there are negative bins

    # QCD shape is obtained, either MC or Data

    #for histo in (histo_diff_os, histo_diff_ss):
    #  for bini in range(histo.GetSize()):
    #        print histo.GetBinContent(bini)

    if not args.plot and not args.ratio:
        # save the distribution and exit
        fout = TFile(outname + '.root', 'RECREATE')
        fout.cd()

        histo_diff_os.SetName("qcd_os")
        histo_diff_ss.SetName("qcd_ss")

        histo_diff_os.Write()
        histo_diff_ss.Write()

        fout.Write()
        fout.Close()
        exit(0)

    if not args.osss_pl:
        histo_diff_os.Divide(histo_diff_ss)

    pad1 = TPad("pad1","This is pad1", 0., 0.,  1., 1.)
    pad1.Draw()

    pad1.cd()

    if args.title == 'default':
        title_plot = "%s %s" % (channel, sys_name)
    else:
        title_plot = args.title

    title_x = args.title_x if args.title_x else args.distr_name
    title_y = args.title_y if args.title_y else "OS/SS"

    min_y, max_y = 0., max(histo_diff_os.GetMaximum(), histo_diff_ss.GetMaximum())
    if args.y_range:
        min_y, max_y = [float(x) for x in args.y_range.split(',')]
        histo_diff_os.SetMaximum(max_y)
        histo_diff_os.SetMinimum(min_y)

    # tune the typography: font sizes, title offsets etc
    histo_diff_os.GetXaxis().SetTitleFont(63)
    histo_diff_os.GetXaxis().SetTitleSize(20)
    histo_diff_os.GetYaxis().SetTitleFont(63)
    histo_diff_os.GetYaxis().SetTitleSize(20)
    histo_diff_os.GetYaxis().SetTitleOffset(1.4)

    # set titles and draw
    histo_diff_os.SetXTitle(title_x)
    histo_diff_os.SetYTitle(title_y)
    histo_diff_os.SetTitle(title_plot)
    histo_diff_os.Draw()

    if args.osss_pl:
        histo_diff_ss.SetLineColor(kRed)
        histo_diff_ss.Draw("same")

        # add the legend
        leg = TLegend(0.6, 0.7, 0.9, 0.89)
        leg.SetBorderSize(0)
        leg.AddEntry(histo_diff_os, "opposite sign", "l")
        leg.AddEntry(histo_diff_ss, "same sign", "l")
        leg.Draw("same")

    if args.vert_lines:
        x_positions = [float(x) for x in args.vert_lines.split(',')]
        for x in x_positions:
            l = TLine(x, min_y, x, max_y)
            l.SetLineStyle(7)
            l.Draw("same")

    #left_title = TPaveText(0.1, 0.9, 0.4, 0.94, "brNDC")
    #left_title.AddText("CMS preliminary at 13 TeV")
    left_title = TPaveText(0.12, 0.82, 0.2, 0.89, "brNDC")
    left_title.AddText("CMS")
    left_title.SetTextFont(1)
    left_title.SetFillColor(0)

    #right_title = TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC")
    #right_title.AddText("L = %s fb^{-1}" % (args.lumi / 1000. if args.lumi else args.lumi_label))
    right_title = TPaveText(0.65, 0.9, 0.9, 0.95, "brNDC")
    right_title.AddText("%s fb^{-1} (13 TeV)" % (args.lumi_label / 1000. if args.lumi_label else args.lumi / 1000.))
    right_title.SetTextFont(132)
    right_title.SetFillColor(0)

    left_title .Draw("same")
    right_title.Draw("same")

    if not args.no_label_osss:
        label = TPaveText(0.5, 0.2, 0.8, 0.3, "brNDC")
        label.AddText("anti-iso region")
        label.SetTextFont(132)
        label.SetFillColor(0)
        label.Draw("same")

    cst.SaveAs(outname + '.png')

else:
    '''
    plot stack plot and/or ratio plot
    '''

    from ROOT import gStyle, gROOT, TCanvas, TPad
    gROOT.SetBatch()
    if args.exp_legend:
        cst = TCanvas("cst","stacked hists",10,10,1000,700)
    else:
        cst = TCanvas("cst","stacked hists",10,10,700,700)
    gStyle.SetOptStat(0)

    if args.title == 'default':
        title_plot = "%s %s" % (channel, sys_name)
    else:
        title_plot = args.title

    title_x = args.title_x if args.title_x else args.distr_name
    title_y = args.title_y if args.title_y else "events/bin"

    pad_right_edge = 0.8 if args.exp_legend else 1.

    if args.ratio and args.plot:
        pad1 = TPad("pad1","This is pad1", 0., 0.3,   pad_right_edge, 1.)
        pad2 = TPad("pad2","This is pad2", 0., 0.05,  pad_right_edge, 0.3)

        # trying to draw the axis labels on top of histograms
        pad1.GetFrame().SetFillColor(42)
        pad1.GetFrame().SetBorderMode(1)
        pad1.GetFrame().SetBorderSize(5)

        pad2.GetFrame().SetFillColor(42)
        pad2.GetFrame().SetBorderMode(1)
        pad2.GetFrame().SetBorderSize(5)

        #pad2.SetTopMargin(0.01) # doesn't work
        #gStyle.SetPadTopMargin(0.05) # nope
        #ROOT.gPad.SetTopMargin(0.01) # nope
        if args.logy:
            logging.info("setting Pad logy")
            pad1.SetLogy()
        pad1.Draw()
        pad2.Draw() # these have to be before cd()
        # now set margins:
        pad1.cd()
        ROOT.gPad.SetBottomMargin(0.02)
        if args.exp_legend:
            ROOT.gPad.SetRightMargin(0.02)
        #ROOT.gPad.SetTopMargin(0.01)
        pad2.cd()
        #ROOT.gPad.SetBottomMargin(0.001)
        ROOT.gPad.SetTopMargin(0.02)
        if args.exp_legend:
            ROOT.gPad.SetRightMargin(0.02)

    elif args.ratio:
        pad2 = TPad("pad2","This is pad2", 0., 0.05,  pad_right_edge, 1.)
        pad2.Draw()
        pad2.cd()
        if args.exp_legend:
            ROOT.gPad.SetRightMargin(0.02)
    else:
        pad1 = TPad("pad1","This is pad1", 0., 0.,  pad_right_edge, 1.)
        pad1.Draw()
        pad1.cd()
        if args.exp_legend:
            ROOT.gPad.SetRightMargin(0.02)

    # if fake rate then the sums should be divided by second distr sum
    if args.fake_rate:
        histos_data_sum.Divide(histos_data_sum2)
        hs_sum2.Divide(hs_sum1_2)

    if args.x_range:
        x_min, x_max = (float(x) for x in args.x_range.split(','))
        histos_data_sum .SetAxisRange(x_min, x_max, "X")
        hs_sum2         .SetAxisRange(x_min, x_max, "X")

    # drop bin
    if args.drop_bin:
        bn_to_drop = histos_data_sum.FindBin(args.drop_bin)
        histos_data_sum.SetBinContent (bn_to_drop, 0)
        histos_data_sum.SetBinError   (bn_to_drop, 0)

    print "MC sum error", hs_sum2.GetBinError(2), hs_sum2.GetBinContent(2)

    # plotting
    if args.ratio:
        pad2.cd()

        ratio_max = 1. + args.ratio_range
        ratio_min = 1. - args.ratio_range
        # calculate ratios
        if not (args.no_data or args.no_data_plot):
            histo_data_relative = histos_data_sum.Clone()
            histo_data_relative.SetName("rel_data")
            histo_data_relative.SetStats(False)
            histo_data_relative.SetMaximum(ratio_max)
            histo_data_relative.SetMinimum(ratio_min)
            histo_data_relative.Divide(hs_sum2)

            histo_data_relative.GetXaxis().SetLabelFont(63)
            histo_data_relative.GetXaxis().SetLabelSize(14) # labels will be 14 pixels
            histo_data_relative.GetYaxis().SetLabelFont(63)
            histo_data_relative.GetYaxis().SetLabelSize(14) # labels will be 14 pixels
            histo_data_relative.GetXaxis().SetTitleFont(63)
            histo_data_relative.GetXaxis().SetTitleSize(20)
            histo_data_relative.GetYaxis().SetTitleFont(63)
            histo_data_relative.GetYaxis().SetTitleSize(20)

        hs_sum1_relative = hs_sum2.Clone()
        hs_sum1_relative.SetName("rel_mc")
        hs_sum1_relative.SetStats(False)

        #histo_data_relative.GetYaxis().SetRange(0.5, 1.5)
        #histo_data_relative.GetYaxis().SetUserRange(0.5, 1.5)
        #hs_sum1_relative.GetYaxis().SetRange(0.5, 1.5)
        #hs_sum1_relative.GetYaxis().SetUserRange(0.5, 1.5)

        hs_sum1_relative.SetMaximum(ratio_max)
        hs_sum1_relative.SetMinimum(ratio_min)

        hs_sum1_relative.Divide(hs_sum2)

        #h2.GetYaxis()->SetLabelOffset(0.01)

        hs_sum1_relative.GetXaxis().SetLabelFont(63)
        hs_sum1_relative.GetXaxis().SetLabelSize(14) # labels will be 14 pixels
        hs_sum1_relative.GetYaxis().SetLabelFont(63)
        hs_sum1_relative.GetYaxis().SetLabelSize(14) # labels will be 14 pixels

        hs_sum1_relative.GetXaxis().SetTitleFont(63)
        hs_sum1_relative.GetXaxis().SetTitleSize(20)

        hs_sum1_relative.GetYaxis().SetTitleFont(63)
        hs_sum1_relative.GetYaxis().SetTitleSize(20)

        # if there is stack plot
        # removing the margin space to stack plot
        # and add the X label to the ratio
        if args.plot:
            #hs_sum1_relative.GetYaxis().SetLabelOffset(0.01)
            #histo_data_relative.GetYaxis().SetLabelOffset(0.01)
            hs_sum1_relative   .SetXTitle(title_x)
            hs_sum1_relative   .GetXaxis().SetTitleOffset(4.) # place the title not overlapping with labels...
            if not (args.no_data or args.no_data_plot):
                histo_data_relative.SetXTitle(title_x)
                histo_data_relative.GetXaxis().SetTitleOffset(4.)

        hs_sum1_relative   .GetYaxis().SetTitleOffset(1.4) # place the title not overlapping with labels...

        hs_sum1_relative   .SetYTitle("Data/MC")

        hs_sum1_relative.Draw("e2")
        if not (args.no_data or args.no_data_plot):
            histo_data_relative.SetYTitle("Data/MC")
            histo_data_relative.GetYaxis().SetTitleOffset(1.4)
            histo_data_relative.Draw("e p same")

    if args.plot:
        pad1.cd()

        max_y = hs_sum2.GetMaximum() if (args.no_data or args.no_data_plot) else max([h.GetMaximum() for h in (hs_sum2, histos_data_sum)])
        if args.y_max:
            max_y = args.y_max
        min_y = hs_sum2.GetMinimum() if (args.no_data or args.no_data_plot) else max([h.GetMinimum() for h in (hs_sum2, histos_data_sum)])

        if args.y_range:
            min_y, max_y = [float(x) for x in args.y_range.split(',')]
            hs_sum2        .SetMaximum(max_y)
            hs_sum2        .SetMinimum(min_y)
            #hs_sum2.GetYaxis().SetRange(min_y, max_y)
            #hs_sum2.GetYaxis().SetRangeUser(min_y, max_y)
            hs_sum2         .SetAxisRange(min_y, max_y, "Y")
            if not (args.no_data or args.no_data_plot):
                #histos_data_sum.GetYaxis().SetRange(min_y, max_y)
                #histos_data_sum.GetYaxis().SetRangeUser(min_y, max_y)
                histos_data_sum.SetMaximum(max_y)
                histos_data_sum.SetMinimum(min_y)
                histos_data_sum .SetAxisRange(min_y, max_y, "Y")

        # remove the label on stack plot if ratio is there
        if args.ratio:
            hs_sum2.GetXaxis().SetLabelOffset(999)
            hs_sum2.GetXaxis().SetLabelSize(0)
            if not (args.no_data or args.no_data_plot):
                histos_data_sum.GetXaxis().SetLabelOffset(999)
                histos_data_sum.GetXaxis().SetLabelSize(0)
            #hs.GetXaxis().SetLabelOffset(999)
            #hs.GetXaxis().SetLabelSize(0)
        else:
            if not (args.no_data or args.no_data_plot):
                histos_data_sum.SetXTitle(title_x)
            #hs            .SetXTitle(title_x)
            hs_sum2        .SetXTitle(title_x)

        if not (args.no_data or args.no_data_plot):
            histos_data_sum.SetMaximum(max_y * 1.1)
        if args.logy: # and args.fake_rate:
            logging.info("setting histos logy") # some bug 
            hs_sum2        .SetMaximum(max_y * 10.)
            hs_sum2        .SetMinimum(min_y / 10.)
            if not (args.no_data or args.no_data_plot):
                histos_data_sum.SetMaximum(max_y * 10.)
                histos_data_sum.SetMinimum(min_y / 10.)
        else: # if not args.logy:
            hs_sum2   .SetMinimum(0)
            if not (args.no_data or args.no_data_plot):
                histos_data_sum.SetMinimum(0)

        #hs.GetYaxis().SetTitleFont(63)
        #hs.GetYaxis().SetTitleSize(20)
        hs_sum2.GetYaxis().SetTitleFont(63)
        hs_sum2.GetYaxis().SetTitleSize(20)
        if not (args.no_data or args.no_data_plot):
            histos_data_sum.GetYaxis().SetTitleFont(63)
            histos_data_sum.GetYaxis().SetTitleSize(20)

        #histos_data_sum.GetYaxis().SetLabelFont(63)
        #histos_data_sum.GetYaxis().SetLabelSize(14) # labels will be 14 pixels

        # axis labels
        #hs_sum1.GetYaxis().SetLabelSize(0.02)
        #hs_sum1.GetXaxis().SetLabelSize(0.02)
        #histos_data_sum.GetYaxis().SetLabelSize(0.02)
        #histos_data_sum.GetXaxis().SetLabelSize(0.02)

        hs_sum2.GetYaxis().SetLabelFont(63)
        hs_sum2.GetXaxis().SetLabelFont(63)
        hs_sum2.GetYaxis().SetLabelSize(14)
        hs_sum2.GetXaxis().SetLabelSize(14)
        if not (args.no_data or args.no_data_plot):
            histos_data_sum.GetYaxis().SetLabelFont(63)
            histos_data_sum.GetXaxis().SetLabelFont(63)
            histos_data_sum.GetYaxis().SetLabelSize(14)
            histos_data_sum.GetXaxis().SetLabelSize(14)

        #hs             .GetYaxis().SetTitleOffset(1.4)
        hs_sum2        .GetYaxis().SetTitleOffset(1.5) # place the title not overlapping with labels...

        #hs             .SetYTitle(title_y)
        hs_sum2        .SetYTitle(title_y)

        #hs             .SetTitle(title_plot)
        hs_sum2        .SetTitle(title_plot)

        if not (args.no_data or args.no_data_plot):
            histos_data_sum.GetYaxis().SetTitleOffset(1.5)
            histos_data_sum.SetYTitle(title_y)
            histos_data_sum.SetTitle(title_plot)

        if args.y_range:
            min_y, max_y = [float(x) for x in args.y_range.split(',')]
            hs_sum2         .SetAxisRange(min_y, max_y, "Y")
            if not (args.no_data or args.no_data_plot):
                histos_data_sum .SetAxisRange(min_y, max_y, "Y")

        # damn root's inability to adjust maxima and all these workarounds...
        if not (args.no_data or args.no_data_plot):
            histos_data_sum.Draw("e1 p")
            if not args.fake_rate:
                hs.Draw("same")

            # MC sum plot
            if args.fake_rate:
                hs_sum2.SetFillStyle(3004);
                hs_sum2.SetFillColor(1);
                #hs_sum2.SetMarkerColorAlpha(0, 0.1);
                hs_sum2.SetMarkerStyle(25);
                hs_sum2.SetMarkerColor(kRed);
                hs_sum2.SetLineColor(kRed);
                #hs_sum2.Draw("same e2")
                hs_sum2.Draw("same e")
            elif args.cumulative_fractions:
                #for bini in range(hs_sum2.GetSize()):
                #    hs_sum2.SetBinError(bini, 0.)
                #hs_sum2.Draw("same")
                pass
            else:
                # only error band in usual case
                hs_sum2.Draw("same e2")

            histos_data_sum.Draw("same e1p")
        else:
            # the histogramStack cannot have title in root... therefore it cannot be plotted first..
            # thus I have to plot sum of MC first to get the titles right..
            if args.fake_rate:
                hs_sum2.Draw("e")
            elif args.cumulative_fractions:
                hs_sum2.Draw()
            else:
                # only error band in usual case
                hs_sum2.Draw("e2")

            if not args.fake_rate:
                hs.Draw("same")

            # MC sum plot
            if args.fake_rate:
                hs_sum2.SetFillStyle(3004);
                hs_sum2.SetFillColor(1);
                #hs_sum2.SetMarkerColorAlpha(0, 0.1);
                hs_sum2.SetMarkerStyle(25);
                hs_sum2.SetMarkerColor(kRed);
                hs_sum2.SetLineColor(kRed);
                #hs_sum2.Draw("same e2")
                hs_sum2.Draw("same e")
            elif args.cumulative_fractions:
                hs_sum2.Draw("same")
            else:
                # only error band in usual case
                hs_sum2.Draw("same e2")

    '''
    cout << "setting title" << endl;
    // title for the plot
    TPaveText* left_title = new TPaveText(0.1, 0.9, 0.4, 0.94, "brNDC");
    left_title->AddText("CMS preliminary at 13 TeV");
    left_title->SetTextFont(1);
    left_title->SetFillColor(0);
    cout << "drawing left title" << endl;
    left_title->Draw("same");
    '''

    if args.exp_legend:
        #left_title = TPaveText(0.05, 0.9, 0.4, 0.94, "brNDC")
        left_title = TPaveText(0.12, 0.8, 0.35, 0.88, "brNDC")
    else:
        #left_title = TPaveText(0.1, 0.9, 0.4, 0.94, "brNDC")
        #left_title = TPaveText(0.12, 0.8, 0.35, 0.88, "brNDC")
        left_title = TPaveText(0.1, 0.8, 0.25, 0.88, "brNDC")

    if args.no_data or args.no_data_plot:
        left_title.AddText("CMS simulation")
    else:
        left_title.AddText("CMS")
    left_title.SetTextFont(1)
    left_title.SetFillColor(0)

    '''
    TPaveText* right_title = new TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC");
    TString s_title(""); s_title.Form("L = %.1f fb^{-1}", lumi/1000);
    right_title->AddText(s_title);
    right_title->SetTextFont(132);
    right_title->SetFillColor(0);
    cout << "drawing right title" << endl;
    right_title->Draw("same");
    '''

    if args.exp_legend:
        #right_title = TPaveText(0.85, 0.9, 1.0, 0.96, "brNDC")
        right_title = TPaveText(0.75, 0.9, 0.98, 0.98, "brNDC")
    else:
        right_title = TPaveText(0.65, 0.9, 0.9,  0.95, "brNDC")
    right_title.AddText("%s fb^{-1} (13 TeV)" % (args.lumi / 1000. if args.lumi else args.lumi_label))
    right_title.SetTextFont(132)
    right_title.SetFillColor(0)

    left_title .Draw("same")
    right_title.Draw("same")

    # trying to get the ticks on the pad axes to be in the foreground
    #pad1.Draw("same")
    #cst.Update()
    #cst.Draw("same")
    #pad1.DrawFrame()
    ROOT.gPad.RedrawAxis()

    cst.cd()
    if not args.fake_rate and not args.skip_legend:
       if args.exp_legend: leg.SetBorderSize(0)
       leg.Draw("same")

    if args.output_name:
        filename = out_dir + args.output_name + '.png'
    else:
        stack_or_ratio = ('_stack' if args.plot else '') + ('_ratio' if args.ratio else '')
        shape_chan = ('_x_' + args.shape) if args.shape else ''
        filename = out_dir + '_'.join((args.mc_file.replace('/', ',').split('.root')[0], args.data_file.replace('/', ',').split('.root')[0], distr_names[0], channel, sys_name)) + stack_or_ratio + shape_chan + ('_dataqcd' if args.qcd > 0. else '') + ('_fakerate' if args.fake_rate else '') + ('_cumulative' if args.cumulative else '') + ('_cumulative-fractions' if args.cumulative_fractions else '') + ('_logy' if args.logy else '') + ('_normalize' if args.normalize else '') + ('_nolegend' if args.skip_legend else '') + ('_noQCD' if args.skip_QCD else '') + ('_nodata' if args.no_data else '') + ('_nodataplot' if args.no_data_plot else '') + args.output_suffix + ".png"

    if isfile(filename) and not args.overwrite:
        print filename, "exists, nothing written"
    else:
        cst.SaveAs(filename)


logging.info("segfault after here")

