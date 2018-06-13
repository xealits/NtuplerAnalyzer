import argparse
import logging
from os.path import isfile

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

parser.add_argument("--ratio-range", type=float, default=0.5, help="range of ratio plot (1-range 1+range)")

parser.add_argument("--lumi", type=float, help="use to skip the final normalizing step")

parser.add_argument("--y-max",   type=float, help="set the maximum on Y axis")
parser.add_argument("--x-range", type=str,   help="set the range on X axis")
parser.add_argument("--y-range",   type=str, help="set the range on Y axis `0,1.5`")

#parser.add_argument("--Vloose-shapes",   action='store_true', help="get TT proc shapes from selVloose")

parser.add_argument("--exp-legend",   action='store_true', help="experimentary legend drawing")

parser.add_argument("--sort-dy",   action='store_true', help="sort DY up")

parser.add_argument("--no-data",   action='store_true', help="don't draw data")

parser.add_argument("--bin-norm",  action='store_true', help="normalize per bin width")

parser.add_argument("--qcd",   type=float, default=0., help="get QCD from corresponding _ss channel and transfer it with the given factor (try --qcd 1)")
parser.add_argument("--qcd-factor",   type=float, default=1., help="factor for MC QCD")
parser.add_argument("--osss",     action='store_true', help="plot the ratio in OS/SS of data-other bkcg")
parser.add_argument("--osss-mc",  action='store_true', help="plot the ratio in OS/SS of MC QCD")

parser.add_argument("--wjets", type=float, default=0., help="scale factor for wjets (from the control region)")
parser.add_argument("--factor-dy", type=float, help="scale factor for dy (from the control region)")

parser.add_argument("--draw-overflows", type=float, default=0., help="draw the overflow bin with the specified width")

parser.add_argument("--fake-rate", action='store_true', help='(ad-hocish) the fake rate in MC and data, i.e. ratio of whatever two distributions given as "nominator/denominator"')

parser.add_argument("--cumulative", action='store_true', help="plot stack of cumulative distribution (from right)")
parser.add_argument("--cumulative-fractions", action='store_true', help="cumulative distributions, each bin normalized to 1")

#parser.add_argument("-f", "--form-shapes", type=str, default='usual', help="plot the shapes of distributions normalized to 1")
parser.add_argument("-f", "--form-shapes", action='store_true', help="plot the shapes of distributions normalized to 1")
#parser.add_argument("-e", "--shape-evolution", type=str, default='usual', help="plot the same distr in different channels normalized to 1")
parser.add_argument("--processes", type=str, default='all', help="set processes to consider (all by default)")

parser.add_argument("--skip-legend", action='store_true', help="don't plot the legend (full view of distribution)")
parser.add_argument("--skip-QCD", action='store_true', help="skip MC QCD")

parser.add_argument("--legend-shift", type=float, help="shift legend on X axis")

parser.add_argument("--drop-bin", type=float, help="drop values in the bin")

parser.add_argument("--lumi-label", type=float, default=35.8, help="set lumi label on the plot")
parser.add_argument("--title-x", type=str, default="", help="set title of X axis on the plot")
parser.add_argument("--title-y", type=str, default="", help="set title of Y axis on the plot")
parser.add_argument("--title",   type=str, default="default", help="set title of the plot")

parser.add_argument("--output-name",   type=str, help="name for the output file")

args = parser.parse_args()

assert isfile(args.mc_file)
assert isfile(args.data_file)

logging.info("import ROOT")

import ROOT
from ROOT import TFile, THStack, TLegend, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
import plotting_root
from draw_overflows import DrawOverflow

#channel = "mu_presel"
#channel = "mu_sel"
#channel = "mu_lj"
#channel = "mu_lj_out"

#sys_name = "NOMINAL"
sys_name = args.systematic
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
for distr_name in distr_names:
    data_distr_name = distr_name[:]
    if data_distr_name in ('Mt_lep_met_f_w_mu_trk_b', 'Mt_lep_met_f_w_mu_trk_h'):
        data_distr_name = 'Mt_lep_met_f'
    data_distr_name.replace('_w_pu_sum', '_w_pu').replace('_w_pu_b', '_w_pu').replace('_w_pu_h2', '_w_pu').replace('_pu_h2', '').replace('_pu_b', '').replace('_pu_sum', '')
    # get overflows if requested
    if args.draw_overflows:
        #h_overflows = draw_overflows.DrawOverflow(histo, args.draw_overflows)
        histos_data_distrs.append((data_distr_name, [(DrawOverflow(fdata.Get(channel + '/data/NOMINAL/' + '_'.join([channel, 'data', 'NOMINAL', data_distr_name])), args.draw_overflows), 'data', channel) for channel in channels]))
    else:
        histos_data_distrs.append((data_distr_name, [(fdata.Get(channel + '/data/NOMINAL/' + '_'.join([channel, 'data', 'NOMINAL', data_distr_name])), 'data', channel) for channel in channels]))

    if args.qcd > 0. or args.osss or args.osss_mc:
        histos_data_distrs_ss.append((data_distr_name, [(fdata.Get(channel + '_ss' + '/data/NOMINAL/' + '_'.join([channel + '_ss', 'data', 'NOMINAL', data_distr_name])), 'data', channel) for channel in channels]))

if args.cumulative:
    histos_data_per_distr    = [(name, [(histo.GetCumulative(False), n, c) for histo, n, c in distrs]) for name, distrs in histos_data_distrs]
    histos_data_per_distr_ss = [(name, [(histo.GetCumulative(False), n, c) for histo, n, c in distrs]) for name, distrs in histos_data_distrs_ss]
else:
    histos_data_per_distr    = histos_data_distrs
    histos_data_per_distr_ss = histos_data_distrs_ss


# normalize data distrs to different bin width:
if args.bin_norm:
  for _, distrs in histos_data_per_distr + histos_data_per_distr_ss:
      for histo, _, _ in distrs:
          for bini in range(histo.GetSize()):
              content = histo.GetBinContent(bini)
              error   = histo.GetBinError(bini)
              width   = histo.GetXaxis().GetBinUpEdge(bini) - histo.GetXaxis().GetBinLowEdge(bini)
              histo.SetBinContent(bini, content/width)
              histo.SetBinError(bini, error/width)

logging.info("# data histograms = %d" % len(histos_data_per_distr))

syst_factors = {
              "FragUp":             1.001123,
            "FragDown":             0.999304,
         "SemilepBRUp":             1.008565,
       "SemilepBRDown":             0.987050,
          "PetersonUp":             1.000032,
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

def get_histos(infile, channels, shape_channel, sys_name, distr_name, skip_QCD=False):
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
           if skip_QCD and nick == 'qcd':
               continue
           if processes_requirement and nick not in processes_requirement:
               continue
           #logging.info(nick)

           if nick == 'data':
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
               tau_selection = any(ch in channel for ch in ('ctr_old_mu_sel', 'ctr_old_mu_sel', 'ctr_old_el_sel', 'optel_tight_2L1M', 'optmu_tight_2L1M'))
               tau_process   = nick in ('tt_mutau3ch', 'tt_eltau3ch', 'tt_mutau', 'tt_eltau', 'tt_taultauh', 'dy_tautau', 's_top_eltau', 's_top_mutau')
               if tau_selection and tau_process:
                   tauIDSF_factor = 0.95

               pu_factor = 1.
               if 'PUUp' in fixed_sys_name:
                   pu_factor = 1. / (1.01995 if '_el_' in channel else 0.9979) # 0.97 # 1./ 0.9979
               elif 'PUDown' in fixed_sys_name:
                   pu_factor = 1. / (1.07766 if '_el_' in channel else 1.0485) # 1.17 # 1./ 1.485
               else:
                   pu_factor = 1. / (1.04678 if '_el_' in channel else 1.022)  # 1.06 # 1./ 1.02135 an 1/1.014 with weight counter..
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

               final_factor = args.lumi * tauIDSF_factor * pu_factor * th_factor
               logging.info("final factor %20s %20f %5f %5f %5f  %20f" % (nick, args.lumi, tauIDSF_factor, pu_factor, th_factor, final_factor))
               histo.Scale(final_factor)

           # wjets normalization
           if args.wjets > 0 and 'wjets' in nick:
               histo.Scale(args.wjets)

           # dy normalization
           if args.factor_dy and 'dy' in nick:
               histo.Scale(args.factor_dy)

           # mc qcd normalization
           if nick == 'qcd':
               histo.Scale(args.qcd_factor)

           # rename the histo for the correct systematic name in case of TT systematics
           histo_name = '_'.join([channel, nick, sys_name, distr_name])
           histo.SetName(histo_name)

           #histo = histo_key.ReadObj()
           logging.info("%s   %s   %x = %f %f" % (histo.GetName(), histo_name, histo_name == '_'.join([channel, nick, fixed_sys_name, distr_name]), histo.GetEntries(), histo.Integral()))

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

           if args.cumulative:
               used_histos.append((histo.GetCumulative(False), nick, channel)) # hopefully root wont screw this up
           else:
               used_histos.append((histo, nick, channel)) # hopefully root wont screw this up

    return used_histos


f = TFile(args.mc_file)
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
    # get SS data and MC
    used_histos_per_distr_ss = get_histos(f, [channel + '_ss'], args.shape, sys, distr)
    data_hist = fdata.Get(channel + '_ss' + '/data/NOMINAL/' + '_'.join([channel + '_ss', 'data', 'NOMINAL', distr]))

    # find QCD shape in SS
    mc_hist = used_histos_per_distr_ss[0][0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
    mc_hist.SetName('mc_sum2_ss')
    for h, nick, channel in used_histos_per_distr_ss[1:]:
        if nick == 'qcd': continue # don't include QCD MC -- its sum of no-QCD MC
        mc_hist.Add(h)

    qcd_hist = data_hist - mc_hist
    # scale by the given factor and return
    qcd_hist.Scale(args.qcd)

    return qcd_hist


if args.qcd > 0. or args.osss or args.osss_mc:
    # get all the same distributions for _ss channel
    #used_histos_per_distr_ss = [(distr_name, get_histos(f, [c + '_ss' for c in channels], args.shape, sys_name, distr_name)) for distr_name in distr_names]
    used_histos_per_distr_ss = [(distr_name, get_histos(f, [c + '_ss' for c in channels], args.shape, 'NOMINAL', distr_name)) for distr_name in distr_names]

    for distr, histos in used_histos_per_distr_ss:
        # loop through normal list
        hs_sum2 = histos[0][0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
        hs_sum2.SetName('mc_sum2_ss')

        for h, nick, channel in histos[1:]: # I sum different channels?... the old hacked not fixed yet
            if nick == 'qcd': continue # don't include QCD MC -- its sum of no-QCD MC
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
            qcd_hist = data_hist - mc_hist
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
                shape_qcd = datadriven_qcd(args.shape, distr, 'NOMINAL') # nominal only for now
                shape_qcd.Scale(qcd_hist.Integral() / shape_qcd.Integral())
                qcd_hist = shape_qcd
            
            for bini in range(qcd_hist.GetSize()):
                if qcd_hist.GetBinContent(bini) < 0:
                    qcd_hist.SetBinContent(bini, 0)

            #datadriven_qcd_name = "qcd"
            # optmu_tight_2L1M_wjets_NOMINAL_Mt_lep_met
            # I remove _ss from channel name
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
        if args.qcd > 0. and nick == 'qcd':
            # then take the qcd from the differences in ss region
            h = qcd_hists[(distr, channel + '_ss')]
            # and substitute it in the list (if it won't brak everything)
            histos[i] = (h, nick, channel) # h is new here!
        hs_sum2.Add(h)
        if nick != 'qcd':
            hs_sum_noqcd.Add(h)

    hs_sum2.SetFillStyle(3004);
    hs_sum2.SetFillColor(1);
    hs_sum2.SetMarkerStyle(1)
    hs_sum2.SetMarkerColor(0)

    hs_sums2.append(hs_sum2)
    hs_sums_noqcd.append(hs_sum_noqcd)

hs_sum2 = hs_sums2[0]
if args.qcd > 0. or args.osss or args.osss_mc:
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
histos_data_sum = histos_data_sums_per_distr[0]

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

# get MC stack and legend for it
#hs, leg = plotting_root.stack_n_legend(used_histos)
if args.legend_shift:
    hs_legs_per_distr = [(distr_name, plotting_root.stack_n_legend(histos, args.legend_shift, args.exp_legend, sort_dy=args.sort_dy)) for distr_name, histos in used_histos_per_distr]
else:
    hs_legs_per_distr = [(distr_name, plotting_root.stack_n_legend(histos, exp_legend=args.exp_legend, sort_dy=args.sort_dy)) for distr_name, histos in used_histos_per_distr]

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

logging.info("data   = %f" % histos_data_sum.Integral())
logging.info("mc sum = %f %f" % (hs_sum1.Integral(), hs_sum2.Integral()))

# set data histo styles and
# add data to the legend
#for histo_data in histos_data:
#    histo_data.SetMarkerStyle(21)
#    leg.AddEntry(histo_data, "data", "e1 p")
# no, add just the data_sum entry
histos_data_sum.SetMarkerStyle(21)
if not args.no_data:
    leg.AddEntry(histos_data_sum, "data", "e1 p")

out_dir = args.output_directory + '/' if args.output_directory else './'

histos_data = histos_data_per_distr[0][1]

histos_data[0][0].GetXaxis().SetLabelFont(63)
histos_data[0][0].GetXaxis().SetLabelSize(14) # labels will be 14 pixels
histos_data[0][0].GetYaxis().SetLabelFont(63)
histos_data[0][0].GetYaxis().SetLabelSize(14) # labels will be 14 pixels


if not args.plot and not args.ratio:
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

    filename = out_dir + args.mc_file.split('.root')[0] + '_HISTOSEL_%s_%s_%s%s.root' % (distr_name, channel, sys_name, options)
    logging.info('saving root %s' % filename)
    #fout = TFile(filename, 'RECREATE')
    fout = TFile(filename, 'CREATE')
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

    logging.info("no segfault here")

elif args.form_shapes:
    from ROOT import gStyle, gROOT, TCanvas, TPad
    gROOT.SetBatch()
    gStyle.SetOptStat(0)
    cst = TCanvas("cst","stacked hists",10,10,700,700)
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

    histos_data[0][0].SetMaximum(max_y * 1.1)
    histos_data[0][0].SetXTitle(distr_name)

    if not args.no_data:
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

elif args.osss or args.osss_mc:
    from ROOT import gStyle, gROOT, TCanvas, TPad
    gROOT.SetBatch()
    gStyle.SetOptStat(0)
    cst = TCanvas("cst","stacked hists",10,10,700,700)

    # the no-qcd MC sum is in
    # hs_sum_noqcd and hs_sum_noqcd_ss
    # the data is in
    # histos_data_per_distr and histos_data_per_distr_ss
    #histos_data_distrs = [(distr_name, [(histo, 'data', channel)]]

    if args.osss:
        histo_diff_os = histos_data_per_distr   [0][1][0][0] - hs_sum_noqcd
        histo_diff_ss = histos_data_per_distr_ss[0][1][0][0] - hs_sum_noqcd_ss
    elif args.osss_mc:
        #used_histos_per_distr = [(distr_name, get_histos(f, channels, args.shape, sys_name, distr_name, skip_QCD = args.skip_QCD)) for distr_name in distr_names]
        # [(distr_name, [h, nick, channel in histos])]
        histo_diff_os = [h for h, nick, _ in used_histos_per_distr   [0][1] if nick == 'qcd'][0]
        histo_diff_ss = [h for h, nick, _ in used_histos_per_distr_ss[0][1] if nick == 'qcd'][0]

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

    if args.y_range:
        min_y, max_y = [float(x) for x in args.y_range.split(',')]
        histo_diff_os.SetMaximum(max_y)
        histo_diff_os.SetMinimum(min_y)

    histo_diff_os.SetXTitle(title_x)
    histo_diff_os.SetYTitle(title_y)
    histo_diff_os.SetTitle(title_plot)
    histo_diff_os.Draw()

    left_title = TPaveText(0.1, 0.9, 0.4, 0.94, "brNDC")
    left_title.AddText("CMS preliminary at 13 TeV")
    left_title.SetTextFont(1)

    right_title = TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC")
    right_title.AddText("L = %s fb^{-1}" % (args.lumi / 1000. if args.lumi else args.lumi_label))
    right_title.SetTextFont(132)
    right_title.SetFillColor(0)

    left_title .Draw("same")
    right_title.Draw("same")

    label = TPaveText(0.5, 0.2, 0.8, 0.3, "brNDC")
    label.AddText("anti-iso region")
    label.SetTextFont(132)
    label.SetFillColor(0)
    label.Draw("same")

    cst.SaveAs(out_dir + '_'.join((args.mc_file.replace('/', ',').split('.root')[0], args.data_file.replace('/', ',').split('.root')[0], distr_name, channel, sys_name)) + ('_osss-factor' if args.osss else '_osss-mc-factor') + '.png')

else:
    '''
    plot stack plot and/or ratio plot
    '''

    from ROOT import gStyle, gROOT, TCanvas, TPad
    gROOT.SetBatch()
    gStyle.SetOptStat(0)
    if args.exp_legend:
        cst = TCanvas("cst","stacked hists",10,10,1000,700)
    else:
        cst = TCanvas("cst","stacked hists",10,10,700,700)

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
        hs_sum1.Divide(hs_sum1_2)

    if args.x_range:
        x_min, x_max = (float(x) for x in args.x_range.split(','))
        histos_data_sum .SetAxisRange(x_min, x_max, "X")
        hs_sum1         .SetAxisRange(x_min, x_max, "X")

    # drop bin
    if args.drop_bin:
        bn_to_drop = histos_data_sum.FindBin(args.drop_bin)
        histos_data_sum.SetBinContent (bn_to_drop, 0)
        histos_data_sum.SetBinError   (bn_to_drop, 0)

    # plotting
    if args.ratio:
        pad2.cd()

        # calculate ratios
        histo_data_relative = histos_data_sum.Clone()
        hs_sum1_relative = hs_sum1.Clone()
        histo_data_relative.SetName("rel_data")
        hs_sum1_relative.SetName("rel_mc")

        histo_data_relative.SetStats(False)
        hs_sum1_relative.SetStats(False)

        #histo_data_relative.GetYaxis().SetRange(0.5, 1.5)
        #histo_data_relative.GetYaxis().SetUserRange(0.5, 1.5)
        #hs_sum1_relative.GetYaxis().SetRange(0.5, 1.5)
        #hs_sum1_relative.GetYaxis().SetUserRange(0.5, 1.5)

        ratio_max = 1. + args.ratio_range
        ratio_min = 1. - args.ratio_range
        histo_data_relative.SetMaximum(ratio_max)
        histo_data_relative.SetMinimum(ratio_min)
        hs_sum1_relative.SetMaximum(ratio_max)
        hs_sum1_relative.SetMinimum(ratio_min)

        histo_data_relative.Divide(hs_sum1)
        hs_sum1_relative.Divide(hs_sum1)

        #h2.GetYaxis()->SetLabelOffset(0.01)

        histo_data_relative.GetXaxis().SetLabelFont(63)
        histo_data_relative.GetXaxis().SetLabelSize(14) # labels will be 14 pixels
        histo_data_relative.GetYaxis().SetLabelFont(63)
        histo_data_relative.GetYaxis().SetLabelSize(14) # labels will be 14 pixels

        hs_sum1_relative.GetXaxis().SetLabelFont(63)
        hs_sum1_relative.GetXaxis().SetLabelSize(14) # labels will be 14 pixels
        hs_sum1_relative.GetYaxis().SetLabelFont(63)
        hs_sum1_relative.GetYaxis().SetLabelSize(14) # labels will be 14 pixels

        hs_sum1_relative.GetXaxis().SetTitleFont(63)
        hs_sum1_relative.GetXaxis().SetTitleSize(20)
        histo_data_relative.GetXaxis().SetTitleFont(63)
        histo_data_relative.GetXaxis().SetTitleSize(20)

        hs_sum1_relative.GetYaxis().SetTitleFont(63)
        hs_sum1_relative.GetYaxis().SetTitleSize(20)
        histo_data_relative.GetYaxis().SetTitleFont(63)
        histo_data_relative.GetYaxis().SetTitleSize(20)

        # if there is stack plot
        # removing the margin space to stack plot
        # and add the X label to the ratio
        if args.plot:
            #hs_sum1_relative.GetYaxis().SetLabelOffset(0.01)
            #histo_data_relative.GetYaxis().SetLabelOffset(0.01)
            hs_sum1_relative   .SetXTitle(title_x)
            histo_data_relative.SetXTitle(title_x)
            hs_sum1_relative   .GetXaxis().SetTitleOffset(4.) # place the title not overlapping with labels...
            histo_data_relative.GetXaxis().SetTitleOffset(4.)

        hs_sum1_relative   .GetYaxis().SetTitleOffset(1.4) # place the title not overlapping with labels...
        histo_data_relative.GetYaxis().SetTitleOffset(1.4)

        hs_sum1_relative   .SetYTitle("Data/MC")
        histo_data_relative.SetYTitle("Data/MC")

        hs_sum1_relative.Draw("e2")
        if not args.no_data:
            histo_data_relative.Draw("e p same")

    if args.plot:
        pad1.cd()

        max_y = max([h.GetMaximum() for h in (hs_sum1, histos_data_sum)])
        if args.y_max:
            max_y = args.y_max
        min_y = max([h.GetMinimum() for h in (hs_sum1, histos_data_sum)])

        if args.y_range:
            min_y, max_y = [float(x) for x in args.y_range.split(',')]
            histos_data_sum.SetMaximum(max_y)
            hs_sum1        .SetMaximum(max_y)
            histos_data_sum.SetMinimum(min_y)
            hs_sum1        .SetMinimum(min_y)
            #hs_sum1.GetYaxis().SetRange(min_y, max_y)
            #hs_sum1.GetYaxis().SetRangeUser(min_y, max_y)
            #histos_data_sum.GetYaxis().SetRange(min_y, max_y)
            #histos_data_sum.GetYaxis().SetRangeUser(min_y, max_y)
            hs_sum1         .SetAxisRange(min_y, max_y, "Y")
            histos_data_sum .SetAxisRange(min_y, max_y, "Y")

        # remove the label on stack plot if ratio is there
        if args.ratio:
            hs_sum1.GetXaxis().SetLabelOffset(999)
            hs_sum1.GetXaxis().SetLabelSize(0)
            histos_data_sum.GetXaxis().SetLabelOffset(999)
            histos_data_sum.GetXaxis().SetLabelSize(0)
            #hs.GetXaxis().SetLabelOffset(999)
            #hs.GetXaxis().SetLabelSize(0)
        else:
            histos_data_sum.SetXTitle(title_x)
            #hs            .SetXTitle(title_x)
            hs_sum1        .SetXTitle(title_x)

        histos_data_sum.SetMaximum(max_y * 1.1)
        if args.logy: # and args.fake_rate:
            logging.info("setting histos logy") # some bug 
            histos_data_sum.SetMaximum(max_y * 10.)
            hs_sum1        .SetMaximum(max_y * 10.)
            histos_data_sum.SetMinimum(min_y / 10.)
            hs_sum1        .SetMinimum(min_y / 10.)
        else: # if not args.logy:
            histos_data_sum.SetMinimum(0)
            hs_sum1   .SetMinimum(0)

        #hs.GetYaxis().SetTitleFont(63)
        #hs.GetYaxis().SetTitleSize(20)
        hs_sum1.GetYaxis().SetTitleFont(63)
        hs_sum1.GetYaxis().SetTitleSize(20)
        histos_data_sum.GetYaxis().SetTitleFont(63)
        histos_data_sum.GetYaxis().SetTitleSize(20)

        # axis labels
        hs_sum1.GetYaxis().SetLabelSize(0.02)
        hs_sum1.GetXaxis().SetLabelSize(0.02)
        histos_data_sum.GetYaxis().SetLabelSize(0.02)
        histos_data_sum.GetXaxis().SetLabelSize(0.02)

        #hs             .GetYaxis().SetTitleOffset(1.4)
        hs_sum1        .GetYaxis().SetTitleOffset(1.5) # place the title not overlapping with labels...
        histos_data_sum.GetYaxis().SetTitleOffset(1.5)

        #hs             .SetYTitle(title_y)
        hs_sum1        .SetYTitle(title_y)
        histos_data_sum.SetYTitle(title_y)

        #hs             .SetTitle(title_plot)
        hs_sum1        .SetTitle(title_plot)
        histos_data_sum.SetTitle(title_plot)

        if args.y_range:
            min_y, max_y = [float(x) for x in args.y_range.split(',')]
            hs_sum1         .SetAxisRange(min_y, max_y, "Y")
            histos_data_sum .SetAxisRange(min_y, max_y, "Y")

        # damn root's inability to adjust maxima and all these workarounds...
        if not args.no_data:
            histos_data_sum.Draw("e1 p")
            if not args.fake_rate:
                hs.Draw("same")

            # MC sum plot
            if args.fake_rate:
                hs_sum1.SetFillStyle(3004);
                hs_sum1.SetFillColor(1);
                #hs_sum1.SetMarkerColorAlpha(0, 0.1);
                hs_sum1.SetMarkerStyle(25);
                hs_sum1.SetMarkerColor(kRed);
                hs_sum1.SetLineColor(kRed);
                #hs_sum1.Draw("same e2")
                hs_sum1.Draw("same e")
            else:
                # only error band in usual case
                hs_sum1.Draw("same e2")

            histos_data_sum.Draw("same e1p")
        else:
            # the histogramStack cannot have title in root... therefore it cannot be plotted first..
            # thus I have to plot sum of MC first to get the titles right..
            if args.fake_rate:
                hs_sum1.Draw("e")
            else:
                # only error band in usual case
                hs_sum1.Draw("e2")

            if not args.fake_rate:
                hs.Draw("same")

            # MC sum plot
            if args.fake_rate:
                hs_sum1.SetFillStyle(3004);
                hs_sum1.SetFillColor(1);
                #hs_sum1.SetMarkerColorAlpha(0, 0.1);
                hs_sum1.SetMarkerStyle(25);
                hs_sum1.SetMarkerColor(kRed);
                hs_sum1.SetLineColor(kRed);
                #hs_sum1.Draw("same e2")
                hs_sum1.Draw("same e")
            else:
                # only error band in usual case
                hs_sum1.Draw("same e2")

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
        left_title = TPaveText(0.05, 0.9, 0.4, 0.94, "brNDC")
    else:
        left_title = TPaveText(0.1, 0.9, 0.4, 0.94, "brNDC")
    left_title.AddText("CMS preliminary at 13 TeV")
    left_title.SetTextFont(1)

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
        right_title = TPaveText(0.85, 0.9, 1.0, 0.94, "brNDC")
    else:
        right_title = TPaveText(0.75, 0.9, 0.9, 0.94, "brNDC")
    right_title.AddText("L = %s fb^{-1}" % (args.lumi / 1000. if args.lumi else args.lumi_label))
    right_title.SetTextFont(132)
    right_title.SetFillColor(0)

    left_title .Draw("same")
    right_title.Draw("same")

    cst.cd()
    if not args.fake_rate and not args.skip_legend:
       if args.exp_legend: leg.SetBorderSize(0)
       leg.Draw("same")

    if args.output_name:
        cst.SaveAs(out_dir + args.output_name + '.png')
    else:
        stack_or_ratio = ('_stack' if args.plot else '') + ('_ratio' if args.ratio else '')
        shape_chan = ('_x_' + args.shape) if args.shape else ''
        cst.SaveAs(out_dir + '_'.join((args.mc_file.replace('/', ',').split('.root')[0], args.data_file.replace('/', ',').split('.root')[0], distr_names[0], channel, sys_name)) + stack_or_ratio + shape_chan + ('_dataqcd' if args.qcd > 0. else '') + ('_fakerate' if args.fake_rate else '') + ('_cumulative' if args.cumulative else '') + ('_cumulative-fractions' if args.cumulative_fractions else '') + ('_logy' if args.logy else '') + ('_normalize' if args.normalize else '') + ('_nolegend' if args.skip_legend else '') + ('_noQCD' if args.skip_QCD else '') + ".png")


logging.info("segfault after here")

