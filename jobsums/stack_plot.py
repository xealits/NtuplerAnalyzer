import argparse
import logging
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

parser.add_argument("--fake-rate", action='store_true', help='(ad-hocish) the fake rate in MC and data, i.e. ratio of whatever two distributions given as "nominator/denominator"')

parser.add_argument("--cumulative", action='store_true', help="plot stack of cumulative distribution (from right)")
parser.add_argument("--cumulative-fractions", action='store_true', help="not implemented yet (cumulative distributions, each bin normalized to 1)")

#parser.add_argument("-f", "--form-shapes", type=str, default='usual', help="plot the shapes of distributions normalized to 1")
parser.add_argument("-f", "--form-shapes", action='store_true', help="plot the shapes of distributions normalized to 1")
#parser.add_argument("-e", "--shape-evolution", type=str, default='usual', help="plot the same distr in different channels normalized to 1")
parser.add_argument("--processes", type=str, default='all', help="set processes to consider (all by default)")

args = parser.parse_args()

logging.info("import ROOT")

import ROOT
from ROOT import TFile, THStack, TLegend, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
import plotting_root

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

if args.processes == 'all':
    processes_requirement = None
else:
    processes_requirement = args.processes.split(',')


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

histos_data_distrs = [(distr_name, [(fdata.Get(channel + '/data/NOMINAL/' + '_'.join([channel, 'data', 'NOMINAL', distr_name])), 'data', channel) for channel in channels]) for distr_name in distr_names]
if args.cumulative:
    histos_data_per_distr = [(name, [(histo.GetCumulative(False), n, c) for histo, n, c in distrs]) for name, distrs in histos_data_distrs]
else:
    histos_data_per_distr = histos_data_distrs

logging.info("# data histograms = %d" % len(histos_data_per_distr))

histos_data_sums_per_distr = []
for distr_name, histos_per_channel in histos_data_per_distr:
    histos_data_sum = histos_per_channel[0][0].Clone()
    histos_data_sum.Sumw2(ROOT.kTRUE) # are errors saved here?
    for h, _, _ in histos_per_channel[1:]:
        histos_data_sum.Add(h)
    histos_data_sums_per_distr.append(histos_data_sum)
histos_data_sum = histos_data_sums_per_distr[0]

def get_histos(infile, channels, shape_channel, sys_name, distr_name):
    """get_histos(infile)

    the file contains usual structure
    <channels>/<processes>/<systematics>/<histograms>

    returns [(histo, nick)]
    """
    used_histos = [] # (histo, nick)

    for channel in channels:
        chan = infile.Get(channel)
        processes_keys = list(chan.GetListOfKeys())
        #sorted_pkeys = sorted(processes_keys, key=lambda pkey: pkey.GetName() if pkey.GetName() not in ('qcd') else 'z_' + pkey.GetName())
        #for process in sorted_pkeys:
        for process in processes_keys:
           nick = process.GetName()
           if processes_requirement and nick not in processes_requirement:
               continue
           #logging.info(nick)

           if nick == 'data':
               continue

           fixed_sys_name = sys_name
           if any(sn in sys_name for sn in ['TOPPT', 'FSR', 'ISR', 'HDAMP', 'TuneCUETP8M2T4', 'QCDbasedCRTune', 'GluonMoveCRTune']):
               fixed_sys_name = sys_name if 'tt' in nick else 'NOMINAL'

           histo_name = '_'.join([channel, nick, fixed_sys_name, distr_name])
           logging.debug("%s" % histo_name)

           h_init = process.ReadObj().Get(fixed_sys_name + '/' + histo_name)

           if not h_init:
               logging.debug("absent %s" % histo_name)
               continue

           if shape_channel and nick in ('dy_other', 'dy_tautau', 'wjets'):
               histo_name = '_'.join([shape_channel, nick, fixed_sys_name, distr_name])
               h_shape_path = shape_channel + '/' + nick + '/' + fixed_sys_name + '/' + histo_name
               logging.info("shape from %s" % h_shape_path)

               h_shape = infile.Get(h_shape_path)
               histo = h_shape.Clone()
               histo.Scale(h_init.Integral() / h_shape.Integral())
           else:
               histo = h_init

           #histo = histo_key.ReadObj()
           logging.info("%s   %s   %x = %f %f" % (histo_name, histo_name, histo_name == '_'.join([channel, nick, fixed_sys_name, distr_name]), histo.GetEntries(), histo.Integral()))

           if args.cumulative:
               used_histos.append((histo.GetCumulative(False), nick, channel)) # hopefully root wont screw this up
           else:
               used_histos.append((histo, nick, channel)) # hopefully root wont screw this up

    return used_histos


f = TFile(args.mc_file)
used_histos_per_distr = [(distr_name, get_histos(f, channels, args.shape, sys_name, distr_name)) for distr_name in distr_names]

hs_sums2 = []

for _, histos in used_histos_per_distr:
    # loop through normal list
    hs_sum2 = histos[0][0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
    hs_sum2.SetName('mc_sum2')

    for h, _, _ in histos[1:]:
        hs_sum2.Add(h)

    hs_sum2.SetFillStyle(3004);
    hs_sum2.SetFillColor(1);
    hs_sum2.SetMarkerStyle(1)
    hs_sum2.SetMarkerColor(0)

    hs_sums2.append(hs_sum2)
hs_sum2 = hs_sums2[0]

if args.normalize:
    ratio = histos_data_sum.Integral() / hs_sum2.Integral()
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
hs_legs_per_distr = [(distr_name, plotting_root.stack_n_legend(histos)) for distr_name, histos in used_histos_per_distr]

# sum of MC to get the sum of errors

# loop through TList
hs_sums1 = []

for _, (hs, _) in hs_legs_per_distr:
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
leg.AddEntry(histos_data_sum, "data", "e1 p")

out_dir = args.output_directory + '/' if args.output_directory else './'

histos_data = histos_data_per_distr[0][1]

histos_data[0][0].GetXaxis().SetLabelFont(63)
histos_data[0][0].GetXaxis().SetLabelSize(14) # labels will be 14 pixels
histos_data[0][0].GetYaxis().SetLabelFont(63)
histos_data[0][0].GetYaxis().SetLabelSize(14) # labels will be 14 pixels


if not args.plot and not args.ratio:
    filename = out_dir + args.mc_file.split('.root')[0] + '_%s_%s_%s.root' % (distr_name, channel, sys_name)
    logging.info('saving root %s' % filename)
    fout = TFile(filename, 'RECREATE')
    fout.cd()
    #histo_data.Write() #'data')
    for h, nick, _ in histos_data + used_histos:
        logging.info(nick)
        h.Write()
    hs.Write()
    hs_sum1.Write()
    hs_sum2.Write()
    leg.Write('leg')
    fout.Write()
    fout.Close()

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
        histo.Scale(1/histo.Integral())

    gStyle.SetHatchesSpacing(2)

    # plot stuff
    pad.cd()

    # find max Y for plot
    max_y = max([h.GetMaximum() for h, _, _ in used_histos + histos_data])
    histos_data[0][0].SetMaximum(max_y * 1.1)
    histos_data[0][0].SetXTitle(distr_name)
    histos_data[0][0].Draw('e1 p')
    for i, (histo, nick, _) in enumerate(histos_data[1:] + used_histos):
        #histo.SetFillColor( # it shouldn't be needed with hist drawing option
        # nope, it's needed...
        histo.SetLineColor(plotting_root.nick_colour[nick])
        histo.SetLineWidth(4)
        if i < 2:
            histo.SetFillStyle([3354, 3345][i])
        else:
            histo.SetFillColorAlpha(0, 0.0)
        histo.Draw("hist same")
    histos_data[0][0].Draw('same e1 p')

    leg.Draw("same")

    cst.SaveAs(out_dir + '_'.join((args.mc_file.replace('/', ',').split('.root')[0], args.data_file.replace('/', ',').split('.root')[0], distr_name, channel, sys_name)) + '_shapes_%d-channels.png' % len(channels))

else:
    from ROOT import gStyle, gROOT, TCanvas, TPad
    gROOT.SetBatch()
    gStyle.SetOptStat(0)
    cst = TCanvas("cst","stacked hists",10,10,700,700)

    if args.ratio and args.plot:
        pad1 = TPad("pad1","This is pad1", 0., 0.3,  1., 1.)
        pad2 = TPad("pad2","This is pad2", 0., 0.05,  1., 0.3)
        #pad2.SetTopMargin(0.01) # doesn't work
        #gStyle.SetPadTopMargin(0.05) # nope
        #ROOT.gPad.SetTopMargin(0.01) # nope
        if args.logy:
            pad1.SetLogy()
        pad1.Draw()
        pad2.Draw() # these have to be before cd()
        # now set margins:
        pad1.cd()
        ROOT.gPad.SetBottomMargin(0.02)
        #ROOT.gPad.SetTopMargin(0.01)
        pad2.cd()
        #ROOT.gPad.SetBottomMargin(0.001)
        ROOT.gPad.SetTopMargin(0.02)
    elif args.ratio:
        pad2 = TPad("pad2","This is pad2", 0., 0.05,  1., 1.)
        pad2.Draw()
    else:
        pad1 = TPad("pad1","This is pad1", 0., 0.,  1., 1.)
        pad1.Draw()

    # if fake rate then the sums should be divided by second distr sum
    if args.fake_rate:
        histos_data_sum.Divide(histos_data_sum2)
        hs_sum1.Divide(hs_sum1_2)

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

        ratio_max = 2. if args.fake_rate else 1.5
        ratio_min = 0. if args.fake_rate else 0.5
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

        # if there is stack plot
        # removing the margin space to stack plot
        # and add the X label to the ratio
        if args.plot:
            #hs_sum1_relative.GetYaxis().SetLabelOffset(0.01)
            #histo_data_relative.GetYaxis().SetLabelOffset(0.01)
            hs_sum1_relative.SetXTitle(args.distr_name)
            histo_data_relative.SetXTitle(args.distr_name)
            hs_sum1_relative.GetXaxis().SetTitleOffset(4.) # place the title not overlapping with labels...
            histo_data_relative.GetXaxis().SetTitleOffset(4.)

        hs_sum1_relative.Draw("e2")
        histo_data_relative.Draw("e p same")

    if args.plot:
        pad1.cd()

        max_y = max([h.GetMaximum() for h in (hs_sum1, histos_data_sum)])
        min_y = max([h.GetMinimum() for h in (hs_sum1, histos_data_sum)])

        # remove the label on stack plot if ratio is there
        if args.ratio:
            hs_sum1.GetXaxis().SetLabelOffset(999)
            hs_sum1.GetXaxis().SetLabelSize(0)
            histos_data_sum.GetXaxis().SetLabelOffset(999)
            histos_data_sum.GetXaxis().SetLabelSize(0)
            #hs.GetXaxis().SetLabelOffset(999)
            #hs.GetXaxis().SetLabelSize(0)
        else:
            histos_data_sum.SetXTitle(args.distr_name)
            #hs        .SetXTitle(distr_name)
            hs_sum1   .SetXTitle(args.distr_name)

        histos_data_sum.SetMaximum(max_y * 1.1)
        if args.logy: # and args.fake_rate:
            histos_data_sum.SetMaximum(max_y * 10.)
            hs_sum1        .SetMaximum(max_y * 10.)
            histos_data_sum.SetMinimum(min_y / 10.)
            hs_sum1        .SetMinimum(min_y / 10.)
        else: # if not args.logy:
            histos_data_sum.SetMinimum(0)
            hs_sum1   .SetMinimum(0)

        histos_data_sum.SetTitle("%s %s" % (channel, sys_name))
        hs_sum1   .SetTitle("%s %s" % (channel, sys_name))

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
        if not args.fake_rate:
            leg.Draw("same")

    stack_or_ratio = ('_stack' if args.plot else '') + ('_ratio' if args.ratio else '')
    shape_chan = ('_x_' + args.shape) if args.shape else ''
    cst.SaveAs(out_dir + '_'.join((args.mc_file.replace('/', ',').split('.root')[0], args.data_file.replace('/', ',').split('.root')[0], distr_names[0], channel, sys_name)) + stack_or_ratio + shape_chan + ('_fakerate' if args.fake_rate else '') + ('_cumulative' if args.cumulative else '') + ('_logy' if args.logy else '') + ('_normalize' if args.normalize else '') + ".png")


