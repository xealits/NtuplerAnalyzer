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
parser.add_argument("-c", "--channel",  type=str, default='mu_sel', help="selection of events")
parser.add_argument("-s", "--systematic",  type=str, default='NOMINAL', help="systematic variation")
parser.add_argument("-d", "--distr",  type=str, default='Mt_lep_met', help="recorded distribution")
parser.add_argument("-x", "--shape",  type=str, default='', help="selection for shape distributions of wjets and dy")
parser.add_argument("-p", "--plot",  action = "store_true", help="don't save the root file, plot stacked histogram in png")
parser.add_argument("-r", "--ratio", action = "store_true", help="don't save the root file, make ratio plot (in addition to stack or alone)")
parser.add_argument("-f", "--form-shapes", type=str, default='', help="plot the shapes of distributions normalized to 1")
parser.add_argument("-o", "--output-directory",  type=str, default='', help="optional output directory")

args = parser.parse_args()

logging.info("import ROOT")

import ROOT
from ROOT import TFile, THStack, TLegend, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
import plotting_root

#channel = "mu_presel"
#channel = "mu_sel"
#channel = "mu_lj"
#channel = "mu_lj_out"
channel = args.channel

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
distr_name = args.distr

fdata = TFile(args.data_file)

# only nominal in data
#histo_data = fdata.Get(channel + '/data/' + sys_name + '/' + '_'.join([channel, sys_name, 'data', distr_name]))
histo_data = fdata.Get(channel + '/data/NOMINAL/' + '_'.join([channel, 'data', 'NOMINAL', distr_name]))
histo_data.Print()


def get_histos(infile, channel, shape_channel, sys_name, distr_name):
    """get_histos(infile)

    the file contains usual structure
    <channels>/<processes>/<systematics>/<histograms>

    returns [(histo, nick)]
    """
    used_histos = [] # (histo, nick)

    chan = infile.Get(channel)
    processes_keys = list(chan.GetListOfKeys())
    sorted_pkeys = sorted(processes_keys, key=lambda pkey: pkey.GetName() if pkey.GetName() not in ('qcd') else 'z_' + pkey.GetName())
    for process in sorted_pkeys:
       nick = process.GetName()
       #logging.info(nick)

       if nick == 'data':
           continue

       fixed_sys_name = sys_name
       if 'TOPPT' in sys_name:
           fixed_sys_name = sys_name if 'tt' in nick else 'NOMINAL'

       histo_name = '_'.join([channel, nick, fixed_sys_name, distr_name])

       h_init = process.ReadObj().Get(fixed_sys_name + '/' + histo_name)

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

       used_histos.append((histo, nick)) # hopefully root wont screw this up

    return used_histos


f = TFile(args.mc_file)
used_histos = get_histos(f, channel, args.shape, sys_name, distr_name)

if args.form_shapes:
    if args.form_shapes == 'usual':
        shape_nicks = ['wjets', 'tt_mutau', 'tt_eltau', 'tt_lj']
    else:
        shape_nicks = args.form_shapes.split(',')
    used_histos = [i for i in used_histos if i[1] in shape_nicks]

# get MC stack and legend for it
hs, leg = plotting_root.stack_n_legend(used_histos)

# set data histo styles and
# add data to the legend
histo_data.SetMarkerStyle(21)
leg.AddEntry(histo_data, "data", "e1 p")

# sum of MC to get the sum of errors

# loop through TList
histos = hs.GetHists() # TList
hs_sum1 = histos[0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
hs_sum1.SetName('mc_sum1')

for h in histos[1:]:
    hs_sum1.Add(h)

hs_sum1.SetFillStyle(3004);
hs_sum1.SetFillColor(1);
hs_sum1.SetMarkerStyle(1)
hs_sum1.SetMarkerColor(0)

# loop through normal list
hs_sum2 = used_histos[0][0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
hs_sum2.SetName('mc_sum2')

for h, _ in used_histos[1:]:
    hs_sum2.Add(h)

hs_sum2.SetFillStyle(3004);
hs_sum2.SetFillColor(1);
hs_sum2.SetMarkerStyle(1)
hs_sum2.SetMarkerColor(0)

logging.info("data   = %f" % histo_data.Integral())
logging.info("mc sum = %f %f" % (hs_sum1.Integral(), hs_sum2.Integral()))

out_dir = args.output_directory + '/' if args.output_directory else './'

if not args.plot and not args.ratio:
    fout = TFile(out_dir + args.mc_file.split('.root')[0] + '_%s_%s_%s.root' % (distr_name, channel, sys_name), 'RECREATE')
    fout.cd()
    histo_data.Write() #'data')
    for h in histos:
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
    histo_data.Sumw2(ROOT.kTRUE) # to correctly save the errors after scaling
    histo_data.Scale(1/histo_data.Integral()) # but errors are not scaled here?
    for histo, nick in used_histos:
        histo.Scale(1/histo.Integral())

    # find max Y for plot
    max_y = max([h.GetMaximum() for h, _ in used_histos + [(histo_data, 'data')]])

    # plot stuff
    pad.cd()
    histo_data.SetMaximum(max_y * 1.1)
    histo_data.Draw('e1 p')
    for histo, nick in used_histos:
	#histo.SetFillColor( # it shouldn't be needed with hist drawing option
        # nope, it's needed...
        histo.SetLineColor(plotting_root.nick_colour[nick])
        histo.SetLineWidth(2)
        histo.SetFillColorAlpha(0, 0.0)
        histo.Draw("hist same")

    leg.Draw("same")

    cst.SaveAs(out_dir + '_'.join((args.mc_file.replace('/', ',').split('.root')[0], args.data_file.replace('/', ',').split('.root')[0], distr_name, channel, sys_name)) + '_shapes.png')

else:
    from ROOT import gStyle, gROOT, TCanvas, TPad
    gROOT.SetBatch()
    gStyle.SetOptStat(0)
    cst = TCanvas("cst","stacked hists",10,10,700,700)

    if args.ratio and args.plot:
        pad1 = TPad("pad1","This is pad1", 0., 0.3,  1., 1.)
        pad2 = TPad("pad2","This is pad2", 0., 0.05,  1., 0.3)
        pad1.Draw()
        pad2.Draw()
    elif args.ratio:
        pad2 = TPad("pad2","This is pad2", 0., 0.05,  1., 1.)
        pad2.Draw()
    else:
        pad1 = TPad("pad1","This is pad1", 0., 0.,  1., 1.)
        pad1.Draw()

    # plotting
    if args.ratio:
        pad2.cd()

        # calculate ratios
        histo_data_relative = histo_data.Clone()
        hs_sum1_relative = hs_sum1.Clone()
        histo_data_relative.SetName("rel_data")
        hs_sum1_relative.SetName("rel_mc")

        histo_data_relative.SetStats(False)
        hs_sum1_relative.SetStats(False)

        #histo_data_relative.GetYaxis().SetRange(0.5, 1.5)
        #histo_data_relative.GetYaxis().SetUserRange(0.5, 1.5)
        #hs_sum1_relative.GetYaxis().SetRange(0.5, 1.5)
        #hs_sum1_relative.GetYaxis().SetUserRange(0.5, 1.5)

        histo_data_relative.SetMaximum(1.5)
        histo_data_relative.SetMinimum(0.5)
        hs_sum1_relative.SetMaximum(1.5)
        hs_sum1_relative.SetMinimum(0.5)

        histo_data_relative.Divide(hs_sum1)
        hs_sum1_relative.Divide(hs_sum1)

        hs_sum1_relative.Draw("e2")
        histo_data_relative.Draw("e p same")

    if args.plot:
        pad1.cd()

        histo_data.SetMinimum(0)
        hs_sum1   .SetMinimum(0)

        histo_data.SetTitle("%s %s" % (channel, sys_name))
        hs_sum1   .SetTitle("%s %s" % (channel, sys_name))

        histo_data.SetXTitle(distr_name)
        #hs        .SetXTitle(distr_name)
        hs_sum1   .SetXTitle(distr_name)

        histo_data.Draw("e1 p")
        hs.Draw("same")
        hs_sum1.Draw("same e2")
        histo_data.Draw("same e1p")
        leg.Draw("same")

    stack_or_ratio = ('_stack' if args.plot else '') + ('_ratio' if args.ratio else '')
    shape_chan = ('_x_' + args.shape) if args.shape else ''
    cst.SaveAs(out_dir + '_'.join((args.mc_file.replace('/', ',').split('.root')[0], args.data_file.replace('/', ',').split('.root')[0], distr_name, channel, sys_name)) + stack_or_ratio + shape_chan + ".png")


