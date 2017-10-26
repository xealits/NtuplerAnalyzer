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

args = parser.parse_args()

logging.info("import ROOT")

import ROOT
from ROOT import TFile, THStack, TLegend, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan

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
logging.debug("data: %f" % histo_data.Integral())

f = TFile(args.mc_file)

hs = THStack("mc_stack", "mc_stack")

nick_colour = {
"data": kWhite,
"dy": kGray,
"dy_other": kGray,
"dy_tautau": kGray+2,
"wjets": kRed+1,
"dibosons": kCyan,
"singletop": kAzure,
"s_top_eltau": kAzure,
"s_top_lj": kAzure+1,
"s_top_other": kAzure+2,

"tt_taultauhtt_other" :  kCyan-5,
"tt_jj": kGreen+4,
"tt_lj": kGreen+3,
"tt_em": kYellow-7,
"tt_{em}": kYellow-7,
"tt_ee": kAzure-9,
"tt_mm": kGreen-9,
"tt_eltau": kOrange+2,
"tt_mutau": kOrange+1,
"tt_{l\\tau-l}": kOrange+3,
"tt_other": kCyan-5,

"qcd": kViolet,
}

leg = TLegend(0.7, 0.7, 0.89, 0.89)

used_histos = []

for chan in list(f.GetListOfKeys()):
    if chan.GetName() == 'weight_counter' or chan.GetName() == 'events_counter':
        continue
    if chan.GetName() != channel and chan.GetName() != args.shape:
        continue

    processes_keys = list(chan.ReadObj().GetListOfKeys())
    sorted_pkeys = sorted(processes_keys, key=lambda pkey: pkey.GetName() if pkey.GetName() not in ('qcd') else 'z_' + pkey.GetName())
    for process in sorted_pkeys:
        nick = process.GetName()
        logging.info(nick)

        if nick == 'data':
            continue

        # if the channel with shapes for dy and wjets are given
        # don't use other channels
        if args.shape and nick in ('dy', 'wjets') and chan.GetName() != args.shape:
            continue

        for sys in list(process.ReadObj().GetListOfKeys()):
            fixed_sys_name = sys_name
            if 'TOPPT' in sys_name:
                fixed_sys_name = sys_name if 'tt' in nick else 'NOMINAL'
            if fixed_sys_name != sys.GetName():
                continue

            try:
                for histo_key in list(sys.ReadObj().GetListOfKeys()):
                    # rename for given shapes of dy and wjets
                    if args.shape and nick in ('dy', 'wjets'):
                        histo_name = '_'.join([args.shape, nick, fixed_sys_name, distr_name])
                    else:
                        histo_name = '_'.join([channel, nick, fixed_sys_name, distr_name])

                    if histo_key.GetName() != histo_name:
                        continue

                    histo = histo_key.ReadObj()
                    logging.info("%s   %s   %x = %f %f" % (histo_name, histo_key.GetName(), histo_key.GetName() == '_'.join([channel, nick, fixed_sys_name, distr_name]), histo.GetEntries(), histo.Integral()))

                    col = nick_colour[nick]
                    histo.SetFillColor( col );
                    histo.SetMarkerStyle(20);
                    histo.SetLineStyle(0);
                    histo.SetMarkerColor(col);
                    used_histos.append(histo) # hopefully root wont screw this up
                    hs.Add(histo, "HIST")
                    leg.AddEntry(histo, nick, "F")
            except Exception as e:
                print "failed", nick, e.__class__, e.__doc__

# get the sum of histograms
histos = hs.GetHists() # TList
hs_sum1 = histos[0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
hs_sum1.SetName('mc_sum1')

for h in histos[1:]:
    hs_sum1.Add(h)

hs_sum1.SetFillStyle(3004);
hs_sum1.SetFillColor(1);
hs_sum1.SetMarkerStyle(1)
hs_sum1.SetMarkerColor(0)

hs_sum2 = used_histos[0].Clone() #TH1F("sum","sum of histograms",100,-4,4);
hs_sum2.SetName('mc_sum2')

for h in used_histos[1:]:
    hs_sum2.Add(h)

hs_sum2.SetFillStyle(3004);
hs_sum2.SetFillColor(1);
hs_sum2.SetMarkerStyle(1)
hs_sum2.SetMarkerColor(0)

fout = TFile(args.mc_file.split('.root')[0] + '_%s_%s_%s.root' % (distr_name, channel, sys_name), 'RECREATE')
fout.cd()
histo_data.Write() #'data')
hs.Write()
hs_sum1.Write()
hs_sum2.Write()
leg.Write('leg')
fout.Write()
fout.Close()


