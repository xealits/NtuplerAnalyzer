import os
import argparse
import logging
from os.path import isfile


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "just stack bunch of histograms and save together with others (standard histonames)",
    epilog = """Example:\npython stack_quick.py quick-test/v40/MC2016_Summer16_DYJetsToLL_50toInf_madgraph_all.root --chan dy_mutau_sel --sys NOMINAL --distr geom_tau_sv_sign --stack ttbar,wjets,dibosons,single_top,dy_other,dy_tautau --also data"""
    )

parser.add_argument("in_file",    help="file with histos")
parser.add_argument("out_file",   help="output filename")
parser.add_argument("--chan",  type=str, default='mu_sel', help="selection channel")
parser.add_argument("--sys",   type=str, default='NOMINAL', help="systematic")
parser.add_argument("--distr", type=str, default='geom_tau_sv_sign', help='distr')
parser.add_argument("--stack", type=str, help='processes to build a stack of')
parser.add_argument("--scale-stack", type=float, help='scale factor for processes in the stack')
parser.add_argument("--also",  type=str, help='processes to also save in output')
parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

assert isfile(args.in_file)

logging.info("import ROOT")

import ROOT
from ROOT import TFile, THStack, TLegend, TPaveText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan
import plotting_root
#from draw_overflows import DrawOverflow

tfile = TFile(args.in_file)

def std_histo_name(chan, proc, sys, distr):
    return '{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}'.format(chan=chan, proc=proc, sys=sys, distr=distr)

mc_stack_n = 0
def make_stack(tfile, stack_procs, chan, sys, distr):
    global mc_stack_n
    hs = THStack("mc_stack_%d" % mc_stack_n, "mc_stack_%d" % mc_stack_n)
    for proc in stack_procs.split(','):
        histo_name = std_histo_name(args.chan, proc, args.sys, args.distr)
        logging.debug(histo_name)

        if not tfile.Get(histo_name):
            logging.warning(histo_name + ' not in the file ' + args.in_file)
            continue

        histo = tfile.Get(histo_name)

        if args.scale_stack:
            histo.Scale(args.scale_stack)

        logging.debug(proc + ' ' + str(histo.Integral()))

        # histo styling
        col = plotting_root.nick_info[proc]['color']
        histo.SetFillColor(col);
        #histo.SetLineColor( col ); # it's needed for shapes
        histo.SetMarkerStyle(20);
        #histo.SetLineStyle(proc_ocurance);
        histo.SetMarkerColor(col);

        hs.Add(histo, "HIST")

    mc_stack_n += 1
    return hs

stacks = []
if args.stack:
   stacks = [make_stack(tfile, s, args.chan, args.sys, args.distr) for s in args.stack.split(';')]

histos = []
if args.also:
   histos = []
   for proc in args.also.split(','):
        histo_name = std_histo_name(args.chan, proc, args.sys, args.distr)
        logging.debug(histo_name)

        if not tfile.Get(histo_name):
            logging.warning(histo_name + ' (also) not in the file ' + args.in_file)
            continue

        histo = tfile.Get(histo_name)
        histos.append(histo)

if stacks + histos:
    fout = TFile(args.out_file, "RECREATE")
    fout.Write()
    fout.cd()

    for h in stacks+histos:
        h.Write()

    fout.Write()
    fout.Close()

else:
    print "nothing to output"



