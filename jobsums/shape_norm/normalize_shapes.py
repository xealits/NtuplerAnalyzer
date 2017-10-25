import argparse
import logging
logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "normalize the shapes of a selection to integrals of another one (presel to sel etc)",
    #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
    )

parser.add_argument('infile',     help="the input distributions")
#parser.add_argument('outfile',    help="where to store the output")
parser.add_argument('-s', '--selection', default='mu_presel', help="the selection with shapes (mu_presel)")
parser.add_argument('-i', '--integrals', default='mu_lj',     help="the integrals (mu_lj)")
parser.add_argument('-n', '--name',      default='mu_presel_to_mu_lj',     help="the name of new normalized selection (mu_presel_to_mu_lj)")

args = parser.parse_args()

logging.info("import ROOT")

import ROOT
from ROOT import TFile


fin  = TFile(args.infile)
fout = TFile('test1_' + args.name + '.root', "RECREATE")

for channel in list(fin.GetListOfKeys()):
    chan = channel.GetName()
    if chan == 'weight_counter' or chan == 'events_counter' or chan != args.selection:
        continue

    logging.info("shapes in %s" % chan)
    # for shapes selection loop throught all channels, sys and histo
    # get the integral of the same channels, sys, histo of normalization selection
    # scale and save to output file
    for process in list(channel.ReadObj().GetListOfKeys()):
        proc = process.GetName()
        if proc == 'data':
           continue
        logging.info(proc)
        for systematics in list(process.ReadObj().GetListOfKeys()):
            sys = systematics.GetName()

            #out_dir_name = '%s/%s/%s/' % (chan, proc, sys)
            out_dir_name = '%s/%s/%s/' % (args.name, proc, sys)
            # create these directories in output file
            if fout.Get(out_dir_name):
                #logger.debug('found ' + out_dir_name)
                out_dir = fout.Get(out_dir_name)
            else:
                #logger.debug('made  ' + out_dir_name)
                out_dir_c = fout.Get(args.name) if fout.Get(args.name) else fout.mkdir(args.name)
                #out_dir_c.cd() # this cd is needed?
                out_dir_p = out_dir_c.Get(proc) if out_dir_c.Get(proc) else out_dir_c.mkdir(proc)
                #out_dir_p.cd()
                out_dir   = out_dir_p.mkdir(sys)

            for histo_key in list(systematics.ReadObj().GetListOfKeys()):
                # so, this is the shape distr
                h = histo_key.ReadObj()

                h_out_name = histo_key.GetName() #'_'.join([args.name, sys, proc, 'Mt_lep_met'])
                # and substitute the channel name for the new one
                h_out_name = ''.join([args.name] + h_out_name.split(chan)[1:])
                h_out = h.Clone()
                h_out.SetDirectory(out_dir)
                h_out.SetName(h_out_name)
                try:
                    h_out.Sumw2(ROOT.kTRUE) # to correctly save the errors after scaling
                except Warning:
                    # it will raise Warnin for histograms which already have this -- have to filter it
                    pass

                h_norm_name = histo_key.GetName()
                # and substitute the channel name for the integrals name
                h_norm_name = ''.join([args.integrals] + h_norm_name.split(chan)[1:])
                h_normalization = fin.Get('%s/%s/%s/%s' % (args.integrals, proc, sys, h_norm_name))

                if proc == 'dy' and sys == 'NOMINAL' and histo_key.GetName() == '_'.join([chan, sys, proc, 'Mt_lep_met']):
                    logging.info("%s   %s   %s" % (histo_key.GetName(), h_norm_name, h_out_name))

                scale = 1
                init_integral = h_out.Integral()
                try:
                    #h_out.Scale(h_normalization.Integral() / h_out.Integral())
                    scale = h_normalization.Integral() / init_integral
                except ZeroDivisionError:
                    logging.error("zero division in %s %s %s %s" % (chan, proc, sys, histo_key.GetName()))
                    h_normalization.Print()
                    h_out.Print()

                h_out.Scale(scale)

                # for control
                if proc == 'dy' and sys == 'NOMINAL' and histo_key.GetName() == '_'.join([chan, sys, proc, 'Mt_lep_met']):
                    logging.info("%s %s :   %f    %f    %f" % (h_norm_name, h_out_name, init_integral, h_normalization.Integral(), h_out.Integral()))

                out_dir.cd() # it seems this is needed despite SetDirectory
                h_out.Write()

fout.Write()
fout.Close()

fin.Close()


