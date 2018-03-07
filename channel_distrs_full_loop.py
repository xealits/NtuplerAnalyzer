from sys import argv
import os
from os.path import isfile
import logging
import threading

#logging.basicConfig(level=logging.DEBUG)
log_common = logging.getLogger("common")


if __name__ == '__main__':
    #main(argv)
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = "Select the distrs with systematics from processed jobs.",
        #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
        )
    #def main(input_dir, dtag, outdir, range_min, range_max):

    #parser.add_argument("input_file", help="file with results from job")
    #parser.add_argument("dtag",      help="basically it's the filename with the TTree from jobs")
    parser.add_argument("outdir",    help="where to store the output")
    parser.add_argument("-l", "--log-file",  type=str, default='', help="if given the log of the process will go to this file under outdir/log/<log_file>")
    #parser.add_argument("-s", "--range-min", type=int, default=0,    help="number of event to start processing from")
    #parser.add_argument("-e", "--range-max",   type=int, default=None, help="number of event to end processing")
    parser.add_argument("-i", "--input-files", nargs='+', help='files to process (each is run in a thread)')
    parser.add_argument("--old-miniaod-jets",  action='store_true', help="the option to run on pre-v20 incorectly saved jets (initial empty, save in nominal jets)")
    parser.add_argument("--do-W-stitching",  action='store_true', help="turn ON skipping NUP events of inclusive sample")

    args = parser.parse_args()

    # configure log and common file for threads
    try:
        os.makedirs(args.outdir + '/logs/')
        
    except OSError as e:
        print e.errno, e.strerror

    if args.log_file:
        logger_file = args.outdir + '/logs/' + args.log_file.split('/')[-1].split('.root')[0] + '.log'
        hdlr = logging.FileHandler(logger_file)
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        log_common.addHandler(hdlr)
        log_common.setLevel(logging.DEBUG)

    else:
        log_common.basicConfig(level=logging.DEBUG)

    log_common.info('running threads on %d files' % len(args.input_files))


    for input_filename in args.input_files:
        #main(args.input_file, args.outdir, args.range_min, args.range_max)

        # check if output file already exists to not import ROOT
        #input_tfile = TFile(input_filename)
        #tree = input_tfile.Get('ntupler/reduced_ttree')
        #if not range_max: range_max = tree.GetEntries()

        fout_name = input_filename.split('/')[-1].split('.root')[0] + ".root" # no ranges anymore

        if isfile(args.outdir + '/' + fout_name):
            print "output file exists: %s" % (args.outdir + '/' + fout_name)
            continue

        import support_channel_distrs
        from support_channel_distrs import main # ROOT stuff is loaded here
        support_channel_distrs.OLD_MINIAOD_JETS = args.old_miniaod_jets
        support_channel_distrs.W_STITCHING = args.do_W_stitching
        t = threading.Thread(target=main, args=(input_filename, fout_name, args.outdir))
        t.start()
        log_common.info('started thread on %s' % input_filename)


