from sys import argv
import os
import logging
import threading

logging.basicConfig(level=logging.DEBUG)


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
    parser.add_argument("-s", "--range-min", type=int, default=0,    help="number of event to start processing from")
    parser.add_argument("-e", "--range-max",   type=int, default=None, help="number of event to end processing")
    parser.add_argument("-i", "--input-files", nargs='+', help='files to process (each is run in a thread)')

    args = parser.parse_args()

    logging.info('running threads on %d files' % len(args.input_files))

    from support_channel_distrs import main # ROOT stuff is loaded here

    if not os.path.exists(args.outdir + '/logs/'):
        os.makedirs(args.outdir + '/logs/')

    for input_file in args.input_files:
        #main(args.input_file, args.outdir, args.range_min, args.range_max)
        t = threading.Thread(target=main, args=(input_file, args.outdir, args.range_min, args.range_max))
        t.start()
        logging.info('started thread on %s' % input_file)


