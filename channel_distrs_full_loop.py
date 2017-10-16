from sys import argv
import logging



if __name__ == '__main__':
    #main(argv)
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = "Select the distrs with systematics from processed jobs.",
        #epilog = "Example:\n$ python job_submit.py ttbarleps80_eventSelection jobing/my_runAnalysis_cfg_NEWSUBMIT.templ.py bin/ttbar-leptons-80X/analysis/dsets_testing_noHLT_TTbar.json test/tests/outdir_test_v11_ttbar_v8.40/"
        )
    #def main(input_dir, dtag, outdir, range_min, range_max):

    parser.add_argument("input_dir", help="merged-sets directory of results from jobs")
    parser.add_argument("dtag",      help="basically it's the filename with the TTree from jobs")
    parser.add_argument("outdir",    help="where to store the output")
    parser.add_argument("-s", "--range-min", type=int, default=0,    help="number of event to start processing from")
    parser.add_argument("-e", "--range-max",   type=int, default=None, help="number of event to end processing")

    args = parser.parse_args()

    from support_channel_distrs import main

    main(args.input_dir, args.dtag, args.outdir, args.range_min, args.range_max)


