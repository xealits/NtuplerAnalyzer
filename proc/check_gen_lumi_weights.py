import argparse
import logging
from os.path import isfile
from jobsums.mc_norm.dtag_xsecs import usual_gen_lumi_weights


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "check event weight in the files and sumup all of them",
    epilog = """Example:\npython check_gen_lumi_weights.py lstore_outdirs/v37/test1/MC2016_Summer16_DYJetsToLL_50toInf_madgraph*/*root"""
    )

parser.add_argument('inp_files', nargs='+', help='files to check')

args = parser.parse_args()

#"lstore_outdirs/merge-sets/v37/test1/MC2016_Summer16_DYJetsToLL_10to50_amcatnlo.root"
#directory = "lstore_outdirs/merge-sets/v37/test1/"

#for dtag in usual_gen_lumi_weights.keys():
#    dtag_file = directory + dtag + '.root'
#    if isfile(dtag_file):
#        from ROOT import TFile
#        tfile = TFile(dtag_file)
#        weight_histo = tfile.Get('weight_counter')
#        gen_lumi_weight = weight_histo.GetBinContent(2)
#        print '"%20s" : %f,' % (dtag, gen_lumi_weight)

def match_dtag(fname):
    matched_dtag = None
    for dtag in usual_gen_lumi_weights.keys():
        if dtag in fname:
            matched_dtag = dtag
            break
    return matched_dtag

def get_count(fname):
    from ROOT import TFile
    tfile = TFile(fname)
    weight_histo = tfile.Get('weight_counter')
    gen_lumi_weight = weight_histo.GetBinContent(2)
    return gen_lumi_weight

dtag_counts = {}

for fname in args.inp_files:
    dtag = match_dtag(fname)
    if not dtag:
        print "no known dtag for %s" % fname

    file_count = get_count(fname)
    count = dtag_counts.get(dtag, 0)
    dtag_counts[dtag] = count + file_count

for dtag, gen_lumi_weight in dtag_counts.items():
    print '"%20s" : %f,' % (dtag, gen_lumi_weight)

