import argparse
import logging
import textwrap
from yaml import load
import os


logging.basicConfig(level=logging.DEBUG)

# TEMPLATE_cfg_for_crab.py  TEMPLATE_crabConf.py
with open('TEMPLATE_crabConf.py') as f:
    template_crab = f.read()

with open('TEMPLATE_cfg_for_crab.py') as f:
    template_cfg = f.read()

conf_dir = 'python/crab_cfgs/'

# why __main__ if the whole file is script?
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = "make crab config pair for given dataset",
        epilog = "Example:\n    $ python make_crab_cfgs.py '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM' v13.0 Dilep_tauID"
        )

    parser.add_argument("dataset",    help="logical name of the dataset")
    parser.add_argument("version",    help="version of the Ntupler job")
    parser.add_argument("record_scheme", help="ntupler's parameter record_scheme, sets thresholds wich events to save")
    parser.add_argument("-s", "--suffix", type=str, default='', help="suffix for dataset, used for ext MC datasets")
    parser.add_argument("-d", "--dsets-info", type=str, default='dsets_info.yaml', help="file with info on datasets: xsec, dtag, isMC, possible LumiMask")
    parser.add_argument("--without-HLT", action="store_true", default=False, help="turn off HLT in events (to support noHLT MC in 2015)")
    #parser.add_argument("-o", "--output-dir", type=str, default='python/crab_cfgs/', help="directory where config files are created for this vertsion (in version/ subdir)")

    args = parser.parse_args()

    version = args.version
    dset    = args.dataset
    suffix  = args.suffix
    record_scheme  = args.record_scheme
    #out_dir = args.output_dir

    with open(args.dsets_info) as f:
        dsets_info = load(f)

    dtag = dsets_info[dset]['dtag']
    isMC = dsets_info[dset]['isMC']
    LumiMask = dsets_info[dset]['LumiMask']
    withHLT = not args.without_HLT
    is2017rereco = 'is2017rereco' in dsets_info[dset]

    logging.debug('dtag = ' + dtag)
    logging.debug('suffix = ' + suffix)
    logging.debug('LumiMask = ' + LumiMask)
    logging.debug('isMC = ' + str(isMC))
    logging.debug('withHLT = ' + str(withHLT))
    logging.debug('is2017rereco = ' + str(is2017rereco))

    config_file = conf_dir + version + '/%s%s_cfg.py' % (dtag, suffix)

    template_crab = template_crab.format(LumiMask=LumiMask, dtag=dtag, suffix=suffix, version=version, dset=dset,
        config_file=config_file)
    template_cfg = template_cfg .format(isMC=isMC, dtag=dtag, record_scheme=record_scheme, withHLT=withHLT, is2017rereco=is2017rereco)

    if not os.path.exists(conf_dir + version):
        os.makedirs(conf_dir + version)

    crab_cfg_file = conf_dir + version + '/%s%s_crab.py' % (dtag, suffix)
    print(crab_cfg_file)
    with open(crab_cfg_file, 'w') as f:
        f.write(template_crab)

    with open(config_file, 'w') as f:
        f.write(template_cfg)

