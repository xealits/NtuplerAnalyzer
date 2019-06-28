import argparse
import logging

from os.path import isfile


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "generate job queues",
    epilog = "Example:\n$ python genjobs.py v17 oldrun5 --dtags std"
    )

parser.add_argument("csv_file",    help="the csv filw with electron energy scales of EGM")
parser.add_argument("-d", "--debug",    action='store_true', help="DEBUG level of logging")

"""documentation of the csv files with scales

abs eta range,bad|gold,Et range,gainEle,literal runNumber,runNumber min,runNumber max,corrections,unc1,unc2,full_unc,unc3
absEta_0_1,bad,Et_0_20,gainEle_12,runNumber,273158,273445,1.0036,0.0002,0.0000,0.0013,0.0000

at the bottom, moved it to files with "_per_eta":
absEta_2_2.5,gold,runNumber,283946,283963,0.9974,0.0009,0.0001,0.0011,0.0000
-- these are endcaps, the other file does not have the endcap eta region

info is obtained from ploting them in R ggplot2

there are not clear parameters:
bad|gold
gainEle

each gainEle category (gainEle_12, gainEle_6, gainEle_1) has equal amount of entries
bad|gold has 2/3 bad and 1/3 gold

> summary(scales[scales$bad_good == 'gold', "gain_ele"])
 gainEle_1 gainEle_12  gainEle_6 
      1824       1824       1824 
> summary(scales[scales$bad_good == 'bad', "gain_ele"])
 gainEle_1 gainEle_12  gainEle_6 
      2432       2432       2432 

some runs do not have full coverage of et and eta:
> scales[scales$runNumber_min == 273158 & scales$et == 'Et_35_43', "abseta"]
[1] absEta_0_1 absEta_0_1 absEta_0_1
Levels: absEta_0_1 absEta_1_1.4442

???

ggplot(scales[scales$gain_ele == 'gainEle_12' & scales$abseta == 'absEta_0_1',], aes(x=et, y=correction, color=bad_good)) + geom_jitter()

-- only Et 0-20 and 100-inf have both bad and gold!
   from 20 to 100 they have different binning in Et
   and "bad" has 1 more bin then "gold"

check if every run has both "bad" and "gold"
table(scales[scales$abseta == 'absEta_0_1' & scales$et == 'Et_0_20', "runNumber_min"])

--- yes, all = 6 (3 gains and 2 bad-gold)


ggplot(scales[scales$bad_good == 'gold' & scales$abseta == 'absEta_0_1',], aes(x=et, y=correction, color=gain_ele)) + geom_jitter()

-- practically flat dep-s
   gain el 1   0.990  - 0.995
           6   1.000  - 1.005
          12   1.0025 - 1.0075

ggplot(scales[scales$bad_good == 'gold' & scales$abseta == 'absEta_1_1.4442',], aes(x=et, y=correction, color=gain_ele)) + geom_jitter()
-- similar picture, all lower by a percent

ggplot(scales[scales$gain_ele == 'gainEle_12' & scales$abseta == 'absEta_0_1',], aes(x=et, y=correction, color=bad_good)) + geom_jitter()
ggplot(scales[scales$gain_ele == 'gainEle_6'  & scales$abseta == 'absEta_0_1',], aes(x=et, y=correction, color=bad_good)) + geom_jitter()
ggplot(scales[scales$gain_ele == 'gainEle_1'  & scales$abseta == 'absEta_0_1',], aes(x=et, y=correction, color=bad_good)) + geom_jitter()

it seems gain might be a parameter of the electron...

these gains are hardware 3 signals read for each ecal deposit: amplified with these gains
https://indico.cern.ch/event/11994/contributions/84499/attachments/63981/91926/Argiro.pdf
-- some study of that

but endcaps are uniform..

https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer#How_to_apply_corrections_directl
-- add these params to ntuples and test on 1 data and tt files

scale_corr=eScaler.ScaleCorrection(unsigned int runNumber, bool isEBEle, double R9Ele, double etaSCEle, double EtEle);

http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_8_0_25/doc/html/d8/dac/GsfElectron_8h_source.html#l00335
  335       bool isEB ;        // true if particle is in ECAL Barrel

eta < 1.479

https://github.com/ECALELFS/ECALELF/blob/master/ZNtupleDumper/plugins/ZNtupleDumper.cc
 - take the R9 from the PAT electron: electron->r9()

also
REle[index] = electron.full5x5_r9();
R9Ele[index] = energy_3x3SC[index] / photon.superCluster()->rawEnergy(); //original

-- indeed I don't save these

https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEgammaShowerShape
R9: R9 = E3x3/Esupercluster (the ratio E3x3/E5x5 is also sometimes referred to as R9, e.g. in ECAL testbeam analysis). This variable makes a good separation between unconverted photons (energy not spread in tracker) and converted photons (energy spread by B-field before reaching ECAL) and is used to make this separation in the Photon object. The variable has also, in the past, been used to make a tuning corrections ("local containment") to unconverted photons in the barrel (correcting for the position of the shower with respect to the crystal matrix == "where the photon hit in the central crystal"). 

in data r9 peaks at 0.86 with a tail to low values and maximum around 1.0

data needs scale correction, which requires r9
damn smearings of MC also needs r9
"""


