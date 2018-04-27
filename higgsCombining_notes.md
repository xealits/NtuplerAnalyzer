# higgsCombine documentation

usually it is either not to the point or lacking or outdated

main source is the book:
https://www.gitbook.com/book/cms-hcomb/combine/details
(the twiki is outdated! use this book)

it has "Useful links" chapter:
https://cms-hcomb.gitbooks.io/combine/content/part4/index.html

there is the link to the latest tutorial:
https://indico.cern.ch/event/677948/#day-2017-11-29
-- but they actually used 1 year old slides for the tutorial,
command names changed since then and some additional options are required to run those procedures

there is link to how to cite higgsCombine -- it includes statistics information on math implemented there

and main contact is hypernews:
https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination.html







# higgs combine installation

currently used release is `CMSSW_8_1_0`

from the [book](https://cms-hcomb.gitbooks.io/combine/content/part1/):

    export SCRAM_ARCH=slc6_amd64_gcc530
    cmsrel CMSSW_8_1_0
    cd CMSSW_8_1_0/src 
    cmsenv
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd HiggsAnalysis/CombinedLimit

    cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
    git fetch origin
    git checkout v7.0.6
    scramv1 b clean; scramv1 b # always make a clean build


The script to check out CombineHarvester scripts:

    #!/bin/sh
    set -x
    set -e

    if [ -z ${CMSSW_BASE+x} ]; then
      echo "CMSSW environment not set!";
    else
      pushd $CMSSW_BASE/src
      mkdir CombineHarvester; cd CombineHarvester
      git init
      git remote add origin git@github.com:cms-analysis/CombineHarvester.git
      git config core.sparsecheckout true; echo CombineTools/ >> .git/info/sparse-checkout
      git pull origin master
      popd
    fi

just run it and then `scram b` everything

I also have link to my fit-stuff directory:

    cd HiggsAnalysis/CombinedLimit/
    $ ls ttxsec/fit-stuff -l
    lrwxrwxrwx 1 olek t3cms 92 Nov  1 19:16 ttxsec/fit-stuff -> ..../CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/fit-stuff/






# hig comb commands

go to HiggsCombine CMSSW release:

    cd CMSSW_8_1_0/src
    cmsenv
    cd HiggsAnalysis/CombinedLimit/



## join merge datacards

general command is `combineCards.py Name1=card1.txt Name2=card2.txt .... > card.txt`

    combineCards.py mu=ttxsec/fit-stuff/latest_datacard_mu.txt el=ttxsec/fit-stuff/latest_datacard_el.txt > ttxsec/fit-stuff/latest_datacard_both.txt






## NLL scan

Fit MC to Data with all nuisances, storing NLL of each variation:

    combine -M MultiDimFit ttxsec/fit-stuff/latest_datacard_mu.txt --algo grid --points 100 --rMin 0.5 --rMax 1.5 --name MuShapes


if needed

    --freezeNuisanceGroup tt_th_frag,tt_th_match,tt_th_pdf
    --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU"
    --setParameters TOPPT=1

-- `--algo grid` runs the fit over #N points (given by `--points N`) in the space of nuisance parameters.
2000 is many, usually 100 is ok to check.
Options `--rMax N` and `--rMin N` limit the range of signal strength `r` parameter.
HiggsCombine is aimed at Higgs discovery. So they always consider signal strength = 0 and call it "background-only" fit.
In cross-section measurement it is useless and usually the fit breaks somewhere.

The command will output couple files, the main one has the `--name Name` in the filename:

    root -l higgsCombineMuShapes.MultiDimFit.mH120.root
    limit->Draw("2*deltaNLL:r", "deltaNLL>0 && 2*deltaNLL<10", "L")

A horizontal line at NLL = 1 for 1 sigma variation:

    TLine* l = new TLine(0.7, 1, 1.15, 1)
    l->Draw("same")

Same for electrons:

    combine -M MultiDimFit ttxsec/fit-stuff/latest_datacard_el.txt --algo grid --points 2000 --rMax 2 --name ElShapes


(ERROR, check later) trying to get some output, some nuisances included

    combine toy-hgg-125.root -M MultiDimFit --algo grid --points 2000 --setPhysicsModelParameterRanges JES=-2,2 -m 125 --fastScan
    Invalid options: unrecognised option '--setPhysicsModelParameterRanges'


save some useful info and maybe useful RecreateNLL:

    combine -M MultiDimFit ttxsec/fit-stuff/fresh_el.txt --algo grid --points 2000 --rMax 2 --name ElShapes --forceRecreateNLL --saveSpecifiedNuis TOPPT,dy_norm
    combine -M MultiDimFit ttxsec/fit-stuff/fresh_mu.txt --algo grid --points 2000 --rMax 2 --name MuShapes --forceRecreateNLL --saveSpecifiedNuis TOPPT,dy_norm --parameters TOPPT

this crashes: `--setParameterRanges TOPPT=-1.0,2.0`


    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_mu.txt --algo grid --points 2000 --rMax 2 --name MuShapes
    root -l higgsCombineMuShapes.MultiDimFit.mH120.root
    limit->Draw("2*deltaNLL:r", "deltaNLL>0 && 2*deltaNLL<10", "L")

    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_freefrate_mu.txt --algo grid --points 2000 --rMax 2 --name MuShapesFreeFrate
    root -l higgsCombineMuShapesFreeFrate.MultiDimFit.mH120.root
    limit->Draw("2*deltaNLL:r", "deltaNLL>0 && 2*deltaNLL<10", "L")

    combine -M MultiDimFit ttxsec/fit-stuff/oldrun2_freefrate_mu.txt --algo grid --points 2000 --rMax 2 --name MuShapesOld
    root -l higgsCombineMuShapesOld.MultiDimFit.mH120.root
    limit->Draw("2*deltaNLL:r", "deltaNLL>0 && 2*deltaNLL<10", "L")


    combine -M MultiDimFit ttxsec/fit-stuff/oldrun2_freefrate_el.root --algo grid --points 1000 --rMin 0.5 --rMax 1.5 --name ElShapesOld
    root -l higgsCombineElShapesOld.MultiDimFit.mH120.root
    limit->Draw("2*deltaNLL:r", "deltaNLL>0 && 2*deltaNLL<10", "L")

    combine -M MultiDimFit ttxsec/fit-stuff/oldrun2_freefrate_el.root --algo grid --points 1000 --rMin 0.5 --rMax 1.5 --name ElShapesOldNoTTSys --freezeNuisanceGroup th_tt
    root -l higgsCombineElShapesOldNoTTSys.MultiDimFit.mH120.root
    limit->Draw("2*deltaNLL:r", "deltaNLL>0 && 2*deltaNLL<10", "L")

also useful:

--saveSpecifiedNuis




#### Breakdown of the uncertainties

From [tutorial](https://indico.cern.ch/event/677948/contributions/2776352/attachments/1550599/2468832/HComb-Tutorial-FitDiagnostics.pdf) (slide 18)
 and a hypernews question.


muons

they recommend to convert the text file of the datacard to RooFit workplace `.root` file:

    ./scripts/text2workspace.py ttxsec/fit-stuff/latest_datacard_mu.txt

it will make a file with `.txt -> .root`: `ttxsec/fit-stuff/latest_datacard_mu.root`

Then run the preliminary fit of the best signal strength:

    combine -M MultiDimFit --algo none --rMin 0.1 --rMax 2.0 ttxsec/fit-stuff/oldrun2_freefrate_mu.root --saveWorkspace -n MuBestfit

-- it saves the resulting workspace with `--saveWorkspace`, which is used in the following.

Some examples follow. But full run of uncertainty breakdown is below **in bold**.

Then freeze parameters like this:

    combine -M MultiDimFit --algo grid --points 1000 --rMin 0.1 --rMax 1.5 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" -n MuNoSysButTOPPT

Or a group of them (the group is defined in the datacard):

    combine -M MultiDimFit --algo grid --points 1000 --rMin 0.1 --rMax 1.5 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeNuisanceGroup "tau" -n MuNoTau

To fit **only stat uncertainty** freeze all nuisances and run `--fastScan`:

    combine -M MultiDimFit --algo grid --points 100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan  -n MuNoSys

The full uncertainty file has nothing frozen:

    combine -M MultiDimFit --algo grid --points 1000 --rMin 0.1 --rMax 1.5 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n MuFullUncertainty

also some runs:

    combine -M MultiDimFit --algo grid --points 2000 --rMin 0.1 --rMax 2.0 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,bSF,PU,TOPPT" -n MuNoSysButTauES
    combine -M MultiDimFit --algo grid --points 2000 --rMin 0.1 --rMax 2.0 -n MuNoSysButTauFakes higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU,TOPPT"
    combine -M MultiDimFit --algo grid --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeNuisanceGroup "tau" -n MuNoTau



**Full breakdown for muons**:

    ./scripts/text2workspace.py ttxsec/fit-stuff/latest_datacard_mu.txt
    combine -M MultiDimFit --algo none               --rMin 0.7 --rMax 1.3 ttxsec/fit-stuff/latest_datacard_mu.root --saveWorkspace -n MuBestfit

    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan  -n MuNoSys
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" --freezeNuisanceGroup "th_tt" -n MuNoSysButTOPPT
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff" -n MuNoTau
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV" -n MuNoLumi
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n MuFullUncertainty

no the weights if needed --freezeNuisanceGroup tt_th_pdf,tt_th_frag,tt_th_match

To plot the breakdown with the script from Harvester repository:

    plot1DScan.py higgsCombineMuFullUncertainty.MultiDimFit.mH120.root --others 'higgsCombineMuNoLumi.MultiDimFit.mH120.root:lumi:4' 'higgsCombineMuNoSys.MultiDimFit.mH120.root:stat:2' --breakdown lumi,syst,stat

It outputs 2 files, renaming them:

    mv scan.png mu_uncert_break_scan.png
    mv scan.pdf mu_uncert_break_scan.pdf

Old plot, when proper stat scan was not available:

    plot1DScan.py higgsCombineMuFullUncertainty.MultiDimFit.mH120.root --others 'higgsCombineMuNoLumi.MultiDimFit.mH120.root:lumi:4' 'higgsCombineMuNoSysButTOPPT.MultiDimFit.mH120.root:stat:2' --breakdown lumi,syst,stat



muons for v25v26:

    ./scripts/text2workspace.py ttxsec/fit-stuff/latest_datacard_mu.txt
    combine -M MultiDimFit --algo none               --rMin 0.7 --rMax 1.3 ttxsec/fit-stuff/latest_datacard_mu.root --saveWorkspace -n MuBestfit --freezeParameters AlphaS

    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan  -n MuNoSys
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "AlphaS,lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" --freezeNuisanceGroup "th_tt" -n MuNoSysButTOPPT
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "AlphaS,tauID_eff" -n MuNoTau
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "AlphaS,lumi_13TeV" -n MuNoLumi
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root -m 120 --snapshotName MultiDimFit --freezeParameters all --fastScan -n MuNoSys
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n MuFullUncertainty --freezeParameters AlphaS

and

    plot1DScan.py higgsCombineMuFullUncertainty.MultiDimFit.mH120.root --others 'higgsCombineMuNoLumi.MultiDimFit.mH120.root:lumi:4' 'higgsCombineMuNoSys.MultiDimFit.mH120.root:stat:2' --breakdown lumi,syst,stat




**electrons**

    ./scripts/text2workspace.py ttxsec/fit-stuff/latest_datacard_el.txt -m 120
    combine -M MultiDimFit --algo none --rMin 0.1 --rMax 2.0 ttxsec/fit-stuff/latest_datacard_el.root -m 120 --saveWorkspace -n ElBestfit

electrons seriously break at about 0.5 signal strength -- fit scan and uncertainty breakdown doesnt work

    ./scripts/text2workspace.py ttxsec/fit-stuff/latest_datacard_el.txt
    combine -M MultiDimFit --algo none              --rMin 0.7 --rMax 1.3 ttxsec/fit-stuff/latest_datacard_el.root -m 120 --saveWorkspace -n ElBestfit
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.7 --rMax 1.3 higgsCombineElBestfit.MultiDimFit.mH120.root -m 120 --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" -n ElNoSysButTOPPT
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.7 --rMax 1.3 higgsCombineElBestfit.MultiDimFit.mH120.root -m 120 --snapshotName MultiDimFit --freezeNuisanceGroup "tau" -n ElNoTau
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.7 --rMax 1.3 higgsCombineElBestfit.MultiDimFit.mH120.root -m 120 --snapshotName MultiDimFit --freezeParameters "lumi_13TeV" -n ElNoLumi
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.7 --rMax 1.3 higgsCombineElBestfit.MultiDimFit.mH120.root -m 120 --snapshotName MultiDimFit --freezeParameters all --fastScan -n ElNoSys
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.7 --rMax 1.3 higgsCombineElBestfit.MultiDimFit.mH120.root -m 120 --snapshotName MultiDimFit -n ElFullUncertainty



And plots:

    plot1DScan.py higgsCombineElFullUncertainty.MultiDimFit.mH120.root

    plot1DScan.py higgsCombineElFullUncertainty.MultiDimFit.mH120.root --others 'higgsCombineElNoLumi.MultiDimFit.mH120.root:lumi:4' 'higgsCombineElNoSysButTOPPT.MultiDimFit.mH120.root:stat:2' --breakdown lumi,syst,stat
    plot1DScan.py higgsCombineElFullUncertainty.MultiDimFit.mH120.root --others 'higgsCombineElNoLumi.MultiDimFit.mH120.root:lumi:4' 'higgsCombineElNoSys.MultiDimFit.mH120.root:stat:2' --breakdown lumi,syst,stat

    mv scan.png el_uncert_break_scan.png
    mv scan.pdf el_uncert_break_scan.pdf



Some old notes on fitting processing oldrun2 fine

    ./scripts/text2workspace.py ttxsec/fit-stuff/oldrun2_freefrate_el_fine.txt
    combine -M MultiDimFit --algo none               --rMin 0.5 --rMax 1.5 ttxsec/fit-stuff/oldrun2_freefrate_el_fine.root --saveWorkspace -n ElBestfit

    combine -M MultiDimFit --algo grid --points 1000 --rMin 0.5 --rMax 1.5 higgsCombineElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" --freezeNuisanceGroup "th_tt" -n ElNoSysButTOPPT
    combine -M MultiDimFit --algo grid --points 1000 --rMin 0.5 --rMax 1.5 higgsCombineElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff" -n ElNoTau
    combine -M MultiDimFit --algo grid --points 1000 --rMin 0.5 --rMax 1.5 higgsCombineElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n ElFullUncertainty





**both**

    ./scripts/text2workspace.py ttxsec/fit-stuff/latest_datacard_both.txt
    combine -M MultiDimFit --algo none               --rMin 0.7 --rMax 1.3 ttxsec/fit-stuff/latest_datacard_both.root --saveWorkspace -n BothBestfit
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" --freezeNuisanceGroup "th_tt" -n BothNoSysButTOPPT
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan -n BothNoSys
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff" -n BothNoTau
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV" -n BothNoLumi
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n BothFullUncertainty

or in v25v26:

    ./scripts/text2workspace.py ttxsec/fit-stuff/latest_datacard_both.txt
    combine -M MultiDimFit --algo none               --rMin 0.7 --rMax 1.3 ttxsec/fit-stuff/latest_datacard_both.root --saveWorkspace -n BothBestfit --freezeParameters AlphaS
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "AlphaS,lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" --freezeNuisanceGroup "th_tt" -n BothNoSysButTOPPT
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan -n BothNoSys
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "AlphaS,tauID_eff" -n BothNoTau
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "AlphaS,lumi_13TeV" -n BothNoLumi
    combine -M MultiDimFit --algo grid --points  100 --rMin 0.7 --rMax 1.3 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n BothFullUncertainty --freezeParameters AlphaS

plotting and saving:

    plot1DScan.py higgsCombineBothFullUncertainty.MultiDimFit.mH120.root --others 'higgsCombineBothNoLumi.MultiDimFit.mH120.root:lumi:4' 'higgsCombineBothNoSys.MultiDimFit.mH120.root:stat:2' --breakdown lumi,syst,stat

    mv scan.png both_uncert_break_scan.png
    mv scan.pdf both_uncert_break_scan.pdf







### Expected NLL

Its basically running the same fit on generated toy datasets instead of data. Or on Asimov dataset.
To run on toys instread of data add option `-t Ntoys` or `-t -1` to run on Asimov dataset.

Examples:

    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_mu.txt           --algo grid --points 2000 --rMax 2 --name ExpectedMu -t -1 --expectSignal 1
    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_mu_inclusive.txt --algo grid --points 2000 --rMax 2 --name ExpectedInclusiveMu -t -1 --expectSignal 1

    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_freefrate_mu.txt           --algo grid --points 2000 --rMax 2 --name FFrateExpectedMu -t -1 --expectSignal 1
    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_freefrate_mu_inclusive.txt --algo grid --points 2000 --rMax 2 --name FFrateExpectedInclusiveMu -t -1 --expectSignal 1

freeze tau uncert:

    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_mu.txt --algo grid --points 2000 --rMax 2 --name ExpectedMu -t -1 --expectSignal 1 --freezeNuisanceGroup=tau

to freeze all systs:
(none run must find only the best fit point)

    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_freefrate_mu.txt --algo none --rMax 2 --name ExpectedMuBest -t -1 --expectSignal 1 --saveWorkspace
    combine -M MultiDimFit higgsCombineExpectedMuBest.MultiDimFit.mH120.root --algo grid --points 200 --rMax 2 --name ExpectedMuStats -t -1 --expectSignal 1 --snapshotName MultiDimFit --freezeNuisances all

-- old syntax doesnt work

still nothing works


Expected uncertainties

I keep 1 week systematic on to emulate statistical uncertainty.

    ./scripts/text2workspace.py ttxsec/fit-stuff/optimized5_freefrate_mu.txt
    combine -M MultiDimFit --algo none -t -1 --expectSignal 1 --rMin 0.1 --rMax 2.0 ttxsec/fit-stuff/optimized5_freefrate_mu.root --saveWorkspace -n ExpectedMuBestfit
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.1 --rMax 2.0 higgsCombineExpectedMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" -n ExpectedMuNoSysButTOPPT
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.1 --rMax 2.0 higgsCombineExpectedMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeNuisanceGroup "tau" -n ExpectedMuNoTau
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.1 --rMax 2.0 higgsCombineExpectedMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n ExpectedMuFullUncertainty


    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeNuisanceGroup "tau" -n ExpectedMuNoTau


    ./scripts/text2workspace.py ttxsec/fit-stuff/oldrun2_freefrate_mu.txt
    combine -M MultiDimFit --algo none -t -1 --expectSignal 1               --rMin 0.5 --rMax 1.5 ttxsec/fit-stuff/oldrun2_freefrate_mu.root --saveWorkspace -n ExpectedMuBestfit
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" --freezeNuisanceGroup "th_tt" -n ExpectedMuNoSysButTOPPT
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff" -n ExpectedMuNoTau
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan -n ExpectedMuNoSys
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n ExpectedMuFullUncertainty




electrons:

    ./scripts/text2workspace.py ttxsec/fit-stuff/optimized5_freefrate_el.txt
    combine -M MultiDimFit --algo none -t -1 --expectSignal 1 --rMin 0.1 --rMax 2.0 ttxsec/fit-stuff/optimized5_freefrate_el.root --saveWorkspace -n ExpectedElBestfit
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.1 --rMax 2.0 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" -n ExpectedElNoSysButTOPPT
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.1 --rMax 2.0 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeNuisanceGroup "tau" -n ExpectedElNoTau
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.1 --rMax 2.0 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n ExpectedElFullUncertainty


    ./scripts/text2workspace.py ttxsec/fit-stuff/oldrun2_freefrate_el_fine.txt
    combine -M MultiDimFit --algo none -t -1 --expectSignal 1 --rMin 0.5 --rMax 1.5 ttxsec/fit-stuff/oldrun2_freefrate_el_fine.root --saveWorkspace -n ExpectedElBestfit
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" --freezeNuisanceGroup "th_tt" -n ExpectedElNoSysButTOPPT
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff" -n ExpectedElNoTau
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan -n ExpectedElNoSys
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n ExpectedElFullUncertainty


    ./scripts/text2workspace.py ttxsec/fit-stuff/oldrun2_freefrate_both.txt
    combine -M MultiDimFit --algo none -t -1 --expectSignal 1               --rMin 0.5 --rMax 1.5 ttxsec/fit-stuff/oldrun2_freefrate_both.root --saveWorkspace -n ExpectedBothBestfit
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" --freezeNuisanceGroup "th_tt" -n ExpectedBothNoSysButTOPPT
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff" -n ExpectedBothNoTau
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan -n ExpectedBothNoSys
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 2000 --rMin 0.5 --rMax 1.5 higgsCombineExpectedBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n ExpectedBothFullUncertainty


new

    ./scripts/text2workspace.py ttxsec/fit-stuff/latest_datacard_el.txt

el

expected

    combine -M MultiDimFit --algo none -t -1 --expectSignal 1 --rMin 0.5 --rMax 1.5 ttxsec/fit-stuff/latest_datacard_el.root --saveWorkspace -n ExpectedElBestfit
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 100 --rMin 0.5 --rMax 1.5 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,tauID_eff,tau_fakes,dy_norm,wjets_norm,qcd_norm,JES,JER,TauES,bSF,PU" --freezeNuisanceGroup "th_tt" -n ExpectedElNoSysButTOPPT
    
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 100 --rMin 0.5 --rMax 1.5 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff,tau_fakes" -n ExpectedElNoTau
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 100 --rMin 0.5 --rMax 1.5 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan -n ExpectedElNoSys
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 100 --rMin 0.5 --rMax 1.5 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n ExpectedElFullUncertainty

to separate lumi:

    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 100 --rMin 0.5 --rMax 1.5 higgsCombineExpectedElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV" -n ExpectedElNoLumi

data fit

    combine -M MultiDimFit --algo none --rMin 0.5 --rMax 1.5 ttxsec/fit-stuff/latest_datacard_el.root --saveWorkspace -n ElBestfit
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.5 --rMax 1.5 higgsCombineElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff,tau_fakes" -n ElNoTau
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.5 --rMax 1.5 higgsCombineElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan -n ElNoSys
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.5 --rMax 1.5 higgsCombineElBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n ElFullUncertainty

mu

    combine -M MultiDimFit --algo none -t -1 --expectSignal 1 --rMin 0.5 --rMax 1.5 ttxsec/fit-stuff/latest_datacard_mu.root --saveWorkspace -n ExpectedMuBestfit --freezeParameters AlphaS
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 100 --rMin 0.5 --rMax 1.5 higgsCombineExpectedMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff,tau_fakes,AlphaS" -n ExpectedMuNoTau
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 100 --rMin 0.5 --rMax 1.5 higgsCombineExpectedMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan -n ExpectedMuNoSys
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 100 --rMin 0.5 --rMax 1.5 higgsCombineExpectedMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n ExpectedMuFullUncertainty --freezeParameters AlphaS

    combine -M MultiDimFit --algo none --rMin 0.5 --rMax 1.5 ttxsec/fit-stuff/latest_datacard_mu.root --saveWorkspace -n MuBestfit --freezeParameters AlphaS
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.5 --rMax 1.5 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff,tau_fakes,AlphaS" -n MuNoTau
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.5 --rMax 1.5 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan -n MuNoSys
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.5 --rMax 1.5 higgsCombineMuBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n MuFullUncertainty --freezeParameters AlphaS

both

    combine -M MultiDimFit --algo none -t -1 --expectSignal 1 --rMin 0.5 --rMax 1.5 ttxsec/fit-stuff/latest_datacard_both.root --saveWorkspace -n ExpectedBothBestfit --freezeParameters AlphaS
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 100 --rMin 0.5 --rMax 1.5 higgsCombineExpectedBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff,tau_fakes,AlphaS" -n ExpectedBothNoTau
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 100 --rMin 0.5 --rMax 1.5 higgsCombineExpectedBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan -n ExpectedBothNoSys
    combine -M MultiDimFit --algo grid -t -1 --expectSignal 1 --points 100 --rMin 0.5 --rMax 1.5 higgsCombineExpectedBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n ExpectedBothFullUncertainty --freezeParameters AlphaS

    combine -M MultiDimFit --algo none --rMin 0.5 --rMax 1.5 ttxsec/fit-stuff/latest_datacard_both.root --saveWorkspace -n BothBestfit --freezeParameters AlphaS
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.5 --rMax 1.5 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff,tau_fakes,AlphaS" -n BothNoTau
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.5 --rMax 1.5 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan -n BothNoSys
    combine -M MultiDimFit --algo grid --points 100 --rMin 0.5 --rMax 1.5 higgsCombineBothBestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit -n BothFullUncertainty --freezeParameters AlphaS



plot both channels:

    python plot_uncertainty_scans.py





### Bias band plot

Could not get Combine/Harvester to make this yet.

therefore, separate script with hard-codes

    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_freefrate_el.txt           --algo grid --points 2000 --rMax 2 --name ExpectedElForSig09 -t -1 --expectSignal 0.9 &
    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_freefrate_el.txt           --algo grid --points 2000 --rMax 2 --name ExpectedElForSig10 -t -1 --expectSignal 1.0 &
    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_freefrate_el.txt           --algo grid --points 2000 --rMax 2 --name ExpectedElForSig11 -t -1 --expectSignal 1.1 &

    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_freefrate_mu.txt           --algo grid --points 2000 --rMax 2 --name ExpectedMuForSig09 -t -1 --expectSignal 0.9 &
    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_freefrate_mu.txt           --algo grid --points 2000 --rMax 2 --name ExpectedMuForSig10 -t -1 --expectSignal 1.0 &
    combine -M MultiDimFit ttxsec/fit-stuff/optimized5_freefrate_mu.txt           --algo grid --points 2000 --rMax 2 --name ExpectedMuForSig11 -t -1 --expectSignal 1.1 &


    ((TTree*) _file0->Get("limit"))->Scan("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL<0.01")
    ((TTree*) _file1->Get("limit"))->Scan("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL<0.01")
    ((TTree*) _file2->Get("limit"))->Scan("2*deltaNLL:r", "2*deltaNLL>0 && 2*deltaNLL<0.01")

    ((TTree*) _file0->Get("limit"))->Scan("2*deltaNLL:r", "2*deltaNLL>0.9 && 2*deltaNLL<1.1")
    ((TTree*) _file1->Get("limit"))->Scan("2*deltaNLL:r", "2*deltaNLL>0.9 && 2*deltaNLL<1.1")
    ((TTree*) _file2->Get("limit"))->Scan("2*deltaNLL:r", "2*deltaNLL>0.9 && 2*deltaNLL<1.1")

    ((TTree*) _file0->Get("limit"))->Scan("2*deltaNLL:r", "2*deltaNLL>3.9 && 2*deltaNLL<4.1")
    ((TTree*) _file1->Get("limit"))->Scan("2*deltaNLL:r", "2*deltaNLL>3.9 && 2*deltaNLL<4.1")
    ((TTree*) _file2->Get("limit"))->Scan("2*deltaNLL:r", "2*deltaNLL>3.9 && 2*deltaNLL<4.1")








## post-fit pulls in fit

Some examples:

    combine -M FitDiagnostics ttxsec/fit-stuff/tt_el_2.txt --name El
    combine -M FitDiagnostics ttxsec/fit-stuff/fresh_mu.txt --name MuShape
    combine -M FitDiagnostics ttxsec/fit-stuff/fresh_el.txt --name ElShape

    combine -M FitDiagnostics ttxsec/fit-stuff/oldrun2_freefrate_mu.txt --name MuShape  --noMCbonly 0  &
    combine -M FitDiagnostics ttxsec/fit-stuff/oldrun2_freefrate_el_fine.txt --name ElShape  --noMCbonly 0  &

    python test/diffNuisances.py fitDiagnostics.root -g pulls_of_nuisances.root -o

-- the scripts break on background-obly fits here. The opption `--noMCbonly 0/1` doesnt work.

Therefore Icorrected diffNuisances to output only s+b -- option "-o".
No need now with --noMCbonly -- no, it seems to not work.

    combine -M FitDiagnostics ttxsec/fit-stuff/latest_datacard_mu.txt --name MuShape  --noMCbonly 0  &
    combine -M FitDiagnostics ttxsec/fit-stuff/latest_datacard_mu.txt --name MuShape  --noMCbonly 0 --freezeParameters AlphaS &
    combine -M FitDiagnostics ttxsec/fit-stuff/latest_datacard_both.txt --name BothShape  --noMCbonly 0 --freezeParameters AlphaS &
    combine -M FitDiagnostics ttxsec/fit-stuff/latest_datacard_el.txt --name ElShape  --noMCbonly 0  &

Produce nuisances with custom script which skips background-only step:

    python ttxsec/fit-stuff/diffNuisances.py fitDiagnosticsMuShape.root -g pulls_of_nuisances_mu_shape.root -o
    python ttxsec/fit-stuff/diffNuisances.py fitDiagnosticsBothShape.root -g pulls_of_nuisances_both_shape.root -o

Plot result:

    root -l pulls_of_nuisances_mu_shape.root
    gStyle->SetOptStat(0)
    nuisancs->Draw()

    TPaveText *pt = new TPaveText(.5, 1.5, 3.5, 2.5)
    pt->AddText("muon-tau channel")
    pt->Draw("same")

    root -l pulls_of_nuisances_both_shape.root
    gStyle->SetOptStat(0)
    nuisancs->Draw()

    TPaveText *pt = new TPaveText(.5, 1.5, 3.5, 2.5)
    pt->AddText("both channels")
    pt->Draw("same")

And electrons:

    python ttxsec/fit-stuff/diffNuisances.py fitDiagnosticsElShape.root -g pulls_of_nuisances_el_shape.root -o
    root -l pulls_of_nuisances_el_shape.root

    TPaveText *pt = new TPaveText(.5, 1.5, 3.5, 2.5)
    pt->AddText("electron-tau channel")
    pt->Draw("same")


At some point I got this error:

> now it crashes on electrons (!!!! worked yesterday fully all was done with electrons!!) with:
>    File "ttxsec/fit-stuff/diffNuisances.py", line 50, in <module>
>        if fit_b == None or fit_b.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the background fit 'fit_b'" % args[0]
>    RuntimeError: File fitDiagnostics.root does not contain the output of the background fit 'fit_b'

but fixed it later




## post-fit distributions

TODO: need to finish this, add a script with propper plotting of the post-fit distrs?

Examples:

    combine -M MaxLikelihoodFit ttxsec/play_shapes_tt_mu.txt --saveShapes --saveWithUncertainties
    combine -M FitDiagnostics ttxsec/fit-stuff/tt_mu_2_spec.txt --saveShapes --saveWithUncertainties --name Mu
    combine -M FitDiagnostics ttxsec/fit-stuff/fresh_mu.txt --saveShapes --saveWithUncertainties --name MuShapes
    combine -M FitDiagnostics ttxsec/fit-stuff/latest_datacard_mu.txt  --saveShapes --saveWithUncertainties --name MuShapes
    combine -M FitDiagnostics ttxsec/fit-stuff/latest_datacard_mu.txt  --saveShapes --saveWithUncertainties --name MuShapes --freezeParameters AlphaS
    combine -M FitDiagnostics ttxsec/fit-stuff/latest_datacard_el.txt  --saveShapes --saveWithUncertainties --name ElShapes
    ls fitDiagnosticsMu.root

-- saves post-fit uncertainty, and then additional stack-plot should put the distr-s together

also notice:
 WARNING --  From combine v7, method MaxLikelihoodFit has been renamed to FitDiagnostics

stacking script:

    python stack_postfit.py fitDiagnosticsMuShapes.root ctr_old_mu_sel_lj --mu




## covariance matrices

they are done in simple

    combine -M FitDiagnostics ttxsec/fit-stuff/latest_datacard_mu.txt --plots

-- the png files are saved for fit with b or b+s

more plots with:

    combine -M FitDiagnostics ttxsec/fit-stuff/latest_datacard_mu.txt --plots --saveShapes


Overall:

    combine -M FitDiagnostics ttxsec/fit-stuff/latest_datacard_mu.txt --plots --saveShapes --saveWithUncertainties --name MuShapes





## datacard to "workspace"

required here and there by stuff and also is recommended step before running all commands

    ./scripts/text2workspace.py ttxsec/tt_el_2.txt
    ./scripts/text2workspace.py ttxsec/fit-stuff/fresh_mu.txt

    ./scripts/text2workspace.py ttxsec/fit-stuff/optimized5_mu.txt
    ./scripts/text2workspace.py ttxsec/fit-stuff/optimized5_freefrate_mu.txt



## toys

THIS DIDNT WORK, use the section "expected NLL".

    combine -M FitDiagnostics ttxsec/tt_el_2.txt --name Data
    combine -M FitDiagnostics ttxsec/tt_el_2.txt -t 20000 --toysFrequentist --noErrors --minos none --name Toys --expectSignal 1

-- (without --expectSignal 1 it seems this generates toys with signal = 0)

    cd test
    root -l
    .L plotParametersFromToys.C+

    plotParametersFromToys("../fitDiagnosticsToys.root","../fitDiagnosticsData.root", "../tt_el_2.root", "r>0 && r<2")

-- weirdness

    --expectSignal arg (=0)       If set to non-zero, generate *signal* 
                                  toys instead of background ones, with 
                                  the specified signal strength.


returns pdf `tree_fit_sb.pdf`


## toys impacts

    ./scripts/text2workspace.py ttxsec/fit-stuff/datacard_AN_v2_mu

Asimov toy

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts  -t -1 --expectSignal 1 -d ttxsec/fit-stuff/datacard_AN_v2_mu.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts  -t -1 --expectSignal 1 -d ttxsec/fit-stuff/datacard_AN_v2_mu.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts  -t -1 --expectSignal 1 -d ttxsec/fit-stuff/datacard_AN_v2_mu.root -m 125 -o latest_mu_toys_asimov_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_mu_toys_asimov_impacts.json -o latest_mu_toys_asimov_impacts

1000 toys

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts  -t 100 --expectSignal 1 -d ttxsec/fit-stuff/datacard_AN_v2_mu.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts  -t 100 --expectSignal 1 -d ttxsec/fit-stuff/datacard_AN_v2_mu.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts  -t 100 --expectSignal 1 -d ttxsec/fit-stuff/datacard_AN_v2_mu.root -m 125 -o latest_mu_toys_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_mu_toys_impacts.json -o latest_mu_toys_impacts


mu

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts --name MuImpactsToys  -t -1 --expectSignal 1 --freezeParameters AlphaS -d ttxsec/fit-stuff/latest_datacard_mu.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts --name MuImpactsToys  -t -1 --expectSignal 1 --freezeParameters AlphaS -d ttxsec/fit-stuff/latest_datacard_mu.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts --name MuImpactsToys  -t -1 --expectSignal 1 --freezeParameters AlphaS -d ttxsec/fit-stuff/latest_datacard_mu.root -m 125 -o latest_mu_toys_asimov_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_mu_toys_asimov_impacts.json -o latest_mu_toys_asimov_impacts

el

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts --name ElImpactsToys  -t -1 --expectSignal 1 -d ttxsec/fit-stuff/latest_datacard_el.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts --name ElImpactsToys  -t -1 --expectSignal 1 -d ttxsec/fit-stuff/latest_datacard_el.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts --name ElImpactsToys  -t -1 --expectSignal 1 -d ttxsec/fit-stuff/latest_datacard_el.root -m 125 -o latest_el_toys_asimov_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_el_toys_asimov_impacts.json -o latest_el_toys_asimov_impacts

both

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts --name BothImpactsToys  -t -1 --expectSignal 1 --freezeParameters AlphaS -d ttxsec/fit-stuff/latest_datacard_both.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts --name BothImpactsToys  -t -1 --expectSignal 1 --freezeParameters AlphaS -d ttxsec/fit-stuff/latest_datacard_both.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts --name BothImpactsToys  -t -1 --expectSignal 1 --freezeParameters AlphaS -d ttxsec/fit-stuff/latest_datacard_both.root -m 125 -o latest_both_toys_asimov_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_both_toys_asimov_impacts.json -o latest_both_toys_asimov_impacts






## impacts

Uses Harvester script to produce impacts diagram

First produce `.root` workspace file:

    ./scripts/text2workspace.py ttxsec/fit-stuff/latest_datacard_mu.txt
    ./scripts/text2workspace.py ttxsec/fit-stuff/latest_datacard_el.txt
    ./scripts/text2workspace.py ttxsec/fit-stuff/latest_datacard_both.txt

electrons (I had some issues with these, but fixed them somehow):

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root -m 125 -o latest_el_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_el_impacts.json -o latest_el_impacts

if needed `--freezeNuisanceGroup tt_th_frag,tt_th_match`


../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root --freezeNuisanceGroup tt_th_pdf --freezeParameters AlphaS -m 125 --robustFit 1 --doInitialFit
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root --freezeNuisanceGroup tt_th_pdf --freezeParameters AlphaS -m 125 --robustFit 1 --doFits --parallel 5
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root --freezeNuisanceGroup tt_th_pdf --freezeParameters AlphaS -m 125 -o latest_el_impacts.json
../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_el_impacts.json -o latest_el_impacts


../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeNuisanceGroup tt_th_pdf --freezeParameters AlphaS -m 125 --robustFit 1 --doInitialFit
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeNuisanceGroup tt_th_pdf --freezeParameters AlphaS -m 125 --robustFit 1 --doFits --parallel 5
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeNuisanceGroup tt_th_pdf --freezeParameters AlphaS -m 125 -o latest_mu_impacts.json
../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_mu_impacts.json -o latest_mu_impacts


with PDFs

../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root --freezeParameters AlphaS -m 125 --robustFit 1 --doInitialFit         --name ElImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root --freezeParameters AlphaS -m 125 --robustFit 1 --doFits --parallel 5  --name ElImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root --freezeParameters AlphaS -m 125 -o latest_el_impacts.json            --name ElImpacts
../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_el_impacts.json -o latest_el_impacts


../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root -m 125 --robustFit 1 --doInitialFit         --name ElImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root -m 125 --robustFit 1 --doFits --parallel 5  --name ElImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root -m 125 -o latest_el_impacts.json            --name ElImpacts
../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_el_impacts.json -o latest_el_impacts


../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeParameters AlphaS -m 125 --robustFit 1 --doInitialFit         --name MuImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeParameters AlphaS -m 125 --robustFit 1 --doFits --parallel 5  --name MuImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeParameters AlphaS -m 125 -o latest_mu_impacts.json            --name MuImpacts
../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_mu_impacts.json -o latest_mu_impacts

both

../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_both.root --freezeParameters AlphaS -m 125 --robustFit 1 --doInitialFit         --name BothImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_both.root --freezeParameters AlphaS -m 125 --robustFit 1 --doFits --parallel 5  --name BothImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_both.root --freezeParameters AlphaS -m 125 -o latest_both_impacts.json          --name BothImpacts
../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_both_impacts.json -o latest_both_impacts




mu without updowns tt_updowns and pdfs:

../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 --robustFit 1 --doInitialFit         --name MuImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 --robustFit 1 --doFits --parallel 5  --name MuImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 -o latest_mu_impacts.json            --name MuImpacts
../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_mu_impacts.json -o latest_mu_impacts

same el:

../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 --robustFit 1 --doInitialFit         --name MuImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 --robustFit 1 --doFits --parallel 5  --name MuImpacts
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_el.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 -o latest_el_impacts.json            --name MuImpacts
../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_el_impacts.json -o latest_el_impacts

mu toys

../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -t -1 --expectSignal 1 -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 --robustFit 1 --doInitialFit         --name MuImpactsToys
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -t -1 --expectSignal 1 -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 --robustFit 1 --doFits --parallel 5  --name MuImpactsToys
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -t -1 --expectSignal 1 -d ttxsec/fit-stuff/latest_datacard_mu.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 -o latest_mu_impacts_toys.json            --name MuImpactsToys
../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_mu_impacts_toys.json -o latest_mu_impacts_toys

el toys

../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -t -1 --expectSignal 1 -d ttxsec/fit-stuff/latest_datacard_el.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 --robustFit 1 --doInitialFit         --name ElImpactsToys
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -t -1 --expectSignal 1 -d ttxsec/fit-stuff/latest_datacard_el.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 --robustFit 1 --doFits --parallel 5  --name ElImpactsToys
../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -t -1 --expectSignal 1 -d ttxsec/fit-stuff/latest_datacard_el.root --freezeParameters AlphaS --freezeNuisanceGroup tt_updowns,tt_th_pdf,tt_th_match,tt_th_frag -m 125 -o latest_el_impacts_toys.json            --name ElImpactsToys
../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_el_impacts_toys.json -o latest_el_impacts_toys



muons:

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_mu.root -m 125 -o latest_mu_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_mu_impacts.json -o latest_mu_impacts

both

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_both.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_both.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/latest_datacard_both.root -m 125 -o latest_both_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_both_impacts.json -o latest_both_impacts



Examples from the past follow

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d tt_el_12-1.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d tt_el_12-1.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d tt_el_12-1.root -m 125 -o el_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts.py -i el_impacts.json -o el_impacts


    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d tt_mu_12-1.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d tt_mu_12-1.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d tt_mu_12-1.root -m 125 -o mu_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts.py -i mu_impacts.json -o mu_impacts

it creates pdf


    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_mu.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_mu.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_mu.root -m 125 -o mu_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts.py -i mu_impacts.json -o fixed_mu_impacts


    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_el.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_el.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_el.root -m 125 -o el_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts.py -i el_impacts.json -o fixed_el_impacts

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_el_mixed.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_el_mixed.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_el_mixed.root -m 125 -o el_impacts.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts.py -i el_impacts.json -o fixed_el_impacts

muon fit investigation:

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_mu_just_fixed.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_mu_just_fixed.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_mu_just_fixed.root -m 125 -o mu_impacts_just.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i mu_impacts_just.json -o fixed_mu_impacts_just


    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_el_just.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_el_just.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/fresh_el_just.root -m 125 -o el_impacts_just.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i el_impacts_just.json -o fixed_el_impacts_just

    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i el_impacts_just.json -o fixed_el_impacts_just



    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_mu.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_mu.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_mu.root -m 125 -o mu_impacts_opt5.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i mu_impacts_opt5.json -o fixed_mu_impacts_opt5

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_mu_inclusive.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_mu_inclusive.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_mu_inclusive.root -m 125 -o mu_impacts_opt5incl.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i mu_impacts_opt5incl.json -o fixed_mu_impacts_opt5incl


    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_freefrate_mu.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_freefrate_mu.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_freefrate_mu.root -m 125 -o mu_impacts_opt5.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i mu_impacts_opt5.json -o mu_impacts_opt5

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_freefrate_el.root -m 125 --robustFit 1 --doInitialFit
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_freefrate_el.root -m 125 --robustFit 1 --doFits --parallel 5
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/optimized5_freefrate_el.root -m 125 -o el_impacts_opt5.json
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i el_impacts_opt5.json -o el_impacts_opt5


    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/oldrun2_freefrate_el.root -m 125 --robustFit 1 --doInitialFit --freezeParameters "FSR"
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/oldrun2_freefrate_el.root -m 125 --robustFit 1 --doFits --parallel 5 --freezeParameters "FSR"
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/oldrun2_freefrate_el.root -m 125 -o el_impacts_old2.json --freezeParameters "FSR"
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i el_impacts_old2.json -o el_impacts_old2

    ttxsec/fit-stuff/oldrun2_freefrate_el_fine.root
    ttxsec/fit-stuff/oldrun2_freefrate_el.root

    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/oldrun2_freefrate_el_fine.root -m 125 --robustFit 1 --doInitialFit  --freezeParameters "TauES"
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/oldrun2_freefrate_el_fine.root -m 125 --robustFit 1 --doFits --parallel 5  --freezeParameters "TauES"
    ../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts -d ttxsec/fit-stuff/oldrun2_freefrate_el_fine.root -m 125 -o el_impacts_old2_fine.json  --freezeParameters "TauES"
    ../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i el_impacts_old2_fine.json -o el_impacts_old2_fine








## table of impacts












