Installation
============

The module works with CMSSW_8_0_26+ (CMSSW_8_0_29 recommended)
and requires [TopQuarkAnalysis/BFragmentation](https://gitlab.cern.ch/CMS-TOPPAG/BFragmentationAnalyzer) wich produces some values for [certain systematics](https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Fragmentation).

Other modules are standard and it should not need anything additional.

An example installation script:

    cmsrel CMSSW_8_0_30
    cd CMSSW_8_0_30/src 
    cmsenv

    # copy of the BFragmentation installation
    mkdir TopQuarkAnalysis
    cd TopQuarkAnalysis
    git clone ssh://git@gitlab.cern.ch:7999/CMS-TOPPAG/BFragmentationAnalyzer.git
    cd -
    scram b

    # installing the NtuplerAnalyzer
    mkdir UserCode
    cd UserCode
    git clone git@github.com:xealits/NtuplerAnalyzer.git
    scram b -j 10

If something crashes during the compilation you can cd into separate directories,
compile them one by one and then compile-link everything together.
Like:

    cd CMSSW_8_0_30/src/UserCode/NtuplerAnalyzer/src/
    scram b -j 5
    cd ../plugins/
    scram b -j 5


Then a test run should work with:

    cd CMSSW_8_0_30/src/UserCode/NtuplerAnalyzer/
    cmsRun python/ConfFile_cfg.py

--- it accesses test files in my CERNBox (EOS) area, for example 1 TT file:

    $ ls -l /eos/user/o/otoldaie/TT_165F54A0-A3BE-E611-B3F7-0025905A606A.root
    -rw-r--r--. 1 otoldaie zh 3816687534 Mar 15  2018 /eos/user/o/otoldaie/TT_165F54A0-A3BE-E611-B3F7-0025905A606A.root

--- so, now it should be readable for everybody.
However the directory is not readable. For some reason I cannot chmod a directory on EOS.




Old notes (keeping just in case)
--------------------------------

In the past the RecoilCorrections for DY were calculated online and the ntupler required [HTT-utilities/RecoilCorrections](https://github.com/CMS-HTT/RecoilCorrections/blob/master/instructions.txt).
Now it is done "offline" in processing the output ttrees from the ntupler.
The ntupler saves the parameters needed for the corrections (pT of visible particles and neutrinos in DY etc).




Structure
=========

The repository contains 2 projects:

* ntupler, the thing that runs over crab on grid and reduces the MINIAODs to TTree with only important parameters
* proc, the scripts for processing resulting TTree

The proc scripts are in `NtuplerAnalyzer/proc/` directory.

The ntupler is a standard EDM module and uses standard directories:

* `plugins/` with couple EDM plugins, but the main one is `plugins/NtuplerAnalyzer.cc`
* `src/` with functions for muon/electron/jet IDs, a bunch of old helper functions in ParUtils and MacroUtils etc
* `python/` with the default config for the ntupler in `python/CfiFile_cfi.py` and a test config in `python/ConfFile_cfg.py`
* also a script for generating crab jobs in `NtuplerAnalyzer/make_crab_cfgs.py`, which uses template files for jobs

   + `TEMPLATE_crabConf.py` the config for the crab, according to which the crab will generate a bunch of jobs and run them on different computers in grid
   + `TEMPLATE_cfg_for_crab.py` the cmsRun cfg.py file which will be run by crab on some computer in the grid

* also the crab jobs submition uses files with dataset names like `dsets_all_usual_mc` and the file with info on datasets `dsets_info.yaml`, which contains "dtag" of the dataset -- it is the dataset name used in the proc scripts




Example of crab jobs submition
==============================

It can be found in `NtuplerAnalyzer/last_jobs`.
For example:

    $ head last_jobs
    # v35 safe run for antiiso mc qcd
    for c in `cat dsets_all_usual_qcd_ext`
    do
    python make_crab_cfgs.py $c v35 tauIDantiIso
    done
    
    for c in `ls python/crab_cfgs/v35/*crab.py`
    do
    crab submit -c $c &
    done

--- the `make_crab_cfgs.py` script makes the jobs in `python/crab_cfgs/v35/` and then crab submits them. Here the for loop launches all submitions in parallel. It is not very safe: it might not pick up your grid credentials or may create too many processes etc.

Don't forget to have your grid credentials before launching crab.
Otherwise it will ask you to initialize.
Example on my lxplus account:

    source /cvmfs/cms.cern.ch/crab3/crab.sh
    voms-proxy-init --voms cms --out ~/x509_proxy && export X509_USER_PROXY=`readlink -e ~/x509_proxy `
    export X509_USER_PROXY=/afs/cern.ch/user/o/otoldaie/x509_proxy





ntupler operation
=================

Ntupler reduces the MINIAOD events in both ways:
1) it saves less events in the TTree according to some conditions like "is there only 1 good lepton in the event?" etc,
2) it reduces the size of each stored event by saving only necessary parameters, these are simple types like int/double/LorentzVector or vectors of those (like a vector of LorentzVectors for momenta of all jets in the event etc), when in MINIAOD you have much richer PAT classes etc

In principle, one can save all events with ntupler. Of course, it will require more space for output.
But I saw that in data (SingleElectron and SingleMuon datasets) most of events pass the requirements anyway.
It should be otherwise in TT MC dataset. But I did not check how large will be the increase in the output.




Typical sizes and parameters of crab jobs
-----------------------------------------

TOMEASURE


