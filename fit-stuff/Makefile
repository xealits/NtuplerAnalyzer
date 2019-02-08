# a bunch of fitting scripts collected in makefile
# instead of keeping them in separate files


# assume .txt ending
# needed for easy matching with workspace name
both: mucard=ttxsec/latest_datacard_mu
both: elcard=ttxsec/latest_datacard_el
both: bothcard=ttxsec/latest_datacard_both
both:
	combineCards.py mu=${mucard}.txt el=${elcasrd}.txt > ${bothcard}.txt

workspace: datacard=ttxsec/latest_datacard_mu
workspace:
	../scripts/text2workspace.py ${datacard}.txt

nll: workspace=ttxsec/latest_datacard_el.root
nll: nn=el
nll: freeze=--freezeParameters lumi_13TeV
nll: points=100
nll: min=0.5
nll: max=1.5
nll:
	combine -M MultiDimFit ${workspace} --algo grid --points ${points} --rMin ${min} --rMax ${max} --name ${nn} ${freeze}

# for toys pass
#  -t -1 --expectSignal 1
impacts: workspace=ttxsec/fit-stuff/latest_datacard_el.root
impacts: nn=el
impacts: freeze=--freezeParameters lumi_13TeV
impacts: toys=
impacts:
	../../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts ${toys} --name ${nn} ${freeze} -d ${workspace} -m 125 --robustFit 1 --doInitialFit
	../../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts ${toys} --name ${nn} ${freeze} -d ${workspace} -m 125 --robustFit 1 --doFits --parallel 5
	../../../CombineHarvester/CombineTools/scripts/combineTool.py -M Impacts ${toys} --name ${nn} ${freeze} -d ${workspace} -m 125 -o latest_impacts_${nn}.json
	../../../CombineHarvester/CombineTools/scripts/plotImpacts_my.py -i latest_impacts_${nn}.json -o latest_impacts_${nn}

impacts_toys: toys=-t -1 --expectSignal 1
impacts_toys: impacts

# for toys pass
#  -t -1 --expectSignal 1
fit_unc: workspace=ttxsec/latest_datacard_mu.root
fit_unc: nn=el
fit_unc: points=100
fit_unc: min=0.5
fit_unc: max=1.5
fit_unc: toys=
fit_unc:
	combine -M MultiDimFit --algo none              ${toys} --rMin ${min} --rMax ${max} ${workspace} --saveWorkspace -n ${nn}Bestfit
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan                   -n ${nn}NoSysWLumi
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff,tau_fakes,lumi_13TeV" -n ${nn}NoTau
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit                                                     -n ${nn}FullUncertaintyWLumi
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV"                     -n ${nn}FullUncertaintyNoLumi
	plot1DScan.py higgsCombine${nn}FullUncertaintyWLumi.MultiDimFit.mH120.root --others 'higgsCombine${nn}FullUncertaintyNoLumi.MultiDimFit.mH120.root:lumi:4' 'higgsCombine${nn}NoSys.MultiDimFit.mH120.root:stat:2' --breakdown lumi,syst,stat
	mv scan.pdf scan_${nn}.pdf
	mv scan.png scan_${nn}.png

fit_unc_toys: fit_unc=-t -1 --expectSignal 1
fit_unc_toys: fit_unc

# pulls are part-manual
pre_pulls: workspace=ttxsec/latest_datacard_mu.root
pre_pulls: nn=mu
pre_pulls: freeze=--freezeParameters lumi_13TeV
pre_pulls:
	combine -M FitDiagnostics ${workspace}   --name ${nn}Shape    --noMCbonly 0 ${freeze}
	python ttxsec/fit-stuff/diffNuisances.py fitDiagnostics${nn}Shape.root -g pulls_of_nuisances_${nn}_shape.root   -o

full_card_fit: workspace=ttxsec/latest_datacard_el.root
full_card_fit: nn=el
full_card_fit: freeze=--freezeParameters lumi_13TeV
full_card_fit: points=100
full_card_fit: min=0.5
full_card_fit: max=1.5
full_card_fit: nll impacts impacts_toys fit_unc fit_unc_toys pre_pulls

full_three_card_fit: mucard=ttxsec/latest_datacard_mu
full_three_card_fit: elcard=ttxsec/latest_datacard_el
full_three_card_fit: bothcard=ttxsec/latest_datacard_both
full_three_card_fit: freeze=--freezeParameters lumi_13TeV
full_three_card_fit: points=100
full_three_card_fit: min=0.5
full_three_card_fit: max=1.5
full_three_card_fit: both
full_three_card_fit: muname=mu
full_three_card_fit: elname=el
full_three_card_fit: bothname=both
full_three_card_fit:
	make workspace datacard=${elcard}
	make workspace datacard=${mucard}
	make workspace datacard=${bothcard}
	make full_card_fit workspace=${elcard}.root   nn=${elname}
	make full_card_fit workspace=${mucard}.root   nn=${muname}
	make full_card_fit workspace=${bothcard}.root nn=${bothname}
	python new_impact_table.py --merge-pdf latest_impacts_${elname}.json latest_impacts_${muname}.json latest_impacts_${bothname}.json > impact_table_${elname}_${muname}_${bothname}

