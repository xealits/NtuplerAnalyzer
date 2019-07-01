# a bunch of fitting scripts collected in makefile
# instead of keeping them in separate files


# assume .txt ending
# needed for easy matching with workspace name
both: mucard=ttxsec/latest_datacard_mu
both: elcard=ttxsec/latest_datacard_el
both: bothcard=ttxsec/latest_datacard_both
both:
	combineCards.py mu=${mucard}.txt el=${elcard}.txt > ${bothcard}.txt

workspace: datacard=ttxsec/latest_datacard_mu
workspace:
	../scripts/text2workspace.py ${datacard}.txt

fit: workspace=ttxsec/datacard_v37_test172_bunchFULLFIT1_mu.root
fit: nn=el
fit: freeze=--freezeParameters lumi_13TeV
fit: options=
fit:
	combine -M FitDiagnostics   ${workspace} --name ${nn} ${freeze} ${options}

nll: workspace=ttxsec/latest_datacard_el.root
nll: nn=el
nll: freeze=--freezeParameters lumi_13TeV
nll: options=
nll: points=100
nll: min=0.5
nll: max=1.5
nll:
	combine -M MultiDimFit ${workspace} --algo grid --points ${points} --rMin ${min} --rMax ${max} --name ${nn} ${freeze} ${options}

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
impacts_toys: freeze=--freezeParameters lumi_13TeV
impacts_toys:
	make impacts nn=${nn}_toys toys="${toys}" workspace=${workspace} freeze="${freeze}"


# breakdown from Yuta's question
#fit_unc_breakdown_yutas_q: workspace_nomc=ttxsec/latest_datacard_mu.root
fit_unc_breakdown_yutas_q: workspace=ttxsec/latest_datacard_mu.root
fit_unc_breakdown_yutas_q: chan=el
fit_unc_breakdown_yutas_q: nn=el
fit_unc_breakdown_yutas_q: points=100
fit_unc_breakdown_yutas_q: min=0.5
fit_unc_breakdown_yutas_q: max=1.5
fit_unc_breakdown_yutas_q: toys=
fit_unc_breakdown_yutas_q:
	combine -M MultiDimFit --algo none              ${toys} --rMin ${min} --rMax ${max} ${workspace} --saveWorkspace -n ${nn}Bestfit
	#combine -M MultiDimFit --algo none              ${toys} --rMin ${min} --rMax ${max} ${workspace_nomc} --saveWorkspace -n ${nn}BestfitNoMC
	# trying to merge MC and no MC
	# stat unc
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan                   -n ${nn}NoSys
	# no tau ID and MC stat
	#combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}BestfitNoMC.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff,lumi_13TeV" -n ${nn}NoTauNoMC
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff,lumi_13TeV,prop_bin${chan}_selTight_ljout_bin8,prop_bin${chan}_selTight_ljout_bin7,prop_bin${chan}_selTight_ljout_bin6,prop_bin${chan}_selTight_ljout_bin5,prop_bin${chan}_selTight_ljout_bin4,prop_bin${chan}_selTight_ljout_bin3,prop_bin${chan}_selTight_ljout_bin2,prop_bin${chan}_selTight_ljout_bin1,prop_bin${chan}_selTight_ljout_bin0,prop_bin${chan}_selTight_lj_bin8,prop_bin${chan}_selTight_lj_bin7,prop_bin${chan}_selTight_lj_bin6,prop_bin${chan}_selTight_lj_bin5,prop_bin${chan}_selTight_lj_bin4,prop_bin${chan}_selTight_lj_bin3,prop_bin${chan}_selTight_lj_bin2,prop_bin${chan}_selTight_lj_bin1,prop_bin${chan}_selTight_lj_bin0" -n ${nn}NoTauNoMC
	# no MC stat
	#combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}BestfitNoMC.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV" -n ${nn}NoMC
	# correct no MC
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,prop_bin${chan}_selTight_ljout_bin8,prop_bin${chan}_selTight_ljout_bin7,prop_bin${chan}_selTight_ljout_bin6,prop_bin${chan}_selTight_ljout_bin5,prop_bin${chan}_selTight_ljout_bin4,prop_bin${chan}_selTight_ljout_bin3,prop_bin${chan}_selTight_ljout_bin2,prop_bin${chan}_selTight_ljout_bin1,prop_bin${chan}_selTight_ljout_bin0,prop_bin${chan}_selTight_lj_bin8,prop_bin${chan}_selTight_lj_bin7,prop_bin${chan}_selTight_lj_bin6,prop_bin${chan}_selTight_lj_bin5,prop_bin${chan}_selTight_lj_bin4,prop_bin${chan}_selTight_lj_bin3,prop_bin${chan}_selTight_lj_bin2,prop_bin${chan}_selTight_lj_bin1,prop_bin${chan}_selTight_lj_bin0" -n ${nn}NoMC
	#combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV,prop_bin_el_selTight_ljout_bin8" -n ${nn}NoMC
	# all, with MC
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV"                     -n ${nn}FullUncertaintyNoLumi
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit                                                     -n ${nn}FullUncertaintyWLumi
	# plot
	plot1DScan.py --output scan_${nn} higgsCombine${nn}FullUncertaintyNoLumi.MultiDimFit.mH120.root --others 'higgsCombine${nn}NoMC.MultiDimFit.mH120.root:No MC:4' 'higgsCombine${nn}NoTauNoMC.MultiDimFit.mH120.root:No MC No Tau:3' 'higgsCombine${nn}NoSys.MultiDimFit.mH120.root:stat:2' --breakdown "MC,Tau,other syst,stat"


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
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters all --fastScan                   -n ${nn}NoSys
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "tauID_eff,tau_fakes,lumi_13TeV" -n ${nn}NoTau
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit                                                     -n ${nn}FullUncertaintyWLumi
	combine -M MultiDimFit --algo grid --points 100 ${toys} --rMin ${min} --rMax ${max} higgsCombine${nn}Bestfit.MultiDimFit.mH120.root --snapshotName MultiDimFit --freezeParameters "lumi_13TeV"                     -n ${nn}FullUncertaintyNoLumi
	plot1DScan.py --output scan_${nn} higgsCombine${nn}FullUncertaintyWLumi.MultiDimFit.mH120.root --others 'higgsCombine${nn}FullUncertaintyNoLumi.MultiDimFit.mH120.root:lumi:4' 'higgsCombine${nn}NoSys.MultiDimFit.mH120.root:stat:2' --breakdown lumi,syst,stat

fit_unc_toys: workspace=ttxsec/latest_datacard_mu.root
fit_unc_toys: nn=el
fit_unc_toys: points=100
fit_unc_toys: min=0.5
fit_unc_toys: max=1.5
fit_unc_toys: toys=-t -1 --expectSignal 1
fit_unc_toys:
	make fit_unc toys="${toys}" min=${min} max=${max} points=${points} nn=${nn}_Expected workspace=${workspace}

# pulls are part-manual
pre_pulls: workspace=ttxsec/latest_datacard_mu.root
pre_pulls: nn=mu
pre_pulls: freeze=--freezeParameters lumi_13TeV
pre_pulls:
	combine -M FitDiagnostics ${workspace}   --name ${nn}Shape    --noMCbonly 0 ${freeze}
	python ttxsec/fit-stuff/diffNuisances.py fitDiagnostics${nn}Shape.root -g pulls_of_nuisances_${nn}_shape.root   -o

# same pulls with less systematics for readability
pre_pulls_short: workspace=ttxsec/latest_datacard_mu.root
pre_pulls_short: nn=mu
pre_pulls_short: freeze=--freezeParameters lumi_13TeV --freezeNuisanceGroup tt_th_pdf_only
pre_pulls_short:
	combine -M FitDiagnostics ${workspace}   --name ${nn}ShapeShort    --noMCbonly 0 ${freeze}
	python ttxsec/fit-stuff/diffNuisances.py fitDiagnostics${nn}ShapeShort.root -g pulls_of_nuisances_${nn}_shapeshort.root   -o

full_card_fit: workspace=ttxsec/latest_datacard_el.root
full_card_fit: nn=el
full_card_fit: freeze=--freezeParameters lumi_13TeV
full_card_fit: points=100
full_card_fit: min=0.5
full_card_fit: max=1.5
full_card_fit: nll impacts impacts_toys fit_unc fit_unc_toys pre_pulls pre_pulls_short

postfit_distrs: fit_tag=v37_test13_bunchFULLFIT2
postfit_distrs: bins=0,20,40,60,80,100,130,160,200,250
postfit_distrs:
	combine -M FitDiagnostics ttxsec/datacard_${fit_tag}_mu.root  --saveShapes --saveWithUncertainties --name ${fit_tag}_MuShapes --freezeParameters lumi_13TeV
	combine -M FitDiagnostics ttxsec/datacard_${fit_tag}_el.root  --saveShapes --saveWithUncertainties --name ${fit_tag}_ElShapes --freezeParameters lumi_13TeV
	# plot results
	python ttxsec/fit-stuff/stack_postfit.py --custom-bins "${bins}" fitDiagnostics${fit_tag}_MuShapes.root mu_sel_lj    --mu --ratio --ratio-range 0.5 --y-max 150
	python ttxsec/fit-stuff/stack_postfit.py --custom-bins "${bins}" fitDiagnostics${fit_tag}_MuShapes.root mu_sel_ljout --mu --ratio --ratio-range 0.5 --y-max 150
	python ttxsec/fit-stuff/stack_postfit.py --custom-bins "${bins}" fitDiagnostics${fit_tag}_ElShapes.root el_sel_lj    --lumi 35.8 --ratio --ratio-range 0.5 --y-max 100
	python ttxsec/fit-stuff/stack_postfit.py --custom-bins "${bins}" fitDiagnostics${fit_tag}_ElShapes.root el_sel_ljout --lumi 35.8 --ratio --ratio-range 0.5 --y-max 100


full_three_card_fit: freeze=--freezeParameters lumi_13TeV
full_three_card_fit: points=100
full_three_card_fit: min=0.5
full_three_card_fit: max=1.5
full_three_card_fit: fit_tag=v37_test13_bunchFULLFIT2
full_three_card_fit:
	# workspaces
	make both mucard=ttxsec/datacard_${fit_tag}_mu elcard=ttxsec/datacard_${fit_tag}_el bothcard=ttxsec/datacard_${fit_tag}_both
	make workspace datacard=ttxsec/datacard_${fit_tag}_el
	make workspace datacard=ttxsec/datacard_${fit_tag}_mu
	make workspace datacard=ttxsec/datacard_${fit_tag}_both
	# also the cards without MC stat for shorter, readable pulls
	head -n-2 ttxsec/datacard_${fit_tag}_el.txt >  ttxsec/datacard_${fit_tag}_el_noMCstat.txt
	head -n-2 ttxsec/datacard_${fit_tag}_mu.txt >  ttxsec/datacard_${fit_tag}_mu_noMCstat.txt
	make both mucard=ttxsec/datacard_${fit_tag}_mu_noMCstat elcard=ttxsec/datacard_${fit_tag}_el_noMCstat bothcard=ttxsec/datacard_${fit_tag}_both_noMCstat
	make workspace datacard=ttxsec/datacard_${fit_tag}_el_noMCstat
	make workspace datacard=ttxsec/datacard_${fit_tag}_mu_noMCstat
	make workspace datacard=ttxsec/datacard_${fit_tag}_both_noMCstat
	# fits
	make full_card_fit workspace=ttxsec/datacard_${fit_tag}_el.root   nn=el_${fit_tag}
	make full_card_fit workspace=ttxsec/datacard_${fit_tag}_mu.root   nn=mu_${fit_tag}
	make full_card_fit workspace=ttxsec/datacard_${fit_tag}_both.root nn=both_${fit_tag}
	python new_impact_table.py --merge-pdf latest_impacts_el_${fit_tag}.json latest_impacts_mu_${fit_tag}.json latest_impacts_both_${fit_tag}.json > impact_table_${fit_tag}
	# also the short pulls
	make pre_pulls_short workspace=ttxsec/datacard_${fit_tag}_el_noMCstat   nn=el_${fit_tag}_noMCstat
	make pre_pulls_short workspace=ttxsec/datacard_${fit_tag}_mu_noMCstat   nn=mu_${fit_tag}_noMCstat
	make pre_pulls_short workspace=ttxsec/datacard_${fit_tag}_both_noMCstat nn=both_${fit_tag}_noMCstat
	# plot NLL
	python ttxsec/fit-stuff/plot_uncertainty_scans.py ${fit_tag}
	# and pulls
	python ttxsec/fit-stuff/pulls_plot.py pulls_of_nuisances_both_${fit_tag}_noMCstat_shapeshort.root
	python ttxsec/fit-stuff/pulls_plot.py pulls_of_nuisances_mu_${fit_tag}_noMCstat_shapeshort.root
	python ttxsec/fit-stuff/pulls_plot.py pulls_of_nuisances_el_${fit_tag}_noMCstat_shapeshort.root
	# all pulls
	python ttxsec/fit-stuff/pulls_plot.py pulls_of_nuisances_both_${fit_tag}_shape.root
	python ttxsec/fit-stuff/pulls_plot.py pulls_of_nuisances_mu_${fit_tag}_shape.root
	python ttxsec/fit-stuff/pulls_plot.py pulls_of_nuisances_el_${fit_tag}_shape.root
	# post-fit distrs
	make postfit_distrs fit_tag=${fit_tag}



#plots_for_export: nt=v37
#plots_for_export: proc=test13
#plots_for_export: distrs=bunchFULLFIT2
#plots_for_export: fit=2
.ONESHELL:
plots_for_export: fit_tag=v37_test13_bunchFULLFIT2
plots_for_export:
	mkdir -p export/${fit_tag}/
	cp *${fit_tag}*.png *${fit_tag}*.pdf *${fit_tag}*.json export/${fit_tag}/
	cp impact_table_*${fit_tag}*                           export/${fit_tag}/
	mkdir -p export/${fit_tag}/renames/
	#
	cd export/${fit_tag}/renames/
	# impacts
	ln -s ../latest_impacts_mu_*_toys.pdf   prefit_mu_impacts_toys.pdf
	ln -s ../latest_impacts_el_*_toys.pdf   prefit_el_impacts_toys.pdf
	ln -s ../latest_impacts_both_*_toys.pdf prefit_both_impacts_toys.pdf
	ln -s ../latest_impacts_mu_*${fit_tag}.pdf    postfit_mu_impacts.pdf
	ln -s ../latest_impacts_el_*${fit_tag}.pdf    postfit_el_impacts.pdf
	ln -s ../latest_impacts_both_*${fit_tag}.pdf  postfit_both_impacts.pdf
	# sys breakdown (the plots are not used in AN, only the numbers)
	ln -s ../scan_*mu*${fit_tag}.png    mu_uncert_break_scan.png
	ln -s ../scan_*el*${fit_tag}.png    el_uncert_break_scan.png
	ln -s ../scan_*both*${fit_tag}.png  both_uncert_break_scan.png
	ln -s ../scan_*mu*Expected.png      mu_uncert_break_scan_exp.png
	ln -s ../scan_*el*Expected.png      el_uncert_break_scan_exp.png
	ln -s ../scan_*both*Expected.png    both_uncert_break_scan_exp.png
	# pdfs
	ln -s ../scan_*mu*${fit_tag}.pdf    mu_uncert_break_scan.pdf
	ln -s ../scan_*el*${fit_tag}.pdf    el_uncert_break_scan.pdf
	ln -s ../scan_*both*${fit_tag}.pdf  both_uncert_break_scan.pdf
	ln -s ../scan_*mu*Expected.pdf      mu_uncert_break_scan_exp.pdf
	ln -s ../scan_*el*Expected.pdf      el_uncert_break_scan_exp.pdf
	ln -s ../scan_*both*Expected.pdf    both_uncert_break_scan_exp.pdf
	# TODO: NLL
	ln -s ../uncertainty_scans_${fit_tag}_el_exp.png    uncertainty_scans_el_exp.png
	ln -s ../uncertainty_scans_${fit_tag}_mu_exp.png    uncertainty_scans_mu_exp.png
	ln -s ../uncertainty_scans_${fit_tag}_both_exp.png  uncertainty_scans_both_exp.png
	ln -s ../uncertainty_scans_${fit_tag}_el_exp_obs.png    uncertainty_scans_el_exp_obs.png
	ln -s ../uncertainty_scans_${fit_tag}_mu_exp_obs.png    uncertainty_scans_mu_exp_obs.png
	ln -s ../uncertainty_scans_${fit_tag}_both_exp_obs.png  uncertainty_scans_both_exp_obs.png
	# pulls
	# 2 sets: all and readable
	ln -s ../pulls_of_nuisances_both_${fit_tag}_noMCstat_shapeshort.png  postfit_both_nuisances_lessparams.png
	ln -s ../pulls_of_nuisances_mu_${fit_tag}_noMCstat_shapeshort.png    postfit_mu_nuisances_lessparams.png
	ln -s ../pulls_of_nuisances_el_${fit_tag}_noMCstat_shapeshort.png    postfit_el_nuisances_lessparams.png
	ln -s ../pulls_of_nuisances_both_${fit_tag}_shape.png  postfit_both_nuisances.png
	ln -s ../pulls_of_nuisances_mu_${fit_tag}_shape.png    postfit_mu_nuisances.png
	ln -s ../pulls_of_nuisances_el_${fit_tag}_shape.png    postfit_el_nuisances.png
	# post-fit distrs
	ln -s ../postfit-distr_fitDiagnostics${fit_tag}_MuShapes_mu_sel_lj.png     postfit_mu_Mt_lj.png
	ln -s ../postfit-distr_fitDiagnostics${fit_tag}_MuShapes_mu_sel_ljout.png  postfit_mu_Mt_ljout.png
	ln -s ../postfit-distr_fitDiagnostics${fit_tag}_ElShapes_el_sel_lj.png     postfit_el_Mt_lj.png
	ln -s ../postfit-distr_fitDiagnostics${fit_tag}_ElShapes_el_sel_ljout.png  postfit_el_Mt_ljout.png



