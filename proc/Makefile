
sub_proc:
	ls /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/proc/queue_dir/${nt}/${proc}/${node}/com
	ssh fermi0${node} 'screen -d -m tcsh -c "source /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/proc/queue_dir/${nt}/${proc}/${node}/com "'

#make merge_proc dtags=merge-sets/jobs_hadd.dtags.large_ab nt=v25 proc=p3
merge_proc:
	for d in `cat ${dtags}`; \
	do \
	hadd merge-sets/${nt}/${proc}/$$d.root /lstore/cms/olek/outdirs/${nt}/${proc}/$$d/*.root & \
	done

merge_tt:
	for fl in `ls tt_files_*`; \
	do \
	hadd merge-sets/${nt}/${proc}/$$fl.root `sed 's,^,/lstore/cms/olek/outdirs/${nt}/${proc}/${dtag}/,' $$fl` & \
	done

#plot_mc: dtags=merge-sets/jobs_hadd.dtags.mc
plot_mc: dir=merge-sets/v25/pStage2Test2/
plot_mc:
	# qcd
	for d in MC2016_Summer16_QCD_HT-100-200 MC2016_Summer16_QCD_HT-1000-1500 MC2016_Summer16_QCD_HT-1500-2000 MC2016_Summer16_QCD_HT-200-300 MC2016_Summer16_QCD_HT-2000-Inf MC2016_Summer16_QCD_HT-300-500 MC2016_Summer16_QCD_HT-500-700 MC2016_Summer16_QCD_HT-700-1000 ; \
	do \
	python sumup_ttree_distrs.py "event_leptons[0].pt()" --ttree ttree_out --cond-com "selection_stage == 5" --histo-range 50,0,200 --output test1_lep_pt_$$d.root --histo-name ctr_old_mu_sel/qcd/NOMINAL/test1_lep_pt --save-weight ${dir}/$$d.root & \
	done
	# dy
	for d in MC2016_Summer16_DYJetsToLL_10to50_amcatnlo MC2016_Summer16_DYJetsToLL_50toInf_madgraph ; \
	do \
	python sumup_ttree_distrs.py "event_leptons[0].pt()" --ttree ttree_out --cond-com "selection_stage == 5 && gen_proc_id == 1" --histo-range 50,0,200 --output test1_lep_pt_$$d_dy_tautau.root --histo-name ctr_old_mu_sel/dy_tautau/NOMINAL/test1_lep_pt --save-weight ${dir}/$$d.root & \
	python sumup_ttree_distrs.py "event_leptons[0].pt()" --ttree ttree_out --cond-com "selection_stage == 5 && gen_proc_id == 0" --histo-range 50,0,200 --output test1_lep_pt_$$d_dy_other.root  --histo-name ctr_old_mu_sel/dy_other/NOMINAL/test1_lep_pt  --save-weight ${dir}/$$d.root & \
	done
	# single-top
	for d in MC2016_Summer16_SingleTbar_tW_5FS_powheg MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg ; \
	do \
	python sumup_ttree_distrs.py "event_leptons[0].pt()" --ttree ttree_out --cond-com "selection_stage == 5 && gen_proc_id == 10" --histo-range 50,0,200 --output test1_lep_pt_$$d_s_top_mutau.root --histo-name ctr_old_mu_sel/dy_tautau/NOMINAL/test1_lep_pt --save-weight ${dir}/$$d.root & \
	python sumup_ttree_distrs.py "event_leptons[0].pt()" --ttree ttree_out --cond-com "selection_stage == 5 && gen_proc_id == 0" --histo-range 50,0,200 --output test1_lep_pt_$$d_dy_other.root  --histo-name ctr_old_mu_sel/dy_other/NOMINAL/test1_lep_pt  --save-weight ${dir}/$$d.root & \
MC2016_Summer16_SingleT_tW_5FS_powheg





distrs_data_nominal: p=pStage2Run1
distrs_data_nominal: nt=v25
distrs_data_nominal: dtag=SingleMuon
distrs_data_nominal: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_data_nominal: cond=selection_stage == 107
distrs_data_nominal: chan=mu_sel
distrs_data_nominal: options=
distrs_data_nominal:
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" --distr Mt_lep_met_f --chan ${chan} ${bins} merge-sets/${nt}/${p}/*${dtag}*.root



distrs_mc_nominal: p=pStage2Run1
distrs_mc_nominal: nt=v25
distrs_mc_nominal: dtag=MC
distrs_mc_nominal: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_nominal: cond=selection_stage == 107
distrs_mc_nominal: chan=mu_sel
distrs_mc_nominal: options=
distrs_mc_nominal: draw="sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))"
distrs_mc_nominal: distr=Mt_lep_met_f
distrs_mc_nominal:
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight merge-sets/${nt}/${p}/*${dtag}*.root

distrs_mc_updowns: p=pStage2Run1
distrs_mc_updowns: nt=v25
distrs_mc_updowns: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_updowns: cond=selection_stage == 107
distrs_mc_updowns: chan=mu_sel
distrs_mc_updowns: options=
distrs_mc_updowns: draw="sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))"
distrs_mc_updowns: distr=Mt_lep_met_f
distrs_mc_updowns:
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys ISRUp              merge-sets/${nt}/${p}/MC*TT*isrup*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys ISRDown            merge-sets/${nt}/${p}/MC*TT*isrdown*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys FSRUp              merge-sets/${nt}/${p}/MC*TT*fsrup*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys FSRDown            merge-sets/${nt}/${p}/MC*TT*fsrdown*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys HDAMPUp            merge-sets/${nt}/${p}/MC*TT*hdampUP*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys HDAMPDown          merge-sets/${nt}/${p}/MC*TT*hdampDOWN*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys TuneCUETP8M2T4Up   merge-sets/${nt}/${p}/MC*TT*CUETP8M2T4up*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys TuneCUETP8M2T4Down merge-sets/${nt}/${p}/MC*TT*CUETP8M2T4down*.root

distrs_mc_common_weights: p=pStage2Run1
distrs_mc_common_weights: nt=v25
distrs_mc_common_weights: dtag=MC
distrs_mc_common_weights: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_common_weights: cond=selection_stage == 107
distrs_mc_common_weights: chan=mu_sel
distrs_mc_common_weights: options=
distrs_mc_common_weights:
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PUUp --distr Mt_lep_met_f    --chan ${chan}    --weight "event_weight*event_weight_PUUp/event_weight_PU"      merge-sets/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PUDown --distr Mt_lep_met_f  --chan ${chan}  --weight "event_weight*event_weight_PUUp/event_weight_PU"      merge-sets/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys bSFUp   --distr Mt_lep_met_f --chan ${chan} --weight "event_weight*event_weight_bSFUp/event_weight_bSF"    merge-sets/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys bSFDown --distr Mt_lep_met_f --chan ${chan} --weight "event_weight*event_weight_bSFDown/event_weight_bSF"  merge-sets/${nt}/${p}/*${dtag}*.root


distrs_mc_lep_weights: p=pStage2Run1
distrs_mc_lep_weights: nt=v25
distrs_mc_lep_weights: dtag=MC
distrs_mc_lep_weights: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_lep_weights: cond=selection_stage == 107
distrs_mc_lep_weights: chan=mu_sel
distrs_mc_lep_weights: options=
distrs_mc_lep_weights: sysname=LEP
distrs_mc_lep_weights:
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys ${sysname}Up --distr Mt_lep_met_f   --chan ${chan} --weight "event_weight*event_weight_LEPUp"         merge-sets/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys ${sysname}Down --distr Mt_lep_met_f --chan ${chan} --weight "event_weight*event_weight_LEPDown"     merge-sets/${nt}/${p}/*${dtag}*.root


distrs_mc_objects: p=pStage2Run1
distrs_mc_objects: nt=v25
distrs_mc_objects: dtag=MC
distrs_mc_objects: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_objects: cond=selection_stage == 107
distrs_mc_objects: chan=mu_sel
distrs_mc_objects: options=
distrs_mc_objects:
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_JERDown.Px()*event_met_JERDown.Px() + event_met_JERDown.Py()*event_met_JERDown.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_JERDown.Px()*event_leptons[0].Px() + event_met_JERDown.Py()*event_leptons[0].Py())))"    --cond-com "${cond}" ${bins} --sys JERDown   --weight event_weight  merge-sets/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_JERUp.Px()*event_met_JERUp.Px() + event_met_JERUp.Py()*event_met_JERUp.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_JERUp.Px()*event_leptons[0].Px() + event_met_JERUp.Py()*event_leptons[0].Py())))"                --cond-com "${cond}" ${bins}             --sys JERUp     --weight event_weight  merge-sets/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_JESDown.Px()*event_met_JESDown.Px() + event_met_JESDown.Py()*event_met_JESDown.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_JESDown.Px()*event_leptons[0].Px() + event_met_JESDown.Py()*event_leptons[0].Py())))"    --cond-com "${cond}" ${bins} --sys JESDown   --weight event_weight  merge-sets/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_JESUp.Px()*event_met_JESUp.Px() + event_met_JESUp.Py()*event_met_JESUp.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_JESUp.Px()*event_leptons[0].Px() + event_met_JESUp.Py()*event_leptons[0].Py())))"                --cond-com "${cond}" ${bins}             --sys JESUp     --weight event_weight  merge-sets/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_TESDown.Px()*event_met_TESDown.Px() + event_met_TESDown.Py()*event_met_TESDown.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_TESDown.Px()*event_leptons[0].Px() + event_met_TESDown.Py()*event_leptons[0].Py())))"    --cond-com "${cond}" ${bins} --sys TESDown   --weight event_weight  merge-sets/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_TESUp.Px()*event_met_TESUp.Px() + event_met_TESUp.Py()*event_met_TESUp.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_TESUp.Px()*event_leptons[0].Px() + event_met_TESUp.Py()*event_leptons[0].Py())))"                --cond-com "${cond}" ${bins}             --sys TESUp     --weight event_weight  merge-sets/${nt}/${p}/*${dtag}*.root



distrs_mc_tt_weights: p=pStage2Run1
distrs_mc_tt_weights: nt=v25
distrs_mc_tt_weights: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_tt_weights: cond=selection_stage == 107
distrs_mc_tt_weights: chan=mu_sel
distrs_mc_tt_weights: options=
distrs_mc_tt_weights:
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys TOPPTDown   --weight event_weight                               merge-sets/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys TOPPTUp     --weight "event_weight*event_weight_toppt"          merge-sets/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys FragUp      --weight "event_weight*event_weight_FragUp"         merge-sets/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys FragDown    --weight "event_weight*event_weight_FragDown"       merge-sets/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys SemilepBRUp --weight "event_weight*event_weight_SemilepBRUp"    merge-sets/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys SemilepBRDown --weight "event_weight*event_weight_SemilepBRDown" merge-sets/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PetersonUp    --weight "event_weight*event_weight_PetersonUp"   merge-sets/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PetersonDown  --weight "event_weight"                           merge-sets/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root



