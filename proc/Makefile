qwatch: n=10
qwatch:
	watch -n ${n} 'qstat | grep " r     " | wc -l ; qstat | wc -l'

clear_gstore_failed: nt=v37
clear_gstore_failed:
	for f in `find gstore_outdirs/${nt}/ -name failed -type d | sed 's,gstore_outdirs/,,'`; do \
	    gfal-rm -r "srm://srm01.ncg.ingrid.pt:8444/srm/managerv2?SFN=/cmst3/store/user/otoldaie/$$f" ; \
	done

check_job_error_not_found:
	cat `wc -l job_*.e*[0-9] | grep '^   12 ' | sed 's/^   12 //'` | grep ' does not exist'

batch_plots: batch_jobs=batch_jobs
batch_plots:
	ls ${batch_jobs}/ | sed -e 's,^\(.*\)$$,qsub -l h_vmem=512M "/lstore/cms/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/proc/'${batch_jobs}'/\1",' | grep -v submit > ${batch_jobs}/submit

find_missing_jobs: nt=v37
find_missing_jobs: proc=test8
find_missing_jobs:
	cat queue_dir/${nt}/${proc}/jobs_dir/job_* | sed -n '/[^/]*\.root/ s,^.*/\([^/][^/]*\.root\).*,\1, p' | sort > job_files
	find lstore_outdirs/v37/test8/ -name "*root" -exec basename '{}' \; 2> /dev/null | sort > job_files_produced
	diff job_files job_files_produced | tee job_files_missing

#sed -e 's,^,qsub -l h_vmem=${job_memory} "/lstore/cms/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/proc/./queue_dir/${nt}/${proc}/jobs_dir/,' -e 's,\.o[0-9]*$$,",' jobs_failed > jobs_failed_submit
find_failed_jobs: job_dir=
find_failed_jobs: job_memory=1G
find_failed_jobs: job_name_tag=job
find_failed_jobs:
	# not always correct:
	#wc -l job_*.e* | grep " [0-9][0-9] " > jobs_failed
	grep --files-with-matches "output written,"    ${job_name_tag}_*.o* | sort > jobs_produced
	grep --files-with-matches "output file exists" ${job_name_tag}_*.o* | sort > jobs_output_exists
	cat jobs_produced jobs_output_exists | sort > jobs_done
	ls ${job_name_tag}_*.o* | sort > jobs_all
	diff jobs_all jobs_done | sed -n '/${job_name_tag}_/ s,^< ,, p' > jobs_failed
	# add the submition line for each job
	sed -e 's,^,qsub -l h_vmem=${job_memory} "${job_dir},' -e 's,\.o[0-9]*$$,",' jobs_failed > jobs_failed_submit
	wc -l jobs_failed_submit

find_failed_jobs_stage2: nt=v37
find_failed_jobs_stage2: proc=test14
find_failed_jobs_stage2: job_name_tag=job
find_failed_jobs_stage2:
	make find_failed_jobs job_dir=/lstore/cms/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/proc/./queue_dir/${nt}/${proc}/jobs_dir/ job_memory=1G job_name_tag=${job_name_tag}

find_failed_jobs_plots: job_dir=batch_jobs/
find_failed_jobs_plots: job_name_tag=j
find_failed_jobs_plots:
	make find_failed_jobs job_dir=${job_dir} job_memory=512M job_name_tag=${job_name_tag}

find_failed_plots: correct_num=8
find_failed_plots:
	wc -l j_*.e* | sort | grep -v "^[ ][ ]*${correct_num} \|^[ ][ ]*[0-9][0-9]* total" > failed_plots
	sed -i -e 's,\.e[0-9][0-9]*$$,",' -e 's,^[ ][ ]*[0-9][0-9]*[ ][ ]*,qsub -l h_vmem=512M "/lstore/cms/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/proc/batch_jobs/,' failed_plots

resub_failed_jobs_stage2: nt=v37
resub_failed_jobs_stage2: proc=test14
resub_failed_jobs_stage2: job_name_tag=job
resub_failed_jobs_stage2: find_failed_jobs_stage2
	rm j_*[0-9]
	sh ./jobs_failed_submit
	make qwatch

resub_failed_jobs_plots: job_dir=batch_jobs/
resub_failed_jobs_plots: job_name_tag=j
resub_failed_jobs_plots: find_failed_jobs_plots
	rm j_*[0-9]
	sh ./jobs_failed_submit
	make qwatch

check_job_output: out=jobsums/distrs/v37_test175_bunchFULLFIT1/all/
check_job_output:
	find ${out} -name "*root" -size -10k | wc -l
	find ${out} -name "*root"            | wc -l

test_pattern: nt=v37
test_pattern: proc=test14
test_pattern: job_name_tag=job
test_pattern: pattern=${nt}_${proc}_${job_name_tag}
test_pattern:
	echo ${pattern}

sumup_weights: nt=v37
sumup_weights: proc=test5
sumup_weights: distrs=bunch1
sumup_weights:
	for dtag in `ls lstore_outdirs/${nt}/${proc}/`; \
	do \
	echo $$dtag ; \
	mkdir -p jobsums/distrs/${nt}_${proc}_${distrs}/all/$$dtag/weights/ ; \
	python sumup_weights.py --overwrite jobsums/distrs/${nt}_${proc}_${distrs}/all/$$dtag/weights/the_sums.root lstore_outdirs/${nt}/${proc}/$$dtag/$$dtag*.root ; \
	done

std_changroups: nt=v37
std_changroups: proc=test5
std_changroups: distrs=bunch1
std_changroups:
	# groups:
	# control
	# presel
	# sel
	python loop_selections.py --options " --wjets-nup-cut " --nt ${nt} --proc ${proc} --distrs ${distrs} --chan-groups ${changroup}

# groups:
# control
# presel
# sel
std_changroups_all: nt=v37
std_changroups_all: proc=test5
std_changroups_all: distrs=bunch1
std_changroups_all: changroup=wjets,wjets_ss,dy_dileptons,tt_dileptons,tt_presel_cands,tt_presel_lj_mu,tt_presel_lj_el,fit_tt_leptau,fit_tt_leptau_lj,fit_tt_leptau_ljout,fit_tt_leptau_Vloose,fit_tt_leptau_Vloose_lj,fit_tt_leptau_Vloose_ljout,tt_leptau,tt_leptau_lj,tt_leptau_ljout,tt_leptau_Vloose,tt_leptau_Vloose_lj,tt_leptau_Vloose_ljout
std_changroups_all: sumup_weights

std_changroups_control: nt=v37
std_changroups_control: proc=test5
std_changroups_control: distrs=bunch1
std_changroups_control: changroup=wjets,wjets_ss,dy_dileptons,tt_dileptons
std_changroups_control: sumup_weights

std_changroups_presel: nt=v37
std_changroups_presel: proc=test5
std_changroups_presel: distrs=bunch1
std_changroups_presel: changroup=tt_presel_cands,tt_presel_lj_mu,tt_presel_lj_el
std_changroups_presel: sumup_weights

std_changroups_sel_fit: nt=v37
std_changroups_sel_fit: proc=test5
std_changroups_sel_fit: distrs=bunch1
std_changroups_sel_fit: changroup=fit_tt_leptau,fit_tt_leptau_lj,fit_tt_leptau_ljout,fit_tt_leptau_Vloose,fit_tt_leptau_Vloose_lj,fit_tt_leptau_Vloose_ljout,tt_leptau,tt_leptau_lj,tt_leptau_ljout,tt_leptau_Vloose,tt_leptau_Vloose_lj,tt_leptau_Vloose_ljout
std_changroups_sel_fit: sumup_weights


cur_gstore_lstore_jobs: nt=v37
cur_gstore_lstore_jobs: proc=test1
cur_gstore_lstore_jobs: pattern=*MC
cur_gstore_lstore_jobs:
	find gstore_outdirs/${nt}/ -name "${pattern}*root" | sed 's,^.*/,,' | sort > cur_jobs_gstore
	find lstore_outdirs/${nt}/${proc} -name "${pattern}*root" | sed 's,^.*/,,' | sort > cur_jobs_lstore

sub_proc:
	ls /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/proc/queue_dir/${nt}/${proc}/${node}/com
	ssh fermi0${node} 'screen -d -m tcsh -c "source /exper-sw/cmst3/cmssw/users/olek/CMSSW_8_0_26_patch1/src/UserCode/NtuplerAnalyzer/proc/queue_dir/${nt}/${proc}/${node}/com "'

#make merge_proc dtags=merge-sets/jobs_hadd.dtags.large_ab nt=v25 proc=p3
merge_proc: merge_dir=merge-sets/
merge_proc:
	echo `date` merge_proc nt=${nt} proc=${proc} dtags=${dtags} >> logs
	for d in `cat ${dtags}`; \
	do \
	hadd ${merge_dir}/${nt}/${proc}/$$d.root /lstore/cms/olek/outdirs/${nt}/${proc}/$$d/*.root & \
	done

merge_proc_parallel: merge_dir=merge-sets/
merge_proc_parallel: n_proc=4
merge_proc_parallel: dtags=jobs_hadd.dtags
merge_proc_parallel:
	echo `date` merge_proc_parallel nt=${nt} proc=${proc} dtags=${dtags} >> logs
	cat ${dtags} | xargs -P ${n_proc} -I DTAG sh -c "hadd ${merge_dir}/${nt}/${proc}/DTAG.root /lstore/cms/olek/outdirs/${nt}/${proc}/DTAG/*.root || true"

merge_tt: merge_dir=merge-sets/
merge_tt:
	echo `date` merge_tt nt=${nt} proc=${proc} dtags=${dtags} >> logs
	for fl in `ls tt_files_*`; \
	do \
	hadd ${merge_dir}/${nt}/${proc}/$$fl.root `sed 's,^,/lstore/cms/olek/outdirs/${nt}/${proc}/${dtag}/,' $$fl` & \
	done

hadds_distrs:
	for dtag in `cat ${dtag_file}`; \
	do \
	hadd ${distr_dir}_hadds/$$dtag.root ${distr_dir}/$$dtag*root & \
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
distrs_data_nominal: merge_dir=merge-sets/
distrs_data_nominal:
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" --distr Mt_lep_met_f --chan ${chan} ${bins} ${merge_dir}/${nt}/${p}/*${dtag}*.root



distrs_mc_nominal: p=pStage2Run1
distrs_mc_nominal: nt=v25
distrs_mc_nominal: dtag=MC
distrs_mc_nominal: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_nominal: cond=selection_stage == 107
distrs_mc_nominal: chan=mu_sel
distrs_mc_nominal: options=
distrs_mc_nominal: draw="sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))"
distrs_mc_nominal: distr=Mt_lep_met_f
distrs_mc_nominal: merge_dir=merge-sets/
distrs_mc_nominal:
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight ${merge_dir}/${nt}/${p}/*${dtag}*.root

distrs_mc_updowns: p=pStage2Run1
distrs_mc_updowns: nt=v25
distrs_mc_updowns: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_updowns: cond=selection_stage == 107
distrs_mc_updowns: chan=mu_sel
distrs_mc_updowns: options=
distrs_mc_updowns: draw="sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))"
distrs_mc_updowns: distr=Mt_lep_met_f
distrs_mc_updowns: merge_dir=merge-sets/
distrs_mc_updowns:
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys ISRUp              ${merge_dir}/${nt}/${p}/MC*TT*isrup*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys ISRDown            ${merge_dir}/${nt}/${p}/MC*TT*isrdown*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys FSRUp              ${merge_dir}/${nt}/${p}/MC*TT*fsrup*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys FSRDown            ${merge_dir}/${nt}/${p}/MC*TT*fsrdown*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys HDAMPUp            ${merge_dir}/${nt}/${p}/MC*TT*hdampUP*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys HDAMPDown          ${merge_dir}/${nt}/${p}/MC*TT*hdampDOWN*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys TuneCUETP8M2T4Up   ${merge_dir}/${nt}/${p}/MC*TT*CUETP8M2T4up*.root
	python selections.py ${options} --out-dir distrs/  ${draw} --cond-com "${cond}" --distr ${distr} --chan ${chan} ${bins} --weight event_weight --sys TuneCUETP8M2T4Down ${merge_dir}/${nt}/${p}/MC*TT*CUETP8M2T4down*.root

distrs_mc_common_weights: p=pStage2Run1
distrs_mc_common_weights: nt=v25
distrs_mc_common_weights: dtag=MC
distrs_mc_common_weights: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_common_weights: cond=selection_stage == 107
distrs_mc_common_weights: chan=mu_sel
distrs_mc_common_weights: options=
distrs_mc_common_weights: merge_dir=merge-sets/
distrs_mc_common_weights:
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PUUp --distr Mt_lep_met_f    --chan ${chan}    --weight "event_weight*event_weight_PUUp/event_weight_PU"   ${merge_dir}/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PUDown --distr Mt_lep_met_f  --chan ${chan}  --weight "event_weight*event_weight_PUUp/event_weight_PU"     ${merge_dir}/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys bSFUp   --distr Mt_lep_met_f --chan ${chan} --weight "event_weight*event_weight_bSFUp/event_weight_bSF"    ${merge_dir}/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys bSFDown --distr Mt_lep_met_f --chan ${chan} --weight "event_weight*event_weight_bSFDown/event_weight_bSF"  ${merge_dir}/${nt}/${p}/*${dtag}*.root


distrs_mc_lep_weights: p=pStage2Run1
distrs_mc_lep_weights: nt=v25
distrs_mc_lep_weights: dtag=MC
distrs_mc_lep_weights: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_lep_weights: cond=selection_stage == 107
distrs_mc_lep_weights: chan=mu_sel
distrs_mc_lep_weights: options=
distrs_mc_lep_weights: sysname=LEP
distrs_mc_lep_weights: merge_dir=merge-sets/
distrs_mc_lep_weights:
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys ${sysname}IDUp    --distr Mt_lep_met_f --chan ${chan} --weight "event_weight*event_weight_LEPidUp"       ${merge_dir}/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys ${sysname}IDDown  --distr Mt_lep_met_f --chan ${chan} --weight "event_weight*event_weight_LEPidDown"     ${merge_dir}/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys ${sysname}TRGUp   --distr Mt_lep_met_f --chan ${chan} --weight "event_weight*event_weight_LEPtrgUp"      ${merge_dir}/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys ${sysname}TRGDown --distr Mt_lep_met_f --chan ${chan} --weight "event_weight*event_weight_LEPtrgDown"    ${merge_dir}/${nt}/${p}/*${dtag}*.root




distrs_mc_objects_mt: p=pStage2Run1
distrs_mc_objects_mt: nt=v25
distrs_mc_objects_mt: dtag=MC
distrs_mc_objects_mt: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_objects_mt: cond=selection_stage == 107
distrs_mc_objects_mt: chan=mu_sel
distrs_mc_objects_mt: options=
distrs_mc_objects_mt: met=event_met
distrs_mc_objects_mt: sys_name=NOMINAL
distrs_mc_objects_mt: weight=event_weight
distrs_mc_objects_mt: distr_name=Mt_lep_met_f
distrs_mc_objects_mt: merge_dir=merge-sets/
distrs_mc_objects_mt:
	python selections.py ${options} --out-dir distrs/ --distr ${distr_name} --chan ${chan} "sqrt(2*(sqrt((${met}.Px()*${met}.Px() + ${met}.Py()*${met}.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (${met}.Px()*event_leptons[0].Px() + ${met}.Py()*event_leptons[0].Py())))"    --cond-com "${cond}" ${bins} --sys ${sys_name}   --weight ${weight}  ${merge_dir}/${nt}/${p}/*${dtag}*.root

distrs_mc_objects: p=pStage2Run1
distrs_mc_objects: nt=v25
distrs_mc_objects: dtag=MC
distrs_mc_objects: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_objects: cond=selection_stage == 107
distrs_mc_objects: chan=mu_sel
distrs_mc_objects: options=
distrs_mc_objects: merge_dir=merge-sets/
distrs_mc_objects:
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_JERDown.Px()*event_met_JERDown.Px() + event_met_JERDown.Py()*event_met_JERDown.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_JERDown.Px()*event_leptons[0].Px() + event_met_JERDown.Py()*event_leptons[0].Py())))"    --cond-com "selection_stage_JERDown ${cond}" ${bins} --sys JERDown   --weight event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_JERUp.Px()*event_met_JERUp.Px() + event_met_JERUp.Py()*event_met_JERUp.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_JERUp.Px()*event_leptons[0].Px() + event_met_JERUp.Py()*event_leptons[0].Py())))"                --cond-com "selection_stage_JERUp   ${cond}" ${bins} --sys JERUp     --weight event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_JESDown.Px()*event_met_JESDown.Px() + event_met_JESDown.Py()*event_met_JESDown.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_JESDown.Px()*event_leptons[0].Px() + event_met_JESDown.Py()*event_leptons[0].Py())))"    --cond-com "selection_stage_JESDown ${cond}" ${bins} --sys JESDown   --weight event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_JESUp.Px()*event_met_JESUp.Px() + event_met_JESUp.Py()*event_met_JESUp.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_JESUp.Px()*event_leptons[0].Px() + event_met_JESUp.Py()*event_leptons[0].Py())))"                --cond-com "selection_stage_JESUp   ${cond}" ${bins} --sys JESUp     --weight event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_TESDown.Px()*event_met_TESDown.Px() + event_met_TESDown.Py()*event_met_TESDown.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_TESDown.Px()*event_leptons[0].Px() + event_met_TESDown.Py()*event_leptons[0].Py())))"    --cond-com "selection_stage_TESDown ${cond}" ${bins} --sys TESDown   --weight event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met_TESUp.Px()*event_met_TESUp.Px() + event_met_TESUp.Py()*event_met_TESUp.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met_TESUp.Px()*event_leptons[0].Px() + event_met_TESUp.Py()*event_leptons[0].Py())))"                --cond-com "selection_stage_TESUp   ${cond}" ${bins} --sys TESUp     --weight event_weight  ${merge_dir}/${nt}/${p}/*${dtag}*.root



distrs_mc_tt_weights: p=pStage2Run1
distrs_mc_tt_weights: nt=v25
distrs_mc_tt_weights: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_tt_weights: cond=selection_stage == 107
distrs_mc_tt_weights: chan=mu_sel
distrs_mc_tt_weights: options=
distrs_mc_tt_weights: merge_dir=merge-sets/
distrs_mc_tt_weights:
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys TOPPTDown   --weight event_weight                                ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys TOPPTUp     --weight "event_weight*event_weight_toppt"           ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys FragUp      --weight "event_weight*event_weight_FragUp"          ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys FragDown    --weight "event_weight*event_weight_FragDown"        ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys SemilepBRUp --weight "event_weight*event_weight_SemilepBRUp"     ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys SemilepBRDown --weight "event_weight*event_weight_SemilepBRDown" ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PetersonUp    --weight "event_weight*event_weight_PetersonUp"    ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PetersonDown  --weight "event_weight"                            ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root



distrs_mc_tt_weights_scale: p=pStage2Run1
distrs_mc_tt_weights_scale: nt=v25
distrs_mc_tt_weights_scale: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_tt_weights_scale: cond=selection_stage == 107
distrs_mc_tt_weights_scale: chan=mu_sel
distrs_mc_tt_weights_scale: options=
distrs_mc_tt_weights_scale: merge_dir=merge-sets/
distrs_mc_tt_weights_scale:
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys MrUp     --weight "event_weight*event_weight_me_f_rUp"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys MrDown   --weight "event_weight*event_weight_me_f_rDn"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys MfUp     --weight "event_weight*event_weight_me_fUp_r"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys MfDown   --weight "event_weight*event_weight_me_fDn_r"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys MfrUp    --weight "event_weight*event_weight_me_frUp"               ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys MfrDown  --weight "event_weight*event_weight_me_frDn"               ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root


distrs_mc_tt_weights_alphas: p=pStage2Run1
distrs_mc_tt_weights_alphas: nt=v25
distrs_mc_tt_weights_alphas: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_tt_weights_alphas: cond=selection_stage == 107
distrs_mc_tt_weights_alphas: chan=mu_sel
distrs_mc_tt_weights_alphas: options=
distrs_mc_tt_weights_alphas: merge_dir=merge-sets
distrs_mc_tt_weights_alphas:
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys AlphaSUp     --weight "event_weight*event_weight_AlphaS_up"          ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys AlphaSDown   --weight "event_weight*event_weight_AlphaS_dn"          ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root


# pdfs only up -- nominal for down
distrs_mc_tt_weights_pdf: p=pStage2Run1
distrs_mc_tt_weights_pdf: nt=v25
distrs_mc_tt_weights_pdf: bins=--custom-range 0,20,40,60,80,100,130,160,200,250
distrs_mc_tt_weights_pdf: cond=selection_stage == 107
distrs_mc_tt_weights_pdf: chan=mu_sel
distrs_mc_tt_weights_pdf: options=
distrs_mc_tt_weights_pdf: merge_dir=merge-sets
distrs_mc_tt_weights_pdf:
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n1Up     --weight "event_weight*event_weight_pdf[0]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n2Up     --weight "event_weight*event_weight_pdf[1]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n3Up     --weight "event_weight*event_weight_pdf[2]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n4Up     --weight "event_weight*event_weight_pdf[3]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n5Up     --weight "event_weight*event_weight_pdf[4]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n6Up     --weight "event_weight*event_weight_pdf[5]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n7Up     --weight "event_weight*event_weight_pdf[6]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n8Up     --weight "event_weight*event_weight_pdf[7]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n9Up     --weight "event_weight*event_weight_pdf[8]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n10Up    --weight "event_weight*event_weight_pdf[9]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n11Up    --weight "event_weight*event_weight_pdf[10]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n12Up    --weight "event_weight*event_weight_pdf[11]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n13Up    --weight "event_weight*event_weight_pdf[12]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n14Up    --weight "event_weight*event_weight_pdf[13]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n15Up    --weight "event_weight*event_weight_pdf[14]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n16Up    --weight "event_weight*event_weight_pdf[15]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n17Up    --weight "event_weight*event_weight_pdf[16]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n18Up    --weight "event_weight*event_weight_pdf[17]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n19Up    --weight "event_weight*event_weight_pdf[18]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n20Up    --weight "event_weight*event_weight_pdf[19]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n21Up    --weight "event_weight*event_weight_pdf[20]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n22Up    --weight "event_weight*event_weight_pdf[21]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n23Up    --weight "event_weight*event_weight_pdf[22]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n24Up    --weight "event_weight*event_weight_pdf[23]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n25Up    --weight "event_weight*event_weight_pdf[24]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n26Up    --weight "event_weight*event_weight_pdf[25]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n27Up    --weight "event_weight*event_weight_pdf[26]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n28Up    --weight "event_weight*event_weight_pdf[27]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n29Up    --weight "event_weight*event_weight_pdf[28]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n30Up    --weight "event_weight*event_weight_pdf[29]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n31Up    --weight "event_weight*event_weight_pdf[30]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n32Up    --weight "event_weight*event_weight_pdf[31]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n33Up    --weight "event_weight*event_weight_pdf[32]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n34Up    --weight "event_weight*event_weight_pdf[33]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n35Up    --weight "event_weight*event_weight_pdf[34]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n36Up    --weight "event_weight*event_weight_pdf[35]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n37Up    --weight "event_weight*event_weight_pdf[36]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n38Up    --weight "event_weight*event_weight_pdf[37]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n39Up    --weight "event_weight*event_weight_pdf[38]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n40Up    --weight "event_weight*event_weight_pdf[39]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n41Up    --weight "event_weight*event_weight_pdf[40]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n42Up    --weight "event_weight*event_weight_pdf[41]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n43Up    --weight "event_weight*event_weight_pdf[42]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n44Up    --weight "event_weight*event_weight_pdf[43]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n45Up    --weight "event_weight*event_weight_pdf[44]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n46Up    --weight "event_weight*event_weight_pdf[45]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n47Up    --weight "event_weight*event_weight_pdf[46]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n48Up    --weight "event_weight*event_weight_pdf[47]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n49Up    --weight "event_weight*event_weight_pdf[48]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n50Up    --weight "event_weight*event_weight_pdf[49]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n51Up    --weight "event_weight*event_weight_pdf[50]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n52Up    --weight "event_weight*event_weight_pdf[51]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n53Up    --weight "event_weight*event_weight_pdf[52]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n54Up    --weight "event_weight*event_weight_pdf[53]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n55Up    --weight "event_weight*event_weight_pdf[54]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root
	python selections.py ${options} --out-dir distrs/ --distr Mt_lep_met_f --chan ${chan} "sqrt(2*(sqrt((event_met.Px()*event_met.Px() + event_met.Py()*event_met.Py())*(event_leptons[0].Px()*event_leptons[0].Px() + event_leptons[0].Py()*event_leptons[0].Py())) - (event_met.Px()*event_leptons[0].Px() + event_met.Py()*event_leptons[0].Py())))" --cond-com "${cond}" ${bins} --sys PDFCT14n56Up    --weight "event_weight*event_weight_pdf[55]"              ${merge_dir}/${nt}/${p}/MC2016_Summer16_TTJets_powheg.root


distrs_usual_data_sel: merge_dir=lstore_outdirs/merge-sets/
distrs_usual_data_sel: nt=v40
distrs_usual_data_sel: p=p1
distrs_usual_data_sel: more=""
distrs_usual_data_sel: suff=""
distrs_usual_data_sel:
	# data
	make distrs_data_nominal merge_dir=${merge_dir} nt=${nt} p=${p} cond="selection_stage==107 ${more}"  chan="mu_sel${suff}"    &
	make distrs_data_nominal merge_dir=${merge_dir} nt=${nt} p=${p} cond="selection_stage==106 ${more}"  chan="mu_sel${suff}_ss" &
	make distrs_data_nominal merge_dir=${merge_dir} nt=${nt} p=${p} cond="(selection_stage==107 || selection_stage==105) ${more}"  chan="mu_selVloose${suff}"    &
	make distrs_data_nominal merge_dir=${merge_dir} nt=${nt} p=${p} cond="(selection_stage==106 || selection_stage==104) ${more}"  chan="mu_selVloose${suff}_ss" &
	#
	make distrs_data_nominal merge_dir=${merge_dir} nt=${nt} p=${p} options="--el-procs" cond="selection_stage==117 ${more}"  chan="el_sel${suff}"    dtag=SingleEle &
	make distrs_data_nominal merge_dir=${merge_dir} nt=${nt} p=${p} options="--el-procs" cond="selection_stage==116 ${more}"  chan="el_sel${suff}_ss" dtag=SingleEle &
	make distrs_data_nominal merge_dir=${merge_dir} nt=${nt} p=${p} options="--el-procs" cond="(selection_stage==117 || selection_stage==115) ${more}"  chan="el_selVloose${suff}"    dtag=SingleEle &
	make distrs_data_nominal merge_dir=${merge_dir} nt=${nt} p=${p} options="--el-procs" cond="(selection_stage==116 || selection_stage==114) ${more}"  chan="el_selVloose${suff}_ss" dtag=SingleEle &


distrs_usual_nominal_sel: merge_dir=lstore_outdirs/merge-sets/
distrs_usual_nominal_sel: nt=v40
distrs_usual_nominal_sel: p=p1
distrs_usual_nominal_sel: more=""
distrs_usual_nominal_sel: suff=""
distrs_usual_nominal_sel: proc=distrs_mc_nominal
distrs_usual_nominal_sel:
	# mc nominal
	make ${proc} merge_dir=${merge_dir} nt=${nt} p=${p} cond="selection_stage==107 ${more}"  chan="mu_sel${suff}"    &
	make ${proc} merge_dir=${merge_dir} nt=${nt} p=${p} cond="selection_stage==106 ${more}"  chan="mu_sel${suff}_ss" &
	make ${proc} merge_dir=${merge_dir} nt=${nt} p=${p} cond="(selection_stage==107 || selection_stage==105) ${more}"  chan="mu_selVloose${suff}"    &
	make ${proc} merge_dir=${merge_dir} nt=${nt} p=${p} cond="(selection_stage==106 || selection_stage==104) ${more}"  chan="mu_selVloose${suff}_ss" &
	
	make ${proc} merge_dir=${merge_dir} nt=${nt} p=${p} options="--el-procs" cond="selection_stage==117 ${more}"  chan="el_sel${suff}"    &
	make ${proc} merge_dir=${merge_dir} nt=${nt} p=${p} options="--el-procs" cond="selection_stage==116 ${more}"  chan="el_sel${suff}_ss" &
	make ${proc} merge_dir=${merge_dir} nt=${nt} p=${p} options="--el-procs" cond="(selection_stage==117 || selection_stage==115) ${more}"  chan="el_selVloose${suff}"    &
	make ${proc} merge_dir=${merge_dir} nt=${nt} p=${p} options="--el-procs" cond="(selection_stage==116 || selection_stage==114) ${more}"  chan="el_selVloose${suff}_ss" &



# some handy scripts
jobmaking: inp_dir=lstore_outdirs/v37/test1/MC2016_Summer16_TTJets_powheg/
jobmaking: template_file=template_file
jobmaking: template_name=INPUTFILE_PATH
jobmaking: job_dir=./
jobmaking:
	ls ${inp_dir} | grep .root | xargs -n 1 -I FOO sh -c 'sed "s,${template_name},${inp_dir}/FOO," ${template_file} > ${job_dir}/FOO.job'

# same with find
jobmaking2: inp_dir=lstore_outdirs/v37/test1/MC2016_Summer16_TTJets_powheg/
jobmaking2: template_file=template_file
jobmaking2: template_name=INPUTFILE_PATH
jobmaking2: job_dir=./
jobmaking2:
	find ${inp_dir} -name "*.root" -exec sh -c 'sed "s,${template_name},${inp_dir}/{}," ${template_file} > ${job_dir}/`basename {}`.job' \;

run_jobs_script: jobs_script=
run_jobs_script: max_procs=12
run_jobs_script:
	cat ${jobs_script} | xargs --delimiter='\n' --max-procs=${max_procs} -I COMMAND sh -c "COMMAND"

run_script_jobs: job_dir=
run_script_jobs: max_procs=12
run_script_jobs:
	ls ${job_dir} | xargs --max-procs=${max_procs} -n 1 -I FOO sh -c 'source ${job_dir}/FOO'

