# doc string
# make histosel
histosel: systs_files=syst-test.syst-list.1
histosel: n_proc=2
histosel: distrs_mc=distrs/distrs_v25v26resub2_mc_noqcd.root
histosel: distrs_data=../quick-test/v25v26-resub2_hadds/SingleMuon.root
histosel: out=histosels/redo-v25-v26_resub2_test-prallel
histosel: lumi=35800
histosel: chan=mu_sel
histosel: distr=Mt_lep_met_c
histosel: opt=
histosel:
	echo `date` histosel distrs_mc=${distrs_mc} distrs_data=${distrs_data} out=${out} chan=${chan} distr=${distr} systs_files=${systs_files} >> logs
	cat ${systs_files} | xargs -P ${n_proc} -I SYS sh -c "python stack_plot.py ${distrs_mc} ${distrs_data} -o ${out} --lumi ${lumi} -c ${chan} -d ${distr} ${opt} --wjet 0.8 --qcd 1.05 --systematic SYS || true"


plot_sys_vars: histosel_file=histosels/redo-v25-v26_1_PUnorm/HISTOSELS_v25v25redo_PUnormFix.root
plot_sys_vars: outdir=an-v7-and-paper-draft-3_answers/sys-vars/
plot_sys_vars:
	# all common systematics, except PDF (and alpha_s)
	for sysdef in PUUp,PUDown TOPPTUp,TOPPTDown bSFUp,bSFDown TESUp,TESDown JESUp,JESDown JERUp,JERDown FragUp,FragDown PetersonUp,PetersonDown SemilepBRUp,SemilepBRDown MfUp,MfDown MfrUp,MfrDown MrUp,MrDown ; \
	python mc_per_sys.py --ratio-range 0.1 --y-max 1500. ${histosel_file} -c el_sel_lj    -d Mt_lep_met_f -s $$sysdef --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " -o ${outdir} --outname el_sel_lj_all_Mt_lep_met_f_$$sysdef.png ; \
	python mc_per_sys.py --ratio-range 0.1 --y-max 1500. ${histosel_file} -c el_sel_ljout -d Mt_lep_met_f -s $$sysdef --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " -o ${outdir} --outname el_sel_ljout_all_Mt_lep_met_f_$$sysdef.png ; \
	python mc_per_sys.py --ratio-range 0.1 --y-max 2500. ${histosel_file} -c mu_sel_lj    -d Mt_lep_met_f -s $$sysdef --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " -o ${outdir} --outname mu_sel_lj_all_Mt_lep_met_f_$$sysdef.png ; \
	python mc_per_sys.py --ratio-range 0.1 --y-max 2500. ${histosel_file} -c mu_sel_ljout -d Mt_lep_met_f -s $$sysdef --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " -o ${outdir} --outname mu_sel_ljout_all_Mt_lep_met_f_$$sysdef.png ; \
	end
	# lep el systematics
	for sysdef in LEPelIDUp,LEPelIDDown LEPelTRGUp,LEPelTRGDown ; \
	python mc_per_sys.py --ratio-range 0.1 --y-max 1500 ${histosel_file} -c el_sel_lj    -d Mt_lep_met_f -s $$sysdef --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " -o an-v7-and-paper-draft-3_answers/sys-vars/ --outname el_sel_lj_all_Mt_lep_met_f_$$sysdef.png ; \
	python mc_per_sys.py --ratio-range 0.1 --y-max 1500 ${histosel_file} -c el_sel_ljout -d Mt_lep_met_f -s $$sysdef --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " -o an-v7-and-paper-draft-3_answers/sys-vars/ --outname el_sel_ljout_all_Mt_lep_met_f_$$sysdef.png ; \
	end
	# lep mu systematics
	for sysdef in LEPmuIDUp,LEPmuIDDown LEPmuTRGUp,LEPmuTRGDown ; \
	python mc_per_sys.py --ratio-range 0.1 --y-max 2500. ${histosel_file} -c mu_sel_lj    -d Mt_lep_met_f -s $$sysdef --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " -o an-v7-and-paper-draft-3_answers/sys-vars/ --outname mu_sel_lj_all_Mt_lep_met_f_$$sysdef.png ; \
	python mc_per_sys.py --ratio-range 0.1 --y-max 2500. ${histosel_file} -c mu_sel_ljout -d Mt_lep_met_f -s $$sysdef --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " -o an-v7-and-paper-draft-3_answers/sys-vars/ --outname mu_sel_ljout_all_Mt_lep_met_f_$$sysdef.png ; \
	end

plot_sys_vars_tau_ID_misID: distr_mc=distrs/distrs_v25v26resub2_mc_noqcd.root
plot_sys_vars_tau_ID_misID: distr_data=../quick-test/v25v26-resub2_hadds/SingleMuon.root
plot_sys_vars_tau_ID_misID: output=plots/distrs_v25v26resub2_sys_var_tauID_misID/
plot_sys_vars_tau_ID_misID: data_nick=data
plot_sys_vars_tau_ID_misID: lumi_mu=35800
plot_sys_vars_tau_ID_misID: lumi_el=31300
plot_sys_vars_tau_ID_misID:
	# plots with the variations
	# tauID
	# UP
	python stack_plot.py ${distr_mc} ${distr_data}  -o ${output}  -c mu_sel_lj    --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_mu} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauIDDown --factor-rate-systematic tt_mutau,tt_eltau,tt_taultauh,s_top_mutau,s_top_eltau,s_top_other,dy_tautau:0.95
	python stack_plot.py ${distr_mc} ${distr_data}  -o ${output}  -c mu_sel_ljout --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_mu} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauIDDown --factor-rate-systematic tt_mutau,tt_eltau,tt_taultauh,s_top_mutau,s_top_eltau,s_top_other,dy_tautau:0.95
	# UP
	python stack_plot.py ${distr_mc} ${distr_data}  -o ${output}  -c mu_sel_lj    --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_mu} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauIDUp   --factor-rate-systematic tt_mutau,tt_eltau,tt_taultauh,s_top_mutau,s_top_eltau,s_top_other,dy_tautau:1.05
	python stack_plot.py ${distr_mc} ${distr_data}  -o ${output}  -c mu_sel_ljout --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_mu} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauIDUp   --factor-rate-systematic tt_mutau,tt_eltau,tt_taultauh,s_top_mutau,s_top_eltau,s_top_other,dy_tautau:1.05
	#
	# same for el
	python stack_plot.py ${distr_mc} ${distr_data}  -o ${output}  -c el_sel_lj    --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_el} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauIDDown --factor-rate-systematic tt_mutau,tt_eltau,tt_taultauh,s_top_mutau,s_top_eltau,s_top_other,dy_tautau:0.95
	python stack_plot.py ${distr_mc} ${distr_data}  -o ${output}  -c el_sel_ljout --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_el} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauIDDown --factor-rate-systematic tt_mutau,tt_eltau,tt_taultauh,s_top_mutau,s_top_eltau,s_top_other,dy_tautau:0.95
	#
	python stack_plot.py ${distr_mc} ${distr_data}  -o ${output}  -c el_sel_lj    --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_el} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauIDUp   --factor-rate-systematic tt_mutau,tt_eltau,tt_taultauh,s_top_mutau,s_top_eltau,s_top_other,dy_tautau:1.05
	python stack_plot.py ${distr_mc} ${distr_data}  -o ${output}  -c el_sel_ljout --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_el} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauIDUp   --factor-rate-systematic tt_mutau,tt_eltau,tt_taultauh,s_top_mutau,s_top_eltau,s_top_other,dy_tautau:1.05
	# 
	# misID
	# DOWN
	python stack_plot.py  ${distr_mc} ${distr_data} -o ${output} -c mu_sel_lj    --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_mu} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauMisIDDown --factor-rate-systematic tt_lj,tt_taulj,tt_other,s_top_lj,dibosons:0.5
	python stack_plot.py  ${distr_mc} ${distr_data} -o ${output} -c mu_sel_ljout --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_mu} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauMisIDDown --factor-rate-systematic tt_lj,tt_taulj,tt_other,s_top_lj,dibosons:0.5
	# UP
	python stack_plot.py  ${distr_mc} ${distr_data} -o ${output} -c mu_sel_lj    --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_mu} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauMisIDUp   --factor-rate-systematic tt_lj,tt_taulj,tt_other,s_top_lj,dibosons:1.5
	python stack_plot.py  ${distr_mc} ${distr_data} -o ${output} -c mu_sel_ljout --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_mu} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauMisIDUp   --factor-rate-systematic tt_lj,tt_taulj,tt_other,s_top_lj,dibosons:1.5
	# for el
	python stack_plot.py  ${distr_mc} ${distr_data} -o ${output} -c el_sel_lj    --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_el} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauMisIDDown --factor-rate-systematic tt_lj,tt_taulj,tt_other,s_top_lj,dibosons:0.5
	python stack_plot.py  ${distr_mc} ${distr_data} -o ${output} -c el_sel_ljout --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_el} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauMisIDDown --factor-rate-systematic tt_lj,tt_taulj,tt_other,s_top_lj,dibosons:0.5
	#
	python stack_plot.py  ${distr_mc} ${distr_data} -o ${output} -c el_sel_lj    --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_el} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauMisIDUp   --factor-rate-systematic tt_lj,tt_taulj,tt_other,s_top_lj,dibosons:1.5
	python stack_plot.py  ${distr_mc} ${distr_data} -o ${output} -c el_sel_ljout --data-nick ${data_nick} -d Mt_lep_met_f --bin-norm --lumi ${lumi_el} --qcd 1.05 --overwrite --wjets 0.8  --rename-sys TauMisIDUp   --factor-rate-systematic tt_lj,tt_taulj,tt_other,s_top_lj,dibosons:1.5
	#
	hadd ${output}/distrs/histosel.root ${output}/distrs/distrs*root
	#
	python mc_per_sys.py --y-max 80.  ${output}/distrs/histosel.root -c el_sel_lj    -d Mt_lep_met_f -s TauIDUp,TauIDDown --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " --outname el_sel_lj_all_Mt_lep_met_f_TauIDUp,TauIDDown.png
	python mc_per_sys.py --y-max 80.  ${output}/distrs/histosel.root -c el_sel_ljout -d Mt_lep_met_f -s TauIDUp,TauIDDown --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " --outname el_sel_ljout_all_Mt_lep_met_f_TauIDUp,TauIDDown.png
	python mc_per_sys.py --y-max 150. ${output}/distrs/histosel.root -c mu_sel_lj    -d Mt_lep_met_f -s TauIDUp,TauIDDown --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " --outname mu_sel_lj_all_Mt_lep_met_f_TauIDUp,TauIDDown.png
	python mc_per_sys.py --y-max 150. ${output}/distrs/histosel.root -c mu_sel_ljout -d Mt_lep_met_f -s TauIDUp,TauIDDown --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " --outname mu_sel_ljout_all_Mt_lep_met_f_TauIDUp,TauIDDown.png
	#
	python mc_per_sys.py --y-max 80.  ${output}/distrs/histosel.root -c el_sel_lj    -d Mt_lep_met_f -s TauMisIDUp,TauMisIDDown --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " --outname el_sel_lj_all_Mt_lep_met_f_TauMisIDUp,TauMisIDDown.png
	python mc_per_sys.py --y-max 80.  ${output}/distrs/histosel.root -c el_sel_ljout -d Mt_lep_met_f -s TauMisIDUp,TauMisIDDown --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " --outname el_sel_ljout_all_Mt_lep_met_f_TauMisIDUp,TauMisIDDown.png
	python mc_per_sys.py --y-max 150. ${output}/distrs/histosel.root -c mu_sel_lj    -d Mt_lep_met_f -s TauMisIDUp,TauMisIDDown --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " --outname mu_sel_lj_all_Mt_lep_met_f_TauMisIDUp,TauMisIDDown.png
	python mc_per_sys.py --y-max 150. ${output}/distrs/histosel.root -c mu_sel_ljout -d Mt_lep_met_f -s TauMisIDUp,TauMisIDDown --no-data --ratio  --add-titles --title-x "M_{T} [GeV]" --title " " --outname mu_sel_ljout_all_Mt_lep_met_f_TauMisIDUp,TauMisIDDown.png




hadd_distrs: nt=v37
hadd_distrs: proc=test1
hadd_distrs: distrs=bunch2
hadd_distrs: dtag=MC2016_Summer16_TTJets_powheg
hadd_distrs: max_procs=4
hadd_distrs:
	ls distrs/${nt}_${proc}_${distrs}/all/${dtag}/ | xargs -P ${max_procs} -I CHANS sh -c "hadd distrs/${nt}_${proc}_${distrs}/all/${dtag}/CHANS.root distrs/${nt}_${proc}_${distrs}/all/${dtag}/CHANS/*root"

multihadd2_distrs: nt=v37
multihadd2_distrs: proc=test1
multihadd2_distrs: distrs=bunch2
multihadd2_distrs: dtag=MC2016_Summer16_TTJets_powheg
multihadd2_distrs: size=30
multihadd2_distrs: max_procs=10
multihadd2_distrs:
	python multihadd2.py -i ${size} distrs/${nt}_${proc}_${distrs}/all/${dtag}/all_chans.root `find distrs/${nt}_${proc}_${distrs}/all/${dtag}/ -name "*root"` > temp.multihadd2_coms
	cat temp.multihadd2_coms | xargs -P ${max_procs} -d "\n" -I COMM sh -c "COMM || true"
	hadd distrs/${nt}_${proc}_${distrs}/all/${dtag}.root distrs/${nt}_${proc}_${distrs}/all/${dtag}/all_chans.root.*

multihadd2_distrs_all_dtags: nt=v37
multihadd2_distrs_all_dtags: proc=test1
multihadd2_distrs_all_dtags: distrs=bunch2
multihadd2_distrs_all_dtags: size=30
multihadd2_distrs_all_dtags: max_procs=10
multihadd2_distrs_all_dtags:
	for dtag in `ls distrs/${nt}_${proc}_${distrs}/all/`; \
	do \
	make multihadd2_distrs dtag=$$dtag nt=${nt} proc=${proc} distrs=${distrs} size=${size} max_procs=${max_procs} ; \
	done

# histosels for the full fit

# hadd ${merge_dir}/${nt}/${proc}/$$d.root /lstore/cms/olek/outdirs/${nt}/${proc}/$$d/*.root &
run_histos: systematics=NOMINAL
run_histos: mc=distrs/distrs_v25v26resub2_mc_noqcd.root
run_histos: data=../quick-test/v25v26-resub2_hadds/SingleMuon.root
run_histos: out=histosels/redo-v25-v26_resub2/
run_histos: lumi=35800
run_histos: chan=mu_sel
run_histos: distr=Mt_lep_met_f
run_histos: wjets_norm=0.8
run_histos: qcd_opt=
run_histos: shape_opt=
run_histos: data_nick=data_other
run_histos:
	for sys in ${systematics}; \
	do \
	python stack_plot.py ${mc} ${data} -o ${out} --data-nick ${data_nick} --lumi ${lumi} -c ${chan} -d ${distr} --wjet ${wjets_norm} ${qcd_opt} ${shape_opt} --systematic $$sys ; \
	done

#run_histos: sys_file=syst-test.syst-list.1
run_std_histos: mc=distrs/distrs_v25v26resub2_mc_noqcd.root
run_std_histos: data_mu=../quick-test/v25v26-resub2_hadds/SingleMuon.root
run_std_histos: data_el=../quick-test/v25v26-resub2_hadds/SingleElectron.root
run_std_histos: out=histosels/redo-v25-v26_resub2/
run_std_histos: distr=Mt_lep_met_c
run_std_histos: lumi_mu=35800
run_std_histos: lumi_el=31300
run_std_histos: wjets_norm=0.8
run_std_histos:
	# log parameters
	echo mc=${mc} data_mu=${data_mu} data_el=${data_el} lumi_mu=${lumi_mu} lumi_el=${lumi_el} distr=${distr} out=${out} wjets_norm=${wjets_norm} > ${out}/run_std_histos.params
	# the batch with data-driven QCD
	# the qcd has the experimental systematics, but not the TT modeling
	#make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=mu_sel_lj    systematics="NOMINAL `cat syst-test.syst-list.experimental | tr '\n' ' '`" qcd_opt="--qcd 1.05" shape_opt="-x mu_selVloose_lj" & \
	#make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=el_sel_lj    systematics="NOMINAL `cat syst-test.syst-list.experimental | tr '\n' ' '`" qcd_opt="--qcd 1.05" shape_opt="-x el_selVloose_lj" & \
	#make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=mu_sel_ljout systematics="NOMINAL `cat syst-test.syst-list.experimental | tr '\n' ' '`" qcd_opt="--qcd 1.05" shape_opt="-x mu_selVloose_ljout" & \
	#make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=el_sel_ljout systematics="NOMINAL `cat syst-test.syst-list.experimental | tr '\n' ' '`" qcd_opt="--qcd 1.05" shape_opt="-x el_selVloose_ljout" &
	make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=mu_sel_lj    systematics="NOMINAL" qcd_opt="--qcd 1.05" shape_opt="-x mu_selVloose_lj" & \
	make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=el_sel_lj    systematics="NOMINAL" qcd_opt="--qcd 1.05" shape_opt="-x el_selVloose_lj" & \
	make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=mu_sel_ljout systematics="NOMINAL" qcd_opt="--qcd 1.05" shape_opt="-x mu_selVloose_ljout" & \
	make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=el_sel_ljout systematics="NOMINAL" qcd_opt="--qcd 1.05" shape_opt="-x el_selVloose_ljout" & \
	# tt modeling with no data-driven qcd
	#make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} distr=${distr} out=${out} chan=mu_sel_lj    systematics="`cat syst-test.syst-list.parton_showers syst-test.syst-list.updowns syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x mu_selVloose_lj" &
	#make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} distr=${distr} out=${out} chan=el_sel_lj    systematics="`cat syst-test.syst-list.parton_showers syst-test.syst-list.updowns syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x el_selVloose_lj" & 
	#make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} distr=${distr} out=${out} chan=mu_sel_ljout systematics="`cat syst-test.syst-list.parton_showers syst-test.syst-list.updowns syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x mu_selVloose_ljout" &
	#make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} distr=${distr} out=${out} chan=el_sel_ljout systematics="`cat syst-test.syst-list.parton_showers syst-test.syst-list.updowns syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x el_selVloose_ljout" & 
	# no updowns now
	#make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=mu_sel_lj    systematics="`cat syst-test.syst-list.parton_showers syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x mu_selVloose_lj" & \
	#make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=el_sel_lj    systematics="`cat syst-test.syst-list.parton_showers syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x el_selVloose_lj" & \
	#make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=mu_sel_ljout systematics="`cat syst-test.syst-list.parton_showers syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x mu_selVloose_ljout" & \
	#make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=el_sel_ljout systematics="`cat syst-test.syst-list.parton_showers syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x el_selVloose_ljout" &
	make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=mu_sel_lj    systematics="`cat syst-test.syst-list.experimental syst-test.syst-list.parton_showers syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x mu_selVloose_lj" & \
	make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=el_sel_lj    systematics="`cat syst-test.syst-list.experimental syst-test.syst-list.parton_showers syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x el_selVloose_lj" & \
	make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=mu_sel_ljout systematics="`cat syst-test.syst-list.experimental syst-test.syst-list.parton_showers syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x mu_selVloose_ljout" & \
	make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=el_sel_ljout systematics="`cat syst-test.syst-list.experimental syst-test.syst-list.parton_showers syst-test.syst-list.pdfs | tr '\n' ' '`" shape_opt="-x el_selVloose_ljout" & \
	make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=mu_sel & \
	make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} distr=${distr} out=${out} wjets_norm=${wjets_norm} chan=el_sel & 



# all the plotting for the AN and paper

plot_control_dy_mumu:
	echo `date` plot_control_dy_mumu nt=${nt} proc=${proc} distrs=${distrs} chan=dy_mumu >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	# lep_pt, Mt, dilep mass
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c dy_mumu -d lep_pt         --wjet 0.8  --ratio-range 0.3 --title "" --title-x "leading muon p_{T} [GeV]"             --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c dy_mumu -d Mt_lep_met_c   --wjet 0.8  --ratio-range 0.3 --title "" --title-x "M_{T}(#mu, MET) [GeV]"    --bin-norm  --title-y "[Events/GeV]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c dy_mumu -d dilep_mass_dy  --wjet 0.8  --ratio-range 0.3 --title "" --title-x "mass(#mu + #mu) [GeV]"                --title-y "[Events]"

plot_control_dy_elel:
	echo `date` plot_control_dy_elel nt=${nt} proc=${proc} distrs=${distrs} chan=dy_elel >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	# lep_pt, Mt, dilep mass
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c dy_elel -d lep_pt         --wjet 0.8  --ratio-range 0.3 --title "" --title-x "leading electron p_{T} [GeV]"      --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c dy_elel -d Mt_lep_met_c   --wjet 0.8  --ratio-range 0.3 --title "" --title-x "M_{T}(e, MET) [GeV]"    --bin-norm --title-y "[Events/GeV]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c dy_elel -d dilep_mass_dy  --wjet 0.8  --ratio-range 0.3 --title "" --title-x "mass(e + e) [GeV]"                 --title-y "[Events]"

plot_control_wjets_mu:
	echo `date` plot_control_wjets_mu nt=${nt} proc=${proc} distrs=${distrs} chan=wjets_mu >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c wjets_mu -d lep_pt         --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --title "" --title-x "leading muon p_{T} [GeV]"            --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c wjets_mu -d Mt_lep_met_c   --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --title "" --title-x "M_{T}(#mu, MET) [GeV]"    --bin-norm --title-y "[Events/GeV]"

plot_control_wjets_el:
	echo `date` plot_control_wjets_el nt=${nt} proc=${proc} distrs=${distrs} chan=wjets_el >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c wjets_el -d lep_pt         --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --title "" --title-x "leading electron p_{T} [GeV]"      --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c wjets_el -d Mt_lep_met_c   --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --title "" --title-x "M_{T}(e, MET) [GeV]"    --bin-norm --title-y "[Events/GeV]"

# special wjets with amcatnlo
plot_control_wjets_amcatnlo:
	echo `date` plot_control_wjets_amcatnlo nt=${nt} proc=${proc} distrs=${distrs} chan=wjets_el wjets_mu >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	#
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc_wjetsamcatnlo.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c wjets_mu -d lep_pt         --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --title "" --title-x "leading muon p_{T} [GeV]"            --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc_wjetsamcatnlo.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c wjets_mu -d Mt_lep_met_c   --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --title "" --title-x "M_{T}(#mu, MET) [GeV]"    --bin-norm --title-y "[Events/GeV]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc_wjetsamcatnlo.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c wjets_el -d lep_pt         --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --title "" --title-x "leading electron p_{T} [GeV]"      --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc_wjetsamcatnlo.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c wjets_el -d Mt_lep_met_c   --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --title "" --title-x "M_{T}(e, MET) [GeV]"    --bin-norm --title-y "[Events/GeV]"


plot_lj_1: chan_tag=mu
plot_lj_1: mW=80.385
plot_lj_1: mt=172.5
plot_lj_1:
	echo `date` plot_lj_1 nt=${nt} proc=${proc} distrs=${distrs} >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	# Comparison qcd data-driven VS MC in wjets and presel control regions
	# wjets
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_wjets -d Mt_lep_met_c --wjet 0.8 --qcd        1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"   --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_wjets -d Mt_lep_met_c --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"   --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_wjets -d lep_pt --wjet 0.8 --qcd        1.05 --ratio-range 0.3                  --title "" --title-x "lepton p_{T} [GeV]" --title-y "[Events]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_wjets -d lep_pt --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3                  --title "" --title-x "lepton p_{T} [GeV]" --title-y "[Events]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_wjets -d met_pt --wjet 0.8 --qcd        1.05 --ratio-range 0.3       --bin-norm --title "" --title-x "E^{miss}_{T} [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_wjets -d met_pt --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3       --bin-norm --title "" --title-x "E^{miss}_{T} [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	# presel
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d Mt_lep_met_c --wjet 0.8 --qcd        1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"   --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d Mt_lep_met_c --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"   --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d lep_pt --wjet 0.8 --qcd        1.05 --ratio-range 0.3                  --title "" --title-x "lepton p_{T} [GeV]" --title-y "[Events]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d lep_pt --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3                  --title "" --title-x "lepton p_{T} [GeV]" --title-y "[Events]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d met_pt --wjet 0.8 --qcd        1.05 --ratio-range 0.3       --bin-norm --title "" --title-x "E^{miss}_{T} [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d met_pt --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3       --bin-norm --title "" --title-x "E^{miss}_{T} [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	# lj parameters, with MC and data-driven QCD
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d dijet_mass        --wjet 0.8 --qcd        1.05 --ratio-range 0.3       --bin-norm --title "" --title-x  "di-jet mass [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d trijet_mass       --wjet 0.8 --qcd        1.05 --ratio-range 0.3       --bin-norm --title "" --title-x "tri-jet mass [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d dijet_trijet_mass --wjet 0.8 --qcd        1.05 --ratio-range 0.3       --bin-norm --title "" --title-x "#sqrt{(m_{2j} - ${mW})^{2} + (m_{3j} - ${mt})^{2}} [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d dijet_mass        --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3       --bin-norm --title "" --title-x  "di-jet mass [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d trijet_mass       --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3       --bin-norm --title "" --title-x "tri-jet mass [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_preselCand -d dijet_trijet_mass --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3       --bin-norm --title "" --title-x "#sqrt{(m_{2j} - ${mW})^{2} + (m_{3j} - ${mt})^{2}} [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02

plot_tt_dilep:
	echo `date` plot_tt_dilep nt=${nt} proc=${proc} distrs=${distrs} >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	# plots for control of the shapes of single-top and dibosons
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c tt_elmu -d Mt_lep_met_c --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"   --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c tt_elmu -d lep_pt       --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3            --title "" --title-x "lepton p_{T} [GeV]" --title-y "[Events]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c tt_elmu -d met_pt       --wjet 0.8 --qcd-factor 1.5  --ratio-range 0.3 --bin-norm --title "" --title-x "E^{miss}_{T} [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02

plot_tt_dilep_mu: datafile=SingleMuon
plot_tt_dilep_mu: lumi=35800
plot_tt_dilep_mu: plot_tt_dilep


#plot_sel_1: nt=v25v26
#plot_sel_1: proc=redo
#plot_sel_1: distrs=4
#plot_sel_1: lumi=31300
#plot_sel_1: datafile=SingleElectron
#plot_sel_1: chan_tag=el
#plot_sel_1: title_lep=e

plot_presel_1: presel_lev=
plot_presel_1:
	echo `date` plot_presel_1 nt=${nt} proc=${proc} distrs=${distrs} chan=${chan_tag}_presel${presel_lev} >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	# main preselection
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_presel${presel_lev}      -d lep_pt  --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "${title_lep} p_{T} [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_presel${presel_lev}      -d lep_eta --wjet 0.8 --qcd 1.05 --ratio-range 0.3            --title "" --title-x "${title_lep} #eta" --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_presel${presel_lev}      -d tau_pt  --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "#tau p_{T} [GeV]"  --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_presel${presel_lev}      -d tau_eta --wjet 0.8 --qcd 1.05 --ratio-range 0.3            --title "" --title-x "#tau #eta"         --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_presel${presel_lev}      -d bjet_pt  --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "leading b-tagged jet p_{T} [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_presel${presel_lev}      -d bjet_eta --wjet 0.8 --qcd 1.05 --ratio-range 0.3            --title "" --title-x "leading b-tagged jet #eta"        --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_presel${presel_lev}      -d Mt_lep_met_c --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"         --title-y "[Events/GeV]"     --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_presel${presel_lev}      -d met_pt      --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "MET [GeV]"         --title-y "[Events/GeV]"     --draw-overflows 5 --legend-shift 0.02

# TODO: add the presel lj parameters, 2-j and 3-j masses, with Data and MC QCD
plot_presel_lj_1: presel_lev=
plot_presel_lj_1: the_distr=Mt_lep_met_f
plot_presel_lj_1:
	echo `date` plot_presel_lj_1 nt=${nt} proc=${proc} distrs=${distrs} chan=${chan_tag}_presel${presel_lev} lj-ljout ${the_distr} >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	# event categories pre-fit plots
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_presel${presel_lev}_ljout -d ${the_distr} --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"   --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_presel${presel_lev}_lj    -d ${the_distr} --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"   --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02


# the prefit Mt distrs in categories
plot_sel_lj_1: the_distr=Mt_lep_met_f
plot_sel_lj_1:
	echo `date` plot_sel_lj_1 pre-fit nt=${nt} proc=${proc} distrs=${distrs} chan=${chan_tag}_sel lj-ljout ${the_distr} >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	# event categories pre-fit plots
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel_ljout -x ${chan_tag}_selVloose_ljout -d Mt_lep_met_c --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"   --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel_lj    -x ${chan_tag}_selVloose_lj    -d Mt_lep_met_c --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"   --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02



# main selection plots without categories
plot_sel_1:
	echo `date` plot_sel_1 nt=${nt} proc=${proc} distrs=${distrs} chan=${chan_tag}_sel >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	# main selection
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d lep_pt  --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "${title_lep} p_{T} [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d lep_eta --wjet 0.8 --qcd 1.05 --ratio-range 0.3            --title "" --title-x "${title_lep} #eta" --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d tau_pt  --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "#tau p_{T} [GeV]"  --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d tau_eta --wjet 0.8 --qcd 1.05 --ratio-range 0.3            --title "" --title-x "#tau #eta"         --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d bjet_pt  --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "leading b-tagged jet p_{T} [GeV]" --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d bjet_eta --wjet 0.8 --qcd 1.05 --ratio-range 0.3            --title "" --title-x "leading b-tagged jet #eta"        --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d Mt_lep_met_c --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"         --title-y "[Events/GeV]"     --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d met_pt      --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "MET [GeV]"         --title-y "[Events/GeV]"     --draw-overflows 5 --legend-shift 0.02




plot_sel_2:
	echo `date` plot_sel_2 nt=${nt} proc=${proc} distrs=${distrs} chan=${chan_tag}_sel >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	# N jets
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d nbjets      --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --title "" --title-x "N b-tagged jets"      --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d nrjets      --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --title "" --title-x "N not b-tagged jets"  --title-y "[Events]"
	# N vertices
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d nvtx      --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --title "Nominal" --title-x "N reco vertices"  --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d nvtx      --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --title "PU Up"   --title-x "N reco vertices"  --title-y "[Events]" --systematic PUUp
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d nvtx      --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --title "PU Down" --title-x "N reco vertices"  --title-y "[Events]" --systematic PUDown



plot_sel_mt_1:
	echo `date` plot_sel_mt_1 nt=${nt} proc=${proc} distrs=${distrs} chan=${chan_tag}_sel >> logs
	mkdir -p plots/${nt}_${proc}_${distrs}/
	# Mt investigation
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d cos_phi_met_lep  --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --title "" --title-x "cos(${title_lep}, MET)"         --title-y "[Events]"
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel    -x ${chan_tag}_selVloose    -d phi_met_lep      --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --title "" --title-x "#phi(${title_lep}, MET)"        --title-y "[Events]"


plot_sel_std_1_mu: nt=v25v26
plot_sel_std_1_mu: proc=redo
plot_sel_std_1_mu: distrs=4
plot_sel_std_1_mu: lumi=35800
plot_sel_std_1_mu: datafile=SingleMuon
plot_sel_std_1_mu: chan_tag=mu
plot_sel_std_1_mu: title_lep=#mu
plot_sel_std_1_mu: plot_sel_1 plot_sel_lj_1

plot_sel_std_1_el: nt=v25v26
plot_sel_std_1_el: proc=redo
plot_sel_std_1_el: distrs=4
plot_sel_std_1_el: lumi=31300
plot_sel_std_1_el: datafile=SingleElectron
plot_sel_std_1_el: chan_tag=el
plot_sel_std_1_el: title_lep=e
plot_sel_std_1_el: plot_sel_1 plot_sel_lj_1

# this does not work on its own with defaults and dependencies
# # it runs only 1 dependency
# hence the explicit re-make
plot_sel: nt=v25v26
plot_sel: proc=redo
plot_sel: distrs=4
plot_sel:
	make plot_sel_std_1_mu nt=${nt} proc=${proc} distrs=${distrs}
	make plot_sel_std_1_el nt=${nt} proc=${proc} distrs=${distrs}


test_controls_v37_test5_issue1: tests=tests1
test_controls_v37_test5_issue1: mc=distrs/v37_test5_bunch2_mc_wjets_wNjets.root
test_controls_v37_test5_issue1: data_dir=distrs/v37_test5_bunch2/
test_controls_v37_test5_issue1:
	mkdir -p plots/v37_test5_bunch2_mc_${tests}/
	# presel
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_presel       --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.8 --qcd-factor 1.5
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_presel       --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.8 --qcd-factor 1.5
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_presel_lj    --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.8 --qcd-factor 1.5
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_presel_lj    --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.8 --qcd-factor 1.5
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_presel_ljout --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.8 --qcd-factor 1.5
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_presel_ljout --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.8 --qcd-factor 1.5
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_preselCand   --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.8 --qcd-factor 1.5
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_preselCand   --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.8 --qcd-factor 1.5
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_preselCand   --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.8 --qcd 1.05
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_preselCand   --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.8 --qcd 1.05
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_preselCand   --data-nick data_other -d Mt_lep_met_c2 --lumi 35800 --overwrite --wjets 0.8 --qcd-factor 1.5
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_preselCand   --data-nick data_other -d Mt_lep_met_c2 --lumi 35800 --overwrite --wjets 0.8 --qcd-factor 1.5
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_preselCand   --data-nick data_other -d Mt_lep_met_c2 --lumi 35800 --overwrite --wjets 0.8 --qcd 1.05
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_preselCand   --data-nick data_other -d Mt_lep_met_c2 --lumi 35800 --overwrite --wjets 0.8 --qcd 1.05
	# sel
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_sel       --data-nick data_other -d Mt_lep_met_c --lumi 35800 --overwrite --wjets 0.92 --qcd 1.05
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_sel       --data-nick data_other -d Mt_lep_met_c --lumi 35800 --overwrite --wjets 0.92 --qcd 1.05
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_sel       --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd 1.05
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_sel       --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd 1.05
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_sel       --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd 1.05
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_sel       --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd 1.05
	# lj categories
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_sel_lj    --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 0.0001
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_sel_lj    --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 0.0001
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c mu_sel_ljout --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 0.0001
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c el_sel_ljout --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 0.0001
	# other
	# tt_elmu
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c tt_elmu --data-nick data_other -d Mt_lep_met_f --lumi 35800 --overwrite --wjets 0.8
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c tt_elmu --data-nick data_other -d Mt_lep_met_f --lumi 35800 --overwrite --wjets 0.8
	# check the effect of top pT
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c tt_elmu --data-nick data_other -d Mt_lep_met_f --lumi 35800 --overwrite --wjets 0.8 --systematic TOPPTUp
	# the source of the difference in the shape
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c tt_elmu --data-nick data_other -d lep_pt  --lumi 35800 --overwrite --wjets 0.8
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c tt_elmu --data-nick data_other -d lep_pt  --lumi 35800 --overwrite --wjets 0.8
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c tt_elmu --data-nick data_other -d lep_eta --lumi 35800 --overwrite --wjets 0.8
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c tt_elmu --data-nick data_other -d lep_eta --lumi 35800 --overwrite --wjets 0.8
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c tt_elmu --data-nick data_other -d met_c   --lumi 35800 --overwrite --wjets 0.8
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c tt_elmu --data-nick data_other -d met_c   --lumi 35800 --overwrite --wjets 0.8
	# dy
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c dy_mumu --data-nick data_other -d Mt_lep_met_c --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 2.5
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c dy_elel --data-nick data_other -d Mt_lep_met_c --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 2.5
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c dy_mumu --data-nick data_other -d met_c --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 2.5
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c dy_elel --data-nick data_other -d met_c --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 2.5
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c dy_mumu --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 2.5
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c dy_elel --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 2.5
	# wjets
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c wjets_mu --data-nick data_other -d Mt_lep_met_c --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 2.5
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c wjets_el --data-nick data_other -d Mt_lep_met_c --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 2.5
	python stack_plot.py ${mc} ${data_dir}/SingleMuon.root     -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c wjets_mu --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 2.5
	python stack_plot.py ${mc} ${data_dir}/SingleElectron.root -p -r  -o plots/v37_test5_bunch2_mc_${tests}/  -c wjets_el --data-nick data_other -d lep_pt --lumi 35800 --overwrite --wjets 0.92 --qcd-factor 2.5




