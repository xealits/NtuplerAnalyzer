
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
run_histos:
	for sys in ${systematics}; \
	do \
	python stack_plot.py ${mc} ${data} -o ${out} --lumi ${lumi} -c ${chan} -d ${distr} --wjet ${wjets_norm} ${qcd_opt} ${shape_opt} --systematic $$sys ; \
	done

#run_histos: sys_file=syst-test.syst-list.1
run_std_histos: mc=distrs/distrs_v25v26resub2_mc_noqcd.root
run_std_histos: data_mu=../quick-test/v25v26-resub2_hadds/SingleMuon.root
run_std_histos: data_el=../quick-test/v25v26-resub2_hadds/SingleElectron.root
run_std_histos: out=histosels/redo-v25-v26_resub2/
run_std_histos: lumi_mu=35800
run_std_histos: lumi_el=31300
run_std_histos:
	make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} out=${out} chan=mu_sel_lj    systematics=NOMINAL qcd_opt="--qcd 1.05" shape_opt="-x mu_selVloose_lj" & \
	make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} out=${out} chan=el_sel_lj    systematics=NOMINAL qcd_opt="--qcd 1.05" shape_opt="-x el_selVloose_lj" & \
	make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} out=${out} chan=mu_sel_ljout systematics=NOMINAL qcd_opt="--qcd 1.05" shape_opt="-x mu_selVloose_ljout" & \
	make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} out=${out} chan=el_sel_ljout systematics=NOMINAL qcd_opt="--qcd 1.05" shape_opt="-x el_selVloose_ljout" & \
	make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} out=${out} chan=mu_sel_lj    systematics="`cat syst-test.syst-list.1.no_nominal syst-test.syst-list.2_pdfs | tr '\n' ' '`" shape_opt="-x mu_selVloose_lj" & \
	make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} out=${out} chan=el_sel_lj    systematics="`cat syst-test.syst-list.1.no_nominal syst-test.syst-list.2_pdfs | tr '\n' ' '`" shape_opt="-x el_selVloose_lj" & \
	make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} out=${out} chan=mu_sel_ljout systematics="`cat syst-test.syst-list.1.no_nominal syst-test.syst-list.2_pdfs | tr '\n' ' '`" shape_opt="-x mu_selVloose_ljout" & \
	make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} out=${out} chan=el_sel_ljout systematics="`cat syst-test.syst-list.1.no_nominal syst-test.syst-list.2_pdfs | tr '\n' ' '`" shape_opt="-x el_selVloose_ljout" & \
	make run_histos mc=${mc} data=${data_mu} lumi=${lumi_mu} out=${out} chan=mu_sel & \
	make run_histos mc=${mc} data=${data_el} lumi=${lumi_el} out=${out} chan=el_sel & 

