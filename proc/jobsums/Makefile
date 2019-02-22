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
	cat ${systs_files} | xargs -P ${n_proc} -I SYS sh -c "python stack_plot.py ${distrs_mc} ${distrs_data} -o ${out} --lumi ${lumi} -c ${chan} -d ${distr} ${opt} --wjet 0.8 --qcd 1.05 --systematic SYS || true"


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

#plot_sel_1: nt=v25v26
#plot_sel_1: proc=redo
#plot_sel_1: distrs=4
#plot_sel_1: lumi=31300
#plot_sel_1: datafile=SingleElectron
#plot_sel_1: chan_tag=el
#plot_sel_1: title_lep=e


plot_sel_1:
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel_ljout -x ${chan_tag}_selVloose_ljout -d Mt_lep_met_f --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"   --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02
	python stack_plot.py distrs/${nt}_${proc}_${distrs}_mc.root distrs/${nt}_${proc}_${distrs}/${datafile}.root -o plots/${nt}_${proc}_${distrs}/ -p -r --lumi ${lumi} -c ${chan_tag}_sel_lj    -x ${chan_tag}_selVloose_lj    -d Mt_lep_met_f --wjet 0.8 --qcd 1.05 --ratio-range 0.3 --bin-norm --title "" --title-x "M_{T}(${title_lep}, MET) [GeV]"   --title-y "[Events/GeV]" --draw-overflows 5 --legend-shift 0.02

plot_sel1_mu: nt=v25v26
plot_sel1_mu: proc=redo
plot_sel1_mu: distrs=4
plot_sel1_mu: lumi=35800
plot_sel1_mu: datafile=SingleMuon
plot_sel1_mu: chan_tag=mu
plot_sel1_mu: title_lep=#mu
plot_sel1_mu: plot_sel_1

plot_sel1_el: nt=v25v26
plot_sel1_el: proc=redo
plot_sel1_el: distrs=4
plot_sel1_el: lumi=31300
plot_sel1_el: datafile=SingleElectron
plot_sel1_el: chan_tag=el
plot_sel1_el: title_lep=e
plot_sel1_el: plot_sel_1

# this does not work on its own with defaults and dependencies
# # it runs only 1 dependency
# hence the explicit re-make
plot_sel: nt=v25v26
plot_sel: proc=redo
plot_sel: distrs=4
plot_sel:
	make plot_sel1_mu nt=${nt} proc=${proc} distrs=${distrs}
	make plot_sel1_el nt=${nt} proc=${proc} distrs=${distrs}

