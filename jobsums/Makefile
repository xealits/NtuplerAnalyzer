
plots_sel: data_in=try2_split12/SingleMuon.root
plots_sel: mc_in=mc_mu_try2_split12.root
plots_sel: out_dir=./
plots_sel: sel=mu
plots_sel:
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met -s bSFUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met -s bSFDown
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met -s JESUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met -s JESDown
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met -s bSFUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met -s bSFDown
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met -s JESUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met -s JESDown
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d Mt_lep_met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d Mt_lep_met -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d dijet_mass
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d trijet_mass
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d dijet_trijet_mass
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d lep_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d lep_eta
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d tau_eta
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d tau_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d jet_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d jet_eta
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d bjet_eta
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d bjet_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d nvtx
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d nvtx -s PUUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_sel -d nvtx -s PUDown

plots_sel_all_sys: data_in=try2_split12/SingleMuon.root
plots_sel_all_sys: mc_in=mc_mu_try2_split12.root
plots_sel_all_sys: out_dir=./
plots_sel_all_sys: channel=mu_sel
plots_sel_all_sys:
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s NOMINAL     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s JESUp     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s JESDown   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s JERUp     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s JERDown   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s TauESUp   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s TauESDown 
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s bSFUp     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s bSFDown   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s PUUp      
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s PUDown    
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s TOPPTUp      
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${channel} -d Mt_lep_met -s TOPPTDown    

plots_sel_all_sys_fine: data_in=try2_split12/SingleMuon.root
plots_sel_all_sys_fine: mc_in=mc_mu_try2_split12.root
plots_sel_all_sys_fine: out_dir=./
plots_sel_all_sys_fine: sel=mu
plots_sel_all_sys_fine:
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s NOMINAL     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s JESUp     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s JESDown   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s JERUp     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s JERDown   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s TauESUp   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s TauESDown 
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s bSFUp     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s bSFDown   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s PUUp      
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s PUDown    
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s TOPPTUp      
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj -d Mt_lep_met_f -s TOPPTDown    
	# out
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s NOMINAL     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s JESUp     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s JESDown   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s JERUp     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s JERDown   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s TauESUp   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s TauESDown 
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s bSFUp     
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s bSFDown   
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s PUUp      
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s PUDown    
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s TOPPTUp      
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_lj_out -d Mt_lep_met_f -s TOPPTDown    

plots_presel: data_in=try2_split12/SingleMuon.root
plots_presel: mc_in=mc_mu_try2_split12.root
plots_presel: out_dir=./
plots_presel: sel=mu
plots_presel:
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d Mt_lep_met_f
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d Mt_lep_met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d Mt_lep_met -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d dijet_mass
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d trijet_mass
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d dijet_trijet_mass
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d lep_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d lep_pt -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d lep_eta
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d tau_eta
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d tau_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d tau_pt -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d jet_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d jet_pt -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d jet_eta
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d bjet_eta
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d bjet_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d bjet_pt -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d met -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d nvtx
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d nvtx -s PUUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ${sel}_presel -d nvtx -s PUDown


plots_misc: data_in=try2_split12/SingleMuon.root
plots_misc: mc_in=mc_mu_try2_split12.root
plots_misc: out_dir=./
plots_misc:
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c adv_mu_sel
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c mu_sel
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c mu_lj
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c mu_lj_out
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c mu_sel_ss
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_wjet
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_wjet_old
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_mumu
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_tt_em
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_tt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_tt_ss
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_wjet     -d met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_wjet_old -d met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_wjet     -d lep_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_wjet_old -d lep_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_tt_em -d M_lep_lep
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_tt_em -d met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_tt_em -d lep_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_tt_em -d Mt_lep_met -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_tt_em -d met        -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_tt_em -d lep_pt     -s TOPPTUp
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_mumu  -d met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_mumu  -d lep_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_mumu -d M_lep_lep
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_tt    -d lep_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_tt_ss -d lep_pt
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_tt    -d met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_tt_ss -d met
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_tt    -d M_lep_tau
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} -p -r -c ctr_mu_dy_tt_ss -d M_lep_tau

