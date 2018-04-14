
plots_sel: out_dir=./
plots_sel: sel=mu
plots_sel: tag=mu_sel
plots_sel: qcd=1.
plots_sel: lumi=35.8
plots_sel:
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_Mt-lep-met_f -p -r --qcd ${qcd} --title "" --title-x "M_{T}(l, MET)"  -c ${sel} -d Mt_lep_met_f   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_Mt-lep-met -p -r --qcd ${qcd} --title "" --title-x "M_{T}(l, MET)"  -c ${sel} -d Mt_lep_met   || true
	# Mt in lj out
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_lj_Mt-lep-met    -p -r --qcd ${qcd} --title "" --title-x "bck. cat. M_{T}(l, MET)"  -c ${sel}_lj    -d Mt_lep_met  || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_ljout_Mt-lep-met -p -r --qcd ${qcd} --title "" --title-x "sig. cat. M_{T}(l, MET)"  -c ${sel}_ljout -d Mt_lep_met  || true
	# rest
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_lep-pt     -p -r --qcd ${qcd} --title "" --title-x "lepton p_{T}"  -c ${sel} -d lep_pt   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_lep-eta    -p -r --qcd ${qcd} --title "" --title-x "lepton #eta"  -c ${sel} -d lep_eta   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_lep-pt_toppt  -p -r --qcd ${qcd} --title-x "lepton p_{T}"  -c ${sel} -d lep_pt -s TOPPTUp   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_met        -p -r --qcd ${qcd} --title "" --title-x "MET"  -c ${sel} -d met   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_met_toppt  -p -r --qcd ${qcd} --title-x "MET"  -c ${sel} -d met -s TOPPTUp   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_bjet-pt    -p -r --qcd ${qcd} --title "" --title-x "leading M b-jet p_{T}"  -c ${sel} -d bMjet_pt   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_bjet-eta   -p -r --qcd ${qcd} --title "" --title-x "leading M b-jet #eta"   -c ${sel} -d bMjet_eta   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_bLjet-pt   -p -r --qcd ${qcd} --title "" --title-x "leading L b-jet p_{T}"  -c ${sel} -d bLjet_pt   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_bLjet-eta  -p -r --qcd ${qcd} --title "" --title-x "leading L b-jet #eta"   -c ${sel} -d bLjet_eta   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_Rjet-pt    -p -r --qcd ${qcd} --title "" --title-x "leading non-b jet p_{T}"  -c ${sel} -d Rjet_pt   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_Rjet-eta   -p -r --qcd ${qcd} --title "" --title-x "leading non-b jet #eta"   -c ${sel} -d Rjet_eta   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_jet-pt     -p -r --qcd ${qcd} --title "" --title-x "leading non-b jet p_{T}"  -c ${sel} -d LRjet_pt   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_jet-eta    -p -r --qcd ${qcd} --title "" --title-x "leading non-b jet #eta"   -c ${sel} -d LRjet_eta   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_bMjet-pt_toppt -p -r --qcd ${qcd} --title-x "leading M b-jet p_{T}"  -c ${sel} -d bMjet_pt -s TOPPTUp   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_b_discr_med    -p -r --qcd ${qcd} --title "" --title-x "leading M b-jet CSVv2 b discr"  -c ${sel} -d b_discr_med   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_b_discr_loose  -p -r --qcd ${qcd} --title "" --title-x "leading L b-jet CSVv2 b discr"  -c ${sel} -d b_discr_loose   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_b_discr_rest   -p -r --qcd ${qcd} --title "" --title-x "leading non-b-jet CSVv2 b discr"  -c ${sel} -d b_discr_rest   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_nvtx       -p -r --qcd ${qcd} --title "" --title-x "N vtx"  -c ${sel} -d nvtx   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_nvtx_up    -p -r --qcd ${qcd} --title "" --title-x "PU Up N vtx"   -c ${sel} -d nvtx -s PUUp   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_nvtx_down  -p -r --qcd ${qcd} --title "" --title-x "PU Down N vtx" -c ${sel} -d nvtx -s PUDown   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_nvtx_init  -p -r --qcd ${qcd} --title "" --title-x "init N vtx"  -c ${sel} -d nvtx_init   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_tau-pt     -p -r --qcd ${qcd} --title "" --title-x "tau candidate p_{T}" -c ${sel} -d tau_pt   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_tau-pt     -p -r --qcd ${qcd} --title "" --title-x "tau candidate #eta"  -c ${sel} -d tau_eta   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_M-lep-tau  -p -r --qcd ${qcd} --title "" --title-x "lepton-tau candidate mass"  -c ${sel} -d M_lep_tau   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_Mt-tau-met -p -r --qcd ${qcd} --title "" --title-x "M_{T}(#tau, MET)"  -c ${sel} -d Mt_tau_met   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_njets_m    -p -r --qcd ${qcd} --title "" --title-x "N M b-tags"  -c ${sel} -d nMbjets   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_njets_l    -p -r --qcd ${qcd} --title "" --title-x "N L b-tags"  -c ${sel} -d nLbjets   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_njets_r    -p -r --qcd ${qcd} --title "" --title-x "N non-b-taged jets"  -c ${sel} -d nRjets   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_njets_a    -p -r --qcd ${qcd} --title "" --title-x "N all jets"  -c ${sel} -d njets   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_jet_cats   -p -r --qcd ${qcd} --title "" --title-x "jet categories"  -c ${sel} -d njets_cats   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_ntaus      -p -r --qcd ${qcd} --title "" --title-x "N tau candidates"  -c ${sel} -d ntaus   || true
	#python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} --lumi ${lumi}  --output-name ${tag}_dijet_m    -p -r --qcd ${qcd} --title "" --title-x "di-jet mass"  -c ${sel} -d dijet_mass   || true
	#python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} --lumi ${lumi}  --output-name ${tag}_trijet_m   -p -r --qcd ${qcd} --title "" --title-x "tri-jet mass" -c ${sel} -d trijet_mass   || true
	#python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} --lumi ${lumi}  --output-name ${tag}_dijet_trijet_d -p -r --qcd ${qcd} --title "" --title-x "di-tri-jet mass distance" -c ${sel} -d dijet_trijet_mass   || true
	# met control   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} -p -r --qcd ${qcd} --title-x "miniaod met"   -c ${sel} -d init_met   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} -p -r --qcd ${qcd} --title-x "corrected met" -c ${sel} -d corr_met   || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} -p -r --qcd ${qcd} -c ${sel} -d all_sum_control   || true

plots_sel: out_dir=./
plots_sel: sel=mu
plots_sel: tag=mu_sel
plots_sel: qcd=1.
plots_sel: lumi=35.8
plots_sel:
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_lj_Mt-lep-met    -p -r --qcd ${qcd} --title "" --title-x "bck. cat. M_{T}(l, MET)"  -c ${sel}_lj    -d Mt_lep_met  || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir}  --lumi ${lumi} --output-name ${tag}_ljout_Mt-lep-met -p -r --qcd ${qcd} --title "" --title-x "sig. cat. M_{T}(l, MET)"  -c ${sel}_ljout -d Mt_lep_met  || true

current_wjets:
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name wjets_el_Mt_lep_met_datadriven_qcd --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_el_wjet -d Mt_lep_met_f --qcd 1. --wjet 0.85 --title "" --title-x "M_{T}(e, MET)"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name wjets_el_Mt_lep_met                --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_el_wjet -d Mt_lep_met_f --qcd-factor 2. --wjet 0.85 --title "" --title-x "M_{T}(e, MET)"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name wjets_Mt_lep_met_datadriven_qcd    --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_mu_wjet -d Mt_lep_met_f --qcd 1. --wjet 0.8  --title "" --title-x "M_{T}(#mu, MET)"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name wjets_Mt_lep_met                   --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_mu_wjet -d Mt_lep_met_f --qcd-factor 2. --wjet 0.8  --title "" --title-x "M_{T}(#mu, MET)"
	# mu wjets lep, met
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name wjets_lep_pt_datadriven_qcd --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_mu_wjet -d lep_pt --qcd 1.        --wjet 0.8  --title "" --title-x "#mu p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name wjets_lep_pt                --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_mu_wjet -d lep_pt --qcd-factor 2. --wjet 0.8  --title "" --title-x "#mu p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name wjets_met_datadriven_qcd    --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_mu_wjet -d met --qcd 1.        --wjet 0.8  --title "" --title-x "MET"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name wjets_met                   --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_mu_wjet -d met --qcd-factor 2. --wjet 0.8  --title "" --title-x "MET"

current_presel:
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d Mt_lep_met_f --qcd 1. --wjet 0.85 --title "" --title-x "M_{T}(e, MET)"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d Mt_lep_met_f --qcd-factor 2. --wjet 0.85 --title "" --title-x "M_{T}(e, MET)"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d lep_pt --qcd 1.        --wjet 0.8  --title "" --title-x "e p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d lep_pt --qcd-factor 2. --wjet 0.8  --title "" --title-x "e p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d met --qcd 1.        --wjet 0.8  --title "" --title-x "MET"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d met --qcd-factor 2. --wjet 0.8  --title "" --title-x "MET"
	# mu
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_presel_Mt_lep_met        --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d Mt_lep_met_f --qcd 1. --wjet 0.8  --title "" --title-x "M_{T}(#mu, MET)"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_presel_Mt_lep_met_mc_qcd --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d Mt_lep_met_f --qcd-factor 2. --wjet 0.8  --title "" --title-x "M_{T}(#mu, MET)"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_presel_lep_pt_mc_qcd --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d lep_pt  --qcd-factor 2. --wjet 0.8  --title "" --title-x "#mu p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_presel_lep_pt        --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d lep_pt  --qcd 1.        --wjet 0.8  --title "" --title-x "#mu p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_presel_lep_eta       --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d lep_eta --qcd 1.        --wjet 0.8  --title "" --title-x "#mu #eta"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_presel_met           --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d met --qcd 1.        --wjet 0.8  --title "" --title-x "MET"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_presel_met_mc_qcd    --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d met --qcd-factor 2. --wjet 0.8  --title "" --title-x "MET"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_presel_bjet_pt       --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d bMjet_pt  --qcd 1.        --wjet 0.8  --title "" --title-x "leading b-jet p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_presel_bjet_eta      --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d bMjet_eta --qcd 1.        --wjet 0.8  --title "" --title-x "leading b-jet #eta"
	#python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d bMjet_pt  --qcd-factor 2. --wjet 0.8  --title "" --title-x "leading b-jet p_{T}"
	#python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d bMjet_eta --qcd-factor 2. --wjet 0.8  --title "" --title-x "leading b-jet #eta"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_presel_lep_pt     --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d lep_pt  --qcd 1.      --wjet 0.85  --title "" --title-x "e p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_presel_lep_eta    --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d lep_eta --qcd 1.      --wjet 0.85  --title "" --title-x "e #eta"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_presel_bjet_pt    --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d bMjet_pt  --qcd 1.    --wjet 0.85  --title "" --title-x "leading b-jet p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_presel_bjet_eta   --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d bMjet_eta --qcd 1.    --wjet 0.85  --title "" --title-x "leading b-jet #eta"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_presel_met        --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d met --qcd 1.          --wjet 0.85  --title "" --title-x "MET"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_presel_Mt_lep_met --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d Mt_lep_met_f --qcd 1. --wjet 0.85  --title "" --title-x "M_{T}(e, MET)"

current_presel_lj:
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --bin-norm --output-name mu_presel_dijet_mc_qcd       --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d dijet_mass  --qcd-factor 2. --wjet 0.8  --title "" --title-x "di-jet mass"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --bin-norm --output-name mu_presel_dijet              --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d dijet_mass  --qcd        1. --wjet 0.8  --title "" --title-x "di-jet mass"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --bin-norm --output-name mu_presel_trijet_mc_qcd      --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d trijet_mass --qcd-factor 2. --wjet 0.8  --title "" --title-x "tri-jet mass"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --bin-norm --output-name mu_presel_trijet             --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d trijet_mass --qcd        1. --wjet 0.8  --title "" --title-x "tri-jet mass"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --bin-norm --output-name mu_presel_dijettrijet_mc_qcd --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d dijet_trijet_mass --qcd-factor 2. --wjet 0.8  --title "" --title-x "di-jet tri-jet cut parameter"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --bin-norm --output-name mu_presel_dijettrijet        --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d dijet_trijet_mass --qcd        1. --wjet 0.8  --title "" --title-x "di-jet tri-jet cut parameter"
	#
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --bin-norm --output-name el_presel_dijet_mc_qcd       --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d dijet_mass  --qcd-factor 2. --wjet 0.8  --title "" --title-x "di-jet mass"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --bin-norm --output-name el_presel_dijet              --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d dijet_mass  --qcd        1. --wjet 0.8  --title "" --title-x "di-jet mass"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --bin-norm --output-name el_presel_trijet_mc_qcd      --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d trijet_mass --qcd-factor 2. --wjet 0.8  --title "" --title-x "tri-jet mass"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --bin-norm --output-name el_presel_trijet             --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d trijet_mass --qcd        1. --wjet 0.8  --title "" --title-x "tri-jet mass"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --bin-norm --output-name el_presel_dijettrijet_mc_qcd --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d dijet_trijet_mass --qcd-factor 2. --wjet 0.8  --title "" --title-x "di-jet tri-jet cut parameter"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --bin-norm --output-name el_presel_dijettrijet        --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d dijet_trijet_mass --qcd        1. --wjet 0.8  --title "" --title-x "di-jet tri-jet cut parameter"



current_sel:
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_sel_lep_pt        --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_sel -d lep_pt  --qcd 1.      --wjet 0.8  --title "" --title-x "#mu p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_sel_lep_eta       --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_sel -d lep_eta --qcd 1.      --wjet 0.8  --title "" --title-x "#mu #eta"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_sel_bjet_pt       --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_sel -d bMjet_pt  --qcd 1.    --wjet 0.8  --title "" --title-x "leading b-jet p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_sel_bjet_eta      --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_sel -d bMjet_eta --qcd 1.    --wjet 0.8  --title "" --title-x "leading b-jet #eta"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_sel_met           --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_sel -d met --qcd 1.          --wjet 0.8  --title "" --title-x "MET"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_sel_Mt_lep_met    --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_sel -d Mt_lep_met_f --qcd 1. --wjet 0.8  --title "" --title-x "M_{T}(#mu, MET)" -x ctr_old_mu_selVloose
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_sel_lj_Mt_lep_met    --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_sel_lj    -d Mt_lep_met_f --qcd 1. --wjet 0.8  --title "" --title-x "M_{T}(#mu, MET)" -x ctr_old_mu_selVloose
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root --output-name mu_sel_ljout_Mt_lep_met --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_sel_ljout -d Mt_lep_met_f --qcd 1. --wjet 0.8  --title "" --title-x "M_{T}(#mu, MET)" -x ctr_old_mu_selVloose
	#
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_sel_lep_pt     --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_sel -d lep_pt  --qcd 1.      --wjet 0.85  --title "" --title-x "e p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_sel_lep_eta    --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_sel -d lep_eta --qcd 1.      --wjet 0.85  --title "" --title-x "e #eta"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_sel_bjet_pt    --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_sel -d bMjet_pt  --qcd 1.    --wjet 0.85  --title "" --title-x "leading b-jet p_{T}"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_sel_bjet_eta   --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_sel -d bMjet_eta --qcd 1.    --wjet 0.85  --title "" --title-x "leading b-jet #eta"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_sel_met        --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_sel -d met --qcd 1.          --wjet 0.85  --title "" --title-x "MET"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_sel_Mt_lep_met --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_sel -d Mt_lep_met_f --qcd 1. --wjet 0.85  --title "" --title-x "M_{T}(e, MET)" -x ctr_old_el_selVloose
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_sel_lj_Mt_lep_met    --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_sel_lj    -d Mt_lep_met_f --qcd 1. --wjet 0.85  --title "" --title-x "M_{T}(e, MET)" -x ctr_old_el_selVloose
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_sel_ljout_Mt_lep_met --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_sel_ljout -d Mt_lep_met_f --qcd 1. --wjet 0.85  --title "" --title-x "M_{T}(e, MET)" -x ctr_old_el_selVloose
	#

current_pu:
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_sel_nvtx       --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_sel -d nvtx  --qcd 1. --wjet 0.85  --title "" --title-x "N vtx"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_sel_nvtx_up    --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_sel -d nvtx  --qcd 1. --wjet 0.85  --title "" --title-x "PU Up N vtx"   -s PUUp
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_sel_nvtx_down  --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_sel -d nvtx  --qcd 1. --wjet 0.85  --title "" --title-x "PU Down N vtx" -s PUDown
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name mu_sel_nvtx       --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_sel -d nvtx  --qcd 1. --wjet 0.80  --title "" --title-x "N vtx"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name mu_sel_nvtx_up    --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_sel -d nvtx  --qcd 1. --wjet 0.80  --title "" --title-x "PU Up N vtx"   -s PUUp
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name mu_sel_nvtx_down  --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_sel -d nvtx  --qcd 1. --wjet 0.80  --title "" --title-x "PU Down N vtx" -s PUDown
	#
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_presel_nvtx       --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d nvtx  --qcd 1. --wjet 0.85  --title "" --title-x "N vtx"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_presel_nvtx_up    --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d nvtx  --qcd 1. --wjet 0.85  --title "" --title-x "PU Up N vtx"   -s PUUp
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleElectron.root --output-name el_presel_nvtx_down  --exp-legend -o plots/current/ --lumi 31300 -p -r -c ctr_old_el_presel -d nvtx  --qcd 1. --wjet 0.85  --title "" --title-x "PU Down N vtx" -s PUDown
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name mu_presel_nvtx       --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d nvtx  --qcd 1. --wjet 0.80  --title "" --title-x "N vtx"
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name mu_presel_nvtx_up    --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d nvtx  --qcd 1. --wjet 0.80  --title "" --title-x "PU Up N vtx"   -s PUUp
	python stack_plot.py mc_v25p1.root ../merge-sets/v25/p1/SingleMuon.root     --output-name mu_presel_nvtx_down  --exp-legend -o plots/current/ --lumi 35800 -p -r -c ctr_old_mu_presel -d nvtx  --qcd 1. --wjet 0.80  --title "" --title-x "PU Down N vtx" -s PUDown

current_Mt_dep:
	python distr_ratios.py mc_v25p1.root --custom-filename mu_Mt_tau_pt_lt -o plots/current/ -c ctr_old_mu_sel:Mt_lep_met_f -p tt_mutau --no-data --ratio-channel-distrs ctr_old_mu_sel:Mt_lep_met_t2,ctr_old_mu_sel:Mt_lep_met_t3,ctr_old_mu_sel:Mt_lep_met_t4,ctr_old_mu_sel:Mt_lep_met_t5
	python distr_ratios.py mc_v25p1.root --custom-filename mu_Mt_tau_pt_lj -o plots/current/ -c ctr_old_mu_sel:Mt_lep_met_f -p tt_lj    --no-data --ratio-channel-distrs ctr_old_mu_sel:Mt_lep_met_t2,ctr_old_mu_sel:Mt_lep_met_t3,ctr_old_mu_sel:Mt_lep_met_t4,ctr_old_mu_sel:Mt_lep_met_t5
	python distr_ratios.py mc_v25p1.root --custom-filename el_Mt_tau_pt_lt -o plots/current/ -c ctr_old_el_sel:Mt_lep_met_f -p tt_eltau --no-data --ratio-channel-distrs ctr_old_el_sel:Mt_lep_met_t2,ctr_old_el_sel:Mt_lep_met_t3,ctr_old_el_sel:Mt_lep_met_t4,ctr_old_el_sel:Mt_lep_met_t5
	python distr_ratios.py mc_v25p1.root --custom-filename el_Mt_tau_pt_lj -o plots/current/ -c ctr_old_el_sel:Mt_lep_met_f -p tt_lj    --no-data --ratio-channel-distrs ctr_old_el_sel:Mt_lep_met_t2,ctr_old_el_sel:Mt_lep_met_t3,ctr_old_el_sel:Mt_lep_met_t4,ctr_old_el_sel:Mt_lep_met_t5
	#
	python distr_ratios.py mc_v25p1.root --custom-filename mu_presel_Mt_tau_pt_lt -o plots/current/ -c ctr_old_mu_presel:Mt_lep_met_f -p tt_mutau --no-data --ratio-channel-distrs ctr_old_mu_presel:Mt_lep_met_t2,ctr_old_mu_presel:Mt_lep_met_t3,ctr_old_mu_presel:Mt_lep_met_t4,ctr_old_mu_presel:Mt_lep_met_t5
	python distr_ratios.py mc_v25p1.root --custom-filename mu_presel_Mt_tau_pt_lj -o plots/current/ -c ctr_old_mu_presel:Mt_lep_met_f -p tt_lj    --no-data --ratio-channel-distrs ctr_old_mu_presel:Mt_lep_met_t2,ctr_old_mu_presel:Mt_lep_met_t3,ctr_old_mu_presel:Mt_lep_met_t4,ctr_old_mu_presel:Mt_lep_met_t5
	python distr_ratios.py mc_v25p1.root --custom-filename el_presel_Mt_tau_pt_lt -o plots/current/ -c ctr_old_el_presel:Mt_lep_met_f -p tt_eltau --no-data --ratio-channel-distrs ctr_old_el_presel:Mt_lep_met_t2,ctr_old_el_presel:Mt_lep_met_t3,ctr_old_el_presel:Mt_lep_met_t4,ctr_old_el_presel:Mt_lep_met_t5
	python distr_ratios.py mc_v25p1.root --custom-filename el_presel_Mt_tau_pt_lj -o plots/current/ -c ctr_old_el_presel:Mt_lep_met_f -p tt_lj    --no-data --ratio-channel-distrs ctr_old_el_presel:Mt_lep_met_t2,ctr_old_el_presel:Mt_lep_met_t3,ctr_old_el_presel:Mt_lep_met_t4,ctr_old_el_presel:Mt_lep_met_t5




plots_ctr: out_dir=./
plots_ctr: qcd=1.
plots_ctr: wjets=1.
plots_ctr: lumi=35.8
plots_ctr:
	# wjets
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} --output-name wjets_lep_pt     --lumi ${lumi} -p -r --qcd ${qcd} --wjets ${wjets} --title "" --title-x "muon p_{T}"    -c ctr_mu_wjet -d lep_pt       || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} --output-name wjets_Mt_lep_met --lumi ${lumi} -p -r --qcd ${qcd} --wjets ${wjets} --title "" --title-x "M_{T}(l, MET)" -c ctr_mu_wjet -d Mt_lep_met_f || true
	python stack_plot.py ${mc_in} ${data_in} -o ${out_dir} --output-name wjets_met        --lumi ${lumi} -p -r --qcd ${qcd} --wjets ${wjets} --title "" --title-x "MET"    -c ctr_mu_wjet -d met       || true


plots_sel_old: data_in=try2_split12/SingleMuon.root
plots_sel_old: mc_in=mc_mu_try2_split12.root
plots_sel_old: out_dir=./
plots_sel_old: sel=mu
plots_sel_old:
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



fit_datafile: data_in=../merge-sets/v23/p4/SingleMuon
fit_datafile: mc_in=v23p4_mc_fast_run1W0Jets
fit_datafile: lumi=35800
fit_datafile: wjet=0.85
fit_datafile: qcd=1.
fit_datafile: chan=ctr_old_mu_sel
fit_datafile: out_dir=histosel-roots/
fit_datafile: add_options=""
fit_datafile: distr=Mt_lep_met_f
fit_datafile:
	# histosel-roots//v23p4_mc_fast_run1W0Jets_HISTOSEL_Mt_lep_met_f_ctr_old_mu_sel_NOMINAL_data-qcd.root
	# the NOMINAL sets
	python stack_plot.py ${mc_in}.root ${data_in}.root --lumi ${lumi} --wjet ${wjet} --qcd ${qcd} -c ${chan} -d ${distr} ${add_options} -o ${out_dir} || true
	# the other systematics without data QCD
	#for SYS in TOPPTUp TOPPTDown PUUp PUDown LEPUp LEPDown JESUp JESDown JERUp JERDown TauESUp TauESDown bSFUp bSFDown;
	# object systs
	for SYS in JES JER TauES bSF; \
        do \
	python stack_plot.py ${mc_in}.root ${data_in}.root --lumi ${lumi} --wjet ${wjet} -s $${SYS}Up   -c ${chan} -d ${distr} ${add_options} -o ${out_dir} ; \
	python stack_plot.py ${mc_in}.root ${data_in}.root --lumi ${lumi} --wjet ${wjet} -s $${SYS}Down -c ${chan} -d ${distr} ${add_options} -o ${out_dir} ; \
        done
	# basic weight systs
	for SYS in LEP PU TOPPT AlphaS Frag SemilepBR Peterson Mr Mf Mfr; \
        do \
	python stack_plot.py ${mc_in}.root ${data_in}.root --lumi ${lumi} --wjet ${wjet} -s $${SYS}Up   -c ${chan} -d ${distr} ${add_options} -o ${out_dir} ; \
	python stack_plot.py ${mc_in}.root ${data_in}.root --lumi ${lumi} --wjet ${wjet} -s $${SYS}Down -c ${chan} -d ${distr} ${add_options} -o ${out_dir} ; \
        done
	# pdfs
	# PDFCT14n%dUp
	for NUM in `seq 56`; \
        do \
	python stack_plot.py ${mc_in}.root ${data_in}.root --lumi ${lumi} --wjet ${wjet} -s PDFCT14n$${NUM}Up   -c ${chan} -d ${distr} ${add_options} -o ${out_dir} ; \
	python stack_plot.py ${mc_in}.root ${data_in}.root --lumi ${lumi} --wjet ${wjet} -s PDFCT14n$${NUM}Down -c ${chan} -d ${distr} ${add_options} -o ${out_dir} ; \
        done
	# the hadd of produced files
	hadd ${out_dir}/${mc_in}_histosel_${chan}_${lumi}.root ${out_dir}/${mc_in}_HISTOSEL*_${chan}_*.root

clean_fit_datafile: mc_in=mc_mu_try2_split12
clean_fit_datafile:
	rm histosel-roots/${mc_in}_HISTOSEL*root


