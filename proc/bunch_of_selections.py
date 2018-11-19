
distrs_leptonic = [('std_mt_vars', 'Mt_lep_met_f'), ('event_dilep_mass', 'dilep_mass'), ('event_leptons[0].pt()', 'lep_pt')]
distrs_tauonic  = [('std_mt_vars', 'Mt_lep_met_f'), ('event_dilep_mass', 'dilep_mass'), ('event_leptons[0].pt()', 'lep_pt'),
        ('event_taus[0].pt()', 'tau_pt'), ('event_taus_sv_sign[0]', 'tau_sv_sign')]

# channels and conditions
tt_dileptons = (['tt_elmu'], '', distrs_leptonic)
tt_leptau    = (['el_sel', 'el_selVloose', 'el_sel_ss', 'el_selVloose_ss', 'mu_sel', 'mu_selVloose', 'mu_sel_ss', 'mu_selVloose_ss'], '--cond-com "event_taus_sv_sign[0] > 2.5"', distrs_tauonic)
dy_dileptons = (['dy_mumu', 'dy_elel'],   '--cond-com "std_mt_vars < 40."', distrs_leptonic)
dy_leptau    = (['dy_mutau', 'dy_eltau'], '--cond-com "event_taus_sv_sign[0] > 2.5 && std_mt_vars < 40."', distrs_tauonic)


# dtags and systematics
data = ['SingleMuon', 'SingleElectron'], ['nom']
tt = ['MC2016_Summer16_TTJets_powheg'],  ["nom,common", "obj", "tt_weights", "tt_hard", "tt_pdf1", "tt_pdf10", "tt_pdf20", "tt_pdf30", "tt_pdf40", "tt_pdf50,tt_alpha"]
other_mc = ['MC2016_Summer16_DYJetsToLL_10to50_amcatnlo',
'MC2016_Summer16_DYJetsToLL_50toInf_madgraph',
'MC2016_Summer16_SingleT_tW_5FS_powheg',
'MC2016_Summer16_SingleTbar_tW_5FS_powheg',
'MC2016_Summer16_W1Jets_madgraph',
'MC2016_Summer16_W2Jets_madgraph',
'MC2016_Summer16_W3Jets_madgraph',
'MC2016_Summer16_W4Jets_madgraph',
'MC2016_Summer16_WJets_madgraph',
'MC2016_Summer16_WWTo2L2Nu_powheg',
'MC2016_Summer16_WWToLNuQQ_powheg',
'MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin',
'MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin',
'MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin',
'MC2016_Summer16_WZTo3LNu_powheg',
'MC2016_Summer16_ZZTo2L2Nu_powheg',
'MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin',
'MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo',
'MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg',
'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg'], ['nom', 'common', 'obj']


intro = """
#!/bin/sh
export X509_USER_PROXY=/home/t3cms/olek/x509_proxy
export SCRAM_ARCH=slc6_amd64_gcc530
export BUILD_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
export CMS_PATH=$VO_CMS_SW_DIR
cd /lstore/cms/olek/CMSSW_8_0_26_patch1/src/
cmsenv
cd UserCode/NtuplerAnalyzer/proc/

{draw_commands}
"""

draw_command_template = """python -W ignore sumup_ttree_draw.py --cut-w0jets "{draw}"           --ttree ttree_out --histo-range std  --output quick-test/v27-dilep2/{dtag}_{chan}_{sys}_{distr_name}.root  --std-histos --histo-name {chan}/std_procs/{sys}/{distr_name} {cond} --per-weight --try-xsec --save-weight lstore_outdirs/merge-sets/v27/dilep1/*${dtag}*.root"""

#batch_jobs/j_${dtag}_${chans}_${systs}

for (dtags, systs) in [data, tt, other_mc]:
    for dtag in dtags:
        for sys in systs:
            for chans, cond, distrs  in [tt_dileptons, tt_leptau, dy_dileptons, dy_leptau]:
                for chan in chans:
                    with open('batch_jobs/j_{dtag}_{chan}_{sys}'.format(dtag=dtag, chan=chan, sys=sys), 'w') as f:
                        draw_coms = '\n'.join([draw_command_template.format(draw=draw, distr_name=name, dtag=dtag, chan=chan, sys=sys, cond=cond) for draw, name in distrs])
                        f.write(intro.format(draw_commands=draw_coms) + '\n')


