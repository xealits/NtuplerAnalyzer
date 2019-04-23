import logging

# new selection stages
all_std_channels = {
'mu_selSV':          ('({selection_stage}==  9 || {selection_stage}==  7) && event_taus_sv_sign[0] > 2.5', 'selection_stage'),
'mu_selSV_ss':       ('({selection_stage}==  8 || {selection_stage}==  6) && event_taus_sv_sign[0] > 2.5', 'selection_stage'),
'el_selSV':          ('({selection_stage}== 19 || {selection_stage}== 17) && event_taus_sv_sign[0] > 2.5', 'selection_stage'),
'el_selSV_ss':       ('({selection_stage}== 18 || {selection_stage}== 16) && event_taus_sv_sign[0] > 2.5', 'selection_stage'),
'mu_selSVVloose':    ('({selection_stage}==  9 || {selection_stage}==  7 || {selection_stage}==  5) && event_taus_sv_sign[0] > 2.5', 'selection_stage'),
'mu_selSVVloose_ss': ('({selection_stage}==  8 || {selection_stage}==  6 || {selection_stage}==  4) && event_taus_sv_sign[0] > 2.5', 'selection_stage'),
'el_selSVVloose':    ('({selection_stage}== 19 || {selection_stage}== 17 || {selection_stage}== 15) && event_taus_sv_sign[0] > 2.5', 'selection_stage'),
'el_selSVVloose_ss': ('({selection_stage}== 18 || {selection_stage}== 16 || {selection_stage}== 14) && event_taus_sv_sign[0] > 2.5', 'selection_stage'),

'mu_selTight':     ('({selection_stage}==  9)', 'selection_stage'),
'mu_selTight_ss':  ('({selection_stage}==  8)', 'selection_stage'),
'el_selTight':     ('({selection_stage}== 19)', 'selection_stage'),
'el_selTight_ss':  ('({selection_stage}== 18)', 'selection_stage'),
'mu_sel':          ('({selection_stage}==  9 || {selection_stage}==  7)', 'selection_stage'),
'mu_sel_ss':       ('({selection_stage}==  8 || {selection_stage}==  6)', 'selection_stage'),
'el_sel':          ('({selection_stage}== 19 || {selection_stage}== 17)', 'selection_stage'),
'el_sel_ss':       ('({selection_stage}== 18 || {selection_stage}== 16)', 'selection_stage'),
'mu_selVloose':    ('({selection_stage}==  9 || {selection_stage}==  7 || {selection_stage}==  5)', 'selection_stage'),
'mu_selVloose_ss': ('({selection_stage}==  8 || {selection_stage}==  6 || {selection_stage}==  4)', 'selection_stage'),
'el_selVloose':    ('({selection_stage}== 19 || {selection_stage}== 17 || {selection_stage}== 15)', 'selection_stage'),
'el_selVloose_ss': ('({selection_stage}== 18 || {selection_stage}== 16 || {selection_stage}== 14)', 'selection_stage'),

'mu_selTight_ljout':     ('({selection_stage}==  9) && event_jets_lj_var >  60.', 'selection_stage'),
'mu_selTight_ljout_ss':  ('({selection_stage}==  8) && event_jets_lj_var >  60.', 'selection_stage'),
'el_selTight_ljout':     ('({selection_stage}== 19) && event_jets_lj_var >  60.', 'selection_stage'),
'el_selTight_ljout_ss':  ('({selection_stage}== 18) && event_jets_lj_var >  60.', 'selection_stage'),
'mu_sel_ljout':          ('({selection_stage}==  9 || {selection_stage}==  7) && event_jets_lj_var >  60.', 'selection_stage'),
'mu_sel_ljout_ss':       ('({selection_stage}==  8 || {selection_stage}==  6) && event_jets_lj_var >  60.', 'selection_stage'),
'el_sel_ljout':          ('({selection_stage}== 19 || {selection_stage}== 17) && event_jets_lj_var >  60.', 'selection_stage'),
'el_sel_ljout_ss':       ('({selection_stage}== 18 || {selection_stage}== 16) && event_jets_lj_var >  60.', 'selection_stage'),
'mu_selVloose_ljout':    ('({selection_stage}==  9 || {selection_stage}==  7 || {selection_stage}==  5) && event_jets_lj_var >  60.', 'selection_stage'),
'mu_selVloose_ljout_ss': ('({selection_stage}==  8 || {selection_stage}==  6 || {selection_stage}==  4) && event_jets_lj_var >  60.', 'selection_stage'),
'el_selVloose_ljout':    ('({selection_stage}== 19 || {selection_stage}== 17 || {selection_stage}== 15) && event_jets_lj_var >  60.', 'selection_stage'),
'el_selVloose_ljout_ss': ('({selection_stage}== 18 || {selection_stage}== 16 || {selection_stage}== 14) && event_jets_lj_var >  60.', 'selection_stage'),

'mu_selTight_lj':     ('({selection_stage}==  9) && event_jets_lj_var <= 60.', 'selection_stage'),
'mu_selTight_lj_ss':  ('({selection_stage}==  8) && event_jets_lj_var <= 60.', 'selection_stage'),
'el_selTight_lj':     ('({selection_stage}== 19) && event_jets_lj_var <= 60.', 'selection_stage'),
'el_selTight_lj_ss':  ('({selection_stage}== 18) && event_jets_lj_var <= 60.', 'selection_stage'),
'mu_sel_lj':          ('({selection_stage}==  9 || {selection_stage}==  7) && event_jets_lj_var <= 60.', 'selection_stage'),
'mu_sel_lj_ss':       ('({selection_stage}==  8 || {selection_stage}==  6) && event_jets_lj_var <= 60.', 'selection_stage'),
'el_sel_lj':          ('({selection_stage}== 19 || {selection_stage}== 17) && event_jets_lj_var <= 60.', 'selection_stage'),
'el_sel_lj_ss':       ('({selection_stage}== 18 || {selection_stage}== 16) && event_jets_lj_var <= 60.', 'selection_stage'),
'mu_selVloose_lj':    ('({selection_stage}==  9 || {selection_stage}==  7 || {selection_stage}==  5) && event_jets_lj_var <= 60.', 'selection_stage'),
'mu_selVloose_lj_ss': ('({selection_stage}==  8 || {selection_stage}==  6 || {selection_stage}==  4) && event_jets_lj_var <= 60.', 'selection_stage'),
'el_selVloose_lj':    ('({selection_stage}== 19 || {selection_stage}== 17 || {selection_stage}== 15) && event_jets_lj_var <= 60.', 'selection_stage'),
'el_selVloose_lj_ss': ('({selection_stage}== 18 || {selection_stage}== 16 || {selection_stage}== 14) && event_jets_lj_var <= 60.', 'selection_stage'),

# additional channels
'dy_mutau': ('({selection_stage}== 102 || {selection_stage}== 103)', 'selection_stage_dy'),
'dy_eltau': ('({selection_stage}== 112 || {selection_stage}== 113)', 'selection_stage_dy'),
'dy_mumu':  ('({selection_stage}== 102 || {selection_stage}== 103 || {selection_stage}== 105)', 'selection_stage_dy_mumu'),
'dy_elel':  ('({selection_stage}== 112 || {selection_stage}== 113 || {selection_stage}== 115)', 'selection_stage_dy_mumu'),

'tt_elmu':  ('({selection_stage}== 205)', 'selection_stage_em'),
'tt_alliso_presel_el'    :  ('({selection_stage}== 513)', 'selection_stage_tt_alliso'),
'tt_alliso_presel_el_ss' :  ('({selection_stage}== 512)', 'selection_stage_tt_alliso'),
'tt_alliso_presel_mu'    :  ('({selection_stage}== 503)', 'selection_stage_tt_alliso'),
'tt_alliso_presel_mu_ss' :  ('({selection_stage}== 502)', 'selection_stage_tt_alliso'),
}


systs_weights_nominal = {
'NOMINAL': 'event_weight',
}

systs_weights_common = {
'bSFUp'  : "event_weight*event_weight_bSFUp/event_weight_bSF"  ,
'bSFDown': "event_weight*event_weight_bSFDown/event_weight_bSF",

'PUUp'        : "event_weight*event_weight_PUUp"    ,
'PUDown'      : "event_weight*event_weight_PUDown"  ,
'LEPelIDUp'   : "event_weight*event_weight_LEPelIDUp"   ,
'LEPelIDDown' : "event_weight*event_weight_LEPelIDDown" ,
'LEPelTRGUp'  : "event_weight*event_weight_LEPelTRGUp"  ,
'LEPelTRGDown': "event_weight*event_weight_LEPelTRGDown",
'LEPmuIDUp'   : "event_weight*event_weight_LEPmuIDUp"   ,
'LEPmuIDDown' : "event_weight*event_weight_LEPmuIDDown" ,
'LEPmuTRGUp'  : "event_weight*event_weight_LEPmuTRGUp"  ,
'LEPmuTRGDown': "event_weight*event_weight_LEPmuTRGDown",
}

systs_weights_tt = {
'TOPPTDown'     : 'event_weight'                           ,
'TOPPTUp'       : 'event_weight*event_weight_toppt'        ,
'FragUp'        : 'event_weight*event_weight_FragUp'       ,
'FragDown'      : 'event_weight*event_weight_FragDown'     ,
'SemilepBRUp'   : 'event_weight*event_weight_SemilepBRUp'  ,
'SemilepBRDown' : 'event_weight*event_weight_SemilepBRDown',
'PetersonUp'    : 'event_weight*event_weight_PetersonUp'   ,
'PetersonDown'  : 'event_weight'                           ,
}

systs_weights_tt_hard = {
"MrUp"    : "event_weight*event_weight_me_f_rUp",
"MrDown"  : "event_weight*event_weight_me_f_rDn",
"MfUp"    : "event_weight*event_weight_me_fUp_r",
"MfDown"  : "event_weight*event_weight_me_fDn_r",
"MfrUp"   : "event_weight*event_weight_me_frUp" ,
"MfrDown" : "event_weight*event_weight_me_frDn" ,
}

systs_weights_tt_alpha = {
'AlphaSUp'  : "event_weight*event_weight_AlphaS_up",
'AlphaSDown': "event_weight*event_weight_AlphaS_dn",
}

systs_weights_tt_pdf1 = {
'PDFCT14n1Up'     : "event_weight*event_weight_pdf[0]" ,
'PDFCT14n2Up'     : "event_weight*event_weight_pdf[1]" ,
'PDFCT14n3Up'     : "event_weight*event_weight_pdf[2]" ,
'PDFCT14n4Up'     : "event_weight*event_weight_pdf[3]" ,
'PDFCT14n5Up'     : "event_weight*event_weight_pdf[4]" ,
'PDFCT14n6Up'     : "event_weight*event_weight_pdf[5]" ,
'PDFCT14n7Up'     : "event_weight*event_weight_pdf[6]" ,
'PDFCT14n8Up'     : "event_weight*event_weight_pdf[7]" ,
'PDFCT14n9Up'     : "event_weight*event_weight_pdf[8]" ,
'PDFCT14n10Up'    : "event_weight*event_weight_pdf[9]" ,
}

systs_weights_tt_pdf10 = {
'PDFCT14n11Up'    : "event_weight*event_weight_pdf[10]",
'PDFCT14n12Up'    : "event_weight*event_weight_pdf[11]",
'PDFCT14n13Up'    : "event_weight*event_weight_pdf[12]",
'PDFCT14n14Up'    : "event_weight*event_weight_pdf[13]",
'PDFCT14n15Up'    : "event_weight*event_weight_pdf[14]",
'PDFCT14n16Up'    : "event_weight*event_weight_pdf[15]",
'PDFCT14n17Up'    : "event_weight*event_weight_pdf[16]",
'PDFCT14n18Up'    : "event_weight*event_weight_pdf[17]",
'PDFCT14n19Up'    : "event_weight*event_weight_pdf[18]",
'PDFCT14n20Up'    : "event_weight*event_weight_pdf[19]",
}

systs_weights_tt_pdf20 = {
'PDFCT14n21Up'    : "event_weight*event_weight_pdf[20]",
'PDFCT14n22Up'    : "event_weight*event_weight_pdf[21]",
'PDFCT14n23Up'    : "event_weight*event_weight_pdf[22]",
'PDFCT14n24Up'    : "event_weight*event_weight_pdf[23]",
'PDFCT14n25Up'    : "event_weight*event_weight_pdf[24]",
'PDFCT14n26Up'    : "event_weight*event_weight_pdf[25]",
'PDFCT14n27Up'    : "event_weight*event_weight_pdf[26]",
'PDFCT14n28Up'    : "event_weight*event_weight_pdf[27]",
'PDFCT14n29Up'    : "event_weight*event_weight_pdf[28]",
'PDFCT14n30Up'    : "event_weight*event_weight_pdf[29]",
}

systs_weights_tt_pdf30 = {
'PDFCT14n31Up'    : "event_weight*event_weight_pdf[30]",
'PDFCT14n32Up'    : "event_weight*event_weight_pdf[31]",
'PDFCT14n33Up'    : "event_weight*event_weight_pdf[32]",
'PDFCT14n34Up'    : "event_weight*event_weight_pdf[33]",
'PDFCT14n35Up'    : "event_weight*event_weight_pdf[34]",
'PDFCT14n36Up'    : "event_weight*event_weight_pdf[35]",
'PDFCT14n37Up'    : "event_weight*event_weight_pdf[36]",
'PDFCT14n38Up'    : "event_weight*event_weight_pdf[37]",
'PDFCT14n39Up'    : "event_weight*event_weight_pdf[38]",
'PDFCT14n40Up'    : "event_weight*event_weight_pdf[39]",
}

systs_weights_tt_pdf40 = {
'PDFCT14n41Up'    : "event_weight*event_weight_pdf[40]",
'PDFCT14n42Up'    : "event_weight*event_weight_pdf[41]",
'PDFCT14n43Up'    : "event_weight*event_weight_pdf[42]",
'PDFCT14n44Up'    : "event_weight*event_weight_pdf[43]",
'PDFCT14n45Up'    : "event_weight*event_weight_pdf[44]",
'PDFCT14n46Up'    : "event_weight*event_weight_pdf[45]",
'PDFCT14n47Up'    : "event_weight*event_weight_pdf[46]",
'PDFCT14n48Up'    : "event_weight*event_weight_pdf[47]",
'PDFCT14n49Up'    : "event_weight*event_weight_pdf[48]",
'PDFCT14n50Up'    : "event_weight*event_weight_pdf[49]",
}

systs_weights_tt_pdf50 = {
'PDFCT14n51Up'    : "event_weight*event_weight_pdf[50]",
'PDFCT14n52Up'    : "event_weight*event_weight_pdf[51]",
'PDFCT14n53Up'    : "event_weight*event_weight_pdf[52]",
'PDFCT14n54Up'    : "event_weight*event_weight_pdf[53]",
'PDFCT14n55Up'    : "event_weight*event_weight_pdf[54]",
'PDFCT14n56Up'    : "event_weight*event_weight_pdf[55]",
}


systs_weights_tt_alpha_pdf = {
'AlphaSUp'  : "event_weight*event_weight_AlphaS_up",
'AlphaSDown': "event_weight*event_weight_AlphaS_dn",

'PDFCT14n1Up'     : "event_weight*event_weight_pdf[0]" ,
'PDFCT14n2Up'     : "event_weight*event_weight_pdf[1]" ,
'PDFCT14n3Up'     : "event_weight*event_weight_pdf[2]" ,
'PDFCT14n4Up'     : "event_weight*event_weight_pdf[3]" ,
'PDFCT14n5Up'     : "event_weight*event_weight_pdf[4]" ,
'PDFCT14n6Up'     : "event_weight*event_weight_pdf[5]" ,
'PDFCT14n7Up'     : "event_weight*event_weight_pdf[6]" ,
'PDFCT14n8Up'     : "event_weight*event_weight_pdf[7]" ,
'PDFCT14n9Up'     : "event_weight*event_weight_pdf[8]" ,
'PDFCT14n10Up'    : "event_weight*event_weight_pdf[9]" ,
'PDFCT14n11Up'    : "event_weight*event_weight_pdf[10]",
'PDFCT14n12Up'    : "event_weight*event_weight_pdf[11]",
'PDFCT14n13Up'    : "event_weight*event_weight_pdf[12]",
'PDFCT14n14Up'    : "event_weight*event_weight_pdf[13]",
'PDFCT14n15Up'    : "event_weight*event_weight_pdf[14]",
'PDFCT14n16Up'    : "event_weight*event_weight_pdf[15]",
'PDFCT14n17Up'    : "event_weight*event_weight_pdf[16]",
'PDFCT14n18Up'    : "event_weight*event_weight_pdf[17]",
'PDFCT14n19Up'    : "event_weight*event_weight_pdf[18]",
'PDFCT14n20Up'    : "event_weight*event_weight_pdf[19]",
'PDFCT14n21Up'    : "event_weight*event_weight_pdf[20]",
'PDFCT14n22Up'    : "event_weight*event_weight_pdf[21]",
'PDFCT14n23Up'    : "event_weight*event_weight_pdf[22]",
'PDFCT14n24Up'    : "event_weight*event_weight_pdf[23]",
'PDFCT14n25Up'    : "event_weight*event_weight_pdf[24]",
'PDFCT14n26Up'    : "event_weight*event_weight_pdf[25]",
'PDFCT14n27Up'    : "event_weight*event_weight_pdf[26]",
'PDFCT14n28Up'    : "event_weight*event_weight_pdf[27]",
'PDFCT14n29Up'    : "event_weight*event_weight_pdf[28]",
'PDFCT14n30Up'    : "event_weight*event_weight_pdf[29]",
'PDFCT14n31Up'    : "event_weight*event_weight_pdf[30]",
'PDFCT14n32Up'    : "event_weight*event_weight_pdf[31]",
'PDFCT14n33Up'    : "event_weight*event_weight_pdf[32]",
'PDFCT14n34Up'    : "event_weight*event_weight_pdf[33]",
'PDFCT14n35Up'    : "event_weight*event_weight_pdf[34]",
'PDFCT14n36Up'    : "event_weight*event_weight_pdf[35]",
'PDFCT14n37Up'    : "event_weight*event_weight_pdf[36]",
'PDFCT14n38Up'    : "event_weight*event_weight_pdf[37]",
'PDFCT14n39Up'    : "event_weight*event_weight_pdf[38]",
'PDFCT14n40Up'    : "event_weight*event_weight_pdf[39]",
'PDFCT14n41Up'    : "event_weight*event_weight_pdf[40]",
'PDFCT14n42Up'    : "event_weight*event_weight_pdf[41]",
'PDFCT14n43Up'    : "event_weight*event_weight_pdf[42]",
'PDFCT14n44Up'    : "event_weight*event_weight_pdf[43]",
'PDFCT14n45Up'    : "event_weight*event_weight_pdf[44]",
'PDFCT14n46Up'    : "event_weight*event_weight_pdf[45]",
'PDFCT14n47Up'    : "event_weight*event_weight_pdf[46]",
'PDFCT14n48Up'    : "event_weight*event_weight_pdf[47]",
'PDFCT14n49Up'    : "event_weight*event_weight_pdf[48]",
'PDFCT14n50Up'    : "event_weight*event_weight_pdf[49]",
'PDFCT14n51Up'    : "event_weight*event_weight_pdf[50]",
'PDFCT14n52Up'    : "event_weight*event_weight_pdf[51]",
'PDFCT14n53Up'    : "event_weight*event_weight_pdf[52]",
'PDFCT14n54Up'    : "event_weight*event_weight_pdf[53]",
'PDFCT14n55Up'    : "event_weight*event_weight_pdf[54]",
'PDFCT14n56Up'    : "event_weight*event_weight_pdf[55]",
}

systs_weights_tt_updowns = {
'TuneCUETP8M2T4Down' : 'event_weight',
'TuneCUETP8M2T4Up'   : 'event_weight',
'FSRDown'            : 'event_weight',
'FSRUp'              : 'event_weight',
'HDAMPDown'          : 'event_weight',
'HDAMPUp'            : 'event_weight',
'ISRDown'            : 'event_weight',
'ISRUp'              : 'event_weight',
}

systs_weights_objects = {
'JERDown':  'event_weight',
'JERUp'  :  'event_weight',
'JESDown':  'event_weight',
'JESUp'  :  'event_weight',
'TESDown':  'event_weight',
'TESUp'  :  'event_weight',
}

updown_dtags_to_sys = {
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4down' : 'TuneCUETP8M2T4Down',
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4up'   : 'TuneCUETP8M2T4Up',
'MC2016_Summer16_TTJets_powheg_fsrdown'        : 'FSRDown',
'MC2016_Summer16_TTJets_powheg_fsrup'          : 'FSRUp',
'MC2016_Summer16_TTJets_powheg_hdampDOWN'      : 'HDAMPDown',
'MC2016_Summer16_TTJets_powheg_hdampUP'        : 'HDAMPUp',
'MC2016_Summer16_TTJets_powheg_isrdown'        : 'ISRDown',
'MC2016_Summer16_TTJets_powheg_isrup'          : 'ISRUp',
}

named_systs_weights_all = {'nom': systs_weights_nominal,
'common': systs_weights_common,
'tt_weights': systs_weights_tt,
'tt_hard':  systs_weights_tt_hard,
'tt_pdf':   systs_weights_tt_alpha_pdf,
'tt_alpha': systs_weights_tt_alpha,
'tt_pdf1':  systs_weights_tt_pdf1,
'tt_pdf10': systs_weights_tt_pdf10,
'tt_pdf20': systs_weights_tt_pdf20,
'tt_pdf30': systs_weights_tt_pdf30,
'tt_pdf40': systs_weights_tt_pdf40,
'tt_pdf50': systs_weights_tt_pdf50,
}

def extend_full_sys_list(sys_list_with_nicknames):
    # the list may contain either known nicknames of groups of systematics
    # or concrete names of systematics
    full_list = []
    for sname in sys_list_with_nicknames:
        if sname in named_systs_weights_all:
            # extend-expand this list of the nickname for the output
            full_list.extend(named_systs_weights_all[sname])
        else:
            # otherwise it's a concrete systematic name
            full_list.append(sname)
    return full_list

systs_weights_all = {}
for s_d in named_systs_weights_all.values():
    systs_weights_all.update(s_d)

systs_weights_all.update(systs_weights_tt_updowns)

# systs2 with event functions

systs2_weights_initial = {
'INITIAL': lambda ev: 1.,
}

systs2_weights_nominal = {
'NOMINAL': lambda ev: ev.event_weight,
}

systs2_weights_common = {
'PUUp'   : lambda ev: ev.event_weight*ev.event_weight_PUUp    / ev.event_weight_PU  if ev.event_weight_PU  > 0. else 0.,
'PUDown' : lambda ev: ev.event_weight*ev.event_weight_PUDown  / ev.event_weight_PU  if ev.event_weight_PU  > 0. else 0.,
'bSFUp'  : lambda ev: ev.event_weight*ev.event_weight_bSFUp   / ev.event_weight_bSF if ev.event_weight_bSF > 0. else 0.,
'bSFDown': lambda ev: ev.event_weight*ev.event_weight_bSFDown / ev.event_weight_bSF if ev.event_weight_bSF > 0. else 0.,

'LEPelIDUp'   : lambda ev: ev.event_weight*ev.event_weight_LEPelIDUp   ,
'LEPelIDDown' : lambda ev: ev.event_weight*ev.event_weight_LEPelIDDown ,
'LEPelTRGUp'  : lambda ev: ev.event_weight*ev.event_weight_LEPelTRGUp  ,
'LEPelTRGDown': lambda ev: ev.event_weight*ev.event_weight_LEPelTRGDown,
'LEPmuIDUp'   : lambda ev: ev.event_weight*ev.event_weight_LEPmuIDUp   ,
'LEPmuIDDown' : lambda ev: ev.event_weight*ev.event_weight_LEPmuIDDown ,
'LEPmuTRGUp'  : lambda ev: ev.event_weight*ev.event_weight_LEPmuTRGUp  ,
'LEPmuTRGDown': lambda ev: ev.event_weight*ev.event_weight_LEPmuTRGDown,
}

systs2_weights_tt = {
'TOPPTDown'     : lambda ev: ev.event_weight                           ,
'TOPPTUp'       : lambda ev: ev.event_weight*ev.event_weight_toppt        ,
'FragUp'        : lambda ev: ev.event_weight*ev.event_weight_FragUp       ,
'FragDown'      : lambda ev: ev.event_weight*ev.event_weight_FragDown     ,
'SemilepBRUp'   : lambda ev: ev.event_weight*ev.event_weight_SemilepBRUp  ,
'SemilepBRDown' : lambda ev: ev.event_weight*ev.event_weight_SemilepBRDown,
'PetersonUp'    : lambda ev: ev.event_weight*ev.event_weight_PetersonUp   ,
'PetersonDown'  : lambda ev: ev.event_weight                           ,
}

systs2_weights_tt_hard = {
'MrUp'    : lambda ev: ev.event_weight*ev.event_weight_me_f_rUp,
'MrDown'  : lambda ev: ev.event_weight*ev.event_weight_me_f_rDn,
'MfUp'    : lambda ev: ev.event_weight*ev.event_weight_me_fUp_r,
'MfDown'  : lambda ev: ev.event_weight*ev.event_weight_me_fDn_r,
'MfrUp'   : lambda ev: ev.event_weight*ev.event_weight_me_frUp ,
'MfrDown' : lambda ev: ev.event_weight*ev.event_weight_me_frDn ,
}

systs2_weights_tt_alpha = {
'AlphaSUp'  : lambda ev: ev.event_weight*ev.ev.event_weight_AlphaS_up,
'AlphaSDown': lambda ev: ev.event_weight*ev.ev.event_weight_AlphaS_dn,
}

systs2_weights_tt_pdf1 = {
'PDFCT14n1Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[0] ,
'PDFCT14n2Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[1] ,
'PDFCT14n3Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[2] ,
'PDFCT14n4Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[3] ,
'PDFCT14n5Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[4] ,
'PDFCT14n6Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[5] ,
'PDFCT14n7Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[6] ,
'PDFCT14n8Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[7] ,
'PDFCT14n9Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[8] ,
'PDFCT14n10Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[9] ,
}

systs2_weights_tt_pdf10 = {
'PDFCT14n11Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[10],
'PDFCT14n12Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[11],
'PDFCT14n13Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[12],
'PDFCT14n14Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[13],
'PDFCT14n15Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[14],
'PDFCT14n16Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[15],
'PDFCT14n17Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[16],
'PDFCT14n18Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[17],
'PDFCT14n19Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[18],
'PDFCT14n20Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[19],
}

systs2_weights_tt_pdf20 = {
'PDFCT14n21Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[20],
'PDFCT14n22Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[21],
'PDFCT14n23Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[22],
'PDFCT14n24Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[23],
'PDFCT14n25Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[24],
'PDFCT14n26Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[25],
'PDFCT14n27Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[26],
'PDFCT14n28Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[27],
'PDFCT14n29Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[28],
'PDFCT14n30Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[29],
}

systs2_weights_tt_pdf30 = {
'PDFCT14n31Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[30],
'PDFCT14n32Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[31],
'PDFCT14n33Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[32],
'PDFCT14n34Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[33],
'PDFCT14n35Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[34],
'PDFCT14n36Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[35],
'PDFCT14n37Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[36],
'PDFCT14n38Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[37],
'PDFCT14n39Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[38],
'PDFCT14n40Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[39],
}

systs2_weights_tt_pdf40 = {
'PDFCT14n41Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[40],
'PDFCT14n42Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[41],
'PDFCT14n43Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[42],
'PDFCT14n44Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[43],
'PDFCT14n45Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[44],
'PDFCT14n46Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[45],
'PDFCT14n47Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[46],
'PDFCT14n48Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[47],
'PDFCT14n49Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[48],
'PDFCT14n50Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[49],
}

systs2_weights_tt_pdf50 = {
'PDFCT14n51Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[50],
'PDFCT14n52Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[51],
'PDFCT14n53Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[52],
'PDFCT14n54Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[53],
'PDFCT14n55Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[54],
'PDFCT14n56Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[55],
}


systs2_weights_tt_alpha_pdf = {
'AlphaSUp'  : lambda ev: ev.event_weight*ev.event_weight_AlphaS_up,
'AlphaSDown': lambda ev: ev.event_weight*ev.event_weight_AlphaS_dn,

'PDFCT14n1Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[0] ,
'PDFCT14n2Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[1] ,
'PDFCT14n3Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[2] ,
'PDFCT14n4Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[3] ,
'PDFCT14n5Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[4] ,
'PDFCT14n6Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[5] ,
'PDFCT14n7Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[6] ,
'PDFCT14n8Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[7] ,
'PDFCT14n9Up'     : lambda ev: ev.event_weight*ev.event_weight_pdf[8] ,
'PDFCT14n10Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[9] ,
'PDFCT14n11Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[10],
'PDFCT14n12Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[11],
'PDFCT14n13Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[12],
'PDFCT14n14Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[13],
'PDFCT14n15Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[14],
'PDFCT14n16Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[15],
'PDFCT14n17Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[16],
'PDFCT14n18Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[17],
'PDFCT14n19Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[18],
'PDFCT14n20Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[19],
'PDFCT14n21Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[20],
'PDFCT14n22Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[21],
'PDFCT14n23Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[22],
'PDFCT14n24Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[23],
'PDFCT14n25Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[24],
'PDFCT14n26Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[25],
'PDFCT14n27Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[26],
'PDFCT14n28Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[27],
'PDFCT14n29Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[28],
'PDFCT14n30Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[29],
'PDFCT14n31Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[30],
'PDFCT14n32Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[31],
'PDFCT14n33Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[32],
'PDFCT14n34Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[33],
'PDFCT14n35Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[34],
'PDFCT14n36Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[35],
'PDFCT14n37Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[36],
'PDFCT14n38Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[37],
'PDFCT14n39Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[38],
'PDFCT14n40Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[39],
'PDFCT14n41Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[40],
'PDFCT14n42Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[41],
'PDFCT14n43Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[42],
'PDFCT14n44Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[43],
'PDFCT14n45Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[44],
'PDFCT14n46Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[45],
'PDFCT14n47Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[46],
'PDFCT14n48Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[47],
'PDFCT14n49Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[48],
'PDFCT14n50Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[49],
'PDFCT14n51Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[50],
'PDFCT14n52Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[51],
'PDFCT14n53Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[52],
'PDFCT14n54Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[53],
'PDFCT14n55Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[54],
'PDFCT14n56Up'    : lambda ev: ev.event_weight*ev.event_weight_pdf[55],
}

systs2_weights_tt_updowns = {
'TuneCUETP8M2T4Down' : lambda ev: ev.event_weight,
'TuneCUETP8M2T4Up'   : lambda ev: ev.event_weight,
'FSRDown'            : lambda ev: ev.event_weight,
'FSRUp'              : lambda ev: ev.event_weight,
'HDAMPDown'          : lambda ev: ev.event_weight,
'HDAMPUp'            : lambda ev: ev.event_weight,
'ISRDown'            : lambda ev: ev.event_weight,
'ISRUp'              : lambda ev: ev.event_weight,
}

systs2_weights_objects = {
'JERDown':  lambda ev: ev.event_weight,
'JERUp'  :  lambda ev: ev.event_weight,
'JESDown':  lambda ev: ev.event_weight,
'JESUp'  :  lambda ev: ev.event_weight,
'TESDown':  lambda ev: ev.event_weight,
'TESUp'  :  lambda ev: ev.event_weight,
}


named_systs2_weights_all = {'nom': systs2_weights_nominal,
'common': systs2_weights_common,
'obj':    systs2_weights_objects,
'tt_weights': systs2_weights_tt,
'tt_hard':  systs2_weights_tt_hard,
'tt_pdf':   systs2_weights_tt_alpha_pdf,
'tt_alpha': systs2_weights_tt_alpha,
'tt_pdf1':  systs2_weights_tt_pdf1,
'tt_pdf10': systs2_weights_tt_pdf10,
'tt_pdf20': systs2_weights_tt_pdf20,
'tt_pdf30': systs2_weights_tt_pdf30,
'tt_pdf40': systs2_weights_tt_pdf40,
'tt_pdf50': systs2_weights_tt_pdf50,
}

systs2_weights_all = {}
for s_d in named_systs2_weights_all.values():
    systs2_weights_all.update(s_d)
systs2_weights_all.update(systs_weights_tt_updowns)
systs2_weights_all.update(systs2_weights_initial)


# ----- event-function weights
def event_function(param_name):
    return lambda ev: getattr(ev, param_name)

# dtags and systematics
std_dtag_systs = {
'data': (['SingleMuon', 'SingleElectron'], ['nom']),
'tt'  : (['MC2016_Summer16_TTJets_powheg'],  ["nom,common", "obj", "tt_weights", "tt_hard", "tt_pdf1", "tt_pdf10", "tt_pdf20", "tt_pdf30", "tt_pdf40", "tt_pdf50,tt_alpha"]), #select_sparse_channels
'other_mc': (['MC2016_Summer16_DYJetsToLL_10to50_amcatnlo',
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
'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg'], ['nom', 'common', 'obj']),

'qcd_mc': (['MC2016_Summer16_QCD_HT-100-200',
'MC2016_Summer16_QCD_HT-200-300',
'MC2016_Summer16_QCD_HT-300-500',
'MC2016_Summer16_QCD_HT-500-700',
'MC2016_Summer16_QCD_HT-700-1000',
'MC2016_Summer16_QCD_HT-1000-1500',
'MC2016_Summer16_QCD_HT-1500-2000',
'MC2016_Summer16_QCD_HT-2000-Inf'], ['nom'])
}




# selection stages and parameters per object sys variation
# object sys affects the thresholds of event, can move it between categories

systs_objects = {
'JERDown':  'selection_stage_JERDown',
'JERUp'  :  'selection_stage_JERUp'  ,
'JESDown':  'selection_stage_JESDown',
'JESUp'  :  'selection_stage_JESUp'  ,
'TESDown':  'selection_stage_TESDown',
'TESUp'  :  'selection_stage_TESUp'  ,
}

systs_objects_mt_variation = {
'JERUp'    : 'event_met_lep_mt_JERUp',   
'JERDown'  : 'event_met_lep_mt_JERDown', 
'JESUp'    : 'event_met_lep_mt_JESUp',   
'JESDown'  : 'event_met_lep_mt_JESDown', 
'TESUp'    : 'event_met_lep_mt_TESUp',   
'TESDown'  : 'event_met_lep_mt_TESDown', 
}

systs_objects_met_variation = {
'JERUp'    : 'event_met_JERUp.pt()',   
'JERDown'  : 'event_met_JERDown.pt()', 
'JESUp'    : 'event_met_JESUp.pt()',   
'JESDown'  : 'event_met_JESDown.pt()', 
'TESUp'    : 'event_met_TESUp.pt()',   
'TESDown'  : 'event_met_TESDown.pt()', 
}



# standard distrs

'''
distr_ranges = {
    'Mt_lep_met_c': '--custom-range 0,20,40,60,80,100,130,160,200,250',
    'Mt_lep_met_c2': '--custom-range 0,20,40,60,80,100,120,140,170,200,250,500',
    'Mt_lep_met_f':  '--histo-range 20,0,250',
    'met_f':         '--histo-range 20,0,300',
    'met_c':         '--custom-range 0,20,40,60,80,100,120,140,200,500',
    'dilep_mass':    '--histo-range 100,0,400',
    'lep_pt':        '--histo-range 20,0,150',
    'lep_eta':       '--histo-range 26,-2.6,2.6',
    'tau_sv_sign':   '--histo-range 42,-1,20',
    'tau_pt':        '--histo-range 20,0,100',
    'tau_eta':       '--histo-range 26,-2.6,2.6',
    'yield':         '--histo-range 3,0.0,3.0'
}

distrs_leptonic = [('std_mt_vars', 'Mt_lep_met_c'), ('std_mt_vars', 'Mt_lep_met_c2'), ('event_met_lep_mt', 'Mt_lep_met_f'), ('event_dilep_mass', 'dilep_mass'), ('event_leptons[0].pt()', 'lep_pt')]
distrs_tauonic_std  = [('std_mt_vars', 'Mt_lep_met_c'), ('std_mt_vars', 'Mt_lep_met_c2'), ('event_met_lep_mt', 'Mt_lep_met_f'), ('event_dilep_mass', 'dilep_mass'), ('event_leptons[0].pt()', 'lep_pt'),
        ('event_taus[0].pt()', 'tau_pt'), ('event_taus[0].eta()', 'tau_eta'),
        ('event_taus_sv_sign[0]', 'tau_sv_sign')]
'''

# a parameter can change due to object systamtics -- variation of objects in the event
# it happens for kinematic distributions and energy corrections
# variation of jet energy leads to change in MET and corresponding MET-related parameters

# calculation of standard distrs

systs_objects_mt_variation = {
'NOMINAL'  : lambda ev: ev.event_met_lep_mt,
'JERUp'    : lambda ev: ev.event_met_lep_mt_JERUp,
'JERDown'  : lambda ev: ev.event_met_lep_mt_JERDown,
'JESUp'    : lambda ev: ev.event_met_lep_mt_JESUp,
'JESDown'  : lambda ev: ev.event_met_lep_mt_JESDown,
'TESUp'    : lambda ev: ev.event_met_lep_mt_TESUp,
'TESDown'  : lambda ev: ev.event_met_lep_mt_TESDown,
}

systs_objects_met_variation = {
'NOMINAL'  : lambda ev: ev.event_met.pt(),
'JERUp'    : lambda ev: ev.event_met_JERUp.pt(),
'JERDown'  : lambda ev: ev.event_met_JERDown.pt(),
'JESUp'    : lambda ev: ev.event_met_JESUp.pt(),
'JESDown'  : lambda ev: ev.event_met_JESDown.pt(),
'TESUp'    : lambda ev: ev.event_met_TESUp.pt(),
'TESDown'  : lambda ev: ev.event_met_TESDown.pt(),
}

from ROOT import TMath

distr_defs = {
    'Mt_lep_met_c':  (systs_objects_mt_variation,  ('custom-range', [0,20,40,60,80,100,130,160,200,250])),
    'Mt_lep_met_c2': (systs_objects_mt_variation,  ('custom-range', [0,20,40,60,80,100,120,140,170,200,250,500])),
    'Mt_lep_met_f':  (systs_objects_mt_variation,  ('histo-range',  [20,0,250])),
    'Mt_lep_met_init_f': ({'NOMINAL': lambda ev: ev.event_met_lep_mt_init2}, ('histo-range',  [20,0,250])),
    'met_init_f':        ({'NOMINAL': lambda ev: ev.event_met_init2.pt()},   ('histo-range',  [30,0,300])),
    'met_f':         (systs_objects_met_variation, ('histo-range',  [30,0,300])),
    'met_c':         (systs_objects_met_variation, ('custom-range', [0,20,40,60,80,100,120,140,200,500])),
    'cos_phi_met_lep':   ({'NOMINAL': lambda ev: TMath.Cos(ev.event_leptons[0].phi() - ev.event_met.phi())},       ('histo-range',  [21,-1.05,1.05])),
    'phi_met_lep':       ({'NOMINAL': lambda ev: ev.event_leptons[0].phi() - ev.event_met.phi()},    ('histo-range',  [21,-2*3.15,2*3.15])),
    'dilep_mass':    ({'NOMINAL': lambda ev: ev.event_dilep_mass},       ('histo-range',  [100,0,400])),
    'dilep_mass_dy': ({'NOMINAL': lambda ev: ev.event_dilep_mass},       ('histo-range',  [40,80,100])),
    'lep_pt_f':      ({'NOMINAL': lambda ev: ev.event_leptons[0].pt()},  ('histo-range',  [20,0,150])),
    'lep_pt':        ({'NOMINAL': lambda ev: ev.event_leptons[0].pt()},  ('histo-range',  [40,0,200])),
    'lep_eta':       ({'NOMINAL': lambda ev: ev.event_leptons[0].eta()}, ('histo-range',  [26,-2.6,2.6])),
    'elmu_el_pt':    ({'NOMINAL': lambda ev: ev.event_leptons[0].pt()  if abs(ev.event_leptons_ids[0]) == 11 else (ev.event_leptons[1].pt()  if len(ev.event_leptons_ids)>1 else -111.)},  ('histo-range',  [40,0,200])),
    'elmu_el_eta':   ({'NOMINAL': lambda ev: ev.event_leptons[0].eta() if abs(ev.event_leptons_ids[0]) == 11 else (ev.event_leptons[1].eta() if len(ev.event_leptons_ids)>1 else -111.)}, ('histo-range',  [26,-2.6,2.6])),
    'elmu_mu_pt':    ({'NOMINAL': lambda ev: ev.event_leptons[0].pt()  if abs(ev.event_leptons_ids[0]) == 13 else (ev.event_leptons[1].pt()  if len(ev.event_leptons_ids)>1 else -111.)},  ('histo-range',  [40,0,200])),
    'elmu_mu_eta':   ({'NOMINAL': lambda ev: ev.event_leptons[0].eta() if abs(ev.event_leptons_ids[0]) == 13 else (ev.event_leptons[1].eta() if len(ev.event_leptons_ids)>1 else -111.)}, ('histo-range',  [26,-2.6,2.6])),
    'tau_sv_sign':   ({'NOMINAL': lambda ev: ev.event_taus_sv_sign[0] if len(ev.event_taus_sv_sign) > 0 else -111.},  ('histo-range',  [42,-1,20])),
    'tau_pt':        ({'NOMINAL': lambda ev: ev.event_taus[0].pt()    if len(ev.event_taus) > 0   else -111.},    ('histo-range',  [20,0,100])),
    'tau_eta':       ({'NOMINAL': lambda ev: ev.event_taus[0].eta()   if len(ev.event_taus) > 0   else -111.},    ('histo-range',  [26,-2.6,2.6])),
    'bjet_pt':       ({'NOMINAL': lambda ev: ev.event_jets_b[0].pt()  if len(ev.event_jets_b) > 0 else -111.},  ('histo-range',  [20,0,300])),
    'bjet_eta':      ({'NOMINAL': lambda ev: ev.event_jets_b[0].eta() if len(ev.event_jets_b) > 0 else -111.}, ('histo-range',  [26,-2.6,2.6])),
    'rjet_lead_pt':  ({'NOMINAL': lambda ev: ev.event_jets_r[0].pt()  if len(ev.event_jets_r) > 0 else -111.},  ('histo-range',  [20,0,300])),
    'rjet_lead_eta': ({'NOMINAL': lambda ev: ev.event_jets_r[0].eta() if len(ev.event_jets_r) > 0 else -111.}, ('histo-range',  [26,-2.6,2.6])),
    # TODO: implement saving all jets pT into 1 histo
    #'rjet_all_pt':   ({'NOMINAL': lambda ev: ev.event_jets_b[0].pt()  if len(ev.event_jets_b) > 0 else -111.},  ('histo-range',  [20,0,300])),
    #'rjet_all_eta':  ({'NOMINAL': lambda ev: ev.event_jets_b[0].eta() if len(ev.event_jets_b) > 0 else -111.}, ('histo-range',  [26,-2.6,2.6])),

    'nbjets':        ({'NOMINAL': lambda ev: ev.event_jets_n_bjets},                        ('histo-range',  [5,0.0,5.0])),
    'nrjets':        ({'NOMINAL': lambda ev: ev.event_jets_n_jets - ev.event_jets_n_bjets}, ('histo-range',  [6,0.0,6.0])),
    'najets':        ({'NOMINAL': lambda ev: ev.event_jets_n_jets},                         ('histo-range',  [10,0.0,10.0])),
    'relIso_el':     ({'NOMINAL': lambda ev: ev.event_leptons_alliso_reliso[0] if len(ev.event_leptons_alliso_reliso) > 0 else -111.},  ('custom-range', [0, 0.059, 0.1, 0.175, 0.3, 0.5])),
    'relIso_el_ext': ({'NOMINAL': lambda ev: ev.event_leptons_alliso_reliso[0] if len(ev.event_leptons_alliso_reliso) > 0 else -111.},  ('custom-range', [0, 0.059, 0.1, 0.175, 0.3, 0.5, 2., 4.])),
    'relIso_mu':     ({'NOMINAL': lambda ev: ev.event_leptons_alliso_reliso[0] if len(ev.event_leptons_alliso_reliso) > 0 else -111.},  ('custom-range', [0, 0.15, 0.25, 0.5, 1., 2., 4.])),
    'relIso_mu_ext': ({'NOMINAL': lambda ev: ev.event_leptons_alliso_reliso[0] if len(ev.event_leptons_alliso_reliso) > 0 else -111.},  ('custom-range', [0, 0.15, 0.25, 0.5, 1., 2., 4., 8., 16.])),
    'yield':         ({'NOMINAL': lambda ev: 1},     ('histo-range',  [3,0.0,3.0])),
}

# CHECK: implement multi-dim histos?

def make_histo(name, histo_range):
    import ctypes
    from ROOT import TH1D
    range_type, bin_edges = histo_range

    if range_type == 'custom-range':
        #bin_edges = [float(b) for b in args.custom_range.split(',')]
        n_bin_edges = len(bin_edges)
        n_bins = n_bin_edges - 1
        logging.debug("making %d bins from %s" % (n_bins, repr(bin_edges)))
        root_bin_edges = (ctypes.c_double * n_bin_edges)(*bin_edges)
        output_histo = TH1D(name, "", n_bins, root_bin_edges)

    elif range_type == 'histo-range':
        n_bins, lbin, rbin = bin_edges
        output_histo = TH1D(name, "", n_bins, lbin, rbin)
    else:
        raise ValueError('UNKNOWN range type %s in %s' % (range_type, name))

    return output_histo

# also a selection stage can change with object variations, events can move between categories
# now it is implemented only for main selection (no control selections have it)

# calculation of standard channels

main_sel_stages = {
'NOMINAL':  lambda ev: ev.selection_stage,
'JERDown':  lambda ev: ev.selection_stage_JERDown,
'JERUp'  :  lambda ev: ev.selection_stage_JERUp  ,
'JESDown':  lambda ev: ev.selection_stage_JESDown,
'JESUp'  :  lambda ev: ev.selection_stage_JESUp  ,
'TESDown':  lambda ev: ev.selection_stage_TESDown,
'TESUp'  :  lambda ev: ev.selection_stage_TESUp  ,
}

main_presel_stages = {'NOMINAL': lambda ev: ev.selection_stage_presel}

#def main_sel_stages(sname, event):
#    if sname in main_selection_stages:
#        return main_selection_stages[sname](ev)
#    else:
#        return main_selection_stages['NOMINAL'](ev)

# new selection stages
std_channels_ev_loop = {
'mu_selTight':     (lambda sel_stage, ev: (sel_stage==  9), main_sel_stages),
'mu_selTight_ss':  (lambda sel_stage, ev: (sel_stage==  8), main_sel_stages),
'el_selTight':     (lambda sel_stage, ev: (sel_stage== 19), main_sel_stages),
'el_selTight_ss':  (lambda sel_stage, ev: (sel_stage== 18), main_sel_stages),
'mu_sel':          (lambda sel_stage, ev: (sel_stage==  9 or sel_stage==  7), main_sel_stages),
'mu_sel_ss':       (lambda sel_stage, ev: (sel_stage==  8 or sel_stage==  6), main_sel_stages),
'el_sel':          (lambda sel_stage, ev: (sel_stage== 19 or sel_stage== 17), main_sel_stages),
'el_sel_ss':       (lambda sel_stage, ev: (sel_stage== 18 or sel_stage== 16), main_sel_stages),
'mu_selVloose':    (lambda sel_stage, ev: (sel_stage==  9 or sel_stage==  7 or sel_stage ==  5), main_sel_stages),
'mu_selVloose_ss': (lambda sel_stage, ev: (sel_stage==  8 or sel_stage==  6 or sel_stage ==  4), main_sel_stages),
'el_selVloose':    (lambda sel_stage, ev: (sel_stage== 19 or sel_stage== 17 or sel_stage == 15), main_sel_stages),
'el_selVloose_ss': (lambda sel_stage, ev: (sel_stage== 18 or sel_stage== 16 or sel_stage == 14), main_sel_stages),

'mu_selTight_ljout':     (lambda sel_stage, ev: (sel_stage==  9) and ev.event_jets_lj_var >  60., main_sel_stages),
'mu_selTight_ljout_ss':  (lambda sel_stage, ev: (sel_stage==  8) and ev.event_jets_lj_var >  60., main_sel_stages),
'el_selTight_ljout':     (lambda sel_stage, ev: (sel_stage== 19) and ev.event_jets_lj_var >  60., main_sel_stages),
'el_selTight_ljout_ss':  (lambda sel_stage, ev: (sel_stage== 18) and ev.event_jets_lj_var >  60., main_sel_stages),
'mu_sel_ljout':          (lambda sel_stage, ev: (sel_stage==  9 or sel_stage==  7) and ev.event_jets_lj_var >  60., main_sel_stages),
'mu_sel_ljout_ss':       (lambda sel_stage, ev: (sel_stage==  8 or sel_stage==  6) and ev.event_jets_lj_var >  60., main_sel_stages),
'el_sel_ljout':          (lambda sel_stage, ev: (sel_stage== 19 or sel_stage== 17) and ev.event_jets_lj_var >  60., main_sel_stages),
'el_sel_ljout_ss':       (lambda sel_stage, ev: (sel_stage== 18 or sel_stage== 16) and ev.event_jets_lj_var >  60., main_sel_stages),
'mu_selVloose_ljout':    (lambda sel_stage, ev: (sel_stage==  9 or sel_stage==  7 or sel_stage==  5) and ev.event_jets_lj_var >  60., main_sel_stages),
'mu_selVloose_ljout_ss': (lambda sel_stage, ev: (sel_stage==  8 or sel_stage==  6 or sel_stage==  4) and ev.event_jets_lj_var >  60., main_sel_stages),
'el_selVloose_ljout':    (lambda sel_stage, ev: (sel_stage== 19 or sel_stage== 17 or sel_stage== 15) and ev.event_jets_lj_var >  60., main_sel_stages),
'el_selVloose_ljout_ss': (lambda sel_stage, ev: (sel_stage== 18 or sel_stage== 16 or sel_stage== 14) and ev.event_jets_lj_var >  60., main_sel_stages),

'mu_selTight_lj':     (lambda sel_stage, ev: (sel_stage==  9) and ev.event_jets_lj_var <= 60., main_sel_stages),
'mu_selTight_lj_ss':  (lambda sel_stage, ev: (sel_stage==  8) and ev.event_jets_lj_var <= 60., main_sel_stages),
'el_selTight_lj':     (lambda sel_stage, ev: (sel_stage== 19) and ev.event_jets_lj_var <= 60., main_sel_stages),
'el_selTight_lj_ss':  (lambda sel_stage, ev: (sel_stage== 18) and ev.event_jets_lj_var <= 60., main_sel_stages),
'mu_sel_lj':          (lambda sel_stage, ev: (sel_stage==  9 or sel_stage==  7) and ev.event_jets_lj_var <= 60., main_sel_stages),
'mu_sel_lj_ss':       (lambda sel_stage, ev: (sel_stage==  8 or sel_stage==  6) and ev.event_jets_lj_var <= 60., main_sel_stages),
'el_sel_lj':          (lambda sel_stage, ev: (sel_stage== 19 or sel_stage== 17) and ev.event_jets_lj_var <= 60., main_sel_stages),
'el_sel_lj_ss':       (lambda sel_stage, ev: (sel_stage== 18 or sel_stage== 16) and ev.event_jets_lj_var <= 60., main_sel_stages),
'mu_selVloose_lj':    (lambda sel_stage, ev: (sel_stage==  9 or sel_stage==  7 or sel_stage==  5) and ev.event_jets_lj_var <= 60., main_sel_stages),
'mu_selVloose_lj_ss': (lambda sel_stage, ev: (sel_stage==  8 or sel_stage==  6 or sel_stage==  4) and ev.event_jets_lj_var <= 60., main_sel_stages),
'el_selVloose_lj':    (lambda sel_stage, ev: (sel_stage== 19 or sel_stage== 17 or sel_stage== 15) and ev.event_jets_lj_var <= 60., main_sel_stages),
'el_selVloose_lj_ss': (lambda sel_stage, ev: (sel_stage== 18 or sel_stage== 16 or sel_stage== 14) and ev.event_jets_lj_var <= 60., main_sel_stages),

'mu_presel':          (lambda sel_stage, ev: (sel_stage >  1 and sel_stage < 10), main_presel_stages),
'mu_presel_lj':       (lambda sel_stage, ev: (sel_stage >  1 and sel_stage < 10) and ev.event_jets_lj_var <= 60., main_presel_stages),
'mu_presel_ljout':    (lambda sel_stage, ev: (sel_stage >  1 and sel_stage < 10) and ev.event_jets_lj_var  > 60., main_presel_stages),
'el_presel':          (lambda sel_stage, ev: (sel_stage > 10 and sel_stage < 20), main_presel_stages),
'el_presel_lj':       (lambda sel_stage, ev: (sel_stage > 10 and sel_stage < 20) and ev.event_jets_lj_var <= 60., main_presel_stages),
'el_presel_ljout':    (lambda sel_stage, ev: (sel_stage > 10 and sel_stage < 20) and ev.event_jets_lj_var  > 60., main_presel_stages),

'mu_preselCand':      (lambda sel_stage, ev: (sel_stage==  9), main_presel_stages),
'mu_preselCand_ss':   (lambda sel_stage, ev: (sel_stage==  8), main_presel_stages),
'el_preselCand':      (lambda sel_stage, ev: (sel_stage== 19), main_presel_stages),
'el_preselCand_ss':   (lambda sel_stage, ev: (sel_stage== 18), main_presel_stages),
'mu_preselCand_lj':      (lambda sel_stage, ev: (sel_stage==  9 and ev.event_jets_lj_var <= 60.), main_presel_stages),
'mu_preselCand_lj_ss':   (lambda sel_stage, ev: (sel_stage==  8 and ev.event_jets_lj_var <= 60.), main_presel_stages),
'el_preselCand_lj':      (lambda sel_stage, ev: (sel_stage== 19 and ev.event_jets_lj_var <= 60.), main_presel_stages),
'el_preselCand_lj_ss':   (lambda sel_stage, ev: (sel_stage== 18 and ev.event_jets_lj_var <= 60.), main_presel_stages),
'mu_preselCand_ljout':      (lambda sel_stage, ev: (sel_stage==  9 and ev.event_jets_lj_var  > 60.), main_presel_stages),
'mu_preselCand_ljout_ss':   (lambda sel_stage, ev: (sel_stage==  8 and ev.event_jets_lj_var  > 60.), main_presel_stages),
'el_preselCand_ljout':      (lambda sel_stage, ev: (sel_stage== 19 and ev.event_jets_lj_var  > 60.), main_presel_stages),
'el_preselCand_ljout_ss':   (lambda sel_stage, ev: (sel_stage== 18 and ev.event_jets_lj_var  > 60.), main_presel_stages),

# additional channels
# TODO: made the met_lep_mt variation with sys objects
'dy_mutau': (lambda sel_stage, ev: ((sel_stage== 102 or sel_stage== 103) and ev.event_met_lep_mt < 40.), {'NOMINAL': lambda ev: ev.selection_stage_dy}),
'dy_eltau': (lambda sel_stage, ev: ((sel_stage== 112 or sel_stage== 113) and ev.event_met_lep_mt < 40.), {'NOMINAL': lambda ev: ev.selection_stage_dy}),
'dy_mutau_3j': (lambda sel_stage, ev: ((sel_stage== 102 or sel_stage== 103) and ev.event_met_lep_mt < 40. and ev.event_jets_n_jets > 2), {'NOMINAL': lambda ev: ev.selection_stage_dy}),
'dy_eltau_3j': (lambda sel_stage, ev: ((sel_stage== 112 or sel_stage== 113) and ev.event_met_lep_mt < 40. and ev.event_jets_n_jets > 2), {'NOMINAL': lambda ev: ev.selection_stage_dy}),

'dy_mutau_ss':    (lambda sel_stage, ev: ((sel_stage== 202 or sel_stage== 203) and ev.event_met_lep_mt < 40.), {'NOMINAL': lambda ev: ev.selection_stage_dy}),
'dy_eltau_ss':    (lambda sel_stage, ev: ((sel_stage== 212 or sel_stage== 213) and ev.event_met_lep_mt < 40.), {'NOMINAL': lambda ev: ev.selection_stage_dy}),
'dy_mutau_3j_ss': (lambda sel_stage, ev: ((sel_stage== 202 or sel_stage== 203) and ev.event_met_lep_mt < 40. and ev.event_jets_n_jets > 2), {'NOMINAL': lambda ev: ev.selection_stage_dy}),
'dy_eltau_3j_ss': (lambda sel_stage, ev: ((sel_stage== 212 or sel_stage== 213) and ev.event_met_lep_mt < 40. and ev.event_jets_n_jets > 2), {'NOMINAL': lambda ev: ev.selection_stage_dy}),

'dy_mumu':  (lambda sel_stage, ev: (sel_stage== 102 or sel_stage== 103 or sel_stage== 105), {'NOMINAL': lambda ev: ev.selection_stage_dy_mumu}),
'dy_elel':  (lambda sel_stage, ev: (sel_stage== 112 or sel_stage== 113 or sel_stage== 115), {'NOMINAL': lambda ev: ev.selection_stage_dy_mumu}),

'wjets_mu'    : (lambda sel_stage, ev: (sel_stage==  7), {'NOMINAL': lambda ev: ev.selection_stage_wjets}),
'wjets_mu_ss' : (lambda sel_stage, ev: (sel_stage==  6), {'NOMINAL': lambda ev: ev.selection_stage_wjets}),
'wjets_el'    : (lambda sel_stage, ev: (sel_stage== 17), {'NOMINAL': lambda ev: ev.selection_stage_wjets}),
'wjets_el_ss' : (lambda sel_stage, ev: (sel_stage== 16), {'NOMINAL': lambda ev: ev.selection_stage_wjets}),

'tt_elmu':  (lambda sel_stage, ev: (sel_stage> 200 and sel_stage < 210 and ev.event_leptons[0].pt() > 32. and ev.event_leptons[1].pt() > 32.), {'NOMINAL': lambda ev: ev.selection_stage_em}),
'tt_alliso_presel_el'    :  (lambda sel_stage, ev: sel_stage==513, {'NOMINAL': lambda ev: ev.selection_stage_tt_alliso}),
'tt_alliso_presel_el_ss' :  (lambda sel_stage, ev: sel_stage==512, {'NOMINAL': lambda ev: ev.selection_stage_tt_alliso}),
'tt_alliso_presel_mu'    :  (lambda sel_stage, ev: sel_stage==503, {'NOMINAL': lambda ev: ev.selection_stage_tt_alliso}),
'tt_alliso_presel_mu_ss' :  (lambda sel_stage, ev: sel_stage==502, {'NOMINAL': lambda ev: ev.selection_stage_tt_alliso}),
}
# calculation of standard systematic weights

more_channels_ev_loop = {
'mu_selSV':          (lambda sel_stage, ev: (sel_stage==  9 or sel_stage==  7) and ev.event_taus_sv_sign[0] > 2.5, main_sel_stages),
'mu_selSV_ss':       (lambda sel_stage, ev: (sel_stage==  8 or sel_stage==  6) and ev.event_taus_sv_sign[0] > 2.5, main_sel_stages),
'el_selSV':          (lambda sel_stage, ev: (sel_stage== 19 or sel_stage== 17) and ev.event_taus_sv_sign[0] > 2.5, main_sel_stages),
'el_selSV_ss':       (lambda sel_stage, ev: (sel_stage== 18 or sel_stage== 16) and ev.event_taus_sv_sign[0] > 2.5, main_sel_stages),
'mu_selSVVloose':    (lambda sel_stage, ev: (sel_stage==  9 or sel_stage==  7 or sel_stage==  5) and ev.event_taus_sv_sign[0] > 2.5, main_sel_stages),
'mu_selSVVloose_ss': (lambda sel_stage, ev: (sel_stage==  8 or sel_stage==  6 or sel_stage==  4) and ev.event_taus_sv_sign[0] > 2.5, main_sel_stages),
'el_selSVVloose':    (lambda sel_stage, ev: (sel_stage== 19 or sel_stage== 17 or sel_stage== 15) and ev.event_taus_sv_sign[0] > 2.5, main_sel_stages),
'el_selSVVloose_ss': (lambda sel_stage, ev: (sel_stage== 18 or sel_stage== 16 or sel_stage== 14) and ev.event_taus_sv_sign[0] > 2.5, main_sel_stages),
}

all_channels_ev_loop = {}
all_channels_ev_loop.update(more_channels_ev_loop)
all_channels_ev_loop.update(std_channels_ev_loop)




# the standard distrs in channels

distrs_mt_fit   = {'Mt_lep_met_c'}
distrs_mt       = {'Mt_lep_met_c', 'Mt_lep_met_c2', 'Mt_lep_met_f'}
distrs_mt_calc  = {'lep_pt', 'met_c', 'phi_met_lep', 'cos_phi_met_lep'}
distrs_leptonic = {'Mt_lep_met_c', 'Mt_lep_met_c2', 'Mt_lep_met_f', 'dilep_mass', 'lep_pt', 'lep_eta', 'met_c'}
distrs_lep      = {'Mt_lep_met_c', 'lep_pt', 'lep_eta', 'met_c'}
distrs_relIso   = {'relIso_mu', 'relIso_mu_ext', 'relIso_el', 'relIso_el_ext'}
distrs_dy       = {'met_c',          'Mt_lep_met_c', 'Mt_lep_met_f', 'dilep_mass', 'dilep_mass_dy', 'lep_pt', 'lep_eta'}
distrs_wjets    = {'met_c', 'met_f', 'met_init_f', 'Mt_lep_met_c', 'Mt_lep_met_f', 'Mt_lep_met_init_f', 'dilep_mass',                  'lep_pt', 'lep_eta', 'phi_met_lep', 'cos_phi_met_lep'}

distrs_on_jets = {'nbjets', 'nrjets', 'bjet_pt', 'bjet_eta', 'rjet_lead_pt', 'rjet_lead_eta'}

distrs_tauonic_std  = {'tau_pt', 'tau_eta', 'tau_sv_sign'}
#distrs_tauonic_nomt  = distrs_tauonic_std - {'Mt_lep_met_c', 'Mt_lep_met_c2', 'Mt_lep_met_f'}
#distrs_tauonic = distrs_tauonic_nomt


main_sys = ['nom', 'common', 'obj']
full_sys = ['nom', 'common', 'obj', 'tt_weights', "tt_hard", 'tt_pdf']
presel_sys = ['NOMINAL', 'TOPPTDown', 'TOPPTUp', 'PUUp', 'PUDown']

channels_distrs = {
'tt_alliso_presel'    : (['tt_alliso_presel_el', 'tt_alliso_presel_el_ss', 'tt_alliso_presel_mu', 'tt_alliso_presel_mu_ss', ], sorted(distrs_lep.union(distrs_relIso)), main_sys),
'tt_dileptons'    : (['tt_elmu'], sorted(distrs_leptonic.union({'phi_met_lep'}).union({'elmu_el_pt', 'elmu_el_eta', 'elmu_mu_pt', 'elmu_mu_eta'})), main_sys + ['TOPPTUp', 'TOPPTDown']),
'tt_leptauSV'     : (['el_selSV', 'el_selSVVloose', 'el_selSV_ss', 'el_selSVVloose_ss', 'mu_selSV', 'mu_selSVVloose', 'mu_selSV_ss', 'mu_selSVVloose_ss'], sorted(distrs_tauonic_std.union(distrs_leptonic) - distrs_mt_fit), ['nom']),

'fit_tt_leptau'       : (['el_sel',       'el_sel_ss',       'mu_sel',       'mu_sel_ss'],       sorted(distrs_mt), full_sys),
'fit_tt_leptau_lj'    : (['el_sel_lj',    'el_sel_lj_ss',    'mu_sel_lj',    'mu_sel_lj_ss'],    sorted(distrs_mt), full_sys),
'fit_tt_leptau_ljout' : (['el_sel_ljout', 'el_sel_ljout_ss', 'mu_sel_ljout', 'mu_sel_ljout_ss'], sorted(distrs_mt), full_sys),

'fit_tt_leptau_Vloose'       : (['el_selVloose',       'el_selVloose_ss',       'mu_selVloose',       'mu_selVloose_ss'],        sorted(distrs_mt), full_sys),
'fit_tt_leptau_Vloose_lj'    : (['el_selVloose_lj',    'el_selVloose_lj_ss',    'mu_selVloose_lj',    'mu_selVloose_lj_ss'],     sorted(distrs_mt), full_sys),
'fit_tt_leptau_Vloose_ljout' : (['el_selVloose_ljout', 'el_selVloose_ljout_ss', 'mu_selVloose_ljout', 'mu_selVloose_ljout_ss'],  sorted(distrs_mt), full_sys),

'fit_tt_leptau_Tight'       : (['el_selTight',       'el_selTight_ss',       'mu_selTight',       'mu_selTight_ss'],        sorted(distrs_mt), full_sys),
'fit_tt_leptau_Tight_lj'    : (['el_selTight_lj',    'el_selTight_lj_ss',    'mu_selTight_lj',    'mu_selTight_lj_ss'],     sorted(distrs_mt), full_sys),
'fit_tt_leptau_Tight_ljout' : (['el_selTight_ljout', 'el_selTight_ljout_ss', 'mu_selTight_ljout', 'mu_selTight_ljout_ss'],  sorted(distrs_mt), full_sys),

'tt_leptau'       : (['el_sel',       'el_sel_ss',       'mu_sel',       'mu_sel_ss'],       sorted(distrs_tauonic_std.union(distrs_leptonic, distrs_on_jets, distrs_mt_calc) - distrs_mt_fit), main_sys),
'tt_leptau_lj'    : (['el_sel_lj',    'el_sel_lj_ss',    'mu_sel_lj',    'mu_sel_lj_ss'],    sorted(distrs_tauonic_std.union(distrs_leptonic, distrs_mt_calc) - distrs_mt_fit), main_sys),
'tt_leptau_ljout' : (['el_sel_ljout', 'el_sel_ljout_ss', 'mu_sel_ljout', 'mu_sel_ljout_ss'], sorted(distrs_tauonic_std.union(distrs_leptonic, distrs_mt_calc) - distrs_mt_fit), main_sys),

'tt_leptau_Vloose'       : (['el_selVloose',       'el_selVloose_ss',       'mu_selVloose',       'mu_selVloose_ss'],        sorted(distrs_tauonic_std.union(distrs_leptonic) - distrs_mt_fit), main_sys),
'tt_leptau_Vloose_lj'    : (['el_selVloose_lj',    'el_selVloose_lj_ss',    'mu_selVloose_lj',    'mu_selVloose_lj_ss'],     sorted(distrs_tauonic_std.union(distrs_leptonic) - distrs_mt_fit), main_sys),
'tt_leptau_Vloose_ljout' : (['el_selVloose_ljout', 'el_selVloose_ljout_ss', 'mu_selVloose_ljout', 'mu_selVloose_ljout_ss'],  sorted(distrs_tauonic_std.union(distrs_leptonic) - distrs_mt_fit), main_sys),

'tt_leptau_Tight'       : (['el_selTight',       'el_selTight_ss',       'mu_selTight',       'mu_selTight_ss'],        sorted(distrs_tauonic_std.union(distrs_leptonic) - distrs_mt_fit), main_sys),
'tt_leptau_Tight_lj'    : (['el_selTight_lj',    'el_selTight_lj_ss',    'mu_selTight_lj',    'mu_selTight_lj_ss'],     sorted(distrs_tauonic_std.union(distrs_leptonic) - distrs_mt_fit), main_sys),
'tt_leptau_Tight_ljout' : (['el_selTight_ljout', 'el_selTight_ljout_ss', 'mu_selTight_ljout', 'mu_selTight_ljout_ss'],  sorted(distrs_tauonic_std.union(distrs_leptonic) - distrs_mt_fit), main_sys),

'tt_presel_cands'       : (['el_preselCand', 'el_preselCand_ss', 'mu_preselCand', 'mu_preselCand_ss',], sorted(distrs_tauonic_std.union(distrs_leptonic, distrs_on_jets) - distrs_mt_fit), presel_sys),
'tt_presel_cands_lj'    : (['el_preselCand_lj', 'el_preselCand_lj_ss', 'mu_preselCand_lj', 'mu_preselCand_lj_ss',], sorted(distrs_tauonic_std.union(distrs_leptonic, distrs_on_jets) - distrs_mt_fit), presel_sys),
'tt_presel_cands_ljout' : (['el_preselCand_ljout', 'el_preselCand_ljout_ss', 'mu_preselCand_ljout', 'mu_preselCand_ljout_ss',], sorted(distrs_tauonic_std.union(distrs_leptonic, distrs_on_jets) - distrs_mt_fit), presel_sys),
'tt_presel_lj_mu' : (['mu_presel',    'mu_presel_lj',    'mu_presel_ljout',], sorted(distrs_on_jets.union(distrs_lep)), presel_sys),
'tt_presel_lj_el' : (['el_presel',    'el_presel_lj',    'el_presel_ljout',], sorted(distrs_on_jets.union(distrs_lep)), presel_sys),

'dy_dileptons'    : (['dy_mumu',  'dy_elel'],   sorted(distrs_dy     ), main_sys),
'dy_leptau'       : (['dy_mutau', 'dy_eltau', 'dy_mutau_3j', 'dy_eltau_3j'],  sorted(distrs_tauonic_std.union(distrs_leptonic).union(distrs_on_jets).union(distrs_mt_calc)), main_sys),
'dy_leptau_ss'    : (['dy_mutau_ss', 'dy_eltau_ss', 'dy_mutau_3j_ss', 'dy_eltau_3j_ss'],  sorted(distrs_tauonic_std.union(distrs_leptonic).union(distrs_on_jets).union(distrs_mt_calc)), main_sys),
'wjets'           : (['wjets_mu', 'wjets_el'],  sorted(distrs_wjets  ), main_sys),
'wjets_ss'        : (['wjets_mu_ss', 'wjets_el_ss'],  sorted(distrs_wjets  ), main_sys),

'yields_el_presel' : (['el_preselCand', 'el_preselCand_ss', 'el_presel',], sorted({'yield'}), ['NOMINAL']),
'yields_el_sel'    : (['el_sel', 'el_sel_ss', 'el_sel_lj', 'el_sel_lj_ss', 'el_sel_ljout', 'el_sel_ljout_ss'], sorted({'yield'}), ['NOMINAL']),
'yields_mu_presel' : (['mu_preselCand', 'mu_preselCand_ss', 'mu_presel',], sorted({'yield'}), ['NOMINAL']),
'yields_mu_sel'    : (['mu_sel', 'mu_sel_ss', 'mu_sel_lj', 'mu_sel_lj_ss', 'mu_sel_ljout', 'mu_sel_ljout_ss'], sorted({'yield'}), ['NOMINAL']),

'yields_control_wjets'    : (['wjets_mu', 'wjets_el', 'wjets_mu_ss', 'wjets_el_ss'], sorted({'yield'}), ['NOMINAL']),
'yields_control_dyjets'   : (['dy_mumu',  'dy_elel'], sorted({'yield'}), ['NOMINAL']),
'yields_control_tt_em'    : (['tt_elmu',], sorted({'yield'}), ['NOMINAL']),
}


# the standard systematics per dtag

sample_info = {
'data': ([
'Data13TeV_SingleElectron2016B_03Feb2017_ver2',
'Data13TeV_SingleElectron2016C_03Feb2017_v1',
'Data13TeV_SingleElectron2016D_03Feb2017_v1',
'Data13TeV_SingleElectron2016E_03Feb2017_v1',
'Data13TeV_SingleElectron2016F_03Feb2017_v1',
'Data13TeV_SingleElectron2016G_03Feb2017_v1',
'Data13TeV_SingleElectron2016H_03Feb2017_ver2',
'Data13TeV_SingleElectron2016H_03Feb2017_ver3',
'Data13TeV_SingleMuon2016B_03Feb2017_ver2',
'Data13TeV_SingleMuon2016C_03Feb2017_v1',
'Data13TeV_SingleMuon2016D_03Feb2017_v1',
'Data13TeV_SingleMuon2016E_03Feb2017_v1',
'Data13TeV_SingleMuon2016F_03Feb2017_v1',
'Data13TeV_SingleMuon2016G_03Feb2017_v1',
'Data13TeV_SingleMuon2016H_03Feb2017_ver2',
'Data13TeV_SingleMuon2016H_03Feb2017_ver3'], ['nom']),

'tt'  : (['MC2016_Summer16_TTJets_powheg'],
  full_sys), # "tt_pdf1", "tt_pdf10", "tt_pdf20", "tt_pdf30", "tt_pdf40", "tt_pdf50,tt_alpha"]),

'tt_syst': ([
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4down',
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4down_ext1',
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4up',
'MC2016_Summer16_TTJets_powheg_CUETP8M2T4up_ext1',
'MC2016_Summer16_TTJets_powheg_fsrdown',
'MC2016_Summer16_TTJets_powheg_fsrdown_ext1',
'MC2016_Summer16_TTJets_powheg_fsrdown_ext2',
'MC2016_Summer16_TTJets_powheg_fsrup',
'MC2016_Summer16_TTJets_powheg_fsrup_ext1',
'MC2016_Summer16_TTJets_powheg_fsrup_ext2',
'MC2016_Summer16_TTJets_powheg_hdampDOWN',
'MC2016_Summer16_TTJets_powheg_hdampDOWN_ext1',
'MC2016_Summer16_TTJets_powheg_hdampUP',
'MC2016_Summer16_TTJets_powheg_hdampUP_ext1',
'MC2016_Summer16_TTJets_powheg_isrdown',
'MC2016_Summer16_TTJets_powheg_isrdown_ext1',
'MC2016_Summer16_TTJets_powheg_isrdown_ext2',
'MC2016_Summer16_TTJets_powheg_isrup'],
  ['nom']),

'other_mc': ([
'MC2016_Summer16_DYJetsToLL_10to50_amcatnlo',
'MC2016_Summer16_DYJetsToLL_10to50_amcatnlo_v1_ext1',
'MC2016_Summer16_DYJetsToLL_10to50_amcatnlo_v2',
'MC2016_Summer16_DYJetsToLL_50toInf_madgraph',
'MC2016_Summer16_DYJetsToLL_50toInf_madgraph_ext2_v1',

'MC2016_Summer16_W1Jets_madgraph',
'MC2016_Summer16_W2Jets_madgraph',
'MC2016_Summer16_W2Jets_madgraph_ext1',
'MC2016_Summer16_W3Jets_madgraph',
'MC2016_Summer16_W3Jets_madgraph_ext1',
'MC2016_Summer16_W4Jets_madgraph',
'MC2016_Summer16_W4Jets_madgraph_ext1',
'MC2016_Summer16_W4Jets_madgraph_ext2',
'MC2016_Summer16_WJets_madgraph',
'MC2016_Summer16_WJets_madgraph_ext2_v1',

'MC2016_Summer16_WJets_amcatnlo',
'MC2016_Summer16_WJets_amcatnlo_ext2_v2',

'MC2016_Summer16_WWTo2L2Nu_powheg',
'MC2016_Summer16_WWToLNuQQ_powheg',
'MC2016_Summer16_WZTo1L1Nu2Q_amcatnlo_madspin',
'MC2016_Summer16_WZTo1L3Nu_amcatnlo_madspin',
'MC2016_Summer16_WZTo2L2Q_amcatnlo_madspin',
'MC2016_Summer16_WZTo3LNu_powheg',
'MC2016_Summer16_ZZTo2L2Nu_powheg',
'MC2016_Summer16_ZZTo2L2Q_amcatnlo_madspin',
'MC2016_Summer16_SingleT_tW_5FS_powheg',
'MC2016_Summer16_SingleTbar_tW_5FS_powheg',
'MC2016_Summer16_schannel_4FS_leptonicDecays_amcatnlo',
'MC2016_Summer16_tchannel_antitop_4f_leptonicDecays_powheg',
'MC2016_Summer16_tchannel_top_4f_leptonicDecays_powheg'], main_sys),

'qcd_mc': (['MC2016_Summer16_QCD_HT-100-200',
'MC2016_Summer16_QCD_HT-200-300',
'MC2016_Summer16_QCD_HT-300-500',
'MC2016_Summer16_QCD_HT-500-700',
'MC2016_Summer16_QCD_HT-700-1000',
'MC2016_Summer16_QCD_HT-1000-1500',
'MC2016_Summer16_QCD_HT-1500-2000',
'MC2016_Summer16_QCD_HT-2000-Inf'], ['nom'])
}




