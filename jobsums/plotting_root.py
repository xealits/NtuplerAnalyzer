from ROOT import THStack, TLegend, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor

def rgb(r, g, b):
    '''rgb(r, g, b):

    from:
    http://webhome.phy.duke.edu/~dmb60/the-guide/
    TColor* color = gROOT->GetColor(TColor::GetColor(red,green,blue));//Use ints from 0 to 255 
    color->SetAlpha(0.5);//0 is fully transparent, 1 fully opaque
    hist->SetFillColor(color->GetNumber());
    '''
    return TColor.GetColor(r, g, b)

nick_info = {
"data":        {'color': kWhite,   'legend': 'data'},
"dy":          {'color': kGray,    'legend': 'DY'},
"dy_other":    {'color': kGray,    'legend': 'DY #rightarrow other'},
"dy_tautau":   {'color': kGray+2,  'legend': 'DY #rightarrow #tau#tau'},
"wjets":       {'color': kCyan+2,  'legend': 'W+jets#rightarrow l'},
"wjets_tauh":  {'color': kCyan+3,  'legend': 'W+jets#rightarrow #tau_{h}'},
"wjets_taul":  {'color': kCyan+4,  'legend': 'W+jets#rightarrow #tau_{l}'},
"dibosons":    {'color': kCyan,    'legend': 'dibosons'},
"singletop":   {'color': kAzure,   'legend': 'singletop'},
"s_top_eltau": {'color': kAzure,   'legend': 's.top#rightarrow e#tau'},
"s_top_mutau": {'color': kAzure,   'legend': 's.top#rightarrow #mu#tau'},
"s_top_elmu":  {'color': kAzure+1, 'legend': 's.top#rightarrow e#mu'},
"s_top_lj":    {'color': kAzure+2, 'legend': 's.top#rightarrow lj'},
"s_top_other": {'color': kAzure+3, 'legend': 's.top#rightarrow other'},

"tt_jj":         {'color': kGreen+4,  'legend': 'tt_jj'},
"tt_em":         {'color': kYellow-7, 'legend': 'tt_em'},
"tt_{em}":       {'color': kYellow-7, 'legend': 'tt_{em}'},
"tt_ee":         {'color': kAzure-9,  'legend': 'tt_ee'},
"tt_mm":         {'color': kGreen-9,  'legend': 'tt_mm'},
"tt_{l\\tau-l}": {'color': kOrange+3, 'legend': 'tt_{l\\tau-l}'},

"tt_mutau3ch":   {'color': rgb(255, 255, 178), 'legend': 't#bar{t}#rightarrow #mu#tau_{3h}'},
"tt_eltau3ch":   {'color': rgb(254,204,92),    'legend': 't#bar{t}#rightarrow e#tau_{3h}'},

"tt_mutau":      {'color': rgb(255, 255, 178), 'legend': 't#bar{t}#rightarrow #mu#tau_{h}'},
"tt_eltau":      {'color': rgb(254,204,92),    'legend': 't#bar{t}#rightarrow e#tau_{h}'},
"tt_elmu":       {'color': rgb(254,224,210),   'legend': 't#bar{t}#rightarrow e#mu'},
"tt_taultauh":   {'color': rgb(253,141,60),    'legend': 't#bar{t}#rightarrow #tau_{l}#tau_{h}'},
"tt_ltaul":      {'color': rgb(252,146,114),   'legend': 't#bar{t}#rightarrow l#tau_{l}'},
"tt_taueltaumu": {'color': rgb(252,146,114),   'legend': 't#bar{t}#rightarrow #tau){e}#tau_{#mu}'},
"tt_ljw":        {'color': rgb(240,59,32),     'legend': 't#bar{t}#rightarrow ljw'},
"tt_ljb":        {'color': rgb(240,59,32),     'legend': 't#bar{t}#rightarrow ljb'},
"tt_ljo":        {'color': rgb(240,59,32),     'legend': 't#bar{t}#rightarrow ljo'},
"tt_lj":         {'color': rgb(240,59,32),     'legend': 't#bar{t}#rightarrow lj'},
"tt_taulj":      {'color': rgb(210,29,32),     'legend': 't#bar{t}#rightarrow #tau_{l}j'},
"tt_other":      {'color': rgb(189,0,38),      'legend': 't#bar{t}#rightarrow other'},

"qcd": {'color': kViolet, 'legend': 'qcd'},
}

nick_colour = {nick: nick_info[nick]['color'] for nick in nick_info}

nick_order = {
"qcd": 0,

"tt_mutau3ch":  -1,
"tt_eltau3ch":  -1,
"tt_mutau": -2,
"tt_eltau": -2,
"tt_elmu":  -3,
"tt_ltaul": -4,
"tt_taueltaumu": -4,
"tt_taultauh": -4,
"tt_ljw": -5,
"tt_ljb": -6,
"tt_ljo": -7,
"tt_lj":  -8,
"tt_taulj": -9,
"tt_other": -11,

"wjets": -20,
"wjets_tauh": -21,
"wjets_taul": -22,

"dy_tautau": -25,
"dy_other":  -26,

"s_top_eltau": -30,
"s_top_mutau": -31,
"s_top_elmu":  -32,
"s_top_lj":    -33,
"s_top_other": -34,

"dibosons": -40,
}

def stack_n_legend(used_histos, shift=0., exp_legend=False):
    '''stack_n_legend(used_histos)

    used_histos = [(histo, nick of process, channel), ...]
    first check if there are different channels -- then add the channel name to the legend
    also the line style changes for different channels of the same process for histograms
    (not sure if it is a good hack)
    '''

    channels = [c for _, _, c in used_histos]
    homogeneous_channels = all(c == channels[0] for c in channels)

    # build Stach and legend
    hs = THStack("mc_stack", "mc_stack")
    #leg = TLegend(0.7 - (0.15 if not homogeneous_channels else 0.), 0.4, 0.89, 0.89)
    if exp_legend:
        leg = TLegend(0.8 - shift, 0.4, 1. - shift, 0.89)
    else:
        leg = TLegend(0.7 - shift, 0.4, 0.89 - shift, 0.89)

    process_counter = {} # for distinguishing processes in different channels

    for histo, nick, channel in sorted(used_histos, key=lambda h_n: nick_order.get(h_n[1], 1)):
        proc_ocurance = process_counter.setdefault(nick, 1)
        process_counter[nick] += 1

        col = nick_colour[nick]
        histo.SetFillColor( col );
        #histo.SetLineColor( col ); # it's needed for shapes
        histo.SetMarkerStyle(20);
        histo.SetLineStyle(proc_ocurance);
        histo.SetMarkerColor(col);
        #used_histos.append(histo) # hopefully root wont screw this up
        hs.Add(histo, "HIST")

    # to have legend in the same order
    for histo, nick, channel in sorted(used_histos, key=lambda h_n: -nick_order.get(h_n[1], 1)):
        if homogeneous_channels:
            leg.AddEntry(histo, nick_info[nick]['legend'], "F")
        else:
            leg.AddEntry(histo, "%s %s" % (nick, channel), "F")

    return hs, leg


