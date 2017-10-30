from ROOT import THStack, TLegend, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan

nick_colour = {
"data": kWhite,
"dy": kGray,
"dy_other": kGray,
"dy_tautau": kGray+2,
"wjets": kRed+1,
"dibosons": kCyan,
"singletop": kAzure,
"s_top_eltau": kAzure,
"s_top_lj": kAzure+1,
"s_top_other": kAzure+2,

"tt_taultauh" : kCyan-3,
"tt_jj": kGreen+4,
"tt_lj": kGreen+3,
"tt_em": kYellow-7,
"tt_{em}": kYellow-7,
"tt_ee": kAzure-9,
"tt_mm": kGreen-9,
"tt_eltau": kOrange+2,
"tt_mutau": kOrange+1,
"tt_{l\\tau-l}": kOrange+3,
"tt_other": kCyan-5,

"qcd": kViolet,
}


def stack_n_legend(used_histos):
    # build Stach and legend
    hs = THStack("mc_stack", "mc_stack")
    leg = TLegend(0.7, 0.7, 0.89, 0.89)
    for histo, nick in used_histos:
        col = nick_colour[nick]
        histo.SetFillColor( col );
        histo.SetLineColor( col );
        histo.SetMarkerStyle(20);
        histo.SetLineStyle(0);
        histo.SetMarkerColor(col);
        #used_histos.append(histo) # hopefully root wont screw this up
        hs.Add(histo, "HIST")
        leg.AddEntry(histo, nick, "F")
    return hs, leg

