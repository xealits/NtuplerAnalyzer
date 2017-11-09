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

nick_colour = {
"data": kWhite,
"dy": kGray,
"dy_other": kGray,
"dy_tautau": kGray+2,
#"wjets": kRed+1,
"wjets": kCyan+2,
"dibosons": kCyan,
"singletop": kAzure,
"s_top_eltau": kAzure,
"s_top_mutau": kAzure,
"s_top_lj": kAzure+1,
"s_top_other": kAzure+2,

#"tt_taultauh" : kCyan-3,
#"tt_jj": kGreen+4,
#"tt_lj": kGreen+3,
#"tt_em": kYellow-7,
#"tt_{em}": kYellow-7,
#"tt_ee": kAzure-9,
#"tt_mm": kGreen-9,
#"tt_eltau": kOrange+2,
#"tt_mutau": kOrange+1,
#"tt_{l\\tau-l}": kOrange+3,
#"tt_other": kCyan-5,

"tt_jj": kGreen+4,
"tt_em": kYellow-7,
"tt_{em}": kYellow-7,
"tt_ee": kAzure-9,
"tt_mm": kGreen-9,
"tt_{l\\tau-l}": kOrange+3,
#"tt_mutau": kOrange+0,
#"tt_eltau": kOrange+1,
#"tt_taultauh" : kOrange+7,
#"tt_lj": kRed+1, #kOrange+4,
#"tt_other": kOrange+5,

"tt_mutau":    rgb(255, 255, 178),
"tt_eltau":    rgb(254,204,92),
"tt_taultauh": rgb(253,141,60),
"tt_lj":       rgb(240,59,32),
"tt_other":    rgb(189,0,38),

"qcd": kViolet,
}

nick_order = {
"qcd": 0,

"tt_mutau": -1,
"tt_eltau": -2,
"tt_taultauh": -3,
"tt_lj": -4,
"tt_other": -5,

"wjets": -10,
"dy_tautau": -11,
"dy_other": -12,

"s_top_eltau": -20,
"s_top_mutau": -21,
"s_top_lj": -22,
"s_top_other": -23,
"dibosons": -30,
}

def stack_n_legend(used_histos):
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
    leg = TLegend(0.75 - (0.15 if not homogeneous_channels else 0.), 0.5, 0.89, 0.89)

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
            leg.AddEntry(histo, nick, "F")
        else:
            leg.AddEntry(histo, "%s %s" % (nick, channel), "F")

    return hs, leg


