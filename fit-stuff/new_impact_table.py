import argparse
import logging
import json
import logging
from os.path import isfile


logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "print out the impacts in table format",
    epilog = "Example:\npython impact_table.py latest_el_impacts.json\npython impact_table.py latest_both_impacts.json --merge-pdfs\npython new_impact_table.py --merge-pdfs latest_el_impacts.json latest_mu_impacts.json latest_both_impacts.json"
    )

parser.add_argument("--merge-pdfs", action='store_true', help="add PDFs in quadrature")
parser.add_argument("--quote-dileps", action='store_true', help="quote dilep systematics with correlations and extrapolation")
parser.add_argument("--precision",  type=int, default=1, help="number of digits after coma in output")
parser.add_argument("impacts_file", nargs='+', help="json file with impacts")

args = parser.parse_args()

assert all(isfile(f) for f in args.impacts_file)

impact_name_exp = {
#"lumi_13TeV"  : ("integrated luminosity",       2.5,  0),
"tauID_eff"   : ("$\\tauh$ jet identification", 0.0,  0),
#"tau_fakes"  : ("$\\tauh$ misidentification background", ),
"tau_fakes"   : ("\\ttbar background",          0.2,  0),
"dy_norm"     : ("DY background",               0.9,  1),
"tt_norm"     : ("tt background",               0.2,  0),
"stop_norm"   : ("tW background",               1.1,  1),
"dibos_norm"  : ("dibosons background",         0.2,  1),
"wjets_norm"  : ("W+jets background",           0.2,  0),
"qcd_norm"    : ("multijet background",     '$<0.1$', 0),

"LEP"         : ("lepton ID efficiency",        0, 0),
"LEPel"       : ("electron ID efficiency",      0, 0),
"LEPmu"       : ("muon ID efficiency",          0, 0),
"LEPelID"     : ("electron ID efficiency",      0, 0),
"LEPmuID"     : ("muon ID efficiency",          0, 0),
"LEPelTRG"    : ("electron trigger efficiency", 0, 0),
"LEPmuTRG"    : ("muon trigger efficiency",     0, 0),

"bSF"         : ("\\cPqb tagging efficiency",   0.4,  1),
"JES"         : ("jet energy scale",        0.4, 1), #"JES",
"JER"         : ("jet energy resolution",   0.4, 1), #"JER",
"TauES"       : ("$\\tauh$ energy scale",   0.0, 0), #"TES",
"TES"         : ("$\\tauh$ energy scale",   0.0, 0), #"TES",
"PU"          : ("pileup (PU) uncertainty", 0.1, 1), #"pileup uncertainty",
}

impact_name_th = {
"TOPPT"         : ("top quark $p_{T}$ modeling", 0.5, 1),
"Frag"          : ("b fragmentation",            0.7, 1),
"Peterson"      : ("b frag., ), Peterson variation", -0.999, 1),
"SemilepBR"     : ("semi-leptonic b-hadron branching fraction", 0.1, 1), # "semi-leptonic b-hadron BRs",

"Mr"            : ("ME Renormalization scale",    -0.999, 1),
"Mf"            : ("ME Factorization scale",      -0.999, 1),
"Mfr"           : ("ME combined scale variation", -0.999, 1),
"ME"            : ("\\ttbar ME scale",             0.2  , 1),

"FSR"           : ("\\ttbar FSR scale", 0.8, 1), #"final   state radiation",
"ISR"           : ("\\ttbar ISR scale", 0.4, 1), #"initial state radiation",
"HDAMP"         : ("ME-PS matching",    0.2, 1), # "ME-PS matching (hdamp)",
"TuneCUETP8M2T4": ("underlying event",  0.3, 1), # "underlying event (Pythia tune)",
}

common_uncertainties = [
("integrated luminosity", 2.5, 2.5, 2.5, 2.5, 1),
("statistics", 1.4, 1.1, 0.9, 0.2, 0),
]

extrapolation_uncertainties = [
("\\ttbar ME scales", 0.3, 0.4, 0.3, 0.2, 0),
("PDF",               1.2, 1.4, 1.3, 0.2, 0),
("underlying event",  0.3, 0.2, 0.2, 0.3, 0),
("ME-PS matching",    0.8, 1.0, 0.9, 0.2, 0),
("\\ttbar ISR scale", 0.5, 0.3, 0.3, 0.4, 0),
("\\ttbar FSR scale", 1.9, 2.0, 1.9, 0.8, 0),
]

dilepton_only_impacts = [
("trigger efficiency"                   , 0.3, 0),
#("lepton identification and isolation"  , 2.0, 1),
("electron efficiencies"            , 1.5, 1),
("muon efficiencies"                , 1.3, 1),
("muon momentum scale"              , 0.1, 1),
("electron momentum scale"          , 0.1, 1),
("tW ME scale"                      , 0.2, 1),
("DY ME scale"                      , 0.1, 1),
("tW ISR scale"                     , 0.1, 1),
("tW FSR scale"                     , 0.1, 1),
("colour reconnection"              , 0.3, 1),
# also some leptau systs are calculated on the fly -- attach them!
("PDF dilep"    , 1.1, 1),
("MC statistics", 1.1, 0),
]


impact_name_dictionary = {}

impact_name_dictionary.update(impact_name_exp)
impact_name_dictionary.update(impact_name_th)

#print impact_name_dictionary

impacts_exp = {}
impacts_th  = {}
impacts_mc_stat  = []

n_files = 0
for fname in args.impacts_file:
    f = open(fname)

    impacts_js = (json.loads(f.read()))['params']

    #print fname

    for param in impacts_js:
        if param['name'] not in impact_name_exp.keys(): continue
        if 'prop_bin' in  param['name']: continue
        #print param['name']
        impacts_exp.setdefault(param['name'], []).append(param['impact_r'])

    # mc stat
    impacts_mc_stat.append(0.)
    for param in impacts_js:
        if 'prop_bin' not in  param['name']: continue
        impacts_mc_stat[-1] += param['impact_r']**2

    for param in impacts_js:
        name = param['name']
        if name in impact_name_exp.keys(): continue
        impacts_th.setdefault(name, []).append(param['impact_r'])

    n_files += 1
    #impacts_exp = {param['name']: param['impact_r'] for param in impacts_js if param['name'] in impact_name_exp.keys()}
    #impacts_th  = {param['name']: param['impact_r'] for param in impacts_js if param['name'] not in impact_name_exp.keys()}

if 'Frag' in impacts_th and 'Peterson' in impacts_th:
    # these should be added as 'Frag'
    #print "adding frags"
    frags    = impacts_th.pop('Frag')
    peterson = impacts_th.pop('Peterson')
    impacts_th['Frag'] = [(p**2 + f **2)**0.5 for f, p in zip(frags, peterson)]

if 'Mfr' in impacts_th and 'Mf' in impacts_th and 'Mr' in impacts_th:
    # these should be added as 'Frag'
    #print "adding frags"
    m_fr = impacts_th.pop('Mfr')
    m_f  = impacts_th.pop('Mf')
    m_r  = impacts_th.pop('Mr')
    impacts_th['ME'] = [(fr**2 + f**2 + r**2)**0.5 for fr, r, f in zip(m_fr, m_f, m_r)]


#table_row_format = "\\lstinline!%25s!  & %4.2f \\\\"
#table_row_format = "\\lstinline!%25s!  & %s \\\\"
table_row_format = "%25s  & %s \\\\"

imp_form = '%4.{}f'.format(args.precision)

#print impacts_th
#print impacts_exp

if not args.merge_pdfs:
    for name, imp in sorted(impacts.items(), key=lambda k: k[1], reverse=True):
        print table_row_format % (name, imp*100)

else:
    #mcstat_quadrature_sum = [0. for _ in range(n_files)]
    pdf_quadrature_sum    = [0. for _ in range(n_files)]

    final_leptau_impacts = []

    for nick, impacts in sorted(impacts_exp.items(), key=lambda k: k[1][-1], reverse=True):
        if 'PDF' in nick or 'AlphaS' in nick:
            for i, imp in enumerate(impacts):
                pdf_quadrature_sum[i] += imp**2
            continue

        name = impact_name_dictionary.get(nick, (nick,))[0]
        leptau_impacts = [imp_form % (imp*100) for imp in impacts]
        if args.quote_dileps and nick in impact_name_dictionary:
            name, dilep_unc, correlation_with_dilep = impact_name_dictionary[nick]
            #dilep_sys = ' & %s ' % (' & '.join(str(i) for i in [dilep_unc, correlation_with_dilep]))
            leptau_impacts.extend(str(i) for i in [dilep_unc, correlation_with_dilep])

        final_leptau_impacts.append((name, leptau_impacts))

        #leptau_impacts = ' & '.join(imp_form % (imp*100) for imp in impacts)
        #if args.quote_dileps and nick in impact_name_dictionary:
        #    name, dilep_unc, correlation_with_dilep = impact_name_dictionary[nick]
        #    dilep_sys = ' & %s ' % (' & '.join(str(i) for i in [dilep_unc, correlation_with_dilep]))
        #    all_systematics = leptau_impacts + dilep_sys
        #    print table_row_format % (name, all_systematics)
        #else:
        #    name = impact_name_dictionary.get(nick, (nick,))[0]
        #    print table_row_format % (name, leptau_impacts)

    for nick, impacts in sorted(impacts_th.items(), key=lambda k: k[1][-1], reverse=True):
        if 'PDF' in nick or 'AlphaS' in nick:
            for i, imp in enumerate(impacts):
                pdf_quadrature_sum[i] += imp**2
            continue
        # the MC stat systematics
        if "_bin" in nick:
            continue # done separately
            #for i, imp in enumerate(impacts):
            #    #pdf_quadrature_sum[i] += imp**2
            #    mcstat_quadrature_sum[i] += imp**2

        name = impact_name_dictionary.get(nick, (nick,))[0]
        leptau_impacts = [imp_form % (imp*100) for imp in impacts]
        if args.quote_dileps and nick in impact_name_dictionary:
            name, dilep_unc, correlation_with_dilep = impact_name_dictionary[nick]
            #dilep_sys = ' & %s ' % (' & '.join(str(i) for i in [dilep_unc, correlation_with_dilep]))
            leptau_impacts.extend(str(i) for i in [dilep_unc, correlation_with_dilep])

        final_leptau_impacts.append((name, leptau_impacts))

        ##print name, impacts
        #print table_row_format % (name, ' & '.join(imp_form % (imp*100) for imp in impacts))

    for name, uncertainties in final_leptau_impacts:
        print table_row_format % (name, ' & '.join(uncertainties))

    pdf_impacts = [pdf**0.5 for pdf in pdf_quadrature_sum]
    print table_row_format % ('PDF', ' & '.join(imp_form % (imp*100) for imp in pdf_impacts))

    #print mcstat_quadrature_sum
    mcstat_impacts = [mcstat**0.5 for mcstat in impacts_mc_stat]
    print table_row_format % ('MC statistics', ' & '.join(imp_form % (imp*100) for imp in mcstat_impacts))

# print dilep-only uncertainties
# and acceptance uncertainty
if args.quote_dileps:

    print
    for name, uncertainty, correlation in dilepton_only_impacts:
        print ' & '.join(str(i) for i in [name] + ['$<0.1$']*3 + [uncertainty, correlation]) + ' \\\\'

    for static_unc_group in [common_uncertainties, extrapolation_uncertainties]:
        print
        for extr in static_unc_group:
            print ' & '.join(str(i) for i in extr) + ' \\\\'

