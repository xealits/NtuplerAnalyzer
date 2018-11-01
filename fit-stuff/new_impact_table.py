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
parser.add_argument("--precision",  type=int, default=1, help="number of digits after coma in output")
parser.add_argument("impacts_file", nargs='+', help="json file with impacts")

args = parser.parse_args()

assert all(isfile(f) for f in args.impacts_file)

impact_name_exp = {
"lumi_13TeV"  : "Luminosity",
"tauID_eff"   : "$\\tauh$ jet identification",
"tau_fakes"   : "$\\tauh$ misidentification background",
"dy_norm"     : "DY normalization",
"stop_norm"   : "Single Top normalization",
"dibos_norm"  : "Dibosons normalization",
"wjets_norm"  : "W+jets normalization",
"qcd_norm"    : "QCD normalization",
"bSF"         : "b-tagging efficiency",
"LEP"         : "lepton ID efficiency",
"LEPel"       : "electron ID efficiency",
"LEPmu"       : "muon ID efficiency",
"LEPelID"     : "electron ID efficiency",
"LEPmuID"     : "muon ID efficiency",
"LEPelTRG"    : "electron trigger efficiency",
"LEPmuTRG"    : "muon trigger efficiency",
"JES"         : "jet energy scale",
"JER"         : "jet energy resolution",
"TauES"       : "$\\tauh$ Energy Scale",
"TES"         : "$\\tauh$ Energy Scale",
"PU"          : "Pile-Up uncertainty",
"TOPPT"       : "top quark $p_{T}$",
}

impact_name_th = {
"Frag"        : "b fragmentation",
"Peterson"    : "b frag., Peterson variation",
"SemilepBR"   : "b decay tables, Semileptonic BR",
"Mr"          : "ME Renormalization scale",
"Mf"          : "ME Factorization scale",
"Mfr"         : "ME combined scale variation",
"FSR" : "Final State Radiation",
"ISR" : "Initial State Radiation",
"HDAMP" : "ME-PS matching (hdamp)",
"TuneCUETP8M2T4" : "Underlying event (Pythia tune)",
}

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
    for nick, impacts in sorted(impacts_exp.items(), key=lambda k: k[1][-1], reverse=True):
        if 'PDF' in nick or 'AlphaS' in nick:
            for i, imp in enumerate(impacts):
                pdf_quadrature_sum[i] += imp**2
            continue
        name = impact_name_dictionary.get(nick, nick)
        print table_row_format % (name, ' & '.join(imp_form % (imp*100) for imp in impacts))

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

        name = impact_name_dictionary.get(nick, nick)
        #print name, impacts
        print table_row_format % (name, ' & '.join(imp_form % (imp*100) for imp in impacts))

    pdf_impacts = [pdf**0.5 for pdf in pdf_quadrature_sum]
    print table_row_format % ('PDF', ' & '.join(imp_form % (imp*100) for imp in pdf_impacts))

    #print mcstat_quadrature_sum
    mcstat_impacts = [mcstat**0.5 for mcstat in impacts_mc_stat]
    print table_row_format % ('MC statistics', ' & '.join(imp_form % (imp*100) for imp in mcstat_impacts))


