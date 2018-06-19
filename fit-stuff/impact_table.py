import argparse
import logging
import json
import logging
from os.path import isfile


logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "print out the impacts in table format",
    epilog = "Example:\npython impact_table.py latest_el_impacts.json\npython impact_table.py latest_both_impacts.json --merge-pdfs"
    )

parser.add_argument("impacts_file", help="json file with impacts")
parser.add_argument("--merge-pdfs", action='store_true', help="add PDFs in quadrature")

args = parser.parse_args()

assert isfile(args.impacts_file)

impact_name_exp = {
"lumi_13TeV"  : "Luminosity",
"tauID_eff"   : "Tau ID efficiency",
"tau_fakes"   : "Tau mis-ID probability",
"dy_norm"     : "DY normalization",
"stop_norm"   : "Single Top normalization",
"dibos_norm"  : "Dibosons normalization",
"wjets_norm"  : "W+jets normalization",
"qcd_norm"    : "QCD normalization",
"bSF"         : "b-tagging efficiency",
"LEP"         : "lepton ID efficiency",
"JES"         : "Jet Energy Scale",
"JER"         : "Jet Energy Resolution",
"TauES"       : "Tau Energy Scale",
"PU"          : "Pile-Up uncertainty",
"TOPPT"       : "Top $p_{T}$",
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

f = open(args.impacts_file)

impacts_js = (json.loads(f.read()))['params']

impacts_exp = {param['name']: param['impact_r'] for param in impacts_js if param['name'] in impact_name_exp.keys()}
impacts_th  = {param['name']: param['impact_r'] for param in impacts_js if param['name'] not in impact_name_exp.keys()}

table_row_format = "\\lstinline!%25s!  & %4.2f \\\\"

if not args.merge_pdfs:
    for name, imp in sorted(impacts.items(), key=lambda k: k[1], reverse=True):
        print table_row_format % (name, imp*100)

else:
    pdf_quadrature_sum = 0.
    for nick, imp in sorted(impacts_exp.items(), key=lambda k: k[1], reverse=True):
        if 'PDF' in nick or 'AlphaS' in nick:
            pdf_quadrature_sum += imp**2
            continue
        name = impact_name_dictionary.get(nick, nick)
        print table_row_format % (name, imp*100)

    for nick, imp in sorted(impacts_th.items(), key=lambda k: k[1], reverse=True):
        if 'PDF' in nick or 'AlphaS' in nick:
            pdf_quadrature_sum += imp**2
            continue
        name = impact_name_dictionary.get(nick, nick)
        print table_row_format % (name, imp*100)

    pdf_impact = pdf_quadrature_sum**0.5
    print table_row_format % ('PDF', pdf_impact*100)

