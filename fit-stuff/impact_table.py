import argparse
import logging
import json
import logging
from os.path import isfile


logging.basicConfig(level=logging.DEBUG)


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "print out the impacts in table format",
    epilog = "Example:\npython impact_table.py latest_el_impacts.json"
    )

parser.add_argument("impacts_file", help="json file with impacts")
parser.add_argument("--merge-pdfs", action='store_true', help="add PDFs in quadrature")

args = parser.parse_args()

assert isfile(args.impacts_file)

f = open(args.impacts_file)

impacts_js = (json.loads(f.read()))['params']

impacts = {param['name']: param['impact_r'] for param in impacts_js}

table_row_format = "\\lstinline!%10s!  & %4.2f \\\\"

if not args.merge_pdfs:
    for name, imp in sorted(impacts.items(), key=lambda k: k[1], reverse=True):
        print table_row_format % (name, imp*100)

else:
    pdf_quadrature_sum = 0.
    for name, imp in sorted(impacts.items(), key=lambda k: k[1], reverse=True):
        if 'PDF' in name:
            pdf_quadrature_sum += imp**2
            continue
        print table_row_format % (name, imp*100)
    pdf_impact = pdf_quadrature_sum**0.5
    print table_row_format % ('PDF', pdf_impact*100)

