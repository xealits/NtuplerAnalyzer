import json
import logging


#f = open("latest_mu_impacts.json")
f = open("latest_el_impacts.json")

impacts_js = (json.loads(f.read()))['params']

impacts = {param['name']: param['impact_r'] for param in impacts_js}

for name, imp in sorted(impacts.items(), key=lambda k: k[1], reverse=True):
    print "\\lstinline!%10s!  & %4.1f \\\\" % (name, imp*100)

