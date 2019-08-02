import logging
from sys import argv

import sumup_ttree_events

def round2(fl):
    return round(fl, 2)

def round3(fl):
    return round(fl, 3)


output_histos_dict = {}

global prefix
prefix = 'tau'
def output_histos(histoname):
    return output_histos_dict[prefix + '_' + histoname]

"""
root [12] ttree_out->Draw("event_taus_pat_sv_leng[0]:event_taus[0].energy()>>h(20,30,150,20,0,1)", "", "col")
(Long64_t) 2752
root [13] 
root [13] 
root [13] ttree_out->Draw("event_taus_sv_leng[0]:event_taus[0].energy()>>h(20,30,150,20,(Long64_t) 2752")

gen_proc_id

Output:
process,energy,pat_sign,pat_len,geom_sign,geom_len

"""

def process_event_for_fake_taus(ev):
    global prefix
    # pass conditions

    process = ev.gen_proc_id

    # tt gen = lepton + jets
    prefix, is_lep_jets = 'lj',  (20 < process < 30)
    prefix, is_lep_tau  = 'tau', (process == 32 or process == 42)
    if not (is_lep_jets or is_lep_tau):
        return None

    # main with Vloose taus lep-tau selection, both electrons and muons
    #sel_el = ev.selection_stage > 13
    #sel_mu = 3 < ev.selection_stage < 10
    # medium:
    sel_el = ev.selection_stage > 15
    sel_mu = 5 < ev.selection_stage < 10
    if not (sel_el or sel_mu): return None

    # only SV refitted taus
    sv_sign = ev.event_taus_sv_sign[0]
    if sv_sign < 0:
        return None

    # SV > 2.
    if sv_sign < 2.: return None

    energy  = ev.event_taus[0].energy()
    sv_sign_pat = ev.event_taus_pat_sv_sign[0]
    sv_leng_pat = ev.event_taus_pat_sv_leng[0]
    sv_leng = ev.event_taus_sv_leng[0]

    output_histos('tau_sv_leng_VS_energy')    .Fill(sv_leng,     energy)
    output_histos('tau_sv_pat_leng_VS_energy').Fill(sv_leng_pat, energy)

    # printout
    print ('%d,' + ','.join(["%f"]*5)) % (process, energy, sv_sign_pat, sv_leng_pat, sv_sign, sv_leng)

if __name__ == '__main__':

    args = sumup_ttree_events.parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    from ROOT import TLorentzVector, TH1D, TH2D

    for pref in ['tau', 'ljw', 'ljb', 'ljo']:
        #output_histos_dict[pref + '_tau_sv_leng_VS_energy'] = TH2D(pref + "_tau_sv_leng_VS_energy", "", 20, -0.0, 1.,  15, 30, 150)
        output_histos_dict[pref + '_tau_sv_leng_VS_energy'] = TH2D(pref + "_tau_sv_leng_VS_energy", "", 21, -0.1, 2., 15, 0, 150)
        #output_histos('tau_sv_leng_VS_energy').Fill(main_tau_SV_leng, main_tau_p4.energy())
        output_histos_dict[pref + '_tau_sv_pat_leng_VS_energy'] = TH2D(pref + "_tau_sv_pat_leng_VS_energy", "", 20, -0.0, 0.1, 15, 30, 150)


    sumup_ttree_events.event_process = process_event_for_fake_taus
    sumup_ttree_events.output_histos = output_histos_dict

    print "process,energy,pat_sign,pat_len,geom_sign,geom_len"
    sumup_ttree_events.main(args)

