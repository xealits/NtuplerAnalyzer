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


def process_event_for_fake_taus(ev):
    global prefix
    # pass conditions

    # tt gen = lepton + jets
    t_wid  = abs(ev.gen_t_w_decay_id)
    tb_wid = abs(ev.gen_tb_w_decay_id)
    prefix, is_lep_jets = 'lj',  (t_wid == 1 and (tb_wid == 11 or tb_wid == 13)) or (tb_wid == 1 and (t_wid == 11 or t_wid == 13))
    prefix, is_lep_tau  = 'tau', (t_wid > 15*15 and (tb_wid == 11 or tb_wid == 13)) or (tb_wid > 15*15 and (t_wid == 11 or t_wid == 13))
    if not (is_lep_jets or is_lep_tau):
        return None

    # reco conditions: 1 good 3ch high SV tau
    #if not (abs(ev.leps_ID) == 13 and ev.HLT_mu and ev.lep_matched_HLT[0]): continue
    if not (abs(ev.leps_ID) == 13 or abs(ev.leps_ID) == 11): return None
    if len(ev.tau_p4) < 1: return None

    # medium tau and SV > 2.
    taus_main = zip(ev.tau_IDlev, ev.tau_p4, ev.tau_id, ev.tau_refited_index)
    taus_main .sort(key=lambda it: (it[0], it[1].pt()), reverse=True)
    # use tau_refited_index to access teh geometric SV
    if taus_main[0][0] < 3: return None

    # SV data
    main_tau = taus_main[0]
    main_tau_refit_index = main_tau[3]
    if main_tau_refit_index > -1:
        main_tau_SV_sign = ev.tau_SV_geom_flightLenSign [main_tau_refit_index]
        main_tau_SV_leng = ev.tau_SV_geom_flightLen     [main_tau_refit_index]

        # dalitz parameters
        tau_track_lengs = len(ev.tau_SV_fit_track_OS_p4) == len(ev.tau_SV_fit_track_SS1_p4) == len(ev.tau_SV_fit_track_SS2_p4)
        if not tau_track_lengs:
            print "the problem in tau track lists tau_track_lengs"
            return None

        out_of_range  = main_tau_refit_index > len(ev.tau_SV_fit_track_OS_p4)
        out_of_range |= main_tau_refit_index > len(ev.tau_SV_fit_track_SS1_p4)
        out_of_range |= main_tau_refit_index > len(ev.tau_SV_fit_track_SS2_p4)
        if out_of_range:
            print "the problem in tau track lists main_tau_refit_index"
            return None

        dalitz_m1 = (ev.tau_SV_fit_track_OS_p4[main_tau_refit_index] + ev.tau_SV_fit_track_SS1_p4[main_tau_refit_index]).mass()
        dalitz_m2 = (ev.tau_SV_fit_track_OS_p4[main_tau_refit_index] + ev.tau_SV_fit_track_SS2_p4[main_tau_refit_index]).mass()

        tau_track_mass = (ev.tau_SV_fit_track_OS_p4[main_tau_refit_index] + ev.tau_SV_fit_track_SS1_p4[main_tau_refit_index] + ev.tau_SV_fit_track_SS2_p4[main_tau_refit_index]).mass()
    else:
        main_tau_SV_sign = -11.
        main_tau_SV_leng = -11.
        dalitz_m1 = -11.
        dalitz_m2 = -11.
        tau_track_mass = -11.

    # SV > 2.
    if main_tau_SV_sign < 2.: return None

    # basic data on the main tau
    main_tau_p4 = main_tau[1]
    if main_tau_p4.pt() < 21.: return None
    #TLorentzVector::TLorentzVector(double x, double y, double z, double t) =>
    tau_p4_tlor = TLorentzVector(main_tau_p4.X(), main_tau_p4.Y(), main_tau_p4.Z(), main_tau_p4.T())

    # for fake tous figure out is the origin b or W
    # print out the final products around tau, sorted by dR
    print 'tau', round2(main_tau_p4.eta()), round2(main_tau_p4.phi()), round3(main_tau_SV_sign), round3(main_tau_SV_leng)
    fake_taus_t_w  = []
    fake_taus_t_b  = []
    fake_taus_tb_w = []
    fake_taus_tb_b = []
    final_products = [
        ('tw1',   fake_taus_t_w,  ev.gen_t_w1_final_pdgIds,  ev.gen_t_w1_final_statuses,  ev.gen_t_w1_final_p4s),
        #('tw2',  fake_taus_t_w,  ev.gen_t_w2_final_pdgIds,  ev.gen_t_w2_final_statuses,  ev.gen_t_w2_final_p4s), # repetition
        ('tb',    fake_taus_t_b,  ev.gen_t_b_final_pdgIds,   ev.gen_t_b_final_statuses,   ev.gen_t_b_final_p4s),
        ('tbw1',  fake_taus_tb_w, ev.gen_tb_w1_final_pdgIds, ev.gen_tb_w1_final_statuses, ev.gen_tb_w1_final_p4s),
        #('tbw2', fake_taus_tb_w, ev.gen_tb_w2_final_pdgIds, ev.gen_tb_w2_final_statuses, ev.gen_tb_w2_final_p4s),
        ('tbb',   fake_taus_tb_b, ev.gen_tb_b_final_pdgIds,  ev.gen_tb_b_final_statuses,  ev.gen_tb_b_final_p4s),
        ]

    # process the products per origin
    if is_lep_jets:
      for origin, fakes_close_in_dR, pdgIDs, statuses, p4s in final_products:
        #tau_p4_tlor = TLorentzVector(main_tau_p4)
        products_per_dR = []
        for prod in zip(pdgIDs, statuses, p4s):
            # skip neutrinos (did not remove them in b decays)
            # actually should have kept them everywhere
            if abs(prod[0]) in (12, 14, 16): continue

            #dR = prod[2].DeltaR(main_tau_p4)
            #prod_p4_tlor = TLorentzVector(prod[2])
            prod_p4_tlor = TLorentzVector(prod[2].X(), prod[2].Y(), prod[2].Z(), prod[2].T())
            dR = prod_p4_tlor.DeltaR(tau_p4_tlor)

            if dR > 0.4: continue
            else:
                a_close_fake = (dR, prod[0], prod[1], prod[2])
                products_per_dR   .append(a_close_fake)
                fakes_close_in_dR .append(a_close_fake)

        if not products_per_dR: continue

        # sort by dR
        products_per_dR   .sort(key=lambda p: p[0])
        fakes_close_in_dR .sort(key=lambda p: p[0])

        # print in 1 line
        #products_per_dR.append((round3(dR), prod[0], prod[1], round(prod[2].eta(),2), round(prod[2].phi(), 2)))
        print origin, [(round3(prod[0]), prod[1], prod[2], round2(prod[3].eta()), round2(prod[3].eta())) for prod in products_per_dR]

    if is_lep_jets:
        #is_w = (ev.gen_t_w1_final_pdgIds.size() + ev.gen_tb_w1_final_pdgIds.size()) > 0
        #is_b = (ev.gen_t_b_final_pdgIds.size()  + ev.gen_tb_b_final_pdgIds.size() ) > 0
        is_w = (len(fake_taus_t_w) + len(fake_taus_tb_w)) > 0
        is_b = (len(fake_taus_t_b) + len(fake_taus_tb_b)) > 0
        if is_w and not (is_b):
            prefix = 'ljw'
        elif is_b and not (is_w):
            prefix = 'ljb'
        else:
            prefix = 'ljo'

    for origin, products_per_dR in [('tw', fake_taus_t_w), ('tb', fake_taus_t_b), ('tbw', fake_taus_tb_w), ('tbb', fake_taus_tb_b)]:
        # save info about products
        if origin in ('tb', 'tbb'):
            for dR, pdgID, _, _ in products_per_dR:
                output_histos('tau_fake_tb').Fill(abs(pdgID))
                if dR < 0.1:
                    output_histos('tau_fake_tb_close').Fill(abs(pdgID))
        elif 'w' in origin:
            for dR, pdgID, _, _ in products_per_dR:
                output_histos('tau_fake_tw').Fill(abs(pdgID))
                if dR < 0.1:
                    output_histos('tau_fake_tw_close').Fill(abs(pdgID))

    output_histos('tau_IDlev').Fill(main_tau[0])
    output_histos('tau_pt').Fill(main_tau_p4.pt())
    output_histos('tau_eta').Fill(main_tau_p4.eta())

    output_histos('tau_mass')       .Fill(main_tau_p4.mass())
    output_histos('tau_track_mass') .Fill(tau_track_mass)

    output_histos('tau_R1').Fill(ev.tau_leadChargedHadrCand_pt[0] / main_tau_p4.pt())
    output_histos('tau_R2').Fill(ev.tau_leadCand_pt[0]            / main_tau_p4.pt())

    # SV related params
    output_histos('tau_sv_sign').Fill(main_tau_SV_sign)
    output_histos('tau_sv_leng').Fill(main_tau_SV_leng)

    output_histos('tau_sv_leng_VS_energy').Fill(main_tau_SV_leng, main_tau_p4.energy())

    output_histos('tau_m1').Fill(dalitz_m1)
    output_histos('tau_m2').Fill(dalitz_m2)

    output_histos('tau_m1_m2').Fill(dalitz_m1, dalitz_m2)

if __name__ == '__main__':

    args = sumup_ttree_events.parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    from ROOT import TLorentzVector, TH1D, TH2D

    for pref in ['tau', 'ljw', 'ljb', 'ljo']:
        output_histos_dict[pref + '_tau_IDlev'] =   TH1D(pref + "_tau_IDlev", "", 10,  0, 10)
        output_histos_dict[pref + '_tau_pt'] =      TH1D(pref + "_tau_pt", "", 20,  0, 200)
        output_histos_dict[pref + '_tau_eta'] =     TH1D(pref + "_tau_eta", "", 20,  -2.5, 2.5)

        output_histos_dict[pref + '_tau_R1'] =     TH1D(pref + "_tau_R1", "", 61,  0., 1.22)
        output_histos_dict[pref + '_tau_R2'] =     TH1D(pref + "_tau_R2", "", 61,  0., 1.22)

        output_histos_dict[pref + '_tau_mass']       = TH1D(pref + "_tau_mass", "", 51,  0., 2.55)
        output_histos_dict[pref + '_tau_track_mass'] = TH1D(pref + "_tau_track_mass", "", 51,  0., 2.55)

        output_histos_dict[pref + '_tau_m1'] =     TH1D(pref + "_tau_m1", "", 51,  0., 1.53)
        output_histos_dict[pref + '_tau_m2'] =     TH1D(pref + "_tau_m2", "", 51,  0., 1.53)
        output_histos_dict[pref + '_tau_m1_m2'] =  TH2D(pref + "_tau_m1_m2", "", 51,  0., 1.53, 51,  0., 1.53)

        output_histos_dict[pref + '_tau_sv_sign'] = TH1D(pref + "_tau_sv_sign", "", 21,  -1., 20.)
        output_histos_dict[pref + '_tau_sv_leng'] = TH1D(pref + "_tau_sv_leng", "", 21,  -0.1, 2.)
        output_histos_dict[pref + '_tau_sv_leng_VS_energy'] = TH2D(pref + "_tau_sv_leng_VS_energy", "", 21, -0.1, 2., 15, 0, 150)
        output_histos_dict[pref + '_tau_fake_tb'] = TH1D(pref + "_tau_fake_tb", "", 1000,  0., 1000.)
        output_histos_dict[pref + '_tau_fake_tw'] = TH1D(pref + "_tau_fake_tw", "", 1000,  0., 1000.)
        output_histos_dict[pref + '_tau_fake_tb_close'] = TH1D(pref + "_tau_fake_tb_close", "", 1000,  0., 1000.)
        output_histos_dict[pref + '_tau_fake_tw_close'] = TH1D(pref + "_tau_fake_tw_close", "", 1000,  0., 1000.)


    sumup_ttree_events.event_process = process_event_for_fake_taus
    sumup_ttree_events.output_histos = output_histos_dict

    sumup_ttree_events.main(args)

