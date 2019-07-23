import argparse
import logging
from os.path import isfile, basename
import ctypes
from sys import argv


parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "sumup loops over TTree",
    epilog = """Example:\npython sumup_ttree_events.py --output foobar.root gstore_outdirs/v28/SingleMuon/Ntupler_v28_Data13TeV_SingleMuon2016*/1*/*/*root"""
    )

parser.add_argument('--ttree',    type=str, default="ntupler/reduced_ttree", help='path to the TTree in the file')

parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")
parser.add_argument("--output", type=str, default="output.root", help="filename for output")
parser.add_argument("--n-events", type=int, default=-1, help="maximum events per file (-1 means no limit, default)")

parser.add_argument("--save-weight", action='store_true', help="save the weight counters in the output file")
parser.add_argument("--per-weight",  action='store_true', help="normalize by event weight of datasets")

parser.add_argument('input_files', nargs='+', help="""the files to sum up, passed by shell, as:

/gstore/t3cms/store/user/otoldaie/v19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v19_MC2016_Summer16_TTJets_powheg/180226_022336/0000/MC2016_Summer16_TTJets_powheg_*.root""")



'''
out_histo  = TH1D("os_tauall", "", 44, -2, 20)
tau_b      = TH1D("os_taub", "", 44, -2, 20)
tau_w      = TH1D("os_tauw", "", 44, -2, 20)
tau_w_c    = TH1D("os_tauw_c", "", 44, -2, 20)
tau_w_notc = TH1D("os_tauw_notc", "", 44, -2, 20)

ss_tau_all    = TH1D("ss_tauall", "", 44, -2, 20)
ss_tau_b      = TH1D("ss_taub", "", 44, -2, 20)
ss_tau_w      = TH1D("ss_tauw", "", 44, -2, 20)
ss_tau_w_c    = TH1D("ss_tauw_c", "", 44, -2, 20)
ss_tau_w_notc = TH1D("ss_tauw_notc", "", 44, -2, 20)
'''


event_process = None
output_histos = None


def process_files(input_files, ttree_name, output_filename, max_events, per_weight=False, save_weight=False):
    global event_process, output_histos
    weight_counter = None

    from ROOT import TFile

    for filename in input_files:
        if not isfile(filename):
            logging.info("missing: " + filename)
            continue

        logging.debug(filename)
        logging.debug(ttree_name)

        tfile = TFile(filename)
        ttree = tfile.Get(ttree_name)

        for iev, ev in enumerate(ttree):
            if max_events > 0 and iev > max_events: break
            #event_process(ev, output_histos)
            event_process(ev)

        # calculate weights to scale MC
        if per_weight or save_weight:
            wcounter = tfile.Get('ntupler/weight_counter')
            if not wcounter:
                # try top level
                wcounter = tfile.Get('weight_counter')
            assert bool(wcounter)
            if not weight_counter:
                weight_counter = wcounter
                weight_counter.SetDirectory(0)
                weight_counter.SetName("sumup_weight_counter")
            else:
                weight_counter.Add(wcounter)

        tfile.Close()

    if per_weight:
        #weight_counter = tfile.Get('ntupler/weight_counter')
        for out_histo in [n_tau_cands, tau_pt, tau_sv_sign]:
            out_histo.Scale(1./weight_counter.GetBinContent(2))


    fout = TFile(output_filename, "RECREATE")
    fout.Write()

    fout.cd()

    #n_tau_cands.Write()
    #tau_pt.Write()
    #tau_sv_sign.Write()
    for histo in output_histos.values():
        histo.Write()

    # save weight counters etc in the top directory in the file
    if save_weight:
        fout.cd()
        weight_counter.SetName('weight_counter')
        weight_counter.Write()

    fout.Write()
    fout.Close()


'''
for extension and command line support:
set up the event processing and run main with command line parameters
'''

def main(args):
    #def process_files(input_files, ttree_name, output_filename, per_weight=False, save_weight=False):
    process_files(args.input_files, args.ttree, args.output, args.n_events, args.per_weight, args.save_weight)

if __name__ == '__main__':
    #global event_process, output_histos

    args = parser.parse_args(command_line_args)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.info("import ROOT")

    import ROOT
    from ROOT import TH1D, TFile, TTree, gROOT, TLorentzVector
    from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
    from ROOT import TLegend
    #from ROOT import THStack

    logging.info("done")

    from plotting_root import rgb

    gROOT.SetBatch()

    # The histograms to collect
    # 
    name_prefix = "dy_"

    std_output = {
    'n_tau_cands':    TH1D(name_prefix + "n_tau_cands", "", 10,  0, 10),
    'tau_pt':         TH1D(name_prefix + "tau_pt",      "", 100, 0, 200),
    'n_tau_sv_cands': TH1D(name_prefix + "n_tau_sv_cands", "", 10,  0, 10),
    'tau_sv_pt':      TH1D(name_prefix + "tau_sv_pt",      "", 100, 0, 200),
    'tau_sv_sign':    TH1D(name_prefix + "tau_sv_sign",    "", 21, -1, 20),
    }

    output_histos = std_output


    def std_process(ev):
        # pass conditions
        if not (abs(ev.leps_ID) == 13 and ev.HLT_mu and ev.lep_matched_HLT[0]): return None
        pass_lep = ev.lep_p4[0].pt() > 30 and abs(ev.lep_p4[0].eta()) < 2.4
        if not pass_lep: return None

        if len(ev.tau_p4) < 1: return None
        # sv_ind = -1 if there was no valid SV refited
        # otherwise it is an index into vectors on SVs
        taus_main = [(IDlev, p4, PDGid, sv_ind) for IDlev, p4, PDGid, sv_ind in zip(ev.tau_IDlev, ev.tau_p4, ev.tau_id, ev.tau_refited_index) if IDlev > 0 and p4.pt() > 30. and abs(p4.eta()) < 2.4]
        taus_main .sort(key=lambda it: (it[0], it[1].pt()), reverse=True)

        #pass_tau = ev.tau_p4[0].pt() > 30 and abs(ev.tau_p4[0].eta()) < 2.4 and len(ev.tau_IDlev) > 0 and ev.tau_IDlev[0] > 2
        pass_tau = len(taus_main) > 0 and taus_main[0][0] > 2
        if not pass_tau: return None
        # pass OS
        if not (taus_main[0][2] * ev.lep_id[0] < 0): return None

        # check there no b-jets among not-tau jets
        main_tau = taus_main[0]
        mtau = main_tau[1]
        main_tau_tlor = TLorentzVector(mtau.X(), mtau.Y(), mtau.Z(), mtau.T())

        n_bjets = 0
        for i in xrange(ev.jet_p4.size()):
            pfid = ev.jet_PFID[i]
            p4 = ev.jet_initial_p4[i]
            if pfid < 1 or abs(p4.eta()) > 2.5: continue # Loose PFID and eta
            if ev.jet_matching_lep[i]: continue

            tj_p4 = TLorentzVector(p4.X(), p4.Y(), p4.Z(), p4.T())
            # if matches to tau
            if tj_p4.DeltaR(main_tau_tlor) < 0.4: continue

            # else count the b-tagged jets
            jet_b_discr = ev.jet_b_discr[i]
            if jet_b_discr > 0.8484: # medium b WP
                n_bjets += 1

        if n_bjets > 0: return None


        # save stuff
        output_histos['n_tau_cands'].Fill(len(taus_main))
        output_histos['tau_pt'].Fill(mtau.pt())

        # if the main tau has valid SV
        #
        if main_tau[3] < 0: return None
        output_histos['n_tau_sv_cands'].Fill(len(taus_main))
        output_histos['tau_sv_pt'].Fill(mtau.pt())
        # TODO: WARNING! THE PT 0 PROBABLY DOES NOT CORRESPOND TO SV 0
        output_histos['tau_sv_sign'].Fill(ev.tau_SV_geom_flightLenSign[main_tau[3]])


    event_process = std_process



    main(args)

