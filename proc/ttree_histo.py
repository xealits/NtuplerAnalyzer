import argparse
import logging
from os.path import isfile, basename




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "process the ttree in the file",
    epilog = """Example:\npython ttree_histo.py --input MC2016_Summer16_TTJets_powheg.root:ttree_out --output test1.root"""
    )

parser.add_argument("--input",  type=str, default="input.root:ttree",  help="filename and ttree name of the input")
parser.add_argument("--output", type=str, default="output.root", help="filename for output")
parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")


args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)


logging.info("import ROOT")

import ROOT
from ROOT import TH1D, TFile, TTree, gROOT
#from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
#from ROOT import TLegend
#from ROOT import THStack

logging.info("done")

from plotting_root import rgb

gROOT.SetBatch()



input_filename, ttree_path = args.input.split(':')

input_file = TFile(input_filename)
ttree      = input_file.Get(ttree_path)

#                                  n bins  first  last
histo_sig = TH1D("tau_pt_sig", "", 200,    0.,    250.)
histo_bck = TH1D("tau_pt_bck", "", 200,    0.,    250.)

for event in ttree:
    # you can find what the ttree contains by running in ROOT interpreter:
    # root -l MC2016_Summer16_TTJets_powheg.root
    # ttree_out->Print()

    # the events are stored with an index number
    # which encodes the selection requirements the event passes
    # the index is called selection_stage
    # selection_stage == 5 is our final selection
    # i.e. there is more signal and less bacground
    # 3 is for preselection, i.e. much more background than signal
    # there is also numbers 2 and 4, which are supporting selections and not needed now

    # processing only the final selection events
    # you can try selection 3 to check the difference
    if not event.selection_stage == 5:
        continue

    '''
    for the simulated events (called "Monte-Carlo" events, from the method they are generated)
    we know not only the response in the detector but also the internal decay chain
    from this information (called "generator information") we know exactly which event this was, signal or background

    it is encoded in gen_proc_id
    TT bar is main process for us, it is the signal process and also the main background
    but there are many subprocesses and many details about them
    so, gen_proc_id stores all these subprocesses for TTbar:

    genproc_tt_eltau3ch  = 42
    genproc_tt_eltau     = 41
    genproc_tt_mutau3ch  = 32
    genproc_tt_mutau     = 31
    
    genproc_tt_ljb       = 24
    genproc_tt_ljw       = 23
    genproc_tt_ljo       = 22
    genproc_tt_lj        = 21

    genproc_tt_taultauh  = 12
    genproc_tt_taulj     = 11
    
    genproc_tt_elmu       = 3
    genproc_tt_ltaul      = 2
    genproc_tt_taueltaumu = 1
    
    genproc_tt_other = 0

    -- our signal is tt_mutau and tt_mutau3ch (tau3ch is special type of tau)
       and main background is tt_lj with all the variations tt_ljb tt_ljw tt_ljo
       it is tt->lepton+jets process
       where one of W decays into lepton,
       but the other one instead of decaying to tau decays into quarks which produce jets in the detector
       so there is no tau
       but the detector mis-identifies one of the jets as product of decay of a tau
       those 4 variations lj, ljb, ljw, ljo show which exactly jet "faked" a tau
    '''

    if event.gen_proc_id in (31, 32):
        # it is signal event
        histo_sig.Fill(event.event_taus[0].pt())

    elif 20 < event.gen_proc_id < 25:
        # background
        histo_bck.Fill(event.event_taus[0].pt())

    # -- here I could test whether the tau exists in the event,
    #    but I know that at stage 5 events are required to have tau
    # otherwise, at stage 3, an event might not have tau and the program will break trying to access event_taus[0]

'''
Try comparing pt of tau -- maybe fake taus and true taus have different distribution?
Try plotting other distributions, like eta of tau:
event_taus[0].eta()

-- this parameter takes values within (-2.5, 2.5)

Or try making 2D distributions.
For example, couple interesting ones are special parameters of taus:
 event_taus_sv_dalitz_m1[0] and event_taus_sv_dalitz_m2[0].
In root interpreter you can see it as:
ttree_out->Draw("event_taus_sv_dalitz_m1[0]:event_taus_sv_dalitz_m2[0]>>h(20,0,2,20,0,2)", "(gen_proc_id == 32 || gen_proc_id == 31)", "col")
'''

input_file.Close()


# the option RECREATE makes it overwrite the existing file
fout = TFile(args.output, "RECREATE")
fout.Write()

histo_sig.Write()
histo_bck.Write()

fout.Write()
fout.Close()

