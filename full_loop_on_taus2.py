import logging
from numpy import sqrt, sign, log
from array import array
from numpy.linalg import norm
from sys import argv


ttree_name = "ntupler/reduced_ttree"
ttree_name = "tauMatch/tauMatch_ttree"

gen_match  = '-g' in argv
gen_match2 = '-g2' in argv
match_quality = '-m' in argv
run_test  = '-t' in argv
run_million  = '-1000000' in argv
logging.basicConfig(level=logging.DEBUG)

logging.debug("loading ROOT...")

# ROOT is heavy, init it after comline
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, gROOT, gSystem, TCanvas, kRed
from ROOT.TMath import Sqrt
gROOT.SetBatch(True)

if '-w' in argv:
    #filename, nick, suff = "/eos/user/o/otoldaie/ttbar-leptons-80X_data/v12.2_merged-sets/MC2016_Summer16_WJets_madgraph.root", 'wjets', ''
    filename, nick, suff = '/afs/cern.ch/user/o/otoldaie/work/private/16/CMSSW_8_0_25/src/UserCode/ttbar-leptons-80X/outdir/v12.6/merged-sets/MC2016_Summer16_WJets_amcatnlo.root', 'wjets', ''
elif '-q' in argv:
    filenames = ['MC2016_Summer16_QCD_HT-100-200.root',
                 'MC2016_Summer16_QCD_HT-1500-2000.root',
                 'MC2016_Summer16_QCD_HT-2000-Inf.root',
                 'MC2016_Summer16_QCD_HT-500-700.root',
                 'MC2016_Summer16_QCD_HT-1000-1500.root',
                 'MC2016_Summer16_QCD_HT-200-300.root',
                 'MC2016_Summer16_QCD_HT-300-500.root',
                 'MC2016_Summer16_QCD_HT-700-1000.root']
    filename = '/afs/cern.ch/user/o/otoldaie/work/private/16/CMSSW_8_0_25/src/UserCode/ttbar-leptons-80X/outdir/v12.6/merged-sets/MC2016_Summer16_QCD_HT-300-500.root'
    nick = 'qcd'
    suff = ''
elif '-d' in argv:
    filename, nick, suff = '/afs/cern.ch/user/o/otoldaie/work/private/16/CMSSW_8_0_25/src/UserCode/ttbar-leptons-80X/outdir/v12.6/merged-sets/Data13TeV_SingleMuon2016B_03Feb2017_ver2.root', 'data', ''
else:
    #filename, nick = "/eos/user/o/otoldaie/ttbar-leptons-80X_data/v12.2_merged-sets/MC2016_Summer16_TTJets_powheg.root" , 'tt'
    #filename, nick = "/afs/cern.ch/user/o/otoldaie/work/private/16/CMSSW_8_0_25/src/UserCode/ttbar-leptons-80X/outdir/v12.4/merged-sets/MC2016_Summer16_TTJets_powheg.root" , 'tt'
    #filename, nick = '/eos/user/o/otoldaie/ttbar-leptons-80X_data/v12.4_merged-sets/MC2016_Summer16_TTJets_powheg.root', 'tt'
    filename, nick, suff = '/afs/cern.ch/user/o/otoldaie/work/private/16/CMSSW_8_0_25/src/UserCode/ttbar-leptons-80X/outdir/v12.6/merged-sets/MC2016_Summer16_TTJets_powheg_1.root', 'tt', '1'
    filename, nick, suff = '~/work/private/16/CMSSW_8_0_26_patch1/src/UserCode/ttbar-leptons-80X/outdir/tauMatch-1/merged-sets/MC2016_Summer16_TTJets_powheg.root', 'tt', ''
    #filename, nick = '/afs/cern.ch/user/o/otoldaie/work/private/16/CMSSW_8_0_25/src/UserCode/ttbar-leptons-80X/outdir/v12.6/merged-sets/MC2016_Summer16_TTJets_powheg.root', 'tt'

isTT = 'TT' in filename

logging.debug("done")
logging.debug("options: %s %s %s %s" % (run_test, isTT, gen_match, gen_match2))

#def PFTau_FlightLength_significance(TVector3 pv,TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov ){
def PFTau_FlightLength_significance(pv,  PVcov, sv, SVcov):
   SVPV = sv - pv
   FD = ROOT.TVectorF()
   FD.ResizeTo(3);
   #FD(0) = SVPV.X();
   #FD(1) = SVPV.Y();
   #FD(2) = SVPV.Z();
   FD.SetElements(array('f', [SVPV.X(), SVPV.Y(), SVPV.Z()]))

   # input covs are
   # ROOT.Math.SMatrix(float, 3, 3, ROOT.Math.MatRepSym(float, 3) )

   #TMatrixT<double> PVcv;
   PVcv = ROOT.TMatrixT(float)()
   PVcv.ResizeTo(3,3);
   #TMatrixT<double> SVcv;
   SVcv = ROOT.TMatrixT(float)()
   SVcv.ResizeTo(3,3);
   #for(int nr =0; nr<PVcov.GetNrows(); nr++){
   #  for(int nc =0; nc<PVcov.GetNcols(); nc++){
   #    PVcv(nr,nc) = PVcov(nr,nc);
   #  }
   #}
   #for(int nr =0; nr<SVcov.GetNrows(); nr++){
   #  for(int nc =0; nc<SVcov.GetNcols(); nc++){
   #    SVcv(nr,nc) = SVcov(nr,nc);
   #  }
   #}
   # assume rows -- first index, coloumns -- second
   for nr in range(3):
      for nc in range(3):
         PVcv[nr][nc] = PVcov(nr, nc)
         SVcv[nr][nc] = SVcov(nr, nc)

   #TMatrixT<double> SVPVMatrix(3,1);
   #SVPVMatrix = ROOT.TMatrixT(float)(3,1) # Error in <operator*=(const TMatrixT &)>: source matrix has wrong shape
   SVPVMatrix = ROOT.TMatrixT(float)(3,3)
   #for(int i=0; i<SVPVMatrix.GetNrows();i++){
   #  SVPVMatrix(i,0)=FD(i);
   #}
   for i in range(SVPVMatrix.GetNrows()):
       SVPVMatrix[i][0] = FD(i)
       SVPVMatrix[i][1] = 0
       SVPVMatrix[i][2] = 0

   #TMatrixT<double> SVPVMatrixT=SVPVMatrix;
   SVPVMatrixT = SVPVMatrix.Clone()
   SVPVMatrixT.T()

   if run_test:
      SVcv.Print()
      PVcv.Print()
   SVcv += PVcv
   PVSVcv = SVcv
   if run_test:
      SVcv.Print()
      PVSVcv.Print()

   if run_test:
       SVPVMatrixT.Print()
       PVSVcv     .Print()
       SVPVMatrix .Print()

   #TMatrixT<double> lambda2 = SVPVMatrixT*(SVcv + PVcv)*SVPVMatrix;
   #lambda2 = SVPVMatrixT*(SVcv + PVcv)*SVPVMatrix
   #lambda2 = SVPVMatrixT*PVSVcv*SVPVMatrix
   # doc says "Compute target = target * source inplace"
   # https://root.cern.ch/doc/v608/classTMatrixT.html
   lambda2 = SVPVMatrixT
   lambda2 *= PVSVcv
   lambda2 *= SVPVMatrix
   if run_test:
       lambda2.Print()

   sigmaabs = Sqrt(lambda2(0,0))/SVPV.Mag()
   sign = SVPV.Mag()/sigmaabs

   return sign, sigmaabs


in_file = TFile(filename, "read")
ntuple  = in_file.Get(ttree_name)

logging.debug("opened file and ntuple")

out_suffix = ''
if match_quality:  out_suffix += 'matchQuality_'
if gen_match:  out_suffix += 'genMatch_'
if gen_match2: out_suffix += 'genMatch2_'
if run_million: out_suffix += 'runMillion_'
out_suffix += nick + '_' + suff
out_filename = "test_full_loop_on_taus2_%s.root" % out_suffix

logging.debug("will write to " + out_filename)
out_file = TFile(out_filename, "recreate")
h_refit_SV_cov00 = TH1D("refit_SV_cov00_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_SV_cov01 = TH1D("refit_SV_cov01_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_SV_cov02 = TH1D("refit_SV_cov02_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_SV_cov10 = TH1D("refit_SV_cov10_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_SV_cov11 = TH1D("refit_SV_cov11_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_SV_cov12 = TH1D("refit_SV_cov12_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_SV_cov20 = TH1D("refit_SV_cov20_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_SV_cov21 = TH1D("refit_SV_cov21_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_SV_cov22 = TH1D("refit_SV_cov22_%s" % nick, "", 100, -0.00005, 0.00005)

h_refit_SV_cov00.SetDirectory(out_file)
h_refit_SV_cov01.SetDirectory(out_file)
h_refit_SV_cov02.SetDirectory(out_file)
h_refit_SV_cov11.SetDirectory(out_file)
h_refit_SV_cov11.SetDirectory(out_file)
h_refit_SV_cov12.SetDirectory(out_file)
h_refit_SV_cov22.SetDirectory(out_file)
h_refit_SV_cov21.SetDirectory(out_file)
h_refit_SV_cov22.SetDirectory(out_file)

h_refit_flightLen_to_other_PV = TH1D("refit_flightLen_to_other_PV_%s" % nick, "", 100, 0, 0.5)
h_refit_flightLen_to_other_PV_large = TH1D("refit_flightLen_to_other_PV_large_%s" % nick, "", 100, 0, 10)

h_refit_PV_to_other_goodPV_small = TH1D("refit_PV_to_other_goodPV_small_%s" % nick, "", 100, 0, 0.02)
h_refit_PV_to_other_goodPV       = TH1D("refit_PV_to_other_goodPV_%s"       % nick, "", 100, 0, 0.5)
h_refit_PV_to_other_goodPV_large = TH1D("refit_PV_to_other_goodPV_large_%s" % nick, "", 100, 0, 10)
h_refit_PV_to_all_goodPV_small   = TH1D("refit_PV_to_all_goodPV_small_%s" % nick, "", 100, 0, 0.02)
h_refit_PV_to_all_goodPV         = TH1D("refit_PV_to_all_goodPV_%s"       % nick, "", 100, 0, 0.5)
h_refit_PV_to_all_goodPV_large   = TH1D("refit_PV_to_all_goodPV_large_%s" % nick, "", 100, 0, 10)

h_refit_flightLen_small = TH1D("refit_flightLen_small_%s" % nick, "", 100, 0, 0.02)
h_refit_flightLen     = TH1D("refit_flightLen_%s" % nick, "", 100, 0, 0.1)
h_refit_flightSign    = TH1D("refit_flightSign_%s" % nick, "", 100, 0, 10)
h_refit_flightLen_lt  = TH1D("refit_flightLen_lt", "", 100, 0, 0.1)
h_refit_flightLen_lj  = TH1D("refit_flightLen_lj", "", 100, 0, 0.1)
h_refit_flightSign_lt = TH1D("refit_flightSign_lt", "", 100, 0, 10)
h_refit_flightSign_lj = TH1D("refit_flightSign_lj", "", 100, 0, 10)
h_refit_flightSign_lj.SetLineColor(kRed)

h_refit_flightErr    = TH1D("refit_flightErr_%s" % nick, "", 100, 0, 0.02)
h_refit_flightErr_lt = TH1D("refit_flightErr_lt", "", 100, 0, 0.02)
h_refit_flightErr_lj = TH1D("refit_flightErr_lj", "", 100, 0, 0.02)

h_pat_flightLen = TH1D("pat_flightLen_%s" % nick, "", 100, 0, 0.1)
h_pat_flightLen_lt = TH1D("pat_flightLen_lt", "", 100, 0, 0.1)
h_pat_flightLen_lj = TH1D("pat_flightLen_lj", "", 100, 0, 0.1)
h_pat_flightSign = TH1D("pat_flightSign_%s" % nick, "", 100, 0, 25)
h_pat_flightSign_lt = TH1D("pat_flightSign_lt", "", 100, 0, 25)
h_pat_flightSign_lj = TH1D("pat_flightSign_lj", "", 100, 0, 25)
h_pat_flightSign_lj.SetLineColor(kRed)

h_pat_flightLen_refit_flightLen    = TH2D("pat_flightLen_refit_flightLen_%s" % nick, "", 100, 0, 0.01, 100, 0, 0.01)
h_pat_flightSign_refit_flightSign  = TH2D("pat_flightSign_refit_flightSign_%s" % nick, "", 100, 0, 10, 100, 0, 10)

h_pat_flightLen_flightSign    = TH2D("pat_flightLen_flightSign_%s" % nick, "", 100, 0, 0.05, 100, 0, 10)
h_pat_flightLen_flightSign_lt = TH2D("pat_flightLen_flightSign_lt",        "", 100, 0, 0.05, 100, 0, 10)
h_pat_flightLen_flightSign_lj = TH2D("pat_flightLen_flightSign_lj",        "", 100, 0, 0.05, 100, 0, 10)
h_pat_bTag_flightSign    = TH2D("pat_bTag_flightSign_%s" % nick, "", 100, 0, 50, 100, 0, 1)
h_pat_bTag_flightSign_lt = TH2D("pat_bTag_flightSign_lt",        "", 100, 0, 50, 100, 0, 1)
h_pat_bTag_flightSign_lj = TH2D("pat_bTag_flightSign_lj",        "", 100, 0, 50, 100, 0, 1)
#tau_dR_matched_jet

h_refit_m1m2    = TH2D("refit_m1m2_%s" % nick, "", 100, 0, 2, 100, 0, 2)
h_refit_m1m2_lt = TH2D("refit_m1m2_lt", "", 100, 0, 2, 100, 0, 2)
h_refit_m1m2_lj = TH2D("refit_m1m2_lj", "", 100, 0, 2, 100, 0, 2)

h_refit_flightLen_flightSign    = TH2D("refit_flightLen_flightSign_%s" % nick, "", 100, 0, 0.05, 100, 0, 10)
h_refit_flightLen_flightSign_lt = TH2D("refit_flightLen_flightSign_lt",        "", 100, 0, 0.05, 100, 0, 10)
h_refit_flightLen_flightSign_lj = TH2D("refit_flightLen_flightSign_lj",        "", 100, 0, 0.05, 100, 0, 10)
h_refit_bTag_flightSign    = TH2D("refit_bTag_flightSign_%s" % nick, "", 100, 0, 10, 100, 0, 1.01)
h_refit_bTag_flightSign_lt = TH2D("refit_bTag_flightSign_lt",        "", 100, 0, 10, 100, 0, 1.01)
h_refit_bTag_flightSign_lj = TH2D("refit_bTag_flightSign_lj",        "", 100, 0, 10, 100, 0, 1.01)
#tau_dR_matched_jet

h_refit_flightSign_lt_Btag = TH1D("refit_flightSign_lt_Btag", "", 100, 0, 10)
h_refit_flightSign_lj_Btag = TH1D("refit_flightSign_lj_Btag", "", 100, 0, 10)
h_pat_flightSign_lt_Btag   = TH1D("pat_flightSign_lt_Btag",   "", 100, 0, 25)
h_pat_flightSign_lj_Btag   = TH1D("pat_flightSign_lj_Btag",   "", 100, 0, 25)
h_pat_flightSign_lj_Btag.SetLineColor(kRed)

h_refit_flightSign_lt_noBtag = TH1D("refit_flightSign_lt_noBtag", "", 100, 0, 10)
h_refit_flightSign_lj_noBtag = TH1D("refit_flightSign_lj_noBtag", "", 100, 0, 10)
h_pat_flightSign_lt_noBtag   = TH1D("pat_flightSign_lt_noBtag",   "", 100, 0, 25)
h_pat_flightSign_lj_noBtag   = TH1D("pat_flightSign_lj_noBtag",   "", 100, 0, 25)
h_pat_flightSign_lj_noBtag.SetLineColor(kRed)

# flightLen vs tau Energy
h_refit_flightLen_Energy    = TH2D("refit_flightLen_Energy_%s" % nick, "", 10, 0, 0.02, 10, 0, 150)
h_refit_flightLen_Energy_lt = TH2D("refit_flightLen_Energy_lt",        "", 10, 0, 0.02, 10, 0, 150)
h_refit_flightLen_Energy_lt_large = TH2D("refit_flightLen_Energy_lt_large",        "", 10, 0, 0.1, 10, 0, 150)
h_refit_flightLen_Energy_lj = TH2D("refit_flightLen_Energy_lj",        "", 10, 0, 0.02, 10, 0, 150)
h_pat_flightLen_Energy      = TH2D("pat_flightLen_Energy_%s" % nick,   "", 10, 0, 0.02, 10, 0, 150)
h_pat_flightLen_Energy_lt   = TH2D("pat_flightLen_Energy_lt",          "", 10, 0, 0.02, 10, 0, 150)
h_pat_flightLen_Energy_lj   = TH2D("pat_flightLen_Energy_lj",          "", 10, 0, 0.02, 10, 0, 150)
#tau_dR_matched_jet

# PV XY Z cov
h_refit_PV_xy  = TH2D("refit_PV_xy_%s"  % nick, "", 100, 0.1048-0.01, 0.1048+0.01, 100, 0.1686-0.01, 0.1686+0.01)
h_refit_PV_z   = TH1D("refit_PV_z_%s"   % nick, "", 100, -50, 50)
h_refit_PV_cov00 = TH1D("refit_PV_cov00_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_PV_cov01 = TH1D("refit_PV_cov01_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_PV_cov02 = TH1D("refit_PV_cov02_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_PV_cov10 = TH1D("refit_PV_cov10_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_PV_cov11 = TH1D("refit_PV_cov11_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_PV_cov12 = TH1D("refit_PV_cov12_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_PV_cov20 = TH1D("refit_PV_cov20_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_PV_cov21 = TH1D("refit_PV_cov21_%s" % nick, "", 100, -0.00005, 0.00005)
h_refit_PV_cov22 = TH1D("refit_PV_cov22_%s" % nick, "", 100, -0.00005, 0.00005)

h_n_goodPV = TH1D("n_goodPV_%s" % nick, "", 50, 0, 50)

logging.debug("set histograms")

N_events = 1000000 # ntuple.GetEntries()
n_events = 0

tau_ID_lev = 1 # >3 = Tight

# looping:
for i, event in enumerate(ntuple):
    if run_test and i > 1: break
    #if not (event.tau_SV_fit_isOk.size() > 0 and event.tau_SV_fit_isOk[0]>0 and event.PV_fit_isOk>0 and (event.tau_SV_fit_matchingQuality[0] < 0.001 or not match_quality)):
    if not (event.tau_IDlev[0] > tau_ID_lev and event.tau_decayMode[0] > 9 and event.tau_refited_index.size() > 0 and event.tau_refited_index[0]>-1 and event.PV_fit_isOk>0 and (event.tau_SV_fit_matchingQuality[0] < 0.001 or not match_quality)):
        continue
    # check tracks are present:
    #if not (event.tau_SV_fit_track_SS_p4.size() > 0 and event.tau_SV_fit_track_OS1_p4.size() > 0 and event.tau_SV_fit_track_OS2_p4.size() > 0):
    if not (event.tau_SV_fit_track_OS_p4.size() > 0 and event.tau_SV_fit_track_SS1_p4.size() > 0 and event.tau_SV_fit_track_SS2_p4.size() > 0):
        continue
    n_events += 1
    if run_million and n_events > N_events: break

    # match the 0 tau to any of gen taus
    if gen_match and not any(sqrt((gen.Eta() - event.tau_p4[0].Eta())**2 + (gen.Phi() - event.tau_p4[0].Phi())**2) < 0.1 for gen in event.gen_tau_p4): continue
    if gen_match2:
       #and not any(sqrt((gen.Eta() - event.tau_p4[0].Eta())**2 + (gen.Phi() - event.tau_p4[0].Phi())**2) < 0.1 for gen in event.gen_tau_p4): continue
       matched = False
       if abs(event.gen_t_pt) > 15*15: # hadronic tau
           matched = ((event.gen_t_w_final_p4.Eta() - event.tau_p4[0].Eta())**2 + (event.gen_t_w_final_p4.Phi() - event.tau_p4[0].Phi())**2) < 0.04  # 0.2 dR
       if abs(event.gen_tb_pt) > 15*15: # hadronic tau
           matched = ((event.gen_tb_w_final_p4.Eta() - event.tau_p4[0].Eta())**2 + (event.gen_tb_w_final_p4.Phi() - event.tau_p4[0].Phi())**2) < 0.04
       if not matched: continue
    # ref:
    # http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_5_3_9/doc/html/d5/d6b/DataFormats_2Math_2interface_2deltaR_8h_source.html#l00019

    h_n_goodPV.Fill(event.PV_x.size())
    #h_refit_SV_cov00.Fill(sign(event.tau_SV_cov[0](0, 0)) *log(abs(event.tau_SV_cov[0](0, 0))))
    #h_refit_SV_cov01.Fill(sign(event.tau_SV_cov[0](0, 1)) *log(abs(event.tau_SV_cov[0](0, 1))))
    #h_refit_SV_cov02.Fill(sign(event.tau_SV_cov[0](0, 2)) *log(abs(event.tau_SV_cov[0](0, 2))))
    #h_refit_SV_cov10.Fill(sign(event.tau_SV_cov[0](1, 1)) *log(abs(event.tau_SV_cov[0](1, 1))))
    #h_refit_SV_cov11.Fill(sign(event.tau_SV_cov[0](1, 1)) *log(abs(event.tau_SV_cov[0](1, 1))))
    #h_refit_SV_cov12.Fill(sign(event.tau_SV_cov[0](1, 2)) *log(abs(event.tau_SV_cov[0](1, 2))))
    #h_refit_SV_cov20.Fill(sign(event.tau_SV_cov[0](2, 2)) *log(abs(event.tau_SV_cov[0](2, 2))))
    #h_refit_SV_cov21.Fill(sign(event.tau_SV_cov[0](2, 1)) *log(abs(event.tau_SV_cov[0](2, 1))))
    #h_refit_SV_cov22.Fill(sign(event.tau_SV_cov[0](2, 2)) *log(abs(event.tau_SV_cov[0](2, 2))))
    #
    h_refit_SV_cov00.Fill(event.tau_SV_cov[0](0, 0))
    h_refit_SV_cov01.Fill(event.tau_SV_cov[0](0, 1))
    h_refit_SV_cov02.Fill(event.tau_SV_cov[0](0, 2))
    h_refit_SV_cov10.Fill(event.tau_SV_cov[0](1, 0))
    h_refit_SV_cov11.Fill(event.tau_SV_cov[0](1, 1))
    h_refit_SV_cov12.Fill(event.tau_SV_cov[0](1, 2))
    h_refit_SV_cov20.Fill(event.tau_SV_cov[0](2, 0))
    h_refit_SV_cov21.Fill(event.tau_SV_cov[0](2, 1))
    h_refit_SV_cov22.Fill(event.tau_SV_cov[0](2, 2))

    '''
    pv = array((event.PV_x[0], event.PV_y[0], event.PV_z[0]))
    pv_err = array((event.PV_x_err[0], event.PV_y_err[0], event.PV_z_err[0]))
    sv = array((event.tau_SV_fit_x[0], event.tau_SV_fit_y[0], event.tau_SV_fit_z[0]))
    sv_err = array((event.tau_SV_cov[0](0,0), event.tau_SV_cov[0](1,1), event.tau_SV_cov[0](2,2)))

    flight_vec = pv - sv
    flight_len = norm(flight_vec)
    flight_vec_err = sqrt(sv_err**2 + pv_err**2)
    flight_sign = norm(flight_vec / flight_vec_err)
    '''

    h_refit_PV_xy.Fill(event.PV_fit_x, event.PV_fit_y)
    h_refit_PV_z .Fill(event.PV_fit_z)
    #h_refit_PV_cov00.Fill(sign(event.PV_cov(0, 0))*log(abs(event.PV_cov(0, 0))))
    #h_refit_PV_cov01.Fill(sign(event.PV_cov(0, 1))*log(abs(event.PV_cov(0, 1))))
    #h_refit_PV_cov02.Fill(sign(event.PV_cov(0, 2))*log(abs(event.PV_cov(0, 2))))
    #h_refit_PV_cov10.Fill(sign(event.PV_cov(1, 0))*log(abs(event.PV_cov(1, 0))))
    #h_refit_PV_cov11.Fill(sign(event.PV_cov(1, 1))*log(abs(event.PV_cov(1, 1))))
    #h_refit_PV_cov12.Fill(sign(event.PV_cov(1, 2))*log(abs(event.PV_cov(1, 2))))
    #h_refit_PV_cov20.Fill(sign(event.PV_cov(2, 0))*log(abs(event.PV_cov(2, 0))))
    #h_refit_PV_cov21.Fill(sign(event.PV_cov(2, 1))*log(abs(event.PV_cov(2, 1))))
    #h_refit_PV_cov22.Fill(sign(event.PV_cov(2, 2))*log(abs(event.PV_cov(2, 2))))
    #

    h_refit_PV_cov00.Fill(event.PV_cov(0, 0))
    h_refit_PV_cov01.Fill(event.PV_cov(0, 1))
    h_refit_PV_cov02.Fill(event.PV_cov(0, 2))
    h_refit_PV_cov10.Fill(event.PV_cov(1, 0))
    h_refit_PV_cov11.Fill(event.PV_cov(1, 1))
    h_refit_PV_cov12.Fill(event.PV_cov(1, 2))
    h_refit_PV_cov20.Fill(event.PV_cov(2, 0))
    h_refit_PV_cov21.Fill(event.PV_cov(2, 1))
    h_refit_PV_cov22.Fill(event.PV_cov(2, 2))

    pv = ROOT.TVector3(event.PV_fit_x, event.PV_fit_y, event.PV_fit_z)
    sv = ROOT.TVector3(event.tau_SV_fit_x[0], event.tau_SV_fit_y[0], event.tau_SV_fit_z[0])
    flight_len = (sv-pv).Mag()
    if run_test: logging.debug("flight len  = %f" % flight_len)

    ##def PFTau_FlightLength_significance(TVector3 pv,TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov ){
    try:
        flight_sign, flight_err = PFTau_FlightLength_significance(pv, event.PV_cov, sv, event.tau_SV_cov[0])
    except ZeroDivisionError:
        print "division by zero in PFTau_FlightLength_significance"

    if run_test: logging.debug("flight sign = %f" % flight_sign)

    pat_flightLen  = event.tau_flightLength[0]
    pat_flightSign = event.tau_flightLengthSignificance[0]

    if event.PV_x.size() > 1:
      # check distances to other PVs
      for i in range(event.PV_x.size()):
        pv_i = ROOT.TVector3(event.PV_x[i], event.PV_y[i], event.PV_z[i])
        pv_dist = (pv_i - pv).Mag()
        h_refit_PV_to_all_goodPV_small.Fill(pv_dist)
        h_refit_PV_to_all_goodPV      .Fill(pv_dist)
        h_refit_PV_to_all_goodPV_large.Fill(pv_dist)

        if i == 0: continue
        h_refit_PV_to_other_goodPV_small.Fill(pv_dist)
        h_refit_PV_to_other_goodPV      .Fill(pv_dist)
        h_refit_PV_to_other_goodPV_large.Fill(pv_dist)
        #pv_i = array('f', (event.PV_x[i], event.PV_y[i], event.PV_z[i]))
        flight_len_i = (pv_i - sv).Mag()
        h_refit_flightLen_to_other_PV.Fill(flight_len_i)
        h_refit_flightLen_to_other_PV_large.Fill(flight_len_i)

    # masses of tracks
    #mass1 = (event.tau_SV_fit_track_SS_p4[0] + event.tau_SV_fit_track_OS1_p4[0]).mass()
    #mass2 = (event.tau_SV_fit_track_SS_p4[0] + event.tau_SV_fit_track_OS2_p4[0]).mass()
    mass1 = (event.tau_SV_fit_track_OS_p4[0] + event.tau_SV_fit_track_SS1_p4[0]).mass()
    mass2 = (event.tau_SV_fit_track_OS_p4[0] + event.tau_SV_fit_track_SS2_p4[0]).mass()

    h_refit_m1m2.Fill(mass1, mass2)

    #h_pat_bTag_flightSign_lt = TH2D("pat_bTag_flightSign_lt",        "", 100, 0, 1, 100, 0, 50)
    #h_pat_bTag_flightSign_lj = TH2D("pat_bTag_flightSign_lj",        "", 100, 0, 1, 100, 0, 50)
    ##tau_dR_matched_jet
    h_pat_flightLen .Fill(pat_flightLen)
    h_pat_flightSign.Fill(pat_flightSign)
    h_refit_flightLen_small .Fill(flight_len)
    h_refit_flightLen .Fill(flight_len)
    h_refit_flightSign.Fill(flight_sign)
    h_refit_flightErr.Fill(flight_err)

    h_pat_flightLen_flightSign  .Fill(pat_flightLen, pat_flightSign)
    h_refit_flightLen_flightSign.Fill(flight_len, flight_sign)

    h_pat_flightLen_refit_flightLen  .Fill(pat_flightLen, flight_len)
    h_pat_flightSign_refit_flightSign.Fill(pat_flightSign, flight_sign)

    h_refit_flightLen_Energy .Fill(flight_len,    event.tau_p4[0].energy())
    h_pat_flightLen_Energy   .Fill(pat_flightLen, event.tau_p4[0].energy())

    i_matched_tau_jet = event.tau_dR_matched_jet[0]
    if i_matched_tau_jet > -1:
        h_pat_bTag_flightSign.Fill(pat_flightSign, event.jet_b_discr[i_matched_tau_jet])
        h_refit_bTag_flightSign.Fill(flight_sign, event.jet_b_discr[i_matched_tau_jet])

    if isTT:
        if abs(event.gen_t_w_decay_id * event.gen_tb_w_decay_id) == 13: # lj
            h_refit_flightLen_lj .Fill(flight_len)
            h_refit_flightSign_lj.Fill(flight_sign)
            h_refit_flightErr_lj.Fill(flight_err)
            h_pat_flightLen_lj .Fill(pat_flightLen)
            h_pat_flightSign_lj.Fill(pat_flightSign)

            h_refit_flightLen_Energy_lj .Fill(flight_len,    event.tau_p4[0].energy())
            h_pat_flightLen_Energy_lj   .Fill(pat_flightLen, event.tau_p4[0].energy())

            h_pat_flightLen_flightSign_lj  .Fill(pat_flightLen, pat_flightSign)
            h_refit_flightLen_flightSign_lj.Fill(flight_len, flight_sign)

            h_refit_m1m2_lj.Fill(mass1, mass2)

            if i_matched_tau_jet > -1:
                h_pat_bTag_flightSign_lj.Fill(pat_flightSign, event.jet_b_discr[i_matched_tau_jet])
                h_refit_bTag_flightSign_lj.Fill(flight_sign, event.jet_b_discr[i_matched_tau_jet])
                if event.jet_b_discr[i_matched_tau_jet] < 0.85:
                    h_refit_flightSign_lj_noBtag .Fill(flight_sign)
                    h_pat_flightSign_lj_noBtag   .Fill(pat_flightSign)
                else:
                    h_refit_flightSign_lj_Btag .Fill(flight_sign)
                    h_pat_flightSign_lj_Btag   .Fill(pat_flightSign)

        if (abs(event.gen_t_w_decay_id) > 15*15 and abs(event.gen_tb_w_decay_id) == 13) or (abs(event.gen_t_w_decay_id) == 13 and abs(event.gen_tb_w_decay_id) > 15*15): # lt
            h_refit_flightLen_lt .Fill(flight_len)
            h_refit_flightSign_lt.Fill(flight_sign)
            h_refit_flightErr_lt.Fill(flight_err)
            h_pat_flightLen_lt .Fill(pat_flightLen)
            h_pat_flightSign_lt.Fill(pat_flightSign)

            h_refit_flightLen_Energy_lt_large .Fill(flight_len,    event.tau_p4[0].energy())
            h_refit_flightLen_Energy_lt .Fill(flight_len,    event.tau_p4[0].energy())
            h_pat_flightLen_Energy_lt   .Fill(pat_flightLen, event.tau_p4[0].energy())

            h_pat_flightLen_flightSign_lt  .Fill(pat_flightLen, pat_flightSign)
            h_refit_flightLen_flightSign_lt.Fill(flight_len, flight_sign)

            h_refit_m1m2_lt.Fill(mass1, mass2)

            if i_matched_tau_jet > -1:
                h_pat_bTag_flightSign_lt.Fill(pat_flightSign, event.jet_b_discr[i_matched_tau_jet])
                h_refit_bTag_flightSign_lt.Fill(flight_sign, event.jet_b_discr[i_matched_tau_jet])
                if event.jet_b_discr[i_matched_tau_jet] < 0.85:
                    h_refit_flightSign_lt_noBtag .Fill(flight_sign)
                    h_pat_flightSign_lt_noBtag   .Fill(pat_flightSign)
                else:
                    h_refit_flightSign_lt_Btag .Fill(flight_sign)
                    h_pat_flightSign_lt_Btag   .Fill(pat_flightSign)



logging.debug("writing output")

out_file.Write()

