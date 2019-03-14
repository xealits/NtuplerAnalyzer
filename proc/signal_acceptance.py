import argparse
import logging
import os
from os.path import isfile, basename
import ctypes
from array import array
from sys import exit




parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "signal acceptance study",
    epilog = """Example:
python signal_acceptance.py temp/ gstore_outdirs/v34/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v34_MC2016_Summer16_TTJets_powheg/180717_001103/0000/MC2016_Summer16_TTJets_powheg_90.root
python proc/signal_acceptance.py temp/ NtuplerAnalyzer_test_METfiltersOFF_TTJets2_signal_recordAll.root
"""
    )

parser.add_argument("output_dir",  type=str, default="temp/", help="the path of output directory")
parser.add_argument("--debug",  action='store_true', help="DEBUG level of logging")
parser.add_argument("--no-tau-cut",  action='store_true', help="don't apply tau cut")

parser.add_argument('input_file', help="""the job files to process: 
/gstore/t3cms/store/user/otoldaie/v19/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/Ntupler_v19_MC2016_Summer16_TTJets_powheg/180226_022336/0000/MC2016_Summer16_TTJets_powheg_0.root""")


args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

input_filename = args.input_file
if not isfile(input_filename):
    logging.info("missing: " + input_filename)
    exit(1)

job_filename = input_filename.split('/')[-1]
output_path = args.output_dir + '/' + job_filename
logging.debug(output_path)

assert not isfile(output_path)

logging.info("import ROOT")

import ROOT
from ROOT import TH1D, TFile, TTree, gROOT
from ROOT import kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, TColor
from ROOT import TLegend
#from ROOT import THStack

logging.info("done")

from gen_proc_defs import MUf_nom_MUr_nom  , MUf_up_MUr_nom   , MUf_down_MUr_nom , MUf_nom_MUr_up   , MUf_up_MUr_up    , MUf_down_MUr_up  , MUf_nom_MUr_down , MUf_up_MUr_down  , MUf_down_MUr_down


# pu weights

pileup_ratio = array('d', [
0.360609416811339, 0.910848525427002, 1.20629960507795, 0.965997726573782, 1.10708082813183, 1.14843491548622, 0.786526251164482, 0.490577792661333, 0.740680941110478,
0.884048630953726,
0.964813189764159, 1.07045369167689, 1.12497267309738, 1.17367530613108, 1.20239808206413, 1.20815108390021, 1.20049333094509, 1.18284686347315, 1.14408796655615,
1.0962284704313, 1.06549162803223, 1.05151011089581, 1.05159666626121, 1.05064452078328, 1.0491726301522, 1.05772537082991, 1.07279673875566, 1.0837536468865, 1.09536667397119,
1.10934472980173, 1.09375894592864, 1.08263679568271, 1.04289345879947, 0.9851490341672, 0.909983816540809, 0.821346330143864, 0.71704523475871, 0.609800913869359, 0.502935245638477,
0.405579825620816, 0.309696044611377, 0.228191137503131, 0.163380359253309, 0.113368437957202, 0.0772279997453792, 0.0508111733313502, 0.0319007262683943, 0.0200879459309245, 0.0122753366005436,
0.00739933885813127, 0.00437426967257811, 0.00260473545284139, 0.00157047254226743, 0.000969500595715493, 0.000733193118123283, 0.000669817107713128, 0.000728548958604492, 0.000934559691182011, 0.00133719688378802,
0.00186652283903214, 0.00314422244976771, 0.00406954793369611, 0.00467888840511915, 0.00505224284441512, 0.00562827194936864, 0.0055889504870752, 0.00522867039470319, 0.00450752163476433, 0.00395300774604375,
0.00330577167682956, 0.00308353042577215, 0.00277846504893301, 0.00223943190687725, 0.00196650068765464, 0.00184742734258922,])

pileup_ratio_up = array('d', [0.351377216124927, 0.717199649125846, 1.14121536968772, 0.84885826611733, 1.00700929402897, 1.03428595270903, 0.717444379696992, 0.344078389355127, 0.499570875027422,
0.606614916257104, 0.632584599390169, 0.731450949466174, 0.827511723989754, 0.910682115553867, 0.960170981598162, 0.988896170761361, 1.02468865580207, 1.05296667126403, 1.05112033565679,
1.0269129153969, 1.00548641752714, 0.998316130432865, 1.01492587998551, 1.03753749807849, 1.05742218946485, 1.08503978097083, 1.12134132247053, 1.15585339474274, 1.19214399856171,
1.23308400947467, 1.24528804633732, 1.26786364716917, 1.26101551498967, 1.23297806722714, 1.18042533075471, 1.10534683838101, 1.00275591661645, 0.889094305531985, 0.768791254270252,
0.655054015673457, 0.533361034358457, 0.423095146361996, 0.329177839117034, 0.250352385505809, 0.188377378855567, 0.137852651411779, 0.0968577167707531, 0.0686240187247059, 0.0473889635126706,
0.0323695027438475, 0.0216752397914914, 0.0145119352923332, 0.00961177893634792, 0.00615582219138384, 0.00430085627914427, 0.00305735512896403, 0.00223567790438986, 0.00189369737638594, 0.00199585978316291,
0.00236236592656064, 0.00372472999463276, 0.00474687312579969, 0.00549508151576102, 0.00603023110946686, 0.0068545111910253, 0.00695838760530896, 0.00666224781277046, 0.00588243140681038, 0.00528714370892014,
0.00453424615273565, 0.00433985030329723, 0.00401493171035719, 0.00332436608713241, 0.00300063798808221, 0.00289925128977536,])

pileup_ratio_down = array('d', [0.37361294640242, 1.1627791004568, 1.26890787896295, 1.10266790442705, 1.23456697093644, 1.26278991594152, 0.909648777562084, 0.759569490571151, 1.09035651921682,
1.34530547603283, 1.48713160105, 1.52535976889483, 1.49730550773404, 1.49792998045778, 1.49767851097519, 1.44431045398336, 1.3681909492045, 1.29912252494785, 1.2274279217797,
1.16525969099909, 1.12531044676724, 1.09094501417685, 1.06405434433422, 1.03997120824565, 1.0185716022098, 1.00560949501652, 0.997570939806059, 0.985543761409897, 0.972557804582185,
0.957832827239337, 0.9139572640153, 0.872252387173971, 0.808388185417578, 0.733817960498049, 0.650440963845892, 0.561688505024782, 0.466564380334112, 0.374428618658619, 0.28845274688129,
0.214909665968644, 0.149991974352384, 0.100014138338029, 0.0642260884603397, 0.0396553405911344, 0.0238687936736627, 0.0137921542898078, 0.00756854010632403, 0.00415483516246187, 0.00221776872027937,
0.00118249725637452, 0.000641889697310868, 0.000383647166012176, 0.000273637590071334, 0.000242902582071058, 0.000291239677209452, 0.000394091114279828, 0.000542541231466254, 0.000771067920964491, 0.00113596447675107,
0.00158061353194779, 0.00261959852500539, 0.00331800452823827, 0.00372426930370732, 0.00392086545082614, 0.00425479965493548, 0.00411256966391362, 0.00374240422174387, 0.00313603438166934, 0.00267155793176928,
0.00216878198028599, 0.00196249821290853, 0.00171433839159669, 0.00133866519755926, 0.00113810604240254, 0.00103447940224886,])



pileup_ratio_ele = array('d', [
   0.413231   ,    1.01701    ,    1.19502     ,   0.883906   ,    1.05852    ,    1.11823    ,    0.789439   ,    0.515477   ,    0.81338    ,    0.990148   ,
   1.0919     ,    1.21784    ,    1.28268     ,   1.33936    ,    1.37267    ,    1.38001    ,    1.37224    ,    1.35253    ,    1.30805    ,    1.25303    ,
   1.21761    ,    1.20085    ,    1.1987      ,   1.19257    ,    1.1807     ,    1.17079    ,    1.15238    ,    1.10667    ,    1.03375    ,    0.935086   ,
   0.793376   ,    0.65125    ,    0.502727    ,   0.369298   ,    0.25859    ,    0.173207   ,    0.110361   ,    0.0677957  ,    0.0403186  ,    0.0236369  ,
   0.0133546  ,    0.00746494 ,    0.00417626  ,   0.00233773 ,    0.0013288  ,    0.000757718,    0.000432788,    0.000266239,    0.000177605,    0.000137241,
   0.000125696,    0.000137018,    0.000167806 ,   0.000215108,    0.000313214,    0.000464376,    0.000669855,    0.000981399,    0.00148275 ,    0.00211313 ,
   0.0035872  ,    0.00465614 ,    0.005359    ,   0.00578897 ,    0.00645001 ,    0.00640537 ,    0.00599263 ,    0.00516618 ,    0.00453067 ,    0.00378886 ,
   0.00353415 ,    0.00318451 ,    0.0025667   ,   0.00225388 ,    0.00211741 ,])

pileup_ratio_up_ele = array('d', [
   0.402668    ,   0.803377   ,    1.15963    ,    0.764147   ,    0.966328    ,   0.995159   ,    0.71563    ,    0.354304   ,   0.541943   ,    0.674778   ,
   0.713035    ,   0.830366   ,    0.942616   ,    1.03882    ,    1.09589     ,   1.12909    ,    1.17068    ,    1.20376    ,   1.20191    ,    1.17404    ,
   1.14928     ,   1.14083    ,    1.15906    ,    1.18279    ,    1.20082     ,   1.22277    ,    1.24559    ,    1.25129    ,   1.23607    ,    1.19534    ,
   1.09539     ,   0.978694   ,    0.82546    ,    0.662451   ,    0.50547     ,   0.367764   ,    0.25369    ,    0.168007   ,   0.10706    ,    0.0667404  ,
   0.0397874   ,   0.0233291  ,    0.0136533  ,    0.00799737 ,    0.00476279  ,   0.00284044 ,    0.00167744 ,    0.00103389 ,   0.000648432,    0.000427764,
   0.000303899 ,   0.000247672,    0.000236803,    0.000258026,    0.000345092 ,   0.000494341,    0.000708128,    0.00104444 ,   0.00159927 ,    0.00231779 ,
   0.00400894  ,   0.00530831 ,    0.00623822 ,    0.00688571 ,    0.00784455  ,   0.00797042 ,    0.00763388 ,    0.00674129 ,   0.00605947 ,    0.00519674 ,
   0.00497399  ,   0.00460162 ,    0.00381017 ,    0.00343914 ,    0.00332295  ,])

pileup_ratio_down_ele = array('d', [
   0.428107   ,   1.29388    ,    1.22078   ,    1.02596    ,    1.1769     ,    1.24377    ,    0.921862   ,    0.814769   ,    1.20901   ,    1.51527   ,
   1.68838    ,   1.73792    ,    1.70826   ,    1.70984    ,    1.71038    ,    1.65067    ,    1.56442    ,    1.48535    ,    1.40303   ,    1.33164   ,
   1.28514    ,   1.24342    ,    1.20714   ,    1.16839    ,    1.12262    ,    1.06993    ,    0.999693   ,    0.900043   ,    0.778486  ,    0.644942  ,
   0.497564   ,   0.37052    ,    0.259917  ,    0.174109   ,    0.111585   ,    0.0687061  ,    0.0404941  ,    0.0232033  ,    0.0129877 ,    0.00721863,
   0.00388086 ,   0.00206418 ,    0.00109735,    0.000585006,    0.000321242,    0.000184087,    0.000114353,    8.57172e-5 ,    7.68135e-5,    8.09537e-5,
   9.42381e-5 ,   0.000118387,    0.0001548 ,    0.000202628,    0.000294603,    0.000431519,    0.00061181 ,    0.000878668,    0.00129927,    0.0018102 ,
   0.00300153 ,   0.00380243 ,    0.0042683 ,    0.00449375 ,    0.00487654 ,    0.00471355 ,    0.0042893  ,    0.00359432 ,    0.00306198,    0.00248572,
   0.0022493  ,   0.00196487 ,    0.0015343 ,    0.00130443 ,    0.00118566 ,])


'''
systematic weights to implement

# renorm refact scales
"M_NOM"    : 10,
"MrUp"     : 11,
"MrDown"   : 12,
"MfUp"     : 13,
"MfDown"   : 14,
"MfrUp"    : 15,
"MfrDown"  : 16,

# PDF and alphaS
"PDF_NOM"       : 27,
"AlphaSUp"      : 28,
"AlphaSDown"    : 29,

and the pdf sets

and the separate datasets
'''





# for all systematic weights
histos = {
'eltau':     (TH1D("gen_tbw_eltau",     "", 400, -200, 200), TH1D("gen_tw_eltau",     "", 400, -200, 200), TH1D("all_ev_eltau",     "", 200, 0, 200), TH1D("cut_ev_eltau",     "", 200, 0, 200)),
'elj':       (TH1D("gen_tbw_elj",       "", 400, -200, 200), TH1D("gen_tw_elj",       "", 400, -200, 200), TH1D("all_ev_elj",       "", 200, 0, 200), TH1D("cut_ev_elj",       "", 200, 0, 200)),
'tauelj':    (TH1D("gen_tbw_tauelj",    "", 400, -200, 200), TH1D("gen_tw_tauelj",    "", 400, -200, 200), TH1D("all_ev_tauelj",    "", 200, 0, 200), TH1D("cut_ev_tauelj",    "", 200, 0, 200)),
'taueltauh': (TH1D("gen_tbw_taueltauh", "", 400, -200, 200), TH1D("gen_tw_taueltauh", "", 400, -200, 200), TH1D("all_ev_taueltauh", "", 200, 0, 200), TH1D("cut_ev_taueltauh", "", 200, 0, 200)),
'mutau':     (TH1D("gen_tbw_mutau",     "", 400, -200, 200), TH1D("gen_tw_mutau",     "", 400, -200, 200), TH1D("all_ev_mutau",     "", 200, 0, 200), TH1D("cut_ev_mutau",     "", 200, 0, 200)),
'muj':       (TH1D("gen_tbw_muj",       "", 400, -200, 200), TH1D("gen_tw_muj",       "", 400, -200, 200), TH1D("all_ev_muj",       "", 200, 0, 200), TH1D("cut_ev_muj",       "", 200, 0, 200)),
'taumuj':    (TH1D("gen_tbw_taumuj",    "", 400, -200, 200), TH1D("gen_tw_taumuj",    "", 400, -200, 200), TH1D("all_ev_taumuj",    "", 200, 0, 200), TH1D("cut_ev_taumuj",    "", 200, 0, 200)),
'taumutauh': (TH1D("gen_tbw_taumutauh", "", 400, -200, 200), TH1D("gen_tw_taumutauh", "", 400, -200, 200), TH1D("all_ev_taumutauh", "", 200, 0, 200), TH1D("cut_ev_taumutauh", "", 200, 0, 200)),
'other':     (TH1D("gen_tbw_other",     "", 400, -200, 200), TH1D("gen_tw_other",     "", 400, -200, 200), TH1D("all_ev_other",     "", 200, 0, 200), TH1D("cut_ev_other",     "", 200, 0, 200))}

# some control distrs
control_histos = {
'pre_gen_lep_pt': TH1D("pre_gen_lep_pt",    "", 90, 20, 200),
'pre_gen_el_pt':  TH1D("pre_gen_el_pt",     "", 90, 20, 200),
'pre_gen_mu_pt':  TH1D("pre_gen_mu_pt",     "", 90, 20, 200),
'gen_lep_pt':     TH1D("gen_lep_pt",    "", 90, 20, 200),
'gen_el_pt':      TH1D("gen_el_pt",     "", 90, 20, 200),
'gen_mu_pt':      TH1D("gen_mu_pt",     "", 90, 20, 200),
}

#el_histos = (histo_all_ev_eltau, histo_cut_ev_eltau)
#mu_histos = (histo_all_ev_mutau, histo_cut_ev_mutau)

def round_pdf_weight(pdf_w):
    if pdf_w > 2.:
        pdf_w = 1.
    return pdf_w

logging.debug(input_filename)

tfile = TFile(input_filename)
ttree = tfile.Get('ntupler/reduced_ttree')

for iev, event in enumerate(ttree):
    is_electron_process = abs(event.gen_t_w_decay_id) == 11 or abs(event.gen_tb_w_decay_id) == 11 or abs(event.gen_t_w_decay_id) == 15*11 or abs(event.gen_tb_w_decay_id) == 15*11
    is_eltau, is_mutau = False, False
    if   (abs(event.gen_t_w_decay_id) == 11    and abs(event.gen_tb_w_decay_id) > 15*15) or (abs(event.gen_tb_w_decay_id) == 11    and abs(event.gen_t_w_decay_id) > 15*15):
        process = 'eltau'
        is_eltau = True
    elif (abs(event.gen_t_w_decay_id) == 11*15 and abs(event.gen_tb_w_decay_id) > 15*15) or (abs(event.gen_tb_w_decay_id) == 11*15 and abs(event.gen_t_w_decay_id) > 15*15):
        process = 'taueltauh'
    elif (abs(event.gen_t_w_decay_id) == 11    and abs(event.gen_tb_w_decay_id) == 1)    or (abs(event.gen_tb_w_decay_id) == 11    and abs(event.gen_t_w_decay_id) == 1):
        process = 'elj'
    elif (abs(event.gen_t_w_decay_id) == 15*11 and abs(event.gen_tb_w_decay_id) == 1)    or (abs(event.gen_tb_w_decay_id) == 15*11 and abs(event.gen_t_w_decay_id) == 1):
        process = 'tauelj'

    elif (abs(event.gen_t_w_decay_id) == 13    and abs(event.gen_tb_w_decay_id) > 15*15) or (abs(event.gen_tb_w_decay_id) == 13    and abs(event.gen_t_w_decay_id) > 15*15):
        process = 'mutau'
        is_mutau = True
    elif (abs(event.gen_t_w_decay_id) == 13*15 and abs(event.gen_tb_w_decay_id) > 15*15) or (abs(event.gen_tb_w_decay_id) == 13*15 and abs(event.gen_t_w_decay_id) > 15*15):
        process = 'taumutauh'
    elif (abs(event.gen_t_w_decay_id) == 13    and abs(event.gen_tb_w_decay_id) == 1)    or (abs(event.gen_tb_w_decay_id) == 13    and abs(event.gen_t_w_decay_id) == 1):
        process = 'muj'
    elif (abs(event.gen_t_w_decay_id) == 15*13 and abs(event.gen_tb_w_decay_id) == 1)    or (abs(event.gen_tb_w_decay_id) == 15*13 and abs(event.gen_t_w_decay_id) == 1):
        process = 'taumuj'

    else:
        process = 'other'

    count_gen_tw, count_gen_tbw, all_histo, cut_histo = histos[process]

    count_gen_tw  .Fill(event.gen_t_w_decay_id)
    count_gen_tbw .Fill(event.gen_tb_w_decay_id)

    # the nominal event weight
    # consists of only PU (no amcatnlo, z or recoil for ttbar)
    if is_electron_process:
        the_mc_weight = pileup_ratio_ele[event.nvtx_gen]
        #weight_pu_el_up = pileup_ratio_up_ele  [ev.nvtx_gen]
        #weight_pu_el_dn = pileup_ratio_down_ele[ev.nvtx_gen]
    else:
        the_mc_weight = pileup_ratio[event.nvtx_gen]
        #weight_pu_mu_up = pileup_ratio_up[ev.nvtx_gen]
        #weight_pu_mu_dn = pileup_ratio_down[ev.nvtx_gen]

    # get systematic weights
    weights_gen_weight_nom   = event.gen_weights_renorm_fact[MUf_nom_MUr_nom] if event.gen_weights_renorm_fact[MUf_nom_MUr_nom] > 0. else 0.00001
    # the nominal is found to be always == 1
    weights_gen_weight_f_rUp = event.gen_weights_renorm_fact[MUf_nom_MUr_up]
    weights_gen_weight_f_rDn = event.gen_weights_renorm_fact[MUf_nom_MUr_down]
    weights_gen_weight_fUp_r = event.gen_weights_renorm_fact[MUf_up_MUr_nom]
    weights_gen_weight_fDn_r = event.gen_weights_renorm_fact[MUf_down_MUr_nom]
    weights_gen_weight_frUp  = event.gen_weights_renorm_fact[MUf_up_MUr_up]
    weights_gen_weight_frDn  = event.gen_weights_renorm_fact[MUf_down_MUr_down]
    w_scale = [weights_gen_weight_f_rUp, weights_gen_weight_f_rDn, weights_gen_weight_fUp_r, weights_gen_weight_fDn_r, weights_gen_weight_frUp, weights_gen_weight_frDn]

    # pdf CT14 nominal
    weights_gen_weight_norm = event.gen_weights_pdf_hessians[0]
    # norm is the nominal PDF
    if 0. <= weights_gen_weight_norm < 0.00001:
        weights_gen_weight_norm = 0.00001
    elif 0. > weights_gen_weight_norm > -0.00001:
        weights_gen_weight_norm = -0.00001

    w_pdf = [round_pdf_weight(event.gen_weights_pdf_hessians[i] / weights_gen_weight_norm) for i in range(1,57)]

    # alphaS
    #weights_gen_weight_alphas = (event.gen_weight_alphas_1, event..gen_weight_alphas_2)
    weights_gen_weight_alphas_up = event.gen_weight_alphas_1  / weights_gen_weight_norm
    weights_gen_weight_alphas_dn = event.gen_weight_alphas_2  / weights_gen_weight_norm

    # sane high weight in tt->taumu tauh as in pdfs
    if weights_gen_weight_alphas_up > 2.:
        weights_gen_weight_alphas_up = 1.
    if weights_gen_weight_alphas_dn > 2.:
        weights_gen_weight_alphas_dn = 1.

    w_alphas = [weights_gen_weight_alphas_up, weights_gen_weight_alphas_dn]

    # save the overall number of events
    # and events passing cuts
    # events must be weighted with their weights
    for sys_i, sys_w in enumerate([the_mc_weight] + w_scale + w_alphas + w_pdf):
        all_histo.Fill(sys_i, sys_w)

    # there is always 1 lepton and 1 tau in signal:
    if not len(event.gen2_leptons_p4) > 0: continue
    if not args.no_tau_cut and len(event.gen_tt_tau_vis_p4) > 0: continue

    # lepton cuts
    if       is_electron_process and not (event.gen2_leptons_p4[0].pt() > 30 and abs(event.gen2_leptons_p4[0].eta()) < 2.4 and (abs(event.gen2_leptons_p4[0].eta()) < 1.442 or abs(event.gen2_leptons_p4[0].eta()) > 1.566)): continue
    elif not is_electron_process and not (event.gen2_leptons_p4[0].pt() > 26 and abs(event.gen2_leptons_p4[0].eta()) < 2.4): continue

    control_histos['pre_gen_lep_pt'].Fill(event.gen2_leptons_p4[0].pt())
    if is_eltau:
        control_histos['pre_gen_el_pt'].Fill(event.gen2_leptons_p4[0].pt())
    elif is_mutau:
        control_histos['pre_gen_mu_pt'].Fill(event.gen2_leptons_p4[0].pt())

    # tau cuts
    if not args.no_tau_cut and not (event.gen_tt_tau_vis_p4[0].pt() > 30 and abs(event.gen_tt_tau_vis_p4[0].eta()) < 2.4): continue

    # jet cuts
    n_jets_pass = 0
    n_jets_taucands_pass = 0
    n_b_jets_pass = 0
    for jet_p4, jet_pdgId in zip(event.gen2_jets_p4, event.gen2_jets_pdgId):
        if jet_p4.pt() > 30 and abs(jet_p4.eta()) < 2.5:
            n_jets_pass += 1
            if abs(jet_p4.eta()) < 2.4:
                n_jets_taucands_pass += 1
            if jet_pdgId == 5:
                n_b_jets_pass += 1

    if not n_jets_pass > 2: continue
    if not n_jets_taucands_pass > 0: continue

    # b-jet
    if not n_b_jets_pass > 0: continue

    control_histos['gen_lep_pt'].Fill(event.gen2_leptons_p4[0].pt())
    if is_eltau:
        control_histos['gen_el_pt'].Fill(event.gen2_leptons_p4[0].pt())
    elif is_mutau:
        control_histos['gen_mu_pt'].Fill(event.gen2_leptons_p4[0].pt())

    for sys_i, sys_w in enumerate([the_mc_weight] + w_scale + w_alphas + w_pdf):
        cut_histo.Fill(sys_i, sys_w)






# make job dir if needed
if not os.path.isdir(args.output_dir):
    try:
        os.makedirs(args.output_dir)

    except OSError as e:
        print e.errno, e.strerror

fout = TFile(output_path, "RECREATE")
fout.Write()

fout.cd()

for count1, count2, out_histo1, out_histo2 in histos.values():
    count1.SetDirectory(fout)
    count1.Write()
    count2.SetDirectory(fout)
    count2.Write()

    out_histo1.SetDirectory(fout)
    out_histo1.Write()
    out_histo2.SetDirectory(fout)
    out_histo2.Write()

for out_histo in control_histos.values():
    out_histo.SetDirectory(fout)
    out_histo.Write()

fout.Write()
fout.Close()

