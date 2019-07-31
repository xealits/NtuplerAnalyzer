from ROOT import *
import sys

gROOT.SetBatch(True)

#  KEY: TH1F	h
#  KEY: TH1F	h_0
#  KEY: TH1F	dyjets
#  KEY: TH1F	qcd
#  KEY: TH1F	single top
#  KEY: TH1F	tt_lj
#  KEY: TH1F	tt_{\mu\tau}
#  KEY: TH1F	tt_{other}
#  KEY: TH1F	w-jets

# MODEL FACTORS
LUMI_FACTOR      = RooRealVar("LUMI_FACTOR",      "The lumi factor",  1.0, 0.9, 1.1)
lf = LUMI_FACTOR
tt_factor        = RooRealVar("tt_factor",        "The tt xsec factor", 1.0, 0.5, 1.5)
#tauID_factor     = RooRealVar("tauID_factor",     "The tt tau ID SF",   0.97, 0.85, 1.1)
tauID_factor     = RooRealVar("tauID_factor",     "The tt tau ID SF",   0.91, 0.80, 1.1) # new, from Tau POG, from ttbar events
tt_taufakefactor = RooRealVar("tt_taufakefactor", "The tt tau fake SF", 1.0, 0.8, 1.2)

pu_shift = RooRealVar("pu_shift", "The PU reweighting shift", 0.0, -1.0, 1.0)
top_pt_shift = RooRealVar("top_pt_shift", "The Top pT shift", 0.0, -1.0, 1.0)
tau_es_shift = RooRealVar("tau_es_shift", "The tau ES shift", 0.0, -1.0, 1.0)

par_0    = RooFormulaVar('shape_SYS_0', '(1 - abs(@0) - @1 - @2)', RooArgList(pu_shift, top_pt_shift, tau_es_shift))
par_pu_1 = RooFormulaVar('shape_SYS_pu1', '(@0 > 0 ? @0 :  0)', RooArgList(pu_shift))
par_pu_2 = RooFormulaVar('shape_SYS_pu2', '(@0 > 0 ?  0 :-@0)', RooArgList(pu_shift))
par_top_pt = RooFormulaVar('shape_SYS_top_pt', '(@0)', RooArgList(top_pt_shift))
par_tau_es = RooFormulaVar('shape_SYS_tau_es', '(@0)', RooArgList(tau_es_shift))

xsec_dyjets_factor = RooRealVar("xsec_dyjets_factor", "The DY x-sec factor",    1.0, 0.7, 1.3)
xsec_wjets_factor  = RooRealVar("xsec_wjets_factor",  "The WJets x-sec factor", 1.0, 0.7, 1.3)
xsec_singletop_factor  = RooRealVar("xsec_singletop_factor",  "The Single Top x-sec factor", 1.0, 0.7, 1.3)
#model_factors = 

# constraints (added externally at the fit)
# https://root.cern.ch/root/html530/tutorials/roofit/rf604_constraints.C.html
lumi_constr  = RooGaussian("lumi_constr",  "lumi_constr",  LUMI_FACTOR, RooFit.RooConst(1.0), RooFit.RooConst(0.026))
#tauID_constr   = RooGaussian("tauID_constr", "tauID_constr", tauID_factor, RooFit.RooConst(0.97), RooFit.RooConst(0.05))
tauID_constr   = RooGaussian("tauID_constr", "tauID_constr", tauID_factor, RooFit.RooConst(0.91), RooFit.RooConst(0.07)) # the new
tauFake_constr = RooGaussian("tauFake_constr", "tauFake_constr", tt_taufakefactor, RooFit.RooConst(1.0), RooFit.RooConst(0.1))

# the shifts are supposedly done by 1 sigma
pu_shift_constr = RooGaussian("pu_shift_constr", "pu_shift_constr", pu_shift, RooFit.RooConst(0.0), RooFit.RooConst(1.0))
# pu shift without constrain:
#RooRealVar::pu_shift = -2.39713 +/- 2.39713  L(-5 - 5) 
# -- correlates with PU troubles?
top_pt_shift_constr = RooGaussian("top_pt_shift_constr", "top_pt_shift_constr", top_pt_shift, RooFit.RooConst(0.0), RooFit.RooConst(2.0))
tau_es_shift_constr = RooGaussian("tau_es_shift_constr", "tau_es_shift_constr", tau_es_shift, RooFit.RooConst(0.0), RooFit.RooConst(1.0))


xsec_dyjets_constr = RooGaussian("xsec_dyjets_constr", "xsec_dyjets_constr", xsec_dyjets_factor, RooFit.RooConst(1.0), RooFit.RooConst(0.3))
xsec_wjets_constr  = RooGaussian("xsec_wjets_constr",  "xsec_wjets_constr",  xsec_wjets_factor , RooFit.RooConst(1.0), RooFit.RooConst(0.3))
xsec_singletop_constr  = RooGaussian("xsec_singletop_constr",  "xsec_singletop_constr",  xsec_singletop_factor , RooFit.RooConst(1.0), RooFit.RooConst(0.1))

#contraints_set = RooArgSet(lumi_constr, tauID_constr, tauFake_constr, xsec_dyjets_constr, xsec_wjets_constr, xsec_singletop_constr)
contraints_set = RooArgSet(lumi_constr, tauID_constr, tauFake_constr, pu_shift_constr, top_pt_shift_constr, tau_es_shift_constr,
    xsec_dyjets_constr, xsec_wjets_constr, xsec_singletop_constr)


LUMI_FACTOR = lf
# CHANNELS
data_names = ("data"      , ('h', 'h_d'))
channels = {
    "dyjets"   : ("dyjets"       , (0.7, 1.3, "(@0*@1*@2)",       lambda par: RooArgList(par, lf, xsec_dyjets_factor)             ) ),
    #"dyjets"   : ("dyjets"       , (0.5, 1.5, "(@0*@1)",       lambda par: RooArgList(par, lf)             ) ),
    "wjets"    : ("w-jets"       , (0.7, 1.3, "(@0*@1*@2)",       lambda par: RooArgList(par, lf, xsec_wjets_factor)                     ) ),
    #"wjets"    : ("w-jets"       , (0.5, 1.5, "(@0*@1)",       lambda par: RooArgList(par, lf)                     ) ),
    #"qcd"      : ("qcd"          , (0.1, 1.9, "(@0*@1)",       lambda par: RooArgList(par, lf)                                ) ),
    "singletop": ("single top"   , (0.7, 1.3, "(@0*@1*@2*@3)",    lambda par: RooArgList(par, lf, tauID_factor, xsec_singletop_factor)   ) ),
    #"singletop": ("single top"   , (0.7, 1.3, "(@0*@1*@2)",    lambda par: RooArgList(par, lf, tauID_factor)   ) ),
    "tt_lj"    : ("tt_lj"        , (0.7, 1.3, "(@0*@1*@2*@3)", lambda par: RooArgList(par, lf, tt_factor, tt_taufakefactor)   ) ),
    "tt_mutau" : ("tt_{\mu\\tau}", (0.7, 1.3, "(@0*@1*@2*@3)", lambda par: RooArgList(par, lf, tt_factor, tauID_factor)       ) ),
    #"tt_eltau" : ("tt_{el\\tau}",  (0.7, 1.3, "(@0*@1*@2*@3)", lambda par: RooArgList(par, lf, tt_factor, tauID_factor)       ) ),
    "tt_eltau" : ("tt_{e\\tau}", (0.7, 1.3, "(@0*@1*@2*@3)", lambda par: RooArgList(par, lf, tt_factor, tauID_factor)       ) ),
    "tt_other" : ("tt_{other}"   , (0.7, 1.3, "(@0*@1*@2)",    lambda par: RooArgList(par, lf, tt_factor)                     ) )}
channel_names = channels.keys() # for consistency later, when creating vectors of ci*pdfi

LUMI_FACTOR.Print()

# CATEGORIES

#filename = "QuickNtupleDistr_sel-mu_sel1lep3jets40met1b1t_allWeights_mu-pt.root"
#category_files = [('incl', 'QuickNtupleDistr_sel-mu_sel1lep3jets40met1b1t_allWeights_MET-lep-Mt.root'),
#    ('second', 'QuickNtupleDistr_sel-mu_sel1lep3jets40met1b1t_allWeights_MET-lep-Mt.root')]

#category_files = [
#    ('muDparameter', 'QuickNtupleDistr_deep-distr_sel-mu_sel1lep3jets40met1b1t_allWeights_D-without-tau_to1_split_fixed_0-800-n10.root'),
#    ('tauMETouLJ', 'QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b_allWeights_tau-MET-MT_passTauOutLJ_WaMC.root')]

#category_files = [
#    ('el1', 'QuickNtupleDistr_deep-distr_sel-el_sel1lep3jets40met1b1t_allWeights_D-without-tau_to1_split_fixed_0-800-n10.root'),
#    ('el2', 'QuickNtupleDistr_deep-distr_sel-el_sel1lep3jets40met1b1t_allWeights_D-without-tau_to1_split_fixed_0-800-n10.root')]

#category_files = [
#    ('mu1', 'QuickNtupleDistr_deep-distr_sel-mu_sel1lep3jets40met1b1t_allWeights_D-without-tau_to1_split_fixed_0-800-n10.root'),
#    ('mu2', 'QuickNtupleDistr_deep-distr_sel-mu_sel1lep3jets40met1b1t_allWeights_D-without-tau_to1_split_fixed_0-800-n10.root')]

#category_files = [
#    ('el', 'QuickNtupleDistr_deep-distr_sel-el_sel1lep3jets40met1b1t_allWeights_D-without-tau_to1_split_fixed_0-800-n10.root'),
#    ('mu', 'QuickNtupleDistr_deep-distr_sel-mu_sel1lep3jets40met1b1t_allWeights_D-without-tau_to1_split_fixed_0-800-n10.root')]

#category_files = [
#    ('mu-MET', 'QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_allWeights_mu-MET-MT_All_WaMC.root'),
#    ('el-MET', 'QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_allWeights_el-MET-MT_All_WaMC.root')]
#    #('el', 'QuickNtupleDistr_deep-distr_sel-el_sel1lep3jets40met1b1t_allWeights_D-without-tau_to1_split_fixed_0-800-n10.root')]

#categories_name, category_files = "el", [
#    ('el-inLJ',   'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_allWeights_el-MET-MT_All_inLJ1600-Scale0984_WaMC.root'),
#    ('el-out-3j', 'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_allWeights_el-MET-MT_All_outLJ1600-3j-Scale0984_WaMC.root'),
#    ('el-out-4j', 'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_allWeights_el-MET-MT_All_outLJ1600-4j-Scale0984_WaMC.root')]
#
#categories_name, category_files = "mu", [
#    ('mu-inLJ',   'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_allWeights_mu-MET-MT_inLJ1600-Scale-1036_WaMC.root'),
#    ('mu-out-3j', 'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_allWeights_mu-MET-MT_outLJ1600-3j-Scale-1036_WaMC.root'),
#    ('mu-out-4j', 'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_allWeights_mu-MET-MT_outLJ1600-4j-Scale-1036_WaMC.root')]
#
#categories_name, category_files = 'combined', [
#    ('el-inLJ',   'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_allWeights_el-MET-MT_All_inLJ1600-Scale0984_WaMC.root'),
#    ('el-out-3j', 'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_allWeights_el-MET-MT_All_outLJ1600-3j-Scale0984_WaMC.root'),
#    ('el-out-4j', 'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_allWeights_el-MET-MT_All_outLJ1600-4j-Scale0984_WaMC.root'),
#    ('mu-inLJ',   'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_allWeights_mu-MET-MT_inLJ1600-Scale-1036_WaMC.root'),
#    ('mu-out-3j', 'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_allWeights_mu-MET-MT_outLJ1600-3j-Scale-1036_WaMC.root'),
#    ('mu-out-4j', 'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_allWeights_mu-MET-MT_outLJ1600-4j-Scale-1036_WaMC.root')]


categories_name, category_files = 'mu', [
    ('mu-inLJ',   'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_inLJ-Scale-1.036_WaMC.root'),
    ('mu-out-3j', 'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_outLJ-3j-Scale-1.036_WaMC.root'),
    ('mu-out-4j', 'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_outLJ-4j-Scale-1.036_WaMC.root')]

categories_name, category_files = 'el', [
    ('el-inLJ',   'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_inLJ-Scale-0.983_WaMC.root'    ),
    ('el-out-3j', 'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_outLJ-3j-Scale-0.983_WaMC.root'),
    ('el-out-4j', 'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_outLJ-4j-Scale-0.983_WaMC.root')]

categories_name, category_files = 'combined', [
    ('el-inLJ',   'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_inLJ-Scale-0.983_WaMC.root'    ),
    ('el-out-3j', 'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_outLJ-3j-Scale-0.983_WaMC.root'),
    ('el-out-4j', 'histos/QuickNtupleDistr_likelihood_sel-el_sel1lep3jets40met1b1t_el-MET-MT_outLJ-4j-Scale-0.983_WaMC.root'),
    ('mu-inLJ',   'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_inLJ-Scale-1.036_WaMC.root'),
    ('mu-out-3j', 'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_outLJ-3j-Scale-1.036_WaMC.root'),
    ('mu-out-4j', 'histos/QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b1t_mu-MET-MT_outLJ-4j-Scale-1.036_WaMC.root')]






#category_files = [
#    ('cat3j1b', 'QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b_allWeights_mu-MET-MT_cat3j1b_WaMC.root'),
#    ('cat4j1b', 'QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b_allWeights_mu-MET-MT_cat4j1b_WaMC.root'),
#    ('cat3j2b', 'QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b_allWeights_mu-MET-MT_cat3j2b_WaMC.root'),
#    ('cat4j2b', 'QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b_allWeights_mu-MET-MT_cat4j2b_WaMC.root')]
    #('tauMETouLJ', 'QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b_allWeights_tau-MET-MT_passTauOutLJ_WaMC.root')]
    #('notPassTau', 'QuickNtupleDistr_likelihood_sel-mu_sel1lep3jets40met1b_allWeights_mu-MET-MT_NOTpassTau_WaMC.root')]
    #('elmu',    'QuickNtupleDistr_likelihood_sel-elmu_selElMu3jets40met1b_allWeights_lep-MET-MT_ElMu_WaMC.root')]

def get_allsys_files(filename):
    files = {'NOMINAL': {'tau': TFile(filename), 'pretau': TFile(filename.replace('.root', '_pretau.root'))}} # nominal file
    #files['PU_UP']    = {'tau': TFile(filename.replace('.root', '_PU_UP.root'))  , 'pretau': TFile(filename.replace('.root', '_pretau_PU_UP.root'))  }
    #files['PU_DOWN']  = {'tau': TFile(filename.replace('.root', '_PU_DOWN.root')), 'pretau': TFile(filename.replace('.root', '_pretau_PU_DOWN.root'))}
    #files['TOP_PT']   = {'tau': TFile(filename.replace('.root', '_TOP_PT.root')), 'pretau': TFile(filename.replace('.root', '_pretau_TOP_PT.root'))}
    #files['TAU_ES']   = {'tau': TFile(filename.replace('.root', '_TAU_ES_DOWN.root')), 'pretau': TFile(filename.replace('.root', '_pretau_TAU_ES_DOWN.root'))}
    return files

filenames = [i for _, i in category_files]
#files = [TFile(f) for f in filenames]
#files = [{'tau': get_allsys_files(f), 'pretau': get_allsys_files(f.replace(".root", "_template.root"))} for f in filenames]
files = [get_allsys_files(f) for f in filenames]
# files is list of catogries,
# each category is dictionary with systematic files,
# in each file we have bunch of channels

#pretau_files = [TFile(f.replace(".root", "_template.root")) for f in filenames]
#pretau_files = [get_allsys_files(f.replace(".root", "_template.root")) for f in filenames]

def get_data_histo(f, data_names):
    for n in data_names:
        if f.FindKey(n):
            return f.Get(n)
    raise ValueError("no appropriate data histogram name in the file")

data_histos = [(category_files[i][0], get_data_histo(f_sys['NOMINAL']['tau'], data_names[1])) for i, f_sys in enumerate(files)]
histos      = [{SYS: {name: f['tau'].Get(c) for name, (c, _) in channels.items() if f['tau'].FindKey(c)}
                     for SYS, f in f_sys.items()}
                         for f_sys in files]
# histos now are a list of categories,
# each category is dictionary of systematics...

template_histos = [{SYS: {name: f['pretau' if name in ('wjets', 'dyjets', 'qcd') else 'tau'].Get(c) for name, (c, _) in channels.items() if f['tau'].FindKey(c) or f['pretau'].FindKey(c)}
                     for SYS, f in f_sys.items()}
                         for f_sys in files]
# template histos are needed for templates, histos -- for normalization integrals
# poor statistics in some MC requires using pre-tau templates here

for c in histos:
    print c.keys()

print 'FOOOOOOOOOOOOOOOOOO', channels['tt_lj']

# the N PARAMETERS for the fit
initial_n_events = [{SYS: {chan_name: chan_h.Integral() for chan_name, chan_h in channels.items()} for SYS, channels in category.items()}
    for category in histos] # thus it corresponds to MC xsec and nominal lumi

data_int_el, mc_int_el = 0, 0
data_int_mu, mc_int_mu = 0, 0
for (cat, data_histo), mc_cat in zip(data_histos, initial_n_events):
    data_int, mc_int = data_histo.Integral(), sum(n for _, n in mc_cat['NOMINAL'].items())
    print cat, data_int, mc_int, mc_int/data_int
    if 'el' in cat:
        data_int_el += data_int
        mc_int_el += mc_int
    else:
        data_int_mu += data_int
        mc_int_mu += mc_int

print 'el', data_int_el, mc_int_el
print 'mu', data_int_mu, mc_int_mu

el_scale, mu_scale = 0.986, 1.036
scales_inclusive = {
 True:  (el_scale, data_int_el, mc_int_el),
 False: (mu_scale, data_int_mu, mc_int_mu)}

nominal_scales_per_channel = []

for i, mc_cat in enumerate(initial_n_events):
    cat_name = category_files[i][0]
    scale, data_int, mc_int = scales_inclusive['el' in cat_name]

    nominal_scales_per_channel = {}
    for chan_name, n_ev in mc_cat['NOMINAL'].items():
        nominal_scales_per_channel[chan_name] =  (1 / mc_int) * data_int * scale

    for SYS, ch_s in mc_cat.items():
        for chan_name, _ in ch_s.items():
            #mc_cat[SYS][chan_name] = (n_ev / mc_int) * data_int * scale
            mc_cat[SYS][chan_name] *= nominal_scales_per_channel[chan_name]

for (cat, data_histo), mc_cat in zip(data_histos, initial_n_events):
    data_int, mc_int = data_histo.Integral(), sum(n for _, n in mc_cat['NOMINAL'].items())
    print cat, data_int, mc_int, mc_int/data_int

print 'FOOOOOOOOOOOOOOOOOO', channels['tt_lj']

#n_events_parameters = [{chn: RooRealVar('n' + chn, 'n' + chn, chan_n, chan_n*channels[chn][1][0], chan_n*channels[chn][1][1]) for chn, chan_n in category.items()}
#n_events_parameters = [{chn: RooConstVar('n' + chn, 'n' + chn, chan_n, chan_n*channels[chn][1][0], chan_n*channels[chn][1][1]) for chn, chan_n in category.items()}
n_events_parameters = [{SYS: {chn: RooConstVar('n_%s_%s' % (SYS, chn), 'n' + chn, chan_n) for chn, chan_n in channels.items()} for SYS, channels in category.items()}
    for category in initial_n_events]

neventsparameters = n_events_parameters

for i, c in enumerate(initial_n_events):
    print category_files[i][0]
    for n, N in c.items():
        print n, N

LUMI_FACTOR.Print()

#sys.exit(1)
#    "tt_other" : ("tt_{other}"   , (0.7, 1.3, "(@0*@1*@2)",    lambda par: RooArgList(par, LUMI_FACTOR, tt_factor)                     ) )}
# for each category, each channel should have a model of N events for systematic shifts,
# in the following it's multiplied by factors
# N events for systematic shifts are set independent from each other, as a quadratic or linear factor

def get_sys_N_events(category):
    '''
    category is dictionary of SYS shifts with channels inside
    category['NOMINAL'] are the central values

    returns dictionary of {channel_name: RooFormulaVar()} with the model for systematic shifts on this channel
    '''
    channels_models = dict()
    for chan_name, chan_n in category['NOMINAL'].items():
        #chan_shift_model = RooFormulaVar('NSYS_' + chan_name, '((@2 + @3 - 2*@1)/2 *pow(@4 - (), 2))', )
        #chan_model = RooFormulaVar('NSYS_' + chan_name, '(@1+())', )
        # parabola is complex -- linear interpolation should be fine
        if 'PU_UP' in category and 'PU_DOWN' in category:
            #chan_model = RooFormulaVar('NSYS_' + chan_name, '(@1 + @0 * (@0 > 0? @2 - @1 : @1 - @3))', RooArgList(pu_shift, chan_n, category['PU_UP'][chan_name], category['PU_DOWN'][chan_name]))
            sys_params = RooArgList(category['NOMINAL'][chan_name], par_0,
                                                  category['PU_UP']  [chan_name], par_pu_1,
                                                  category['PU_DOWN'][chan_name], par_pu_2,
                                                  category['TOP_PT'] [chan_name], par_top_pt)
            #sys_params.add(RooArgList(category['TAU_ES'] [chan_name], par_tau_es))
            chan_model = RooFormulaVar('NSYS_' + chan_name,
                                       '(@0*@1 + @2*@3 + @4*@5 + @6*@7)', sys_params)
                                       #'(@0*@1 + @2*@3 + @4*@5 + @6*@7 + @8*@9)', sys_params)
        else:
            chan_model = chan_n # nominal constant
        print chan_model
        channels_models[chan_name] = chan_model

    return channels_models
    

n_events_parameters = neventsparameters

N_parameters_sys = [get_sys_N_events(category) for category in n_events_parameters]

print 'FOOOOOOOOOOOOOOOOOO', channels['tt_lj']

LUMI_FACTOR = lf

N_parameters_fit = []
for cat in N_parameters_sys:
    params = {}
    for chn, p in cat.items():
        LUMI_FACTOR = lf
        print chn, p
        print channels[chn]
        params[chn] = RooFormulaVar('N_' + chn, channels[chn][1][2], channels[chn][1][3](p))
    N_parameters_fit.append(params)

# final parameters for PDFs, with factors included
#N_parameters_fit = [{chn: RooFormulaVar('N_' + chn,
#                                        channels[chn][1][2],
#                                        channels[chn][1][3](chan_par))
#                     for chn, chan_par in category.items()}
#    for category in N_parameters_sys]

#N_dyjets    = RooFormulaVar("N_dyjets"   , "(@0*@1)",        RooArgList(ndyjets   , LUMI_FACTOR)                             )
#N_qcd       = RooFormulaVar("N_qcd"      , "(@0*@1)",        RooArgList(nqcd      , LUMI_FACTOR)                             )
#N_singletop = RooFormulaVar("N_singletop", "(@0*@1)",        RooArgList(nsingletop, LUMI_FACTOR)                             )
#N_wjets     = RooFormulaVar("N_wjets"    , "(@0*@1)",        RooArgList(nwjets    , LUMI_FACTOR)                             )
#N_tt_lj     = RooFormulaVar("N_tt_lj"    , "(@0*@1*@2*@3)",  RooArgList(ntt_lj,     LUMI_FACTOR, tt_factor, tt_taufakefactor))
#N_tt_mutau  = RooFormulaVar("N_tt_mutau" , "(@0*@1*@2*@3)",  RooArgList(ntt_mutau,  LUMI_FACTOR, tt_factor, tauID_factor)    )
#N_tt_other  = RooFormulaVar("N_tt_other" , "(@0*@1*@2)",     RooArgList(ntt_other,  LUMI_FACTOR, tt_factor)                  )

##dy_xsec = ROOT.RooRealVar("tt_xsec", "The tt x-sce",    837., 750., 900)
#n_dyjets    = dyjets   .Integral()
#n_qcd       = qcd      .Integral()
#n_singletop = singletop.Integral()
#n_wjets     = wjets    .Integral()
#n_tt_lj     = tt_lj    .Integral()
#n_tt_mutau  = tt_mutau .Integral()
#n_tt_other  = tt_other .Integral()

#ndyjets    = RooRealVar("ndyjets"   , "ndyjets"   , n_dyjets   , n_dyjets   * 0.7, n_dyjets   *1.3)
#nqcd       = RooRealVar("nqcd"      , "nqcd"      , n_qcd      , n_qcd      * 0., n_qcd       *2.)
#nsingletop = RooRealVar("nsingletop", "nsingletop", n_singletop, n_singletop* 0.7, n_singletop*1.3)
#nwjets     = RooRealVar("nwjets"    , "nwjets"    , n_wjets    , n_wjets    * 0.7, n_wjets    *1.3)
#ntt_lj     = RooRealVar("ntt_lj"    , "ntt_lj"    , n_tt_lj    , n_tt_lj    * 0.7, n_tt_lj    *1.3)
#ntt_mutau  = RooRealVar("ntt_mutau" , "ntt_mutau" , n_tt_mutau , n_tt_mutau * 0.7, n_tt_mutau *1.3)
#ntt_other  = RooRealVar("ntt_other" , "ntt_other" , n_tt_other , n_tt_other * 0.7, n_tt_other *1.3)

#h_data     = f.Get("h")
#dyjets     = f.Get("dyjets")
#qcd        = f.Get("qcd")
#singletop  = f.Get("single top")
#tt_lj      = f.Get("tt_lj")
#tt_mutau   = f.Get("tt_{\mu\\tau}")
#tt_other   = f.Get("tt_{other}")
#wjets      = f.Get("w-jets")

# check 1st category histograms only
for o in histos[0]:
    print type(o)


#PDFS

LUMI_FACTOR = lf
LUMI_FACTOR.Print()

# Create variable sets
#mupt = RooRealVar("mupt", "#mu p_{T}", 0., 200., "GeV")
distr = RooRealVar("distr", "l E^{miss}_{T} M_{T}", 0., 250., "GeV")
#distr = RooRealVar("distr", "D", 0., 800., "GeV")
#mass = ROOT.RooRealVar("mass", "#mu^{+}#mu^{-} invariant mass", 2., 6., "GeV")
vars = RooArgList()
vars.add(distr)
vars_set = RooArgSet()
vars_set.add(distr)

# Create RooDataHists
# out of these
#data_histos = [f.Get(data_names[1]) for f in files]
#histos      = [{name: f.Get(c) for name, (c, _) in channels.items()} for f in files]
data_rh   = [RooDataHist("data_" + n, "dataset with distr", vars, h_data) for n, h_data in data_histos]
histos_rh = [{SYS: {name: RooDataHist('rh_' + name, name, vars, h) for name, h in channels.items()} for SYS, channels in category.items()}
             for category in template_histos]

#rh_dyjets     = RooDataHist("rh_dyjets",    "dyjets",    vars, dyjets   )
#rh_qcd        = RooDataHist("rh_qcd",       "qcd",       vars, qcd      )
#rh_singletop  = RooDataHist("rh_singletop", "singletop", vars, singletop)
#rh_tt_lj      = RooDataHist("rh_tt_lj",     "tt_lj",     vars, tt_lj    )
#rh_tt_mutau   = RooDataHist("rh_tt_mutau",  "tt_mutau",  vars, tt_mutau )
#rh_tt_other   = RooDataHist("rh_tt_other",  "tt_other",  vars, tt_other )
#rh_wjets      = RooDataHist("rh_wjets",     "wjets",     vars, wjets    )

# Create PDFs
#pdfs = [{name: RooHistPdf('pdf_'+name, 'pdf '+name, vars_set, h, 0) for name, h in category.items()} in histos_rh]
# root complains about argument 4, the passed histogram
#  RooHistPdf::RooHistPdf(const char* name, const char* title, const RooArgSet& vars, const RooDataHist& dhist, int intOrder = 0) =>
#    could not convert argument 4
# it is of class  <class '__main__.RooDataHist'>
#
# "combined fit" example uses this (in c++):
# RooAddPdf  model("model","g1+g2+a",RooArgList(bkg,sig),bkgfrac) 
# but here the error doesn't show such method

pdfs = []
#check this
for c in histos_rh:
    pdfs.append(dict())
    for SYS, channels in c.items():
        #pdf[-1][SYS] = {n} # one-liner went into some troubles with ROOT
        pdfs[-1][SYS] = dict()
        for n, h in channels.items():
            pdf = RooHistPdf('pdf_%s_%s' % (SYS, n), 'pdf %s %s' % (SYS, n), vars_set, h, 0)
            print n, type(h), type(pdf)
            pdfs[-1][SYS][n] = pdf
            # doesn't break...
# again the same structure (need a class for this already):
# a list categories,
# category = dict SYS: channels
# channels = name, RooHistPdf

# pdfs for fit are list of categories,
# where each category is dictionary channel_name: Pdf model

# old pdfs for fit -- NOMINAL templates
nominal_pdfs = [cat['NOMINAL'] for cat in pdfs]

# new PDF models, systematic shapes, linearly interpolated
sys_pdfs = []
for category in pdfs:
    channels_models = dict()
    for chan_name, chan_nominal in category['NOMINAL'].items():
        if 'PU_UP' in category and 'PU_DOWN' in category:
            # trying creating the heaviside coeficients on the fly and adding all shift pdfs with this coefficients
            chan_model = RooAddPdf('pdf_SYS_' + chan_name, '(1 - |e|)*S0 + e*S1 - eS-1',
                                   RooArgList(category['NOMINAL'][chan_name], category['PU_UP'][chan_name], category['PU_DOWN'][chan_name], category['TOP_PT'][chan_name]),
                                   RooArgList(par_0, par_pu_1, par_pu_2, par_top_pt))
                                   #RooArgList(category['NOMINAL'][chan_name], category['PU_UP'][chan_name], category['PU_DOWN'][chan_name], category['TOP_PT'][chan_name], category['TAU_ES'][chan_name]),
                                   #RooArgList(par_0, par_pu_1, par_pu_2, par_top_pt, par_tau_es))

            # creating coefficients on the fly ran into errors
            # the errors made no sense:
            #
            #Traceback (most recent call last):
            #  File "my_combined_histofit.py", line 410, in <module>
            #    RooFormulaVar('shape_SYS2_'+chan_name, '(@0 > 0 ?  0 :-@0)', RooArgList(pu_shift))))
            #TypeError: none of the 6 overloaded methods succeeded. Full details:
            #  RooAddPdf::RooAddPdf(const char* name, const char* title, const RooArgList& pdfList, const RooArgList& coefList, bool recursiveFraction = kFALSE) =>
            #    problem in C++; program state has been reset
            #  RooAddPdf::RooAddPdf() =>
            #    takes at most 0 arguments (4 given)
            #  RooAddPdf::RooAddPdf(const char* name, const char* title = 0) =>
            #    takes at most 2 arguments (4 given)
            #  RooAddPdf::RooAddPdf(const char* name, const char* title, RooAbsPdf& pdf1, RooAbsPdf& pdf2, RooAbsReal& coef1) =>
            #    takes at least 5 arguments (4 given)
            #  RooAddPdf::RooAddPdf(const char* name, const char* title, const RooArgList& pdfList) =>
            #    takes at most 3 arguments (4 given)
            #  RooAddPdf::RooAddPdf(const RooAddPdf& other, const char* name = 0) =>
            #    takes at most 2 arguments (4 given)
            # -- the definition itself iplies option with 4 arguments, and this option is used successfully in the final model!

            # how to do it with many systematics?
        else:
            chan_model = chan_nominal # nominal constant
        print chan_model
        channels_models[chan_name] = chan_model
    sys_pdfs.append(channels_models)


fit_pdfs, pdf_prefix = sys_pdfs, 'pdf_SYS_'
fit_pdfs, pdf_prefix = nominal_pdfs, 'pdf_NOMINAL_'




# systematic shape variations with complex AddPdf models


LUMI_FACTOR = lf
LUMI_FACTOR.Print()

#pdf_dyjets     = RooHistPdf("pdf_dyjets",     "pdf dyjets",    vars_set, rh_dyjets   , 0)
#pdf_qcd        = RooHistPdf("pdf_qcd",        "pdf qcd",       vars_set, rh_qcd      , 0)
#pdf_singletop  = RooHistPdf("pdf_singletop",  "pdf singletop", vars_set, rh_singletop, 0)
#pdf_tt_lj      = RooHistPdf("pdf_tt_lj",      "pdf tt_lj",     vars_set, rh_tt_lj    , 0)
#pdf_tt_mutau   = RooHistPdf("pdf_tt_mutau",   "pdf tt_mutau",  vars_set, rh_tt_mutau , 0)
#pdf_tt_other   = RooHistPdf("pdf_tt_other",   "pdf tt_other",  vars_set, rh_tt_other , 0)
#pdf_wjets      = RooHistPdf("pdf_wjets",      "pdf wjets",     vars_set, rh_wjets    , 0)



#h_data = FileReader.getHistogramFromFile("h", filename)
#h_datadriven_bkg = FileReader.getHistogramFromFile("myDistribution2", filename)
#h_signal = FileReader.getHistogramFromFile("myDistribution", filename)
#h_bk1 = FileReader.getHistogramFromFile("myDistribution", filename)
#h_bk2 = FileReader.getHistogramFromFile("myDistribution", filename)

#dyjets     = FileReader.getHistogramFromFile("dyjets",       filename)
#qcd        = FileReader.getHistogramFromFile("qcd",          filename)
#singletop  = FileReader.getHistogramFromFile("single top",   filename)
#tt_lj      = FileReader.getHistogramFromFile("tt_lj",        filename)
#tt_mutau   = FileReader.getHistogramFromFile("tt_{\mu\tau}", filename)
#tt_other   = FileReader.getHistogramFromFile("tt_{other}",   filename)
#wjets      = FileReader.getHistogramFromFile("w-jets",       filename)



# get normalisation
#n_event_obs = h_data.GetEntries();
#n_signal = h_signal
#n_bkg1 = h_bkg1.Integral()
#n_bkg2 = h_bkg2.Integral()
#n_datadriven_bkg = 800#data-driven estimate

#lumi    = ROOT.RooRealVar("lumi",    "The recorded luminosity pb", 35200., 35200. * 0.97, 35200 * 1.03)
#tt_xsec = ROOT.RooRealVar("tt_xsec", "The tt x-sce",    837., 750., 900)


## Create fit variables
#lowerBound = -10 * sqrt(n_event_obs); 
#upperBound = n_event_obs + 10 * sqrt(n_event_obs);
#nSignal = RooRealVar  ("nSignal", "number of signal events", n_signal, lowerBound, upperBound, "event");
#nBkg1   = RooRealVar  ("nbkg1", "number of bkg1 events", n_bkg1, lowerBound, upperBound, "event");
#nBkg2   = RooRealVar  ("nbkg2", "number of bkg1 events", n_bkg2, lowerBound, upperBound, "event");
#ndatadriven_bkg = RooRealVar("ndatadriven_bkg", "number of datadriven_bkg events", n_datadriven_bkg, lowerBound, upperBound, "event");


# 0 at the end is "intOrder" https://root.cern.ch/doc/v608/classRooHistPdf.html



## Create test model of 1st category
#model1 = RooAddPdf("model1", "sig+bkg1+bkg2+datadriven_bkg",
#    RooArgList(*[pdfs[0][n] for n in channel_names]),
#    RooArgList(*[N_parameters_fit[0][n] for n in channel_names]))
#ifitResult = model1.fitTo(data_rh[0], RooFit.SumW2Error(kTRUE), RooFit.PrintEvalErrors(-1))




# Create combined model of all categories
# Define category to distinguish physics and control samples events

LUMI_FACTOR.Print()
print LUMI_FACTOR, type(LUMI_FACTOR)
lf = LUMI_FACTOR
sample = RooCategory("sample", "sample") 

print lf
LUMI_FACTOR = lf
LUMI_FACTOR.Print()

#sample = RooCategory(*("sample" for _ in category_files))
for cat_name, _ in category_files:
    sample.defineType(cat_name) ;
#sample.defineType("control") ;

LUMI_FACTOR.Print()

# Construct combined dataset in (x,sample)
#RooDataSet combData("combData","combined data",x,Index(sample),Import("physics",*data),Import("control",*data_ctl))
# x is observable
#RooDataSet combData("combData", "combined data", distr, Index(sample), Import("physics",*data), Import("control",*data_ctl))
#combData  = RooDataSet("combData", "combined data", distr, RooFit.Index(sample), *(RooFit.Import(cat, d) for cat, d in zip([i for i, _  in category_files], data_rh)))

LUMI_FACTOR.Print()

print "FOOOOOOOOOOOOOOOOOOOOOO"
x = RooFit.Import(category_files[0][0], data_rh[0])
print x, type(x)
y = RooFit.Index(sample)
print y, type(y)

# https://root.cern.ch/root/html/tutorials/roofit/rf401_importttreethx.C.html
#// Create a binned dataset that imports contents of all TH1 mapped by index category c
#  RooDataHist* dh = new RooDataHist("dh","dh",x,Index(c),Import("SampleA",*hh_1),Import("SampleB",*hh_2),Import("SampleC",*hh_3)) ;
#  dh->Print() ;

#data_hist_combined = RooDataHist('dh', 'dh', vars, RooFit.Index(sample),
#    RooFit.Import(category_files[0][0], data_rh[0]),
#    RooFit.Import(category_files[1][0], data_rh[1]))
#
#data_hist_combined = RooDataHist('dh', 'dh', vars, RooFit.Index(sample),
#    RooFit.Import(category_files[0][0], data_rh[0]))#,
#    #RooFit.Import(category_files[1][0], data_rh[1]))

data_hist_combined = RooDataHist('dh', 'dh', vars, RooFit.Index(sample),
    *(RooFit.Import(category_files[i][0], data_rh[i]) for i, _ in enumerate(category_files)))

#combData  = RooDataSet("combData", "combined data", vars_set, RooFit.Index(sample), RooFit.Import(category_files[0][0], data_rh[0]))
#combData  = RooDataSet("combData", "combined data", vars_set, y)
#combData  = RooDataSet("combData", "combined data", vars_set, RooFit.Index(sample),
#    RooFit.Import(category_files[0][0], data_rh[0]),
#    RooFit.Import(category_files[1][0], data_rh[1]))

#// Construct a simultaneous pdf using category sample as index
#RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
#
#// Associate model with the physics state and model_ctl with the control state
#simPdf.addPdf(model,"physics") ;
#simPdf.addPdf(model_ctl,"control")

# construct models for each category
print 'constructing models'
models = [(cat, RooAddPdf("model_"+cat, "sig+bkg1+bkg2+datadriven_bkg",
    RooArgList(*[fit_pdfs[i][n] for n in channel_names if n in fit_pdfs[i]]),
    RooArgList(*[N_parameters_fit[i][n] for n in channel_names if n in N_parameters_fit[i]]))) for i, (cat, _) in enumerate(category_files)]

simModel = RooSimultaneous("simModel", "simultaneous pdf", sample)

for cat, m in models:
    simModel.addPdf(m, cat) ;

LUMI_FACTOR.Print()

print 'fitting'

ifitResult = simModel.fitTo(data_hist_combined, RooFit.ExternalConstraints(contraints_set), RooFit.SumW2Error(kTRUE))

LUMI_FACTOR.Print()

##get results
#nSignal_fit = nSignal.getVal();
#nbkg1_fit = nwj.getVal();
#nbkg2_fit = nzj.getVal();
#ndatadriven_bkg_fit = nqcd.getVal();
#nSignal_fiterr = nSignal.getError();
#nbkg1_fiterr = nbkg1.getError();
#nbkg2_fiterr = nbkg2.getError();
#ndatadriven_bkg_fiterr = ndatadriven_bkg.getError();

# make the plots with fit templates
c_templates = TCanvas("templates", "templates", 10, 10, 1200, 800)
c_templates.Divide(3,2)

templates_frames = []
for i, (name, _) in enumerate(category_files):
    cur_frame = distr.frame()
    cur_frame.SetTitle(name)
    templates_frames.append(cur_frame)
    #data_rh[i].plotOn(cur_frame)
    if 'tt_mutau' in fit_pdfs[i]: fit_pdfs [i]['tt_mutau'] .plotOn(cur_frame, RooFit.LineColor(kOrange + 1))
    else:  fit_pdfs [i]['tt_eltau'] .plotOn(cur_frame, RooFit.LineColor(kOrange+1))
    fit_pdfs [i]['tt_lj']    .plotOn(cur_frame, RooFit.LineStyle(kSolid), RooFit.LineColor(kGreen+3))
    fit_pdfs [i]['singletop'].plotOn(cur_frame, RooFit.LineStyle(kSolid), RooFit.LineColor(kAzure))
    fit_pdfs [i]['wjets']    .plotOn(cur_frame, RooFit.LineStyle(kDashed), RooFit.LineColor(kRed))
    fit_pdfs [i]['dyjets']   .plotOn(cur_frame, RooFit.LineStyle(kDashed), RooFit.LineColor(kGray))
    #fit_pdfs [i]['qcd']   .plotOn(cur_frame, RooFit.LineStyle(kDashed), RooFit.LineColor(kViolet))
    c_templates.cd(i+1)
    cur_frame.Draw()

c_templates.SaveAs("my_combined_histofit_templates.png")

# make the plots with post-fit distributions
c_postfit = TCanvas("post_fit_distr", "post_fit_distr", 10, 10, 1200, 800)
c_postfit.Divide(3,2)

for i, (name, _) in enumerate(category_files):
    cur_frame = distr.frame()
    cur_frame.SetTitle(name)
    data_rh[i].plotOn(cur_frame)

    model_sig = 'tt_mutau' if 'tt_mutau' in fit_pdfs[i].keys() else "tt_eltau"
    print 'model sig', model_sig
    
    models [i][1].plotOn(cur_frame, RooFit.Components("%ssingletop" % pdf_prefix),           RooFit.LineColor(kAzure), RooFit.FillColor(kAzure), RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())
    models [i][1].plotOn(cur_frame, RooFit.Components("%ssingletop,%s" % (pdf_prefix, pdf_prefix+model_sig) ),           RooFit.LineColor(kOrange+1), RooFit.FillColor(kOrange+1), RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())
    models [i][1].plotOn(cur_frame, RooFit.Components("%ssingletop,%s,%stt_lj" % (pdf_prefix, pdf_prefix+model_sig, pdf_prefix)), RooFit.LineColor(kGreen+3),  RooFit.FillColor(kGreen+3),   RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())
    models [i][1].plotOn(cur_frame, RooFit.Components("%ssingletop,%s,%stt_lj,%swjets" % (pdf_prefix, pdf_prefix+model_sig, pdf_prefix, pdf_prefix)), RooFit.LineColor(kRed),   RooFit.FillColor(kRed),   RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())
    models [i][1].plotOn(cur_frame, RooFit.Components("%ssingletop,%s,%stt_lj,%swjets,%sdyjets" % (pdf_prefix, pdf_prefix+model_sig, pdf_prefix, pdf_prefix, pdf_prefix)), RooFit.LineColor(kGray),   RooFit.FillColor(kGray),   RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())
    models [i][1].plotOn(cur_frame, RooFit.Components("%ssingletop,%s,%stt_lj,%swjets,%sdyjets,%sqcd" % (pdf_prefix, pdf_prefix+model_sig, pdf_prefix, pdf_prefix, pdf_prefix, pdf_prefix)), RooFit.LineColor(kViolet),   RooFit.FillColor(kViolet),   RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())

    c_postfit.cd(i+1)
    cur_frame.Draw()

c_postfit.SaveAs("my_combined_histofit_postfit.png")


#Draw the results

xframe1 = distr.frame()
data_rh[0].plotOn(xframe1)
# this one is for errors of MC:
#models [0][1].plotOn(xframe1, RooFit.LineWidth(0))
#models [0][1].plotOn(xframe1)
#models [0][1].plotOn(xframe1, RooFit.Components("pdf_tt_mutau"), RooFit.LineStyle(kDashed), RooFit.LineColor(kRed))
#models [0][1].plotOn(xframe1, RooFit.Components(RooArgSet("pdf_tt_mutau", "pdf_tt_lj")), RooArgSet(RooFit.LineStyle(kDashed), RooFit.LineStyle(kSolid)), RooArgSet(RooFit.LineColor(kRed), RooFit.LineColor(kGreen)))
#models [0][1].plotOn(xframe1, RooFit.Components("pdf_tt_mutau,pdf_tt_lj"))
#    modelPdf_->plotOn(frame, RooFit::Components(compStack), RooFit::LineColor(lineColors_[N]),RooFit::FillColor(fillColors_[N]),RooFit::FillStyle(1001), RooFit::LineWidth(2), RooFit::Slice(*cat_, slice.c_str()),RooFit::DrawOption("F"),RooFit::ProjWData(*data_));
#    modelPdf_->plotOn(frame, RooFit::Components(compStack), RooFit::LineColor(lineColors_[N]),RooFit::FillColor(fillColors_[N]),RooFit::FillStyle(1001), RooFit::LineWidth(2),RooFit::DrawOption("FL"));
#    modelPdf_->plotOn(frame, RooFit::Components(compStack), RooFit::LineColor(lineColors_[N]),RooFit::FillColor(fillColors_[N]),RooFit::FillStyle(1001), RooFit::LineWidth(2));
# std::vector<int> fillColors_;
# https://github.com/iross/UWAnalysis/blob/master/StatTools/interface/plotter.h
#models [0][1].plotOn(xframe1, RooFit.Components("pdf_tt_mutau,pdf_tt_lj"), RooFit.FillColor([kRed, kGreen]))
# he hacks:
#//now create th stack
#	for(int N=components_.size()-1;N>=0 ;N=N-1) {
#	  //create component stack
#	TString compStack;
# draws allcomp, then allcomp-1... until all is drawn

print "PDFs"
print fit_pdfs[0].keys()
print fit_pdfs[1].keys()

model0_sig = 'pdf_tt_mutau' if 'tt_mutau' in fit_pdfs[0].keys() else "pdf_tt_eltau"
model1_sig = 'pdf_tt_mutau' if 'tt_mutau' in fit_pdfs[1].keys() else "pdf_tt_eltau"

models [0][1].plotOn(xframe1, RooFit.Components("pdf_singletop"),           RooFit.LineColor(kAzure), RooFit.FillColor(kAzure), RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())
models [0][1].plotOn(xframe1, RooFit.Components("pdf_singletop," +model0_sig ),           RooFit.LineColor(kOrange+1), RooFit.FillColor(kOrange+1), RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())
models [0][1].plotOn(xframe1, RooFit.Components("pdf_singletop," +model0_sig+ ",pdf_tt_lj"), RooFit.LineColor(kGreen+3),  RooFit.FillColor(kGreen+3),   RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())
models [0][1].plotOn(xframe1, RooFit.Components("pdf_singletop," +model0_sig+ ",pdf_tt_lj,pdf_wjets"), RooFit.LineColor(kRed),   RooFit.FillColor(kRed),   RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())
models [0][1].plotOn(xframe1, RooFit.Components("pdf_singletop," +model0_sig+ ",pdf_tt_lj,pdf_wjets,pdf_dyjets"), RooFit.LineColor(kGray),   RooFit.FillColor(kGray),   RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())
models [0][1].plotOn(xframe1, RooFit.Components("pdf_singletop," +model0_sig+ ",pdf_tt_lj,pdf_wjets,pdf_dyjets,pdf_qcd"), RooFit.LineColor(kViolet),   RooFit.FillColor(kViolet),   RooFit.FillStyle(1001), RooFit.DrawOption('FL'), RooFit.MoveToBack())
#data_rh[0].plotOn(xframe1, RooFit.DrawOption('S'))
# not clear how to overlay data but this came along:
#MoveToBack

xframe2 = distr.frame()
data_rh[1].plotOn(xframe2)
models [1][1].plotOn(xframe2)
models [1][1].plotOn(xframe2, RooFit.Components(model1_sig), RooFit.LineStyle(kDashed), RooFit.LineColor(kRed))

xframe3 = distr.frame()
#data_rh[2].plotOn(xframe3)
#models [2][1].plotOn(xframe3)
#models [2][1].plotOn(xframe3, RooFit.Components("pdf_tt_mutau"), RooFit.LineStyle(kDashed), RooFit.LineColor(kRed))

if 'tt_mutau' in fit_pdfs[0]: fit_pdfs [0]['tt_mutau'] .plotOn(xframe3, RooFit.LineColor(kRed))
else:  fit_pdfs [0]['tt_eltau'] .plotOn(xframe3, RooFit.LineColor(kOrange+1))
fit_pdfs [0]['tt_lj']    .plotOn(xframe3, RooFit.LineStyle(kDashed), RooFit.LineColor(kGreen+3))
fit_pdfs [0]['singletop'].plotOn(xframe3, RooFit.LineStyle(kDashed), RooFit.LineColor(kAzure))
fit_pdfs [0]['wjets']    .plotOn(xframe3, RooFit.LineStyle(kDashed), RooFit.LineColor(kRed))
fit_pdfs [0]['dyjets']   .plotOn(xframe3, RooFit.LineStyle(kDashed), RooFit.LineColor(kGray))
#fit_pdfs [0]['qcd']   .plotOn(xframe3, RooFit.LineStyle(kDashed), RooFit.LineColor(kViolet))

xframe4 = distr.frame()
#data_rh[3].plotOn(xframe4)
#models [3][1].plotOn(xframe4)
#models [3][1].plotOn(xframe4, RooFit.Components("pdf_tt_mutau"), RooFit.LineStyle(kDashed), RooFit.LineColor(kRed))
#models [3][1].plotOn(xframe4, RooFit.Components("pdf_tt_lj"),  RooFit.LineStyle(kDashed), RooFit.LineColor(kGreen))
#models [3][1].plotOn(xframe4, RooFit.Components("pdf_wjets"),  RooFit.LineStyle(kDashed), RooFit.LineColor(kYellow))
#models [3][1].plotOn(xframe4, RooFit.Components("pdf_dyjets"), RooFit.LineStyle(kDashed), RooFit.LineColor(kOrange))
# draw the fit_pdfs, the templates
if 'tt_mutau' in fit_pdfs[1]: fit_pdfs [1]['tt_mutau'] .plotOn(xframe4, RooFit.LineColor(kRed))
else:  fit_pdfs [1]['tt_eltau'] .plotOn(xframe4, RooFit.LineColor(kRed))
fit_pdfs [1]['tt_lj']    .plotOn(xframe4, RooFit.LineStyle(kDashed), RooFit.LineColor(kGreen))
fit_pdfs [1]['wjets']    .plotOn(xframe4, RooFit.LineStyle(kDashed), RooFit.LineColor(kYellow))
fit_pdfs [1]['dyjets']   .plotOn(xframe4, RooFit.LineStyle(kDashed), RooFit.LineColor(kAzure))
fit_pdfs [1]['singletop'].plotOn(xframe4, RooFit.LineStyle(kDashed), RooFit.LineColor(kBlue))

##model.plotOn(xframe, RooFit.Components(RooArgSet(fit_pdfs[0]["dyjets"], fit_pdfs[0]["qcd"], fit_pdfs[0]["singletop"], fit_pdfs[0]["wjets"], fit_pdfs[0]["tt_lj"], fit_pdfs[0]["tt_mutau"], fit_pdfs[0]["tt_other"])))
#simModel.plotOn(xframe)
#simModel.plotOn(xframe, RooFit.Components(RooArgSet(fit_pdfs[0]["dyjets"], fit_pdfs[0]["qcd"], fit_pdfs[0]["singletop"], fit_pdfs[0]["wjets"], fit_pdfs[0]["tt_lj"], fit_pdfs[0]["tt_mutau"], fit_pdfs[0]["tt_other"])))
#
#simModel.plotOn(xframe, RooFit.Components("pdf_tt_mutau"), RooFit.LineStyle(kDashed), RooFit.LineColor(kRed))
#simModel.plotOn(xframe, RooFit.Components("pdf_qcd"), RooFit.LineStyle(kDashed), RooFit.LineColor(kGreen))
##model.plotOn(xframe, RooFit.Components("pdf_tt_mutau"), RooFit.LineStyle(kDashed), RooFit.LineColor(kRed))
##model.plotOn(xframe, RooFit.Components("pdf_qcd"), RooFit.LineStyle(kDashed), RooFit.LineColor(kGreen))
##  // Plot multiple background components specified by object reference
##  // Note that specified components may occur at any level in object tree
##  // (e.g bkg is component of 'model' and 'sig2' is component 'sig')
##  model.plotOn(xframe,Components(RooArgSet(bkg,sig2)),LineStyle(kDotted)) ;

# they say this should work for "whole dataset"
# but they don't work with histograms
#data_hist_combined.plotOn(xframe)
#simModel          .plotOn(xframe, RooFit.ProjWData(vars, data_hist_combined))
#
## plotting 1 category
#data_hist_combined.plotOn(xframe, RooFit.Cut("sample==sample::"+category_files[0][0]))
#simModel.plotOn(xframe, RooFit.Slice(sample, category_files[0][0]), RooFit.ProjWData(sample, data_hist_combined))


#mcstudy_Canvas.Divide(2,2)
#mcstudy_Canvas.cd(1)
#frame_cross_par.Draw()
#mcstudy_Canvas.cd(2)
#frame_cross_err.Draw()
#mcstudy_Canvas.cd(3)
#frame_cross_pul.Draw()
#mcstudy_Canvas.cd(4)
#frame_nll.Draw()

#Draw the results
c1 = TCanvas()
c1.Divide(2,2)

c1.cd(1)
xframe1.Draw()

c1.cd(2)
xframe2.Draw()

c1.cd(3)
xframe3.Draw()

c1.cd(4)
xframe4.Draw()

c1.SaveAs("my_combined_histofit.png")

# // P l o t   m o d e l   s l i c e s   o n   d a t a    s l i c e s 
# // ----------------------------------------------------------------
# 
# // Make a frame for the physics sample
# RooPlot* frame1 = x.frame(Bins(30),Title("Physics sample")) ;
# 
# // Plot all data tagged as physics sample
# combData.plotOn(frame1,Cut("sample==sample::physics")) ;
# 
# // Plot "physics" slice of simultaneous pdf. 
# // NBL You _must_ project the sample index category with data using ProjWData 
# // as a RooSimultaneous makes no prediction on the shape in the index category 
# // and can thus not be integrated
# simPdf.plotOn(frame1,Slice(sample,"physics"),ProjWData(sample,combData)) ;
# simPdf.plotOn(frame1,Slice(sample,"physics"),Components("px"),ProjWData(sample,combData),LineStyle(kDashed)) ;
# 
# // The same plot for the control sample slice
# RooPlot* frame2 = x.frame(Bins(30),Title("Control sample")) ;
# combData.plotOn(frame2,Cut("sample==sample::control")) ;
# simPdf.plotOn(frame2,Slice(sample,"control"),ProjWData(sample,combData)) ;
# simPdf.plotOn(frame2,Slice(sample,"control"),Components("px_ctl"),ProjWData(sample,combData),LineStyle(kDashed)) ;
# 
# 
# 
# TCanvas* c = new TCanvas("rf501_simultaneouspdf","rf403_simultaneouspdf",800,400) ;
# c->Divide(2) ;
# c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;
# c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;

print "##############"

for i, category in enumerate(N_parameters_fit):
    print category_files[i][0]
    for n, par in category.items():
        par.Print()

#tt_xsec         .Print()

LUMI_FACTOR     .Print()
tt_factor       .Print()
tauID_factor    .Print()
tt_taufakefactor.Print()

pu_shift.Print()
top_pt_shift.Print()
tau_es_shift.Print()

xsec_dyjets_factor   .Print()
xsec_wjets_factor    .Print()
xsec_singletop_factor.Print()
#model_factors = 

#print models [0][1].plotOn.__doc__

#Now save the data and the PDF into a Workspace, for later use for statistical analysis
#ws = RooWorkspace("ws")
#getattr(ws,'import')(dataset)
#getattr(ws,'import')(totPDF)
#
#fOutput = TFile("Workspace_my_fit.root","RECREATE")
#ws.Write()
#fOutput.Write()
#fOutput.Close()

# nope, just import the module to rerun everything each time

print distr

