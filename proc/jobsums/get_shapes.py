import argparse
import logging
from os.path import isfile
from sys import exit
import pdb

from draw_overflows import DrawOverflow



parser = argparse.ArgumentParser(
    formatter_class = argparse.RawDescriptionHelpFormatter,
    description = "to compare shapes from different files",
    epilog = 'Example:\npython jobsums/get_shapes.py --formula "h1 : h2 " -c ctr_old_mu_sel_lj -d Mt_lep_met_f -p tt_lj h1:jobsums/AN-v1_mu_file-orig.root h2:merge-sets/v23/p10/MC2016_Summer16_TTJets_powheg.root'
    )

parser.add_argument("-c", "--channel", type=str, default='mu_sel', help="channel")
parser.add_argument("-p", "--process", type=str, default='tt_lj', help="process")
parser.add_argument('-s', '--systematic',  type=str, default='NOMINAL', help="systematic")
parser.add_argument('-d', '--distr',  type=str, default='Mt_lep_met_f', help="distr")
parser.add_argument('--debug',    action='store_true', help="logging the debug messages")
parser.add_argument('--logy',     action='store_true', help="log Y")
parser.add_argument('--no-norm',  action='store_true', help="don't normalize each histo to 1")
parser.add_argument('--norm-binwidth',  action='store_true', help="normalize each bin to width")
parser.add_argument('--norm-formulas',  action='store_true', help="normalize the formulas of histos to 1")
parser.add_argument("--y-range",     type=str,      help="set Y range as `ymin,ymax`")
parser.add_argument("--x-title",     type=str,      help="title of X axis")
parser.add_argument("--y-title",     type=str,      help="title of Y axis")
parser.add_argument("--fonts1",      type=int, default=43, help="axis font")
parser.add_argument("--fonts1-size", type=int, default=25, help="axis and legend font size")
parser.add_argument("--font-size-axes-titles", type=int, default=40, help="axis and legend font size")
parser.add_argument("--title",       type=str,      help="title of the pad")
parser.add_argument('--left-title',  type=str, default="CMS #font[52]{Simulation}", help="the text of the left label title")
parser.add_argument('--left-title-inside',  action='store_true', help="put the left label inside the figure box")

parser.add_argument("--draw-overflows", type=float, help="draw the bin with overflows in the given width")

parser.add_argument("--vertical-line",  type=str,   help="format X,Ymin,Ymax draw a line from Ymin to Ymax at the given X")

parser.add_argument("--add-outlier-plot",  type=str,   help="format S,B with Signal and Background outlier events")


parser.add_argument("--outname",     type=str,      help="set the output filename")

parser.add_argument('-e', '--draw-error-bars',  action='store_true', help="root's option e")

parser.add_argument('--no-legend',  action='store_true', help="do not plot the legend")

parser.add_argument('--formula',  type=str, help="to plot h1 overlayed with 123*h2+982*h3: `h1 : 123 h2 , 982h3`")

parser.add_argument('--set-file',  type=str, help=" histos come from this file, skip filename in the input_files definition, but : must be there")

parser.add_argument('--no-std-path',  action='store_true', help="the histo is directly given by the path, it is not the standard `chan_proc_sys_distr`")

parser.add_argument("--output-type",       type=str, default='png', help="output type (png, eps etc)")
#parser.add_argument('input_files', nargs='+', help="""the files process, `histo_name:filename[:chan,proc,sys[,distr]]`""")
parser.add_argument('input_files', nargs='+', help="""the files process, `histo_name:filename[:chan/proc/sys[/distr]]` with any `p/a/th` to a TH histo""")

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)

if args.y_range:
    y_min, y_max = [float(x) for x in args.y_range.split(',')]
else:
    y_min = y_max = None


# process title to get rid of escaping
x_title = args.x_title
y_title = args.y_title


logging.info("import ROOT")

import ROOT
from ROOT import gStyle, gROOT, gPad, TFile, TCanvas, TPad, THStack, TH1D, TLegend, TLine, TArrow, TPaveText, TText, kGreen, kYellow, kOrange, kViolet, kAzure, kWhite, kGray, kRed, kCyan, kBlue, kBlack
gROOT.SetBatch()
gStyle.SetOptStat(0)
#gStyle.SetLineWidth(2)

from plotting_root import rgb, nick_colour

default_colors = [kRed, kBlue-4, kViolet-2, kBlack, kBlack]
color_i = 0
def new_color():
    global color_i
    if default_colors:
        return default_colors.pop()
    else:
        color_i += 1
        return color_i

colors = {'red': kRed, 'black': kBlack, 'blue': kBlue}
styles = {}

if args.left_title_inside:
    leg = TLegend(0.6, 0.55, 0.89, 0.89)
else:
    leg = TLegend(0.7, 0.55, 0.99, 0.89)
leg.SetBorderSize(0)
leg.SetTextFont(args.fonts1)
leg.SetTextSize(args.fonts1_size)

#leg.SetEntrySeparation(0.1)

histos = {}
legend_names = {
'tau40': '30<p_{T}<40',
'tau60': '40<p_{T}<60',
'tauMore60': '60<p_{T}',
'tau40/inclusive2': '30<p_{T}<40',
'tau60/inclusive2': '40<p_{T}<60',
'tauMore60/inclusive2': '60<p_{T}',

'signal':     '#splitline{Signal}{t#bar{t}#rightarrow l#nu #tau_{h}#nu b#bar{b}}',
'background': "#splitline{Background}{t#bar{t}#rightarrow l#nu q#bar{q}' b#bar{b}}",
}

for i, fileparameter in enumerate(args.input_files):
    pars = fileparameter.split(':')
    assert len(pars) in (2, 3)

    #if len(pars) == 4:
    #    legend, nick, filename, opts = pars
    #    logging.debug('legend name for  %s  is  %s' % (nick, legend))
    #    legend_names[nick] = legend

    if len(pars) == 3:
        nick, filename, path_to_h = pars
    else:
        nick, filename = pars
        path_to_h = None

    if args.set_file:
        filename = args.set_file

    if not isfile(filename):
        logging.info("missing: " + filename)
        continue
    logging.debug("%s" % filename)

    if args.no_std_path and path_to_h:
        histo_path = path_to_h

    else:
        channel = process = systematic = distr = ''
        opts = path_to_h
        if opts and len(opts.split('/')) == 4:
            channel, process, systematic, distr = opts.split('/')
        elif opts and len(opts.split('/')) == 3:
            channel, process, systematic = opts.split('/')
        elif opts and len(opts.split('/')) == 2:
            channel, process = opts.split('/')

        if not channel:
            channel    = args.channel
        if not process:
            process    = args.process
        if not systematic:
            systematic = args.systematic
        if not distr:
            distr = args.distr

        histo_path = "{chan}/{proc}/{sys}/{chan}_{proc}_{sys}_{distr}".format(chan=channel, proc=process, sys=systematic, distr=distr)

    tfile = TFile(filename)
    logging.debug("%s" % histo_path)

    # if requested add the overflow bin
    if args.draw_overflows:
        histo = DrawOverflow(tfile.Get(histo_path), args.draw_overflows)
    else:
        histo = tfile.Get(histo_path)

    histo.Sumw2()

    if args.title:
        histo.SetTitle(args.title)

    logging.debug("histo    : %20s" % nick)
    logging.debug("integral : %f" % histo.Integral())
    histo.SetDirectory(0)
    histo.SetFillColor(0)
    #histo.SetLineColor(1 + i)
    color = new_color()
    histo.SetLineColor(color)
    histo.SetLineWidth(3)

    histo.SetMarkerStyle(19)
    histo.SetMarkerColor(color)
    histo.SetMarkerColorAlpha(color, 0.)

    if not args.no_norm:
        histo.Scale(1./histo.Integral())
    histo.SetName(nick)
    if y_max is not None and y_min is not None:
        histo.SetMaximum(y_max)
        histo.SetMinimum(y_min)
    histos[nick] = histo


cst = TCanvas("cst","stacked hists",10,10,700,700)
cst.Draw("Y+")

y_axis_tick_len = 0.03
bottom_margin = 0.17

if args.add_outlier_plot:
    padsplit_x = 0.89
    pad1 = TPad("pad1","This is pad1", 0.,  0., padsplit_x, 1.)
    pad2 = TPad("pad2","This is pad2", padsplit_x, 0., 1., 1.)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.)
    pad1.SetBottomMargin(bottom_margin)
    pad2.SetTitle("")
    pad2.SetLeftMargin(0.)
    pad2.SetBottomMargin(bottom_margin)
    pad2.SetRightMargin(0.5)
    pad1.Draw()
    pad2.Draw()

else:
    padsplit_x = 1.
    pad1 = TPad("pad1","This is pad1", 0., 0.,  1., 1.)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.01)
    pad1.SetBottomMargin(bottom_margin)
    pad1.Draw()

if args.title:
    pad1.SetTitle(args.title)


pad1.SetTicks()
if args.logy:
    pad1.SetLogy()

if args.add_outlier_plot:
    pad2.SetTicks()
    pad1.SetTicks(1, 0) 
    pad2.SetTicks(0, 0)
    if args.logy:
        pad2.SetLogy()

pad1.cd()

plots_to_legend = []

drawn = False

if not args.formula:
    for i, (nick, h) in enumerate(histos.items()):
        if i == 0:
            if y_max is not None and y_min is not None:
                h.SetMaximum(y_max)
                h.SetMinimum(y_min)
            if args.x_title:
                #nom_MC.SetTitle(args.process)
                h.SetXTitle(x_title)
            if args.y_title:
                h.SetYTitle(y_title)
            if args.title:
                h.SetTitle(args.title)
            h.Draw()
        else:
            h.Draw("same")
        logging.debug('adding from legend %s nick %s' % (repr(legend_names), nick))
        leg.AddEntry(histo, legend_names.get(nick, nick), "l")

else:
    # parse the formula
    import re
    form_histos = []

    # a sub-form is a distribution to plot
    # it can be a linear combination of distrs
    # or a ratio of linear combinations
    for form_i, subform in enumerate(args.formula.split(':')): # sub-formulas
        # ' number  nick ' pairs
        logging.debug("subform " + subform)
        nicks = [nick for nick in histos.keys() if nick in subform]
        logging.debug("nicks %s" % str(nicks))

        # parse styles out of subform
        color = style = None
        if   subform.count(',') == 2:
            subform, color, style = subform.split(',')
        elif subform.count(',') == 1:
            subform, color = subform.split(',')
        elif subform.count(',') == 0:
            pass
        else:
            logging.error('wrong subform, too many "," %s' % subform)
            continue

        # a subform might be a ratio of two distributions
        sums_histos = []
        for sum_i, linear_sum in enumerate(subform.split('/')):

            # get the nicks in the subform's linear_sum and scale them by their factors
            nick_pairs = []
            for nick_pair in linear_sum.split('+'):
                for nick in nicks:
                    if nick in nick_pair:
                        logging.debug("nick_pair " + nick_pair)
                        number = nick_pair.replace(nick, '')
                        logging.debug("number " + number)
                        factor = float(number) if number.strip() else 1.
                        logging.debug("factor %f" % factor)
                        nick_pairs.append((factor, nick))
                        break

            logging.debug(nick_pairs)

            ## normalize the histo
            #if not args.norm:
            #    histo.Scale(1./histo.Integral())

            # construct the histo of the subform
            factor, nick = nick_pairs[0]
            histo = histos[nick].Clone()

            if args.title:
                histo.SetTitle(args.title)

            formula_name = "formula_%d_%d" % (form_i, sum_i)
            histo.SetName(formula_name)
            histo.SetDirectory(0)
            histo.Scale(factor)
            logging.debug("%15s %15s %10.1f" % (formula_name, nick, histo.Integral()))

            # add up all the histos in the subform
            for factor, nick in nick_pairs[1:]:
                # not sure if pyroot can do h += h2*factor correctly, without mutation
                h = histos[nick].Clone()
                h.Scale(factor)
                histo.Add(h)
                logging.debug("%15s %15s %10.1f" % (formula_name, nick, histo.Integral()))

            logging.debug("histo %20s %f" % (nick, histo.Integral()))

            #here histo is the linear combination in the formula
            # save them
            sums_histos.append(histo)

        # how to divide this stuff? a / b / c = (a / b) / c
        histo = sums_histos[0].Clone()
        histo.SetName("formula_%d" % form_i)
        histo.SetDirectory(0)
        for operand_histo in sums_histos[1:]:
            histo.Divide(operand_histo)

        # set histo styles, if given
        if color:
            if color in colors:
                col = colors[color]
            else:
                col = int(color)
            histo.SetLineColor(col)
            histo.SetMarkerColor(col)
            histo.SetMarkerColorAlpha(col, 0.)

        if style:
            stl = styles.get(style, int(style))
            histo.SetLineStyle(stl)

        form_histos.append(histo)

        plots_to_legend.append((histo, subform))

    # and draw
    for histo in form_histos:
        if args.title:
            histo.SetTitle(args.title)
        # normalize bin width if needed
        if args.norm_binwidth:
            #
            for bini in range(histo.GetSize()):
                content = histo.GetBinContent(bini)
                error   = histo.GetBinError(bini)
                width   = histo.GetXaxis().GetBinUpEdge(bini) - histo.GetXaxis().GetBinLowEdge(bini)
                histo.SetBinContent(bini, content/width)
                histo.SetBinError(bini, error/width)

            # rescale to 1 just in case
            #histo.Scale(1./histo.Integral())

        # normalize the linear sum
        if args.norm_formulas:
            histo.Scale(1./histo.Integral())


        if drawn:
            if args.title:
                histo.SetTitle(args.title)
            if y_max is not None and y_min is not None:
                logging.debug("setting min-max")
                histo.SetMaximum(y_max)
                histo.SetMinimum(y_min)
            if args.draw_error_bars:
                histo.Draw("hist e same")
            else:
                histo.Draw("hist same")
        else:
            if args.x_title:
                #nom_MC.SetTitle(args.process)
                histo.SetXTitle(x_title)
                #pdb.set_trace()
            if args.y_title:
                histo.SetYTitle(y_title)
            if args.title:
                histo.SetTitle(args.title)

            #histo.GetYaxis().SetTitleOffset(1.4)
            histo.GetXaxis().SetTitleOffset(0.9)

            histo.GetYaxis().SetLabelFont(args.fonts1)
            histo.GetXaxis().SetLabelFont(args.fonts1)
            #histo.GetYaxis().SetLabelSize(args.fonts1_size*0.75)
            #histo.GetXaxis().SetLabelSize(args.fonts1_size*0.75)
            histo.GetYaxis().SetLabelSize(args.fonts1_size)
            histo.GetXaxis().SetLabelSize(args.fonts1_size)
            histo.GetYaxis().SetTitleOffset(1.15)

            histo.GetYaxis().SetTitleFont(args.fonts1)
            histo.GetXaxis().SetTitleFont(args.fonts1)
            histo.GetYaxis().SetTitleSize(args.font_size_axes_titles)
            histo.GetXaxis().SetTitleSize(args.font_size_axes_titles)
            if y_max is not None and y_min is not None:
                logging.debug("setting min-max")
                histo.SetMaximum(y_max)
                histo.SetMinimum(y_min)

            # set the length of ticks on Y axis
            # it is done consistently with the possible additional plot
            histo.GetYaxis().SetTickLength(y_axis_tick_len/padsplit_x)
            # draw the first one with line
            histo.SetFillColor(0)
            histo.Draw("hist e" if args.draw_error_bars else "hist")
            drawn = True

            # draw the verticle line in the background
            if args.vertical_line:
                x, ymin, ymax = (float(i) for i in args.vertical_line.split(','))
                #l = TLine(x, ymin, x, ymax)
                l = TArrow(x, ymin, x, ymax, 0.05, "<|")
                l.SetAngle(40)
                l.SetLineWidth(1)
                l.Draw()
                #l.Draw("same")

                #drawn = True
                #histo.Draw("hist same")
                histo.Draw("hist e same" if args.draw_error_bars else "hist same")

for histo, subform in reversed(plots_to_legend):
    if subform.strip() in legend_names:
        leg.AddEntry(histo, legend_names[subform.strip()]) #, "l")
    else:
        leg.AddEntry(histo, subform) #, "l")


if not args.no_legend:
    leg.Draw("same")

if args.add_outlier_plot:
    sig, bkg = (float(i) for i in args.add_outlier_plot.split(','))
    pad2.cd()
    pad2.Draw("Y+")
    cst.Draw("Y+")

    #l = TLine(0., 0.5, 1., 0.5)
    #l.Draw("same")

    sig_hist = TH1D("untagged_events_signal", "", 1, 0., 1.)
    sig_hist.SetBinContent(1, sig)
    sig_hist.SetMaximum(y_max)
    sig_hist.SetMinimum(y_min)
    sig_hist.SetLineColor(colors['red'])
    sig_hist.SetLineWidth(3)

    # the Y axis title
    sig_hist.GetYaxis().SetTitle("Fraction of events")
    sig_hist.GetYaxis().SetTitleFont(args.fonts1)
    sig_hist.GetYaxis().SetTitleOffset(2.5)
    sig_hist.GetYaxis().SetTitleSize(args.font_size_axes_titles)
    sig_hist.GetYaxis().SetLabelSize(0.)

    # remove ticks
    sig_hist.GetXaxis().SetTickLength(0.)

    # Set the correct length of Y axis ticks
    sig_hist.GetYaxis().SetTickLength(y_axis_tick_len/(1.-padsplit_x))

    ## remove the labels on ticks
    #sig_hist.GetXaxis().SetLabelOffset(999)
    #sig_hist.GetXaxis().SetLabelSize(0)

    # set "untagged" as a label
    #sig_hist.GetXaxis().SetLabelOffset(999)
    sig_hist.GetXaxis().SetLabelFont(args.fonts1)
    sig_hist.GetXaxis().SetLabelSize(args.fonts1_size) # does not fit
    #sig_hist.GetXaxis().SetLabelSize(25)

    # set the title "untagged"
    #sig_hist.GetXaxis().SetTitle("untagged")
    # set alphanumeric label
    sig_hist.GetXaxis().SetBinLabel(1, "untagged")
    sig_hist.LabelsOption("v")

    sig_hist.Draw("same Y+") # the Y axis ticks on the right side!!
    sig_hist.Draw("same")

    bkg_hist = TH1D("untagged_events_background", "", 1, 0., 1.)
    bkg_hist.SetBinContent(1, bkg)
    bkg_hist.SetLineColor(colors['black'])
    bkg_hist.SetLineWidth(3)
    bkg_hist.SetLineStyle(2)
    bkg_hist.Draw("same")


cst.cd()

if args.left_title:
    #left_title = TPaveText(0.12, 0.8, 0.36, 0.88, "brNDC")
    if args.left_title_inside:
        left_title = TPaveText(0.19, 0.82, 0.5, 0.87, "brNDC")
    else:
        left_title = TPaveText(0.15,  0.92, 0.5, 0.97, "brNDC")
    #left_title.SetBorderSize(1)
    #left_title.SetTextAlign(11) # this adjust to left (+ margin inside the box and top)
    left_title.SetTextAlign(13) # this adjust to left (+ margin inside the box and top)
    left_title.SetMargin(0)
    #left_title.AddText("CMS simulation")
    left_title.AddText(args.left_title)
    #left_title.SetTextFont(1) # italic/cursive
    #left_title.SetTextFont(1)
    left_title.SetFillColor(0)

    ## this does not show up on the plot
    #left_title = TText(0.05, 0.905, "CMS simulation")
    #left_title.SetTextAlign(11)
    #left_title.SetTextSize(0.12)

    left_title .Draw("same")

#right_title = TPaveText(0.65, 0.91, 0.89, 0.96, "brNDC")
#right_title = TPaveText(0.65, 0.92, 0.9 + (1 - padsplit_x), 0.97, "brNDC")
#right_title = TPaveText(0.65, 0.92, 0.89 + (1 - padsplit_x), 0.97, "brNDC")
if args.add_outlier_plot:
    right_title = TPaveText(0.65, 0.92, 0.945, 0.97, "brNDC")
else:
    right_title = TPaveText(0.65, 0.92, 0.99, 0.97, "brNDC")
right_title.SetTextAlign(33)
right_title.SetMargin(0)
right_title.AddText("(13 TeV)")
right_title.SetTextFont(132)
right_title.SetFillColor(0)

right_title.Draw("same")


if args.outname:
    cst.SaveAs(args.outname)
else:
    outname = './compare_%s_%s_%s_%s_%s%s%s.' % (args.channel, args.process, args.systematic, args.distr, ','.join(histos.keys()), '_formula' if args.formula else '', '_nonorm' if args.no_norm else '')
    outname += args.output_type
    cst.SaveAs(outname)

