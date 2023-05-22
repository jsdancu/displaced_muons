import sys
import os
import ROOT
ROOT.gROOT.SetBatch(1)
import style
ROOT.gStyle.SetPalette(ROOT.kCherry)
#style.makeColorTable()
#ROOT.gStyle.SetPalette(256, "acton.txt")
# ROOT.TColor.CreateColorTableFromFile("acton.txt")
# ROOT.gStyle.SetNumberContours(256)
ROOT.TColor.InvertPalette()
from array import array
import random
import numpy as np
#import root_numpy

import stats_functions

def define_plots_1D(binning1, xaxis_title, yaxis_title, title, bin_edges=0):

    if bin_edges:
        nbins = len(binning1) - 1
        hist1 = ROOT.TH1D(title, title, nbins, binning1)
    else:
        hist1 = ROOT.TH1D(title, title, *binning1)

    hist1.GetXaxis().SetTitle(xaxis_title)
    hist1.GetYaxis().SetTitle(yaxis_title)
    hist1.GetXaxis().SetTitleOffset(0.8)
    hist1.GetYaxis().SetTitleOffset(0.8)

    return hist1

def define_plots_2D(binning1, binning2, xaxis_title, yaxis_title, title, bin_edges1=0, bin_edges2=0):

    if bin_edges1 and bin_edges2:
        nbins1 = len(binning1) - 1
        nbins2 = len(binning2) - 1
        hist1 = ROOT.TH2D(title, title, nbins1, binning1, nbins2, binning2)
    elif  (not bin_edges1) and bin_edges2:
        nbins2 = len(binning2) - 1
        hist1 = ROOT.TH2D(title, title, *binning1, nbins2, binning2)
    elif bin_edges1 and (not bin_edges2):
        nbins1 = len(binning1) - 1
        hist1 = ROOT.TH2D(title, title, nbins1, binning1, *binning2)
    else:
        hist1 = ROOT.TH2D(title, title, *binning1, *binning2)

    #hist1.SetTitle(title)
    hist1.GetXaxis().SetTitle(xaxis_title)
    hist1.GetYaxis().SetTitle(yaxis_title)
    hist1.GetYaxis().SetTitleOffset(0.9)

    return hist1

def plotting_1D(array1, hist1, puWeight=None):

    ntimes = len(array1)

    if puWeight is None:
        puWeight = np.ones(ntimes)

    hist1.FillN(ntimes, array('d', array1), array('d', puWeight))

    return hist1

def plotting_2D(array1, array2, hist1, puWeight=None):

    ntimes = len(array1)

    if puWeight is None:
        puWeight = np.ones(ntimes)

    hist1.FillN(ntimes, array('d', array1), array('d', array2), array('d', puWeight))

    return hist1

def fill_hist_1D(binning, xtitle, title, tree, weights):

    hist = define_plots_1D(binning, xtitle, "Entries", title+str(random.randint(1, 10000000)))
    hist = plotting_1D(tree, hist, weights)
    
    return hist

def fill_hist_2D(binning1, binning2, xtitle, ytitle, title, tree1, tree2, weights):

    hist = define_plots_2D(binning1, binning2, xtitle, ytitle, title+str(random.randint(1, 10000000)))
    hist = plotting_2D(tree1, tree2, hist, weights)
    
    return hist

def plotting(hist, plotdir, plot_name, var, CR, corr_coeff, category, charge, year,luminosity, process, xaxis_log_scale, yaxis_log_scale, drawstyle="HIST"):

    cv = ROOT.TCanvas()
    cv.Draw()
    cv.SetBottomMargin(0.2)
    cv.SetLeftMargin(0.13)
    cv.SetRightMargin(0.17)
   
    if xaxis_log_scale>0:
        cv.SetLogx()
    if yaxis_log_scale>0:
        cv.SetLogy()

    hist.Draw(drawstyle)
    
    category_name_dict = {}
    category_name_dict["muonmuon"] = "#mu#mu"
    category_name_dict["muonelectron"] = "#mue"
    category_name_dict["electronmuon"] = "e#mu"
    category_name_dict["electronelectron"] = "ee"
    category_name_dict["all_cat"] = "all categories"
    category_name = category_name_dict[category]

    if CR == "CR_DY":
        control_reg = "DY+jets CR"
    elif CR == "CR_resolved":
        control_reg = "resolved CR"
    elif CR == "CR_3lep":
        control_reg = "3 lepton CR"
    else:
        control_reg = CR

    style.makeLumiText(0.75, 0.97, lumi=str(luminosity), year=str(year))
    #style.makeText(0.75, 0.93, 0.95, 0.97, "("+str(year)+")")
    style.makeCMSText(0.13, 0.97,additionalText="Preliminary", dx=0.1)

    if process == "":
        style.makeText(0.2, 0.80, 0.3, 0.90, category_name+", "+charge+", "+control_reg)
    else:    
        style.makeText(0.2, 0.80, 0.3, 0.90, category_name+", "+charge+", "+control_reg+", "+str(process))
    if corr_coeff:
        style.makeText(0.08, 0.07, 0.4, 0.08, "mutual info. = {0:.3g}".format(corr_coeff))

    cv.Modified()

    # if process == "":
    #     cv.SaveAs(os.path.join(plotdir, str(var)+"_"+str(CR)+"_"+str(category)+"_"+str(charge)+"_"+str(year)+".pdf"))
    #     cv.SaveAs(os.path.join(plotdir, str(var)+"_"+str(CR)+"_"+str(category)+"_"+str(charge)+"_"+str(year)+".png"))
    # else:
    #     cv.SaveAs(os.path.join(plotdir, str(var)+"_"+str(CR)+"_"+str(category)+"_"+str(charge)+"_"+str(year)+"_"+str(process)+".pdf"))
    #     cv.SaveAs(os.path.join(plotdir, str(var)+"_"+str(CR)+"_"+str(category)+"_"+str(charge)+"_"+str(year)+"_"+str(process)+".png"))

    cv.SaveAs(os.path.join(plotdir, plot_name+".pdf"))
    cv.SaveAs(os.path.join(plotdir, plot_name+".png"))

    cv.Close()

def plotting_joint(hist1, hist2, hist3, hist4, plotdir, type, CR, category, year, drawstyle="HIST"):

    cv = ROOT.TCanvas()
    cv.Draw()
    cv.SetBottomMargin(0.2)
    cv.SetLeftMargin(0.13)
    cv.SetRightMargin(0.05)

    hist1.SetMaximum(hist1.GetMaximum()*25)

    hist1.Scale(1./hist1.Integral())
    hist2.Scale(1./hist2.Integral())
    hist3.Scale(1./hist3.Integral())
    hist4.Scale(1./hist4.Integral())

    hist1.SetLineColor(ROOT.kBlack)
    hist1.Draw(drawstyle)
    hist2.SetLineColor(ROOT.kBlue)
    hist2.Draw(drawstyle+" SAME")
    hist3.SetLineColor(ROOT.kRed)
    hist3.Draw(drawstyle+" SAME")
    hist4.SetLineColor(ROOT.kGreen)
    hist4.Draw(drawstyle+" SAME")

    legend = style.makeLegend(0.7,0.45,0.85,0.75)
    legend.AddEntry(hist1, "sig", "f")
    legend.AddEntry(hist2, "sig ev. w.", "f")
    legend.AddEntry(hist3, "bkg", "f")
    legend.AddEntry(hist4, "bkg ev. w.", "f")
    legend.Draw("SAME")

    if CR == "CR_DY":
        control_reg = "DY+jets CR"
    elif CR == "CR_resolved":
        control_reg = "resolved CR"
    else:
        control_reg = CR

    #style.makeLumiText(0.75, 0.97, lumi=str(args.luminosity), year=str(args.year))
    style.makeText(0.75, 0.93, 0.95, 0.97, "("+str(year)+")")
    style.makeCMSText(0.13, 0.97,additionalText="Preliminary", dx=0.1)
    style.makeText(0.2, 0.80, 0.3, 0.90, category+", SS+OS, "+control_reg)

    cv.Modified()

    cv.SaveAs(os.path.join(plotdir, str(type)+"_"+str(CR)+"_"+str(category)+"_SS+OS_"+str(year)+".pdf"))
    cv.SaveAs(os.path.join(plotdir, str(type)+"_"+str(CR)+"_"+str(category)+"_SS+OS_"+str(year)+".png"))

    cv.Close()

def plotting_2in1(hist1, hist2, plotdir, var, CR, category, charge, year, process, drawstyle="HIST"):
        
    cv = ROOT.TCanvas()
    cv.Draw()
    cv.SetBottomMargin(0.2)
    cv.SetLeftMargin(0.13)
    cv.SetRightMargin(0.17)

    hist1.Scale(1./hist1.Integral())
    hist2.Scale(1./hist2.Integral())
    hist1.SetMaximum(hist1.GetMaximum()*2)

    hist1.SetLineColor(ROOT.TColor.GetColor("#ef5350"))
    hist1.Draw(drawstyle)
    hist2.SetLineColor(ROOT.TColor.GetColor("#388e3c"))
    hist2.Draw(drawstyle+" SAME")
    
    legend = style.makeLegend(0.85,0.4,0.98,0.7)
    legend.AddEntry(hist1, "sig","l")
    legend.AddEntry(hist2, "bkg","l")
    legend.Draw("SAME")

    category_name_dict = {}
    category_name_dict["muonmuon"] = "#mu#mu"
    category_name_dict["muonelectron"] = "#mue"
    category_name_dict["electronmuon"] = "e#mu"
    category_name_dict["electronelectron"] = "ee"
    category_name_dict["all_cat"] = "all categories"
    category_name = category_name_dict[category]

    if CR == "CR_DY":
        control_reg = "DY+jets CR"
    elif CR == "CR_resolved":
        control_reg = "resolved CR"
    else:
        control_reg = CR

    #style.makeLumiText(0.75, 0.97, lumi=str(args.luminosity), year=str(args.year))
    style.makeText(0.75, 0.93, 0.95, 0.97, "("+str(year)+")")
    style.makeCMSText(0.13, 0.97,additionalText="Preliminary", dx=0.1)
    style.makeText(0.2, 0.80, 0.3, 0.90, category_name+", "+charge+", "+control_reg+", "+process)

    cv.Modified()

    cv.SaveAs(os.path.join(plotdir, str(var)+"_"+str(CR)+"_"+str(category)+"_"+charge+"_"+str(year)+"_"+str(process)+".pdf"))
    cv.SaveAs(os.path.join(plotdir, str(var)+"_"+str(CR)+"_"+str(category)+"_"+charge+"_"+str(year)+"_"+str(process)+".png"))

    cv.Close()

def plotting_scale_factor_2D(hist, hist_val, hist_err_down, hist_err_up, edges1, edges2, plotdir, plot_name, var, CR, category, category_name, charge, year,luminosity, xaxis_log_scale, yaxis_log_scale):
        
    cv = ROOT.TCanvas()
    cv.Draw()
    cv.SetBottomMargin(0.2)
    cv.SetLeftMargin(0.13)
    cv.SetRightMargin(0.17)
   
    if xaxis_log_scale>0:
        cv.SetLogx()
    if yaxis_log_scale>0:
        cv.SetLogy()

    #hist.SetMaximum(np.ceil(hist.GetMaximum()))
    # hist.SetMaximum(1.05)
    # hist.SetMinimum(0.8)
    if hist.GetMaximum() < 1.5:
        hist.SetMaximum(hist.GetMaximum()+0.3)
    else:
        hist.SetMaximum(hist.GetMaximum()+1.0)
    
    if hist.GetMinimum() < 0.5:
        hist.SetMinimum(0.0)
    else:
        hist.SetMinimum(hist.GetMinimum()-0.1)

    # hist.SetMaximum(0.1)
    # hist.SetMinimum(0.0)
    
    hist.Draw("COLZ")

    width1 = []
    for i in range(len(edges1)-1):
        width1.append(edges1[i+1] - edges1[i])
    x = edges1[:-1] + 0.1 * np.array(width1)

    width2 = []
    for i in range(len(edges2)-1):
        width2.append(edges2[i+1] - edges2[i])
    y = edges2[:-1] + 0.5 * np.array(width2)

    text = ROOT.TLatex()
    text.SetTextSize(.022)
    text.SetIndiceSize(2.0)
    # if hist_err_up==[]:
    #     for i, x_val in enumerate(x):
    #         for j, y_val in enumerate(y):
    #             text.DrawLatex(x_val,y_val,f'{hist_val[i,j]:.2f}#pm{hist_err_down[i,j]:.2f}')
    #             #change bin error value
    #             #hist.SetBinError(i+1,j+1,hist_err_down[i,j])
    # else:
    #     for i, x_val in enumerate(x):
    #         for j, y_val in enumerate(y):
    #             text.DrawLatex(x_val,y_val,f'{hist_val[i,j]:.2f}^{{+{hist_err_up[i,j]:.2f}}}_{{-{hist_err_down[i,j]:.2f}}}')
    #             #change bin error value
    #             #hist.SetBinError(i+1,j+1,(hist_err_down[i,j]+hist_err_up[i,j])/2.0)

    for i in range(hist.GetNbinsX()):
        for j in range(hist.GetNbinsY()):
            # print("i: ", i)
            # print("j: ", j)
            sf = hist.GetBinContent(i+1,j+1)
            sf_err_stats = hist_err_down[i, j]
            sf_err_syst = hist_err_up[i, j]

            text.DrawLatex(x[i],y[j],f'{sf:.2f}^{{#pm{sf_err_stats:.2f} (stats)}}_{{#pm{sf_err_syst:.2f} (syst)}}')

    # for i in range(hist.GetNbinsX()):
    #     for j in range(hist.GetNbinsY()):
    #         # print("i: ", i)
    #         # print("j: ", j)
    #         sf = hist.GetBinContent(i+1,j+1)
    #         sf_err = hist.GetBinError(i+1,j+1)
    #         if hist_err_down == []:
    #             text.DrawLatex(x[i],y[j],f'{sf:.2f}')
    #         else:
    #             text.DrawLatex(x[i],y[j],f'{sf:.2f}#pm{sf_err:.2f}')

    # for i in range(hist.GetNbinsX()):
    #     for j in range(hist.GetNbinsY()):
    #         # print("i: ", i)
    #         # print("j: ", j)
    #         sf = hist.GetBinContent(i+1,j+1)
    #         if hist_err_down == []:
    #             text.DrawLatex(x[i],y[j],f'{sf:.10f}')

    if CR == "CR_DY":
        control_reg = "DY+jets CR"
    elif CR == "CR_resolved":
        control_reg = "resolved CR"
    elif CR == "CR_3lep":
        control_reg = "3 lepton CR"
    else:
        control_reg = CR

    style.makeLumiText(0.75, 0.97, lumi=str(luminosity), year=str(year))
    #style.makeText(0.75, 0.93, 0.95, 0.97, "("+str(year)+")")
    style.makeCMSText(0.13, 0.97,additionalText="Preliminary", dx=0.1)

    #style.makeText(0.2, 0.80, 0.3, 0.90, category_name+", "+charge+", "+control_reg)

    cv.Modified()

    cv.SaveAs(os.path.join(plotdir, plot_name+".pdf"))
    cv.SaveAs(os.path.join(plotdir, plot_name+".png"))

    ####### Save scale factor historgram in root file ###############
    #cv.SaveAs(os.path.join(plotdir, plot_name+".root"))

    cv.Close()

    # file = ROOT.TFile.Open(os.path.join(plotdir, plot_name+".root"), "RECREATE")   
    # hist.Write()
    # file.Close()

"""
def plotting_scale_factor_1D(MC_hists_noID, data_hist_noID, MC_hists_ID, data_hist_ID, MC_hist_err_noID, MC_hist_err_ID, plotdir, var, var_name, MC_types, var_min, var_max, CR, category, category_name, dilep_charge, year, luminosity, xaxis_log_scale):

    hist_stack_noID = ROOT.THStack("mumueOS"+var+"MC", "mumueOS"+var+"MC")
    hist_stack_ID = ROOT.THStack("mumueIDOS"+var+"MC", "mumueIDOS"+var+"MC")

    legend = style.makeLegend(0.72,0.35,0.98,0.7)
    legend.SetHeader("#mu#mue(CustomID)")

    for mc in reversed(list(MC_hists_noID.keys())):
        hist_stack_noID.Add(MC_hists_noID[mc])

    for mc in reversed(list(MC_hists_ID.keys())):
        hist_stack_ID.Add(MC_hists_ID[mc])
        
    cv = ROOT.TCanvas()
    cv.Draw()
    cv.SetBottomMargin(0.2)
    cv.SetLeftMargin(0.13)
    #cv.SetRightMargin(0.25)

    upperPad = ROOT.TPad("upperPad", "upperPad", 0, 0.33, 1, 1)
    lowerPad = ROOT.TPad("lowerPad", "lowerPad", 0, 0, 1, 0.33)

    upperPad.SetBottomMargin(0.00001)
    upperPad.SetLeftMargin(0.12)
    upperPad.SetRightMargin(0.28)
    upperPad.SetBorderMode(0)
    upperPad.SetTopMargin(0.15)
    lowerPad.SetTopMargin(0.00001)
    lowerPad.SetBottomMargin(0.4)
    lowerPad.SetLeftMargin(0.12)
    lowerPad.SetRightMargin(0.28)
    lowerPad.SetBorderMode(0)
    if xaxis_log_scale>0:
        upperPad.SetLogx()
        lowerPad.SetLogx()
    upperPad.Draw()
    lowerPad.Draw()

    upperPad.cd()
    
    array_data_ID, edges = root_numpy.hist2array(data_hist_ID, return_edges=True)
    array_data_noID = root_numpy.hist2array(data_hist_noID)
    edges = edges[0]

    width = []
    for i in range(len(edges)-1):
        width.append(edges[i+1] - edges[i])

    xerr = 0.5 * np.array(width)
    x = edges[:-1] + 0.5 * np.array(width)

    array_eff_data, (array_eff_data_down, array_eff_data_up) = stats_functions.divide_efficiency(array_data_ID, array_data_noID)

    hist_eff_data = ROOT.TGraphAsymmErrors(len(x), x, array_eff_data, xerr, xerr, array_eff_data_down, array_eff_data_up)

    hist_eff_data.SetMarkerColor(ROOT.TColor.GetColor("#cf1f1f"))
    hist_eff_data.SetLineColor(ROOT.TColor.GetColor("#cf1f1f"))
    hist_eff_data.SetMaximum(1.5)
    hist_eff_data.SetMinimum(0.0)
    hist_eff_data.GetYaxis().SetTitleOffset(0.8)
    hist_eff_data.GetYaxis().SetTitle("Efficiency")

    # hist_eff_data.Draw("AP")
    # upperPad.Update()
    
    legend.AddEntry(hist_eff_data, "#varepsilon(data)","p")

    array_stack_ID = root_numpy.hist2array(hist_stack_ID.GetStack().Last())
    array_stack_noID = root_numpy.hist2array(hist_stack_noID.GetStack().Last())
    array_MC_err_ID = np.sqrt(root_numpy.hist2array(MC_hist_err_ID))
    array_MC_err_noID = np.sqrt(root_numpy.hist2array(MC_hist_err_noID))

    #array_eff_MC, (array_eff_MC_down, array_eff_MC_up) = stats_functions.divide_efficiency(array_stack_ID, array_stack_noID)

    array_eff_MC = array_stack_ID/array_stack_noID
    array_eff_MC_err = stats_functions.error_propagation(array_eff_MC, array_stack_ID, array_MC_err_ID, array_stack_noID, array_MC_err_noID)

    hist_eff_MC = ROOT.TGraphAsymmErrors(len(x), x, array_eff_MC, xerr, xerr, array_eff_MC_err, array_eff_MC_err)

    hist_eff_MC.SetMarkerColor(ROOT.TColor.GetColor("#2d67cc"))
    hist_eff_MC.SetLineColor(ROOT.TColor.GetColor("#2d67cc"))
    hist_eff_MC.SetMaximum(1.5)
    hist_eff_MC.SetMinimum(0.0)
    hist_eff_MC.GetYaxis().SetTitleOffset(0.8)
    hist_eff_MC.GetYaxis().SetTitle("Efficiency")

    # hist_eff_MC.Draw("SAME AP")
    # upperPad.Update()

    multigraph = ROOT.TMultiGraph("Efficiency", "Efficiency")
    multigraph.Add(hist_eff_data)
    multigraph.Add(hist_eff_MC)

    multigraph.Draw("AP")

    multigraph.SetMaximum(1.5)
    multigraph.SetMinimum(0.0)
    multigraph.GetXaxis().SetLimits(var_min, var_max)
    multigraph.GetYaxis().SetTitleOffset(0.8)
    multigraph.GetYaxis().SetTitle("Efficiency")

    upperPad.Update()

    legend.AddEntry(hist_eff_MC, "#varepsilon(simulation)","p")
    legend.Draw("SAME")

    line = ROOT.TLine(var_min, 1, var_max, 1)
    line.Draw("SAME")
    upperPad.Update()
    

    ########## Moving to lower pad -> scale factors ###########
    lowerPad.cd()

    print("array_data_ID: ", array_data_ID)
    print("array_data_noID: ", array_data_noID)
    print("array_stack_ID: ", array_stack_ID)
    print("array_stack_noID: ", array_stack_noID)

    array_sf, (array_sf_down, array_sf_up) = stats_functions.binomial_ratio(array_data_ID, array_data_noID, array_stack_ID, array_stack_noID)

    print("array_sf: ", array_sf)
    print("array_sf_down: ", array_sf_down)
    print("array_sf_up: ", array_sf_up)

    hist_sf = ROOT.TGraphAsymmErrors(len(x), x, array_sf, xerr, xerr, array_sf_down, array_sf_up)

    hist_sf.SetMarkerColor(ROOT.kBlack)
    hist_sf.SetLineColor(ROOT.kBlack)
    hist_sf.SetMaximum(2.0)
    hist_sf.SetMinimum(0.0)
    hist_sf.GetXaxis().SetTitle(var_name)
    hist_sf.GetXaxis().SetTitleOffset(2.5)
    hist_sf.GetXaxis().SetLimits(var_min, var_max)
    hist_sf.GetYaxis().SetTitleOffset(0.8)
    hist_sf.GetYaxis().SetTitle("Scale factor")

    hist_sf.Draw("AP")
    lowerPad.Update()

    line.Draw("SAME")

    cv.cd()

    style.makeLumiText(0.85, 0.97, lumi=str(luminosity), year=str(year))
    style.makeCMSText(0.13, 0.97,additionalText="Preliminary", dx=0.1)

    if CR == "CR_3lep":
        control_reg = "3 lepton CR"
    else:
        control_reg = CR

    #style.makeText(0.2, 0.80, 0.3, 0.90, category_name+", "+dilep_charge+", "+control_reg)

    cv.Modified()

    cv.SaveAs(os.path.join(plotdir, "scale_factor_"+str(var)+"_"+str(CR)+"_"+str(category)+"_"+str(dilep_charge)+"_"+str(year)+".pdf"))
    cv.SaveAs(os.path.join(plotdir, "scale_factor_"+str(var)+"_"+str(CR)+"_"+str(category)+"_"+str(dilep_charge)+"_"+str(year)+".png"))

    cv.Close()
"""

def plotting_data_to_MC(MC_hists, hist_data, plotdir, var, var_name, MC_types, var_min, var_max, CR, category, category_name, dilep_charge, year, luminosity, log_scale, plotting, plotting_word, plotting_word_latex, xaxis_log_scale):

    hist_stack = ROOT.THStack(category+dilep_charge+var+"MC", category+dilep_charge+var+"MC")

    legend = style.makeLegend(0.8,0.01,0.98,0.7)

    for mc in reversed(list(MC_hists.keys())):

        MC_hists[mc].SetLineColor(ROOT.TColor.GetColor(MC_types[mc][2]))
        MC_hists[mc].SetFillColor(ROOT.TColor.GetColor(MC_types[mc][1]))
            
        hist_stack.Add(MC_hists[mc])
        legend.AddEntry(MC_hists[mc], MC_types[mc][3],"f")

    cv = ROOT.TCanvas()
    cv.Draw()
    cv.SetBottomMargin(0.2)
    cv.SetLeftMargin(0.13)
    #cv.SetRightMargin(0.25)

    upperPad = ROOT.TPad("upperPad", "upperPad", 0, 0.33, 1, 1)
    lowerPad = ROOT.TPad("lowerPad", "lowerPad", 0, 0, 1, 0.33)
    upperPad.SetBottomMargin(0.00001)
    upperPad.SetLeftMargin(0.15)
    upperPad.SetRightMargin(0.2)
    upperPad.SetBorderMode(0)
    upperPad.SetTopMargin(0.15)
    if log_scale:
        upperPad.SetLogy()
    if xaxis_log_scale==1:
        upperPad.SetLogx()
        lowerPad.SetLogx()
    lowerPad.SetTopMargin(0.00001)
    lowerPad.SetBottomMargin(0.4)
    lowerPad.SetLeftMargin(0.15)
    lowerPad.SetRightMargin(0.2)
    lowerPad.SetBorderMode(0)
    upperPad.Draw()
    lowerPad.Draw()

    upperPad.cd()

    hist_stack.SetMaximum(hist_stack.GetMaximum()*100)#
    hist_stack.SetMinimum(1)
    hist_stack.Draw("HIST")

    hist_stack.GetYaxis().SetTitle("Events")
    hist_stack.GetYaxis().SetTitleOffset(0.8)

    hist_data.SetLineColor(ROOT.kBlack)
    hist_data.Draw("SAME P")
    legend.AddEntry(hist_data, "data","p")

    legend.Draw("SAME")

    upperPad.Modified()

    lowerPad.cd()

    axis = hist_stack.GetStack().Last().Clone("axis")
    axis.SetMinimum(0.0)
    axis.SetMaximum(2.0)
    axis.GetXaxis().SetTitle(var_name)
    axis.GetXaxis().SetTitleOffset(2.5)
    axis.GetYaxis().SetTitle("Data/MC")
    axis.Draw("AXIS")

    line = ROOT.TLine(var_min, 1, var_max, 1)
    line.Draw("SAME")

    for ibin in range(hist_stack.GetHistogram().GetNbinsX()):
        c = hist_stack.GetStack().Last().GetBinCenter(ibin+1)
        w = hist_stack.GetStack().Last().GetBinWidth(ibin+1)
        m = hist_stack.GetStack().Last().GetBinContent(ibin+1)
        #print c, w, m
        if m > 0.0:
            #h = min(hist_stack.GetStack().Last().GetBinError(ibin+1)/m, 0.399)
            h1 = hist_stack.GetStack().Last().GetBinError(ibin+1)/m
            print(f'bin {ibin}: {h1}')
            h = min(h1, 1.0)
            #print h
            box = ROOT.TBox(c-0.5*w, 1-h, c+0.5*w, 1+h)
            box.SetFillStyle(3345)
            box.SetLineColor(ROOT.kGray+1)
            box.SetFillColor(ROOT.kGray)
            style.rootObj.append(box)
            box.Draw("SameF")
            box2 = ROOT.TBox(c-0.5*w, 1-h, c+0.5*w, 1+h)
            box2.SetFillStyle(0)
            box2.SetLineColor(ROOT.kGray+1)
            box2.SetFillColor(ROOT.kGray)
            style.rootObj.append(box2)
            box2.Draw("SameL")

    hist_ratio = hist_data.Clone("ratio histogram"+category+dilep_charge+var)
    hist_ratio.Divide(hist_stack.GetStack().Last())
    hist_ratio.Draw("SAME P")

    cv.cd()

    style.makeLumiText(0.85, 0.97, lumi=str(luminosity), year=str(year))
    style.makeCMSText(0.13, 0.97,additionalText="Preliminary", dx=0.1)

    if CR == "CR_DY":
        control_reg = "DY+jets CR"
    elif CR == "CR_resolved":
        control_reg = "resolved CR"
    elif CR == "CR_mass_sideband_high":
        control_reg = "high mass sideband CR"
    elif CR == "CR_mass_sideband_low":
        control_reg = "low mass sideband CR"
    elif CR == "CR_3lep":
        control_reg = "3 lepton CR"
    else:
        control_reg = CR

    if plotting_word_latex!="":
        #style.makeText(0.2, 0.80, 0.3, 0.90, category_name+", "+dilep_charge+", "+control_reg+", "+plotting_word_latex)
        style.makeText(0.2, 0.80, 0.3, 0.90, category_name+", "+dilep_charge+", "+control_reg)
        style.makeText(0.2, 0.72, 0.3, 0.82, plotting_word_latex)
    else:
        style.makeText(0.2, 0.80, 0.3, 0.90, category_name+", "+dilep_charge+", "+control_reg)

    if hist_stack.GetStack().Last().Integral() > 0:
        style.makeText(0.15, 0.03, 0.4, 0.08, "data/MC = {0:.3g}".format(hist_data.Integral()/hist_stack.GetStack().Last().Integral()))
    else:
        style.makeText(0.15, 0.03, 0.4, 0.08, "data/MC = N/A")

    cv.Modified()

    if plotting_word=="":
        cv.SaveAs(os.path.join(plotdir, str(var)+"_"+str(CR)+"_"+str(category)+"_"+str(dilep_charge)+"_"+str(year)+".pdf"))
        cv.SaveAs(os.path.join(plotdir, str(var)+"_"+str(CR)+"_"+str(category)+"_"+str(dilep_charge)+"_"+str(year)+".png"))
    else:
        cv.SaveAs(os.path.join(plotdir, str(var)+"_"+str(CR)+"_"+str(category)+"_"+str(dilep_charge)+"_"+str(year)+"_"+str(plotting_word)+".pdf"))
        cv.SaveAs(os.path.join(plotdir, str(var)+"_"+str(CR)+"_"+str(category)+"_"+str(dilep_charge)+"_"+str(year)+"_"+str(plotting_word)+".png"))

    cv.Close()
