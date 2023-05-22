import sys
import os
import ROOT
ROOT.gROOT.SetBatch(1)
import style
#ROOT.gStyle.SetPalette(ROOT.kCherry)
style.makeColorTable()
#ROOT.gStyle.SetPalette(256, "acton.txt")
# ROOT.TColor.CreateColorTableFromFile("acton.txt")
# ROOT.gStyle.SetNumberContours(256)
ROOT.TColor.InvertPalette()
from array import array
import random
import numpy as np

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

def plotting_2D(array1, array2, hist1, puWeight=None):

    ntimes = len(array1)

    if puWeight is None:
        puWeight = np.ones(ntimes)

    hist1.FillN(ntimes, array('d', array1), array('d', array2), array('d', puWeight))

    return hist1

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
    hist.SetMaximum(1.5)
    hist.SetMinimum(0.0)
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
    if hist_err_up==[]:
        for i, x_val in enumerate(x):
            for j, y_val in enumerate(y):
                text.DrawLatex(x_val,y_val,f'{hist_val[i,j]:.2f}#pm{hist_err_down[i,j]:.2f}')
                #change bin error value
                #hist.SetBinError(i+1,j+1,hist_err_down[i,j])
    else:
        for i, x_val in enumerate(x):
            for j, y_val in enumerate(y):
                text.DrawLatex(x_val,y_val,f'{hist_val[i,j]:.2f}^{{+{hist_err_up[i,j]:.2f}}}_{{-{hist_err_down[i,j]:.2f}}}')
                #change bin error value
                #hist.SetBinError(i+1,j+1,(hist_err_down[i,j]+hist_err_up[i,j])/2.0)

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

    cv.Close()

    