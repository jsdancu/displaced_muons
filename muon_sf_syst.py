import os
import uproot
import pandas as pd
import numpy as np
from matplotlib import pyplot
import pickle
import json
import argparse
import time
from datetime import date
from array import array
import ROOT

import histogramming

parser = argparse.ArgumentParser()

parser.add_argument('--year', dest='year', action='store', default="2016")
parser.add_argument('--var1', dest='var1', action='store', default="absdxy")
parser.add_argument('--var2', dest='var2', action='store', default="pt")
parser.add_argument('--var_latex1', dest='var_latex1', action='store', default="|d_{xy}| (cm)")
parser.add_argument('--var_latex2', dest='var_latex2', action='store', default="p_{T} (GeV)")
parser.add_argument('--logscale1', dest='logscale1', action='store', type=int, default=1)
parser.add_argument('--logscale2', dest='logscale2', action='store', type=int, default=0)
#parser.add_argument('--input_file', dest='input_file', action='store',default="displaced_muons/peakfit/Run2016/Efficiency/NUM_LooseID_DEN_TrackerMuons_absdxy_1_pt_1.pkl")

args = parser.parse_args()

path_plots = "icenet/output_02Dec22/peakfit/Run"+str(args.year)+"/Efficiency"

binning = {
    "pt": [3.0, 4.0, 6.0, 10.0, 16.0, 30.0],
    "abseta": [0, 0.4, 0.6, 0.9, 1.2, 2.4],
    "eta": [-2.4, -1.2, -0.9, -0.6, -0.4, 0.0, 0.4, 0.6, 0.9, 1.2, 2.4],
    "absdxy": [0.001, 0.1, 1.0, 10.0],
    "absdxy_hack": [0.00001, 0.001, 0.01, 0.1, 10.0],
    "absdxy_sig": [0.001, 3.0, 10.0, 100.0]
    }

luminosities = {
        "2016": 3.9,
        "2017": 2.4,
        "2018": 8.6
        }

run = "Run"+str(args.year)

nbins1 = len(binning[args.var1])-1
nbins2 = len(binning[args.var2])-1
bins1 = array('d', binning[args.var1])
bins2 = array('d', binning[args.var2])

scale_factors = histogramming.define_plots_2D(bins1, bins2, args.var_latex1, args.var_latex2, "scale_factors_"+args.year, 1, 1)
scale_factors.GetZaxis().SetTitle("Scale factor syst")
scale_factors.GetZaxis().SetTitleOffset(0.8)
scale_factors_syst = histogramming.define_plots_2D(bins1, bins2, args.var_latex1, args.var_latex2, "scale_factors_syst_"+args.year, 1, 1)
scale_factors_syst.GetZaxis().SetTitle("Scale factor syst")
scale_factors_syst.GetZaxis().SetTitleOffset(0.8)
sf_vals = np.zeros((nbins1, nbins2))
sf_err_vals = np.zeros((nbins1, nbins2))
sf_syst = np.zeros((nbins1, nbins2))

#for SYST in ['Nominal', 'nVertices_Up', 'nVertices_Down']:
for i in range(1, int(nbins1)+1):
    for j in range(1, int(nbins2)+1):
        
        sf = pd.read_pickle(os.path.join(path_plots, 'Nominal', "NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(i)+"_"+str(args.var2)+"_"+str(j)+".pkl"))
        sf_up = pd.read_pickle(os.path.join(path_plots, 'nVertices_Up', "NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(i)+"_"+str(args.var2)+"_"+str(j)+".pkl"))
        sf_down = pd.read_pickle(os.path.join(path_plots, 'nVertices_Down', "NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(i)+"_"+str(args.var2)+"_"+str(j)+".pkl"))
        # sf = pd.read_pickle(os.path.join(path_plots, 'Nominal', "NUM_LooseID_DEN_TrackerMuons_"+str(args.var1)+"_"+str(i)+"_"+str(args.var2)+"_"+str(j)+".pkl"))
        # sf_up = pd.read_pickle(os.path.join(path_plots, 'nVertices_Up', "NUM_LooseID_DEN_TrackerMuons_"+str(args.var1)+"_"+str(i)+"_"+str(args.var2)+"_"+str(j)+".pkl"))
        # sf_down = pd.read_pickle(os.path.join(path_plots, 'nVertices_Down', "NUM_LooseID_DEN_TrackerMuons_"+str(args.var1)+"_"+str(i)+"_"+str(args.var2)+"_"+str(j)+".pkl"))
        print("i: ", i, ", j: ", j)
        # print("sf_up['scale']: ", sf_up['scale'])
        # print("sf_down['scale']: ", sf_down['scale'])

        scale_factors.SetBinContent(i, j, sf['scale'])

        sf_vals[i-1, j-1] = sf['scale']

        sf_syst[i-1, j-1] = np.abs(sf['scale'] - sf_up['scale'])/2.0 + np.abs(sf['scale'] - sf_down['scale'])/2.0
        scale_factors_syst.SetBinContent(i, j, sf_syst[i-1, j-1])

        sf_err_vals[i-1, j-1] = np.sqrt(sf['scale_err']**2 + sf_syst[i-1, j-1]**2)
        scale_factors.SetBinError(i, j, sf_err_vals[i-1, j-1])

print("sf_syst: ", sf_syst)

################### Extrapolating manually low0-stats bins  ############# BSSSSSSSSSSS ###############
if str(args.year) == "2016":

    x1 = scale_factors.GetBinContent(3, 5)
    x2 = scale_factors.GetBinContent(3, 4)
    x3 = scale_factors.GetBinContent(4, 4)
    x1_err = scale_factors.GetBinError(3, 5)
    x2_err = scale_factors.GetBinError(3, 4)
    x3_err = scale_factors.GetBinError(4, 4)
    sf_vals[3, 4] = 1./3. * (x1 + x2 + x3)
    sf_err_vals[3, 4] = sf_vals[3, 4]/3. * np.sqrt(x1_err**2 + x2_err**2 + x3_err**2)
    scale_factors.SetBinContent(4, 5, sf_vals[3, 4])
    scale_factors.SetBinError(4, 5, sf_err_vals[3, 4])
    print("sf_vals: ", sf_vals[3, 4])
    print("sf_err_vals: ", sf_err_vals[3, 4])

elif str(args.year) == "2017":

    x1 = scale_factors.GetBinContent(3, 5)
    x2 = scale_factors.GetBinContent(3, 4)
    x3 = scale_factors.GetBinContent(3, 3)
    x4 = scale_factors.GetBinContent(4, 3)
    x1_err = scale_factors.GetBinError(3, 5)
    x2_err = scale_factors.GetBinError(3, 4)
    x3_err = scale_factors.GetBinError(3, 3)
    x4_err = scale_factors.GetBinError(4, 3)
    sf_vals[3, 4] = 1./4. * (x1 + x2 + x3 + x4)
    sf_err_vals[3, 4] = sf_vals[3, 4]/4. * np.sqrt(x1_err**2 + x2_err**2 + x3_err**2 + x4_err**2)
    sf_vals[3, 3] = sf_vals[3, 4]
    sf_err_vals[3, 3] = sf_err_vals[3, 4]
    scale_factors.SetBinContent(4, 5, sf_vals[3, 4])
    scale_factors.SetBinError(4, 5, sf_err_vals[3, 4])
    scale_factors.SetBinContent(4, 4, sf_vals[3, 3])
    scale_factors.SetBinError(4, 4, sf_err_vals[3, 3])

elif str(args.year) == "2018":    
    x1 = scale_factors.GetBinContent(1, 5)
    x2 = scale_factors.GetBinContent(1, 4)
    x3 = scale_factors.GetBinContent(2, 4)
    x4 = scale_factors.GetBinContent(3, 4)
    x5 = scale_factors.GetBinContent(4, 4)
    x6 = scale_factors.GetBinContent(4, 5)
    x1_err = scale_factors.GetBinError(1, 5)
    x2_err = scale_factors.GetBinError(1, 4)
    x3_err = scale_factors.GetBinError(2, 4)
    x4_err = scale_factors.GetBinError(3, 4)
    x5_err = scale_factors.GetBinError(4, 4)
    x6_err = scale_factors.GetBinError(4, 5)
    sf_vals[1, 4] = 1./6. * (x1 + x2 + x3 + x4 + x5 + x6)
    sf_err_vals[1, 4] = sf_vals[1, 4]/6. * np.sqrt(x1_err**2 + x2_err**2 + x3_err**2 + x4_err**2 + x5_err**2 + x6_err**2)
    sf_vals[2, 4] = sf_vals[3, 4]
    sf_err_vals[2, 4] = sf_err_vals[3, 4]
    scale_factors.SetBinContent(2, 5, sf_vals[1, 4])
    scale_factors.SetBinError(2, 5, sf_err_vals[1, 4])
    scale_factors.SetBinContent(3, 5, sf_vals[2, 4])
    scale_factors.SetBinError(3, 5, sf_err_vals[2, 4])

histogramming.plotting_scale_factor_2D(scale_factors_syst, sf_syst, [], [], bins1, bins2, path_plots, 
"scale_factor2D_NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year)+"_syst_unc", 
str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

histogramming.plotting_scale_factor_2D(scale_factors, sf_vals, sf_err_vals, [], bins1, bins2, path_plots, 
"scale_factor2D_NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year)+"_syst", 
str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

####### Save scale factor historgram in root file ###############
f = ROOT.TFile(os.path.join(path_plots, "scale_factor2D_NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year)+"_syst.root"),"RECREATE")
scale_factors.Write()
f.Close()

# histogramming.plotting_scale_factor_2D(scale_factors_syst, sf_syst, [], [], bins1, bins2, path_plots, 
# "scale_factor2D_NUM_LooseID_DEN_TrackerMuons_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year)+"_syst_unc", 
# str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

####### Save scale factor historgram in root file ###############
# f = ROOT.TFile(os.path.join(path_plots, "scale_factor2D_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year)+".root"),"RECREATE")
# scale_factors.Write()
# f.Close()



