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

#path_plots = "displaced_muons/peakfit/Run"+str(args.year)+"/Efficiency"
#path_plots = "icenet/output_hack/peakfit/Run"+str(args.year)+"/Efficiency"
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

for SYST in ['Nominal', 'nVertices_Up', 'nVertices_Down']:

    scale_factors = histogramming.define_plots_2D(bins1, bins2, args.var_latex1, args.var_latex2, "scale_factors_"+args.year, 1, 1)
    scale_factors.GetZaxis().SetTitle("Scale factor")
    scale_factors.GetZaxis().SetTitleOffset(0.8)
    sf_vals = np.zeros((nbins1, nbins2))
    sf_err_vals = np.zeros((nbins1, nbins2))

    efficiency_data = histogramming.define_plots_2D(bins1, bins2, args.var_latex1, args.var_latex2, "efficiency data "+args.year, 1, 1)
    efficiency_data.GetZaxis().SetTitle("#varepsilon(data)")
    efficiency_data.GetZaxis().SetTitleOffset(0.8)
    efficiency_MC = histogramming.define_plots_2D(bins1, bins2, args.var_latex1, args.var_latex2, "efficiency MC "+args.year, 1, 1)
    efficiency_MC.GetZaxis().SetTitle("#varepsilon(MC)")
    efficiency_MC.GetZaxis().SetTitleOffset(0.8)

    eff_data = np.zeros((nbins1, nbins2))
    eff_data_err = np.zeros((nbins1, nbins2))
    eff_MC = np.zeros((nbins1, nbins2))
    eff_MC_err = np.zeros((nbins1, nbins2))

    for i in range(1, int(nbins1)+1):
        for j in range(1, int(nbins2)+1):

            #sf = pd.read_pickle(os.path.join(path_plots, "NUM_LooseID_DEN_TrackerMuons_"+str(args.var1)+"_"+str(i)+"_"+str(args.var2)+"_"+str(j)+".pkl"))
            #sf = pd.read_pickle(os.path.join(path_plots, "NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(i)+"_"+str(args.var2)+"_"+str(j)+".pkl"))
            #sf = pd.read_pickle(os.path.join(path_plots, SYST, "NUM_LooseID_DEN_TrackerMuons_"+str(args.var1)+"_"+str(i)+"_"+str(args.var2)+"_"+str(j)+".pkl"))
            sf = pd.read_pickle(os.path.join(path_plots, SYST, "NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(i)+"_"+str(args.var2)+"_"+str(j)+".pkl"))
            print("i: ", i, ", j: ", j, sf)
            sf_vals[i-1, j-1] = sf['scale']
            sf_err_vals[i-1, j-1] = sf['scale_err']

            eff_data[i-1, j-1] = sf['eff'][run]
            eff_data_err[i-1, j-1] = sf['eff_err'][run]
            eff_MC[i-1, j-1] = sf['eff']['JPsi_pythia8']
            eff_MC_err[i-1, j-1] = sf['eff_err']['JPsi_pythia8']

            scale_factors.SetBinContent(i, j, sf['scale'])
            efficiency_data.SetBinContent(i, j, sf['eff'][run])
            efficiency_MC.SetBinContent(i, j, sf['eff']['JPsi_pythia8'])

            efficiency_data.SetBinError(i, j, sf['eff_err'][run])
            efficiency_MC.SetBinError(i, j, sf['eff_err']['JPsi_pythia8'])
            scale_factors.SetBinError(i, j, sf['scale_err'])

    print("sf_vals: ", sf_vals)
    print("sf_err_vals: ", sf_err_vals)

    ################### Extrapolating manually low0-stats bins  ############# BSSSSSSSSSSS ###############
    if str(args.year) == "2016":

        x1 = scale_factors.GetBinContent(3, 5)
        x2 = scale_factors.GetBinContent(3, 4)
        x3 = scale_factors.GetBinContent(4, 4)
        x1_err = scale_factors.GetBinError(3, 5)
        x2_err = scale_factors.GetBinError(3, 4)
        x3_err = scale_factors.GetBinError(4, 4)
        sf_val = 1./3. * (x1 + x2 + x3)
        sf_err = sf_val/3. * np.sqrt(x1_err**2 + x2_err**2 + x3_err**2)
        scale_factors.SetBinContent(4, 5, sf_val)
        scale_factors.SetBinError(4, 5, sf_err)
        print("sf_val: ", sf_val)
        print("sf_err: ", sf_err)

    elif str(args.year) == "2017":

        x1 = scale_factors.GetBinContent(3, 5)
        x2 = scale_factors.GetBinContent(3, 4)
        x3 = scale_factors.GetBinContent(3, 3)
        x4 = scale_factors.GetBinContent(4, 3)
        x1_err = scale_factors.GetBinError(3, 5)
        x2_err = scale_factors.GetBinError(3, 4)
        x3_err = scale_factors.GetBinError(3, 3)
        x4_err = scale_factors.GetBinError(4, 3)
        sf_val = 1./4. * (x1 + x2 + x3 + x4)
        sf_err = sf_val/4. * np.sqrt(x1_err**2 + x2_err**2 + x3_err**2 + x4_err**2)
        scale_factors.SetBinContent(4, 5, sf_val)
        scale_factors.SetBinError(4, 5, sf_err)
        scale_factors.SetBinContent(4, 4, sf_val)
        scale_factors.SetBinError(4, 4, sf_err)

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
        sf_val = 1./6. * (x1 + x2 + x3 + x4 + x5 + x6)
        sf_err = sf_val/6. * np.sqrt(x1_err**2 + x2_err**2 + x3_err**2 + x4_err**2 + x5_err**2 + x6_err**2)
        scale_factors.SetBinContent(2, 5, sf_val)
        scale_factors.SetBinError(2, 5, sf_err)
        scale_factors.SetBinContent(3, 5, sf_val)
        scale_factors.SetBinError(3, 5, sf_err)

    # histogramming.plotting_scale_factor_2D(scale_factors, sf_vals, sf_err_vals, [], bins1, bins2, path_plots+str("/")+SYST, 
    # "scale_factor2D_NUM_LooseID_DEN_TrackerMuons_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year), 
    # str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

    # histogramming.plotting_scale_factor_2D(efficiency_data, eff_data, eff_data_err, [], bins1, bins2, path_plots+str("/")+SYST, 
    # "efficiency_data2D_NUM_LooseID_DEN_TrackerMuons_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year), 
    # str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

    # histogramming.plotting_scale_factor_2D(efficiency_MC, eff_MC, eff_MC_err, [], bins1, bins2, path_plots+str("/")+SYST, 
    # "efficiency_MC2D_NUM_LooseID_DEN_TrackerMuons_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year), 
    # str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

    # ####### Save scale factor historgram in root file ###############
    # f = ROOT.TFile(os.path.join(path_plots, SYST, "scale_factor2D_NUM_LooseID_DEN_TrackerMuons"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year)+".root"),"RECREATE")
    # scale_factors.Write()
    # f.Close()

    histogramming.plotting_scale_factor_2D(scale_factors, sf_vals, sf_err_vals, [], bins1, bins2, path_plots+str("/")+SYST, 
    "scale_factor2D_NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year), 
    str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

    histogramming.plotting_scale_factor_2D(efficiency_data, eff_data, eff_data_err, [], bins1, bins2, path_plots+str("/")+SYST, 
    "efficiency_data2D_NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year), 
    str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

    histogramming.plotting_scale_factor_2D(efficiency_MC, eff_MC, eff_MC_err, [], bins1, bins2, path_plots+str("/")+SYST, 
    "efficiency_MC2D_NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year), 
    str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

    ####### Save scale factor historgram in root file ###############
    f = ROOT.TFile(os.path.join(path_plots, SYST, "scale_factor2D_NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year)+".root"),"RECREATE")
    scale_factors.Write()
    f.Close()



