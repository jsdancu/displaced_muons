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
import re

import histogramming

parser = argparse.ArgumentParser()

parser.add_argument('--year', dest='year', action='store', default="2016")
parser.add_argument('--var1', dest='var1', action='store', default="absdxy_hack")
parser.add_argument('--var2', dest='var2', action='store', default="pt")
parser.add_argument('--var_latex1', dest='var_latex1', action='store', default="|d_{xy}| (cm)")
parser.add_argument('--var_latex2', dest='var_latex2', action='store', default="p_{T} (GeV)")
parser.add_argument('--logscale1', dest='logscale1', action='store', type=int, default=1)
parser.add_argument('--logscale2', dest='logscale2', action='store', type=int, default=0)
#parser.add_argument('--input_file', dest='input_file', action='store',default="displaced_muons/peakfit/Run2016/Efficiency/NUM_LooseID_DEN_TrackerMuons_absdxy_1_pt_1.pkl")

args = parser.parse_args()

path_plots = "displaced_muons/muon_reco_sf_syst"

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

sf_dict = pd.read_pickle(os.path.join(path_plots, "peakfit_systematics_YEAR_"+args.year+".pkl"))
#print("sf: ", sf)

scale_factors = histogramming.define_plots_2D(bins1, bins2, args.var_latex1, args.var_latex2, "scale_factors_"+args.year, 1, 1)
scale_factors.GetZaxis().SetTitle("Scale factor")
scale_factors.GetZaxis().SetTitleOffset(0.8)
sf_vals = np.zeros((nbins1, nbins2))
sf_err_vals = np.zeros((nbins1, nbins2))
sf_err_stats = np.zeros((nbins1, nbins2))
sf_err_syst = np.zeros((nbins1, nbins2))

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

syst_unc_dict = {}

for sf_bin in sf_dict.keys():

    pattern = r"\d+"
    bins = re.findall(pattern, sf_bin)

    sf = sf_dict[sf_bin]['DEFAULT']['scale']
    sf_stats = sf_dict[sf_bin]['DEFAULT']['scale_err']

    sf_syst_array = {}

    eff_data[int(bins[0])-1, int(bins[1])-1] = sf_dict[sf_bin]['DEFAULT']['eff'][run]
    eff_data_err[int(bins[0])-1, int(bins[1])-1] = sf_dict[sf_bin]['DEFAULT']['eff_err'][run]
    eff_MC[int(bins[0])-1, int(bins[1])-1] = sf_dict[sf_bin]['DEFAULT']['eff']['JPsi_pythia8']
    eff_MC_err[int(bins[0])-1, int(bins[1])-1] = sf_dict[sf_bin]['DEFAULT']['eff_err']['JPsi_pythia8']

    efficiency_data.SetBinContent(int(bins[0]), int(bins[1]), sf_dict[sf_bin]['DEFAULT']['eff'][run])
    efficiency_MC.SetBinContent(int(bins[0]), int(bins[1]), sf_dict[sf_bin]['DEFAULT']['eff']['JPsi_pythia8'])
    efficiency_data.SetBinError(int(bins[0]), int(bins[1]), sf_dict[sf_bin]['DEFAULT']['eff_err'][run])
    efficiency_MC.SetBinError(int(bins[0]), int(bins[1]), sf_dict[sf_bin]['DEFAULT']['eff_err']['JPsi_pythia8'])

    for SYST in ['SYMMETRIC-SIGNAL', 'MASS-WINDOW-DOWN', 'MASS-WINDOW-UP']:

        syst_unc = np.abs(sf_dict[sf_bin][SYST]['scale'] - sf) - sf_stats

        if syst_unc > 0.0:
            sf_syst_array[SYST] = syst_unc
        else:
            sf_syst_array[SYST] = 0.0

    sf_err = np.sqrt(sf_stats**2 + np.sum([syst_val ** 2 for syst_val in sf_syst_array.values()]))

    scale_factors.SetBinContent(int(bins[0]), int(bins[1]), sf)
    scale_factors.SetBinError(int(bins[0]), int(bins[1]), sf_err)

    sf_vals[int(bins[0])-1, int(bins[1])-1] = sf
    sf_err_vals[int(bins[0])-1, int(bins[1])-1] = sf_err
    sf_err_stats[int(bins[0])-1, int(bins[1])-1] = sf_stats
    sf_err_syst[int(bins[0])-1, int(bins[1])-1] = np.sqrt(np.sum([syst_val ** 2 for syst_val in sf_syst_array.values()]))

    syst_unc_dict[sf_bin] = sf_syst_array

    print("sf_bin: ", sf_bin)
    print("sf_syst_array: ", sf_syst_array)

print("sf_vals: ", sf_vals)
print("sf_err_stats: ", sf_err_stats)
print("sf_err_syst: ", sf_err_syst)
  
########### Save separate syst unc values for each bin ##########
#syst_unc_dict.to_pickle(os.path.join(path_plots, str(args.year)+"/systematic_uncertainties_"+str(args.year)+".pkl"))
with open(os.path.join(path_plots, str(args.year)+"/systematic_uncertainties_"+str(args.year)+".pkl"), 'wb') as handle:
    pickle.dump(syst_unc_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

################### Extrapolating manually low0-stats bins  ############# BSSSSSSSSSSS ###############

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
# print("sf_vals: ", sf_vals[3, 4])
# print("sf_err_vals: ", sf_err_vals[3, 4])


histogramming.plotting_scale_factor_2D(scale_factors, sf_vals, sf_err_stats, sf_err_syst, bins1, bins2, path_plots+str("/")+str(args.year), 
"scale_factor2D_NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year)+"_syst_unc", 
str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

# histogramming.plotting_scale_factor_2D(efficiency_data, eff_data, eff_data_err, [], bins1, bins2, path_plots+str("/")+str(args.year), 
# "efficiency_data2D_NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year), 
# str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

# histogramming.plotting_scale_factor_2D(efficiency_MC, eff_MC, eff_MC_err, [], bins1, bins2, path_plots+str("/")+str(args.year), 
# "efficiency_MC2D_NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year), 
# str(args.var1)+"_"+str(args.var2), "TnP_JPsi", "", "", "", args.year, luminosities[args.year], args.logscale1, args.logscale2)

####### Save scale factor historgram in root file ###############
f = ROOT.TFile(os.path.join(path_plots, str(args.year), "scale_factor2D_NUM_MyNum_DEN_MyDen_"+str(args.var1)+"_"+str(args.var2)+"_TnP_"+str(args.year)+"_syst.root"),"RECREATE")
scale_factors.Write()
f.Close()