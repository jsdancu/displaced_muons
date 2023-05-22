# Script containing various methods for statistical calculations 
# Credits go to Mikael Mieskolainen, Henrikas Svidras, Vilius Cepaitis

# j.dancu18@imperial.ac.uk, 2022

import os
import uproot
import pandas as pd
import numpy as np
import scipy.stats
from matplotlib import pyplot
import pickle
import json
import seaborn as sns
import argparse
import time
from datetime import date
from pathlib import Path
from termcolor import colored, cprint

def exact_CI(k, n, conf=0.683):
    """ calculated clopper pearson confidence intervals

    Args:
        k: number of passes ("numerator of proportion")
        n: number of trials ("denominator of proportion")
        conf: (optional) confidence level, by default 0.683

    Returns:
        ratio, [lower ratio error, upper ratio error]
    """

    assert k <= n, f"denominator {n} found to be smaller than numerator {k} when calculating interval."

    k = float(k)
    n = float(n)

    alpha = (1 - conf)
    # if (0,0) is passed, then there is no need to perform any calculation.
    if n == 0:
        p = 0
        up = 0
        down = 0
    else:
        p = k/n
        up = 1 if k == n else 1-scipy.stats.beta.ppf(alpha/2, n-k, k+1)
        down = 0 if k == 0 else 1-scipy.stats.beta.ppf(1-alpha/2, n-k+1, k)

    result = (p, p-down, up-p)

    return result

def clopper_pearson_err(k, n, CL=[0.025, 0.975]):
    """ Clopper-Pearson binomial proportion confidence interval.
    Below, beta.ppf (percent point functions) returns inverse CDF for quantiles.
    
    Args:
        k  : observed success counts
        n  : number of trials
        CL : confidence levels
    Returns:
        corresponding interval points
    """
    # Special case must be treated separately
    if n == 0:
        lower = 0
        upper = 0

    elif   k == 0:
        lower = 0
        upper = 1 - (1-CL[1])**(1/n)

    # Special case must be treated separately
    elif k == n:
        lower = CL[0]**(1/n)
        upper = 1

    # Normal case
    else:
        lower = scipy.stats.beta.ppf(q=CL[0], a=k, b=n-k+1)
        upper = scipy.stats.beta.ppf(q=CL[1], a=k+1, b=n-k)

    return np.array([lower, upper])

def katz_binomial_ratio_err(k1,n1, k2,n2, z=1.0):
    """ Katz el al. ratio confidence interval of two binomial proportions.
    """

    RR      = (k1/n1) / (k2/n2) # mean

    logRR   = np.log(RR)
    seLogRR = np.sqrt(1/k1 + 1/k2 - 1/n1 - 1/n2)

    lower   = np.exp(logRR - z*seLogRR)
    upper   = np.exp(logRR + z*seLogRR)

    return np.array([lower, upper])

def error_propagation(y, x1, x1_err, x2, x2_err):
    """Simple error propagation formula assuming uncorrelated variables"""

    y_err1 = np.divide(x1_err, x1, out=np.zeros_like(x1_err), where=x1!=0)
    y_err2 = np.divide(x2_err, x2, out=np.zeros_like(x2_err), where=x2!=0)
    y_err = y * np.sqrt((y_err1)**2 + (y_err2)**2)
    return y_err

def divide_efficiency(n_nom, n_denom, confidence=0.683):
    """ divides two histograms for an efficiency calculation

    Args:
        n_nom: y values of nominator histogram (1d or 2d)
        n_denom: y values of denominator histogram  (1d or 2d)
        confidence: (optional) confidence level, by default 0.683
        scale: (optional, float) if your n_nom and n_denom have previously been scaled, you specify the scale here.

    Returns:
        ratio, [lower ratio error, upper ratio error]
    """
    shape = np.shape(n_nom)

    flattened_n_nom = n_nom.flatten()
    flattened_n_denom = n_denom.flatten()

    rat = []
    err_down = []
    err_up = []

    alpha = (1.0 - confidence)
    CL = [alpha/2.0, 1.0 - alpha/2.0]

    for passes, counts in zip(flattened_n_nom, flattened_n_denom):
        bin_ratio = clopper_pearson_err(passes, counts, CL=CL)
        if (counts > 0) and (passes <= counts):
            rat.append(passes/counts)
        else:
            rat.append(0)
        err_down.append(bin_ratio[0])
        err_up.append(bin_ratio[1])

    err_down = np.array(rat) - np.array(err_down)
    err_up = np.array(err_up) - np.array(rat)

    err_down = np.reshape(err_down, shape)
    err_up = np.reshape(err_up, shape)
    rat = np.reshape(rat, shape)

    return rat, (err_down, err_up)

def binomial_ratio(k1, n1, k2, n2, confidence=0.683, method="CP"):
    """ ratio between two confidence intervals and the corresponding confidence interval
    
    Args:
        k1, n1 : for the binomial 1
        k2, n2 : for the binomial 2
        confidence: (optional) confidence level, by default 0.683
    
    (k1/n1) / (k2/n2)
    Based on the conditional odds ratio ~ p1 / (p1 + p2)
    """
    shape = np.shape(k1)
    
    flattened_k1 = k1.flatten()
    flattened_n1 = n1.flatten()
    flattened_k2 = k2.flatten()
    flattened_n2 = n2.flatten()

    rat = []
    err_down = []
    err_up = []

    alpha = (1.0 - confidence)
    CL = [alpha/2.0, 1.0 - alpha/2.0]

    print("binomial ratio method: ", method)

    for passes1, passes2, counts1, counts2 in zip(flattened_k1, flattened_k2, flattened_n1, flattened_n2):
        if ((counts1 > 0) and (passes1 <= counts1)) and ((counts2 > 0) and (passes2 <= counts2)):
            eff1 = passes1/counts1 #efficiency for data
            eff2 = passes2/counts2 #efficiency for MC
            sf = eff1 / eff2
            rat.append(sf)
            if method == "CP":
                bin_ratio = clopper_pearson_err(passes1, passes1+passes2, CL=CL)
                err_down.append(counts2/counts1 * bin_ratio[0]/(1.0 - bin_ratio[0]))
                err_up.append(counts2/counts1 * bin_ratio[1]/(1.0 - bin_ratio[1]))
            elif method == "Katz":
                bin_ratio = katz_binomial_ratio_err(passes1, counts1, passes2, counts2, z=1.0)
                err_down.append(bin_ratio[0])
                err_up.append(bin_ratio[1])
            elif method == "error_prop":
                bin_ratio1 = clopper_pearson_err(passes1, counts1, CL=CL)
                bin_ratio2 = clopper_pearson_err(passes2, counts2, CL=CL)
                err_down.append(error_propagation(sf, eff1, eff1-bin_ratio1[0], eff2, eff2-bin_ratio2[0]))
                err_up.append(error_propagation(sf, eff1, bin_ratio1[1]-eff1, eff2, bin_ratio2[1]-eff2))
            else:
                raise Exception(f'binomial_ratio: unknown basic method = {method}')
        else:
            rat.append(0)
            err_down.append(0)
            err_up.append(0)
    if (method == "CP") or (method == "Katz"):
        err_down = np.array(rat) - np.array(err_down)
        err_up = np.array(err_up) - np.array(rat)

    err_down = np.reshape(err_down, shape)
    err_up = np.reshape(err_up, shape)
    rat = np.reshape(rat, shape)

    return rat, (err_down, err_up)