import uproot
import pandas as pd
from ROOT import TF1, TH1, TH1D, TH2D, TFile, gDirectory
import matplotlib.pyplot as plt
import numpy as np
import math
import yaml
import argparse
import os
from array import*


parser = argparse.ArgumentParser()
parser.add_argument("config", help="Path to the YAML configuration file")
parser.add_argument('-split', '--split', help='Run with matter and anti-matter splitted', action='store_true')
parser.add_argument('-s', '--std', help='Run with matter and anti-matter splitted', action='store_true')
args = parser.parse_args()


if args.split:
    SPLIT_LIST = ['_matter', '_antimatter']
else:
    SPLIT_LIST = ['']

with open(os.path.expandvars(args.config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

def h1_mass_res(counts, cent_class, pt_range, ct_range, bins=45, name=''):
    th1 = TH1D(f'ct{ct_range[0]}{ct_range[1]}_pT{pt_range[0]}{pt_range[1]}_cen{cent_class[0]}{cent_class[1]}_{name}', '', 100, -0.03, 0.03)

    for index in range(0, len(counts)):
        th1.SetBinContent(index+1, counts[index])
        th1.SetBinError(index + 1, math.sqrt(counts[index]))

    th1.SetDirectory(0)

    return th1

hyp3mass = 2.99131
df = uproot.open(os.path.expandvars(params['MC_PATH']))["SignalTable"].pandas.df()
df['mass_shift'] = df.apply(lambda row: row.m - hyp3mass, axis = 1)

ct_bins = params['CT_BINS']
resultsSysDir = os.environ['HYPERML_RESULTS_{}'.format(params['NBODY'])]
file_name =  resultsSysDir + '/' + params['FILE_PREFIX'] + '_mass_res.root'
results_file = TFile(file_name,"recreate")

gauss = TF1("gauss","gaus",-0.002,0.002)
standard_selection = 'V0CosPA > 0.9999 and NpidClustersHe3 > 80 and He3ProngPt > 1.8 and pt > 2 and pt < 10 and PiProngPt > 0.15 and He3ProngPvDCA > 0.05 and PiProngPvDCA > 0.2 and TPCnSigmaHe3 < 3.5 and TPCnSigmaHe3 > -3.5 and ProngsDCA < 1'
df = df.query(standard_selection)
binning = array('d',ct_bins)
for split in SPLIT_LIST:
    hist_mean = TH1D('hist_mean'+split, ';ct (cm); mean[m_{rec}-m_{gen}] (MeV/c^{2})',len(ct_bins)-1,binning)
    fit_mean = TH1D('fit_mean'+split, ';ct (cm); #mu (MeV/c^{2})',len(ct_bins)-1,binning)
    hist_mean.SetMarkerStyle(8)
    fit_mean.SetMarkerStyle(8)

    if split == '_antimatter':
        hist_mean.SetTitle('{}^{3}_{#bar{#Lambda}} #bar{H} mass shift')
        fit_mean.SetTitle('{}^{3}_{#bar{#Lambda}} #bar{H} mass shift')
    else:
        hist_mean.SetTitle('{}^{3}_{#Lambda}H mass shift')
        fit_mean.SetTitle('{}^{3}_{#Lambda}H mass shift')

    results_file.mkdir(split[1:])
    results_file.cd(split[1:])
    for bin in range(1,len(ct_bins)):
        if split=='':
            dfq = df.query('ct>@ct_bins[@bin-1] and ct<@ct_bins[@bin]')
        elif split=='_antimatter':
            dfq = df.query('ct>@ct_bins[@bin-1] and ct<@ct_bins[@bin] and Matter < 0.5')
        else:
            dfq = df.query('ct>@ct_bins[@bin-1] and ct<@ct_bins[@bin] and Matter > 0.5')
        counts, _ = np.histogram(dfq['mass_shift'], bins=100, range=[-0.03, 0.03])
        histo = h1_mass_res(counts,[0,90],[0,10],[ct_bins[bin-1],ct_bins[bin]])
        histo.Fit(gauss,"QR")
        fit_mean.SetBinContent(bin, gauss.GetParameter(1)*1000)
        fit_mean.SetBinError(bin, gauss.GetParError(1)*1000)
        hist_mean.SetBinContent(bin, histo.GetMean()*1000)
        hist_mean.SetBinError(bin, histo.GetMeanError()*1000)
        histo.Write()
    results_file.cd()
    fit_mean.Write()
    hist_mean.Write()