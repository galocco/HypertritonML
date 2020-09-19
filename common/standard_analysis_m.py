#!/usr/bin/env python3
import argparse
import os
import time
import warnings

import numpy as np
import yaml

import hyp_analysis_utils as hau
import pandas as pd
import xgboost as xgb
from analysis_classes import (ModelApplication, TrainingAnalysis)
from hipe4ml import analysis_utils, plot_utils
from hipe4ml.model_handler import ModelHandler
from ROOT import TFile, gROOT

# avoid pandas warning
warnings.simplefilter(action='ignore', category=FutureWarning)
gROOT.SetBatch()

###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-split', '--split', help='Run with matter and anti-matter splitted', action='store_true')
parser.add_argument('config', help='Path to the YAML configuration file')
args = parser.parse_args()

with open(os.path.expandvars(args.config), 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
###############################################################################

###############################################################################
# define analysis global variables
N_BODY = params['NBODY']
FILE_PREFIX = params['FILE_PREFIX']
LARGE_DATA = params['LARGE_DATA']
LOAD_LARGE_DATA = params['LOAD_LARGE_DATA']

CENT_CLASSES = params['CENTRALITY_CLASS']
PT_BINS = params['PT_BINS']
CT_BINS = params['CT_BINS']
COLUMNS = params['TRAINING_COLUMNS']

SPLIT_MODE = args.split

if SPLIT_MODE:
    SPLIT_LIST = ['_antimatter','_matter']
else:
    SPLIT_LIST = ['']

###############################################################################
# define paths for loading data
signal_path = os.path.expandvars(params['MC_PATH'])
bkg_path = os.path.expandvars(params['BKG_PATH'])
data_path = os.path.expandvars(params['DATA_PATH'])
analysis_res_path = os.path.expandvars(params['ANALYSIS_RESULTS_PATH'])

BKG_MODELS =  ['expo', 'pol1', 'pol2']
results_dir = os.environ['HYPERML_RESULTS_{}'.format(N_BODY)]

###############################################################################
start_time = time.time()                          # for performances evaluation

file_name = results_dir + f'/{FILE_PREFIX}_std_results.root'
results_file = TFile(file_name, 'recreate')
standard_selection = 'V0CosPA > 0.9999 and NpidClustersHe3 > 80 and He3ProngPt > 1.8 and pt > 2 and pt < 10 and PiProngPt > 0.15 and He3ProngPvDCA > 0.05 and PiProngPvDCA > 0.2 and TPCnSigmaHe3 < 3.5 and TPCnSigmaHe3 > -3.5 and ProngsDCA < 1'
my_selection = "and NitsClustersHe3>1"
if(N_BODY==3):
    application_columns = ['score', 'm', 'ct', 'pt', 'centrality', 'positive', 'mppi_vert', "r", "pz",'ArmenterosAlpha']
else:
    application_columns = ['NitsClustersHe3','score', 'm', 'ct', 'pt', 'centrality','ArmenterosAlpha','V0CosPA','NpidClustersHe3','He3ProngPt','PiProngPt','He3ProngPvDCA','PiProngPvDCA','TPCnSigmaHe3','ProngsDCA']

for split in SPLIT_LIST:

    if LARGE_DATA:
        if LOAD_LARGE_DATA:
            try:
                df_skimmed = pd.read_parquet(os.path.dirname(data_path) + '/skimmed_df.parquet.gzip')
            except:
                print("you need to run the training")
        else:
            df_skimmed = hau.get_skimmed_large_data(data_path, CENT_CLASSES, PT_BINS, CT_BINS, COLUMNS, application_columns, N_BODY, split)
            try:
                df_skimmed.to_parquet(os.path.dirname(data_path) + '/skimmed_df.parquet.gzip', compression='gzip')
            except:
                print("you need to run the training")
        ml_application = ModelApplication(N_BODY, data_path, analysis_res_path, CENT_CLASSES, split, df_skimmed)

    else:
        ml_application = ModelApplication(N_BODY, data_path, analysis_res_path, CENT_CLASSES, split)

    for cclass in CENT_CLASSES:
        cent_dir = results_file.mkdir(f'{cclass[0]}-{cclass[1]}{split}')
        #ml_application = ModelApplication(N_BODY,f data_path, analysis_res_path, CENT_CLASSES, split)
        df_applied = ml_application.df_data.query(standard_selection)
        for ptbin in zip(PT_BINS[:-1], PT_BINS[1:]):
            for ctbin in zip(CT_BINS[:-1], CT_BINS[1:]):
                
                sub_dir = cent_dir.mkdir(f'ct_{ctbin[0]}{ctbin[1]}') if 'ct' in FILE_PREFIX else cent_dir.mkdir(f'pt_{ptbin[0]}{ptbin[1]}')
                sub_dir.cd()
                mass_bins = 40 if ctbin[1] < 16 else 36
                mass_array = np.array(df_applied.query("ct<@ctbin[1] and ct>@ctbin[0] and pt<@ptbin[1] and pt>@ptbin[0]")['m'].values, dtype=np.float64)
                #mass_array = np.array(df_applied['m'].values, dtype=np.float64)
                counts, _ = np.histogram(mass_array, bins=mass_bins, range=[2.96, 3.05])
                print(counts)
                h1_minv = hau.h1_invmass(counts, cclass, ptbin, ctbin, bins=mass_bins, name="")

                for bkgmodel in BKG_MODELS:
                    print(bkgmodel)
                    # create dirs for models
                    fit_dir = sub_dir.mkdir(bkgmodel)
                    fit_dir.cd()

                #h1_minv.Write()
                    rawcounts, err_rawcounts, significance, err_significance, mu, mu_err, _, _ = hau.fit_hist(h1_minv, cclass, ptbin, ctbin, model=bkgmodel, fixsigma=params['SIGMA_MC'] , mode=N_BODY)#, split=split)
                    print(bkgmodel)
                    if bkgmodel == "pol2":
                        print("ct:",ctbin,"cm pT:",ptbin,"GeV/c")
                        print("mu: ",mu*1000,"+-",mu_err*1000,"MeV/c^2")
                        print("B: ",1875.61294257+1115.683-(mu*1000),"+-",mu_err*1000,"MeV")

    


print(f'--- analysis time: {((time.time() - start_time) / 60):.2f} minutes ---')
