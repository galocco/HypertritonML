#!/usr/bin/env python3
import argparse
# import collections.abc
import os
import time
import warnings

import numpy as np
import yaml

import hyp_analysis_utils as hau


from ROOT import TH1D, TFile, gROOT

# avoid pandas warning
warnings.simplefilter(action='ignore', category=FutureWarning)
gROOT.SetBatch()

###############################################################################
parser = argparse.ArgumentParser()
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

CENT_CLASSES = params['CENTRALITY_CLASS']
PT_BINS = params['PT_BINS']
CT_BINS = params['CT_BINS']

EFF_MIN, EFF_MAX, EFF_STEP = params['BDT_EFFICIENCY']
FIX_EFF_ARRAY = np.arange(EFF_MIN, EFF_MAX, EFF_STEP)


###############################################################################
# define paths for loading results
results_dir = os.environ['HYPERML_RESULTS_{}'.format(N_BODY)]

input_file_name = results_dir + f'/{FILE_PREFIX}_results.root'
input_file = TFile(input_file_name, 'read')

output_file_name = results_dir + f'/{FILE_PREFIX}_results_merged.root'
output_file = TFile(output_file_name, 'recreate')


###############################################################################
for cclass in CENT_CLASSES:
    cent_dir_name = f'{cclass[0]}-{cclass[1]}'
    cent_dir = output_file.mkdir(cent_dir_name)
    cent_dir.cd()
    
    PreselEff_matter = input_file.Get(f'{cent_dir_name}_matter/PreselEff')
    #print(f'{cent_dir_name}_matter/PreselEff')
    #BDTeff_matter = input_file.Get(f'{cent_dir_name}_matter/BDTeff')
    PreselEff_antimatter = input_file.Get(f'{cent_dir_name}_antimatter/PreselEff')
    #BDTeff_antimatter = input_file.Get(f'{cent_dir_name}_antimatter/BDTeff')

    #PreselEff_matter.Add(PreselEff_antimatter)
    #BDTeff_matter.Add(BDTeff_antimatter)
    #PreselEff_matter.Scale(1./2)
    #BDTeff_matter.Scale(1./2)

    #BDTeff_matter.Write()
    PreselEff_matter.Write()

    
    for ptbin in zip(PT_BINS[:-1], PT_BINS[1:]):
        for ctbin in zip(CT_BINS[:-1], CT_BINS[1:]):
            # get the dir where the inv mass histo are
            subdir_name = f'ct_{ctbin[0]}{ctbin[1]}' if 'ct' in FILE_PREFIX else f'pt_{ptbin[0]}{ptbin[1]}'
            input_subdir_matter = input_file.Get(f'{cent_dir_name}_matter/{subdir_name}')
            input_subdir_antimatter = input_file.Get(f'{cent_dir_name}_antimatter/{subdir_name}')

            # create the subdir in the output file
            output_subdir = cent_dir.mkdir(subdir_name)
            output_subdir.cd()

            # loop over all the histo in the dir
            hist_counter = 0
            for key_m in input_subdir_matter.GetListOfKeys():
                hist = TH1D(key_m.ReadObj())
                for key_a in input_subdir_antimatter.GetListOfKeys():
                    hist_anti = TH1D(key_a.ReadObj())#TH1D(input_subdir_antimatter.GetListOfKeys()[hist_counter].ReadObj())
                    if hist.GetName()==hist_anti.GetName():
                        hist.Add(hist_anti)
                        hist_counter = hist_counter + 1
                        hist.Write()
                        print(hist.GetName())
                        print(hist_anti.GetName())
                        break
                
output_file.Close()
