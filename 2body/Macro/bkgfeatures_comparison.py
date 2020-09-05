#!/usr/bin/env python3
import os

import pandas as pd
from hipe4ml import plot_utils as pu
import uproot

import matplotlib.pyplot as plt

data_path = os.path.expandvars('$HYPERML_TABLES_2/splines_tables/DataTable_18.root')
sig_path = os.path.expandvars('$HYPERML_TABLES_2/splines_tables/SignalTable.root')

sig_selection = '2<=pt<=10'
data_selection = '(m<2.98 or m>3.005) and 2<=pt<=10'

df_sig = uproot.open(sig_path)['SignalTable'].pandas.df().query(sig_selection)
df_data = uproot.open(data_path)['DataTable'].pandas.df().query(data_selection)

df_data['y'] = 0
df_sig['y'] = 1

columns = ['V0CosPA'
  ,'pt'
  ,'ProngsDCA'
  ,'PiProngPvDCAXY'
  ,'He3ProngPvDCAXY'
  ,'He3ProngPvDCA'
  ,'PiProngPvDCA'
  ,'NpidClustersHe3'
  ,'TPCnSigmaHe3'
  ,'TPCnSigmaPi'
  ,'NitsClustersHe3']


df = pd.concat([df_sig, df_data])
pu.plot_corr([df_sig, df_data],columns,columns)
#plt.matshow(df_sig[columns].corr())
plt.show()