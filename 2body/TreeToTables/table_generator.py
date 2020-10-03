#!/usr/bin/env python3

import os
from ROOT import gROOT

gROOT.SetBatch(True)

gROOT.LoadMacro("GenerateTableFromMC.cc")
gROOT.LoadMacro("GenerateTableFromData.cc")
from ROOT import GenerateTableFromMC, GenerateTableFromData

input_dir = os.environ['HYPERML_DATA_2']
output_dir = os.environ['HYPERML_TABLES_2']+"/splines_tables"

print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate Signal Table")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromMC(True, input_dir, output_dir, True)
print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate Data Table")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromData(False, False , input_dir, output_dir, False)
print("++++++++++++++++++++++++++++++++++++++++++")
print("Generate Like-Sign Backgoundd Table")
print("++++++++++++++++++++++++++++++++++++++++++")
GenerateTableFromData(True, False , input_dir ,output_dir, False)