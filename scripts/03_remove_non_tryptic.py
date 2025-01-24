# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 18:20:27 2025

@author: lafields2
"""

import pandas as pd
import csv

data_file = r"D:\Manuscripts\2024_Haiyan_VitaminA\UniProt_PRIDE_retrival\lip_sites_unfiltered_Lung_w_aa.csv"
output_path = r"D:\Manuscripts\2024_Haiyan_VitaminA\UniProt_PRIDE_retrival"

data = pd.read_csv(data_file)
filtered_data = data[data['Sequence T'].str[-1].isin(['R', 'K'])]

file_out_path = output_path + '\\lip_sites_unfiltered_Lung_w_aa_KRonly.csv'
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        filtered_data.to_csv(filec,index=False)