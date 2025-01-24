# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 09:08:40 2022

@author: lawashburn
"""

import pandas as pd
import csv
import math
import scipy
from scipy.stats import ttest_ind
import numpy as np

protein_path = r"C:\Users\lawashburn\Downloads\Haiyan_LiP-MS-main\Haiyan_LiP-MS-main\Example_Input_Files\proteinGroups_Lung.csv" #Sample names associated with protein will end in T
peptide_path = r"C:\Users\lawashburn\Downloads\Haiyan_LiP-MS-main\Haiyan_LiP-MS-main\Example_Input_Files\peptides_Lung.csv" #Sample names associated with protein will end in P
output_path = r"C:\Users\lawashburn\Downloads\Haiyan_LiP-MS-main\Haiyan_LiP-MS-main\out_pt2"
experimental_sample_name_prefix = 'SRT_L'
control_sample_name_prefix = 'WT_L'

p_cutoff = 0.05
min_FC = 2
min_log_FC = 1

exp_1_name_protein = 'LFQ intensity ' + experimental_sample_name_prefix + '1_T'
exp_2_name_protein = 'LFQ intensity ' + experimental_sample_name_prefix + '2_T'
exp_3_name_protein = 'LFQ intensity ' + experimental_sample_name_prefix + '3_T'

ctrl_1_name_protein = 'LFQ intensity ' + control_sample_name_prefix + '1_T'
ctrl_2_name_protein = 'LFQ intensity ' + control_sample_name_prefix + '2_T'
ctrl_3_name_protein = 'LFQ intensity ' + control_sample_name_prefix + '3_T'

exp_1_name_peptide = 'LFQ intensity ' + experimental_sample_name_prefix + '1_P'
exp_2_name_peptide = 'LFQ intensity ' + experimental_sample_name_prefix + '2_P'
exp_3_name_peptide = 'LFQ intensity ' + experimental_sample_name_prefix + '3_P'

ctrl_1_name_peptide = 'LFQ intensity ' + control_sample_name_prefix + '1_P'
ctrl_2_name_peptide = 'LFQ intensity ' + control_sample_name_prefix + '2_P'
ctrl_3_name_peptide = 'LFQ intensity ' + control_sample_name_prefix + '3_P'

#Determine fold change based on protein level results
protein_report = pd.read_csv(protein_path)
protein_report = protein_report.replace(0, np.NaN) #replace empty values with NaN so mean can be taken without impact from 0s
protein_report['Exp mean protein'] = protein_report[[exp_1_name_protein, exp_2_name_protein,exp_3_name_protein]].mean(axis=1,skipna=True) #the mean is taken without considering 0s
                                                                                                                                        
protein_report['Ctrl mean protein'] = protein_report[[ctrl_1_name_protein, ctrl_2_name_protein,ctrl_3_name_protein]].mean(axis=1,skipna=True) #the mean is taken without considering 0s
                                                                                                                                        
protein_report['protein fold change'] = protein_report['Exp mean protein']/protein_report['Ctrl mean protein']
protein_report['Log2(protein fold change)'] = abs(np.log2(protein_report['protein fold change']))

#Determine p-value between ctrl and exp peptides
peptide_report = pd.read_csv(peptide_path)
merged_pep_prot_report = peptide_report.merge(protein_report,on=['Protein names','Gene names'],how='left')

print('Number of protein entries: ',len(protein_report))
print('Number of peptide entries: ',len(peptide_report))
print('Number of entries post-protein/peptide merge: ',len(merged_pep_prot_report))

merged_pep_prot_report['protein p-value'] = ttest_ind(merged_pep_prot_report[[exp_1_name_protein, exp_2_name_protein,exp_3_name_protein]], 
                                merged_pep_prot_report[[ctrl_1_name_protein, ctrl_2_name_protein,ctrl_3_name_protein]], axis=1, equal_var=True, alternative='two-sided',nan_policy='omit')[1]


merged_pep_prot_report['unaltered fold change (SetA)'] = merged_pep_prot_report['protein fold change'] #in this column the FC will not be changed regardless of p-value
merged_pep_prot_report['p-value filtered protein fold change (SetB)'] = merged_pep_prot_report['protein fold change'] #in this column the FC will be equal to 1 if p-value is > 0.05
merged_pep_prot_report['p-value OR log2(FC) filtered protein fold change (SetC)'] = merged_pep_prot_report['protein fold change'] #in this column the FC will be equal to 1 if p-value is > 0.05 OR the FC is insignificant
merged_pep_prot_report['p-value AND log2(FC) filtered protein fold change (SetD)'] = merged_pep_prot_report['protein fold change'] #in this column the FC will be equal to 1 if p-value is > 0.05 AND the FC is insignificant

merged_pep_prot_report.loc[merged_pep_prot_report['protein p-value'] > p_cutoff, ['p-value filtered protein fold change (SetB)']] = 1
merged_pep_prot_report.loc[merged_pep_prot_report['protein p-value'].isnull(), ['p-value filtered protein fold change (SetB)']] = 1

merged_pep_prot_report.loc[merged_pep_prot_report['protein p-value'] > p_cutoff, ['p-value OR log2(FC) filtered protein fold change (SetC)']] = 1
merged_pep_prot_report.loc[merged_pep_prot_report['protein p-value'].isnull(), ['p-value OR log2(FC) filtered protein fold change (SetC)']] = 1
merged_pep_prot_report.loc[abs(merged_pep_prot_report['Log2(protein fold change)']) < min_log_FC, ['p-value OR log2(FC) filtered protein fold change (SetC)']] = 1


merged_pep_prot_report.loc[abs((merged_pep_prot_report['Log2(protein fold change)'])<min_log_FC) & (merged_pep_prot_report['protein p-value']> p_cutoff),
                    ['p-value AND log2(FC) filtered protein fold change (SetD)']] = 1
merged_pep_prot_report.loc[abs((merged_pep_prot_report['Log2(protein fold change)'])<min_log_FC) & (merged_pep_prot_report['protein p-value'].isnull()),
                    ['p-value AND log2(FC) filtered protein fold change (SetD)']] = 1


merged_pep_prot_report = merged_pep_prot_report.replace(0, np.NaN) #replace empty values with NaN so mean can be taken without impact from 0s
merged_pep_prot_report['Ctrl mean peptide'] = merged_pep_prot_report[[ctrl_1_name_peptide, ctrl_2_name_peptide,ctrl_3_name_peptide]].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s
merged_pep_prot_report = merged_pep_prot_report.replace(0, np.NaN) #replace empty values with NaN so mean can be taken without impact from 0s
merged_pep_prot_report['SetA (LFQ intensity Exp1 normalized-peptide)'] = (merged_pep_prot_report[exp_1_name_peptide] / 
                                                                                   merged_pep_prot_report['unaltered fold change (SetA)'])
merged_pep_prot_report['SetA (LFQ intensity Exp2 normalized-peptide)'] = (merged_pep_prot_report[exp_2_name_peptide] / 
                                                                                   merged_pep_prot_report['unaltered fold change (SetA)'])
merged_pep_prot_report['SetA (LFQ intensity Exp3 normalized-peptide)'] = (merged_pep_prot_report[exp_3_name_peptide] / 
                                                                                   merged_pep_prot_report['unaltered fold change (SetA)'])


merged_pep_prot_report['SetB (LFQ intensity Exp1 normalized-peptide)'] = (merged_pep_prot_report[exp_1_name_peptide] / 
                                                                                   merged_pep_prot_report['p-value filtered protein fold change (SetB)'])
merged_pep_prot_report['SetB (LFQ intensity Exp2 normalized-peptide)'] = (merged_pep_prot_report[exp_2_name_peptide] / 
                                                                                   merged_pep_prot_report['p-value filtered protein fold change (SetB)'])
merged_pep_prot_report['SetB (LFQ intensity Exp3 normalized-peptide)'] = (merged_pep_prot_report[exp_3_name_peptide] / 
                                                                                   merged_pep_prot_report['p-value filtered protein fold change (SetB)'])

merged_pep_prot_report['SetC (LFQ intensity Exp1 normalized-peptide)'] = (merged_pep_prot_report[exp_1_name_peptide] / 
                                                                                   merged_pep_prot_report['p-value OR log2(FC) filtered protein fold change (SetC)'])
merged_pep_prot_report['SetC (LFQ intensity Exp2 normalized-peptide)'] = (merged_pep_prot_report[exp_2_name_peptide] / 
                                                                                   merged_pep_prot_report['p-value OR log2(FC) filtered protein fold change (SetC)'])
merged_pep_prot_report['SetC (LFQ intensity Exp3 normalized-peptide)'] = (merged_pep_prot_report[exp_3_name_peptide] / 
                                                                                   merged_pep_prot_report['p-value OR log2(FC) filtered protein fold change (SetC)'])

merged_pep_prot_report['SetD (LFQ intensity Exp1 normalized-peptide)'] = (merged_pep_prot_report[exp_1_name_peptide] / 
                                                                                   merged_pep_prot_report['p-value AND log2(FC) filtered protein fold change (SetD)'])
merged_pep_prot_report['SetD (LFQ intensity Exp2 normalized-peptide)'] = (merged_pep_prot_report[exp_2_name_peptide] / 
                                                                                   merged_pep_prot_report['p-value AND log2(FC) filtered protein fold change (SetD)'])
merged_pep_prot_report['SetD (LFQ intensity Exp3 normalized-peptide)'] = (merged_pep_prot_report[exp_3_name_peptide] / 
                                                                                   merged_pep_prot_report['p-value AND log2(FC) filtered protein fold change (SetD)'])
merged_pep_prot_report = merged_pep_prot_report.replace(0, np.NaN) #replace empty values with NaN so mean can be taken without impact from 0s
merged_pep_prot_report['SetA exp mean intensity peptide'] = merged_pep_prot_report[['SetA (LFQ intensity Exp1 normalized-peptide)', 
                                                                                                            'SetA (LFQ intensity Exp2 normalized-peptide)',
                                                                                                            'SetA (LFQ intensity Exp3 normalized-peptide)']].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s

merged_pep_prot_report['SetB exp mean intensity peptide'] = merged_pep_prot_report[['SetB (LFQ intensity Exp1 normalized-peptide)', 
                                                                                                            'SetB (LFQ intensity Exp2 normalized-peptide)',
                                                                                                            'SetB (LFQ intensity Exp3 normalized-peptide)']].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s                                                                                                                                                              
    
merged_pep_prot_report['SetC exp mean intensity peptide'] = merged_pep_prot_report[['SetC (LFQ intensity Exp1 normalized-peptide)', 
                                                                                                            'SetC (LFQ intensity Exp2 normalized-peptide)',
                                                                                                            'SetC (LFQ intensity Exp3 normalized-peptide)']].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s    
  
merged_pep_prot_report['SetD exp mean intensity peptide'] = merged_pep_prot_report[['SetD (LFQ intensity Exp1 normalized-peptide)', 
                                                                                                            'SetD (LFQ intensity Exp2 normalized-peptide)',
                                                                                                            'SetD (LFQ intensity Exp3 normalized-peptide)']].mean(axis=1,
                                                                                                                                        skipna=True) #the mean is taken without considering 0s                                                                                                                                                                 
    
merged_pep_prot_report['SetA exp mean fold change peptide normalized'] = merged_pep_prot_report['SetA exp mean intensity peptide']/merged_pep_prot_report['Ctrl mean peptide']                                                                                                                                                             
merged_pep_prot_report['SetB exp mean fold change peptide normalized'] = merged_pep_prot_report['SetB exp mean intensity peptide']/merged_pep_prot_report['Ctrl mean peptide'] 
merged_pep_prot_report['SetC exp mean fold change peptide normalized'] = merged_pep_prot_report['SetC exp mean intensity peptide']/merged_pep_prot_report['Ctrl mean peptide'] 
merged_pep_prot_report['SetD exp mean fold change peptide normalized'] = merged_pep_prot_report['SetD exp mean intensity peptide']/merged_pep_prot_report['Ctrl mean peptide'] 

                                                                                                                                                              
file_out_path = output_path + '\\Updated_report_20221206.csv'
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        merged_pep_prot_report.to_csv(filec,index=False)