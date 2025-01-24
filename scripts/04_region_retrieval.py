# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 14:14:51 2025

@author: lafields2
"""


import csv
import pandas as pd
import re
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
from bioservices import UniProt
import requests, sys

db_path = r"C:\Users\lawashburn\Desktop\ALC50_Mass_Search_Files\mouse_uniprotkb_taxonomy_id_10090_AND_reviewe_2024_12_14.fasta"
data_path = r"D:\Manuscripts\2024_Haiyan_VitaminA\glyco_Uniprot_test\Conformotypic peptides match_Lung_glyco.csv"
output_path = r"D:\Manuscripts\2024_Haiyan_VitaminA\region_test"

def get_nth_residue(prot_query, n):
    """
    Extracts the nth residue of a protein sequence from UniProt.

    Parameters:
    prot_query : str
        The UniProt accession code or protein name for the query.
    n : int
        The residue number (1-based index) to retrieve from the protein sequence.

    Returns:
    str
        The nth residue (amino acid) in the protein sequence. If `n` is out of bounds, return an error message.
    """
    
    # URL to access the UniProt data in FASTA format
    url = f"https://www.uniprot.org/uniprot/{prot_query}.fasta"
    
    # Send a GET request to fetch the protein sequence in FASTA format
    response = requests.get(url)
    
    if not response.ok:
        raise ValueError(f"Error retrieving data for {prot_query}. HTTP status code: {response.status_code}")
    
    # Extract the sequence (the lines after the first line, which contains the header)
    fasta_data = response.text.strip().split("\n")[1:]
    sequence = "".join(fasta_data).replace("\n", "")  # Concatenate all lines to get the full sequence
    
    # Check if n is within the bounds of the sequence
    if n <= 0 or n > len(sequence):
        return f"Error: The residue number {n} is out of bounds for this sequence (length: {len(sequence)})."
    
    # Return the nth residue (1-based indexing)
    return sequence[n-1]

def UniProt_search(prot_query):
    """
   Search UniProt for protein modifications and return a DataFrame containing the modified residue information.
   
   Parameters:
   prot_query : str
       The query string for searching UniProt, representing a protein or a list of proteins.
   
   Returns:
   pandas.DataFrame
       A DataFrame with the columns:
           - 'Res #': The residue number(s) of the modified residue(s).
           - 'type': The type of modification (e.g., phosphorylation).
           - 'evidence': The evidence for the modification, set to 'UniProt'.
       
   If no modified residues are found, an empty DataFrame is returned with the columns: ['Res #', 'type', 'evidence'].
   
   The function performs the following steps:
   1. Initializes a UniProt object and retrieves data based on the protein query.
   2. Extracts the 'Modified residue' column from the returned data.
   3. If no modifications are found (NaN value), an empty DataFrame is returned.
   4. Parses the 'Modified residue' string to isolate modified residue details.
   5. Splits the details into three columns: 'Res #', 'type', and 'evidence'.
   6. Cleans up the 'type' column, removing unnecessary text such as 'note=' and quotes.
   7. Extracts isoform and residue numbers from the 'Res #' column and adjusts the formatting.
   8. Converts residue numbers to integers, where possible.
   9. If isoform information is available, it is incorporated into the 'type' description.
   10. The isoform info is removed from the DataFrame after use.
   
   Example usage:
   >>> UniProt_search('P12345')
   """
    u = UniProt()
    res = u.get_df(prot_query.split())
    
    # Extract the 'Region' value
    ptm = res['Region'].iloc[0]
    
    ptm_df = pd.DataFrame()


    # Attempt to extract the 'region' column (make sure the column name matches exactly)
    if 'Region' not in res.columns:
        # No region column available
        return pd.DataFrame(columns=['region_start', 'region_end', 'region_note', 'region_evidence'])
    
    region_entry = res['Region'].iloc[0]  # In your code, you might loop or handle multiple rows
    
    # If 'region' is NaN or missing, return an empty DataFrame
    if isinstance(region_entry, float) and np.isnan(region_entry):
        return pd.DataFrame(columns=['region_start', 'region_end', 'region_note', 'region_evidence'])
    
    # Split out each separate region by 'REGION '
    # region_entry might look like:
    # "REGION 307..531; /note="Interaction with POU5F1"; /evidence="ECO:0000250|UniProtKB:P14618";
    #  REGION 389..433; /note="Intersubunit contact""
    region_split = region_entry.split('REGION ')[1:]
    if not region_split:
        return pd.DataFrame(columns=['region_start', 'region_end', 'region_note', 'region_evidence'])
    
    # Convert to a DataFrame so we can parse each piece
    region_df = pd.DataFrame({'raw': region_split})
    
    # Extract the start and end positions. Example format: 307..531
    region_df['region_start'] = region_df['raw'].str.extract(r'(\d+)\.\.')
    region_df['region_end']   = region_df['raw'].str.extract(r'\.\.(\d+)')
    
    # Convert start/end to numeric if desired
    region_df['region_start'] = pd.to_numeric(region_df['region_start'], errors='coerce')
    region_df['region_end']   = pd.to_numeric(region_df['region_end'],   errors='coerce')
    
    # Extract note (e.g., /note="Interaction with POU5F1")
    region_df['region_note'] = region_df['raw'].str.extract(r'/note="([^"]+)"')
    
    # Extract evidence (e.g., /evidence="ECO:0000250|UniProtKB:P14618")
    region_df['region_evidence'] = region_df['raw'].str.extract(r'/evidence="([^"]+)"')
    
    # Clean up the DataFrame by keeping only the columns we care about
    region_df = region_df[['region_start', 'region_end', 'region_note', 'region_evidence']]
    
    # You can fill missing evidence with something descriptive or leave it NaN
    region_df['region_evidence'].fillna(value='No_evidence_tag', inplace=True)
    
    return region_df


def pride_search(prot_query):
    """
    Search the PRIDE database for proteomics PTM (post-translational modification) data and return a DataFrame containing the modified residue information.
    
    Parameters:
    prot_query : str
        The query string representing the protein to search in the PRIDE database.
    
    Returns:
    pandas.DataFrame
        A DataFrame with the columns:
            - 'Res #': The absolute position of the modified residue in the peptide.
            - 'type': The type of modification (e.g., phosphorylation).
            - 'evidence': The evidence source, set to 'PRIDE'.
        
    If no data is found for the protein or PTMs, an empty DataFrame with the columns: ['Res #', 'type', 'evidence'] is returned.
    
    The function performs the following steps:
    1. Constructs a request URL to search the PRIDE database for proteomics PTM data based on the protein query.
    2. Sends a GET request to the PRIDE API.
    3. Checks for non-OK responses (e.g., 404 status) and returns an empty DataFrame if no data is found.
    4. Processes the response text to extract PTM data (e.g., modification name, peptide position, modification position).
    5. Stores the extracted data in separate lists for further processing.
    6. If no PTM data is found, returns an empty DataFrame.
    7. Converts the peptide and modification positions to numeric values, and calculates the absolute position of the PTM.
    8. Returns a formatted DataFrame with columns for residue position, modification type, and evidence.
    
    Example usage:
    >>> pride_search('P12345')
    """
    requestURL = "https://www.ebi.ac.uk/proteins/api/proteomics-ptm/" + prot_query
    r = requests.get(requestURL, headers={"Accept": "application/json"})

    # Handle non-OK responses (such as 404)
    if not r.ok:
        print(f"No proteomics-PTM data found for {prot_query}. HTTP status code: {r.status_code}")
        # Return an empty DataFrame with the same structure
        return pd.DataFrame(columns=['Res #','type','evidence'])
        
    responseBody = r.text

    # If there is no data or the response is empty, return an empty DataFrame
    if not responseBody.strip():
        print(f"No data returned for {prot_query}.")
        return pd.DataFrame(columns=['Res #','type','evidence'])

    responseBody2 = list(responseBody.split("type"))

    begin_storage = []
    position_storage = []
    ptm_name_storage = []

    for a in range(1, len(responseBody2)):
        responseBody3 = responseBody2[a]
        responseBody4 = list(responseBody3.split(','))
        
        for b in responseBody4:
            if 'begin' in b:
                b_new = b.replace('"begin":"','').replace('"','')
                begin_storage.append(b_new)
            elif 'position' in b:
                b_new = b.replace('"position":','').replace('"','')
                position_storage.append(b_new)
            elif 'ptms' in b:
                b_new = b.replace('"ptms":[{"name":"','').replace('"','')
                ptm_name_storage.append(b_new)

    # If no PTM data was found
    if not ptm_name_storage:
        print(f"No PTM data found for {prot_query}.")
        return pd.DataFrame(columns=['Res #','type','evidence'])

    pride_mods_df = pd.DataFrame({'Mod name':ptm_name_storage, 
                                  'Peptide begin position':begin_storage, 
                                  'Mod position':position_storage})

    # Convert to numeric if possible
    pride_mods_df[['Peptide begin position', 'Mod position']] = pride_mods_df[['Peptide begin position', 'Mod position']].apply(pd.to_numeric, errors='coerce')

    # Drop rows with missing numeric conversion if needed
    pride_mods_df.dropna(subset=['Peptide begin position', 'Mod position'], inplace=True)

    # Calculate the absolute PTM location
    pride_mods_df['Absolute PTM location'] = (pride_mods_df['Peptide begin position'] + pride_mods_df['Mod position']) - 1
    
    res_no = pride_mods_df['Absolute PTM location'].values.tolist()
    type_list =  pride_mods_df['Mod name'].values.tolist()

    pride_mods_df_format = pd.DataFrame({'Res #': res_no, 'type': type_list})
    pride_mods_df_format['evidence'] = 'PRIDE'

    return pride_mods_df_format

def process_dataframe(df, prot_column, n_column, result_column="LiP site AA"):
    """
    Processes an entire DataFrame to retrieve nth residues for each row.

    Parameters:
    df : pd.DataFrame
        Input DataFrame containing protein queries and residue numbers.
    prot_column : str
        Column name containing UniProt accession codes or protein names.
    n_column : str
        Column name containing residue numbers.
    result_column : str
        Column name to store the results. Default is 'nth_residue'.

    Returns:
    pd.DataFrame
        Updated DataFrame with an additional column containing the nth residues or error messages.
    """
    def safe_get_residue(row):
        try:
            return get_nth_residue(row[prot_column], row[n_column])
        except Exception as e:
            return f"Error: {str(e)}"
    
    # Apply the function row-wise and store results in the new column
    df[result_column] = df.apply(safe_get_residue, axis=1)
    return df

fasta_to_df = []
title_to_df = []
accession_to_df = []
with open(db_path) as fasta_file:  # Will close handle cleanly
    for title, sequence in SimpleFastaParser(fasta_file):
        fasta_to_df.append(sequence)
        title_to_df.append(title)


fasta_df = pd.DataFrame()
fasta_df['Record ID'] = title_to_df
fasta_df['Protein Sequence'] = fasta_to_df
fasta_df['Record ID']= fasta_df['Record ID'].str.split('|').str[1].str.strip('|')

data = pd.read_csv(data_path)

merged = data.merge(fasta_df, left_on=('Leading razor protein_P'),right_on='Record ID', how='inner')
queries = len(merged)
protein_seqs = merged['Protein Sequence'].values.tolist()
peptide_seqs = merged['Sequence S'].values.tolist()

pep_storage = []
prot_storage = []
start_index = []
end_index = []
peptide_length = []
protein_cleavage_before = []
peptide_first_cleavage = []
peptide_last_cleavage = []
protein_cleavage_after = []
first_cleavage_type = []
second_cleavage_type = []

for a in range(0,queries):
    protein_query = protein_seqs[a]
    peptide_query = peptide_seqs[a]
    
    if peptide_query in protein_query:
        prot_index = protein_query.index(peptide_query)
        start_index.append(prot_index+1)
        pep_len = len(peptide_query)
        peptide_length.append(pep_len)
        end_index.append(prot_index+pep_len)
        
        if prot_index == 0:
            prot_res_before = 'None'
            protein_cleavage_before.append(prot_res_before)
            first_cleavage_type.append('Non-tryptic')
        else:
            prot_res_before = protein_query[prot_index-1]
            protein_cleavage_before.append(prot_res_before)
            
            if prot_res_before == 'K':
                first_cleavage_type.append('Tryptic')
            elif prot_res_before == 'R':
                first_cleavage_type.append('Tryptic')
            else:
                first_cleavage_type.append('Non-tryptic')
            
        first_pep_res = peptide_query[0]
        last_pep_res = peptide_query[pep_len-1]
        
        if last_pep_res == 'K':
            second_cleavage_type.append('Tryptic')
        elif last_pep_res == 'R':
            second_cleavage_type.append('Tryptic')
        else:
            second_cleavage_type.append('Non-tryptic')
        
        if (prot_index+pep_len)==len(protein_query):
            prot_res_after = 'None'
            protein_cleavage_after.append(prot_res_after)
        else:
            prot_res_after = protein_query[prot_index+pep_len]
            protein_cleavage_after.append(prot_res_after)

        peptide_first_cleavage.append(first_pep_res)
        peptide_last_cleavage.append(last_pep_res)
        pep_storage.append(peptide_query)
        prot_storage.append(protein_query)

index_report = pd.DataFrame()
index_report['Peptide sequence'] = pep_storage
index_report['Protein sequence'] = prot_storage
index_report['Peptide start index'] = start_index
index_report['Peptide end index'] = end_index
index_report['Peptide length'] = peptide_length
index_report['Protein cleavage before'] = protein_cleavage_before
index_report['Peptide first residue'] = peptide_first_cleavage
index_report['Peptide last residue'] = peptide_last_cleavage
index_report['Protein cleavage after'] = protein_cleavage_after
index_report['First cleavage site'] = first_cleavage_type
index_report['Second cleavage site'] = second_cleavage_type

conditions = [
    index_report['First cleavage site'].eq('Tryptic') & index_report['Second cleavage site'].eq('Tryptic'),
    index_report['First cleavage site'].eq('Tryptic') & index_report['Protein cleavage after'].eq('None'),
    index_report['Second cleavage site'].eq('Tryptic') & index_report['Protein cleavage before'].eq('None')
]

choices = ['Fully tryptic','Fully tryptic','Fully tryptic']

index_report['Peptide status'] = np.select(conditions, choices, default='Half tryptic')

merged_final = merged.merge(index_report,left_on=['Sequence S','Protein Sequence'],right_on=['Peptide sequence','Protein sequence'])

half_tryp_filter = merged_final.loc[merged_final['Peptide status'] == 'Half tryptic']

begin_nontryp_filter = half_tryp_filter.loc[half_tryp_filter['First cleavage site'] == 'Non-tryptic']
begin_nontryp_filter['LiP site'] = begin_nontryp_filter['Peptide start index']

end_nontryp_filter = half_tryp_filter.loc[half_tryp_filter['Second cleavage site'] == 'Non-tryptic']
end_nontryp_filter['LiP site'] = end_nontryp_filter['Peptide end index'] +1

lip_site_table_draft = pd.concat([begin_nontryp_filter,end_nontryp_filter])

lip_site_table = process_dataframe(lip_site_table_draft, prot_column="Proteins_P", n_column="LiP site")

file_out_path = output_path + '\\lip_sites_unfiltered.csv'
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        lip_site_table.to_csv(filec,index=False)


concat_ptm_df_filtered = pd.DataFrame()

for index, row in lip_site_table.iterrows():
    prot_query = row['Leading razor protein_P']

    index_query = row['LiP site']
    
    uniprot_ptm_df = UniProt_search(prot_query)

#%%
    #pride_ptm_df = pride_search(prot_query)

    residue_lip = get_nth_residue(prot_query,index_query)
    
    #ptm_df = pd.concat([uniprot_ptm_df, pride_ptm_df],ignore_index=True)
    ptm_df = uniprot_ptm_df.copy()
    
    # print(f'Current index query:{index_query}')
    # print(f'Current region start window:{index_query - PTM_site_variance}')
    # print(f'Current region end window:{index_query + PTM_site_variance}')
    # print(f'Current length of PTM dataframe:{len(ptm_df)}')
    
    ptm_df_filtered = ptm_df[(ptm_df['region_start'] <= index_query) & (ptm_df['region_end'] >= index_query)]
    
    # print(f'Filtered length of PTM dataframe:{len(ptm_df_filtered)}')
    # print('\n\n\n\n\n')
    ptm_df_filtered['Peptide status'] = row['Peptide status']
    ptm_df_filtered['LiP site'] = row['LiP site']
    ptm_df_filtered['LiP Site Residue'] = residue_lip
    ptm_df_filtered['Second cleavage site'] = row['Second cleavage site']
    ptm_df_filtered['First cleavage site'] = row['First cleavage site']
    ptm_df_filtered['Record ID'] = row['Record ID']
    ptm_df_filtered['End position_P'] = row['End position_P']
    ptm_df_filtered['Start position_P'] = row['Start position_P']
    ptm_df_filtered['Sequence S'] = row['Sequence S']
    ptm_df_filtered['Sequence T'] = row['Sequence T']
    ptm_df_filtered['Sequence P'] = row['Sequence P']
       
    concat_ptm_df_filtered = pd.concat([concat_ptm_df_filtered, ptm_df_filtered],ignore_index=True)
    
    # applicable_mods = ptm_df_filtered['IDed sites'].values.tolist()
    
    # if len(applicable_mods)>0:
    #     applicable_mod_storage.append(applicable_mods)
    # else:
    #     applicable_mod_storage.append('No modifications found')

# lip_site_table['Modifications'] = applicable_mod_storage

# lip_site_table = lip_site_table.drop(columns=['Protein sequence','Peptide sequence','Peptide start index','Peptide end index','Peptide length',
#                               'Protein cleavage before','Peptide first residue','Peptide last residue','Protein cleavage after','Protein Sequence'])

concat_ptm_df_filtered = concat_ptm_df_filtered.rename(columns={'IDed sites':'Modification type',
                                                                'Res #':'Modified Res #',
                                                                'evidence':'Database',
                                                                'Record ID':'Protein ID'})
#concat_ptm_df_filtered.drop(columns='type')
concat_ptm_df_filtered_reordered = concat_ptm_df_filtered.loc[:, ['Sequence S', 'Sequence T','Sequence P', 'Start position_P','End position_P',
                                                                  'Protein ID','First cleavage site','Second cleavage site','Peptide status',
                                                                  'LiP site','LiP Site Residue', 'region_start','region_end','region_note','region_evidence']]

#%%

#concat_ptm_df_filtered_reordered['ProtID-LiPsite'] = concat_ptm_df_filtered_reordered['Protein ID'].astype(str) + '-' + concat_ptm_df_filtered_reordered['LiP site'].astype(str)
#concat_ptm_df_filtered_reordered['ProtID-ModSite (PRIDE/UniProt)'] = concat_ptm_df_filtered_reordered['Protein ID'].astype(str) + '-' + concat_ptm_df_filtered_reordered['Modified Res #'].astype(str)

#final_df= pd.merge(concat_ptm_df_filtered_reordered,data,on=['Sequence S','Sequence T','Sequence P','Start position_P','End position_P'], how='left')

file_out_path = output_path + '\\lip_site_region_finder_results_Output_full.csv'
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        concat_ptm_df_filtered_reordered.to_csv(filec,index=False)

        