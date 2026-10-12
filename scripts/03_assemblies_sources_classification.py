#!/usr/bin/env python
# coding: utf-8

"""
This script classifies genome assemblies based on the isolation
source and host metadata. 
"""

import pandas as pd
import argparse
import sys
import os
from tqdm import tqdm

def prepare_input_files(dataset_file, hosts_file, iso_sources_file):
    """
    Import files with the assemblies metadata, the hosts terms classification and the
    isolation sources classification. 
    H = human, EN = natural environment, EB = built environment, O = other host, U = unknown
    """
    # genome asse,blies metadata
    dataset = pd.read_csv(dataset_file, sep='\t', keep_default_na = False, na_values=[' '], index_col=0)
    dataset = dataset.fillna('NA')
    dataset['host_label'] = len(dataset)*['NA']
    dataset['isolation_source_label'] = len(dataset)*['NA']
    dataset['label'] = len(dataset)*['NA']

    # hosts and isoilation sources labels
    hosts_df = pd.read_csv(hosts_file, sep='\t',keep_default_na = False, na_values=[' '],header=None, names=['Host','classification'])
    hosts_dict = dict(zip(hosts_df['Host'], hosts_df['classification']))
    iso_sources_df = pd.read_csv(iso_sources_file, sep='\t',keep_default_na = False, na_values=[' '], header=None, names=['Isolation_source','classification'])
    isolation_sources_dict = dict(zip(iso_sources_df['Isolation_source'], iso_sources_df['classification']))
    print('Files imported successfully')
    return dataset, isolation_sources_dict, hosts_dict


def add_labels_to_dataset(dataset, isolation_sources_dict, hosts_dict):
    """
    Fill the 'host label' and 'isolation source lable' columns 
    """
    # add isolation sources labels
    for i in tqdm(range(len(dataset['isolation_source'])), desc='Adding labels to isolation sources'):
        current_iso_source = dataset.loc[i,'isolation_source']
        if current_iso_source in list(isolation_sources_dict.keys()):
            dataset.loc[i,'isolation_source_label'] = isolation_sources_dict[current_iso_source]
    print('Isolation source labels added successfully')

    # add hosts labels
    for i in tqdm(range(len(dataset['host'])), desc='Adding labels to hosts'):
        current_host = dataset.loc[i,'host']
        if current_host in list(hosts_dict.keys()):
            dataset.loc[i,'host_label'] = hosts_dict[current_host]
    print('Host labels added successfully')
    return dataset


def checkpoint(dataset):
    """
    Checking if all the assemblies have a lable both for the host and the isolation source
    """
    host_problem = False
    iso_source_problem = False
    
    # checkpoint for isolation source column
    if (dataset['isolation_source_label']== 'NA').sum() != 0:
        print('There are unclassified isolation sources')
        iso_source_problem = True  
        # save unclassified iso sources to txt
        output_iso_sources = "./classification_files/isolation_sources_not_classified.txt"
        unclassified_iso_sources_uniq = list(set(dataset.loc[dataset['isolation_source_label'] == 'NA', 'isolation_source'].tolist()))
        with open(output_iso_sources, 'w') as file_2:
            for string in unclassified_iso_sources_uniq:
                file_2.write(f"{string}\n")
    
    # checkpoint for host column
    if (dataset['host_label']== 'NA').sum() != 0:
        print('There are unclassified hosts')
        host_problem = True
        # save unclassified hosts to txt
        output_host = "./classification_files/hosts_not_classified.txt"
        unclassified_hosts_uniq = list(set(dataset.loc[dataset['host_label'] == 'NA', 'host'].tolist()))
        with open(output_host, 'w') as file_1:
            for string in unclassified_hosts_uniq:
                file_1.write(f"{string}\n")
    
    # Exit if host or iso source problem or both are True     
    if iso_source_problem or host_problem or (iso_source_problem and host_problem):
        error_message =  'Check the files in ./classification_files/ directory, update the files with classification and re-run this script.'
        print(error_message)
        sys.exit(1)

def match_labels(dataset, out_dir):
    """
    compare the label assigned to the host with the one assigned to the isolation source
    to decide on a unique label for each assembly
    """
    for index, row in tqdm(dataset.iterrows(), desc='Matching labels', total=dataset.shape[0]):
        if (row['host_label'] in ['H', 'U']) and (row['isolation_source_label'] in ['H', 'U']):
            dataset.at[index, 'label'] = 'H'

        elif (row['host_label'] in ['EN', 'U']) and (row['isolation_source_label'] in ['EN', 'U']):
            dataset.at[index, 'label'] = 'EN'
        
        elif (row['host_label'] in ['EB', 'U']) and (row['isolation_source_label'] in ['EB', 'U']):
            dataset.at[index, 'label'] = 'EB'
        
        elif (row['host_label'] in ['O', 'U']) and (row['isolation_source_label'] in ['O', 'U']):
            dataset.at[index, 'label'] = 'O'
        
        else:
            dataset.at[index, 'label'] = 'U'
        
        if (row['host_label']=='U') and (row['isolation_source_label']=='U'):
            dataset.at[index, 'label'] = 'U'

    # Check if the output directory exists, if it doesn't exist it is created
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    dataset.to_csv(f'{out_dir}labelled_genomes.tsv', sep="\t") 
    dataset[['assembly_accession', 'label']].to_csv(f"{out_dir}acc_number_and_label_only.tsv", sep= "\t")

def main():
    # initiate argument parsing
    parser = argparse.ArgumentParser()

    # Add flags for command line arguments
    parser.add_argument('--metadata', type=str, help='Path to the tsv file with assemblies metadata', default='./results/assemblies_metadata.tsv')
    parser.add_argument('--hosts', type=str, help='Path to txt file with hosts classification', default='./classification_files/hosts_classified.tsv')
    parser.add_argument('--iso_sources', type=str, help='Path to txt file with isolation sources classification', default='./classification_files/isolation_sources_classified.tsv')
    parser.add_argument('--out_dir', type=str, help='Path to the output directory', default='./results/')

    # Add flags input to a variable
    args = parser.parse_args()

    # Access the values of the arguments
    dataset_file = args.metadata
    hosts_file = args.hosts
    iso_sources_file = args.iso_sources
    out_dir = args.out_dir

    # Call the function to prepare the dataset
    dataset, isolation_sources_dict, hosts_dict = prepare_input_files(dataset_file, hosts_file, iso_sources_file)

    # Call the function to add lables
    dataset = add_labels_to_dataset(dataset, isolation_sources_dict, hosts_dict)

    # Checkpoint function
    checkpoint(dataset)

    # Call the function to compare the lables of host and isolation source
    match_labels(dataset, out_dir)


if __name__ == "__main__":
    main()
