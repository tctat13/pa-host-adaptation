#!/usr/bin/env python
# coding: utf-8

import random
import pandas as pd
import argparse
import sys
import os
import shutil

def infected_humnas_filter(dataset_file, filter):
    """
    Filter the human isolates so that the dataset only contains assemblies of 
    bacterial genomes isolated from (potentially) infected humans.
    """
    dataset = pd.read_csv(dataset_file, sep = '\t', keep_default_na = False, na_values=[' '])
    host_terms= []
    isolation_source_terms = []
    divider = '###'
    with open(filter) as f:
        # Iterate over each line in the file
        for line in f:
            term = line.strip()
            # If the divider is found, stop adding lines to the host list
            if divider in term:
                break
            host_terms.append(term)
        
        # Countinue reading lines and add them to the isolation source list
        for line in f:
            term = line.strip()
            isolation_source_terms.append(term)
    
    human_isolates_filtered = dataset[dataset['host'].isin(host_terms) & dataset['isolation_source'].isin(isolation_source_terms)]

    return human_isolates_filtered, dataset


def random_subsample(human_isolates_filtered, random_seed, dataset, out_dir):
    """
    Take a randon subsample of the human isolates so that size_group_H = size_group_EN+EB+O 
    """
    # Check if the output directory exists, if it doesn't exist it is created
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # If random_seed is not None, set random_state to random_seed, else set it to 5 
    random_state = random_seed if random_seed is not None else 5
    
    # Check if the number of H isolates is bigger that the number of EN isolates
    sum_EN_EB_O = (dataset['label'].isin(['EN', 'EB', 'O'])).sum()
    num_H = (dataset['label'] == 'H').sum()
    if num_H > sum_EN_EB_O:
        # Determine the subsample size = number of EN isolates
        subsample_size = sum_EN_EB_O
        # Take a random subsample of the H isolates
        human_isolates_subsample = human_isolates_filtered.sample(n=subsample_size, replace=False, random_state=random_state) 
        # The final dataset contains the H subsample and the assemblies with label EN, EB and O
        final_dataset = dataset.loc[~dataset['label'].isin(['H', 'U'])].append(human_isolates_subsample, ignore_index=True)
        final_dataset.to_csv(f'{out_dir}subsampled_genomes.tsv', sep="\t", index = False) 

        final_dataset_accession_numbers = final_dataset['assembly_accession'].tolist()

    else:
        error_message = '''The number of human isolates is smaller than the number of the natural environmnet isolates so this step 
                            is unnecessary. 
                            Read the instruction in the README.md file on how to run the analysis without performing this step.
                            Alternatively, it is possible to modify the script to take a customized subsample.'''
        print(error_message)
        sys.exit(1)
    
    return final_dataset_accession_numbers

def copy_fasta_files(accession_numbers, source_dir, target_dir):
    # Check if the output directory exists, if it doesn't exist it is created
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Loop through all the genomes in the source directory
    for filename in os.listdir(source_dir):
        temp_filename = '_'.join(filename.split('_')[1:3])
        source_path = os.path.join(source_dir, filename)

        # Check if file exists 
        if os.path.isfile(source_path) and temp_filename in accession_numbers:
            target_path = os.path.join(target_dir, filename)
            # command to copy file from source directory to target directory
            shutil.copy(source_path, target_path)
        
def main():
    # command lie arguments parser
    parser = argparse.ArgumentParser()

    # Add flags for command-line arguments
    parser.add_argument('--in_dir', help = 'Path to directory with fasta files with new header', default='./data/new_fasta_headers/')
    parser.add_argument('--dataset', help='Path to the tsv file with the seze-filtered genomes', default='./results/size_filtered_genomes.tsv') 
    parser.add_argument('--filtering_terms', type=str, help='Path to the file with the terms to filter the human isolates', default='./classification_files/filter_human_isolates.txt')
    parser.add_argument('--seed', type=int, help='Set seed for reproducibility', default=None)
    parser.add_argument('--out_dir', type=str, help='Path to the output directory', default='./results/')
    parser.add_argument('--subsample_dir', type=str, help='Path to the directory where the fasta files will be copied', default='./data/dataset_subsample/')
    
    # Save command line arguments to a variable
    args = parser.parse_args()

    # Access the values of the arguments
    assemblies_dir = args.in_dir
    dataset_file = args.dataset
    filter_txt = args.filtering_terms
    random_seed = args.seed
    out_dir = args.out_dir
    dataset_subsample_dir = args.subsample_dir

    # Call the function to filter the human isolates
    human_isolates_filtered, dataset = infected_humnas_filter(dataset_file, filter_txt)

    # Call the function to take a random subsample of the human isolates
    accession_numbers = random_subsample(human_isolates_filtered, random_seed, dataset, out_dir)

    # Copy the selected genomes to a different directory
    copy_fasta_files(accession_numbers, assemblies_dir, dataset_subsample_dir)

if __name__ == "__main__":
    main()