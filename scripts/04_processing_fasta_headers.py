#!/usr/bin/env python3
# coding: utf-8

"""
This script modifies the header of fasta files by adding a lable
that indicates the source of the assembly (H = human, EN = natural environment, 
EB = built environment, O = other host).
"""

import os
import subprocess
import pandas as pd
import argparse

def map_accession_to_label(dataset_file):
    """
    Read the tsv files containing the metadata and make a dictionary that matches
    the accession number of a genome to its isolation source.
    """
    dataset = pd.read_csv(dataset_file, sep='\t', keep_default_na = False, na_values=[' '], index_col=0)

    accession_label_dict = {}

    for i in range(len(dataset)):
        current_assembly = dataset.loc[i,'assembly_accession']
        if current_assembly not in accession_label_dict.keys():
            accession_label_dict[current_assembly] = str(dataset.loc[i, 'label'])
    
    return accession_label_dict

def edit_fasta_headers(fasta_dir, output_dir, accession_label_dict):
    """
    Add the H, EN, EB, O labels to the fasta header and in the file name.
    """
    # Check if the output directory exists, if it doesn't exist it is created
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate over assembly directories in the fasta directory
    for root, directories, other_files in os.walk(fasta_dir):

        for assembly_dir in directories:
            assembly_path = os.path.join(fasta_dir, assembly_dir)

            # Check if it's a directory and if it's in the accession_label_dict
            if os.path.isdir(assembly_path) and assembly_dir in accession_label_dict:
                label = accession_label_dict[assembly_dir]
        
            # Iterate over the files in the assembly direcctory
            for file in os.listdir(assembly_path):
                # Check file extension
                if file.endswith('.fna'):
                    # Path to the original fasta file
                    original_file_path = os.path.join(assembly_path, file)
                    # Path to the new fasta file with label in the header
                    new_fastafile = os.path.join(output_dir, f"{label}_{assembly_dir}_new_header.fna")

                    # Construct and execute the sed command     
                    command = (f"sed 's/>/>{label}_/' {original_file_path} > {new_fastafile}")
                    subprocess.run(command, shell=True)

def main():
    # initiate argument parsing
    parser = argparse.ArgumentParser()
    
    # add flags for command line argumnets
    parser.add_argument('--metadata', help='Path to the csv file with the labelled genomes', default='./results/labelled_genomes.tsv')
    parser.add_argument('--fasta_dir', help = 'path to the directory with the original fasta files', default='./data/raw_data/assemblies_fasta/ncbi_dataset/data/')
    parser.add_argument('--out_dir', help = 'Path to the directory where the fasta files with the new header will be stored', default='./data/new_fasta_headers/')
    
    # add flags input to varaible
    args = parser.parse_args()

    # access the variables from the command line
    dataset_file = args.metadata
    fasta_dir = args.fasta_dir
    output_dir = args.out_dir
    
    # Call the function to load metadata 
    accession_label_dict = map_accession_to_label(dataset_file)

    # Call the function to edit the fasta headers
    edit_fasta_headers(fasta_dir, output_dir, accession_label_dict)

if __name__ == "__main__":
    main()
