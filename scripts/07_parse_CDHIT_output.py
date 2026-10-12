#!/usr/bin/env python
# coding: utf-8

"""
This scripts parses the output of the CD-HIT clustering algorithm and creates 
a csv file with the clusters composition.
"""

import pandas as pd
import argparse
from tqdm import tqdm

def parse_clusters_file(clusters_file):
    """
    Create a dictionary with cluster name as key and 
    ids of sequences in that cluster as list of values
    """
    clusters_dict = {}
    current_cluster = None

    with open(clusters_file, 'r') as file:
        for line in tqdm(file, desc = "Processing .clstr file. Step 1/2"):
            line = line.strip()
            # New cluster
            if line.startswith('>'):
                line = line.replace('>','').replace(' ', '_')
                current_cluster = line
                clusters_dict[current_cluster] = []
            # Sequences in the current_cluster
            elif current_cluster is not None and line: 
                line = line.split(' ')
                seq = (line[1].replace('>','') + '_' + line[-1]).replace('...','')
                clusters_dict[current_cluster].append(seq)
    return clusters_dict

def create_clusters_csv(clusters_dict, output_dir):
    clusters_df = pd.DataFrame(columns=['cluster', 'repr_seq', 'tot_H', 'tot_EN','tot_EB', 'tot_O'])
    for key, value in tqdm(clusters_dict.items(), desc= "Processing .clstr file. Step 2/2"):
        new_row = []
        H_count = 0
        EN_count = 0
        EB_count = 0
        O_count = 0
        new_row.append(key)
        for sequence in value:
            if sequence.startswith('H'):
                H_count += 1
            elif sequence.startswith('EN'):
                EN_count += 1
            elif sequence.startswith('EB'):
                EB_count += 1
            elif sequence.startswith('O'):
                O_count += 1
        for sequence in value:
            if sequence.endswith('*'):
                new_row.append(sequence)
        new_row.append(H_count)
        new_row.append(EN_count)
        new_row.append(EB_count)
        new_row.append(O_count)
        clusters_df.loc[len(clusters_df)] = new_row   
    clusters_df.to_csv(f"{output_dir}clusters_composition.csv")

def main():

    # command lie arguments parser
    parser = argparse.ArgumentParser()

    # Add flags for cl arguments
    parser.add_argument('--input_file', help='CD-HIT clusters output file', default='./results/CD_HIT_output.clstr')
    parser.add_argument('--out_dir', help = "Path to output directory", default='./results/')

    # Save arguments to a variable
    args = parser.parse_args()

    # Access the arguments
    clusters_file = args.input_file
    output_dir = args.out_dir

    # Call the functions
    clusters_dictionary = parse_clusters_file(clusters_file)
    create_clusters_csv(clusters_dictionary, output_dir)

if __name__ == "__main__":
    main()