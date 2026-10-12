#!/usr/bin/env python
# coding: utf-8

"""
This scripts preprocessed the header in the input multifasta and then 
runs CD-HIT to cluster protein sequences.
"""

import subprocess
import argparse
#import re
"""
def preprocess_headers(input_file, output_file):
    
   # Removes special characters from the headers in the multifatsa with protein sequences
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Replace special characters in the heaer with spaces using regular expression
                header = re.sub(r'[^a-zA-Z0-9_]', ' ', line)
                outfile.write(header)
            else:
                # Write sequence lines unchanged
                outfile.write(line)
"""

def run_CD_HIT(input_sequences, output_dir, identity_threshold, threads):
    """
    Cluster the protein sequences in the multifasta file
    """
    command_list = ['cd-hit', '-i', f"{input_sequences}", '-o', f"{output_dir}CD_HIT_output", '-c', f"{identity_threshold}", '-T', f"{threads}", '-d 0' , '-M 0']
    command = ' '.join(command_list)

    subprocess.run(command, shell=True)

def main():

    # command lie arguments parser
    parser = argparse.ArgumentParser()

    # Add flags for cl arguments
    parser.add_argument('--input_sequences', help='Protein sequences multifasta, Prodigal output', default='./results/predicted_proteins.faa')
    #parser.add_argument('--processed_seq', help = 'Path to the multifasta file with predicted protein sequences with special characters removed from the headers', default='./results/predicted_proteins_processed_header.faa')
    parser.add_argument('--out_dir', help = "Path to output directory", default='./results/')
    parser.add_argument('--id', help = "identity level for the clustering [0.1 - 1.0]", default=0.9)
    parser.add_argument('--threads', help = "inumber of threads used by cd-hit", default=10)

    # Save arguments to a variable
    args = parser.parse_args()

    # Access the arguments
    proteins_multifasta = args.input_sequences
    #proteins_multifasta_processed = args.processed_seq
    output_dir = args.out_dir
    identity_threshold = args.id
    threads = args.threads

    # Call the function the prepares the input file
    #preprocess_headers(proteins_multifasta,proteins_multifasta_processed )
    # Call the function that runs CD-HIT
    run_CD_HIT(proteins_multifasta, output_dir, identity_threshold, threads)

if __name__ == "__main__":
    main()