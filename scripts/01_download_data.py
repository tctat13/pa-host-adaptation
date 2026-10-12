#!/usr/bin/env python
# coding: utf-8

'''
Download the assemblies metadata and fasta files 
using NCBI's dataset command-line tool (v16.x)
'''
import argparse
import subprocess
import os

def download_metadata(organism):
    print(organism)
    # Print a data report containing genome metadata
    command = f"datasets summary genome taxon \"{organism}\" > assemblies_metadata.json"
    subprocess.run(command, shell=True, check=True)


def download_genomes(organism):
    # Download a dehydrated zip archive including the data report and locations of data files (use the rehydrate command to retrieve data files).
    command = f"datasets download genome taxon '{organism}' --dehydrated"
    subprocess.run(command, shell=True, check=True)

def unzip_genomes():
    command = f"unzip ncbi_dataset.zip -d assemblies_fasta"
    subprocess.run(command, shell=True, check=True)
    
def rehydrate_genomes():
    # Download data files for a dehydrated genome data package.
    command = f"datasets rehydrate --directory assemblies_fasta/"
    subprocess.run(command, shell=True, check=True)

def main():
    # Command line arguments parser
    parser = argparse.ArgumentParser()

    # Add flags for cl argumnets
    parser.add_argument("--organims_name", help = "Name organis for which assemblies will be downloaded", default= "pseudomonas aeruginosa", type=str)
    parser.add_argument("--data_dir", help = "Path to directory where the assemblies will be downloaded", default='./data/raw_data/')
    
    # Parse the arguments
    args = parser.parse_args()

    # Access the values of the arguments
    organism = args.organims_name
    data_dir = args.data_dir
    
    # Check if the output directory exists, if it doesn't exist it is created
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    
    # Move into the ./data/raw_data/ directory
    os.chdir(f"{data_dir}")

    # call the functions to download metadata and fasta files
    download_metadata(organism)
    download_genomes(organism)
    unzip_genomes()
    rehydrate_genomes()
    
    # move back to the project directory
    os.chdir('../..')

if __name__ == "__main__":
    main()
