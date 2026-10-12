#!/usr/bin/env python
# coding: utf-8

"""
This script converts the json file with the metadata 
for genome assemblies downloaded by NCBI into a csv file.
Metadata in the output file: accession number, organism name,
strain, isolation source, host, assembly level. 
"""

import json
import pandas as pd
import argparse
import os

def json_to_csv_parser(input_json, output_dir):
	"""
	Parse NCBI's json metadata to tsv
	"""
	with open(input_json, 'r') as f: 
		metadata = json.load(f)

	# definig the metadata (columns) that will be saved for each genome assembly (rows)
	col_names = ['assembly_accession', 'organism_name', 'strain', 'isolation_source', 'host', 'assembly_level']
	metadata_df = pd.DataFrame(columns = col_names)
	print(len(metadata['reports']))
	for i in range(len(metadata['reports'])):

		current_assembly = metadata['reports'][i] # access the asssembly
		current_df_row = ['NA']*len(col_names) # create a list of NAs to save the metadata for the current assembly

		# replace each NA with the corresponsing metadata, if available 
		current_df_row[0] = (current_assembly['accession'])
		current_df_row[1] = (current_assembly['assembly_info']['biosample']['description']['organism']['organism_name'])

		for attribute in current_assembly['assembly_info']['biosample']['attributes']:

			if 'strain' in list(attribute.values()):
				current_df_row[2] = attribute['value']

			if 'isolation_source' in list(attribute.values()):
				current_df_row[3] = attribute['value'].lower()

			if 'host' in list(attribute.values()):
				current_df_row[4] = attribute['value'].lower()
		
		current_df_row[5] = (current_assembly['assembly_info']['assembly_level'])
		
		# once the list has been filled, add it to the df
		metadata_df.loc[len(metadata_df)] = current_df_row
	
	# Check if the output directory exists, if it doesn't exist it is created
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	metadata_df.to_csv(f'{output_dir}assemblies_metadata.tsv', sep='\t')

def main():
	
	# initiate argument parsing
	parser = argparse.ArgumentParser()
	# add flags
	parser.add_argument('--metadata_json', help = "Path to .json file with assemblies metadata", default='./data/raw_data/assemblies_metadata.json')
	parser.add_argument('--out_dir', help = " Path to output directory", default='./results/')
	
	# add flag input to varaible
	args = parser.parse_args()

	# access the variables from the command line
	input_json = args.metadata_json
	output_dir = args.out_dir

	# parse the json file to csv
	json_to_csv_parser(input_json, output_dir)

if __name__ == "__main__":
	main()
