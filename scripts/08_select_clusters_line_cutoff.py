#!/usr/bin/env python
# coding: utf-8

"""
Description:
This script selects a subset of significant clusters with an over-representation of H sequences
based on the linear equation provided with the --line_equation flag.
"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import argparse
import os
from Bio import SeqIO
from tqdm import tqdm

def set_theme():
    """
    Set the theme for the plots.
    """
    sns.set_style("ticks") # white background with ticks
    sns.set_context("paper") # control the size of the plot
    my_colors = ['#8dd3c7','#f2d140','#bebada','#fb8072','#80b1d3','#fdb462']
    sns.set_palette(my_colors)

def parse_line_equation(equation_str):
    """
    Check if the equation is in the correct format and extract the m and q values.
    """
    print('Parsing the line equation...\n')

    # convert the equation string to a lowercase
    equation_str = equation_str.lower()

    # Check if the equation is in the correct format
    if equation_str.count('x') != 1:
        print(f"The equation must contain one x variable. The equation provided is: {equation_str}\n")
        exit()
    
    if equation_str.count('=') != 1:
        print(f"The equation must contain one equal sign. The equation provided is: {equation_str}\n")
        exit()
    
    if equation_str.count('y') != 1:
        print(f"The equation must contain one y variable. The equation provided is: {equation_str}\n")
        exit()

    error_message = f"""The format of the equation {equation_str} is not correct.\n
              Please provide a linear equation in the form: y = mx + q where m and q are numbers.\n
              Make that the values are separated by spaces.\n
              Fractions are not allowed.\n
              If there are decimal values, make sure to use a dot as a decimal separator.\n"""
    
    # Remove y, = and spaces from the equation
    processed_eq= equation_str.replace('y', '').replace('=', '').replace(' ', '')
    
    # Extract the m and q values
    m = processed_eq.split('x')[0].strip(' ') # split the equation at x and take the part before x
    q = processed_eq.split('x')[1] # split the equation at x and take the part after x
    
    # Check if the q value is in the correct format
    if q.count('+') > 1 or q.count('-') > 1:
        print("The equation can contain only one sign before the q value.")
        print('The accepted format is: y = mx + q \n The equation provided is: {equation_str}\n')
        exit()
    
    # Check if the m value is in the correct format
    if m.count('+') > 1 or m.count('-') > 1:
        print("The equation can contain only one sign before the m value.")
        print('The accepted format is y = mx + q \n The equation provided is: {equation_str}\n')
        exit()

    # Check if the m value is a number
    try:
        m = float(m)
    except ValueError:
        print(f"m: {m}")
        print(error_message)
        exit()
    
    # Extract the sign of q and check if it is a number
    if q[0] == '+':
        try:
            q = float(q[1:]) 
        except ValueError:
            print(f"m: {m}")
            print(f"q: {q}")
            print(error_message)
            exit()

    elif q[0] == '-':
        try:
            q = -float(q[1:])
        except ValueError:
            print(f"m: {m}")
            print(f"q: {q}")
            print(error_message)
            exit()

    else:
        print(error_message)
        exit()
    
    print(f"The equation provided is: {equation_str}")
    print(f"m: {m}")
    print(f"q: {q}\n")
    
    return m, q


def extract_subset_list(dataset_csv, m, q, out_dir):
    """
    Save the clusters that fall below the cutoff line.
    """
    print('Extracting the subset of clusters...\n')

    # Read the dataset
    clusters = pd.read_csv(dataset_csv, index_col=0)
    clusters.rename(columns={'repr_seq': 'sequence_id'}, inplace=True) # rename the seq id column
    clusters['sequence_id'] = clusters['sequence_id'].str.replace('_\*$', '', regex=True) # remove the trailing '_*' from the sequence ids
    
    # Extract the subset of clusters that falls below the cutoff line 
    # dataset['tot_EN'] = y
    # dataset['tot_H'] = x
    subset = clusters[clusters['tot_EN'] <= (clusters['tot_H'] * m  +q)]

    # Save the sequence ids of the subset in a text file
    subset_list = subset['sequence_id'].tolist()
    with open(f'./{out_dir}/clusters_line_cutoff_list.txt', 'w') as f:
        for element in subset_list:
            f.write(element + '\n')
    
    # Plot the scatterplot of the dataset with the cut-off line
    set_theme()
    plt.figure(figsize=(12,7))
    sns.scatterplot(data = clusters, x = 'tot_H', y = 'tot_EN', alpha = 0.8, hue = 'reject_H0',  edgecolor = 'black')
    plt.xlabel('H count')
    plt.ylabel('EN count')
    plt.legend(title='p-value < 0.05', bbox_to_anchor=(1.05, 1), loc='upper left')
    x = np.linspace(0, 12000, 10000)
    y = m * x + q
    plt.plot(x, y, color = 'black', linestyle = '--')
    plt.savefig(f'./{out_dir}/clusters_line_cutoff.png',dpi=300, bbox_inches ="tight")

    return subset_list

def extract_subset_sequences(fasta_file, subset_list, out_dir):
    """
    Extract the sequences of the subset from the fasta file.
    """

    # Read the sequences fasta file
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    #check if the output file exists. If it does, remove it to avoid appending to the old file
    output_file = f'{out_dir}/selected_sequences.fasta'
    if os.path.exists(output_file):
        os.remove(output_file)
    print(f'Saving the subset of sequences in {output_file}\n')
    # Write the sequences to the output file
    with open(output_file, 'a') as f:
        for sequence_id in tqdm(subset_list):
            if sequence_id in sequences:
                SeqIO.write(sequences[sequence_id], f, "fasta")
            else:
                print(f"Sequence ID '{sequence_id}' not found in the input fasta file.\n")
        # Add one empty line at the end of the fasta file
        f.write('\n')
    # Remove all the * from the end of the lines
    os.system(f"sed -i 's/\*$//' {output_file}")

def main():

    # Command line arguments parser
    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--output_dir',
                        type=str, 
                        required = False,
                        default='.',
                        metavar = 'OUTPUT_DIRECTORY',
                        help='Path to the output directory.\n', )
    
    parser.add_argument('-c','--clusters_csv', 
                        type=str, 
                        metavar = 'CLUSTERS_CSV_FILE',
                        required = True,
                        help='Input csv file containing the clusters data.\n',)
   
    parser.add_argument('-l','--line_equation', 
                        type=str, 
                        required = True,
                        metavar = 'LINE_EQUATION',
                        help="""Equation of the cutoff line in the form y=mx+q\n
                                Note that the equation must NOT contain spaces.
                                Alternatively, you can provide the equation enclosed in quotes "y = mx + q\n" """, 
                        default='y=0.25x-125')
    
    parser.add_argument('-f','--fasta_file', 
                        type=str, 
                        help='Fasta file from where to extract the sequences of the selected subset of clusters.\n', 
                        metavar = 'FASTA_FILE',
                        default = '/storage/agata/pa_virulence/results/CD_HIT_output',
                        required = False)

    # Save CL arguments to a variable
    args = parser.parse_args()

    equation_str = args.line_equation
    out_dir = args.output_dir
    dataset_csv = args.clusters_csv
    fasta_file = args.fasta_file

    # In the path to the output directory, remove the trailing slash if it exists
    out_dir = out_dir.rstrip('/')

    # Create the output directory if it does not exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # Call the functions
    m, q = parse_line_equation(equation_str)
    subset_id_list = extract_subset_list(dataset_csv, m, q, out_dir)
    
    if fasta_file:
        extract_subset_sequences(fasta_file, subset_id_list, out_dir)

    print(f'\nDone!\nThe output files are saved in: {out_dir}/\n')


if __name__ == "__main__":
    main()
