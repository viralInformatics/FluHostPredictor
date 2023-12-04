#!/data/penglab3-20T/lihuiru/miniconda3/bin/python3
# -*coding: utf-8 -*-
# Author:lihuiru
# Created on 2023/11/7 10:58

import argparse
import re
import sys
from collections import defaultdict
from itertools import product
import numpy as np
import pandas as pd
import subprocess
import os
from Bio import SeqIO
from Bio import pairwise2
from pathlib import Path
from . import predict_host
import pkg_resources

pd.set_option('display.max_columns', None)

HA_TYPES = [f"H{i}" for i in range(1, 17) if i != 3]
NA_TYPES = [f"N{i}" for i in range(1, 10) if i != 2]

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB_PATH = os.path.join(base_dir, 'data', 'flu_db.dmnd')
STD_PATH = os.path.join(base_dir, 'data', 'std.fasta')
STRUCTURE_PATH = os.path.join(base_dir, 'data', 'structure')
STANDARD_PATH = os.path.join(base_dir, 'data', 'standard_seq_protein')
MODEL_PATH = os.path.join(base_dir, 'model')
DATA_PATH = os.path.join(base_dir, 'data')
MARKER_PATH = os.path.join(base_dir, 'data', 'markers_for_extract')


def run_diamond_blast(input_fasta_file, output_path, threads, evalue = 1e-5, suppress_output = False):
    cmd = [
        "diamond", "blastp",
        "-d", DB_PATH,
        "-q", input_fasta_file,
        "-o", output_path,
        "-f", '6',
        "-p", str(threads),
        "-e", str(evalue),
        "--sallseqid",
        "--salltitles"
    ]
    if suppress_output:
        subprocess.call(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    else:
        subprocess.call(cmd)
        print("\nDIAMOND BLAST completed.\n")


def read_annotation_results(output_path, threshold):
    # Read in the alignment file
    data = pd.read_csv(output_path, sep = "\t", header = None,
                       names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

    # Filter based on the evalue threshold
    data = data[data["evalue"] <= threshold]
    # Extract the highest bitscore hits
    best_hits_idx = data.groupby("qseqid")["bitscore"].idxmax()
    best_hits = data.loc[best_hits_idx, ["qseqid", "sseqid"]]
    return best_hits


def map_accession_to_protein(best_hits):
    # Load the mapping dictionary
    protein_sequences = SeqIO.parse(STD_PATH, "fasta")
    id_pro_dic = {}
    for record in protein_sequences:
        protein_type = record.id.split("_")[0]
        id_pro_dic[record.id] = protein_type
    # Merge the best hits with the mapping
    best_hits = best_hits.merge(pd.DataFrame(list(id_pro_dic.items()), columns = ["sseqid", "Protein Abbreviation"]),
                                on = "sseqid", how = "left")
    return best_hits[['qseqid', 'Protein Abbreviation']]


def update_fasta_csv_with_annotations(input_fasta_path, annotations, output_directory, prefix, suppress_output):
    if not suppress_output:
        # Read the FASTA file
        records = list(SeqIO.parse(input_fasta_path, "fasta"))
        input_fasta_filename = os.path.split(input_fasta_path)[1]

        # Create a mapping dictionary from annotations DataFrame
        annotations_dict = annotations.set_index('qseqid')['Protein Abbreviation'].to_dict()

        # Update the record IDs with annotations
        for record in records:
            protein_info = annotations_dict.get(record.id)
            if protein_info:
                # Append the protein information only to the description field
                record.id = f"{record.id}_{protein_info} {record.description.split(' ', 1)[-1]}"
                record.description = ""

        # Write the annotated sequences to a new FASTA file
        output_fasta_filename = f"{prefix}{input_fasta_filename.split('.')[0]}_annotated.fasta"

        output_fasta_path = f"{output_directory}/{output_fasta_filename}"

        with open(output_fasta_path, "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")
        print("\nFASTA file updated with annotations.")


def annotate_fasta_file(input_fasta_path, output_directory = ".", prefix = "", evalue = 1e-5,
                        update_file = True, threads = 10, suppress_output = False):
    """
    Annotate a FASTA file using DIAMOND BLAST against a flu database.

    Parameters:
        input_fasta_path(str): Path to the input FASTA file or a directory.
        output_directory (str): Directory where the output files will be saved.
        prefix (str): Prefix to be added to the output filenames.
        evalue (float): E-value threshold for filtering BLAST hits.
        update_file (bool): If True, update the FASTA file with annotations.
        threads (int): Number of parallel threads.
        suppress_output (bool): If True, suppress output.
    """
    os.makedirs(output_directory, exist_ok = True)
    input_fasta_filename = os.path.split(input_fasta_path)[1]
    add_prefix = prefix + "_" if prefix else ""
    output_filename = f"{add_prefix}{input_fasta_filename.split('.')[0]}.aln"
    output_path = f"{output_directory}/{output_filename}"
    print(output_path)

    # Run DIAMOND BLAST
    run_diamond_blast(input_fasta_path, output_path, threads, evalue, suppress_output)

    # Read and process the BLAST results
    best_hits = read_annotation_results(output_path, evalue)
    annotations = map_accession_to_protein(best_hits)

    # Update the CSV file with annotations
    output_csv_filename = f"{add_prefix}{input_fasta_filename.split('.')[0]}_annotated.csv"
    output_csv_path = f"{output_directory}/{output_csv_filename}"
    annotations.to_csv(output_csv_path, index = False)
    print("CSV file updated with annotations.\n")

    if update_file:
        # Update the FASTA file with annotations
        update_fasta_csv_with_annotations(input_fasta_path, annotations, output_directory, add_prefix, suppress_output)

    return annotations


def load_markers(filepath):
    """
    Load and process markers from an input file.

    Parameters:
        filepath: Path to the Excel file containing host markers.

    Returns:
        Dictionary with protein types as keys and lists of host markers as values.
    """
    column_names = ['Protein Type', 'Amino acid site']
    data = pd.read_csv(filepath)
    data = data.dropna(how = "all", axis = 1)
    data.columns = column_names + data.columns[len(column_names):].tolist()
    data["Amino acid site"] = data["Amino acid site"].str.split('(', expand = True)[0]
    # data["Specific Type"] = data["Protein Type"].str.rsplit("_", n = 1).str[-1]
    # data['Protein Type'] = data['Protein Type'].str.replace(r'_\d+$', '', regex = True)
    return data.groupby('Protein')['Amino acid site'].apply(lambda x: list(set(x))).to_dict(), data


def map_residues_to_h3(protein, marker_dict, convert_to_h3_dict):
    """
    Maps the residue numbers for a given protein to the H3/N2 numbering system.

    Parameters:
        protein (str): The protein identifier.
        marker_dict (dict): Dictionary containing markers for various proteins.
        convert_to_h3_dict (dict): Dictionary that maps residue numbers to H3.

    Returns:
        list: A list of residues mapped to H3 numbering system.
    """
    markers = marker_dict[protein]

    # 如果markers是字符串，将其转换为只含一个元素的列表
    if isinstance(markers, str):
        markers = [markers]

    mapped_residues = []
    for marker in markers:
        # Ensure the marker is in the expected format (e.g., "12A")
        marker_match = re.fullmatch(r"(\d+)([A-Z])", marker)
        if not marker_match:
            # print(f"Warning: marker '{marker}' is not in the correct format and will be skipped.")
            continue

        position, amino_acid = marker_match.groups()
        h3_position = convert_to_h3_dict.get(position)
        if h3_position is None:
            # print(f"Warning: Position {position} does not have an H3 mapping "
            #       f"in the structure comparison file and will be skipped.")
            continue

        mapped_residues.append(h3_position + amino_acid)

    return mapped_residues


def convert_HA_residues(marker_dict, structure_folder):
    """
    Converts HA/NA residues to H3/N2 numbering.

    Parameters:
        marker_dict: Dictionary with protein types as keys and marker lists as values.
        structure_folder: Folder path where the structure mapping files are located.

    Returns:
        Updated marker_dict with HA/NA types converted to H3/N2 numbering.
    """
    updated_marker_dict = marker_dict.copy()  # Create copy
    for protein in list(marker_dict.keys()):
        if protein in HA_TYPES:
            mapping_data = pd.read_csv(f"{structure_folder}/H3_{protein}.txt", sep = "\t", header = None,
                                       names = ['H3', protein])
            convert_to_h3_dict = dict(zip(mapping_data[protein], mapping_data['H3']))

            residues = map_residues_to_h3(protein, marker_dict, convert_to_h3_dict)
            # residues = [convert_to_H3_dict.get(re.search(r"\d+", i).group()) +
            # re.search(r"[A-Z]", i).group() for i in
            #             marker_dict[protein] if convert_to_H3_dict.get(re.search(r"\d+", i).group())]
            if "H3" in updated_marker_dict:
                updated_marker_dict["H3"].extend(residues)
            else:
                updated_marker_dict["H3"] = residues
            del updated_marker_dict[protein]  # del key
        elif protein in NA_TYPES:
            if os.path.isfile(f"{structure_folder}/N2_{protein}.txt"):
                mapping_data = pd.read_csv(f"{structure_folder}/N2_{protein}.txt", sep = "\t", header = None,
                                           names = ['N2', protein])
            else:
                mapping_data = pd.read_csv(f"{structure_folder}/{protein}_N2.txt", sep = "\t", header = None,
                                           names = [protein, 'N2'])
            convert_to_n2_dict = dict(zip(mapping_data[protein], mapping_data['N2']))

            residues = map_residues_to_h3(protein, marker_dict, convert_to_n2_dict)
            if "N2" in updated_marker_dict:
                updated_marker_dict["N2"].extend(residues)
            else:
                updated_marker_dict["N2"] = residues
            del updated_marker_dict[protein]  # del key

    return updated_marker_dict


def annotate_markers(adaptability_path, receptor_path, use_receptor = True):
    """
    Annotate markers by loading adaptability and receptor data, then converting HA types.

    Parameters:
        adaptability_path: Path to the Excel file with adaptability markers.
        receptor_path: Path to the Excel file with receptor markers.
        structure_folder: Folder path where the structure mapping files are located.
        use_receptor:Use receptor markers.

    Returns:
        A dictionary with annotated markers.
    """
    # Load markers from files
    adapt_dict, data1 = load_markers(adaptability_path)
    recep_dict = {}
    if use_receptor:
        recep_dict,data2 = load_markers(receptor_path)

    # Combine the dictionaries
    marker_dict = {**adapt_dict, **recep_dict}

    # Duplicated
    marker_dict = {i: list(set(j)) for i, j in marker_dict.items()}

    # Convert HA/NA residues to H3/N2 numbering and update marker_dict
    marker_dict = convert_HA_residues(marker_dict, STRUCTURE_PATH)

    data = pd.concat([data1,data2])
    # Duplicated
    marker_dict = {i: list(set(j)) for i, j in marker_dict.items()}

    return marker_dict, data


def renumber_sequence(best_alignment):
    """
    Renumber a protein sequence based on the best alignment result.

    Parameters:
        best_alignment (list of tuples): The alignment result between the standard and query sequences.

    Returns:
        list: A list of renumbered positions in the format 'positionamino_acid'.
    """
    # Initialize the list for storing renumbered positions
    renumbered_positions = []
    count = 1  # Start counting positions from 1
    for std_char, query_char in zip(best_alignment[0], best_alignment[1]):
        if std_char != '-':  # Ignore gaps in the standard sequence
            renumbered_positions.append(f"{count}{query_char}")
            count += 1  # Increment the position counter for non-gap characters

    return renumbered_positions


def renumber_proteins(fasta_path, acc_pro_dict, marker_dict):
    """
    Perform global alignment of protein sequences from a FASTA file against standard sequences and renumber them.
    Specifically renumber the positions of proteins corresponding to H1-H18 (excluding H3)
    based on their subtype standard sequences.

    Parameters:
        fasta_path (str): Path to the FASTA file containing protein sequences.
        acc_pro_dict (dict): Dictionary mapping accession IDs to protein abbreviations.
        marker_dict (dict): Dictionary with protein markers.

    Returns:
        dict: A dictionary with protein IDs and their renumbered positions.
    """
    fasta_sequences = SeqIO.parse(fasta_path, 'fasta')
    renumbering_results = {}

    for record in fasta_sequences:
        HA_results = {}
        protein_id = record.id
        protein_abbr = acc_pro_dict.get(protein_id)
        is_hana_type = protein_abbr in HA_TYPES or protein_abbr in NA_TYPES
        if protein_abbr in marker_dict or is_hana_type:
            try:
                # Construct the path to the standard sequence file
                standard_seq_path = os.path.join(STANDARD_PATH, f"{protein_abbr}.fas")
                standard_seq = next(SeqIO.parse(standard_seq_path, 'fasta')).seq

                # Perform global alignment
                alignments = pairwise2.align.globalxx(standard_seq, record.seq)
                best_alignment = max(alignments, key = lambda x: x.score)

                if is_hana_type:
                    # Store the renumbered sequence for HA/NA types
                    HA_results[protein_abbr] = renumber_sequence(best_alignment)
                else:
                    # Store the renumbered sequence for non-HA/NA types
                    renumbering_results[protein_id] = renumber_sequence(best_alignment)

            except Exception as e:
                print(f"An error occurred while processing {protein_id}: {str(e)}")
        else:
            print(f"No markers found for {protein_abbr} in the source data.")

        # Convert other HA subtype numbering to H3
        if is_hana_type:
            renumbered_positions = convert_HA_residues(HA_results, STRUCTURE_PATH)

            pop_num = "H3" if protein_abbr in HA_TYPES else ("N2" if protein_abbr in NA_TYPES else None)
            renumbered_positions[protein_id] = renumbered_positions.pop(pop_num)
            renumbering_results.update(renumbered_positions)

    return renumbering_results


def merge_dicts_with_list(dict_list):
    """
    Function to merge a list of dictionaries. If keys are repeated, values are merged into a list.

    Parameters:
    - dict_list (list): List containing dictionaries.

    Returns:
    - Merged dictionary.
    """
    merged_dict = {}
    for d in dict_list:
        for key, value in d.items():
            if key in merged_dict:
                # Merge values into a list if the key already exists
                if not isinstance(merged_dict[key], list):
                    merged_dict[key] = [merged_dict[key]]
                merged_dict[key].append(value)
            else:
                # Directly add if the key does not exist
                merged_dict[key] = value
    return merged_dict


def generate_combinations(group):
    """
    Groups mutations by specific types and generates all possible combinations.

    Parameters:
    - group (DataFrameGroupBy): Grouped DataFrame.

    Returns:
    - List of all possible combinations.
    """
    # Group mutations by specific type
    spec_type_groups = group.groupby('Specific Type')
    # Create a dictionary to store lists of mutations mapped to proteins by specific type
    mutation_dict = defaultdict(list)
    for spec_type, g in spec_type_groups:
        for _, row in g.iterrows():
            if type(row["Amino acid site"]) == str:
                mutation_dict[spec_type].append({row['Protein']: row['Amino acid site'].strip()})
            else:
                mutation_dict[spec_type].append({row['Protein']: row['Amino acid site']})
    # Extract lists of dictionaries for each key
    values_lists = [mutation_dict[key] for key in mutation_dict]
    # Generate all possible combinations and merge dictionaries
    combinations = [merge_dicts_with_list(comb) for comb in product(*values_lists)]
    return combinations


def generate_protein_dict(grouped_data):
    """Generate protein dictionary"""
    new_protein_dict = defaultdict(list)
    for name, group in grouped_data:
        if 'combination' in name:
            new_protein_dict[name] = generate_combinations(group)
        else:
            new_protein_dict[name].extend(
                {row['Protein']: row['Amino acid site'].strip() if isinstance(row['Amino acid site'], str) else row[
                    'Amino acid site']}
                for _, row in group.iterrows()
            )
    return new_protein_dict


def load_total_markers(data):
    data["Specific Type"] = data["Protein Type"].str.rsplit("_", n = 1).str[-1]
    data['Protein Type'] = data['Protein Type'].str.replace(r'_\d+$', '', regex = True)
    return data.groupby('Protein Type')


def is_subset_complex(dict1, dict2):
    """
    Check if one dictionary is a complex subset of another.

    Parameters:
    - dict1, dict2: Dictionaries to be compared.

    Returns:
    - Boolean: True if dict1 is a subset of dict2, else False.
    """
    for key, value1 in dict1.items():
        if key not in dict2:
            return False

        value2 = dict2[key]

        if isinstance(value1, list) and isinstance(value2, list):
            if not set(value1).issubset(set(value2)):
                return False
        elif isinstance(value1, str) and isinstance(value2, list):
            if value1 not in value2:
                return False
        elif value1 != value2:
            return False

    return True


def format_marker(marker, protein_prefix = ''):
    """
    Formats a single genetic marker. If the marker contains a hyphen ('-'),
    only the part before the hyphen is retained and appended with 'Deletion'.
    If a protein prefix is provided, it's added before the marker.

    Parameters:
        marker (str): The genetic marker to be formatted.
        protein_prefix (str): An optional prefix to be added before the marker.

    Returns:
        str: Formatted genetic marker.
    """
    # Check if the marker contains a hyphen and split accordingly.
    if '-' in marker:
        amino_acid = marker.split('-')[0]
        deletion_suffix = "Deletion"
    else:
        amino_acid = marker
        deletion_suffix = ""

    # Combine the protein prefix, amino acid, and deletion suffix.
    formatted_marker = f"{protein_prefix}-{amino_acid}{deletion_suffix}" \
        if protein_prefix else f"{amino_acid}{deletion_suffix}"
    return formatted_marker


def format_marker_list(markers, protein_prefix = ''):
    """
    Formats a list of markers or a single marker string.
    In case of a list where all elements contain a hyphen, a special formatted string is returned.
    Otherwise, each marker in the list is formatted individually.

    Parameters:
        markers (str or list): A string or list of strings representing genetic markers.
        protein_prefix (str): An optional prefix to be added before each marker.

    Returns:
        str: A single string representing the formatted markers, joined by '&'.
    """
    # Check if the input is a single string and format directly.
    if isinstance(markers, str):
        return format_marker(markers)

    # Determine if all markers in the list contain a hyphen.
    all_contain_dash = all('-' in marker for marker in markers)
    if all_contain_dash:
        # Create a special format string if all markers contain a hyphen.
        start = markers[0].split('-')[0]
        end = markers[-1].split('-')[0]
        return f"{start}-{end}CompleteDeletion"

    # Format each marker individually and join with '&'.
    return '&'.join(format_marker(marker, protein_prefix) for marker in markers)


def process_dictionary(data_dict):
    """
    Processes a dictionary containing genetic markers.
    If the dictionary has a single key-value pair, the value is formatted directly.
    For multiple key-value pairs, each is formatted separately and joined by '&'.

    Parameters:
        data_dict (dict): A dictionary with protein names as keys and genetic markers as values.

    Returns:
        str: A single string representing the formatted contents of the dictionary.
    """
    # Process a single key-value pair directly.
    if len(data_dict) == 1:
        return format_marker_list(next(iter(data_dict.values())))

    # Format each key-value pair separately if there are multiple.
    return '&'.join(format_marker_list(markers, protein) for protein, markers in data_dict.items())


def process_protein_sequence(acc_id, renumbered_position, acc_pro_dic, marker_markers):
    protein_type = acc_pro_dic[acc_id]

    # Skip processing if protein type is unknown
    if protein_type == "Unknown":
        return None, None

    use_protein = "H3" if protein_type in HA_TYPES else protein_type
    expected_markers = marker_markers.get(use_protein, [])
    markers = []

    for marker in expected_markers:
        match = re.match(r"(\d+)([A-Z])", marker)
        if match and match.group() in renumbered_position:
            markers.append(match.group())

    protein = f'H3' if protein_type in HA_TYPES else (
        f'N2' if protein_type in NA_TYPES else protein_type)

    return protein, markers


def check_marker_combinations(total_markers, results_markers, markers_type, input_file_name, data, ha_type, na_type):
    results = []

    # Initialize results with empty/default values
    initial_data = {
        'Strain ID': '',  # or some default value
        'Amino acid site': '',  # or some default value
        'Protein Type': ''  # or some default value
    }
    results.append(initial_data)

    # Sequentially check if each type of marker for the phenotype is present in the identified marker dictionary.
    for marker_protein_type, marker_list in total_markers.items():
        # proba_comb is one of the multiple combinations of markers for each type, which is a dictionary.
        # 'combination-combination_449': [{'PB2': '158G', 'PA': '295P'}, {'PB2': '158A', 'PA': '295P'}]
        # proba_comb is one of these dictionaries, if the dictionary is satisfied,
        # it is considered that there is a combination-combination_449 type of marker.
        for proba_comb in marker_list:
            # If the key-value pair in this dictionary exists in the identified marker dictionary,
            # return a more concise format.
            if is_subset_complex(proba_comb, results_markers):
                if proba_comb and all(proba_comb.values()):
                    markers_formated = process_dictionary(proba_comb)
                    results.append({
                        'Strain ID': input_file_name.split(".")[0],
                        'Amino acid site': markers_formated,
                        'Protein Type': marker_protein_type,
                    })

    results = pd.DataFrame(results)
    final_results = merge_dataframes(results, data, markers_type, ha_type, na_type)
    return final_results


def merge_dataframes(results, data, markers_type, ha_type, na_type):
    """
    Merge two DataFrames based on the 'Protein Type' column, handling rows with and without 'combination' differently.

    Args:
    results (pd.DataFrame): DataFrame containing results data.
    data (pd.DataFrame): DataFrame containing additional data to merge.
    markers_type (str): String indicating the type of markers to be used for renaming columns.

    Returns:
    pd.DataFrame: The merged DataFrame after processing.
    """
    # Pre-compile the regex pattern for performance
    combination_pattern = re.compile(r'combination')

    # Split the 'data' DataFrame into two based on 'combination' presence without using str.contains for each row
    data['HasCombination'] = data['Protein Type'].apply(lambda x: bool(combination_pattern.search(x)))
    data_with_combination = data[data['HasCombination']].drop_duplicates(subset = 'Protein Type')
    data_without_combination = data[~data['HasCombination']]

    # Drop the 'Amino acid site' column and the helper 'HasCombination' column to avoid merging them later
    data_with_combination = data_with_combination.drop(columns = ["Amino acid site", "HasCombination"])

    # Do the same split for the 'results' DataFrame
    results['HasCombination'] = results['Protein Type'].apply(lambda x: bool(combination_pattern.search(x)))
    results_with_combination = results[results['HasCombination']]
    results_without_combination = results[~results['HasCombination']]

    # Merge parts with and without 'combination' separately
    merged_with_combination = pd.merge(results_with_combination, data_with_combination, on = 'Protein Type',
                                       how = 'left')
    merged_without_combination = pd.merge(results_without_combination, data_without_combination,
                                          on = ['Protein Type', 'Amino acid site'], how = 'left')

    # Concatenate the merged parts
    final_results = pd.concat([merged_with_combination, merged_without_combination], ignore_index = True)

    # Rename the 'Amino acid site' column and cleanup 'Protein Type'
    final_results.rename(columns = {'Amino acid site': f'{markers_type.title()} Markers'}, inplace = True)
    final_results['Protein Type'] = final_results['Protein Type'].str.replace('-combination.*', '', regex = True)

    final_results['Protein Type'] = final_results['Protein Type'].apply(lambda x: get_hana_string(x, ha_type, na_type))
    # Drop unnecessary columns and the helper 'HasCombination' column
    final_results.drop(
        columns = ["Specific Type", "Protein", "HasCombination_x", "HasCombination_y", "HasCombination"],
        inplace = True)

    # Replace empty strings with NaN
    final_results.replace('', np.nan, inplace = True)

    # Drop rows where specific columns are NaN
    final_results.dropna(subset = ['Strain ID', 'Adaptation Markers', 'Protein Type'], how = "all", inplace = True)

    final_results.drop_duplicates(inplace = True)
    return final_results


def get_hana_string(protein_type, ha_type, na_type):
    if protein_type in HA_TYPES or protein_type == "H3":
        return f'{ha_type}(H3 numbering)'
    elif protein_type in NA_TYPES or protein_type == "N2":
        return f'{na_type}(N2 numbering)'
    else:
        return protein_type


def identify_markers(input_file_path, renumbering_results, marker_markers, acc_pro_dic, markers_type, data,
                     output_directory = ".", prefix = ""):
    """
    Identifies host markers in protein sequences based on the provided marker markers
    and the renumbered sequences.

    Parameters:
    input_file_path (str): The path to the input file containing sequence data.
    renumbering_results (dict): Dictionary with accession ID as keys and renumbered positions as values.
    marker_markers (dict): Dictionary defining the host markers to be identified.
    acc_pro_dic (dict): Dictionary mapping accession IDs to protein types.
    markers_type (str): Type of markers to be identified (e.g., 'HA', 'NA').
    data (DataFrame): DataFrame containing additional data for marker identification.
    output_directory (str, optional): The directory where output files will be saved. Defaults to the current directory.
    prefix (str, optional): A prefix to be added to the output file names.

    Returns:
    DataFrame: A DataFrame with identified markers for each sequence.
    """
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok = True)
    input_file_name = os.path.split(input_file_path)[1]
    results_markers = defaultdict(list)

    # Process each accession ID and its renumbered position
    for acc_id, renumbered_position in renumbering_results.items():
        protein, markers = process_protein_sequence(acc_id, renumbered_position, acc_pro_dic, marker_markers)
        if protein:
            results_markers[protein] = markers

    # Initialize types
    ha_type = na_type = None

    # Identify HA and NA types based on acc_pro_dic
    for acc, pro in acc_pro_dic.items():
        if pro in HA_TYPES or pro == "H3":
            ha_type = pro
        elif pro in NA_TYPES or pro == "N2":
            na_type = pro

    # This is to handle each HA/NA, including those present in combinations,
    # the file only processed single HA/NA markers
    ori_markers = generate_protein_dict(load_total_markers(data))
    total_markers = defaultdict(list)
    for pro, lst in ori_markers.items():
        for dic in lst:
            if dic and all(dic.values()):
                # After converting through convert_HA_residues, everything will become H3, so there's no impact
                total_markers[pro].append(convert_HA_residues(dic, STRUCTURE_PATH))

    # Check marker combinations and merge results with data
    results_df = check_marker_combinations(total_markers, results_markers, markers_type,
                                           input_file_name, data, ha_type, na_type)

    # Add prefix to filename if provided
    add_prefix = prefix + "_" if prefix else ""
    filename = add_prefix + input_file_name.split(".")[0] + "_markers.csv"
    # Save the results to a CSV file
    results_df.to_csv(f"{output_directory}/{filename}", index = False)

    return results_df


def is_fasta_file(filename):
    # Check if input is a fasta file
    return re.search(r'\.(fasta|faa|fa)$', filename, re.IGNORECASE)


def find_files_with_string(directory, string):
    all_items = os.listdir(directory)
    files_with_string = [item for item in all_items
                         if string in item and os.path.isfile(os.path.join(directory, item))
                         and item.endswith("_annotated.csv")]

    return files_with_string[0]


def parse_args():
    parser = argparse.ArgumentParser(prog = 'fluhp',
                                     description = 'fluhp command line tool for flu marker '
                                                   'extraction, annotation and host prediction.')
    subparsers = parser.add_subparsers(dest = 'subcommand', help = 'Sub-commands')

    # anno subcommand
    anno_parser = subparsers.add_parser('anno',
                                        help = 'Annotate a FASTA file or all FASTA files in a directory '
                                               'using DIAMOND BLAST against a flu database.')
    anno_parser.add_argument('-i', '--input', required = True,
                             help = 'Input FASTA file or directory containing FASTA files.')
    anno_parser.add_argument('-o', '--output_directory', type = str, default = '.',
                             help = 'Directory to save the output files. Defaults to the current directory.')
    anno_parser.add_argument('-p', '--prefix', type = str, default = '', help = 'Prefix for the output filenames.')
    anno_parser.add_argument('-e', '--evalue', type = float, default = 1e-5,
                             help = 'E-value threshold for DIAMOND BLAST hits. Defaults to 1e-5.')
    anno_parser.add_argument('-u', '--update_file', action = 'store_true',
                             help = 'If set, updates the FASTA file with annotations.')
    anno_parser.add_argument('-t', '--threads', type = int, default = 10,
                             help = 'Number of threads for DIAMOND BLAST. Defaults to 10.')

    # extract subcommand
    extract_parser = subparsers.add_parser('extract', help = 'Extract and process protein annotations.')
    extract_parser.add_argument('-i', '--input', required = True,
                                help = 'Input FASTA file or directory containing FASTA files.')
    extract_parser.add_argument('-a', '--anno_path', required = True,
                                help = 'Input annotation CSV file or directory containing annotation CSV files.')
    extract_parser.add_argument('-r', '--receptor_markers', action = 'store_true',
                                help = 'Path to the optional Excel file with receptor markers. If not provided, '
                                       'only adaptability markers are used.')
    extract_parser.add_argument('-o', '--output_directory', type = str, default = '.',
                                help = 'Directory to save the output files. Defaults to the current directory.')
    extract_parser.add_argument('-p', '--prefix', type = str, default = '', help = 'Prefix for the output filenames.')

    # pred command
    pred_parser = subparsers.add_parser('pred', help = 'Predict new data labels using a trained model.')
    pred_parser.add_argument('-i', '--input', required = True, type = str,
                             help = 'Input CSV file with marker data or directory containing such files.')
    pred_parser.add_argument('-m', '--model_path', default = MODEL_PATH+'/gnb_model.joblib', type = str,
                             help = 'Path to the trained model file.')
    pred_parser.add_argument('-th', '--threshold_path', default = MODEL_PATH+'/optimal_threshold.joblib', type = str,
                             help = 'Path to the file containing the optimal threshold value.')
    pred_parser.add_argument('-f', '--top_features_path', default = MODEL_PATH+'/top_features.joblib', type = str,
                             help = 'Path to the file containing the top features.')
    pred_parser.add_argument('-o', '--output_directory', type = str, default = '.',
                             help = 'Directory to save the prediction results. Defaults to the current directory.')
    pred_parser.add_argument('-p', '--prefix', type = str, default = '',
                             help = 'Prefix for the output filenames of the predictions.')

    return parser.parse_args()


def process_anno_cmd(input_file, args):
    """
    Call the appropriate functions to process a single fasta file
    """

    annotate_fasta_file(
        str(input_file),
        args.output_directory,
        args.prefix,
        args.evalue,
        args.update_file,
        args.threads
    )


def process_extract_cmd(input_file, args, is_directory = True):
    """
        Call the appropriate functions to process a single fasta file
    """
    input_filename_pre = os.path.split(input_file)[1].split('.')[0]
    if is_directory:
        anno_filename = find_files_with_string(args.anno_path, input_filename_pre)
        annotations = pd.read_csv(f"{args.anno_path}/{anno_filename}")
    else:
        annotations = pd.read_csv(f"{args.anno_path}")
    acc_pro_dic = dict(zip(annotations.iloc[:, 0], annotations.iloc[:, 1]))
    marker_dict, data = annotate_markers(MARKER_PATH + "/human_adaptation_formated.csv",
                                         MARKER_PATH + "/receptor_binding_alteration_formated.csv",
                                         args.receptor_markers)
    renumbering_results = renumber_proteins(
        fasta_path = str(input_file),
        acc_pro_dict = acc_pro_dic,
        marker_dict = marker_dict
    )
    results_df = identify_markers(
        input_file_path = str(input_file),
        renumbering_results = renumbering_results,
        marker_markers = marker_dict,
        acc_pro_dic = acc_pro_dic,
        output_directory = args.output_directory,
        prefix = args.prefix,
        markers_type = 'adaptation',
        data = data,
    )

    print("\nMarker extracted and saved to file.")


def process_directory(directory, args):
    for file in directory.iterdir():
        if is_fasta_file(str(file)):
            if args.subcommand == 'anno':
                process_anno_cmd(file, args)
            elif args.subcommand == 'extract':
                process_extract_cmd(file, args)


def process_single_file(file, args):
    if args.subcommand == 'anno':
        process_anno_cmd(file, args)
    elif args.subcommand == 'extract' and str(args.anno_path).endswith("_annotated.csv"):
        process_extract_cmd(file, args, is_directory = False)


def run_other_subcommand(args):
    input_path = Path(args.input)
    if input_path.is_dir():
        if (args.subcommand == "extract" and Path(args.anno_path).is_dir()) or (args.subcommand == "anno"):
            process_directory(input_path, args)
        else:
            print(f"Error: {args.anno_path} is not a valid directory")
    elif input_path.is_file():
        process_single_file(input_path, args)
    else:
        print(f"Error: {args.input} is not a valid file or directory", file = sys.stderr)


def main():
    args = parse_args()
    if args.subcommand == 'pred':
        predictions = predict_host.predict_new_data(
            str(Path(args.input)),
            args.model_path,
            args.threshold_path,
            args.top_features_path,
            args.output_directory,
            args.prefix
        )
        print(predictions)
        print(f"Predictions completed.")
    else:
        run_other_subcommand(args)


if __name__ == '__main__':
    main()
