# !/usr/bin/python
# -*- coding: utf-8 -*-
# Author:lhr

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from itertools import product
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import pairwise2, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from pandas.core.common import SettingWithCopyWarning
from .blastSeq import blastSeq
from .changeStandardNameForMPNS import refreshStandardName
from .getBlastMostCommonHitProteinType import getMostCommonHitProtein
from .predict_host import predict_new_data as predh
from .standardize import standardize
from .translateDNA2Protein import translate, makeProteinFileForDownload
import warnings

warnings.filterwarnings('ignore', category=SettingWithCopyWarning)
pd.set_option('display.max_columns', None)

HA_TYPES = [f"H{i}" for i in range(1, 19) if i != 3]
NA_TYPES = [f"N{i}" for i in range(1, 10) if i != 2]
length_diffs = {'H1': 17, 'H10': 17, 'H11': 16, 'H12': 17, 'H13': 18, 'H14': 17, 'H15': 18, 'H16': 19, 'H17': 18,
                'H18': 14, 'H2': 15, 'H3': 16, 'H4': 16, 'H5': 12, 'H6': 16, 'H7': 18, 'H8': 17, 'H9': 18}
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB_PATH = os.path.join(base_dir, 'data', 'flu_db.dmnd')
STD_PATH = os.path.join(base_dir, 'data', 'std.fasta')
COMPLETE_STD_PATH = os.path.join(base_dir, 'data', 'complete_std.fasta')
STRUCTURE_PATH = os.path.join(base_dir, 'data', 'HA_NA_mapdir')
STANDARD_PATH = os.path.join(base_dir, 'data', 'standard_seq_protein')
MODEL_PATH = os.path.join(base_dir, 'model')
DATA_PATH = os.path.join(base_dir, 'data')
MARKER_PATH = os.path.join(base_dir, 'data', 'markers_for_extract')
RESULT_PATH = os.path.join(base_dir, 'result')
TEMP_PATH = 'temp/'
sub_dic = {"human_adaptation": {"479M": "140M", "389N": "46N", "453R": "127R"}}

def calculate_sort_key(line):
    parts = line.split('\t')
    try:
        value = parts[2]
        if "Max" in value:
            return 10000000000 + len(line)
        elif "Min" in value:
            return -100000 + len(line)
        else:
            return float('{:.10f}'.format(float(value)))
    except IndexError:
        return float('inf')


def format_line_stage1(line):
    parts = line.split('\t')
    try:
        value = parts[2]
        if "Max" not in value and "Min" not in value:
            parts[2] = '{0:.2e}'.format(float(value))
        return '\t'.join(parts)

    except IndexError:
        return line


def format_line_stage2(line):
    parts = line.split('\t')
    try:
        value = parts[2]
        if "Max" not in value and "Min" not in value:
            value_float = float(value)
            prefix = "Similar_" if value_float <= 1.0 else "Variant_"
            parts[2] = prefix + value
        elif "Min" in value:
            parts[2] = "similar_" + value
        elif "Max" in value:
            parts[2] = "variant_" + value
        return '\t'.join(parts)
    except IndexError:
        return line


def calculate_reassortment_risk(predictedProteinType, Host, similarVirusIDs = None):

    if similarVirusIDs is None:
        virus_id_dict = {}
    else:
        virus_id_dict = {seq_id: virus_id for seq_id, virus_id in similarVirusIDs}
    excluded_proteins = {'M2', "PA-X", "PB1-F2", "NS1", "NS2", "M1", "Unknown", ''}
    valid_data = []

    for eachProtein in predictedProteinType:
        seq_id = eachProtein[0]
        protein_type = eachProtein[1]

        if protein_type not in excluded_proteins:
            host_found = False
            for eachHost in Host:
                if eachHost[0] == seq_id:
                    host = eachHost[1]
                    host_found = True
                    break
            if host_found:
                virus_id = virus_id_dict.get(seq_id, "Unknown")
                valid_data.append((seq_id, protein_type, host, virus_id))
    hosts = [item[2] for item in valid_data]
    unique_hosts = set(hosts)
    reassortment_risk = 1 if len(unique_hosts) > 1 else 0
    gene_host_virus_info = [(item[1], item[2], item[3]) for item in valid_data]

    return reassortment_risk, gene_host_virus_info


def save_tuples_to_csv(data, file_path):
    sorted_data = sorted(data, key = lambda x: float(x[1]), reverse = True)

    with open(file_path, 'w', newline = '') as file:
        writer = csv.writer(file)

        writer.writerow(['Host', 'Probability'])

        writer.writerows(sorted_data)


def save_reassortment_info_to_csv(reassortment_risk, gene_host_virus_info, file_path):
    with open(file_path, 'w', newline = '') as file:
        writer = csv.writer(file)

        writer.writerow([f'reassortment_risk:{str(reassortment_risk)}'])

        writer.writerow(['Gene', 'Host', 'Most Similar Virus ID'])

        writer.writerows(gene_host_virus_info)




def ivew_task(resultFileDir, tempFileDir, inputFilePath, updateFastaDir):
    DIR = base_dir
    os.makedirs(updateFastaDir, exist_ok = True)
    os.makedirs(tempFileDir, exist_ok = True)
    workName = os.path.basename(inputFilePath)

    inputFile, querySeqType, dic = standardize(inputFilePath, outFileDir = tempFileDir)

    dataBaseName = None
    if querySeqType == "nucleo":
        dataBaseName = "/protType_NA"
    if querySeqType == "protein":
        dataBaseName = "/protType_AA"
    querySeqFile = inputFile
    querySeqFileDir = tempFileDir
    tempFileDir = tempFileDir + "/"
    DBDir = DIR + "/18Mid/standard_seq/allProteinTypeDB"
    queryEValue = "1e-5"
    outfmt = "6"
    print(inputFile)
    outName = "querySeqToProteinTypeDB"
    blastSeq(querySeqFile, querySeqType, DBDir, dataBaseName, queryEValue, outfmt, outName, querySeqFileDir,
             tempFileDir)
    predictedProteinType, similarVirusIDs = getMostCommonHitProtein(outName, tempFileDir)
    CDSList = []
    if querySeqType == "nucleo":
        queryProteinSeqFile = translate(querySeqFile, querySeqFileDir, DIR, tempFileDir)
        querySeqFile = queryProteinSeqFile[0]
        CDSList = queryProteinSeqFile[2]
    proteinPath = makeProteinFileForDownload(tempFileDir, querySeqFile, updateFastaDir, dic, predictedProteinType)
    if querySeqType == "nucleo":
        predictedProteinType = refreshStandardName(predictedProteinType, dic)
    CDSTypeList = []
    for eachCDSType in CDSList:
        for eachType in eachCDSType[1].split('),'):
            if eachType.split('(')[0] in ['NS1', 'NS2', 'M1', 'M2', 'PB1', 'PB1-F2', 'PA-X', 'PA']:
                CDSTypeList.append(eachCDSType[0] + '_' + eachType.split('(')[0])
            else:
                CDSTypeList.append(eachCDSType[0])

    return proteinPath, os.path.join(resultFileDir, workName + ".result"), CDSList, predictedProteinType


def read_annotation_results(output_path, threshold):
    data = pd.read_csv(output_path, sep = "\t", header = None,
                       names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                "qstart", "qend", "sstart", "send", "evalue", "bitscore"])

    data = data[data["evalue"] <= threshold]

    best_hits_idx = data.groupby("qseqid")["bitscore"].idxmax()
    best_hits = data.loc[best_hits_idx, ["qseqid", "sseqid"]]
    return best_hits


def map_accession_to_protein(best_hits):
    protein_sequences = SeqIO.parse(STD_PATH, "fasta")
    id_pro_dic = {}
    for record in protein_sequences:
        protein_type = record.id.split("_")[0]
        id_pro_dic[record.id] = protein_type

    best_hits = best_hits.merge(pd.DataFrame(list(id_pro_dic.items()), columns = ["sseqid", "Protein Abbreviation"]),
                                on = "sseqid", how = "left")
    return best_hits[['qseqid', 'Protein Abbreviation']]


def update_fasta_csv_with_annotations(input_fasta_path, annotations, output_directory, prefix, suppress_output):
    if not suppress_output:

        records = list(SeqIO.parse(input_fasta_path, "fasta"))
        input_fasta_filename = os.path.split(input_fasta_path)[1]

        annotations_dict = annotations.set_index('qseqid')['Protein Abbreviation'].to_dict()

        for record in records:
            protein_info = annotations_dict.get(record.id)
            if protein_info:
                record.id = f"{record.id}_{protein_info} {record.description.split(' ', 1)[-1]}"
                record.description = ""

        output_fasta_filename = f"{prefix}{input_fasta_filename.split('.')[0]}_annotated.fasta"

        output_fasta_path = f"{output_directory}/{output_fasta_filename}"

        with open(output_fasta_path, "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")
        print("\nFASTA file updated with annotations.")


def load_markers(filepath):
    """
    Load and process markers from an input file.
    """
    column_names = ['Protein Type', 'Amino acid site']
    data = pd.read_csv(filepath, na_values = [], keep_default_na = False)
    if 'HA_Type' in data.columns:
        ha_type_column = data['HA_Type']
    else:
        ha_type_column = None

    data = data.dropna(how = 'all', axis = 1)

    if ha_type_column is not None:
        data['HA_Type'] = ha_type_column

    data.columns = column_names + data.columns[len(column_names):].tolist()
    data["Amino acid site"] = data["Amino acid site"].str.split('(', expand = True)[0]
    return data.groupby('Protein')['Amino acid site'].apply(lambda x: list(set(x))).to_dict(), data


def get_h3_dict_and_hatype(protein, marker, convert_to_h3_dict, ha_type_info_dic):
    if protein not in HA_TYPES and protein not in NA_TYPES:
        return None, None

    if protein in NA_TYPES:
        return convert_to_h3_dict, None

    ha_type = "HA2" if "HA2" in marker else "HA1"
    dict_key = ha_type

    if ha_type == "HA1":
        value = ha_type_info_dic.get(f"{protein}-{marker}", None) if ha_type_info_dic else None
        if pd.isna(value) or value is None:
            ha_type = None
    return convert_to_h3_dict[dict_key], ha_type


def adjust_position_type(position, H3_dict):
    if not H3_dict:
        return position

    first_key_type = type(next(iter(H3_dict)))

    if isinstance(position, first_key_type):
        return position
    try:
        return first_key_type(position)
    except ValueError:
        return None


def adjust_position_and_get_h3_position(marker, hatype, H3_dict, protein):
    marker_match = re.search(r"(\d+)([A-Z]|-)", marker)
    converted_sign = None
    if not marker_match:
        return None, None, hatype, converted_sign

    position, amino_acid = marker_match.groups()

    if not hatype and protein in HA_TYPES:
        minus = length_diffs[protein]
        position = str(int(position) - minus)
        hatype = "HA1"

    if H3_dict:
        adjusted_position = adjust_position_type(position, H3_dict)

        res = H3_dict.get(int(adjusted_position), H3_dict.get(adjusted_position))
        if res == "-":
            converted_sign = protein
            adj_position = hatype + "-" + str(adjusted_position) if hatype else str(adjusted_position)
            return adj_position, amino_acid, hatype, converted_sign
        return res, amino_acid, hatype, converted_sign
    else:
        if not hatype:
            hatype = "HA1"
        return f"{hatype}-{position}{amino_acid}"


def map_residues_to_h3(protein, marker_dict, convert_to_h3_dict, hatype = None, ha_type_info_dic = None):
    markers = [marker_dict[protein]] if isinstance(marker_dict[protein], str) else marker_dict[protein]
    mapped_residues = []
    unmapped_residues = []
    for marker in markers:
        if hatype:
            H3_dict = convert_to_h3_dict[hatype]
        else:
            H3_dict, hatype = get_h3_dict_and_hatype(protein, marker, convert_to_h3_dict, ha_type_info_dic)
        if H3_dict is None:
            continue
        marker = marker if marker.endswith("-") else marker.split("-")[-1]
        h3_position, amino_acid, updated_hatype, converted_sign = adjust_position_and_get_h3_position(marker, hatype,
                                                                                                      H3_dict, protein)
        if converted_sign:
            unmapped_residues.append(h3_position)
            continue
        if h3_position is None:
            continue

        hatype_prefix = f"{updated_hatype}-" if updated_hatype else ""
        mapped_residues.append(f"{hatype_prefix}{h3_position}{amino_acid}")
    return mapped_residues, unmapped_residues


def load_mapping_data(filepath, column_names):
    if os.path.isfile(filepath):
        mapping_data = pd.read_csv(filepath, sep = "\t", header = None, names = column_names, na_values = [],
                                   keep_default_na = False)
        return dict(zip(mapping_data[column_names[1]], mapping_data[column_names[0]]))
    return {}


def process_ha_type(protein, marker_dict, structure_folder, hatype, ha_type_info_dic):
    convert_to_h3_dict_ha1 = load_mapping_data(
        f"{structure_folder}/HA1/H3_{protein}.txt", ['H3', protein])
    convert_to_h3_dict_ha2 = load_mapping_data(
        f"{structure_folder}/HA2/H3_{protein}.txt", ['H3', protein])

    combined_dict = {'HA1': convert_to_h3_dict_ha1, 'HA2': convert_to_h3_dict_ha2}
    return map_residues_to_h3(protein, marker_dict, combined_dict, hatype, ha_type_info_dic)


def process_na_type(protein, marker_dict, structure_folder, hatype, ha_type_info_dic):
    filepath = f"{structure_folder}/NA/N2_{protein}.txt"
    alternate_filepath = f"{structure_folder}/NA/{protein}_N2.txt"
    if not os.path.isfile(filepath):
        filepath = alternate_filepath if os.path.isfile(alternate_filepath) else filepath

    convert_to_n2_dict = load_mapping_data(filepath, ['N2', protein])
    return map_residues_to_h3(protein, marker_dict, convert_to_n2_dict, hatype, ha_type_info_dic)


def transform_marker_dict(marker_dict):
    transformed_data = {}
    for key, values in marker_dict.items():
        if key == 'H3':
            sub_dict = {}
            for value in values:
                prefix, suffix = value.split('-', 1)
                if prefix not in sub_dict:
                    sub_dict[prefix] = []
                sub_dict[prefix].append(suffix)
            transformed_data[key] = {k: v[0] if len(v) == 1 else v for k, v in sub_dict.items()}
        else:
            if not isinstance(values, str):
                transformed_data[key] = list(set(values))
            else:
                transformed_data[key] = values
    return transformed_data


def convert_HA_residues(marker_dict, structure_folder, hatype, ha_type_info_dic = None):
    updated_marker_dict = marker_dict.copy()

    for protein in list(marker_dict.keys()):
        if protein == "H3":
            res = []
            if not isinstance(marker_dict["H3"], list):
                marker_dict["H3"] = [marker_dict["H3"]]
            for marker in marker_dict["H3"]:
                if "HA2" not in marker and "HA1" not in marker:
                    marker = adjust_position_and_get_h3_position(marker, hatype = hatype, H3_dict = None,
                                                                 protein = "H3")
                res.append(marker)
            updated_marker_dict["H3"] = res
        if protein in HA_TYPES:
            residues, unmapped_residues = process_ha_type(protein, marker_dict, structure_folder, hatype,
                                                          ha_type_info_dic)
            updated_marker_dict[protein] = unmapped_residues
            updated_marker_dict["H3"] = residues
            del updated_marker_dict[protein]
        elif protein in NA_TYPES:
            residues, unmapped_residues = process_na_type(protein, marker_dict, structure_folder, hatype,
                                                          ha_type_info_dic)
            updated_marker_dict[protein] = unmapped_residues
            updated_marker_dict["N2"] = residues

            del updated_marker_dict[protein]

    return transform_marker_dict(updated_marker_dict)


def read_fasta(file_path):

    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        description = record.description.split('_')[0]
        sequences[description] = str(record.seq)
    return sequences


def compare_sequences(seq_file1, seq_file2):

    seq_dict1 = read_fasta(seq_file1)
    seq_dict2 = read_fasta(seq_file2)

    length_differences = {}
    for key in seq_dict1:
        if key in seq_dict2 and key in [f"H{i}" for i in range(1, 19)]:
            length_differences[key] = abs(len(seq_dict1[key]) - len(seq_dict2[key]))
    return length_differences


def annotate_markers(markers_path, STRUCTURE_PATH, hatype = None):
    """
    Annotate markers by loading virulence data, then converting HA types.
    """

    marker_dict, data = load_markers(markers_path)

    markers_info = data.loc[:, "Protein"] + "-" + data.loc[:, "Amino acid site"]
    if "resistance" in markers_path or "binding" in markers_path:
        ha_type_info_dic = {}
    else:
        ha_type_info_dic = dict(zip(markers_info, data["HA_Type"]))

    marker_dict = {i: list(set(j)) for i, j in marker_dict.items()}

    marker_dict = transform_marker_dict(marker_dict)

    return marker_dict, data, ha_type_info_dic


def renumber_sequence(best_alignment):
    """
    Renumber a protein sequence based on the best alignment result.
    """

    renumbered_positions = []
    count = 1
    for std_char, query_char in zip(best_alignment[0], best_alignment[1]):
        if std_char != '-':
            renumbered_positions.append(f"{count}{query_char}")
            count += 1

    return renumbered_positions


def merge_dictionaries(dict1, dict2):
    merged_dict = {}
    for key in dict1:
        if key == 'H3':

            merged_h3_values = [dict1[key], dict2[key]]
            merged_dict[key] = merged_h3_values
        else:

            merged_dict[key] = dict1[key]
    return merged_dict


def convert_matrix_to_dict(matrix):
    matrix_dict = {}
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if i >= j:
                matrix_dict[(matrix.alphabet[i], matrix.alphabet[j])] = value
                matrix_dict[(matrix.alphabet[j], matrix.alphabet[i])] = value
    return matrix_dict


def perform_alignment_and_renumber(standard_seq_path, query_seq):
    standard_seq = next(SeqIO.parse(standard_seq_path, 'fasta')).seq

    matrix = substitution_matrices.load("BLOSUM62")
    matrix = convert_matrix_to_dict(matrix)

    alignments = pairwise2.align.globalds(standard_seq, query_seq, matrix, -10.0, -0.5)

    best_alignment = max(alignments, key = lambda x: x.score)

    return renumber_sequence(best_alignment)


def process_ha_na(protein_abbr, sequence):
    ha_results = defaultdict(dict)

    if protein_abbr in [f"H{i}" for i in range(1, 19)]:
        ha1_path = f"{STANDARD_PATH}/HA1/{protein_abbr}.fas"
        ha_results["HA1"][protein_abbr] = perform_alignment_and_renumber(ha1_path, sequence)

        ha2_path = f"{STANDARD_PATH}/HA2/{protein_abbr}.fas"
        ha_results["HA2"][protein_abbr] = perform_alignment_and_renumber(ha2_path, sequence)
    elif protein_abbr in [f"N{i}" for i in range(1, 10)]:
        na_path = f"{STANDARD_PATH}/{protein_abbr}.fas"
        ha_results[protein_abbr] = perform_alignment_and_renumber(na_path, sequence)
    return ha_results


def check_fasta_for_non_standard_amino_acids(input_file):

    standard_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    standard_sequences = []

    for record in SeqIO.parse(input_file, "fasta"):
        seq = str(record.seq).upper()
        seq_set = set(seq)

        if seq_set.issubset(standard_amino_acids):

            standard_sequences.append(record)
        else:

            signal_peptide_index = seq.find('*')
            if signal_peptide_index != -1:
                cleaned_seq = seq[:signal_peptide_index]
                cleaned_seq_set = set(cleaned_seq)
                if cleaned_seq_set.issubset(standard_amino_acids):
                    record.seq = Seq(cleaned_seq)
                    standard_sequences.append(record)
                else:
                    print(f"Non-standard amino acids found in sequence {record.id} after removing signal peptide.")
            else:
                print(f"Non-standard amino acids found in sequence {record.id}.")

    return standard_sequences


def renumber_proteins(fasta_path, acc_pro_dict, marker_dict):
    fasta_sequences = check_fasta_for_non_standard_amino_acids(fasta_path)
    renumbering_results = {}

    for record in fasta_sequences:
        protein_id = record.id
        protein_abbr = acc_pro_dict.get(protein_id)
        protein_abbr = "NS2" if protein_abbr == "NEP" else protein_abbr
        is_hana_type = protein_abbr in HA_TYPES or protein_abbr in NA_TYPES

        if protein_abbr in marker_dict or is_hana_type:

            if protein_abbr in [f"H{i}" for i in range(1, 19)]:

                ha_results = process_ha_na(protein_abbr, record.seq)

                renumbered_positions_HA1 = convert_HA_residues(ha_results["HA1"], STRUCTURE_PATH, hatype = "HA1",
                                                               ha_type_info_dic = None)
                renumbered_positions_HA2 = convert_HA_residues(ha_results["HA2"], STRUCTURE_PATH, hatype = "HA2",
                                                               ha_type_info_dic = None)
                renumbered_positions = merge_dictionaries(renumbered_positions_HA1, renumbered_positions_HA2)

                renumbered_positions[protein_id] = renumbered_positions.pop("H3")
                renumbering_results.update(renumbered_positions)
            elif protein_abbr in [f"N{i}" for i in range(1, 10)]:
                ha_results = process_ha_na(protein_abbr, record.seq)
                renumbered_positions = convert_HA_residues(ha_results, STRUCTURE_PATH, hatype = None)
                renumbered_positions[protein_id] = renumbered_positions.pop("N2")
                renumbering_results.update(renumbered_positions)
            else:

                standard_seq_path = os.path.join(STANDARD_PATH, f"{protein_abbr}.fas")
                renumbering_results[protein_id] = perform_alignment_and_renumber(standard_seq_path, record.seq)


        else:
            if pd.isna(protein_abbr):
                protein_abbr = "NA"
            print(f"No markers found for {protein_abbr} in the source data.")
    return renumbering_results


def merge_dicts_with_list(dict_list):
    """
    Function to merge a list of dictionaries. If keys are repeated, values are merged into a list.
    """

    merged_dict = {}
    for d in dict_list:
        for key, value in d.items():
            if key in merged_dict:

                if not isinstance(merged_dict[key], list):
                    merged_dict[key] = [merged_dict[key]]
                merged_dict[key].append(value)
            else:

                merged_dict[key] = value
    return merged_dict


def unique_dicts(dict_list):
    seen = set()
    unique_list = []
    for d in dict_list:

        tuple_repr = tuple(sorted(d.items()))

        if tuple_repr not in seen:
            seen.add(tuple_repr)
            unique_list.append(d)
    return unique_list


def generate_combinations(group):

    spec_type_groups = group.groupby('Specific Type')

    mutation_dict = defaultdict(list)
    for spec_type, g in spec_type_groups:
        for _, row in g.iterrows():
            if type(row["Amino acid site"]) == str:
                mutation_dict[spec_type].append({row['Protein']: row['Amino acid site'].strip()})
            else:
                mutation_dict[spec_type].append({row['Protein']: row['Amino acid site']})

    values_lists = [mutation_dict[key] for key in mutation_dict]

    values_lists = [unique_dicts(sublist) for sublist in values_lists]

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


def is_subset_complex_revised(dict1, dict2):
    """
    Check if one dictionary is a complex subset of another, with revised logic for nested dictionaries.
    """
    for key, value1 in dict1.items():
        if key not in dict2:
            return False

        value2 = dict2[key]

        if isinstance(value1, dict) and isinstance(value2, dict):
            if not is_subset_complex_revised(value1, value2):
                return False

        elif isinstance(value1, list) and isinstance(value2, list):
            if not set(value1).issubset(set(value2)):
                return False
        elif isinstance(value1, str) and isinstance(value2, list):
            if value1 not in value2:
                return False
        elif isinstance(value1, str) and isinstance(value2, str):
            if value1 != value2:
                return False

    return True


def format_marker(marker, protein_prefix = ''):
    if '-' in marker and "HA1" not in marker and "HA2" not in marker:
        amino_acid = marker.split('-')[0]

        deletion_suffix = "-"
    else:
        amino_acid = marker
        deletion_suffix = ""

    formatted_marker = f"{protein_prefix}-{amino_acid}{deletion_suffix}" \
        if protein_prefix else f"{amino_acid}{deletion_suffix}"
    return formatted_marker


def format_marker_list(markers, protein_prefix = ''):
    if isinstance(markers, str):
        return format_marker(markers, protein_prefix)

    not_ha1_ha2 = all('HA' not in marker for marker in markers)
    all_contain_dash = all('-' in marker for marker in markers)
    if all_contain_dash and not_ha1_ha2:
        numbers = [int(marker.split('-')[0]) for marker in markers]

        max_number = max(numbers)
        min_number = min(numbers)
        return f"{protein_prefix}-{min_number}-{max_number}CompleteDeletion"

    return '&'.join(format_marker(marker, protein_prefix) for marker in markers)


def process_dictionary(data_dict):
    formatted_list = []
    for protein, markers in data_dict.items():

        if isinstance(markers, dict):
            for sub_protein, sub_markers in markers.items():
                formatted_marker = format_marker_list(sub_markers, f"{protein}-{sub_protein}")
                formatted_list.append(formatted_marker)
        else:
            formatted_marker = format_marker_list(markers, protein)
            formatted_list.append(formatted_marker)

    return '&'.join(formatted_list)


def process_protein_sequence(acc_id, renumbered_position, acc_pro_dic, marker_markers):
    protein_type = acc_pro_dic[acc_id]

    if protein_type == "Unknown":
        return None, None
    HA_TYPES_ALL = [f"H{i}" for i in range(1, 19)]
    NA_TYPES_ALL = [f"N{i}" for i in range(1, 10)]

    use_protein = "H3" if protein_type in HA_TYPES_ALL else ("N2" if protein_type in NA_TYPES_ALL else protein_type)

    expected_markers = marker_markers.get(use_protein, {})

    protein = f'H3' if protein_type in HA_TYPES_ALL else (
        f'N2' if protein_type in NA_TYPES_ALL else protein_type)
    markers = defaultdict(list)
    if use_protein == "H3":
        for hatype, ha_markers in expected_markers.items():
            for marker in ha_markers:
                match = re.match(r"(\d+)([A-Z]|-)", marker)
                index = 0 if hatype == "HA1" else 1

                if match and match.group() in renumbered_position[index][hatype]:
                    markers[hatype].append(match.group())
        return protein, markers
    markers_oth = []
    for marker in expected_markers:
        match = re.match(r"(\d+)([A-Z]|-)", marker)
        if match and match.group() in renumbered_position:
            markers_oth.append(match.group())
    return protein, markers_oth


def check_marker_combinations(total_markers, results_markers, markers_type, input_file_name, data, ha_type, na_type):
    results = []

    for marker_protein_type, marker_list in total_markers.items():

        for proba_comb in marker_list:

            if proba_comb and is_subset_complex_revised(proba_comb, results_markers):

                if proba_comb and all(proba_comb.values()):
                    markers_formated = process_dictionary(proba_comb)
                    results.append({
                        'Strain ID': input_file_name.split(".")[0],
                        'Amino acid site': markers_formated,
                        'Protein Type': marker_protein_type,
                    })

    results = pd.DataFrame(results)
    final_results = pd.DataFrame()
    if not results.empty:
        final_results = merge_dataframes(results, data, markers_type, ha_type, na_type)
    return final_results


def split_amino_acid_site(site):
    """
    Split and return the part of the amino acid site identifier that is of interest,
    handling special cases.
    """
    site = site.strip()
    if "PB1-F2" in site or "PA-X" in site:
        return site.rsplit("-", 1)[-1]
    elif re.search(r"HA[12]-HA[12]", site):
        parts = re.split(r"HA[12]-", site, maxsplit = 1)
        return parts[1] if len(parts) > 1 else site
    else:
        return site.split("-", 1)[1] if "-" in site else site


def process_results(df):
    """
    Process a DataFrame by modifying the 'Protein Type' column based on conditions derived from the 'Amino acid site'.
    """
    df.loc[:, "Ori Protein Type"] = df.loc[:, "Protein Type"]
    df.loc[:, "Protein Type"] = df.loc[:, "Amino acid site"].apply(
        lambda x: x.split('-', 1)[0] if "PA-X" not in x and "PB1-F2" not in x else x.rsplit('-', 1)[0])
    return df


def merge_dataframes(results, data, markers_type, ha_type, na_type):
    combination_pattern = re.compile(r'combination')

    data['HasCombination'] = data['Protein Type'].apply(lambda x: bool(combination_pattern.search(x)))
    data_with_combination = data[data['HasCombination']].drop_duplicates(subset = 'Protein Type')
    data_without_combination = data[~data['HasCombination']]
    data_without_combination.loc[:, "Protein Type"] = data_without_combination.loc[:, "Protein"]

    results['HasCombination'] = results['Protein Type'].apply(lambda x: bool(combination_pattern.search(x)))

    results_with_combination = results[results['HasCombination']]

    results_with_combination.loc[:, "Ori Protein Type"] = results_with_combination.loc[:, "Protein Type"]
    results_without_combination = process_results(results[~results['HasCombination']])

    final_results = pd.DataFrame()

    if not results_with_combination.empty:
        data_with_combination.drop(columns = ["Amino acid site", "HasCombination"], inplace = True)

        merged_with_combination = pd.merge(results_with_combination, data_with_combination,
                                           on = 'Protein Type', how = 'left')
        final_results = pd.concat([final_results, merged_with_combination], ignore_index = True)

    if not results_without_combination.empty:
        df_copy = results_without_combination.copy()
        df_copy["Amino acid site"] = df_copy["Amino acid site"].apply(split_amino_acid_site)
        merged_without_combination = pd.merge(df_copy, data_without_combination,
                                              on = ['Protein Type', 'Amino acid site'], how = 'left')
        final_results = pd.concat([final_results, merged_without_combination], ignore_index = True)

    if not final_results.empty:
        final_results.rename(columns = {'Amino acid site': f'{markers_type.title()} Markers'}, inplace = True)

        final_results.loc[:, "Protein Type"] = final_results.loc[:, "Ori Protein Type"].apply(
            lambda x: get_hana_string(x, ha_type, na_type))
        final_results.loc[:, "Protein Type"] = final_results.loc[:, "Ori Protein Type"].apply(
            lambda x: get_hana_string(x, ha_type, na_type))

        columns = final_results.columns
        deleted_columns = [i for i in columns if "HasCombination" in i]
        deleted_columns.extend(["Specific Type", "Protein", "Ori Protein Type"])
        final_results.drop(columns = deleted_columns, inplace = True)

        final_results.replace('', np.nan, inplace = True)
        final_results.dropna(subset = ['Strain ID', f'{markers_type.title()} Markers', 'Protein Type'], how = "all",
                             inplace = True)

        subset = [f'{markers_type.title()} Markers', 'Protein Type']

        if markers_type == "resistance":
            subset.extend(['Drug', 'Resistance_level'])

            cols_to_change = ['oseltamivir', 'zanamivir', 'peramivir', 'laninamivir', 'baloxavir', 'adamantane']
            final_results['Drug'] = pd.Categorical(final_results['Drug'], categories = cols_to_change, ordered = True)

            final_results.sort_values(by = 'Drug', inplace = True)

        final_results.drop_duplicates(subset = subset, keep = "first", inplace = True)
        final_results.loc[:, "Protein Type"] = \
            final_results.loc[:, "Protein Type"].str.split("-combination", expand = True)[0]
    return final_results


def get_hana_string(protein_type, ha_type, na_type):
    if any(ha in protein_type for ha in HA_TYPES) or 'H3' in protein_type:

        return 'H3'
    elif any(na in protein_type for na in NA_TYPES) or 'N2' in protein_type:

        return 'N2'
    else:
        return protein_type


def process_protein_sites(df, protein_col, site_col):
    protein_order = ["PB2", "PB1", "PB1-F2", "PA", "PA-X"] + [f"H{i}" for i in range(1, 19)] + ["NP"] + \
                    [f"N{i}" for i in range(1, 11)] + ["NS1", "NS2", "M1", "M2"]

    processed_rows = []
    for index, row in df.iterrows():
        if row[protein_col] == 'combination':

            parts = row[site_col].split('&')

            split_parts = [re.match(r'([A-Za-z0-9-]+)-(\d+)([A-Z-])', part).groups() for part in parts]

            def extract_protein_main(protein):

                if 'HA1' in protein or 'HA2' in protein:
                    return protein.split('-')[0]
                return protein

            sorted_parts = sorted(split_parts, key = lambda x: (
                protein_order.index(extract_protein_main(x[0])) if extract_protein_main(x[0]) in protein_order else len(
                    protein_order),
                int(x[1])
            ))

            new_site = '&'.join([f"{part[0]}-{part[1]}{part[2]}" for part in sorted_parts])
            processed_row = row.to_dict()
            processed_row[site_col] = new_site
            processed_rows.append(processed_row)
        else:
            if ('HA1' in row[site_col] and '&' in row[site_col]) or ('HA2' in row[site_col] and '&' in row[site_col]):

                pattern = r'(\w+-HA[1-2]-\d+[A-Z])'
                matches = re.findall(pattern, row[site_col])

                ha_parts = [match for match in matches if 'HA' in match]

                sorted_ha_parts = sorted(ha_parts, key = lambda x: int(re.findall(r'\d+', x.split('-')[-1])[0]))

                new_site = '&'.join(sorted_ha_parts)
                processed_row = row.to_dict()
                processed_row[site_col] = new_site
                processed_rows.append(processed_row)
            elif "eletion" in row[site_col]:

                processed_rows.append(row.to_dict())
            else:

                site = row[site_col]
                parts = site.split('&')
                numbers = [re.findall(r'\d+', part)[0] for part in parts]
                letters = list(parts[-1][len(numbers[-1]):])
                if any(protein in row[site_col] for protein in protein_order) and '&' in row[site_col]:
                    parts = [i.split("-")[-1] for i in parts]
                    numbers = [re.findall(r'\d+', part)[0] for part in parts]
                    letters = [re.findall(r'[A-Z]', part)[0] for part in parts]

                if len(numbers) == len(letters):
                    combined = list(zip(numbers, letters))

                    sorted_combined = sorted(combined, key = lambda x: int(x[0]))

                    sorted_numbers = [num for num, char in sorted_combined]
                    sorted_letters = [char for num, char in sorted_combined]
                    new_site = '&'.join(sorted_numbers) + ''.join(sorted_letters)
                    processed_row = row.to_dict()
                    processed_row[site_col] = new_site
                    processed_rows.append(processed_row)
                else:

                    processed_rows.append(row.to_dict())

    processed_df = pd.DataFrame(processed_rows)
    return processed_df


def should_continue(pro, dic):
    """
    Determines whether to continue processing based on the presence of "combination" in the
    property string and specific types in a dictionary against expected keys.
    """
    if "combination" not in pro:
        return False

    if any(type in pro for type in HA_TYPES):
        return set(dic.keys()) != {"H3"}

    if any(type in pro for type in NA_TYPES):
        return set(dic.keys()) != {"N2"}

    return False


def process_row(row):
    """
    Processes a row of data to identify virulence markers in protein sequences
    based on the provided original site and original protein type.
    """
    map_dic = {}

    for phenotype_ in sub_dic:
        map_dic.update(sub_dic[phenotype_])

    modified_dict = {key.translate({ord(c): None for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'}): value.translate(
        {ord(c): None for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'}) for key, value in map_dic.items()}

    original_site = row['Original Site']
    original_protein_type = row['Original protein type']

    if original_protein_type == 'combination':

        row['params'] = original_site
    else:
        special_proteins = [f'H{i}' for i in range(1, 19)]
        if original_protein_type in special_proteins:
            patterns = []
            if 'HA' in original_site:

                parts = original_site.split('&')
                for part in parts:
                    ha_patterns = re.search(r'HA[1-2]-\d{1,3}', part)
                    if ha_patterns:
                        patterns.append(ha_patterns.group())
            else:

                markers = re.findall(r'\d+', original_site)
                for marker in markers:
                    if marker in modified_dict:
                        patterns.append(f"HA2-{modified_dict[marker]}")
                    else:
                        patterns.append(f"HA1-{marker}")
            row['params'] = '|'.join(patterns)
        else:

            patterns = re.findall(r'\d+', original_site)
            row['params'] = '|'.join(patterns)

    row['pdb'] = original_protein_type

    return row


def generate_columns_dict(marker_type):
    common_columns = ['Strain ID', 'Protein Type', 'Original Site', 'Original protein type', 'PMID', 'Source',
                      'HA_Type', 'params', 'pdb']

    if marker_type == 'adaptation':
        specific_columns = ['Adaptation Markers', 'Avian residues']
    elif marker_type == 'binding':
        specific_columns = ['Binding Markers']
    elif marker_type == 'resistance':
        specific_columns = ['Resistance Markers', 'Drug', 'Resistance_level', 'Protein Type_y']
    elif marker_type == 'virulence':
        specific_columns = ['Virulence Markers', 'Phenotypic Consequences']
    elif marker_type == 'transmissibility':
        specific_columns = ['Transmissibility Markers']
    else:
        raise ValueError(f"Unknown marker_type: {marker_type}")

    columns = common_columns + specific_columns
    return columns


def sort_results_by_protein_type(results_df):
    """
    Sort the given DataFrame by the 'Protein Type' column according to a specific order.
    """

    protein_order = [
        "PB2", "PB1", "PB1-F2", "PA", "PA-X", "H3", "NP", "N2", "NS1", "NS2", "M1", "M2", "combination"
    ]

    results_df['Protein Type'] = pd.Categorical(results_df['Protein Type'], categories = protein_order, ordered = True)

    results_df = results_df.sort_values(by = 'Protein Type')

    return results_df


def identify_markers(input_file_path, renumbering_results, marker_markers, acc_pro_dic, markers_type, data,
                     output_directory = ".", prefix = "", ha_type_info_dic = None):
    """
    Identifies virulence markers in protein sequences based on the provided marker markers
    and the renumbered sequences.
    """

    os.makedirs(output_directory, exist_ok = True)
    input_file_name = os.path.split(input_file_path)[1]
    results_markers = defaultdict(list)

    for acc_id, renumbered_position in renumbering_results.items():

        protein, markers = process_protein_sequence(acc_id, renumbered_position, acc_pro_dic, marker_markers)
        if protein:
            results_markers[protein] = markers

    ha_type = na_type = None

    for acc, pro in acc_pro_dic.items():
        if pro in HA_TYPES or pro == "H3":
            ha_type = pro
        elif pro in NA_TYPES or pro == "N2":
            na_type = pro

    ori_markers = generate_protein_dict(load_total_markers(data))
    total_markers = defaultdict(list)
    for pro, lst in ori_markers.items():
        for dic in lst:
            if dic and all(dic.values()):

                if should_continue(pro, dic):
                    continue

                total_markers[pro].append(
                    convert_HA_residues(dic, STRUCTURE_PATH, hatype = None, ha_type_info_dic = ha_type_info_dic))

    results_df = check_marker_combinations(total_markers, results_markers, markers_type,
                                           input_file_name, data, ha_type, na_type)

    add_prefix = prefix + "_" if prefix else ""
    filename = add_prefix + input_file_name.split(".")[0] + "_markers.csv"

    if results_df.empty:

        results_df = pd.DataFrame(columns = generate_columns_dict(markers_type.lower()))
        message = "No markers found. Saving an empty file."

        results_df.to_csv(f"{output_directory}/{filename}", index = False, encoding = 'utf-8-sig')

        print(message)
    else:
        if "Phenotypic Consequences" in results_df.columns:
            results_df.loc[results_df['Source'].isnull(), 'Source'] = results_df.loc[
                results_df['Source'].isnull(), 'Phenotypic Consequences']

        results_df = process_protein_sites(results_df, "Protein Type", f"{markers_type.title()} Markers")
        results_df = results_df.apply(process_row, axis = 1)
        results_df = sort_results_by_protein_type(results_df)

        results_df.to_csv(f"{output_directory}/{filename}", index = False, encoding = 'utf-8-sig')

    return results_df


def is_fasta_file(filename):
    try:
        with open(filename, "r") as handle:
            first_record = next(SeqIO.parse(handle, "fasta"))

            return True
    except Exception as e:

        print(f"The file {filename} is not in FASTA format or could not be correctly parsed.")

        return False


def extract_element_fast(query, cds_list):
    """
    Extract the corresponding element from the list based on the query.
    This function first converts the list to a dictionary for faster lookup.

    :param query: The query string (e.g., 'querySeq1').
    :param cds_list: The list of tuples to search through.
    :return: The corresponding element if found, or None if not found.
    """
    cds_dict = {key: value for key, value in cds_list}
    return cds_dict.get(query, None)


def extract_protein_annotations(protein_path, result_path, CDSList):
    base_name = os.path.basename(protein_path).split('.')[0]

    result_dir_name = os.path.dirname(result_path)

    # output_path = os.path.join(result_dir_name, "{}_annotation.csv".format(base_name))
    output_path_for_anno = os.path.join(result_dir_name, "{}_annotation_for_anno.csv".format(base_name))

    sequences = SeqIO.parse(protein_path, "fasta")

    out_file_for_anno = open(output_path_for_anno, 'w')
    # with open(output_path, 'w') as out_file:
    written_instances = set()
    if CDSList:
        # out_file.write("qseqid,oriseqid,newqseqid,GeneType,CDS region,ProteinType,\n")
        out_file_for_anno.write("qseqid,oriseqid,newqseqid,GeneType,CDS region,ProteinType,\n")
    else:
        # out_file.write("qseqid,oriseqid,newqseqid,Protein Abbreviation,CDS region,ProteinType\n")
        out_file_for_anno.write("qseqid,oriseqid,newqseqid,Protein Abbreviation,CDS region,ProteinType\n")
    for record in sequences:
        seq_id = record.id.split("|")[0]
        oriseqid = record.id.split(seq_id)[1].split("|", 1)[1]

        qseqid = re.search("querySeq\d+", seq_id).group()
        protein_nucl = re.search("(querySeq\d+)_(.*)_(.*)", seq_id)

        instance_type = seq_id.rsplit("_", 1)[1]
        proteintype = protein_nucl.group(2) if CDSList else instance_type
        if proteintype == "HA" or proteintype == "NA":
            proteintype = protein_nucl.group(3)
        cds_region = extract_element_fast(qseqid, CDSList) if CDSList else None

        out_file_for_anno.write(
            '"{0}","{1}","{2}","{3}","{4}","{5}"\n'.format(qseqid, oriseqid, record.id, instance_type, cds_region,
                                                           proteintype))

        if CDSList:
            if instance_type in written_instances:
                continue
            written_instances.add(instance_type)

        # out_file.write(
        #     '"{0}","{1}","{2}","{3}","{4}","{5}"\n'.format(qseqid, oriseqid, record.id, instance_type, cds_region,
        #                                                    proteintype))

    print("Annotation written to {}".format(output_path_for_anno))


def parse_args():
    parser = argparse.ArgumentParser(
        prog = 'fluhp',
        description = 'Influenza sequence analysis and host prediction tool'
    )

    # Add a more descriptive help message for the subcommands
    subparsers = parser.add_subparsers(
        dest = 'subcommand',
        help = 'Available analysis modules',
        # required = True  # Enforce selection of a subcommand
    )

    # anno subcommand
    anno_parser = subparsers.add_parser('anno',
                                        help = 'Annotate a FASTA file or all FASTA files in a directory '
                                               'using DIAMOND BLAST against a flu database')
    anno_parser.add_argument('-i', '--input', required = True,
                             help = 'Input FASTA file or directory containing FASTA files')
    anno_parser.add_argument('-o', '--output_directory', type = str, default = "result/",
                             help = 'Directory to save the output files. Defaults to the result directory')
    anno_parser.add_argument('-u', '--updated_directory', type = str, default = 'standardized_fasta/',
                             help = 'Directory to save the standardize fasta files. Defaults to the standardized_fasta directory.')

    # extract subcommand
    extract_parser = subparsers.add_parser('extract', help = 'Extract and process amino acid markers.')
    extract_parser.add_argument('-i', '--input', required = True,
                                help = 'Input FASTA file or directory containing standardized FASTA files')
    extract_parser.add_argument('-a', '--anno_path', required = True,
                                help = 'Input annotation CSV file or directory containing annotation CSV files')
    extract_parser.add_argument('-o', '--output_directory', type = str, default = 'adaptation/',
                                help = 'Directory to save the output files. Defaults to the current directory')
    extract_parser.add_argument('-p', '--prefix', type = str, default = '', help = 'Prefix for the output filenames')

    # Prediction Subcommands (with improved structure)
    def add_prediction_parser(subparsers, parser_name, name, description):
        pred_parser = subparsers.add_parser(
            parser_name,
            help = f'Predict {description} using machine learning models',
            description = f'Advanced predictive modeling for {description} analysis'
        )
        pred_parser.add_argument(
            '-i', '--input',
            required = True,
            help = 'Input CSV file with marker data or directory of marker files'
        )
        pred_parser.add_argument(
            '-t', '--threshold',
            type = float,
            default = 0.5,
            help = 'Probability threshold for prediction (default: 0.5)'
        )
        pred_parser.add_argument(
            '-o', '--output_directory',
            default = f'{name}_predictions/',
            help = f'Output directory for prediction results (default: {name}_predictions/)'
        )
        pred_parser.add_argument(
            '-p', '--prefix',
            default = '',
            help = 'Optional prefix for output prediction files'
        )
        return pred_parser

    # Create specific prediction parsers
    add_prediction_parser(subparsers, 'predh', 'host', 'host adaptation')

    return parser.parse_args()


def process_anno_cmd(input_file, args):
    """
    Call the appropriate functions to process a single fasta file
    """

    os.makedirs(args.output_directory, exist_ok = True)
    proteinPath, resultPath, CDSList, predictedProteinType = ivew_task(args.output_directory, TEMP_PATH,
                                                                       str(input_file),
                                                                       args.updated_directory)

    extract_protein_annotations(proteinPath, resultPath, CDSList)


class FileNotFoundInDirectoryException(Exception):
    """
    Exception raised when no matching file is found in the directory.
    """

    def __init__(self, directory, search_string):
        self.directory = directory
        self.search_string = search_string
        message = f"No file containing '{search_string}' was found in '{directory}'."
        super().__init__(message)


def find_files_with_string(directory, string):
    """
    Find files in a directory that contain the specified string in their name.

    :param directory: The directory to search in.
    :param string: The string to search for in file names.
    :return: The name of the first file found that matches the criteria.
    :raises FileNotFoundInDirectoryException: If no matching file is found.
    """
    all_items = os.listdir(directory)
    string = os.path.splitext(string)[0]
    files_with_string = [item for item in all_items

                         if
                         string == item.rsplit("_annotation_for_anno.csv")[0] and os.path.isfile(
                             os.path.join(directory, item))
                         and item.endswith("_annotation_for_anno.csv")]
    if not files_with_string:
        raise FileNotFoundInDirectoryException(directory, string)
    return files_with_string[0]


def process_extract_cmd(input_file, args, is_directory = True):
    """
    Process a single FASTA file by finding and reading the corresponding annotation file.

    :param input_file: The input FASTA file.
    :param args: Command line arguments.
    :param is_directory: Flag indicating if the input is a directory (default True).
    """
    print(f"Input File:\n{input_file}")

    if ".trans2protein.fasta" in str(input_file):
        input_filename_pre = os.path.split(input_file)[1].split('.trans2protein.fasta')[0]
    else:
        input_filename_pre = os.path.split(input_file)[1].split('.stdName.annotation')[0]
    try:
        if is_directory:
            anno_filename = find_files_with_string(args.anno_path, input_filename_pre)
            annotations = pd.read_csv(os.path.join(args.anno_path, anno_filename))
        else:
            annotations = pd.read_csv(args.anno_path)
    except FileNotFoundError as e:
        print(f"Annotation file not found for {input_file}: {e}")
        return
    except Exception as e:
        print(f"Error processing extraction command for {input_file}: {e}")
        return
    # print(annotations)
    annotations["ProteinType"] = annotations["ProteinType"].apply(lambda x: "NS2" if x.strip() == "NEP" else x)
    acc_pro_dic = dict(zip(annotations.loc[:, "newqseqid"], annotations.loc[:, "ProteinType"]))
    for filename in os.listdir(MARKER_PATH):

        if filename.endswith("_formated.csv") and "human" in filename:
            marker_dict, data, ha_type_info_dic = annotate_markers(os.path.join(MARKER_PATH, f"{filename}"),
                                                                   STRUCTURE_PATH)

            renumbering_results = renumber_proteins(
                fasta_path = str(input_file),
                acc_pro_dict = acc_pro_dic,
                marker_dict = marker_dict,
            )
            markers_type = filename.split("_formated.csv")[0]
            markers_type = markers_type.split('_')[1] if "_" in markers_type else markers_type
            out_dir = args.output_directory if args.output_directory else markers_type
            results_df = identify_markers(
                input_file_path = str(input_file),
                renumbering_results = renumbering_results,
                marker_markers = marker_dict,
                acc_pro_dic = acc_pro_dic,
                # output_directory = markers_type,
                output_directory = out_dir,
                prefix = args.prefix,
                markers_type = markers_type,
                data = data,
                ha_type_info_dic = ha_type_info_dic,
            )
            print("Marker extracted and saved to file.\n")


def is_extractable(file_name):
    return file_name.endswith((".stdName.annotation.fa", ".stdName.annotation.fasta", ".stdName.annotation.fas"))


def process_directory(directory_path, args, max_files = 30000):
    """
    Process files in the specified directory, renaming files that are valid fasta files
    and executing different processing commands based on the subcommand.

    :param directory_path: Path to the directory to be processed (Path object).
    :param args: Command line arguments object.
    :param max_files: Maximum number of files to process.
    """
    directory = Path(directory_path)
    if not directory.is_dir():
        print(f"Error: {directory} is not a directory.")
        return
    file_count = 0
    fasta_count = 0
    files = os.listdir(directory_path)
    files.sort(key = lambda x: [int(s) if s.isdigit() else s.lower() for s in re.split(r'(\d+)', x)])

    for file_name in files:
        if file_count >= max_files:
            print(f"Reached the maximum of {max_files} files to process.")
            break

        file = directory / file_name

        if args.subcommand == 'anno':
            fasta_count += 1

            if is_fasta_file(str(file)):
                process_anno_cmd(file, args)

        elif args.subcommand == 'extract' and is_extractable(str(file)):
            process_extract_cmd(file, args)

        file_count += 1


def process_single_file(file, args):
    file_path = Path(file)

    if not is_fasta_file(file_path):
        print(f"The file {file_path} is not in FASTA format or could not be correctly parsed.")
        return

    if args.subcommand == 'anno':
        if not re.match(r'isolate\d+\.fasta$', file.name):
            new_name = file_path.parent / 'isolate1.fasta'
            file_path.rename(new_name)
            file_path = new_name

        process_anno_cmd(file_path, args)
    elif args.subcommand == 'extract' and str(args.anno_path).endswith("_annotation_for_anno.csv"):
        process_extract_cmd(file_path, args, is_directory = False)


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
    if args.subcommand == "predh":
        predictions = predh(
            str(Path(args.input)),
            f'{MODEL_PATH}/host_ensemble_model.joblib',
            args.threshold,
            args.output_directory,
            args.prefix
        )
        print(f"Predictions completed.")
    else:
        run_other_subcommand(args)


if __name__ == '__main__':
    length_diffs = compare_sequences(STD_PATH, COMPLETE_STD_PATH)
    main()
