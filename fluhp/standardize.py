# !/usr/bin/python
# -*- coding: utf-8 -*-
# Author:lihuiru

from Bio import SeqIO
import os


def is_protein_sequence(seq):
    """
    Determines if the given sequence is a protein sequence.

    Parameters:
        seq (str): The sequence to be analyzed.

    Returns:
        bool: True if the sequence is a protein sequence, False if it is a nucleotide sequence.
    """
    # Convert the sequence to uppercase to ensure case insensitivity
    seq = seq.upper()
    seq = [i for i in seq if i != 'N']
    # Calculate the length of the sequence
    seq_length = len(seq)

    # Calculate the proportion of AGCT characters in the sequence
    known_nucleotide_chars = "ATCG"
    agct_count = sum(seq.count(char) for char in known_nucleotide_chars)
    agct_proportion = agct_count / seq_length

    # If the proportion of AGCT characters is above a certain threshold, it's likely a nucleotide sequence
    # Here we assume a threshold of 85%, but this can be adjusted based on requirements
    threshold = 0.9
    if agct_proportion > threshold:
        return False  # It's a nucleotide sequence

    return True  # It's a protein sequence

def standardize(filePath, outFileDir):
    """
    Standardizes gene sequence files using BioPython and determines if the sequences are DNA or protein.

    Parameters:
        filePath (str): Path to the file to be processed.
        outFileDir (str): Directory where the output file should be stored.

    Returns:
        tuple: A tuple containing the output file name, sequence type ('nucleo' or 'protein'), and a dictionary of standardized sequences.
    """
    file = os.path.basename(filePath)
    outFileName = file + ".stdName"
    is_protein = False
    dic = {}
    n = 0
    # print(f"Reading file: {filePath}")
    # print(f"Output file: {os.path.join(outFileDir, outFileName)}")

    with open(filePath, 'r') as inFile, open(outFileDir + "/" + outFileName, 'w') as outFile:
        for record in SeqIO.parse(inFile, "fasta"):
            n += 1
            standardName = ">querySeq" + str(n)
            sequence = str(record.seq).upper()

            dic[standardName] = f">{record.id}"  # Store the sequence ID and original ID
            outFile.write(f"{standardName}\n{sequence}\n")
            # Check if the sequence is a protein sequence
            if not is_protein and is_protein_sequence(sequence):
                is_protein = True

    sequence_type = "protein" if is_protein else "nucleo"
    return outFileName, sequence_type, dic

if __name__ == '__main__':
    seq = 'GTGGCAAAAACATAATGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    test = is_protein_sequence(seq)
    print(test)
