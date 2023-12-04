# fluhp: Influenza Host Prediction Command Line Tool

`fluhp` is a command-line utility designed for extracting and annotating flu markers from protein sequence data, utilizing DIAMOND BLAST for efficient annotation of influenza A virus sequences. It can handle single files, where each file contains all proteins of an influenza virus strain, or process directories, where each file within represents all proteins of distinct influenza virus strains.

**Note:** The input sequences must be protein sequences. Inputs can either be a single file containing all protein sequences of a flu virus strain or a directory where each file corresponds to all protein sequences of separate flu virus strains.

## Quick Start

```bash
conda create -n fluhp-env python=3.6
conda activate fluhp-env
git clone https://github.com/lihuirull/FluHostPredictor.git
cd FluHostPredictor
pip install .
fluhp anno -i your_data.fasta -o output_directory
```

## Installation

`fluhp` can be installed using the following method:

### From Source

First, clone the repository and install Python dependencies:

```shell
git clone https://github.com/lihuirull/FluHostPredictor.git
cd FluHostPredictor
pip install .
```

This will also attempt to install the `diamond` dependency if it is available from the `bioconda` channel. If not available, you will need to install `diamond` manually.

Please note that installing with `pip` may not install non-Python dependencies such as `diamond`, which you will need to install separately.

## Usage

`fluhp` includes three subcommands: `anno`, `extract`, and `pred`.

### Annotate (`anno`)

Annotate a FASTA file or all FASTA files in a directory using DIAMOND BLAST against a flu database.

**Example:**

```shell
fluhp anno -i tests/test1.fasta -o path/to/output_dir
```

**Arguments:**

- `-i, --input`: Input FASTA file or directory containing FASTA files (required).
- `-o, --output_directory`: Directory to save the output files (default: current directory).
- `-p, --prefix`: Prefix for the output filenames (default: none).
- `-e, --evalue`: E-value threshold for DIAMOND BLAST hits (default: 1e-5).
- `-u, --update_file`: If set, updates the FASTA file with annotations (flag).
- `-t, --threads`: Number of threads for DIAMOND BLAST (default: 10).

### Extract (`extract`)

Extract and process protein annotations from annotated FASTA files.

**Example:**

```shell
fluhp extract -i path/to/input.fasta -a path/to/annotations.csv -o path/to/output_dir
```

**Arguments:**

- `-i, --input`: Input FASTA file or directory containing FASTA files (required).
- `-a, --anno_path`: Input annotation CSV file or directory containing annotation CSV files (required).
- `-r, --receptor_markers`: Include this flag to use the optional Excel file with receptor markers (flag).
- `-o, --output_directory`: Directory to save the output files (default: current directory).
- `-p, --prefix`: Prefix for the output filenames (default: none).

### Predict (`pred`)

Predict new data labels using a trained model.

**Example:**

```shell
fluhp pred -i path/to/marker_data.csv -o path/to/output_dir
```

**Arguments:**

- `-i, --input`: Input CSV file with marker data or directory containing such files (required).
- `-m, --model_path`: Path to the trained model file (default: included model).
- `-th, --threshold_path`: Path to the file containing the optimal threshold value (default: included threshold).
- `-f, --top_features_path`: Path to the file containing the top features (default: included features).
- `-o, --output_directory`: Directory to save the prediction results (default: current directory).
- `-p, --prefix`: Prefix for the output filenames of the predictions (default: none).

## Dependencies

Python dependencies are listed in the `requirements.txt` file and can be installed using `pip`. However, `fluhp` also requires the `diamond` tool, which is not a Python package and needs to be installed separately. Instructions for installing `diamond` can be found at its official [documentation](https://github.com/bbuchfink/diamond/wiki).

For a comprehensive setup, including all dependencies, please follow the installation instructions provided in the sections above.

## Output

The `fluhp` tool will generate the following outputs based on the subcommand executed:

- For `anno`: Annotated FASTA files or a directory of annotated FASTA files with added annotations based on DIAMOND BLAST results.
- For `extract`: A CSV file with extracted markers and annotations formatted as: `Strain ID,Adaptation Markers,Protein Type,Avian residues,Doi,Source,PMID`
- For `pred`: This parameter outputs not only the predicted class labels for the new marker data based on the trained model but also the prediction probabilities. The current study employs a very small probability threshold.

## Notes

Ensure that the DIAMOND database and the trained model are compatible with the sequences you wish to annotate or predict. It is essential to use the same versions of the software and model that were used during the initial training and annotation for consistency and accuracy.

For further assistance or to report issues, please visit the `fluhp` GitHub repository issues section.