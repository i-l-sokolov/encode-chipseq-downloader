# ENCODE ChIP-seq Downloader

A simple tool to download transcription factor ChIP-seq data from the ENCODE Project database.

## Features

- Query ENCODE for ChIP-seq experiments by TF and cell line
- Filter for specific genome assemblies (GRCh38)
- Exclude treated samples
- Generate organized directory structure
- Create download scripts for batch downloading

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### 1. Configure Your Query

Edit `tables/tfs.txt`:
```
TF
CEBPB
ATF4
MYBL2
```

Edit `tables/cells.txt`:
```
Cell_line
A549
HepG2
K562
```

### 2. Run the Downloader

```bash
python code/downloader.py
```

This will:
- Query ENCODE for all TF/cell line combinations
- Download metadata as JSON
- Parse and filter experiments
- Generate download scripts

### 3. Execute Downloads

```bash
chmod +x cell_lines/*/download_files
for file in cell_lines/*/download_files; do ./$file; done
```

## Output Structure

```
cell_lines/
├── A549/
│   ├── json/              # ENCODE metadata
│   ├── files/             # Downloaded BED files
│   ├── files_cut/         # Organized BED files
│   ├── res_dict.json      # Parsed metadata
│   ├── res_dataframe.csv  # Experiment summary
│   └── download_files     # Download script
└── HepG2/
    └── ...
```

## Requirements

- Python 3.8+
- Internet connection
- ~50-100 GB storage (depending on number of cell lines)