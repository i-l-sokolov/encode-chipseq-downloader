"""
ENCODE ChIP-seq Data Downloader
================================
Downloads and organizes ChIP-seq BED files from ENCODE database.
"""

import os
import json
import requests
from glob import glob
from tqdm import tqdm
import pandas as pd
import shutil


def sanitize_name(name):
    """Remove special characters from cell line names."""
    chars_to_replace = [' ', '"', "'", ',', '-', '.', '/']
    for char in chars_to_replace:
        name = name.replace(char, '')
    return name


def get_search_url(cell_line, tf):
    """Construct ENCODE search URL."""
    return (
        f'https://www.encodeproject.org/search/'
        f'?type=Experiment'
        f'&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens'
        f'&assay_title=TF+ChIP-seq'
        f'&status=released'
        f'&target.label={tf}'
        f'&biosample_ontology.classification=cell+line'
        f'&biosample_ontology.term_name={cell_line}'
        f'&frame=embedded'
    )


def create_directories(cell_lines, output_dir='cell_lines'):
    """Create directory structure for all cell lines."""
    rename_dict = {cl: sanitize_name(cl) for cl in cell_lines}
    
    for cell_line in cell_lines:
        cell_folder = rename_dict[cell_line]
        base_path = f'{output_dir}/{cell_folder}'
        
        os.makedirs(f'{base_path}/json', exist_ok=True)
        os.makedirs(f'{base_path}/files', exist_ok=True)
        os.makedirs(f'{base_path}/files_cut', exist_ok=True)
    
    return rename_dict


def download_metadata(cell_lines, tfs, rename_dict, output_dir='cell_lines'):
    """Download experiment metadata from ENCODE."""
    headers = {'accept': 'application/json'}
    
    print("Downloading metadata from ENCODE...")
    for cell_line in tqdm(cell_lines, desc="Cell lines"):
        cell_folder = rename_dict[cell_line]
        
        for tf in tfs:
            url = get_search_url(cell_line, tf)
            response = requests.get(url, headers=headers)
            biosample = response.json()
            
            output_path = f'{output_dir}/{cell_folder}/json/{tf}.json'
            with open(output_path, 'w') as f:
                json.dump(biosample, f)


def parse_metadata(output_dir='cell_lines'):
    """Parse JSON metadata to extract experiment details."""
    print("\nParsing metadata...")
    
    for cell_folder in tqdm(sorted(glob(f'{output_dir}/*')), desc="Processing"):
        results = {}
        
        for json_file in glob(f'{cell_folder}/json/*.json'):
            tf = os.path.basename(json_file).replace('.json', '')
            results[tf] = {}
            
            with open(json_file, 'r') as f:
                biosample = json.load(f)
            
            for exp_idx, experiment in enumerate(biosample.get('@graph', [])):
                default_analysis = experiment.get('default_analysis')
                
                default_idx = None
                for analysis_idx, analysis in enumerate(experiment.get('analyses', [])):
                    if analysis.get('@id') == default_analysis:
                        default_idx = analysis_idx
                        break
                
                if default_idx is not None:
                    results[tf][exp_idx] = {
                        'accession_code': experiment.get('@id'),
                        'n_default': default_idx,
                        'default_analysis': default_analysis,
                        'assembly': experiment['analyses'][default_idx].get('assembly'),
                        'lab': experiment.get('lab', {}).get('@id'),
                        'description': experiment.get('description', ''),
                        'treatment': 'treat' in experiment.get('description', '')
                    }
        
        with open(f'{cell_folder}/res_dict.json', 'w') as f:
            json.dump(results, f)


def create_dataframes(output_dir='cell_lines'):
    """Create summary dataframes from parsed metadata."""
    print("\nCreating summary dataframes...")
    
    for res_dict_path in tqdm(sorted(glob(f'{output_dir}/*/res_dict.json'))):
        with open(res_dict_path, 'r') as f:
            res_dict = json.load(f)
        
        records = []
        for tf, experiments in res_dict.items():
            for exp_data in experiments.values():
                if (exp_data['assembly'] == 'GRCh38' and not exp_data['treatment']):
                    records.append({
                        'prot': tf,
                        'lab': exp_data['lab'],
                        'desc': exp_data['description'],
                        'n_def': exp_data['n_default'],
                        'def_anal': exp_data['default_analysis'],
                        'access_code': exp_data['accession_code']
                    })
        
        df = pd.DataFrame(records)
        output_path = res_dict_path.replace('res_dict.json', 'res_dataframe.csv')
        df.to_csv(output_path, index=False)


def generate_download_scripts(output_dir='cell_lines'):
    """Generate shell scripts to download BED files."""
    print("\nGenerating download scripts...")
    
    def create_download_command(row, output_path):
        exp_id = row['access_code'].split('/')[2]
        analysis_id = row['def_anal'].split('/')[2]
        
        url = (
            f"https://www.encodeproject.org/batch_download/"
            f"?type=Experiment"
            f"&@id=%2Fexperiments%2F{exp_id}%2F"
            f"&files.analyses.@id=%2Fanalyses%2F{analysis_id}%2F"
            f"&files.preferred_default=true"
        )
        
        filename = f"{row['prot']}_{exp_id}_{analysis_id}.txt"
        return f"curl -o {output_path}/{filename} '{url}'"
    
    for df_path in sorted(glob(f'{output_dir}/*/res_dataframe.csv')):
        df = pd.read_csv(df_path)
        output_path = os.path.join(os.path.dirname(df_path), 'files')
        
        commands = df.apply(lambda row: create_download_command(row, output_path), axis=1)
        
        script_path = df_path.replace('res_dataframe.csv', 'download_files')
        commands.to_csv(script_path, index=False, header=False)
        os.chmod(script_path, 0o711)


def organize_files(output_dir='cell_lines'):
    """Organize downloaded files into folders."""
    print("\nOrganizing files...")
    
    for txt_file in tqdm(sorted(glob(f'{output_dir}/*/files/*.txt'))):
        files_folder = txt_file.replace('.txt', '')
        filescut_folder = txt_file.replace('.txt', '').replace('/files/', '/files_cut/')
        
        os.makedirs(files_folder, exist_ok=True)
        os.makedirs(filescut_folder, exist_ok=True)
        
        shutil.copy(txt_file, files_folder)
        shutil.move(txt_file, filescut_folder)


def main():
    """Main execution."""
    # Load configuration
    df_tf = pd.read_csv('config/transcription_factors.txt', sep='\t')
    tfs = df_tf['TF'].tolist()
    
    df_cell = pd.read_csv('config/cell_lines.txt', sep='\t')
    cell_lines = sorted(df_cell['Cell_line'].tolist())
    
    print(f"Processing {len(tfs)} TFs across {len(cell_lines)} cell lines")
    
    # Execute pipeline
    rename_dict = create_directories(cell_lines)
    download_metadata(cell_lines, tfs, rename_dict)
    parse_metadata()
    create_dataframes()
    generate_download_scripts()
    
    print("\n" + "="*60)
    print("Metadata download complete!")
    print("\nNext steps:")
    print("1. Execute: for file in cell_lines/*/download_files; do ./$file; done")
    print("2. Clean files: sed -i '/bigWig$/d' cell_lines/*/files_cut/*/*.txt")
    print("3. Clean files: sed -i '/bigBed$/d' cell_lines/*/files_cut/*/*.txt")
    print("="*60)


if __name__ == "__main__":
    main()
