"""
ENCODE ChIP-seq Data Downloader
================================
Downloads and organizes ChIP-seq BED files from ENCODE database.
"""

import os
import json
import logging
import time
from typing import Dict, List, Optional
import requests
from glob import glob
from tqdm import tqdm
import pandas as pd
import shutil

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
REQUEST_TIMEOUT = 30  # seconds
MAX_RETRIES = 4
RETRY_BACKOFF = 2  # seconds (will be exponentially increased)


def make_request_with_retry(url: str, headers: Dict[str, str]) -> Optional[requests.Response]:
    """
    Make HTTP request with retry logic and proper error handling.

    Args:
        url: The URL to request
        headers: HTTP headers to include

    Returns:
        Response object if successful, None otherwise
    """
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(url, headers=headers, timeout=REQUEST_TIMEOUT)

            # Check for HTTP errors
            if response.status_code == 403:
                logger.error(f"Access denied (403) for URL: {url}")
                logger.error("This may be due to API rate limiting or authentication requirements")
                return None
            elif response.status_code == 404:
                logger.warning(f"Resource not found (404): {url}")
                return None
            elif response.status_code >= 500:
                logger.warning(f"Server error ({response.status_code}), attempt {attempt + 1}/{MAX_RETRIES}")
                if attempt < MAX_RETRIES - 1:
                    time.sleep(RETRY_BACKOFF ** attempt)
                    continue
                return None

            response.raise_for_status()
            return response

        except requests.exceptions.Timeout:
            logger.warning(f"Request timeout (attempt {attempt + 1}/{MAX_RETRIES})")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_BACKOFF ** attempt)
        except requests.exceptions.ConnectionError as e:
            logger.warning(f"Connection error (attempt {attempt + 1}/{MAX_RETRIES}): {e}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_BACKOFF ** attempt)
        except requests.exceptions.RequestException as e:
            logger.error(f"Request failed: {e}")
            return None

    logger.error(f"Failed after {MAX_RETRIES} attempts")
    return None


def sanitize_name(name: str) -> str:
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


def download_metadata(cell_lines: List[str], tfs: List[str],
                      rename_dict: Dict[str, str], output_dir: str = 'cell_lines') -> None:
    """Download experiment metadata from ENCODE with error handling."""
    headers = {
        'accept': 'application/json',
        'User-Agent': 'ENCODE-ChIP-seq-Downloader/1.0'
    }

    logger.info("Downloading metadata from ENCODE...")
    failed_requests = []

    for cell_line in tqdm(cell_lines, desc="Cell lines"):
        cell_folder = rename_dict[cell_line]

        for tf in tfs:
            url = get_search_url(cell_line, tf)
            response = make_request_with_retry(url, headers)

            if response is None:
                failed_requests.append((cell_line, tf))
                logger.warning(f"Failed to download metadata for {cell_line}/{tf}")
                # Create empty result to avoid downstream errors
                biosample = {'@graph': []}
            else:
                try:
                    biosample = response.json()
                except json.JSONDecodeError as e:
                    logger.error(f"Failed to parse JSON for {cell_line}/{tf}: {e}")
                    failed_requests.append((cell_line, tf))
                    biosample = {'@graph': []}

            output_path = f'{output_dir}/{cell_folder}/json/{tf}.json'
            try:
                with open(output_path, 'w') as f:
                    json.dump(biosample, f, indent=2)
            except IOError as e:
                logger.error(f"Failed to write {output_path}: {e}")

    if failed_requests:
        logger.warning(f"Failed to download {len(failed_requests)} metadata files")
        logger.warning(f"Failed: {failed_requests[:5]}...")  # Show first 5


def parse_metadata(output_dir: str = 'cell_lines') -> None:
    """Parse JSON metadata to extract experiment details with error handling."""
    logger.info("Parsing metadata...")

    for cell_folder in tqdm(sorted(glob(f'{output_dir}/*')), desc="Processing"):
        results = {}

        for json_file in glob(f'{cell_folder}/json/*.json'):
            tf = os.path.basename(json_file).replace('.json', '')
            results[tf] = {}

            try:
                with open(json_file, 'r') as f:
                    biosample = json.load(f)
            except (IOError, json.JSONDecodeError) as e:
                logger.error(f"Failed to read/parse {json_file}: {e}")
                continue

            for exp_idx, experiment in enumerate(biosample.get('@graph', [])):
                try:
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
                except (KeyError, IndexError, TypeError) as e:
                    logger.warning(f"Skipping malformed experiment in {tf}: {e}")
                    continue

        try:
            with open(f'{cell_folder}/res_dict.json', 'w') as f:
                json.dump(results, f, indent=2)
        except IOError as e:
            logger.error(f"Failed to write res_dict.json for {cell_folder}: {e}")


def create_dataframes(output_dir: str = 'cell_lines') -> None:
    """Create summary dataframes from parsed metadata with error handling."""
    logger.info("Creating summary dataframes...")

    for res_dict_path in tqdm(sorted(glob(f'{output_dir}/*/res_dict.json'))):
        try:
            with open(res_dict_path, 'r') as f:
                res_dict = json.load(f)
        except (IOError, json.JSONDecodeError) as e:
            logger.error(f"Failed to read {res_dict_path}: {e}")
            continue

        records = []
        for tf, experiments in res_dict.items():
            for exp_data in experiments.values():
                try:
                    if (exp_data.get('assembly') == 'GRCh38' and not exp_data.get('treatment')):
                        records.append({
                            'prot': tf,
                            'lab': exp_data.get('lab'),
                            'desc': exp_data.get('description'),
                            'n_def': exp_data.get('n_default'),
                            'def_anal': exp_data.get('default_analysis'),
                            'access_code': exp_data.get('accession_code')
                        })
                except (KeyError, TypeError) as e:
                    logger.warning(f"Skipping malformed experiment data in {tf}: {e}")
                    continue

        if not records:
            logger.warning(f"No valid records found for {res_dict_path}")
            # Create empty dataframe to avoid downstream errors
            df = pd.DataFrame(columns=['prot', 'lab', 'desc', 'n_def', 'def_anal', 'access_code'])
        else:
            df = pd.DataFrame(records)

        output_path = res_dict_path.replace('res_dict.json', 'res_dataframe.csv')
        try:
            df.to_csv(output_path, index=False)
        except IOError as e:
            logger.error(f"Failed to write {output_path}: {e}")


def generate_download_scripts(output_dir: str = 'cell_lines') -> None:
    """Generate shell scripts to download BED files with error handling."""
    logger.info("Generating download scripts...")

    def create_download_command(row, output_path):
        try:
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
        except (KeyError, IndexError, AttributeError) as e:
            logger.warning(f"Failed to create download command for row: {e}")
            return ""

    for df_path in sorted(glob(f'{output_dir}/*/res_dataframe.csv')):
        try:
            df = pd.read_csv(df_path)
        except (IOError, pd.errors.EmptyDataError) as e:
            logger.error(f"Failed to read {df_path}: {e}")
            continue

        if df.empty:
            logger.warning(f"Empty dataframe at {df_path}, skipping script generation")
            continue

        output_path = os.path.join(os.path.dirname(df_path), 'files')

        commands = df.apply(lambda row: create_download_command(row, output_path), axis=1)
        # Filter out empty commands
        commands = commands[commands != ""]

        if commands.empty:
            logger.warning(f"No valid commands generated for {df_path}")
            continue

        script_path = df_path.replace('res_dataframe.csv', 'download_files')
        try:
            commands.to_csv(script_path, index=False, header=False)
            os.chmod(script_path, 0o711)
        except IOError as e:
            logger.error(f"Failed to write download script {script_path}: {e}")


def organize_files(output_dir: str = 'cell_lines') -> None:
    """Organize downloaded files into folders with error handling."""
    logger.info("Organizing files...")

    txt_files = sorted(glob(f'{output_dir}/*/files/*.txt'))
    if not txt_files:
        logger.warning("No files found to organize")
        return

    for txt_file in tqdm(txt_files, desc="Organizing"):
        try:
            files_folder = txt_file.replace('.txt', '')
            filescut_folder = txt_file.replace('.txt', '').replace('/files/', '/files_cut/')

            os.makedirs(files_folder, exist_ok=True)
            os.makedirs(filescut_folder, exist_ok=True)

            shutil.copy(txt_file, files_folder)
            shutil.move(txt_file, filescut_folder)
        except (IOError, shutil.Error) as e:
            logger.error(f"Failed to organize {txt_file}: {e}")


def main() -> None:
    """Main execution with comprehensive error handling."""
    try:
        # Load configuration
        logger.info("Loading configuration files...")
        try:
            df_tf = pd.read_csv('tables/tfs.txt', sep='\t')
            tfs = df_tf['TF'].tolist()
        except (IOError, pd.errors.EmptyDataError, KeyError) as e:
            logger.error(f"Failed to load transcription factors from tables/tfs.txt: {e}")
            return

        try:
            df_cell = pd.read_csv('tables/cells.txt', sep='\t')
            cell_lines = sorted(df_cell['Cell_line'].tolist())
        except (IOError, pd.errors.EmptyDataError, KeyError) as e:
            logger.error(f"Failed to load cell lines from tables/cells.txt: {e}")
            return

        if not tfs:
            logger.error("No transcription factors found in configuration")
            return
        if not cell_lines:
            logger.error("No cell lines found in configuration")
            return

        logger.info(f"Processing {len(tfs)} TFs across {len(cell_lines)} cell lines")

        # Execute pipeline
        try:
            rename_dict = create_directories(cell_lines)
            download_metadata(cell_lines, tfs, rename_dict)
            parse_metadata()
            create_dataframes()
            generate_download_scripts()
        except Exception as e:
            logger.error(f"Pipeline execution failed: {e}")
            raise

        logger.info("="*60)
        logger.info("Metadata download complete!")
        logger.info("\nNext steps:")
        logger.info("1. Execute: for file in cell_lines/*/download_files; do ./$file; done")
        logger.info("2. Clean files: sed -i '/bigWig$/d' cell_lines/*/files_cut/*/*.txt")
        logger.info("3. Clean files: sed -i '/bigBed$/d' cell_lines/*/files_cut/*/*.txt")
        logger.info("="*60)

    except KeyboardInterrupt:
        logger.warning("\nInterrupted by user")
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()
