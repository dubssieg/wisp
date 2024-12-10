"""Download genomes from NCBI."""
import os
import time
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from argparse import ArgumentParser

log_file = "genome_downloader.log"
if os.path.exists(log_file):
    os.remove(log_file)

logging.basicConfig(
    level=logging.INFO,  # Change to DEBUG for more detailed logs
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%H:%M",  # Format de l'heure : heures et minutes
    handlers=[logging.FileHandler("genome_downloader.log"),  # Logs to a file
              logging.StreamHandler()  # Logs to consol
             ]
)
logger = logging.getLogger(__name__)

def download_genome(genomes_path, https_wget):
    partial_url = https_wget[8:]
    end_url = partial_url.split('/')[-1]
    name_file = f"{partial_url}/{end_url}_genomic.fna.gz"
    command = f"wget -P {genomes_path} {name_file} -nv > /dev/null 2>&1"  # no logs && gzip -d {genomes_path}/{end_url}_genomic.fna.gz")
    process = subprocess.Popen(command, shell=True)
    process.wait()  # Wait for the process to finish before returning
    return process.returncode, partial_url


def download_grouped_url_from_ncbi(grouped_url: list,
                                   genomes_dir: str,
                                   max_parallel: int = 5,
                                   logging_time: int = 100) -> None:
    """Download all genomes from a file, assumming its a standard NCBI summary file

    Args:
        genomes_dir (str): Path to folder to output genomes
    """
    nb = 0
    time_start = time.time()
    futures = []
    with ThreadPoolExecutor(max_workers=max_parallel) as executor:
        active_processes = 0  # Track the number of active processes
        for i, https_wget in enumerate(grouped_url):

            while active_processes >= max_parallel:
                for future in as_completed(futures):
                    return_code, url = future.result()  # Wait for completion of a future                                active_processes -= 1
                    active_processes -= 1

            future = executor.submit(download_genome, genomes_dir, https_wget)
            futures.append(future)
            active_processes += 1
            nb += 1

            if nb % logging_time == 0:
                avg_time = round((time.time() - time_start) / nb, 1)
                logger.info(f"Avg download time per file, at iter {nb}: {avg_time} s")        # Wait for all futures to complete
        for future in as_completed(futures):
            try:
                return_code, url = future.result()
                if return_code != 0:
                    logger.error(f"Error downloading {url}, return code: {return_code}")
            except Exception as exc:
                logger.exception(f"Exception occurred during download: {exc}")

def download_all_url_from_ncbi(summary_file: str,
                               genomes_rootdir: str,
                               chunk_size: int = 5000,
                               max_parallel: int = 5,
                               logging_time: int = 100) -> None:

    all_url_file = get_valid_urls_from_ncbi(summary_file)
    grouped_urls_list = [all_url_file[i:i + chunk_size] for i in range(0, len(all_url_file), chunk_size)]

    for id_group, grouped_urls in enumerate(grouped_urls_list):
        genomes_dir = f"{genomes_rootdir}/group_{id_group}"
        start = time.time()
        logger.info("=" * 60)
        logger.info(f"Download group {id_group} of url to {genomes_dir}")
        download_grouped_url_from_ncbi(grouped_urls, genomes_dir, max_parallel, logging_time)
        duration = time.time() - start
        logger.info(f"Downloaded group {id_group} of size {chunk_size} in {round(duration / 60)} min ")
        logger.info("=" * 60)


def get_valid_urls_from_ncbi(summary_file: str):
    with open(summary_file, "r", encoding='utf-8') as summary_reader:
        next(summary_reader)
        next(summary_reader)
        all_url_file = []
        nb_complete_with_https = 0
        nb_complete = 0
        for idx, line in enumerate(summary_reader):

            if 'complete genome' in line.lower():
                nb_complete += 1
                split = line.split()
                https_urls = [elt for elt in split if elt.startswith('https')]

                if not https_urls:  # No HTTPS URL found
                    # print(f"Warning: No HTTPS URL found in line {idx}: {line.strip()}")
                    continue  # Skip this line
                nb_complete_with_https += 1
                https_wget = https_urls[0]
                all_url_file.append(https_wget)
        nb_lines = idx
        logger.info(f"Extracted from {summary_file} Nb entries: {nb_lines}, "
                    f"nb with 'complete genome': {nb_complete}, and with https: {nb_complete_with_https}")

    return all_url_file


if __name__ == '__main__':
    parser = ArgumentParser(add_help=False)
    parser.add_argument("summary_file", type=str, help="Path to NCBI summary ftp file.")
    parser.add_argument("-o", "--output", help="Specify a output folder", required=True)
    args = parser.parse_args()
    # > python download_refseq_optim.py assembly_summary.txt -o /groups/microtaxo/data/refseq

    # summary_file = "assembly_summary.txt" #args.summary_file
    # genomes_rootdir = args.outputs #"/home/hcourtei/Projects/MicroTaxo/codes/data/refseq"
    # all_url_file = get_valid_urls_from_ncbi("assembly_summary.txt")
    # chunk_size = 5000
    # grouped_urls_list = [all_url_file[i:i + chunk_size] for i in range(0, len(all_url_file), chunk_size)]

    download_all_url_from_ncbi(summary_file=args.summary_file ,
                               genomes_rootdir= args.output,
                               chunk_size=1000,
                               max_parallel=5,
                               logging_time=200)
