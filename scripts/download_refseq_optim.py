"""Download genomes from NCBI."""

import time
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from argparse import ArgumentParser, SUPPRESS

def download_genome(genomes_path, partial_url, end_url):
    command = f"wget -P {genomes_path} {partial_url}/{end_url}_genomic.fna.gz -nv > /dev/null 2>&1"  # && gzip -d {genomes_path}/{end_url}_genomic.fna.gz")
    process = subprocess.Popen(command, shell=True)
    process.wait()  # Wait for the process to finish before returning
    return process.returncode, partial_url

def download_from_ncbi(summary_file: str, genomes_path: str, start: int = 0, max_parallel: int = 5, logging_time: int = 100) -> None:
    """Download all genomes from a file, assumming its a standard NCBI summary file

    Args:
        summary_file (str): Path to a NCBI ftp file
        genomes_path (str): Path to folder to output genomes
        start (int, optional): A line number where to start from. Defaults to 0.
    """
    nb = 0
    time_start = time.time()
    futures = []
    with open(summary_file, "r", encoding='utf-8') as summary_reader:
        # Skipping the two first lines
        next(summary_reader)
        next(summary_reader)
        with ThreadPoolExecutor(max_workers=max_parallel) as executor:
            active_processes = 0  # Track the number of active processes
            for i, line in enumerate(summary_reader):
                if i == 0 :
                    print(line)
                if i > start and 'complete genome' in line.lower():  # Only keeping representative genomes

                    split = line.split()

                    try:
                        https_urls = [elt for elt in split if elt.startswith('https')]
                        if not https_urls:  # No HTTPS URL found
                            print(f"Warning: No HTTPS URL found in line {i}: {line.strip()}")
                            continue  # Skip this line
                        https_wget = https_urls[0]
                        partial_url = https_wget[8:]
                        end_url = partial_url.split('/')[-1]

                        while active_processes >= max_parallel:
                            for future in as_completed(futures):
                                return_code, url = future.result()  # Wait for completion of a future                                active_processes -= 1
                                active_processes -= 1

                        future = executor.submit(download_genome, genomes_path, partial_url, end_url)
                        futures.append(future)
                        active_processes += 1
                        nb += 1
                    except Exception as exc:
                        raise BaseException(f"Job stopped at line {i}.") from exc

                    if nb % logging_time == 0:
                        avg_time = round((time.time() - time_start) / nb, 1)
                        print(f"avg download time per file, at iter {nb} : {avg_time} s")
            # Wait for all futures to complete
            for future in as_completed(futures):
                try:
                    return_code, url = future.result()
                    if return_code != 0:
                        print(f"Error downloading {url}, return code: {return_code}")
                except Exception as exc:
                    print(f"Exception occurred during download: {exc}")

    print(f"Finished downloading {nb} genomes in {(time.time() - time_start)/60} min")



if __name__ == '__main__':
    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "summary_file", type=str, help="Path to NCBI summary ftp file.")
    parser.add_argument('-h', '--help', action='help', default=SUPPRESS,
                        help='Download and rename a set of reference genomes')
    parser.add_argument("-o", "--output", help="Specify a output folder", required=True)
    parser.add_argument("-s", "--start", help="Specify a line number in your summary file to start from",
                        required=False, type=int, default=0)
    parser.add_argument("-n", "--nbThreads", help="nb thread to download",
                        required=False, type=int, default=5)
    args = parser.parse_args()

    download_from_ncbi(args.summary_file, args.output, args.start, args.nbThreads)
    # >> python download_refseq_optim.py  assembly_summary.txt -o /groups/microtaxo/data/refseq/ -s 0 -n 10
    # download_from_ncbi('assembly_summary.txt',
    #                    '/home/hcourtei/Projects/MicroTaxo/codes/data/refseq',
    #                    0,
    #                    5)