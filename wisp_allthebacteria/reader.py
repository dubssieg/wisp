from collections import Counter
from concurrent.futures import (
    FIRST_COMPLETED,
    # ProcessPoolExecutor,
    ThreadPoolExecutor,  # deprecated
    wait,
)

# from concurrent.futures.process import BrokenProcessPool
# import multiprocessing

import gc
import os
import queue
from itertools import product
import logging
from pathlib import Path
import tarfile
import tempfile

# from pympler import asizeof
from loky import get_reusable_executor

from Bio import SeqIO
import psutil
from metadata import Metadata
from utils import format_size, system_stats

LOG = logging.getLogger(__name__)

MAX_RUNNING_TASKS_FASTA = 64
MAX_RUNNING_TASKS_SEQUENCE = 64  # deprecated


class Reader:
    def __init__(
        self,
        metadata: Metadata,
        num_workers: int,
        # sequences_threads: int = 4,
    ):
        LOG.debug(f"Reader({locals()})")
        self._md = metadata
        self._num_workers = num_workers
        # self._sequences_threads = sequences_threads
        self._parent_pid = os.getpid()

    def process_file(
        self,
        file_path: str | Path,
        kmer_size: int,
        window_size: int,
        step: int,
        full: bool = False,
    ) -> dict:
        """Just call process_archive of process_fasta, based on file suffix."""
        file_path = Path(file_path).resolve()
        LOG.debug(f"processing file: {file_path} with {self._num_workers} workers")
        suffix = file_path.suffix
        if suffix == ".xz":
            return self.process_archive(
                file_path,
                kmer_size=kmer_size,
                window_size=window_size,
                step=step,
                full=full,
            )
        elif suffix == ".fa":
            return self.process_fasta(
                file_path,
                kmer_size=kmer_size,
                window_size=window_size,
                step=step,
                full=full,
            )
        LOG.debug(f"processed file: {file_path}")

    def process_archive(
        self,
        archive_path: Path | str,
        kmer_size: int,
        window_size: int,
        step: int,
        full: bool = False,
    ):
        """Extract and process an archive."""
        archive_path = Path(archive_path).resolve()
        archive_name, archive_size = archive_path.name, format_size(
            archive_path.stat().st_size
        )
        LOG.info(f"[{archive_name}] Processing archive: {archive_path}")
        LOG.info(f"[{archive_name}] SIZE: {archive_size}")
        LOG.debug(f"SYSTEM: {system_stats(as_str=True)}")

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir).resolve()
            with tarfile.open(archive_path, "r:xz") as tar:
                tar.extractall(temp_dir)

            extracted_files = list(temp_dir.rglob("*.fa"))
            LOG.debug(f"[{archive_name}] {len(extracted_files)} FASTA files extracted")

            with get_reusable_executor(max_workers=self._num_workers) as executor:
                results = executor.map(
                    Reader.process_fasta,
                    extracted_files,
                    [kmer_size] * len(extracted_files),
                    [window_size] * len(extracted_files),
                    [step] * len(extracted_files),
                    [full] * len(extracted_files),
                )

            sequences_data = list(results)

        merged_data = {}
        for sequence_data in sequences_data:
            for file_id, data in sequence_data.items():
                file_md = self._md[file_id]
                match len(file_md):
                    case 0:
                        tax_id = "no-tax-id"
                    case 1:
                        tax_id = file_md[0]["TaxId"]
                    case _:
                        tax_id = f"multiple-{'|'.join(m['TaxId'] if m else 'no-tax-id' for m in file_md)}"
                merged_data.setdefault(tax_id, {"counters": [], "sources": []})
                merged_data[tax_id]["counters"].extend(data["counters"])
                merged_data[tax_id]["sources"].extend(data["sources"])

        LOG.debug(f"[{archive_path}] Processed {len(merged_data)} tax_id(s)")
        return {"archive": archive_path, "merged_data": merged_data}

    def process_archive_old_version(
        self,
        archive_path: Path | str,
        kmer_size: int,
        window_size: int,
        step: int,
        full: bool = False,
    ):
        """Extract and process an archive."""
        archive_path = Path(archive_path).resolve()
        archive_name, archive_size = archive_path.name, format_size(
            archive_path.stat().st_size
        )
        LOG.info(f"[{archive_path}] Processing archive: {archive_path}")
        LOG.info(f"[{archive_path}] SIZE: {archive_size}")
        LOG.debug(f"SYSTEM: {system_stats(as_str=True)}")

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir).resolve()
            with tarfile.open(archive_path, "r:xz") as tar:
                tar.extractall(temp_dir)

            extracted_files = list(temp_dir.rglob("*.fa"))
            LOG.debug(f"[{archive_name}] Extracted: {len(extracted_files)} FASTA files")

            sequences_data, task_queue = [], queue.Queue()
            for file_path in extracted_files:
                task_queue.put(file_path)

            # with ProcessPoolExecutor(max_workers=self._num_workers) as executor:
            with get_reusable_executor(max_workers=self._num_workers) as executor:
                running_futures, fasta_count = {}, 0

                try:
                    while not task_queue.empty() or running_futures:
                        while (
                            not task_queue.empty()
                            and len(running_futures) < MAX_RUNNING_TASKS_FASTA
                        ):
                            file_path = task_queue.get()
                            future = executor.submit(
                                Reader.process_fasta,
                                file_path=file_path,
                                kmer_size=kmer_size,
                                window_size=window_size,
                                step=step,
                                # num_workers=self._sequences_threads,
                                # parent_pid=self._parent_pid,
                                full=full,
                            )
                            running_futures[future] = file_path

                        done, _ = wait(
                            running_futures.keys(), return_when=FIRST_COMPLETED
                        )

                        LOG.debug(f"SYSTEM: {system_stats(as_str=True)}")
                        for future in done:
                            file_path = running_futures.pop(future)
                            try:
                                sequences_data.append(future.result(timeout=300))

                                # prevent zombies and memory leaks
                                del future
                                gc.collect()
                                done_pids = [
                                    p.pid
                                    for p in psutil.Process().children(recursive=True)
                                    if p.status() == psutil.STATUS_ZOMBIE
                                ]
                                for pid in done_pids:
                                    try:
                                        os.waitpid(pid, 0)
                                    except ChildProcessError:
                                        pass
                                LOG.debug(f"[{file_path.name}] Fasta processed")
                            except TimeoutError:
                                LOG.error(
                                    f"Timeout on {file_path.name}, retrying later..."
                                )
                                task_queue.put(file_path)
                                fasta_count -= 1
                            except Exception:
                                LOG.exception(
                                    f"Error processing FASTA [{archive_name}] {file_path.name}"
                                )
                                raise

                            fasta_count += 1
                            LOG.debug(
                                f"[{file_path.name}] All fasta sequences counted - {fasta_count}/{len(extracted_files)}"
                            )

                finally:
                    LOG.debug(
                        f"[{archive_name}] Forcing ProcessPoolExecutor shutdown + Process cleanup"
                    )
                    for future in running_futures.keys():
                        try:
                            if not future.done():
                                future.cancel()
                        except Exception:
                            pass  # Future already terminated
                    executor.shutdown(wait=True, cancel_futures=True)
                    gc.collect()

                    # Cleanup potential zombie processes
                    current_process = psutil.Process(os.getpid())
                    for child in current_process.children(recursive=True):
                        LOG.warning(
                            f"Killing leftover process: {child.pid} - {child.cmdline()}"
                        )
                        child.terminate()

                    # Wait 5s and force killing if needed
                    _, alive = psutil.wait_procs(current_process.children(), timeout=5)
                    for child in alive:
                        LOG.error(f"Force killing process: {child.pid}")
                        child.kill()
                LOG.debug(f"[{archive_name}] Closing ProcessPoolExecutor")
            LOG.debug(f"{len(sequences_data)} sequences counted - merging by tax_id")

            merged_data = {}
            for sequence_data in sequences_data:
                for file_id, data in sequence_data.items():
                    file_md = self._md[file_id]
                    match len(file_md):
                        case 0:
                            tax_id = "no-tax-id"
                        case 1:
                            tax_id = file_md[0]["TaxId"]
                        case _:
                            tax_id = f"multiple-{'|'.join(m['TaxId'] if m else 'no-tax-id' for m in file_md)}"
                    merged_data.setdefault(tax_id, {"counters": [], "sources": []})
                    merged_data[tax_id]["counters"].extend(data["counters"])
                    merged_data[tax_id]["sources"].extend(data["sources"])

            LOG.debug(f"Sequences counters merged: {len(merged_data)} tax_id(s)")
            LOG.debug(f"[{archive_path}] Extracted files deleted)")

            return {"archive": archive_path, "merged_data": merged_data}

    @staticmethod
    def process_fasta(
        file_path: Path | str,
        kmer_size: int,
        window_size: int,
        step: int,
        full: bool = False,
    ) -> dict:
        """Simpler Fasta counting method, hopefully less bugged"""
        file_path = Path(file_path).resolve()
        file_name, file_size = file_path.name, format_size(file_path.stat().st_size)
        LOG.debug(f"[{file_name}] Processing Fasta - SIZE: {file_size})")

        try:
            sequences = Reader._read_fasta(file_path)
            LOG.debug(f"[{file_name}] {len(sequences)} sequences to count")
        except Exception:
            LOG.exception(f"Error while reading fasta file: {file_path}")
            raise

        source_id_to_data = {}
        for sequence in sequences:
            kmer_count, source = Reader.process_sequence(
                sequence=sequence,
                kmer_size=kmer_size,
                window_size=window_size,
                step=step,
                full=full,
            )
            source_id = source["id"]
            source_id_to_data.setdefault(source_id, {"counters": [], "sources": []})
            source_id_to_data[source_id]["counters"].extend(kmer_count)
            source_id_to_data[source_id]["sources"].append(source)

        LOG.debug(f"[{file_name}] Return all counters & sources")
        return source_id_to_data

    @staticmethod
    def process_fasta_mt(
        file_path: Path | str,
        kmer_size: int,
        window_size: int,
        step: int,
        num_workers: int,
        full: bool = False,
    ) -> dict:
        """Process a FASTA file using controlled multithreading."""
        file_path = Path(file_path).resolve()
        file_name, file_size = file_path.name, format_size(file_path.stat().st_size)
        LOG.debug(f"[{file_name}] Processing Fasta - SIZE: {file_size})")

        try:
            sequences = Reader._read_fasta(file_path)
            LOG.debug(f"[{file_name}] {len(sequences)} sequences to count")
        except Exception:
            LOG.exception(f"Error while reading fasta file: {file_path}")
            raise

        source_id_to_data, task_queue = {}, queue.Queue()
        for seq in sequences:
            task_queue.put(seq)

        # LOG.debug(f"SYSTEM: {system_stats(pid=parent_pid, as_str=True)}")
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            running_futures = set()

            try:
                while not task_queue.empty() or running_futures:
                    while (
                        not task_queue.empty()
                        and len(running_futures) < MAX_RUNNING_TASKS_SEQUENCE
                    ):
                        sequence = task_queue.get()
                        running_futures.add(
                            executor.submit(
                                Reader.process_sequence,
                                sequence=sequence,
                                kmer_size=kmer_size,
                                window_size=window_size,
                                step=step,
                                full=full,
                            )
                        )

                    done, running_futures = wait(
                        running_futures, return_when=FIRST_COMPLETED
                    )

                    for future in list(done):
                        try:
                            if future in running_futures:
                                running_futures.discard(future)
                            kmer_count, source = future.result()
                            source_id = source["id"]
                            source_id_to_data.setdefault(
                                source_id, {"counters": [], "sources": []}
                            )
                            source_id_to_data[source_id]["counters"].extend(kmer_count)
                            source_id_to_data[source_id]["sources"].append(source)
                        except Exception:
                            LOG.exception(file_name)
                            raise
                        finally:
                            running_futures.discard(future)

            finally:
                LOG.debug(f"[{file_name}] Fasta - Closing ThreadPoolExecutor")
                executor.shutdown(wait=True, cancel_futures=True)

        LOG.debug(f"[{file_name}] Fasta - Cleaning task queue and sequences")
        while not task_queue.empty():
            task_queue.get()
            task_queue.task_done()

        del sequences
        gc.collect()

        LOG.debug(f"[{file_name}] Return all counters & sources")
        return source_id_to_data

    @staticmethod
    def process_sequence(
        sequence: str,
        kmer_size: int,
        window_size: int,
        step: int,
        full: bool,
    ) -> tuple[dict, dict]:
        """need to be static for ProcessPoolExecutor"""

        kmer_count = Reader._counter(
            entry=sequence["sequence"],
            kmer_size=kmer_size,
            window_size=window_size,
            step=step,
            full=full,
        )

        source = {k: v for k, v in sequence.items() if k != "sequence"}
        return kmer_count, source

    @staticmethod
    def _counter(
        entry: str,
        kmer_size: int = 4,
        window_size: int = 10000,
        step: int = 5000,
        full: bool = False,
    ) -> dict:
        complements = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "U": "A",
            "R": "Y",
            "Y": "R",
            "K": "M",
            "M": "K",
            "S": "W",
            "W": "S",
            "B": "V",
            "V": "B",
            "D": "H",
            "H": "D",
            "N": "N",
        }

        degenerate_map = {
            "U": ["T"],
            "R": ["G", "A"],
            "Y": ["C", "T"],
            "K": ["G", "T"],
            "M": ["A", "C"],
            "S": ["G", "C"],
            "W": ["A", "T"],
            "B": ["G", "T", "C"],
            "D": ["G", "T", "A"],
            "H": ["A", "T", "C"],
            "V": ["G", "A", "C"],
            "N": ["A", "T", "C", "G"],
        }

        seq_len = len(entry)
        if full or seq_len <= window_size:
            windows = [(0, seq_len)]
        else:
            num_windows = (seq_len - window_size) // step + 1
            step = (
                (seq_len - window_size) // (num_windows - 1)
                if num_windows > 1
                else step
            )
            windows = [
                (i * step, min(i * step + window_size, seq_len))
                for i in range(num_windows)
            ]
            if windows[-1][1] < seq_len:
                windows.append((seq_len - window_size, seq_len))

        window_counters = []

        for start, end in windows:
            kmers = (
                entry[i : i + kmer_size] for i in range(start, end - kmer_size + 1)
            )

            counts = Counter(kmers)
            rev_counts = Counter(
                {Reader._revcomp(k, compl=complements): v for k, v in counts.items()}
            )
            counts += rev_counts

            for filtered_kmer in (alpha * kmer_size for alpha in "ATCG"):
                counts.pop(filtered_kmer, None)

            counts_purged = {}
            for key, count in counts.items():
                list_of_keys = [
                    "".join(item)
                    for item in product(*[degenerate_map.get(x, [x]) for x in key])
                ]
                kmer_number = count // len(list_of_keys)
                for prob_key in list_of_keys:
                    counts_purged[prob_key] = (
                        counts_purged.get(prob_key, 0) + kmer_number
                    )

            window_counters.append(counts_purged)

        return window_counters

    @staticmethod
    def _revcomp(string: str, compl=None) -> str:
        if compl is None:
            compl = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

        try:
            result = "".join(compl[s] for s in reversed(string))
        except KeyError:
            LOG.exception(f"revcom key error: {string}")
            raise
        return result

    @staticmethod
    def _read_fasta(file_path: Path | str) -> list:
        """Just read and parse file, no processing.
        Return a list of:
            'id' = 'SAMD00013333'
            'contig' = 'contig00001'
            'file' = 'SAMEA1561896.fa'
            'sequence' = 'GGAGGGAACAGCGGGGCGGGCGGCGT..."""
        file_path = Path(file_path).resolve()
        if not file_path.exists():
            raise FileNotFoundError(f"File {file_path} not found.")

        sequences = []
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                if ".contig" in record.id:
                    r_id, r_contig = record.id.split(".")[:2]
                else:
                    r_id, r_contig = record.id, ""
                sequences.append(
                    {
                        "id": r_id,
                        "contig": r_contig,
                        "file": file_path.stem,
                        "sequence": str(record.seq),
                    }
                )
        return sequences


if __name__ == "__main__":
    pass
