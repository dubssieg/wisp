from collections import Counter
from concurrent.futures import (
    FIRST_COMPLETED,
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    wait,
)

# from concurrent.futures.process import BrokenProcessPool
import queue
from itertools import product
import logging
from pathlib import Path
import tarfile
import tempfile

from Bio import SeqIO
from metadata import Metadata
from utils import format_size, system_stats

LOG = logging.getLogger(__name__)

MAX_RUNNING_TASKS_FASTA = 64
MAX_RUNNING_TASKS_SEQUENCE = 64


class Reader:
    def __init__(
        self,
        metadata: Metadata,
        num_workers: int = 20,
        sequences_threads: int = 4,
    ):
        self._md = metadata
        self._num_workers = num_workers
        self._sequences_threads = sequences_threads

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
        LOG.debug(
            f"processing file: {file_path} with {self._num_workers} workers and {self._sequences_threads} threads per file"
        )
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
        archive_name = archive_path.name
        archive_size = format_size(archive_path.stat().st_size)
        LOG.info(f"[{archive_path}] Processing archive: {archive_path}")
        LOG.info(f"[{archive_path}] SIZE: {archive_size}")
        LOG.debug(f"SYSTEM: {system_stats(as_str=True)}")

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir).resolve()
            with tarfile.open(archive_path, "r:xz") as tar:
                tar.extractall(path=temp_dir)

            extracted_files = list(temp_dir.rglob("*.fa"))
            LOG.debug(f"[{archive_name}] Extracted: {len(extracted_files)} FASTA files")

            merged_data = {}

            LOG.debug(
                f"[{archive_name}] Processing {self._num_workers} FASTA files in parallel"
            )
            fasta_count = 0

            # all files in queue
            task_queue = queue.Queue()
            for file_path in extracted_files:
                task_queue.put(file_path)

            # LOG.debug(f"SYSTEM: {system_stats(as_str=True)}")
            with ProcessPoolExecutor(max_workers=self._num_workers) as executor:
                running_futures = {}

                while not task_queue.empty() or running_futures:
                    while (
                        not task_queue.empty()
                        and len(running_futures) < MAX_RUNNING_TASKS_FASTA
                    ):
                        file_path = task_queue.get()
                        future = executor.submit(
                            self.process_fasta,
                            file_path=file_path,
                            kmer_size=kmer_size,
                            window_size=window_size,
                            step=step,
                            full=full,
                        )
                        running_futures[future] = file_path

                    done, _ = wait(running_futures.keys(), return_when=FIRST_COMPLETED)

                    for future in done:
                        try:
                            file_data = future.result()
                            file_name = running_futures.pop(future).name
                            LOG.debug(f"[{file_name}] Fasta processed")

                            # merge
                            for tax_id, data in file_data.items():
                                if tax_id not in merged_data:
                                    merged_data[tax_id] = {
                                        "counters": [],
                                        "sources": [],
                                    }
                                merged_data[tax_id]["counters"].extend(data["counters"])
                                merged_data[tax_id]["sources"].extend(data["sources"])

                            fasta_count += 1
                            LOG.debug(
                                f"[{file_name}] Data merged - {fasta_count} / {len(extracted_files)}"
                            )
                        # except BrokenProcessPool:
                        #     LOG.exception(
                        #         f"BrokenProcessPool error processing FASTA [{archive_name}] {file_name}"
                        #     )
                        #     raise
                        except Exception:
                            LOG.exception(
                                f"Error processing FASTA [{archive_name}] {file_name}"
                            )
                            raise
                        LOG.debug(
                            f"[{archive_name}] One future processed (data merged)"
                        )
                    LOG.debug(f"[{archive_name}] All done futures processed")
                LOG.debug(
                    f"[{archive_name}] Queue empty - all future processed - closing ProcessPoolExecutor"
                )
            LOG.debug(f"[{archive_name}] ProcessPoolExecutor closed")
        LOG.debug(f"[{archive_path}] Extracted files deleted")

        return {"archive": archive_path, "merged_data": merged_data}

    def process_fasta(
        self,
        file_path: str | Path,
        kmer_size: int,
        window_size: int,
        step: int,
        full: bool = False,
    ) -> dict:
        """Process a FASTA file using controlled multithreading."""
        file_path = Path(file_path).resolve()
        file_name = file_path.name
        file_size = format_size(file_path.stat().st_size)
        LOG.debug(f"[{file_name}] Processing Fasta - SIZE: {file_size})")
        LOG.debug(f"SYSTEM: {system_stats(as_str=True)}")

        try:
            sequences = self._read_fasta(file_path)
            LOG.debug(f"[{file_name}] (Worker) {len(sequences)} sequences read")
        except Exception:
            LOG.exception(f"Error while reading fasta file: {file_path}")
            raise

        tax_id_to_data = {}

        LOG.debug(f"[{file_name}] Counting - {self._sequences_threads} threads")

        # Queue all sequences
        task_queue = queue.Queue()

        for seq in sequences:
            task_queue.put(seq)

        with ThreadPoolExecutor(max_workers=self._sequences_threads) as executor:
            running_futures = set()

            while not task_queue.empty() or running_futures:
                while (
                    not task_queue.empty()
                    and len(running_futures) < MAX_RUNNING_TASKS_SEQUENCE
                ):
                    seq = task_queue.get()
                    future = executor.submit(
                        self.process_sequence,
                        sequence=seq,
                        kmer_size=kmer_size,
                        window_size=window_size,
                        step=step,
                        full=full,
                        md=self._md[seq["id"]],
                    )
                    running_futures.add(future)

                # wait for at least one task to end before adding more
                done, running_futures = wait(
                    running_futures, return_when=FIRST_COMPLETED
                )

                for future in done:
                    try:
                        tax_id, kmer_count, source = future.result()
                        if tax_id not in tax_id_to_data:
                            tax_id_to_data[tax_id] = {"counters": [], "sources": []}

                        tax_id_to_data[tax_id]["counters"].extend(kmer_count)
                        tax_id_to_data[tax_id]["sources"].extend(
                            [source] * len(kmer_count)
                        )

                    except Exception:
                        LOG.exception(file_name)
                        raise
            LOG.debug(f"[{file_name}] Sequence - Closing ThreadPoolExecutor")
        LOG.debug(f"[{file_name}] Sequence - ThreadPoolExecutor closed")

        return tax_id_to_data

    @staticmethod
    def process_sequence(
        sequence: str,
        kmer_size: int,
        window_size: int,
        step: int,
        full: bool,
        md: dict,
    ):
        """need to be static for ProcessPoolExecutor"""
        match len(md):
            case 0:
                tax_id = "no-tax-id"
            case 1:
                tax_id = md[0]["TaxId"]
            case _:
                tax_id = (
                    f"multiple-{'|'.join(m['TaxId'] if m else 'no-tax-id' for m in md)}"
                )

        kmer_count = Reader._counter(
            entry=sequence["sequence"],
            kmer_size=kmer_size,
            window_size=window_size,
            step=step,
            full=full,
        )

        source = {k: v for k, v in sequence.items() if k != "sequence"}
        return tax_id, kmer_count, source

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
