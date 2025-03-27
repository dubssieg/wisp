import concurrent.futures
import gc
import logging
import os
import sys
import tarfile
import tempfile
from collections import Counter
from itertools import product
from pathlib import Path
import time
from diskcache import Cache

from Bio import SeqIO
import concurrent
from loky import get_reusable_executor
from metadata import Metadata
from utils import format_size, cleanup_zombie_processes, compress

LOG = logging.getLogger(__name__)

MAX_RUNNING_TASKS_FASTA = 64


class Reader:
    def __init__(self, metadata: Metadata, num_workers: int):
        LOG.debug(f"Reader({locals()})")
        self._md = metadata
        self._num_workers = num_workers

    def process_file(
        self,
        file_path: str | Path,
        kmer_sizes: list[int],
        window_size: int,
        step: int,
        full: bool = False,
        batch_size: int | None = None,
        compression: str | None = None,
        merged_data_as_db: bool = False,
    ) -> dict:
        """Just call process_archive of process_fasta, based on file suffix."""
        file_path = Path(file_path).resolve()
        suffix = file_path.suffix
        if suffix == ".xz":
            return self.process_archive(
                file_path,
                kmer_sizes=kmer_sizes,
                window_size=window_size,
                step=step,
                full=full,
                batch_size=batch_size,
                compression=compression,
                return_db=merged_data_as_db,
            )
        elif suffix == ".fa":
            return self.process_fasta(
                file_path,
                kmer_sizes=kmer_sizes,
                window_size=window_size,
                step=step,
                full=full,
                compression=compression,
            )
        LOG.debug(f"processed file: {file_path}")

    def process_archive(
        self,
        archive_path: Path | str,
        kmer_sizes: list[int],
        window_size: int,
        step: int,
        full: bool = False,
        batch_size: int = None,
        compression: str | None = None,
        return_db: bool = True,
    ):
        """Extract and process an archive."""
        archive_path = Path(archive_path).resolve()
        archive_name = archive_path.name

        with tempfile.TemporaryDirectory("_wisp_fasta") as temp_dir:
            temp_dir = Path(temp_dir).resolve()
            archive_size = format_size(archive_path.stat().st_size)
            LOG.info(f"[{archive_name}] Extracting archive")
            LOG.info(f"[{archive_name}] SIZE: {archive_size}, temp dir: {temp_dir}")
            # extraction de l'archive
            with tarfile.open(archive_path, "r:xz") as tar:
                tar.extractall(temp_dir)

            extracted_files = list(temp_dir.rglob("*.fa"))
            LOG.debug(f"[{archive_name}] {len(extracted_files)} FASTA files extracted")

            results = Cache(temp_dir / "results", size_limit=sys.maxsize)
            fasta_count = 0

            if batch_size is None:
                batches = [extracted_files]
            else:
                batches = [
                    extracted_files[i : i + batch_size]
                    for i in range(0, len(extracted_files), batch_size)
                ]
            for batch_i, batch in enumerate(batches):
                LOG.debug(
                    f"[{archive_name}] Batch {batch_i + 1} / {len(batches)} - {len(batch)} fasta files"
                )

                with get_reusable_executor(max_workers=self._num_workers) as executor:
                    # with concurrent.futures.ProcessPoolExecutor(
                    #     max_workers=self._num_workers
                    # ) as executor:
                    futures = {
                        executor.submit(
                            Reader.process_fasta,
                            file_path=file_path,
                            kmer_sizes=kmer_sizes,
                            window_size=window_size,
                            step=step,
                            full=full,
                            compression=compression,
                        ): file_path
                        for file_path in batch
                    }
                    try:
                        for future in concurrent.futures.as_completed(futures):
                            try:
                                result = future.result(timeout=60)  # TODO: conf

                                results[fasta_count] = result
                                fasta_count += 1

                                LOG.debug(
                                    f"[{archive_name}] {futures[future].name}: {fasta_count} / {len(extracted_files)}"
                                )

                            except concurrent.futures.TimeoutError:
                                LOG.error(f"[{archive_name}] Timeout")
                            except Exception:
                                LOG.exception(f"[{archive_name}]")
                                raise
                        LOG.debug(
                            f"[{archive_name}] Batch {batch_i + 1} / {len(batches)} done. Terminating worker pool..."
                        )
                    finally:
                        LOG.debug(f"[{archive_name}] Shutting down executor...")
                        # force workers to stop
                        executor.shutdown(wait=True, kill_workers=True)
                        LOG.debug(f"[{archive_name}] Executor shut down.")

                        # give time for zombies to appear
                        time.sleep(2)
                        LOG.debug(f"[{archive_name}] Checking for zombie processes...")
                        cleanup_zombie_processes(os.getpid())

                        # free mem
                        gc.collect()
                        LOG.debug(
                            f"[{archive_name}] Cleanup complete - Workers pool terminated"
                        )

                LOG.debug(f"[{archive_name}] Workers pool terminated")

            LOG.debug(
                f"[{archive_name}] All {len(extracted_files)} Fasta files done - sorting results..."
            )

            if return_db:
                merged_data_path = Path(tempfile.mkdtemp("_wisp_merged_data"))
            merged_data = {}

            for i in range(fasta_count):
                result = results[i]
                for file_id, data in result.items():
                    file_md = self._md[file_id]
                    match len(file_md):
                        case 0:
                            tax_id = "no-tax-id"
                        case 1:
                            tax_id = file_md[0]["TaxId"]
                        case _:
                            tax_id = f"multiple-{'|'.join(sorted(m['TaxId'] if m else 'no-tax-id' for m in file_md))}"

                    if return_db:
                        if tax_id not in merged_data:
                            merged_data[tax_id] = {
                                "counters": Cache(
                                    merged_data_path / str(tax_id) / "counter",
                                    size_limit=sys.maxsize,
                                ),
                                "sources": Cache(
                                    merged_data_path / str(tax_id) / "source",
                                    size_limit=sys.maxsize,
                                ),
                                "last_id": -1,
                            }

                        Cache(
                            merged_data_path / str(tax_id) / "counter",
                            size_limit=sys.maxsize,
                        )

                        current_id = merged_data[tax_id]["last_id"]
                        db_counter = merged_data[tax_id]["counters"]
                        db_source = merged_data[tax_id]["sources"]
                        for counter, source in zip(data["counters"], data["sources"]):
                            current_id += 1
                            db_counter[current_id] = counter
                            db_source[current_id] = source
                        merged_data[tax_id]["last_id"] = current_id
                    else:
                        merged_data.setdefault(tax_id, {"counters": [], "sources": []})
                        merged_data[tax_id]["counters"].extend(data["counters"])
                        merged_data[tax_id]["sources"].extend(data["sources"])

            LOG.debug(
                f"[{archive_path}] {len(merged_data)} different tax_id(s) found - Deleting temporary directory ({temp_dir}) ..."
            )

        LOG.debug(f"[{archive_name}] Temporary files deleted")

        return {
            "archive": archive_path,
            "merged_data": merged_data,
            "tmp_dir": merged_data_path if return_db else None,
        }

    @staticmethod
    def process_fasta(
        file_path: Path | str,
        kmer_sizes: list[int],
        window_size: int,
        step: int,
        full: bool,
        compression: str | None,
    ) -> dict:
        """Simpler Fasta counting method, hopefully less bugged"""
        file_path = Path(file_path).resolve()
        file_name, file_size = file_path.name, format_size(file_path.stat().st_size)

        try:
            sequences = Reader._read_fasta(file_path)
            LOG.debug(
                f"[{file_name}] SIZE: {file_size} - {len(sequences)} sequences found"
            )
        except Exception:
            LOG.exception(f"Error while reading fasta file: {file_path}")
            raise

        source_id_to_data = {}
        for sequence in sequences:
            kmer_counts = Reader._counter(
                entry=sequence["sequence"],
                kmer_sizes=kmer_sizes,
                window_size=window_size,
                step=step,
                full=full,
            )
            # rearange counters
            num_elements = len(next(iter(kmer_counts.values())))
            kmer_counts_as_list = []
            for i in range(num_elements):
                new_dict = {key: values[i] for key, values in kmer_counts.items()}
                kmer_counts_as_list.append(new_dict)

            source = {"id": sequence["id"], "contig": sequence["contig"]}
            source_id_to_data.setdefault(source["id"], {"counters": [], "sources": []})
            if compression:
                kmer_counts_as_list = compress(
                    kmer_counts_as_list, as_list=True, format=compression
                )
                source = compress(source, format=compression)
            source_id_to_data[sequence["id"]]["counters"].extend(kmer_counts_as_list)
            source_id_to_data[sequence["id"]]["sources"].extend(
                [source] * len(kmer_counts_as_list)
            )

        LOG.debug(f"[{file_name}] Return {len(sequences)} counters & sources")
        return source_id_to_data

    @staticmethod
    def _counter(
        entry: str,
        kmer_sizes: list[int],
        window_size: int,
        step: int,
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

        kmer_counters = {}
        for kmer_size in kmer_sizes:
            window_counters = []
            for start, end in windows:
                kmers = (
                    entry[i : i + kmer_size] for i in range(start, end - kmer_size + 1)
                )

                counts = Counter(kmers)
                rev_counts = Counter(
                    {
                        Reader._revcomp(k, compl=complements): v
                        for k, v in counts.items()
                    }
                )
                counts += rev_counts

                if False:  # TODO: conf=max_size
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
            kmer_counters[kmer_size] = window_counters

        return kmer_counters

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
                if "." in record.id:  # .contig or .something_else
                    r_id, r_contig = record.id.split(".")[:2]
                else:
                    r_id, r_contig = record.id, ""
                sequences.append(
                    {
                        "id": r_id,
                        "contig": r_contig,
                        "sequence": str(record.seq),
                    }
                )
        return sequences


if __name__ == "__main__":
    pass
