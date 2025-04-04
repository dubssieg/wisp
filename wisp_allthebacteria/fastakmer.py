import concurrent.futures
import gc
import logging
import os
import shutil
import sys
import tarfile
import tempfile
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
import time
from typing import Generator
from diskcache import Cache

from Bio import SeqIO
import concurrent
from loky import get_reusable_executor
from metadata import Metadata
from utils import format_size, cleanup_zombie_processes, compress, space_format
from taxdb import TaxDB

LOG = logging.getLogger(__name__)

NO_TAX_ID = "no-tax-id"


class FastaKmer:
    def __init__(
        self,
        metadata: Metadata,
        num_workers: int,
        taxdb: TaxDB,
        tmp_path: str | Path = "/tmp",
    ):
        LOG.debug(f"FastaKmer({locals()})")
        self._md = metadata
        self._num_workers = num_workers
        self._taxdb = taxdb
        self._tmp_path = Path(tmp_path).resolve()

    def process_archive(
        self,
        archive_path: Path | str,
        kmer_sizes: list[int],
        window_size: int,
        step: int,
        full: bool = False,
        batch_size: int = None,
        compression: str | None = None,
        merged_data_as_db: bool = True,
        current_counts: dict = {},
        max_count: int | None = None,
    ) -> Generator[dict, None, None]:
        """Extract and process an archive."""
        archive_path = Path(archive_path).resolve()
        archive_name = archive_path.name

        if self._tmp_path.exists():
            LOG.warning(f"Deleting {self._tmp_path}")
            shutil.rmtree(self._tmp_path)
        self._tmp_path.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory(
            dir=self._tmp_path, suffix="_wisp_fasta"
        ) as tmp_dir:
            tmp_dir = Path(tmp_dir).resolve()
            archive_size = format_size(archive_path.stat().st_size)
            LOG.info(f"[{archive_name}] Extracting archive")
            LOG.info(f"[{archive_name}] SIZE: {archive_size}, temp dir: {tmp_dir}")
            # extraction de l'archive
            with tarfile.open(archive_path, "r:xz") as tar:
                tar.extractall(tmp_dir)
                extracted_files = list(tmp_dir.rglob("*.fa"))
                LOG.debug(
                    f"[{archive_name}] {len(extracted_files)} FASTA files extracted"
                )

            file_selector = None
            if max_count:
                dry_results = next(
                    self._process_archive_workers(
                        archive_name=archive_name,
                        tmp_dir=tmp_dir,
                        batch_size=None,
                        kmer_sizes=kmer_sizes,
                        window_size=window_size,
                        step=step,
                        full=full,
                        compression=False,
                        dry=True,
                        merged_data_as_db=False,
                        file_selector=None,
                    )
                )

                file_selector = self._count_analysis(
                    dry_results=dry_results,
                    max_counts=max_count,
                    current_counts=current_counts,
                )

            yield from self._process_archive_workers(
                archive_name=archive_name,
                tmp_dir=tmp_dir,
                batch_size=batch_size,
                kmer_sizes=kmer_sizes,
                window_size=window_size,
                step=step,
                full=full,
                compression=compression,
                dry=False,
                merged_data_as_db=merged_data_as_db,
                file_selector=file_selector,
            )

            LOG.debug(f"[{archive_path}] Deleting temporary directory ({tmp_dir}) ...")
        LOG.debug(f"[{archive_name}] Temporary files deleted")

    def _process_archive_workers(
        self,
        archive_name: str,
        tmp_dir: Path,
        batch_size: int | None,
        kmer_sizes: list[int],
        window_size: int,
        step: int,
        full: bool,
        compression: str,
        dry: bool,
        merged_data_as_db: bool,
        file_selector: dict | None = None,
    ) -> Generator[dict, None, None]:
        dry_str = "[==DRY==]" if dry else ""

        results = []
        # if dry:
        #     results = [None] * len(extracted_files)
        # else:
        #     results = Cache(tmp_dir / f"results{dry_suffix}", size_limit=sys.maxsize)
        fasta_count = 0

        extracted_files = list(tmp_dir.rglob("*.fa"))
        orig_file_count = len(extracted_files)
        if file_selector is not None:
            extracted_files = [
                p for p in extracted_files if p.name in file_selector.keys()
            ]

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
                futures = {
                    executor.submit(
                        FastaKmer._process_fasta,
                        file_path=file_path,
                        kmer_sizes=kmer_sizes,
                        window_size=window_size,
                        step=step,
                        full=full,
                        compression=compression,
                        dry=dry,
                        file_selector=file_selector,
                    ): file_path
                    for file_path in batch
                }
                try:
                    for future in concurrent.futures.as_completed(futures):
                        try:
                            result = future.result(timeout=60)  # TODO: conf

                            if dry:
                                results.append(result)
                            else:
                                yield self._process_result(
                                    result=result, merged_data_as_db=merged_data_as_db
                                )
                            fasta_count += 1

                            if not dry:
                                LOG.debug(
                                    f"{dry_str}[{archive_name}] {futures[future].name}: {fasta_count} / {len(extracted_files)} ({orig_file_count})"
                                )

                        except concurrent.futures.TimeoutError:
                            LOG.error(f"[{archive_name}] Timeout")
                        except Exception:
                            LOG.exception(f"[{archive_name}]")
                            raise
                    LOG.debug(
                        f"{dry_str}[{archive_name}] Batch {batch_i + 1} / {len(batches)} done. Terminating worker pool..."
                    )
                finally:
                    LOG.debug(f"{dry_str}[{archive_name}] Shutting down executor...")
                    # force workers to stop
                    executor.shutdown(wait=True, kill_workers=True)
                    LOG.debug(f"{dry_str}[{archive_name}] Executor shut down.")

                    # give time for zombies to appear
                    time.sleep(2)
                    LOG.debug(
                        f"{dry_str}[{archive_name}] Checking for zombie processes..."
                    )
                    cleanup_zombie_processes(os.getpid())
                    LOG.debug(f"{dry_str}[{archive_name}] zombie cleanup done")

                    # free mem
                    gc.collect()
                    LOG.debug(
                        f"{dry_str}[{archive_name}] Cleanup complete - Workers pool terminated"
                    )

            LOG.debug(f"{dry_str}[{archive_name}] Workers pool terminated")
        LOG.debug(
            f"{dry_str}[{archive_name}] All {len(extracted_files)} Fasta files processed"
        )
        yield results if dry else None

    @staticmethod
    def _process_fasta(
        file_path: Path | str,
        kmer_sizes: list[int],
        window_size: int,
        step: int,
        full: bool,
        compression: str | None,
        dry: bool = False,
        file_selector: dict | None = None,
    ) -> dict:
        """Fasta counting method"""
        file_path = Path(file_path).resolve()
        file_name, file_size = file_path.name, format_size(file_path.stat().st_size)

        try:
            sequences = FastaKmer._read_fasta(file_path)
            LOG.debug(
                f"[{file_name}] SIZE: {file_size} - {len(sequences)} sequences found"
            )
        except Exception:
            LOG.exception(f"Error while reading fasta file: {file_path}")
            raise

        has_limits = file_selector is not None
        # if file_selector:
        #     file_id_selector =

        if has_limits:
            file_id_limits = file_selector[file_name]

        source_id_to_data = defaultdict(list)
        for sequence in sequences:
            if has_limits:
                file_id = sequence["id"]
                if file_id not in file_id_limits:
                    continue

            if len(sequence["sequence"]) < window_size:
                continue

            kmer_counts = FastaKmer._counter(
                entry=sequence["sequence"],
                kmer_sizes=kmer_sizes,
                window_size=window_size,
                step=step,
                full=full,
                dry=dry,
            )
            # rearange counters
            num_elements = len(next(iter(kmer_counts.values())))

            if has_limits:
                if num_elements > file_id_limits[file_id]:
                    num_elements = file_id_limits[file_id]
                file_id_limits[file_id] -= num_elements
                if file_id_limits[file_id] == 0:
                    del file_id_limits[file_id]

            kmer_counts_as_list = []
            for i in range(num_elements):
                kmer_counts_i = {
                    "counter": {key: values[i] for key, values in kmer_counts.items()},
                    "md": {"id": sequence["id"], "contig": sequence["contig"]},
                }
                if dry:
                    kmer_counts_i["md"]["file_name"] = file_name
                kmer_counts_as_list.append(kmer_counts_i)

            if compression:
                kmer_counts_as_list = compress(
                    kmer_counts_as_list, as_list=True, format=compression
                )
            source_id_to_data[sequence["id"]].extend(kmer_counts_as_list)

            if has_limits:
                if len(file_id_limits) == 0:
                    break

        LOG.debug(
            f"[{file_name}] Return {sum(len(counters) for counters in source_id_to_data.values())} counters ({len(source_id_to_data)} sequence ids)"
        )
        return dict(source_id_to_data)

    def _count_analysis(
        self, dry_results: list, max_counts: int, current_counts: dict
    ) -> dict:
        file_selector = defaultdict(dict)
        tax_id_to_file_ids = defaultdict(lambda: defaultdict(int))

        for result in dry_results:
            for file_id, counter_list in result.items():
                tax_id = self._taxdb.clean_tax_id(self._md[file_id].get("TaxId"))
                for counter in counter_list:
                    tax_id_to_file_ids[tax_id][
                        (file_id, counter["md"]["file_name"])
                    ] += 1

        for tax_id, file_id_counts in tax_id_to_file_ids.items():
            remaining = max_counts - current_counts.get(tax_id, 0)
            if remaining <= 0:
                continue

            for (file_id, file_name), count in sorted(
                file_id_counts.items(), key=lambda x: -x[1]
            ):
                max_insertable = min(remaining, count)
                if max_insertable > 0:
                    file_selector[file_name][file_id] = max_insertable
                    remaining -= max_insertable
                if remaining <= 0:
                    break

        selected = sum(sum(d.values()) for d in file_selector.values())
        total = sum(sum(d.values()) for d in tax_id_to_file_ids.values())
        LOG.debug("Sample analysis done:")
        LOG.debug(f"- Files remaining: {len(file_selector)}")
        LOG.debug(f"- Total samples: {space_format(total)}")
        LOG.debug(f"- Selected: {space_format(selected)}")
        LOG.debug(f"- Skipped:  {space_format(total - selected)}")
        LOG.debug(
            f"- Tax IDs:  {', '.join(map(str, sorted(tax_id_to_file_ids.keys())))}"
        )

        return dict(file_selector)

    def _process_result(self, result: dict, merged_data_as_db: bool = False) -> dict:
        if merged_data_as_db:
            merged_data_path = Path(
                tempfile.mkdtemp(dir=self._tmp_path, suffix="_wisp_merged_data")
            )
        merged_data = {}
        for file_id, counters in result.items():
            md = self._md[file_id]
            tax_id = self._taxdb.clean_tax_id(md.get("TaxId", None))

            if merged_data_as_db:
                if tax_id not in merged_data:
                    merged_data[tax_id] = {
                        "db": Cache(
                            merged_data_path / str(tax_id) / "counter",
                            size_limit=sys.maxsize,
                        ),
                        "last_id": -1,
                    }

                current_id = merged_data[tax_id]["last_id"]
                db = merged_data[tax_id]["db"]
                for counter in counters:
                    current_id += 1
                    db[current_id] = counter
                merged_data[tax_id]["last_id"] = current_id
            else:
                if tax_id not in merged_data:
                    merged_data[tax_id] = []
                merged_data[tax_id].extend(counters)
        return {
            "merged_data": merged_data,
            "tmp_dir": merged_data_path if merged_data_as_db else None,
        }

    @staticmethod
    def _counter(
        entry: str,
        kmer_sizes: list[int],
        window_size: int,
        step: int,
        full: bool = False,
        dry: bool = False,
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

        if dry:
            return {0: windows}

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
                        FastaKmer._revcomp(k, compl=complements): v
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
