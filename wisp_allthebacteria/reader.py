from collections import Counter
from itertools import product
from pathlib import Path
import tarfile
import tempfile
import traceback
from Bio import SeqIO
from tqdm.auto import tqdm
from metadata import Metadata


class Reader:
    def __init__(self, metadata: Metadata):
        self._md = metadata

    def process_file(
        self,
        file_path: str | Path,
        kmer_size: int,
        window_size: int,
        step: int,
        full: bool = False,
    ) -> dict:
        file_path = Path(file_path)
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

    def read_fasta(self, file_path: Path | str) -> list:
        """Just read and parse file, no processing.
        Return a list of:
            'id' = 'SAMD00013333'
            'contig' = 'contig00001'
            'file' = 'SAMEA1561896.fa'
            'sequence' = 'GGAGGGAACAGCGGGGCGGGCGGCGT..."""
        file_path = Path(file_path)
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

    def process_fasta(
        self,
        file_path: str | Path,
        kmer_size,
        window_size: int,
        step: int,
        full: bool = False,
    ) -> dict:
        """Count and get metadata"""
        file_path = Path(file_path)
        sequences = self.read_fasta(file_path)
        tax_id_to_data = {}

        for sequence in tqdm(
            sequences, desc=f"Counting {file_path.name}", leave=False, position=2
        ):
            md = self._md[sequence["id"]]
            match len(md):
                case 0:
                    tax_id = "no-tax-id"
                case 1:
                    tax_id = md[0]["TaxId"]
                case _:
                    tax_id = "|".join([(m["TaxId"] if m else "no-tax-id") for m in md])
                    tax_id = f"multiple-{tax_id}"

            kmer_count = self._counter(
                entry=sequence["sequence"],
                kmer_size=kmer_size,
                window_size=window_size,
                step=step,
                full=full,
            )

            if tax_id not in tax_id_to_data:
                tax_id_to_data[tax_id] = {"metadata": md, "counters": [], "sources": []}

            source = {
                "id": sequence["id"],
                "contig": sequence["contig"],
                "file": sequence["file"],
                "win": 0,
            }
            source = {k: v for k, v in sequence.items() if k != "sequence"}
            tax_id_to_data[tax_id]["counters"].extend(kmer_count)
            tax_id_to_data[tax_id]["sources"].append(source)

        return tax_id_to_data

    def _counter(
        self,
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
                {self._revcomp(k, compl=complements): v for k, v in counts.items()}
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
        except KeyError as exc:
            traceback.print_exc()
            raise IndexError(
                "Complementarity does not include all chars in sequence."
            ) from exc
        return result

    def process_archive(
        self,
        archive_path: Path | str,
        kmer_size: int,
        window_size: int,
        step: int,
        full: bool = False,
    ):
        """Read and process an archive."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with tarfile.open(archive_path, "r:xz") as tar:
                tar.extractall(path=temp_dir)

            extracted_files = list(Path(temp_dir).rglob("*.fa"))
            merged_data = {}

            for file_path in tqdm(
                extracted_files,
                desc=f"Processing {Path(archive_path).name}",
                leave=False,
                position=1,
            ):
                file_data = self.process_fasta(
                    file_path,
                    kmer_size=kmer_size,
                    window_size=window_size,
                    step=step,
                    full=full,
                )

                for tax_id, data in file_data.items():
                    if tax_id not in merged_data:
                        merged_data[tax_id] = {
                            "metadata": data["metadata"],
                            "counters": list(),
                            "sources": list(),
                        }

                    merged_data[tax_id]["counters"].extend(data["counters"])
                    merged_data[tax_id]["sources"].extend(data["sources"])

            return {"archive": archive_path, "merged_data": merged_data}
