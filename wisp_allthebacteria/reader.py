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
        self, file_path: str | Path, kmer_size: int, window_size: int, num_windows: int
    ) -> dict:
        suffix = file_path.suffix
        if suffix == ".xz":
            return self.process_archive(
                file_path,
                kmer_size=kmer_size,
                window_size=window_size,
                num_windows=num_windows,
            )
        elif suffix == ".fa":
            return self.process_fasta(
                file_path,
                kmer_size=kmer_size,
                window_size=window_size,
                num_windows=num_windows,
            )

    def read_fasta(self, file_path: Path | str) -> list:
        """Just read and parse file, no processing.
        Return a list of:
            'id' = 'SAMD00013333.contig00001'
            'description' = 'SAMD00013333.contig00001 len=378640 cov=42.4 ...
            'sequence' = 'GGAGGGAACAGCGGGGCGGGCGGCGT..."""
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"File {file_path} not found.")

        sequences = []
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences.append(
                    {
                        "id": record.id,
                        "description": record.description,
                        "sequence": str(record.seq),
                    }
                )
        return sequences

    def process_fasta(
        self, file_path: str | Path, kmer_size, window_size: int, num_windows: int
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
                num_windows=num_windows,
            )

            if tax_id not in tax_id_to_data:
                tax_id_to_data[tax_id] = {"metadata": md, "counters": [], "sources": []}

            tax_id_to_data[tax_id]["counters"].append(kmer_count)
            tax_id_to_data[tax_id]["sources"].append(
                {
                    "description": sequence["description"],
                    "file": file_path.name,
                }
            )

        return tax_id_to_data

    def _counter(
        self,
        entry: str,
        kmer_size: int = 4,
        window_size: int = 1000,
        num_windows: int = 100,
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

        if len(entry) < window_size * num_windows:
            all_kmers = (
                entry[i : i + kmer_size] for i in range(len(entry) - kmer_size + 1)
            )
        else:
            step = max(1, (len(entry) - window_size) // (num_windows - 1))
            positions = range(0, len(entry) - window_size + 1, step)

            all_kmers = (
                entry[i + j : i + j + kmer_size]
                for i in positions
                for j in range(window_size - kmer_size + 1)
            )

        counts = Counter(all_kmers)
        rev_counts = Counter(
            {self._revcomp(k, compl=complements): v for k, v in counts.items()}
        )
        counts += rev_counts

        for filtered_kmer in (alpha * kmer_size for alpha in "ATCG"):
            counts.pop(filtered_kmer, None)

        counts_purged = {}
        for key, count in counts.items():
            list_of_keys = [degenerate_map.get(x, [x]) for x in key]
            list_of_keys = ["".join(item) for item in product(*list_of_keys)]

            for prob_key in list_of_keys:
                kmer_number = count // len(list_of_keys)
                counts_purged[prob_key] = counts_purged.get(prob_key, 0) + kmer_number

        return counts_purged

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
        num_windows: int,
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
                    num_windows=num_windows,
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
