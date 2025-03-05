import json
import time
import pandas as pd
from pathlib import Path
from urllib.error import HTTPError
from Bio import Entrez, SeqIO
from diskcache import Cache
from tqdm.auto import tqdm
from collections import Counter
import tempfile
import tarfile
from itertools import product


METADATA_PATH = "/data/microtaxo/allthebacteria_sample/metadata"
ASSEMBLY_PATH = "/data/microtaxo/allthebacteria_sample/assembly"
METADATA_FILENAME = "ena_metadata.tsv"
API_CACHE_DIR = "/data/microtaxo/apicache"
EMAIL = "cyrille.leroux@irisa.fr"


class Metadata:
    def __init__(self):
        self._md = None
        self._api = API()

    @property
    def md(self):
        if self._md is None:
            self._md = self._load()
            # self._replace_tax_id()
        return self._md

    def __getitem__(self, seq_id: str) -> dict:
        """Example: SAMD00013333 or SAMD00013333.contig0000"""
        if ".contig" in seq_id:
            seq_id = seq_id.split(".contig")[0]
        data = self.md[seq_id]
        if not isinstance(data, dict):
            self._md[seq_id] = data = self._api[data]
        return data

    def _load(self):
        """Read DataFrame once then free memory (big file)."""
        # read data : key=seq_id from file headers, value=tax_id
        content = Path(METADATA_PATH) / METADATA_FILENAME
        df = pd.read_csv(content, sep="\t", low_memory=False)
        return {
            seq_id: API.clean_tax_id(tax_id)
            for seq_id, tax_id in zip(df["sample_accession"], df["tax_id"])
        }


class API:
    def __init__(self):
        self._api_cache_dir = API_CACHE_DIR
        self._api_cache = Cache(self._api_cache_dir)
        self._df = None
        self._tax_id_errors_ = set()
        Entrez.email = EMAIL

    @property
    def df(self):
        if self._df is None:
            content = Path(METADATA_PATH) / METADATA_FILENAME
            self._df = pd.read_csv(content, sep="\t")
        return self._df

    def __getitem__(self, tax_id):
        if record := self._api_cache.get(tax_id, None):
            return record[0]
        return record

    def tax_ids(self) -> list:
        """ "List of all available tax_ids"""
        return list(self._api_cache)

    def populate_api_cache(self):
        """Run once."""
        errors = set()
        unique_tax_ids = list()
        extra_tax_ids = set()
        # direct tax_id
        for tax_id in tqdm(self.df["tax_id"].unique()):
            if self._get_api_data(tax_id):
                unique_tax_ids.append(tax_id)
            else:
                errors.add(tax_id)

        # extra tax_id
        for tax_id in tqdm(unique_tax_ids):
            if tax_id in self._api_cache:
                if tdata := self._api_cache[tax_id]:
                    extra_tax_ids.update(
                        [int(lin["TaxId"]) for lin in tdata[0].get("LineageEx", {})]
                    )
        for tax_id in tqdm(extra_tax_ids):
            self._get_api_data(tax_id)

        return unique_tax_ids, list(extra_tax_ids), list(errors)

    @staticmethod
    def clean_tax_id(tax_id):
        # error: not a number / None
        if pd.isna(tax_id):
            return None

        # float/string -> int
        try:
            tax_id = int(tax_id)
        except ValueError:
            return None
        return tax_id

    def _get_api_data(self, tax_id):
        tax_id = self.clean_tax_id(tax_id)

        # hit cache
        if tax_id in self._cache:
            return self._cache[tax_id]

        # API call + cache
        try:
            handle = Entrez.efetch(db="taxonomy", id=str(tax_id), retmode="xml")
            records = Entrez.read(handle)
            self._cache[tax_id] = records
            return records
        except HTTPError:
            return None


class Reader:
    def __init__(self, metadata: Metadata):
        self._md = metadata
        self._assembly_path = Path(ASSEMBLY_PATH)

    def process_file(self, filename: str, kmer_size: int) -> dict:
        file_path = self._assembly_path / filename
        suffix = file_path.suffix
        if suffix == ".xz":
            return self.process_archive(file_path, kmer_size=kmer_size)
        elif suffix == ".fa":
            return self.process_fasta(file_path, kmer_size=kmer_size)

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

    def process_fasta(self, file_path: str | Path, kmer_size) -> dict:
        """Count and get metadata"""
        file_path = Path(file_path)
        sequences = self.read_fasta(file_path)
        tax_id_to_data = {}

        for sequence in tqdm(
            sequences, desc=f"Counting {file_path.name}", leave=False, position=2
        ):
            md = self._md[sequence["id"]]
            tax_id = md["TaxId"]

            kmer_count = self._counter(sequence["sequence"], kmer_size)

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

    def _counter(self, entry: str, kmer_size: int = 4) -> dict:
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

        all_kmers = (
            entry[i : i + kmer_size] for i in range(len(entry) - kmer_size + 1)
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
            raise IndexError(
                "Complementarity does not include all chars in sequence."
            ) from exc
        return result

    def process_archive(self, archive_path: Path | str, kmer_size: int = 4):
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
                file_data = self.process_fasta(file_path, kmer_size)

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


class Writer:
    def __init__(self, path: Path | str):
        self._path = Path(path)

    def save_data(self, data: dict):
        """Saves the processed data into a structured directory."""

        self._path.mkdir(parents=True, exist_ok=True)

        archive = Path(data["archive"]).stem.split(".")[0]
        merged_data = data["merged_data"]
        for tax_id, data in tqdm(
            merged_data.items(), desc="write merged data", position=1, leave=False
        ):
            tax_path = self._path / str(tax_id)
            tax_path.mkdir(parents=True, exist_ok=True)

            # Save metadata
            metadata_path = tax_path / "metadata.json"
            if not metadata_path.exists():
                with open(metadata_path, "w", encoding="utf-8") as md_file:
                    json.dump(data["metadata"], md_file, indent=4)

            # Save headers TODO: duplicates if errors?
            headers_path = tax_path / "headers"
            headers_path.mkdir(parents=True, exist_ok=True)

            file_name = Path(archive)
            header_file_path = headers_path / f"{file_name}.txt"
            with open(header_file_path, "a", encoding="utf-8") as header_file:
                for source in tqdm(
                    data["sources"], desc="Save sources", position=2, leave=False
                ):
                    header_file.write(source["description"] + "\n")

            # Save k-mer counts in libSVM format
            libsvm_file_path = tax_path / f"{archive}.libsvm"
            with open(libsvm_file_path, "w", encoding="utf-8") as libsvm_file:
                for source in tqdm(
                    data["sources"],
                    desc="Save k-mer in libSVM format",
                    position=2,
                    leave=False,
                ):

                    for counter in data["counters"]:
                        libsvm_file.write(
                            f"{tax_id} {' '.join([f'{k}:{v}' for k, v in counter.items()])}\n"
                        )


def format_duration(seconds):
    """Convert seconds to a formatted duration string (hh:mm:ss)."""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"


def create_db(input_path: Path | str, output_path: Path, kmer_size: int = 4):
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    json_log_path = output_path / "error_log.json"

    if json_log_path.exists():
        with open(json_log_path, "r", encoding="utf-8") as error_log_file:
            error_log = json.load(error_log_file)
    else:
        error_log = {"complete": [], "incomplete": [], "duration": {}}

    archives = list(input_path.glob("*.xz"))

    md = Metadata()
    reader = Reader(md)
    writer = Writer(output_path)

    for archive_path in tqdm(
        archives, desc=f"Processing {str(input_path)} -> {str(output_path)}", position=0
    ):
        if str(archive_path) in error_log["complete"]:
            continue

        start_time = time.time()

        try:
            if str(archive_path) in error_log["incomplete"]:
                print(f"Retrying {archive_path.name}...")

            merged_data = reader.process_file(archive_path, kmer_size)
            writer.save_data(merged_data)

            error_log["complete"].append(str(archive_path))
            if str(archive_path) in error_log["incomplete"]:
                error_log["incomplete"].remove(str(archive_path))

            # Calculer et enregistrer la durée
            duration = time.time() - start_time
            error_log["duration"][archive_path.name] = format_duration(duration)

        except Exception as e:
            print(f"Error processing {archive_path.name}: {e}")
            if str(archive_path) not in error_log["incomplete"]:
                error_log["incomplete"].append(str(archive_path))

        finally:
            with open(json_log_path, "w", encoding="utf-8") as error_log_file:
                json.dump(error_log, error_log_file, indent=4)


# def encode_dna_to_binary(sequence):
#     nucleotide_to_binary = {"A": 0b00, "T": 0b01, "C": 0b10, "G": 0b11}

#     binary_sequence = bytearray()
#     current_byte = 0
#     bits_added = 0

#     for nucleotide in sequence:
#         current_byte = (current_byte << 2) | nucleotide_to_binary[nucleotide]
#         bits_added += 2

#         if bits_added == 8:
#             binary_sequence.append(current_byte)
#             current_byte = 0
#             bits_added = 0

#     if bits_added > 0:
#         current_byte <<= 8 - bits_added
#         binary_sequence.append(current_byte)

#     return binary_sequence


if __name__ == "__main__":
    md = Metadata()
    # md["SAMD00013333.contig0000"]
    # md.tax_ids()
    # t, e, err = md.populate_api_cache()
    # print(md[367830])
    #
    # content = reader.process_fasta(
    #     Path(ASSEMBLY_PATH) / "achromobacter_xylosoxidans__01/SAMD00013333.fa",
    #     kmer_size=4,
    # )
    # sequence = content[0]["sequence"]

    # reader.process_file("achromobacter_xylosoxidans__01/SAMD00013333.fa")

    # reader = Reader(md)
    # reader.process_file("actinobacillus_lignieresii__01.asm.tar.xz")

    # import pickle

    # with open("/data/microtaxo/merged_data.pkl", "rb") as file:
    #     loaded_data = pickle.load(file)
    #     loaded_data = {
    #         "archive": "actinobacillus_lignieresii__01.asm.tar.xz",
    #         "merged_data": loaded_data,
    #     }
    # writer = Writer("/tmp/microdb1")
    # writer.save_processed_data(loaded_data)

    create_db(ASSEMBLY_PATH, "/data/microtaxo/db_full_4", kmer_size=4)
