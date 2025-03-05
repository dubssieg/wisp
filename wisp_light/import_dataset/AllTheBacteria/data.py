import pandas as pd
from pathlib import Path
from urllib.error import HTTPError
from Bio import Entrez, SeqIO
from diskcache import Cache
from tqdm.auto import tqdm
import tempfile
import tarfile


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
        df = pd.read_csv(content, sep="\t")
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

    def process_file(
        self, filename: str
    ) -> dict:  # attention, pathlib assembly, etc. à revoir
        file_path = self._assembly_path / filename
        suffix = file_path.suffix
        if suffix == ".xz":
            return self.process_archive(file_path)
        elif suffix == ".fa":
            return self.process_fasta(file_path)

    def read_fasta(self, file_path: Path | str) -> list:
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

    def process_fasta(self, file_path: str | Path) -> dict:
        sequences = self.read_fasta(file_path)
        processed_data = []
        tax_id_to_seq = {}

        for sequence in sequences:
            md = self._md[sequence["id"]]

            tax_id = md["TaxId"]

            if tax_id not in tax_id_to_seq:
                tax_id_to_seq[tax_id] = {
                    "description": [],
                    "metadata": md,
                    "sequence": "",
                }

            tax_id_to_seq[tax_id]["description"].append(sequence["description"])
            tax_id_to_seq[tax_id]["sequence"] += sequence["sequence"]

        processed_data = list(tax_id_to_seq.values())

        return processed_data

    def _merge_processed_files(self, file_paths: list[Path]) -> dict[int, dict]:
        """Multiple fasta files.
        Concatenate by tax_id"""
        merged_data = {}

        for file_path in tqdm(file_paths, desc="Process and merge files"):
            processed_data = self.process_fasta(file_path)

            for entry in processed_data:
                tax_id = int(entry["metadata"]["TaxId"])
                if tax_id not in merged_data:
                    merged_data[tax_id] = {
                        "description": [],
                        "metadata": entry["metadata"],
                        "sequence": "",
                    }

                merged_data[tax_id]["description"].extend(entry["description"])
                merged_data[tax_id]["sequence"] += entry["sequence"]

        return merged_data

    def process_archive(self, archive_path: Path | str):
        """Read and process a compressed archive."""
        with tempfile.TemporaryDirectory() as temp_dir:
            with tarfile.open(archive_path, "r:xz") as tar:
                tar.extractall(path=temp_dir)

            extracted_files = list(Path(temp_dir).rglob("*.fa"))

            return self._merge_processed_files(extracted_files)


if __name__ == "__main__":
    md = Metadata()
    # md["SAMD00013333.contig0000"]
    # md.tax_ids()
    # t, e, err = md.populate_api_cache()
    # print(md[367830])
    reader = Reader(md)

    # reader.process_file("achromobacter_xylosoxidans__01/SAMD00013333.fa")
    reader.process_file("actinobacillus_lignieresii__01.asm.tar.xz")
    pass
