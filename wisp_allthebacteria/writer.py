import json
from pathlib import Path
import pickle

from tqdm.auto import tqdm


class Writer:
    def __init__(self, path: Path | str):
        self._path = Path(path)

    def save_data(self, data: dict, as_libsvm: bool = False):
        """Saves the processed data into a structured directory."""

        self._path.mkdir(parents=True, exist_ok=True)

        archive = Path(data["archive"]).stem.split(".")[0]
        merged_data = data["merged_data"]
        # create update files for each tax_id
        for tax_id, tdata in tqdm(
            merged_data.items(), desc="write merged data", position=1, leave=False
        ):
            tax_path = self._path / str(tax_id)
            tax_path.mkdir(parents=True, exist_ok=True)

            # Save metadata
            metadata_path = tax_path / "metadata.json"
            if not metadata_path.exists():
                with open(metadata_path, "w", encoding="utf-8") as md_file:
                    json.dump(tdata["metadata"], md_file, indent=4)

            if as_libsvm:
                # Save headers TODO: duplicates if errors?
                headers_path = tax_path / "headers"
                headers_path.mkdir(parents=True, exist_ok=True)

                file_name = Path(archive)
                header_file_path = headers_path / f"{file_name}.txt"
                with open(header_file_path, "a", encoding="utf-8") as header_file:
                    for source in tqdm(
                        tdata["sources"], desc="Save sources", position=2, leave=False
                    ):
                        header_file.write(source["description"] + "\n")

                # Save k-mer counts in libSVM format
                libsvm_file_path = tax_path / f"{archive}.libsvm"
                lines_to_write = []
                for counter in tdata["counters"]:
                    lines_to_write.append(
                        f"{tax_id} {' '.join([f'{k}:{v}' for k, v in counter.items()])}\n"
                    )

                with open(libsvm_file_path, "w", encoding="utf-8") as libsvm_file:
                    libsvm_file.writelines(lines_to_write)
            else:
                pickle_file_path = tax_path / f"{archive}.pkl"
                with open(pickle_file_path, "wb") as pickle_file:
                    pickle.dump(tdata, pickle_file)
