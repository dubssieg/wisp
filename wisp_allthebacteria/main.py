import json
import time
import traceback
from pathlib import Path


from tqdm.auto import tqdm
from metadata import Metadata
from reader import Reader
from writer import Writer
from api import API


METADATA_PATH = "/data/microtaxo/allthebacteria_sample/metadata"
ASSEMBLY_PATH = "/data/microtaxo/allthebacteria_sample/assembly"

METADATA_FILENAME = "ena_metadata.tsv"


def format_duration(seconds):
    """Convert seconds to a formatted duration string (hh:mm:ss)."""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"


def create_db(
    input_path: Path | str,
    output_path: Path,
    kmer_size: int,
    window_size: int,
    num_windows: int,
):
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

    api = API()
    md = Metadata(csv_path=Path(METADATA_PATH) / METADATA_FILENAME, api=api)
    md.md  # preload
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

            merged_data = reader.process_file(
                archive_path,
                kmer_size=kmer_size,
                window_size=window_size,
                num_windows=num_windows,
            )
            writer.save_data(merged_data)

            error_log["complete"].append(str(archive_path))
            if str(archive_path) in error_log["incomplete"]:
                error_log["incomplete"].remove(str(archive_path))

            duration = time.time() - start_time
            error_log["duration"][archive_path.name] = format_duration(duration)

        except Exception as e:
            print(f"Error processing {archive_path.name}: {e}")
            traceback.print_exc()
            if str(archive_path) not in error_log["incomplete"]:
                error_log["incomplete"].append(str(archive_path))

        finally:
            with open(json_log_path, "w", encoding="utf-8") as error_log_file:
                json.dump(error_log, error_log_file, indent=4)


if __name__ == "__main__":
    api = API()
    # from database import Database

    # db = Database("/data/microtaxo/db_full_4", api)
    # print(db.get_tax_ids_by_rank("phylum"))
    # pass

    import pickle

    writer = Writer(path="/data/microtaxo/db_full_4")
    with open("/data/microtaxo/merged_data.pkl", "rb") as file:
        merged_data = pickle.load(file)
    writer.save_data(merged_data)

    # create_db(
    #     ASSEMBLY_PATH,
    #     "/data/microtaxo/db_full_4",
    #     kmer_size=4,
    #     window_size=100,
    #     num_windows=10,
    # )
