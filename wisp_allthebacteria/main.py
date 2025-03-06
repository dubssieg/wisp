import json
import time
import traceback
from pathlib import Path


from tqdm.auto import tqdm
from metadata import Metadata
from reader import Reader
from writer import Writer


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
    num_window: int,
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

    md = Metadata(Path(METADATA_PATH) / METADATA_FILENAME)
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
                num_window=num_window,
            )
            writer.save_data(merged_data)

            error_log["complete"].append(str(archive_path))
            if str(archive_path) in error_log["incomplete"]:
                error_log["incomplete"].remove(str(archive_path))

            # Calculer et enregistrer la durée
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
    # md = Metadata()
    # md.md
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

    # md = Metadata()

    # md_sup_3 = ["SAMN00189190", "SAMN00189191", "SAMN00189193"]
    # md_none = "SAMEA3400865"
    # md_1 = "SAMD00020420"

    # md[md_1]
    # md[md_sup_3[0]]

    # api = API()
    # print(api.clean_cache())

    create_db(
        ASSEMBLY_PATH,
        "/data/microtaxo/db_full_4",
        kmer_size=4,
        window_size=100,
        num_window=10,
    )
