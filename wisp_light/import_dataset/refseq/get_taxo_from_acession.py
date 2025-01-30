import os

import pandas as pd
import yaml
from Bio import SeqIO, Entrez

"""
> wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -o {output_dir}
voir version info 
enlevé la 1ere ligne et le premier # , enregister en .csv

"""
def filter_assembly_summary(csv_file_in, csv_file_out):
    info_all_genome_df = pd.read_csv(csv_file_in, sep="\t",low_memory=False,
                                     usecols=["assembly_accession", "assembly_level", "ftp_path"])
    print(info_all_genome_df[:10])
    # Filtrer les données
    df_filtered = info_all_genome_df[
        (info_all_genome_df["assembly_level"] == "Complete Genome") &  # Filtre "Complete Genome"
        (info_all_genome_df["ftp_path"].notna()) &                     # Supprime les "na"
        (info_all_genome_df["ftp_path"].str.startswith("https"))       # Garde uniquement les URL HTTPS
    ].drop(columns=["assembly_level"]).reset_index(drop=True)
    print(df_filtered)
    df_filtered.to_csv(csv_file_out, sep="\t", index=False)

# csv_file_in = './assembly_summary.csv'
csv_file_out ="./assembly_summary_filtered_Complete_Genome.csv"
# filter_assembly_summary(csv_file_in, csv_file_out)
# 1ere ligne du fichier: GCF_900128725.1_BCifornacula_v1.0_genomic.fna
# >NZ_LT667500.1 Buchnera aphidicola strain BCifornacula voucher 2912 chromosome 1

info_genome_df = pd.read_csv(csv_file_out, sep="\t")
assembly_accession,https_path = info_genome_df.loc[0]
ftp_path = https_path[8:]
end_url_file = ftp_path.split('/')[-1]
# ftp_path = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/128/725/GCF_900128725.1_BCifornacula_v1.0"
real_ftp_path = f"{ftp_path}/{end_url_file}_genomic.fna.gz"
# 'ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/
genomes_dir = '.'
file_path = f"{genomes_dir}/{end_url_file}_genomic.fna.gz"
print("original filename")
os.system(f"wget -P {genomes_dir} {ftp_path}/{end_url_file}_genomic.fna.gz"
          f" && gzip -d {file_path}")

file_path = file_path.removesuffix(".gz")
accession_map = {}
with open(file_path, "r", encoding='utf-8') as reader:
    first_line = reader.readline()
    print(first_line)
    accession = first_line.split('.')[0][1:]  # Extraction de l'accession from >NZ_LT667500.1

# accession = 'NZ_LT667500'

yaml_file = "NCBI_credentials.yaml"


with open(yaml_file, 'r') as file:
    credentials = yaml.safe_load(file)

Entrez.email = credentials.get("email")  # Remplacez par votre email
Entrez.api_key = credentials.get("api_key") #""
Entrez.max_tries = 5
Entrez.sleep_between_tries = 15

with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as taxo_handle:
    x = SeqIO.read(taxo_handle, 'genbank')
    classif = x.annotations['taxonomy']
    sub = x.annotations['organism']

print(classif)
# From Siegfried
order: int | str = 0
for e in classif:
    if e[-4:] == 'ales':
        order = e
if order:
    group = classif[2] if classif[2][-4:] != 'ales' else classif[1]
    file_name: str = f"{genomes_dir}/{classif[0]}_{classif[1]}_{group}_{order}_{sub.split(' ')[0]}_{sub.split(' ')[1]}.fna"
    # os.system(f"mv {decompressed_path} {file_name}")
    print(f' change to  {file_name}")')
else:
    # we clean genomes we can't retrive classification for
    print(f' rm  {file_path}")')
    # os.system(f"rm {decompressed_path}")

# GET TAXO BY BATCH
# accessions = list(accession_map.keys())
# for i in tqdm(range(0, len(accessions), batch_size), desc="Fetching taxonomy data"):
#     batch = accessions[i:i + batch_size]
#     try:1
#         with Entrez.efetch(db="nucleotide", id=batch, rettype="gb", retmode="text") as taxo_handle:
#             # records = SeqIO.read(taxo_handle, 'genbank')
#             records = SeqI%ùO.par