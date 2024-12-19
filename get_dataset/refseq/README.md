# Downloading RefSeq bacteria

## OLD from Siegfried

1. Download the file at [https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt)
2. Use the script with the command `download_refseq.py assembly_summary.txt -o output_folder/ -e valid.email@mailbox.web`

Mail should be a valid email for Entrez. Output folder will contain the fasta gathered from the refseq.
The script solely downloads genomes flagged as "full genomes". Change [this line](https://github.com/dubssieg/wisp/blob/bd7493f85798d76426cad10148430c9b23383bf6/scripts/download_refseq.py#L22) to modify the behavior.


## NEW from branch optim_wisp

