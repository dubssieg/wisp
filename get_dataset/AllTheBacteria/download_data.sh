wget https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/md5sum.txt
# Télécharger tous les fichiers listés dans md5sum.txt
cat md5sum.txt | awk '{print $2}' | while read file; do
  mkdir -p $(dirname "$file")  # Create the directory structure
  wget -P $(dirname "$file") "https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/$file"
done
