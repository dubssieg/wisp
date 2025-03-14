#!/bin/bash
#SBATCH --job-name=clx_create_db
#SBATCH --chdir=/home/genouest/cnrs_umr6074/cleroux2/microtaxo/wisp/wisp_allthebacteria
#SBATCH --output=slurm/%j.txt
#SBATCH --error=slurm/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=7-00:00:00 
# #SBATCH --mem-per-cpu=2G
#SBATCH --mem=60G
#SBATCH --mail-user=cyrille.leroux@irisa.fr
#SBATCH --mail-type=END,FAIL

echo "=== System Information ==="
uname -a

echo "=== CPU Information ==="
lscpu

echo "=== RAM Information ==="
free -h

echo "=== Disk Space ==="
df -h

echo "=== Environment Variables ==="
printenv | sort

echo "=== NFS Mounts ==="
mount | grep nfs

echo "=== Python Environment ==="
which python
python --version
pip list

echo "=== Running Processes ==="
ps aux --sort=-%mem | head -n 30

echo "=== Slurm Configuration ==="
scontrol show config

echo "=== Logs ==="
if [ -f "slurm.log" ]; then
    tail -n 20 slurm.log
else
    echo "No slurm.log found"
fi

source ~/.bashrc
conda activate wisp

python main.py --create-db

sacct -j $SLURM_JOB_ID --format=JobID,JobName%20,MaxRSS,State,ExitCode
