#!/bin/bash
#SBATCH --job-name=clx_create_db
#SBATCH --chdir=/home/genouest/cnrs_umr6074/cleroux2/microtaxo/wisp/wisp_allthebacteria
#SBATCH --output=slurm/%j.txt
#SBATCH --error=slurm/%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00 
#SBATCH --mem=80G
#SBATCH --mail-user=cyrille.leroux@irisa.fr
#SBATCH --mail-type=END,FAIL

source ~/.bashrc
conda activate wisp


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

echo "=== Overcommit ==="
cat /proc/sys/vm/overcommit_memory
# echo 0 | sudo tee /proc/sys/vm/overcommit_memory

echo "=== Malloc ==="
ldd $(which python)
# LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libc.so.6 python main.py

echo "=== Swapiness ==="
cat /proc/sys/vm/swappiness
# echo 30 | sudo tee /proc/sys/vm/swappiness

export MALLOC_TRIM_THRESHOLD_=0
export MALLOC_MMAP_THRESHOLD_=131072
export MALLOC_ARENA_MAX=2 

# ulimit -Sv 40000000  # Soft limit 40 Go
# ulimit -Hv 42000000  # Hard limit 

# export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libc.so.6

python main.py --create-db

sacct -j $SLURM_JOB_ID --format=JobID,JobName%20,MaxRSS,State,ExitCode
