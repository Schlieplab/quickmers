#!/bin/bash
#SBATCH --cpus-per-task=64           
#SBATCH --time=24:00:00                 
#SBATCH --nodes=1                       
#SBATCH --partition=comp                
#SBATCH --account=root                  
#SBATCH --job-name="test_edit"
#SBATCH --output="slurm_out/%x/slurm-%j.out"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
# ulimit -a
ml load python/3.12.4
source ~/venvs/test/bin/activate 
pip install -e .
ulimit -a
ulimit -t unlimited
ulimit -d unlimited
ulimit -u 8192
ulimit -c unlimited
ulimit -a
# python3 test_hamming.py
python3 test_edit.py