#!/bin/bash
#SBATCH --cpus-per-task=64           
#SBATCH --time=24:00:00                 
#SBATCH --nodes=1                       
#SBATCH --partition=comp                
#SBATCH --account=root                  
#SBATCH --job-name="test_quickmers"
#SBATCH --output="slurm_out/%x/slurm-%j.out"
ml load python/3.12.4
source ~/venvs/test/bin/activate 
pip install -e .
python3 test.py