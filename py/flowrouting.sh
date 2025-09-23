#!/bin/bash
#SBATCH --account=def-gflowers
#SBATCH --time=0-24:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --output=flowrouting.out

source ~/SFU-code/antarctic-glads/venv/bin/activate
python -u flowrouting.py
