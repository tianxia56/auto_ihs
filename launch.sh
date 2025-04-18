#!/bin/bash
#SBATCH --partition=ycga
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000
#SBATCH --job-name=ihs_auto


python run_pipeline.py --pop PJL --aaref arg --rmap pyrho --software hapbin --alt_na na --maf 05
