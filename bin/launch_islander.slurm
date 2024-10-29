#!/bin/bash
#SBATCH --job-name=juxtaposer_job      # Job name
#SBATCH --time=01:00:00                  # Time limit hrs:min:sec #YOU MAY NEED TO CHANGE THIS
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --mem=4G                         # Memory limit. YOU MAY NEED TO CHANGE THIS

source /home/${USER}/.bashrc
source activate prokka_env

cd /projects/wg-islands/Programs/tiger/islander_test; perl /projects/wg-islands/Programs/tiger/bin_ET/islander.pl -verbose $1 &> islander.log
