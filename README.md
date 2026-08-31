# NANOResp
A pipeline for real time analysis of bacterial Nanopore metagenomics of clinical samples. The pipeline was designed for the analysis of respiratory samples and optimized for HPC clusters (slurm manager).

Scripts and configuration files used for manuscript XXXX

## Instructions:

1. Install che conda environments using the txt configuration files
2. Set all variables in the manager.sh script header. Most importantly, you can decide to run the script during live sequencing (basecalling not included) or on after sequencing and basecalling using the ontime tool. Please note: the variables are currently set to work on our cluster, please take some time to adapt it to yours.
4. Usage: `sbatch manager.sh Run1 [Run2 Run3 Run4]`
