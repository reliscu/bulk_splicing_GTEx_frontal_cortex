#!/bin/bash
#SBATCH --job-name=interproscan
#SBATCH --output=interproscan_%j.log
#SBATCH --error=interproscan_%j.err
#SBATCH --account=fc_lareau
#SBATCH --qos=savio_normal
#SBATCH --time=72:00:00
#SBATCH -p savio2
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# ref: https://interproscan6.readthedocs.io/latest/

module load java
module load anaconda3
source activate protein_scan

export NXF_OPTS="-Xms2g -Xmx16g"

NXF_VER=25.10.4 /global/home/users/reliscu/scratch/programs/nextflow run ebi-pf-team/interproscan6 \
    -r 6.0.0 \
    -c licensed_apps.config \
    -profile singularity,slurm \
    --input proteins.fa \
    --datadir data/interpro_data \
    --applications NCBIFam,Pfam,TMbed,AntiFam,PirsR,SFLD,CDD,SMART,SUPERFAMILY,CATH-Gene3D,CATH-FunFam,CDD,MobiDB-lite,Coils,Phobius,PROSITE-patterns,PROSITE-profiles,signalp_euk,DeepTMHMM,PANTHER,PIRSF,PRINTS \
    --outdir data/results \
    --goterms \
    --maxWorkers 10 \
    --cpus 4
