#!/bin/bash

## Script for converting bedgraphs to bigwigs
## Date: 16 June 2023

## Example usage:
## sbatch --array 0-64

# General settings
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=5:00:00
#SBATCH --mem=32GB

# Job name and output
#SBATCH -J bedGraphToBigWig
#SBATCH -o /Users/%u/slurm_out/slurm_%A_%a.out
#SBATCH -e /Users/%u/slurm_err/slurm_%A_%a.err

# Email notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=addz6536@colorado.edu

# Set variables
inDir=/scratch/Users/addz6536/t_kristen_mouseSP140/cutrun/d_macs2/all_treat_pileup.bdg
outDir=/scratch/Users/addz6536/t_kristen_mouseSP140/cutrun/e_bw
chromSizes=/Shares/CL_Shared/db/genomes/mm10/fa/mm10.chrom.sizes

# Define query files
queries=($(ls ${inDir}/*_pileup.bdg | xargs -n 1 basename | sed 's/_treat_pileup.bdg//g'))

# Convert bedgraphs to bigwigs
pwd; hostname; date

echo "Processing file: "${queries[$SLURM_ARRAY_TASK_ID]}
echo $(date +"[%b %d %H:%M:%S] Sorting bedgraph...")

sort \
-k1,1 -k2,2n \
${inDir}/${queries[$SLURM_ARRAY_TASK_ID]}_treat_pileup.bdg \
> ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg

echo $(date +"[%b %d %H:%M:%S] Converting sorted bedgraph to bigwig...")
echo "bedGraphToBigWig v 4"

/Shares/CL_Shared/programs/ucsc/bedGraphToBigWig \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg \
${chromSizes} \
${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.bw

echo $(date +"[%b %d %H:%M:%S] Removing sorted bedgraph...")

rm ${outDir}/${queries[$SLURM_ARRAY_TASK_ID]%.bdg}.sorted.bdg

echo $(date +"[%b %d %H:%M:%S] Done")
